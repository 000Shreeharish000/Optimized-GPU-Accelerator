"""PRAMANA benchmark harness (reproducible, full sets, failure tables).

For every instance of a set it runs the PRAMANA CLI (one process per instance,
JSON output) for each requested engine and, optionally, the reference solver
HiGHS through scipy (LP: scipy.optimize.linprog(method="highs"), MIP:
scipy.optimize.milp). HiGHS is used ONLY here, never linked into the solver.

Outputs (results/<set>/):
  runs.csv        one row per (instance, solver): status, objective, time, certificate, gap, engine
  summary.md      solved counts, shifted geometric means (shift 10 s), agreement with the
                  reference, and the FULL failure table (nothing filtered)

Inclusion rules are fixed in bench/sets/*.txt (committed before runs).

Usage:
  python bench/run_bench.py --set netlib --engines dual --reference --time 300
  python bench/run_bench.py --set gen --engines auto,dual,ipm,pdhg-cpu,pdhg-gpu --time 120
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "python"))
BUILD_EXE = ROOT / "build" / ("pramana.exe" if os.name == "nt" else "pramana")
EXE = BUILD_EXE  # replaced by a private copy in main() so rebuilding never races a benchmark

SETS = {
    "netlib": "data/netlib/*.mps.gz",
    "netlib-infeas": "data/netlib_infeas/*.mps",
    "miplib3": "bench/sets/miplib3_subset.txt",
    "maros": "data/maros/*.qps",
    "gen": "data/gen/*",
    "adversarial": "data/adversarial/*.mps",
}


def instances(set_name: str):
    spec = SETS[set_name]
    if spec.endswith(".txt"):
        names = [l.split()[0] for l in (ROOT / spec).read_text().splitlines() if l.strip() and not l.startswith("#")]
        return [str(ROOT / "data" / "miplib3" / f"{n}.mps.gz") for n in names]
    return sorted(glob.glob(str(ROOT / spec)))


def run_pramana(path: str, engine: str, tlimit: float, workdir: Path, extra=()):
    base = Path(path).name.split(".")[0]
    out = workdir / f"{base}.{engine}.json"
    cmd = [str(EXE), path, "--algo", engine, "--time", str(tlimit), "--log", "0", "--json", str(out)] + list(extra)
    t0 = time.perf_counter()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=tlimit * 2 + 60)
        wall = time.perf_counter() - t0
    except subprocess.TimeoutExpired:
        return {"status": "HARD_TIMEOUT", "time": tlimit * 2 + 60, "engine": engine}
    if not out.exists():
        return {"status": "ERROR", "time": wall, "engine": engine, "detail": (p.stderr or p.stdout)[-300:]}
    j = json.loads(out.read_text())
    cert = j.get("certificate", {})
    tel = j.get("telemetry", {})
    r = {
        "status": j.get("status"),
        "engine_status": j.get("engine_status"),
        "engine": j.get("engine", engine),
        "objective": j.get("objective"),
        "bound": j.get("best_bound"),
        "gap": j.get("relative_gap"),
        "time": j.get("seconds"),
        "wall": wall,
        "certificate": "ACCEPTED" if cert.get("accepted") else "REJECTED",
        "cert_gap": cert.get("certified_relative_gap"),
        "rows": tel.get("model", {}).get("rows"),
        "cols": tel.get("model", {}).get("cols"),
        "nnz": tel.get("model", {}).get("nnz"),
        "router_engine": tel.get("router", {}).get("selected"),
    }
    for k in ("simplex", "ipm", "pdhg", "mip"):
        if k in tel:
            t = tel[k]
            r["iterations"] = t.get("iterations", t.get("nodes"))
            if k == "pdhg":
                r["gpu_kernel_s"] = t.get("kernel_seconds")
                r["gpu_transfer_s"] = t.get("transfer_seconds")
            if k == "mip":
                r["nodes"] = t.get("nodes")
    return r


def run_highs(path: str, tlimit: float):
    import numpy as np
    import scipy.sparse as sp
    from scipy.optimize import Bounds, LinearConstraint, linprog, milp

    from pramana.verify import parse_mps

    M = parse_mps(path, exact=False)
    if M.q:
        return {"status": "N/A", "detail": "QP not supported by scipy.optimize"}
    n, m = M.n, M.m
    sense = -1.0 if M.maximize else 1.0
    c = np.array([sense * float(v) for v in M.cost])
    rows, cols, vals = [], [], []
    for j, col in enumerate(M.cols):
        for i, a in col.items():
            rows.append(i)
            cols.append(j)
            vals.append(float(a))
    A = sp.csr_matrix((vals, (rows, cols)), shape=(m, n))
    inf = np.inf
    rl = np.array([-inf if v is None else float(v) for v in M.rlo])
    ru = np.array([inf if v is None else float(v) for v in M.rup])
    lb = np.array([-inf if v is None else float(v) for v in M.clo])
    ub = np.array([inf if v is None else float(v) for v in M.cup])
    offset = float(M.offset)
    t0 = time.perf_counter()
    if any(M.integer):
        res = milp(c, integrality=np.array([1 if b else 0 for b in M.integer]), bounds=Bounds(lb, ub),
                   constraints=LinearConstraint(A, rl, ru) if m else (), options={"time_limit": tlimit, "disp": False})
        dt = time.perf_counter() - t0
        st = {0: "OPTIMAL", 1: "TIME_LIMIT", 2: "INFEASIBLE", 3: "UNBOUNDED"}.get(res.status, f"STATUS{res.status}")
        obj = None if res.x is None else sense * (float(res.fun) + sense * offset)
        return {"status": st, "objective": obj, "time": dt, "engine": "highs-mip",
                "bound": None if getattr(res, "mip_dual_bound", None) is None else sense * (res.mip_dual_bound + sense * offset)}
    eq = rl == ru
    ubr = np.isfinite(ru) & ~eq
    lbr = np.isfinite(rl) & ~eq
    A_ub = sp.vstack([A[ubr], -A[lbr]]) if (ubr.any() or lbr.any()) else None
    b_ub = np.concatenate([ru[ubr], -rl[lbr]]) if A_ub is not None else None
    A_eq = A[eq] if eq.any() else None
    b_eq = rl[eq] if eq.any() else None
    bounds = list(zip([None if not np.isfinite(v) else v for v in lb], [None if not np.isfinite(v) else v for v in ub]))
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs",
                  options={"time_limit": tlimit, "presolve": True})
    dt = time.perf_counter() - t0
    st = {0: "OPTIMAL", 1: "TIME_LIMIT", 2: "INFEASIBLE", 3: "UNBOUNDED"}.get(res.status, f"STATUS{res.status}")
    obj = None if res.x is None or res.status != 0 else sense * (float(res.fun) + sense * offset)
    return {"status": st, "objective": obj, "time": dt, "engine": "highs"}


def sgm(times, shift=10.0):
    if not times:
        return float("nan")
    return math.exp(sum(math.log(t + shift) for t in times) / len(times)) - shift


def summarize(set_name, rows, solvers, tlimit, outdir: Path):
    by = {}
    for r in rows:
        by.setdefault(r["instance"], {})[r["solver"]] = r
    lines = [f"# Benchmark: {set_name}", "",
             f"Instances: {len(by)}; time limit {tlimit:g}s; machine: see results/machine.json.",
             "Status = certified status (PRAMANA) / solver status (HiGHS). Time excludes file reading.", ""]
    proven = {"OPTIMAL", "INFEASIBLE", "UNBOUNDED"}
    lines += ["| solver | solved (proven) | certified | SGM time (s, shift 10, unsolved at limit) | wrong vs reference |",
              "|---|---|---|---|---|"]
    ref = "highs" if "highs" in solvers else None
    for s in solvers:
        solved = cert = wrong = 0
        ts = []
        for inst, d in by.items():
            r = d.get(s)
            if not r:
                continue
            ok = r.get("status") in proven
            solved += ok
            cert += r.get("certificate") == "ACCEPTED" and ok
            t = r.get("time") if ok and r.get("time") is not None else tlimit
            ts.append(min(float(t), tlimit))
            if ref and s != ref and ok and ref in d and d[ref].get("status") == "OPTIMAL" and r.get("status") == "OPTIMAL":
                a, b = r.get("objective"), d[ref].get("objective")
                tol = 1e-4 if r.get("nodes") not in (None, "") else 1e-6  # MIP: default relative gap
                if a is not None and b is not None and abs(a - b) > tol * max(1, abs(b)):
                    wrong += 1
        lines.append(f"| {s} | {solved}/{len(by)} | {cert if s != 'highs' else '-'} | {sgm(ts):.3f} | {wrong if s != ref else '-'} |")
    lines += ["", "## Full per-instance table (including every failure)", ""]
    hdr = "| instance | " + " | ".join(f"{s} status | {s} obj | {s} time" for s in solvers) + " |"
    lines += [hdr, "|" + "---|" * (1 + 3 * len(solvers))]
    for inst in sorted(by):
        cells = []
        for s in solvers:
            r = by[inst].get(s, {})
            obj = r.get("objective")
            cells += [str(r.get("status", "")), "" if obj is None else f"{obj:.10g}",
                      "" if r.get("time") is None else f"{float(r['time']):.3f}"]
        lines.append(f"| {inst} | " + " | ".join(cells) + " |")
    lines += ["", "## Failures / disagreements", ""]
    nf = 0
    for inst in sorted(by):
        for s in solvers:
            r = by[inst].get(s, {})
            if r.get("status") not in proven:
                lines.append(f"- {inst} / {s}: {r.get('status')} {r.get('detail', '')}")
                nf += 1
    if nf == 0:
        lines.append("- none")
    (outdir / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:12]))


def machine_info():
    import platform
    info = {"platform": platform.platform(), "processor": platform.processor(), "python": platform.python_version(),
            "cpu_count": os.cpu_count()}
    try:
        info["pramana_info"] = subprocess.run([str(EXE), "info"], capture_output=True, text=True, timeout=120).stdout
    except Exception as e:  # noqa: BLE001
        info["pramana_info"] = str(e)
    return info


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", required=True, choices=list(SETS))
    ap.add_argument("--engines", default="auto")
    ap.add_argument("--reference", action="store_true", help="also run HiGHS (scipy) as reference")
    ap.add_argument("--time", type=float, default=60)
    ap.add_argument("--limit", type=int, default=0, help="first N instances only (smoke runs)")
    ap.add_argument("--extra", default="", help="extra CLI flags for PRAMANA, e.g. '--no-cuts'")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    global EXE
    import shutil
    outdir = ROOT / "results" / (a.set + (f"_{a.tag}" if a.tag else ""))
    (ROOT / "results" / "bin").mkdir(parents=True, exist_ok=True)
    EXE = ROOT / "results" / "bin" / f"pramana_{a.set}{('_' + a.tag) if a.tag else ''}{BUILD_EXE.suffix}"
    shutil.copy2(BUILD_EXE, EXE)
    work = outdir / "json"
    work.mkdir(parents=True, exist_ok=True)
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "machine.json").write_text(json.dumps(machine_info(), indent=1))
    insts = instances(a.set)
    if a.limit:
        insts = insts[: a.limit]
    engines = [e for e in a.engines.split(",") if e]
    solvers = engines + (["highs"] if a.reference else [])
    rows = []
    csv_path = outdir / "runs.csv"
    fields = ["instance", "solver", "status", "engine_status", "engine", "objective", "bound", "gap", "time", "wall",
              "certificate", "cert_gap", "rows", "cols", "nnz", "iterations", "nodes", "gpu_kernel_s", "gpu_transfer_s",
              "router_engine", "detail"]
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for path in insts:
            inst = Path(path).name.split(".")[0]
            for e in engines:
                r = run_pramana(path, e, a.time, work, a.extra.split() if a.extra else ())
                r.update({"instance": inst, "solver": e})
                rows.append(r)
                w.writerow(r)
                fh.flush()
                print(f"{inst:24s} {e:9s} {str(r.get('status')):18s} {r.get('objective')!s:>22s} {r.get('time')!s:>10s}",
                      flush=True)
            if a.reference:
                try:
                    r = run_highs(path, a.time)
                except Exception as ex:  # noqa: BLE001
                    r = {"status": "ERROR", "detail": str(ex)[:200]}
                r.update({"instance": inst, "solver": "highs"})
                rows.append(r)
                w.writerow(r)
                fh.flush()
                print(f"{inst:24s} {'highs':9s} {str(r.get('status')):18s} {r.get('objective')!s:>22s} {r.get('time')!s:>10s}",
                      flush=True)
    summarize(a.set, rows, solvers, a.time, outdir)


if __name__ == "__main__":
    main()
