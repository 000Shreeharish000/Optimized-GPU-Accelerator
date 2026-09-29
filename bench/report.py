"""Assemble docs/BENCHMARKS.md (+ figures) from results/*. Nothing is filtered:
every per-instance table (including all failures) stays in results/<set>/summary.md
and is linked from here."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
FIG = ROOT / "docs" / "figures"
PROVEN = {"OPTIMAL", "INFEASIBLE", "UNBOUNDED"}


def rows(set_name):
    p = RES / set_name / "runs.csv"
    if not p.exists():
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh))


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def sgm(ts, shift=10.0):
    return math.exp(sum(math.log(t + shift) for t in ts) / len(ts)) - shift if ts else float("nan")


def is_mip(set_name, inst, solver):
    f = RES / set_name / "json" / f"{inst}.{solver}.json"
    try:
        return json.loads(f.read_text())["telemetry"]["model"].get("integers", 0) > 0
    except Exception:  # noqa: BLE001
        return set_name.startswith("miplib")


def table(set_name, solvers, limit):
    rs = rows(set_name)
    by = {}
    for r in rs:
        by.setdefault(r["instance"], {})[r["solver"]] = r
    out = ["| solver | proven | certified | SGM s (shift 10) | disagreements with HiGHS (1e-6 LP, 1e-4 MIP gap) |", "|---|---|---|---|---|"]
    for s in solvers:
        solved = cert = dis = 0
        ts = []
        for inst, d in by.items():
            r = d.get(s)
            if not r:
                continue
            ok = r["status"] in PROVEN
            solved += ok
            cert += ok and r.get("certificate") == "ACCEPTED"
            ts.append(min(num(r["time"]) or limit, limit) if ok else limit)
            h = d.get("highs")
            if s != "highs" and h and ok and h["status"] == "OPTIMAL" and r["status"] == "OPTIMAL":
                a, b = num(r["objective"]), num(h["objective"])
                tol = 1e-4 if is_mip(set_name, inst, s) else 1e-6  # both solvers' default MIP gap
                if a is not None and b is not None and abs(a - b) > tol * max(1, abs(b)):
                    dis += 1
        if ts:
            out.append(f"| {s} | {solved}/{len(ts)} | {cert if s != 'highs' else '-'} | {sgm(ts):.3f} | "
                       f"{dis if s != 'highs' else '-'} |")
    return out, by


def adversarial_truth(by):
    """Ground truth: scaled_<x> has the optimum of Netlib <x> (exact rescaling)."""
    net = {r["instance"]: r for r in rows("netlib") if r["solver"] == "dual"}
    truth = {}
    for inst in by:
        if inst.startswith("scaled_") and inst[7:] in net:
            truth[inst] = ("OPTIMAL", num(net[inst[7:]]["objective"]))
        elif inst.startswith("kleeminty_"):
            truth[inst] = ("OPTIMAL", 5.0 ** int(inst.split("_")[1]))
        elif inst.startswith("infeasthin"):
            truth[inst] = ("INFEASIBLE", None)
        elif inst.startswith("unbnd"):
            truth[inst] = ("UNBOUNDED", None)
    lines = ["| instance | truth | PRAMANA | HiGHS (scipy) |", "|---|---|---|---|"]
    score = {"auto": [0, 0, 0], "highs": [0, 0, 0]}  # correct, wrong, no-answer
    for inst in sorted(truth):
        st, obj = truth[inst]
        cells = []
        for s in ("auto", "highs"):
            r = by[inst].get(s, {})
            rs, ro = r.get("status"), num(r.get("objective"))
            if rs == st and (obj is None or (ro is not None and abs(ro - obj) <= 1e-6 * max(1, abs(obj)))):
                verdict = "correct"
                score[s][0] += 1
            elif rs in PROVEN:
                verdict = "**WRONG**"
                score[s][1] += 1
            else:
                verdict = "no claim"
                score[s][2] += 1
            cells.append(f"{rs} {'' if ro is None else f'{ro:.8g}'} ({verdict})")
        lines.append(f"| {inst} | {st} {'' if obj is None else f'{obj:.8g}'} | {cells[0]} | {cells[1]} |")
    lines += ["", f"Score (correct / WRONG / no claim): PRAMANA {score['auto']}, HiGHS {score['highs']}."]
    return lines


def perf_profile(by, solvers, limit, path, title):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # noqa: BLE001
        return
    ratios = {s: [] for s in solvers}
    for inst, d in by.items():
        ts = {}
        for s in solvers:
            r = d.get(s)
            ok = r and r["status"] in PROVEN
            ts[s] = max(num(r["time"]) or 1e-4, 1e-4) if ok else None
        best = min([t for t in ts.values() if t is not None], default=None)
        for s in solvers:
            ratios[s].append(ts[s] / best if (ts[s] is not None and best) else float("inf"))
    fig, ax = plt.subplots(figsize=(6.4, 4))
    xs = [1 * (1.25 ** k) for k in range(60)]
    for s in solvers:
        r = ratios[s]
        ax.step(xs, [sum(1 for v in r if v <= x) / len(r) for x in xs], where="post", label=s)
    ax.set_xscale("log")
    ax.set_xlabel("time ratio to best solver (tau)")
    ax.set_ylabel("fraction of instances solved within tau")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=130)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    md = ["# PRAMANA benchmark report", "",
          "Generated by `bench/report.py` from `results/*/runs.csv` (raw per-run JSON in `results/*/json`).",
          "Every per-instance table, including every failure, is in `results/<set>/summary.md`.",
          "Machine: see `results/machine.json`. Reference solver: HiGHS through scipy (same machine, same limits).",
          "PRAMANA times include presolve, postsolve, crossover and certification; file reading excluded for both.", ""]
    sections = [("netlib", ["dual", "highs"], 300, "Netlib LP (all 91 feasible)"),
                ("netlib-infeas", ["dual", "highs"], 300, "Netlib infeasible LPs (all 29)"),
                ("netlib_engines", ["ipm", "pdhg-cpu", "pdhg-gpu"], 30, "Netlib LP: other PRAMANA engines (30 s)"),
                ("gen", ["dual", "ipm", "pdhg-cpu", "pdhg-gpu", "highs"], 120, "Industrial suite (refinery planning / blending / scheduling / UC / dispatch)"),
                ("miplib3", ["auto", "highs"], 60, "MIPLIB 3 pre-registered subset (60 s)"),
                ("miplib3_nocuts", ["auto"], 60, "MIPLIB 3 subset WITHOUT cutting planes (ablation)"),
                ("maros", ["auto"], 60, "Maros-Meszaros convex QP (124 instances, 60 s)"),
                ("adversarial", ["auto", "highs"], 60, "Adversarial suite")]
    for name, solvers, lim, title in sections:
        t, by = table(name, solvers, lim)
        if not by:
            continue
        md += [f"## {title}", ""] + t + ["", f"Full table: [results/{name}/summary.md](../results/{name}/summary.md)", ""]
        if name == "adversarial":
            md += ["### Adversarial ground truth (scaled models have the SAME optimum as the original Netlib model)", ""]
            md += adversarial_truth(by) + [""]
        if name in ("netlib", "miplib3", "gen"):
            perf_profile(by, [s for s in solvers], lim, FIG / f"profile_{name}.png", f"Performance profile: {title}")
            md += [f"![profile](figures/profile_{name}.png)", ""]
    # MIP cuts ablation
    a = {r["instance"]: r for r in rows("miplib3") if r["solver"] == "auto"}
    b = {r["instance"]: r for r in rows("miplib3_nocuts") if r["solver"] == "auto"}
    if a and b:
        md += ["## Ablation: cutting planes (GMI + MIR + cover) on/off", "",
               "| instance | with cuts: status / nodes / s | without cuts: status / nodes / s |", "|---|---|---|"]
        for k in sorted(a):
            if k in b:
                md.append(f"| {k} | {a[k]['status']} / {a[k].get('nodes', '')} / {num(a[k]['time']) or 0:.2f} | "
                          f"{b[k]['status']} / {b[k].get('nodes', '')} / {num(b[k]['time']) or 0:.2f} |")
        md.append("")
    for extra, title in [("gpu", "CPU vs GPU"), ("family", "Repeated-solve family"), ("router", "Router")]:
        p = RES / extra / "summary.md"
        if p.exists():
            md += [p.read_text(), ""]
    cal = RES / "gpu" / "calibrate.json"
    if cal.exists():
        pts = json.loads(cal.read_text()).get("pdhg_iteration_timing", [])
        md += ["## PDHG per-iteration cost: CPU (all threads) vs GPU (own kernels)", "",
               "`pramana calibrate`: one PDHG iteration (2 SpMV + fused vector updates) on synthetic LPs.", "",
               "| nnz | CPU ms/iter | GPU ms/iter | GPU speedup |", "|---|---|---|---|"]
        for p in pts:
            c, gg = p["cpu"]["seconds_per_iteration"], p.get("gpu", {}).get("seconds_per_iteration")
            md.append(f"| {p['nnz']:,} | {c * 1e3:.3f} | {'' if gg is None else f'{gg * 1e3:.3f}'} | "
                      f"{'' if not gg else f'{c / gg:.2f}x'} |")
        md.append("")
    gp = RES / "gpu" / "planning_engines.csv"
    if gp.exists():
        md += ["## CPU vs GPU on the refinery planning family (end-to-end seconds, 300 s limit)", "",
               "Columns dual..pdhg-gpu: exactly certified optimum (PDHG followed by crossover). Raw columns: first-order",
               "answer at 1e-4 relative KKT *without* crossover. Its primal point is NOT feasible to 1e-6 (so it is never",
               "reported as OPTIMAL), but its duals give a rigorous Neumaier-Shcherbina bound; brackets show that bound's",
               "relative distance from the true optimum (reference optimum: any certified run of the same model, else HiGHS).", "",
               "| model | nnz | dual | ipm | pdhg-cpu | pdhg-gpu | raw pdhg-cpu 1e-4 [bound err] | raw pdhg-gpu 1e-4 [bound err] | HiGHS |",
               "|---|---|---|---|---|---|---|---|---|"]
        highs = {r["instance"]: r for r in rows("gen") if r["solver"] == "highs"}
        with open(gp) as fh:
            g = list(csv.DictReader(fh))
        models = sorted({r["model"] for r in g}, key=lambda m: int(next((r["nnz"] for r in g if r["model"] == m and r["nnz"]), 0) or 0))
        for m in models:
            h = highs.get(m, {})
            opt = next((num(r["objective"]) for r in g if r["model"] == m and r.get("status") == "OPTIMAL"), None)
            if opt is None and h.get("status") == "OPTIMAL":
                opt = num(h.get("objective"))
            cell = {}
            nnz = ""
            for r in g:
                if r["model"] == m:
                    nnz = r.get("nnz") or nnz
                    s = r.get("status", "")
                    if "raw" in r["config"]:
                        bound = None
                        try:
                            j = json.loads((RES / "gpu" / f"{m}.{r['config']}.json").read_text())
                            bound = num(j["certificate"]["safe_dual_bound"])
                        except Exception:  # noqa: BLE001
                            pass
                        err = abs(bound - opt) / max(1.0, abs(opt)) if (bound is not None and opt is not None) else None
                        cell[r["config"]] = f"{num(r['seconds']) or 0:.2f}" + (f" [{err:.1e}]" if err is not None else " [-]")
                    else:
                        cell[r["config"]] = f"{num(r['seconds']) or 0:.2f}" + ("" if s == "OPTIMAL" else f" ({s})")
            hc = f"{num(h.get('time')) or 0:.2f}" + ("" if h.get("status") == "OPTIMAL" else f" ({h.get('status')})") if h else ""
            md.append(f"| {m} | {nnz} | {cell.get('dual', '')} | {cell.get('ipm', '')} | {cell.get('pdhg-cpu', '')} | "
                      f"{cell.get('pdhg-gpu', '')} | {cell.get('pdhg-cpu_raw1e-4', '')} | {cell.get('pdhg-gpu_raw1e-4', '')} | {hc} |")
        md += ["", "Note: the HiGHS column comes from the `gen` benchmark run (120 s limit); timings from different",
               "runs on this laptop vary by up to ~2x (thermal/power state), so compare within a column group."]
        md += ["", "![iteration crossover](figures/gpu_iteration_crossover.png)", "",
               "![planning engines](figures/gpu_planning_engines.png)", ""]
    (ROOT / "docs" / "BENCHMARKS.md").write_text("\n".join(md) + "\n")
    print("wrote docs/BENCHMARKS.md")


if __name__ == "__main__":
    main()
