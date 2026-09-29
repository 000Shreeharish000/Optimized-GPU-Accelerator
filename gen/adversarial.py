"""Adversarial LP/MILP suite: models on which weak implementations fail.

  klee_minty_N        exponential for Dantzig pricing (Klee & Minty 1972)
  scaled_<netlib>     Netlib model with rows/cols scaled by 10^U(-6,6) (ill-conditioned data)
  assign_N            N x N assignment LP: massively dual/primal degenerate
  nearsing_N          duplicated rows perturbed by 1e-12 (near-singular bases)
  transp_degen_N      degenerate transportation LP (all supplies/demands equal)
  infeas_box_N        infeasible by a thin margin (Farkas certificate needed)
  unbnd_N             unbounded with a feasible region (primal ray needed)
  bigm_knap_N         weak-relaxation MILP: big-M knapsack with fixed charges

Usage: python gen/adversarial.py --out data/adversarial
"""
from __future__ import annotations

import argparse
import gzip
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpswriter import INF, LinearModel  # noqa: E402


def klee_minty(n: int) -> LinearModel:
    # max sum 2^(n-j) x_j  s.t.  2 sum_{j<i} 2^(i-j) x_j + x_i <= 5^i
    m = LinearModel(f"KLEEMINTY_{n}", maximize=True)
    for j in range(1, n + 1):
        m.var(f"x{j}", 0, INF, cost=2.0 ** (n - j))
    for i in range(1, n + 1):
        coefs = {f"x{j}": 2.0 * 2 ** (i - j) for j in range(1, i)}
        coefs[f"x{i}"] = 1.0
        m.row(f"c{i}", coefs, -INF, 5.0 ** i)
    return m


def assignment(n: int, seed: int) -> LinearModel:
    rng = random.Random(seed)
    m = LinearModel(f"ASSIGN_{n}")
    for i in range(n):
        for j in range(n):
            m.var(f"x{i}_{j}", 0, INF, cost=float(rng.randint(1, 5)))  # many ties -> degeneracy
    for i in range(n):
        m.row(f"r{i}", {f"x{i}_{j}": 1 for j in range(n)}, 1, 1)
    for j in range(n):
        m.row(f"c{j}", {f"x{i}_{j}": 1 for i in range(n)}, 1, 1)
    return m


def transport_degenerate(n: int, seed: int) -> LinearModel:
    rng = random.Random(seed)
    m = LinearModel(f"TRANSPDEG_{n}")
    for i in range(n):
        for j in range(n):
            m.var(f"x{i}_{j}", 0, INF, cost=float(rng.randint(1, 9)))
    for i in range(n):
        m.row(f"s{i}", {f"x{i}_{j}": 1 for j in range(n)}, -INF, 10)
    for j in range(n):
        m.row(f"d{j}", {f"x{i}_{j}": 1 for i in range(n)}, 10, INF)
    return m


def near_singular(n: int, seed: int) -> LinearModel:
    rng = random.Random(seed)
    m = LinearModel(f"NEARSING_{n}")
    for j in range(n):
        m.var(f"x{j}", 0, 10, cost=rng.uniform(-1, 1))
    for i in range(n // 2):
        coefs = {f"x{j}": rng.uniform(-1, 1) for j in range(n) if rng.random() < 0.3}
        if not coefs:
            coefs = {"x0": 1.0}
        m.row(f"r{i}", coefs, -INF, rng.uniform(1, 5))
        dup = {k: v * (1 + 1e-12 * rng.uniform(-1, 1)) for k, v in coefs.items()}
        m.row(f"d{i}", dup, -INF, m.rows[f"r{i}"]["ub"] * (1 + 1e-12))
    return m


def infeasible_thin(n: int) -> LinearModel:
    # sum x >= n + 1e-6 with x in [0,1]: infeasible by a tiny margin
    m = LinearModel(f"INFEASTHIN_{n}")
    for j in range(n):
        m.var(f"x{j}", 0, 1, cost=1.0)
    m.row("cover", {f"x{j}": 1 for j in range(n)}, n + 1e-6, INF)
    for j in range(0, n - 1, 2):
        m.row(f"pair{j}", {f"x{j}": 1, f"x{j + 1}": -1}, -1, 1)
    return m


def unbounded(n: int) -> LinearModel:
    m = LinearModel(f"UNBND_{n}")
    for j in range(n):
        m.var(f"x{j}", 0, INF, cost=-1.0 if j == n - 1 else 1.0)
    for j in range(n - 1):
        m.row(f"r{j}", {f"x{j}": 1, f"x{j + 1}": -1}, -INF, 1)  # x_{j} - x_{j+1} <= 1
    return m


def bigm_knapsack(n: int, seed: int) -> LinearModel:
    rng = random.Random(seed)
    m = LinearModel(f"BIGMKNAP_{n}")
    for j in range(n):
        m.var(f"y{j}", 0, 1, cost=rng.uniform(50, 100), integer=True)
        m.var(f"x{j}", 0, INF, cost=rng.uniform(1, 3))
        m.row(f"link{j}", {f"x{j}": 1, f"y{j}": -1e4}, -INF, 0)  # big-M: weak relaxation
    m.row("demand", {f"x{j}": 1 for j in range(n)}, 37.5 * n, INF)
    for j in range(n):
        m.row(f"cap{j}", {f"x{j}": 1}, -INF, 100)
    return m


def scaled_copy(src_path: str, out_path: str, seed: int) -> None:
    """Scale rows/cols of an MPS model by 10^U(-6,6): same optimum after unscaling, terrible data."""
    rng = random.Random(seed)
    opener = gzip.open if src_path.endswith(".gz") else open
    with opener(src_path, "rt") as fh:
        lines = fh.read().splitlines()
    rscale, cscale = {}, {}

    def rs(r):
        if r not in rscale:
            rscale[r] = 10.0 ** rng.uniform(-6, 6)
        return rscale[r]

    def cs(c):
        if c not in cscale:
            cscale[c] = 10.0 ** rng.uniform(-6, 6)
        return cscale[c]

    out, sec, obj = [], None, None
    for ln in lines:
        if not ln.strip() or ln.startswith("*"):
            continue
        if not ln[0].isspace():
            sec = ln.split()[0]
            out.append(ln)
            continue
        t = ln.split()
        if sec == "ROWS":
            if t[0] == "N" and obj is None:
                obj = t[1]
                rscale[obj] = 1.0
            out.append(ln)
        elif sec == "COLUMNS":
            if "MARKER" in ln or len(t) not in (3, 5):
                out.append(ln)
                continue
            c = t[0]
            parts = [c]
            for k in range(1, len(t), 2):
                r, v = t[k], float(t[k + 1])
                parts += [r, repr(v * rs(r) * cs(c))]
            out.append("    " + "  ".join(parts))
        elif sec in ("RHS", "RANGES"):
            if len(t) not in (2, 3, 4, 5):
                out.append(ln)
                continue
            first = 1 if len(t) in (3, 5) else 0  # set name may be omitted
            parts = [t[0]] if first else ["RHS"]
            for k in range(first, len(t), 2):
                r, v = t[k], float(t[k + 1])
                parts += [r, repr(v * rs(r))]
            out.append("    " + "  ".join(parts))
        elif sec == "BOUNDS":
            if len(t) == 4:
                c = t[2]
                out.append(f" {t[0]} {t[1]} {c} {float(t[3]) / cs(c)!r}")
            else:
                out.append(ln)
        else:
            out.append(ln)
    with open(out_path, "w") as fh:
        fh.write("\n".join(out) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="data/adversarial")
    ap.add_argument("--netlib", default="data/netlib")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    models = [klee_minty(12), klee_minty(20), assignment(30, a.seed), assignment(80, a.seed),
              transport_degenerate(40, a.seed), near_singular(200, a.seed), infeasible_thin(50),
              unbounded(30), bigm_knapsack(25, a.seed)]
    for m in models:
        p = os.path.join(a.out, m.name.lower() + ".mps")
        m.write(p)
        print(f"{m.name}: {m.num_rows} rows, {m.num_cols} cols -> {p}")
    for nm in ["afiro", "sc50a", "adlittle", "blend", "share1b", "scfxm1", "bandm", "e226"]:
        src = os.path.join(a.netlib, nm + ".mps.gz")
        if os.path.exists(src):
            dst = os.path.join(a.out, f"scaled_{nm}.mps")
            scaled_copy(src, dst, a.seed)
            print(f"scaled {nm} -> {dst}")


if __name__ == "__main__":
    main()
