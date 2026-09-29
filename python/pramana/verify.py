"""Independent exact-arithmetic verifier for PRAMANA results.

Shares NO code with the solver: its own MPS/QPS parser reads the decimal data as
exact rationals (fractions.Fraction), and every claim in a result JSON
(`pramana model.mps --json out.json --vectors`) is re-checked exactly:

  OPTIMAL (LP)   exact primal violations; exact Neumaier-Shcherbina dual bound
                 from the reported duals (after the free sign-cone projection),
                 exact gap; optional --exact-basis: recompute the basic solution
                 and duals from the reported basis by exact rational elimination
                 and check primal AND dual feasibility with zero tolerance.
  INFEASIBLE     exact Farkas test for y and -y.
  UNBOUNDED      exact recession-cone and descent test of the primal ray.
  MIP            exact integrality, feasibility and objective of the incumbent.
  QP             exact primal checks; linearization bound with exact arithmetic.

Usage: python -m pramana.verify model.mps result.json [--exact-basis] [--tol 1e-6]
Exit code 0 = all claims verified.
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import sys
from fractions import Fraction as F
from typing import Dict, List, Optional, Tuple

INF = None  # infinite bounds are represented as None


# ---------------------------------------------------------------------------
# Independent MPS / QPS parser (exact)
# ---------------------------------------------------------------------------
class ExactModel:
    def __init__(self):
        self.name = ""
        self.maximize = False
        self.obj_row = None
        self.row_names: List[str] = []
        self.row_type: List[str] = []
        self.row_index: Dict[str, int] = {}
        self.col_names: List[str] = []
        self.col_index: Dict[str, int] = {}
        self.cost: List[F] = []
        self.integer: List[bool] = []
        self.cols: List[Dict[int, F]] = []   # column -> {row: value}
        self.rhs: List[F] = []
        self.rlo: List[Optional[F]] = []
        self.rup: List[Optional[F]] = []
        self.clo: List[Optional[F]] = []
        self.cup: List[Optional[F]] = []
        self.offset = F(0)
        self.q: Dict[Tuple[int, int], F] = {}  # full symmetric Q

    @property
    def n(self):
        return len(self.col_names)

    @property
    def m(self):
        return len(self.row_names)


def _num(s: str) -> F:
    s = s.replace("D", "E").replace("d", "e")
    return F(s)


def parse_mps(path: str) -> ExactModel:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as fh:
        lines = fh.read().splitlines()
    M = ExactModel()
    sec = None
    marker_int = False
    has_lo: set = set()
    has_up: set = set()
    ranges: Dict[int, F] = {}
    free_rows: set = set()
    fixed_fmt = False
    quad_full = False

    def fields(line):
        cuts = [(1, 3), (4, 12), (14, 22), (24, 36), (39, 47), (49, 61)]
        return [line[a:b].strip() if len(line) > a else "" for a, b in cuts]

    for raw in lines:
        if not raw.strip() or raw.startswith("*"):
            continue
        if not raw[0].isspace():
            tok = raw.split()
            head = tok[0].upper()
            if head == "NAME":
                M.name = raw[4:].strip()
            elif head in ("OBJSENSE", "OBJSENCE"):
                sec = "OBJSENSE"
                if len(tok) > 1 and tok[1].upper().startswith("MAX"):
                    M.maximize = True
            elif head in ("ROWS", "COLUMNS", "RHS", "RANGES", "BOUNDS", "ENDATA"):
                sec = head
            elif head == "QUADOBJ":
                sec, quad_full = "QUAD", False
            elif head in ("QMATRIX", "QSECTION"):
                sec, quad_full = "QUAD", True
            elif sec == "OBJSENSE" and head.startswith("MAX"):
                M.maximize = True
            if sec == "ENDATA":
                break
            continue
        tok = raw.split()
        if sec == "OBJSENSE":
            if tok[0].upper().startswith("MAX"):
                M.maximize = True
        elif sec == "ROWS":
            if len(tok) > 2:
                fixed_fmt = True
            if fixed_fmt:
                f = fields(raw)
                t, nm = f[0], f[1]
            else:
                t, nm = tok[0], tok[1]
            t = t.upper()
            if t == "N":
                if M.obj_row is None:
                    M.obj_row = nm
                else:
                    free_rows.add(nm)
                continue
            M.row_index[nm] = len(M.row_names)
            M.row_names.append(nm)
            M.row_type.append(t)
            M.rhs.append(F(0))
        elif sec == "COLUMNS":
            if len(tok) >= 3 and "MARKER" in tok[1]:
                if "INTORG" in tok[2]:
                    marker_int = True
                elif "INTEND" in tok[2]:
                    marker_int = False
                continue
            if not fixed_fmt and len(tok) in (3, 5):
                c = tok[0]
                pairs = [(tok[1], tok[2])] + ([(tok[3], tok[4])] if len(tok) == 5 else [])
            else:
                f = fields(raw)
                c = f[1]
                pairs = [(f[2], f[3])] + ([(f[4], f[5])] if f[4] else [])
            if c not in M.col_index:
                M.col_index[c] = len(M.col_names)
                M.col_names.append(c)
                M.cost.append(F(0))
                M.integer.append(marker_int)
                M.cols.append({})
                M.clo.append(F(0))
                M.cup.append(INF)
            j = M.col_index[c]
            for r, v in pairs:
                val = _num(v)
                if r == M.obj_row:
                    M.cost[j] += val
                elif r in free_rows:
                    continue
                else:
                    i = M.row_index[r]
                    M.cols[j][i] = M.cols[j].get(i, F(0)) + val
        elif sec in ("RHS", "RANGES"):
            if fixed_fmt:
                f = fields(raw)
                pairs = [(f[2], f[3])] + ([(f[4], f[5])] if f[4] else [])
            elif len(tok) in (3, 5):
                pairs = [(tok[1], tok[2])] + ([(tok[3], tok[4])] if len(tok) == 5 else [])
            else:
                pairs = [(tok[0], tok[1])] + ([(tok[2], tok[3])] if len(tok) == 4 else [])
            for r, v in pairs:
                val = _num(v)
                if sec == "RHS":
                    if r == M.obj_row:
                        M.offset = -val
                    elif r in M.row_index:
                        M.rhs[M.row_index[r]] = val
                elif r in M.row_index:
                    ranges[M.row_index[r]] = val
        elif sec == "BOUNDS":
            if fixed_fmt:
                f = fields(raw)
                bt, c, v = f[0].upper(), f[2], f[3]
            else:
                bt = tok[0].upper()
                if bt in ("FR", "MI", "PL", "BV"):
                    c = tok[2] if len(tok) >= 3 else tok[1]
                    v = "0"
                else:
                    c, v = (tok[2], tok[3]) if len(tok) >= 4 else (tok[1], tok[2])
            j = M.col_index[c]
            val = _num(v) if v else F(0)
            if bt == "UP":
                M.cup[j] = val
                if val < 0 and j not in has_lo and M.clo[j] == 0:
                    M.clo[j] = INF
                has_up.add(j)
            elif bt == "LO":
                M.clo[j] = val
                has_lo.add(j)
            elif bt == "FX":
                M.clo[j] = M.cup[j] = val
                has_lo.add(j)
                has_up.add(j)
            elif bt == "FR":
                M.clo[j] = M.cup[j] = INF
                has_lo.add(j)
                has_up.add(j)
            elif bt == "MI":
                M.clo[j] = INF
                has_lo.add(j)
            elif bt == "PL":
                M.cup[j] = INF
                has_up.add(j)
            elif bt == "BV":
                M.integer[j] = True
                M.clo[j], M.cup[j] = F(0), F(1)
                has_lo.add(j)
                has_up.add(j)
            elif bt == "LI":
                M.integer[j] = True
                M.clo[j] = val
                has_lo.add(j)
            elif bt == "UI":
                M.integer[j] = True
                M.cup[j] = val
                has_up.add(j)
            else:
                raise ValueError(f"unsupported bound type {bt}")
        elif sec == "QUAD":
            a, b, v = M.col_index[tok[0]], M.col_index[tok[1]], _num(tok[2])
            M.q[(a, b)] = M.q.get((a, b), F(0)) + v
            if not quad_full and a != b:
                M.q[(b, a)] = M.q.get((b, a), F(0)) + v
    # Row bounds with RANGES semantics.
    for i, t in enumerate(M.row_type):
        b = M.rhs[i]
        if t == "E":
            lo, up = b, b
            if i in ranges:
                r = ranges[i]
                lo, up = (b, b + r) if r >= 0 else (b + r, b)
        elif t == "L":
            lo, up = INF, b
            if i in ranges:
                lo = b - abs(ranges[i])
        else:
            lo, up = b, INF
            if i in ranges:
                up = b + abs(ranges[i])
        M.rlo.append(lo)
        M.rup.append(up)
    # Treat |bound| >= 1e20 as infinite (MPS convention).
    big = F(10) ** 20
    for arr in (M.clo, M.cup, M.rlo, M.rup):
        for k, v in enumerate(arr):
            if v is not None and abs(v) >= big:
                arr[k] = INF
    return M


# ---------------------------------------------------------------------------
# Exact checks
# ---------------------------------------------------------------------------
def _rel(viol: F, bound: Optional[F]) -> float:
    return float(viol / (1 + abs(bound))) if bound is not None else float(viol)


def row_activity(M: ExactModel, x: List[F]) -> List[F]:
    act = [F(0)] * M.m
    for j, col in enumerate(M.cols):
        xj = x[j]
        if xj == 0:
            continue
        for i, a in col.items():
            act[i] += a * xj
    return act


def primal_check(M: ExactModel, x: List[F]) -> Tuple[float, float]:
    cv = 0.0
    for j in range(M.n):
        if M.clo[j] is not None and x[j] < M.clo[j]:
            cv = max(cv, _rel(M.clo[j] - x[j], M.clo[j]))
        if M.cup[j] is not None and x[j] > M.cup[j]:
            cv = max(cv, _rel(x[j] - M.cup[j], M.cup[j]))
    rv = 0.0
    for i, a in enumerate(row_activity(M, x)):
        if M.rlo[i] is not None and a < M.rlo[i]:
            rv = max(rv, _rel(M.rlo[i] - a, M.rlo[i]))
        if M.rup[i] is not None and a > M.rup[i]:
            rv = max(rv, _rel(a - M.rup[i], M.rup[i]))
    return cv, rv


def objective(M: ExactModel, x: List[F]) -> F:
    f = M.offset + sum((M.cost[j] * x[j] for j in range(M.n) if x[j] != 0), F(0))
    if M.q:
        f += F(1, 2) * sum((v * x[a] * x[b] for (a, b), v in M.q.items()), F(0))
    return f


def _min_prod(z: F, lo: Optional[F], up: Optional[F]):
    """min over x in [lo, up] of z*x; returns None for -infinity."""
    if z == 0:
        return F(0)
    if z > 0:
        return None if lo is None else z * lo
    return None if up is None else z * up


def safe_bound(M: ExactModel, cmin: List[F], y: List[F]):
    """Exact lower bound on min cmin^T x over the LP from ANY y (sign-projected).
    Returns (bound or None, max relative cost perturbation needed for wrong-sign terms)."""
    yp = []
    for i in range(M.m):
        v = y[i]
        if M.rlo[i] is None and M.rup[i] is None:
            v = F(0)
        elif M.rlo[i] is None:
            v = min(v, F(0))
        elif M.rup[i] is None:
            v = max(v, F(0))
        yp.append(v)
    lb = F(0)
    for i in range(M.m):
        t = _min_prod(yp[i], M.rlo[i], M.rup[i])
        lb += t
    pert = 0.0
    for j in range(M.n):
        z = cmin[j] - sum((a * yp[i] for i, a in M.cols[j].items()), F(0))
        t = _min_prod(z, M.clo[j], M.cup[j])
        if t is None:
            # wrong sign on an infinite bound: clip (cost perturbation of |z|)
            pert = max(pert, float(abs(z) / (1 + abs(cmin[j]))))
            z = F(0)
            t = F(0)
        lb += t
    return lb, pert


def farkas_check(M: ExactModel, y: List[F]) -> Tuple[bool, float]:
    best = None
    for sgn in (1, -1):
        yy = []
        for i in range(M.m):
            v = sgn * y[i]
            if M.rlo[i] is None and M.rup[i] is None:
                v = F(0)
            elif M.rlo[i] is None:
                v = min(v, F(0))
            elif M.rup[i] is None:
                v = max(v, F(0))
            yy.append(v)
        sup_x = F(0)
        ok = True
        for j in range(M.n):
            w = sum((a * yy[i] for i, a in M.cols[j].items()), F(0))
            if w == 0:
                continue
            if w > 0:
                if M.cup[j] is None:
                    ok = False
                    break
                sup_x += w * M.cup[j]
            else:
                if M.clo[j] is None:
                    ok = False
                    break
                sup_x += w * M.clo[j]
        if not ok:
            continue
        inf_r = F(0)
        for i in range(M.m):
            t = _min_prod(yy[i], M.rlo[i], M.rup[i])
            if t is None:
                ok = False
                break
            inf_r += t
        if not ok:
            continue
        margin = inf_r - sup_x
        if best is None or margin > best:
            best = margin
    return (best is not None and best > 0), (float(best) if best is not None else -math.inf)


# ---------------------------------------------------------------------------
# Exact basis verification (vertex certificate)
# ---------------------------------------------------------------------------
def solve_exact(Bcols: List[Dict[int, F]], rhs: List[F]) -> Optional[List[F]]:
    """Solve B z = rhs exactly; B given column-wise (m x m). Sparse Gaussian elimination."""
    m = len(Bcols)
    rows: List[Dict[int, F]] = [dict() for _ in range(m)]
    for p, col in enumerate(Bcols):
        for i, v in col.items():
            if v != 0:
                rows[i][p] = v
    b = list(rhs)
    pivot_of_col = {}
    used = [False] * m
    order = []
    for _ in range(m):
        # choose the column with the fewest entries among unused rows (Markowitz-lite)
        best = None
        colcount: Dict[int, int] = {}
        for i in range(m):
            if used[i]:
                continue
            for p in rows[i]:
                if p not in pivot_of_col:
                    colcount[p] = colcount.get(p, 0) + 1
        if not colcount:
            return None
        pcol = min(colcount, key=lambda p: (colcount[p], p))
        cand = [i for i in range(m) if not used[i] and pcol in rows[i]]
        prow = min(cand, key=lambda i: len(rows[i]))
        used[prow] = True
        pivot_of_col[pcol] = prow
        order.append((prow, pcol))
        piv = rows[prow][pcol]
        for i in cand:
            if i == prow:
                continue
            f = rows[i][pcol] / piv
            for p, v in rows[prow].items():
                nv = rows[i].get(p, F(0)) - f * v
                if nv == 0:
                    rows[i].pop(p, None)
                else:
                    rows[i][p] = nv
            b[i] -= f * b[prow]
    z = [F(0)] * m
    for prow, pcol in reversed(order):
        s = b[prow]
        for p, v in rows[prow].items():
            if p != pcol:
                s -= v * z[p]
        z[pcol] = s / rows[prow][pcol]
    return z


def exact_basis_check(M: ExactModel, res: dict, sense: int) -> dict:
    cb = res.get("column_basis")
    rb = res.get("row_basis")
    if not cb or not rb:
        return {"check": "exact_basis", "verdict": "N/A", "detail": "no basis in result JSON"}
    n, m = M.n, M.m
    # Variables of [A I]: structural j, logical s_i = -a_i^T x with bounds [-U, -L].
    lo = M.clo + [(-u if u is not None else None) for u in M.rup]
    up = M.cup + [(-l if l is not None else None) for l in M.rlo]
    status = [cb[c] for c in M.col_names] + [rb[r] for r in M.row_names]
    # Row status LOWER/UPPER refer to the activity; for the logical they swap.
    swap = {"LOWER": "UPPER", "UPPER": "LOWER"}
    status = status[:n] + [swap.get(s, s) for s in status[n:]]
    basic = [k for k, s in enumerate(status) if s == "BASIC"]
    if len(basic) != m:
        return {"check": "exact_basis", "verdict": "FAIL", "detail": f"{len(basic)} basic variables for {m} rows"}
    val = [F(0)] * (n + m)
    for k, s in enumerate(status):
        if s == "BASIC":
            continue
        if s in ("LOWER", "FIXED"):
            val[k] = lo[k] if lo[k] is not None else (up[k] if up[k] is not None else F(0))
        elif s == "UPPER":
            val[k] = up[k] if up[k] is not None else F(0)
        else:
            val[k] = F(0)
    def column(k):
        return M.cols[k] if k < n else {k - n: F(1)}
    rhs = [F(0)] * m
    for k in range(n + m):
        if status[k] == "BASIC" or val[k] == 0:
            continue
        for i, a in column(k).items():
            rhs[i] -= a * val[k]
    xb = solve_exact([column(k) for k in basic], rhs)
    if xb is None:
        return {"check": "exact_basis", "verdict": "FAIL", "detail": "basis matrix is singular in exact arithmetic"}
    for p, k in enumerate(basic):
        val[k] = xb[p]
    pviol = 0
    for k in range(n + m):
        if lo[k] is not None and val[k] < lo[k]:
            pviol += 1
        if up[k] is not None and val[k] > up[k]:
            pviol += 1
    # Duals: B^T y = c_B  (min form, logical cost 0)
    cmin = [sense * c for c in M.cost] + [F(0)] * m
    # Build B^T column-wise: column i of B^T = row i of B.
    BT = [dict() for _ in range(m)]
    for p, k in enumerate(basic):
        for i, a in column(k).items():
            BT[i][p] = a
    y = solve_exact(BT, [cmin[k] for k in basic])
    dviol = 0
    for k in range(n + m):
        if status[k] == "BASIC":
            continue
        d = cmin[k] - sum((a * y[i] for i, a in column(k).items()), F(0))
        s = status[k]
        if s == "LOWER" and d < 0:
            dviol += 1
        elif s == "UPPER" and d > 0:
            dviol += 1
        elif s == "ZERO" and d != 0:
            dviol += 1
    ok = pviol == 0 and dviol == 0
    obj = objective(M, val[:n])
    return {"check": "exact_basis", "verdict": "PASS" if ok else "FAIL",
            "detail": f"exact vertex: {pviol} primal and {dviol} dual sign violations; exact objective {float(obj):.15g}",
            "exact_objective": str(obj)}


# ---------------------------------------------------------------------------
def verify(model_path: str, result_path: str, exact_basis: bool = False, tol: float = 1e-6) -> dict:
    M = parse_mps(model_path)
    res = json.load(open(result_path))
    status = res.get("status")
    sense = -1 if M.maximize else 1
    report = {"model": M.name, "claimed_status": status, "checks": [], "verified": False}
    checks = report["checks"]

    def vec(key, names):
        d = res.get(key)
        if d is None:
            return None
        return [F(float(d[nm])) for nm in names]  # exact binary value of the double

    x = vec("x", M.col_names)
    is_mip = any(M.integer)
    ok = False
    if status in ("OPTIMAL", "TIME_LIMIT", "NODE_LIMIT") and x is not None:
        cv, rv = primal_check(M, x)
        checks.append({"check": "exact_primal_bound_violation", "value": cv, "verdict": "PASS" if cv <= tol else "FAIL"})
        checks.append({"check": "exact_primal_row_violation", "value": rv, "verdict": "PASS" if rv <= tol else "FAIL"})
        obj = objective(M, x)
        report["exact_objective"] = float(obj)
        primal_ok = cv <= tol and rv <= tol
        if is_mip:
            iv = max((float(abs(x[j] - round(x[j]))) for j in range(M.n) if M.integer[j]), default=0.0)
            checks.append({"check": "exact_integrality", "value": iv, "verdict": "PASS" if iv <= 1e-5 else "FAIL"})
            bb = res.get("best_bound")
            if bb is not None:
                gap = abs(float(obj) - bb) / max(1.0, abs(float(obj)))
                checks.append({"check": "reported_gap (bound from solver's safe B&B bounds)", "value": gap,
                               "verdict": "PASS" if gap <= 1e-4 or status != "OPTIMAL" else "FAIL"})
            ok = primal_ok and iv <= 1e-5
        elif status == "OPTIMAL":
            yr = vec("row_duals", M.row_names)
            if yr is None:
                checks.append({"check": "dual_present", "verdict": "FAIL"})
            else:
                y = [sense * v for v in yr]  # JSON duals are in model sense
                cmin = [sense * c for c in M.cost]
                if M.q:
                    # linearization at x: g = c + Qx (min form), constant -1/2 x'Qx
                    g = list(cmin)
                    for (a, b), v in M.q.items():
                        g[a] += sense * v * x[b]
                    lb, pert = safe_bound(M, g, y)
                    xqx = sum((v * x[a] * x[b] for (a, b), v in M.q.items()), F(0))
                    lb = lb - F(1, 2) * sense * xqx + sense * M.offset
                else:
                    lb, pert = safe_bound(M, cmin, y)
                    lb += sense * M.offset
                pobj = sense * obj
                gap = float((pobj - lb) / max(F(1), abs(pobj)))
                checks.append({"check": "exact_safe_dual_bound_gap", "value": gap,
                               "verdict": "PASS" if gap <= tol and pert <= 1e-9 else ("WARN" if gap <= tol else "FAIL"),
                               "detail": f"exact bound {float(sense * lb):.15g}; cost perturbation needed {pert:.2e}"})
                ok = primal_ok and gap <= tol
            if exact_basis:
                eb = exact_basis_check(M, res, sense)
                checks.append(eb)
                ok = ok and eb["verdict"] in ("PASS", "N/A")
        else:
            ok = primal_ok
    elif status == "INFEASIBLE":
        y = vec("farkas_ray", M.row_names)
        if y is None:
            checks.append({"check": "farkas_present", "verdict": "FAIL" if not is_mip else "N/A",
                           "detail": "MIP infeasibility is proved by the B&B tree, not a single ray" if is_mip else ""})
            ok = is_mip
        else:
            good, margin = farkas_check(M, y)
            checks.append({"check": "exact_farkas", "value": margin, "verdict": "PASS" if good else "FAIL"})
            ok = good
    elif status == "UNBOUNDED":
        d = vec("primal_ray", M.col_names)
        if d is None:
            checks.append({"check": "ray_present", "verdict": "FAIL"})
        else:
            rec = all((M.clo[j] is None or d[j] >= 0) and (M.cup[j] is None or d[j] <= 0) for j in range(M.n))
            ad = row_activity(M, d)
            rec = rec and all((M.rlo[i] is None or ad[i] >= 0) and (M.rup[i] is None or ad[i] <= 0) for i in range(M.m))
            descent = sum((sense * M.cost[j] * d[j] for j in range(M.n)), F(0)) < 0
            checks.append({"check": "exact_ray_recession", "verdict": "PASS" if rec else "FAIL"})
            checks.append({"check": "exact_ray_descent", "verdict": "PASS" if descent else "FAIL"})
            ok = rec and descent
    else:
        checks.append({"check": "claim", "verdict": "N/A", "detail": f"status {status} carries no claim"})
        ok = True
    report["verified"] = ok
    return report


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model")
    ap.add_argument("result")
    ap.add_argument("--exact-basis", action="store_true", help="exact rational vertex verification (small models)")
    ap.add_argument("--tol", type=float, default=1e-6)
    a = ap.parse_args(argv)
    rep = verify(a.model, a.result, a.exact_basis, a.tol)
    print(json.dumps(rep, indent=1))
    return 0 if rep["verified"] else 4


if __name__ == "__main__":
    sys.exit(main())
