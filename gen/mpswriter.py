"""Tiny algebraic model builder that writes free MPS / QPS (used by all generators).

All data produced by the generators is SYNTHETIC, generated from published model
structures and public ranges (assays, prices, yields). It is not MRPL data.
"""
from __future__ import annotations

import math
from collections import OrderedDict
from typing import Dict, Iterable, Optional

INF = math.inf


class LinearModel:
    def __init__(self, name: str, maximize: bool = False):
        self.name = name
        self.maximize = maximize
        self.cols: "OrderedDict[str, dict]" = OrderedDict()   # name -> {cost, lb, ub, int}
        self.rows: "OrderedDict[str, dict]" = OrderedDict()   # name -> {lb, ub, coefs{col: v}}
        self.quad: Dict[tuple, float] = {}                     # (i, j) i<=j -> v  (1/2 x'Qx)
        self.offset = 0.0

    def var(self, name: str, lb: float = 0.0, ub: float = INF, cost: float = 0.0, integer: bool = False) -> str:
        if name in self.cols:
            raise ValueError(f"duplicate column {name}")
        self.cols[name] = {"cost": cost, "lb": lb, "ub": ub, "int": integer}
        return name

    def cost(self, name: str, c: float) -> None:
        self.cols[name]["cost"] += c

    def row(self, name: str, coefs: Dict[str, float], lb: float = -INF, ub: float = INF) -> str:
        if name in self.rows:
            raise ValueError(f"duplicate row {name}")
        clean = {k: v for k, v in coefs.items() if v != 0}
        for k in clean:
            if k not in self.cols:
                raise KeyError(f"row {name}: unknown column {k}")
        self.rows[name] = {"lb": lb, "ub": ub, "coefs": clean}
        return name

    def qterm(self, a: str, b: str, v: float) -> None:
        key = (a, b) if a <= b else (b, a)
        self.quad[key] = self.quad.get(key, 0.0) + v

    @property
    def num_cols(self):
        return len(self.cols)

    @property
    def num_rows(self):
        return len(self.rows)

    def write(self, path: str) -> None:
        f = lambda v: repr(float(v))  # noqa: E731
        out = [f"NAME {self.name}"]
        if self.maximize:
            out += ["OBJSENSE", "    MAX"]
        out += ["ROWS", " N  OBJ"]
        types = {}
        for r, d in self.rows.items():
            lo, up = d["lb"], d["ub"]
            if lo == up:
                t = "E"
            elif math.isfinite(lo) and math.isfinite(up):
                t = "L"
            elif math.isfinite(lo):
                t = "G"
            elif math.isfinite(up):
                t = "L"
            else:
                t = "N"
            types[r] = t
            out.append(f" {t}  {r}")
        colrows: Dict[str, list] = {c: [] for c in self.cols}
        for r, d in self.rows.items():
            for c, v in d["coefs"].items():
                colrows[c].append((r, v))
        out.append("COLUMNS")
        in_int = False
        for c, d in self.cols.items():
            if d["int"] != in_int:
                tag = "'INTORG'" if d["int"] else "'INTEND'"
                out.append(f"    M{len(out)} 'MARKER' {tag}")
                in_int = d["int"]
            any_entry = False
            if d["cost"] != 0:
                out.append(f"    {c} OBJ {f(d['cost'])}")
                any_entry = True
            for r, v in colrows[c]:
                out.append(f"    {c} {r} {f(v)}")
                any_entry = True
            if not any_entry:
                out.append(f"    {c} OBJ 0")
        if in_int:
            out.append(f"    M{len(out)} 'MARKER' 'INTEND'")
        out.append("RHS")
        if self.offset:
            out.append(f"    RHS OBJ {f(-self.offset)}")
        rng = []
        for r, d in self.rows.items():
            t, lo, up = types[r], d["lb"], d["ub"]
            if t == "E" and lo != 0:
                out.append(f"    RHS {r} {f(lo)}")
            elif t == "L" and up != 0:
                out.append(f"    RHS {r} {f(up)}")
            elif t == "G" and lo != 0:
                out.append(f"    RHS {r} {f(lo)}")
            if t == "L" and math.isfinite(lo) and math.isfinite(up) and lo != up:
                rng.append(f"    RNG {r} {f(up - lo)}")
        if rng:
            out.append("RANGES")
            out += rng
        out.append("BOUNDS")
        for c, d in self.cols.items():
            lo, up = d["lb"], d["ub"]
            if lo == up:
                out.append(f" FX BND {c} {f(lo)}")
                continue
            if not math.isfinite(lo) and not math.isfinite(up):
                out.append(f" FR BND {c}")
                continue
            if not math.isfinite(lo):
                out.append(f" MI BND {c}")
            elif lo != 0:
                out.append(f" LO BND {c} {f(lo)}")
            if math.isfinite(up):
                out.append(f" UP BND {c} {f(up)}")
            elif d["int"]:
                out.append(f" PL BND {c}")
        if self.quad:
            out.append("QUADOBJ")
            for (a, b), v in self.quad.items():
                out.append(f"    {a} {b} {f(v)}")
        out.append("ENDATA")
        with open(path, "w") as fh:
            fh.write("\n".join(out) + "\n")
