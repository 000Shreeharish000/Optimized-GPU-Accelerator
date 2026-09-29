"""PRAMANA Python bindings (ctypes over the C API; no compiled extension needed).

Example
-------
>>> import pramana
>>> m = pramana.Model("demo", maximize=True)
>>> x = m.add_var(lb=0, ub=4, cost=3, name="x")
>>> y = m.add_var(lb=0, cost=5, name="y", integer=True)
>>> m.add_constraint({x: 1, y: 2}, ub=8, name="cap")
>>> r = m.solve(algorithm="auto")
>>> r.status, r.objective, r.certified
('OPTIMAL', 22.0, True)

Every result carries the certificate produced by the independent certifier;
`pramana.verify` (python/pramana/verify.py) re-checks it in exact rational
arithmetic with its own MPS parser.
"""
from __future__ import annotations

import ctypes
import json
import os
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional

__all__ = ["Model", "Result", "solve_file", "parametric", "family", "gpu_info", "load_library"]

_LIB = None


def _candidate_paths():
    env = os.environ.get("PRAMANA_LIB")
    if env:
        yield Path(env)
    here = Path(__file__).resolve().parent
    root = here.parent.parent
    names = ["libpramana.dll", "libpramana.so", "liblibpramana.so", "libpramana.dylib"]
    for d in [here, root / "build", root / "build" / "Release"]:
        for n in names:
            yield d / n


def load_library():
    global _LIB
    if _LIB is not None:
        return _LIB
    for p in _candidate_paths():
        if p.exists():
            _LIB = ctypes.CDLL(str(p))
            break
    if _LIB is None:
        raise OSError("PRAMANA shared library not found; build it (build.bat) or set PRAMANA_LIB")
    L = _LIB
    c_int, c_double, c_char_p, c_void_p = ctypes.c_int, ctypes.c_double, ctypes.c_char_p, ctypes.c_void_p
    dp = ctypes.POINTER(c_double)
    ip = ctypes.POINTER(c_int)
    sig = {
        "pramana_version": ([], c_char_p),
        "pramana_last_error": ([], c_char_p),
        "pramana_model_create": ([c_char_p, c_int], c_void_p),
        "pramana_model_read": ([c_char_p], c_void_p),
        "pramana_model_free": ([c_void_p], None),
        "pramana_add_col": ([c_void_p, c_double, c_double, c_double, c_int, c_char_p], c_int),
        "pramana_add_row": ([c_void_p, c_double, c_double, c_int, ip, dp, c_char_p], c_int),
        "pramana_add_q": ([c_void_p, c_int, c_int, c_double], c_int),
        "pramana_set_offset": ([c_void_p, c_double], None),
        "pramana_num_cols": ([c_void_p], c_int),
        "pramana_num_rows": ([c_void_p], c_int),
        "pramana_write_mps": ([c_void_p, c_char_p], c_int),
        "pramana_solve": ([c_void_p, c_char_p], c_void_p),
        "pramana_result_status_name": ([c_void_p], c_char_p),
        "pramana_result_objective": ([c_void_p], c_double),
        "pramana_result_bound": ([c_void_p], c_double),
        "pramana_result_certified": ([c_void_p], c_int),
        "pramana_result_x": ([c_void_p, dp, c_int], c_int),
        "pramana_result_row_duals": ([c_void_p, dp, c_int], c_int),
        "pramana_result_reduced_costs": ([c_void_p, dp, c_int], c_int),
        "pramana_result_row_activity": ([c_void_p, dp, c_int], c_int),
        "pramana_result_json": ([c_void_p, c_int], c_char_p),
        "pramana_result_free": ([c_void_p], None),
        "pramana_parametric": ([c_void_p, c_char_p], c_void_p),
        "pramana_family": ([c_void_p, c_char_p], c_void_p),
        "pramana_gpu_info": ([], c_void_p),
        "pramana_free_string": ([c_void_p], None),
    }
    for name, (args, res) in sig.items():
        f = getattr(L, name)
        f.argtypes = args
        f.restype = res
    return L


def _take_string(ptr) -> str:
    L = load_library()
    if not ptr:
        raise RuntimeError(L.pramana_last_error().decode())
    s = ctypes.cast(ptr, ctypes.c_char_p).value.decode()
    L.pramana_free_string(ptr)
    return s


INF = float("inf")


class Var:
    __slots__ = ("model", "index", "name")

    def __init__(self, model: "Model", index: int, name: str):
        self.model, self.index, self.name = model, index, name

    def __hash__(self):
        return hash((id(self.model), self.index))

    def __eq__(self, other):
        return isinstance(other, Var) and other.model is self.model and other.index == self.index

    def __repr__(self):
        return f"Var({self.name})"


class Result:
    """Solve result; all values in the model's objective sense."""

    def __init__(self, model: "Model", handle):
        L = load_library()
        self._model = model
        self.status: str = L.pramana_result_status_name(handle).decode()
        self.objective: float = L.pramana_result_objective(handle)
        self.bound: float = L.pramana_result_bound(handle)
        self.certified: bool = bool(L.pramana_result_certified(handle))

        def vec(fn):
            n = fn(handle, None, 0)
            buf = (ctypes.c_double * max(n, 1))()
            fn(handle, buf, n)
            return list(buf[:n])

        self.x = vec(L.pramana_result_x)
        self.row_duals = vec(L.pramana_result_row_duals)
        self.reduced_costs = vec(L.pramana_result_reduced_costs)
        self.row_activity = vec(L.pramana_result_row_activity)
        self.json = json.loads(L.pramana_result_json(handle, 1).decode())
        L.pramana_result_free(handle)

    @property
    def certificate(self) -> dict:
        return self.json.get("certificate", {})

    @property
    def telemetry(self) -> dict:
        return self.json.get("telemetry", {})

    def value(self, v: Var) -> float:
        return self.x[v.index]

    def __repr__(self):
        return f"Result(status={self.status}, objective={self.objective}, certified={self.certified})"


class Model:
    """In-memory model builder:  min/max c^T x + 1/2 x^T Q x,  L <= Ax <= U,  l <= x <= u."""

    def __init__(self, name: str = "model", maximize: bool = False, _handle=None):
        L = load_library()
        self._h = _handle or L.pramana_model_create(name.encode(), 1 if maximize else 0)
        if not self._h:
            raise RuntimeError(L.pramana_last_error().decode())
        self._vars = []

    @classmethod
    def read(cls, path: str) -> "Model":
        L = load_library()
        h = L.pramana_model_read(str(path).encode())
        if not h:
            raise RuntimeError(L.pramana_last_error().decode())
        return cls(_handle=h)

    def __del__(self):
        try:
            if getattr(self, "_h", None):
                load_library().pramana_model_free(self._h)
        except Exception:  # noqa: BLE001
            pass

    def add_var(self, lb: float = 0.0, ub: float = INF, cost: float = 0.0, integer: bool = False,
                name: Optional[str] = None) -> Var:
        L = load_library()
        nm = name or f"x{len(self._vars)}"
        j = L.pramana_add_col(self._h, cost, max(lb, -1e30), min(ub, 1e30), 1 if integer else 0, nm.encode())
        v = Var(self, j, nm)
        self._vars.append(v)
        return v

    def add_constraint(self, coefs: Mapping[Var, float], lb: float = -INF, ub: float = INF,
                       name: Optional[str] = None) -> int:
        L = load_library()
        items = [(v.index, float(a)) for v, a in coefs.items() if a != 0]
        n = len(items)
        cols = (ctypes.c_int * max(n, 1))(*[i for i, _ in items])
        vals = (ctypes.c_double * max(n, 1))(*[a for _, a in items])
        r = L.pramana_add_row(self._h, max(lb, -1e30), min(ub, 1e30), n, cols, vals, (name or "").encode())
        if r < 0:
            raise RuntimeError(L.pramana_last_error().decode())
        return r

    def add_quadratic(self, v1: Var, v2: Var, coef: float) -> None:
        """Adds coef to Q[v1,v2] (and Q[v2,v1]); objective term is 1/2 x^T Q x."""
        load_library().pramana_add_q(self._h, v1.index, v2.index, coef)

    def set_offset(self, offset: float) -> None:
        load_library().pramana_set_offset(self._h, offset)

    @property
    def num_vars(self) -> int:
        return load_library().pramana_num_cols(self._h)

    @property
    def num_constraints(self) -> int:
        return load_library().pramana_num_rows(self._h)

    def write_mps(self, path: str) -> None:
        if load_library().pramana_write_mps(self._h, str(path).encode()) != 0:
            raise RuntimeError(load_library().pramana_last_error().decode())

    def solve(self, **options) -> Result:
        """Options: algorithm (auto|dual|ipm|pdhg|pdhg-cpu|pdhg-gpu|race), time_limit, presolve,
        certify, threads, log_level, tolerance, pdhg_tolerance, ipm_tolerance, mip_gap,
        node_limit, cuts, heuristics, branching, allow_gpu, router_model."""
        L = load_library()
        h = L.pramana_solve(self._h, json.dumps(options).encode())
        if not h:
            raise RuntimeError(L.pramana_last_error().decode())
        return Result(self, h)

    def parametric(self, *, col: Optional[str] = None, row: Optional[str] = None, kind: str = "cost",
                   lo: float = 0.0, hi: float = 1.0) -> dict:
        spec = {"kind": kind, "from": lo, "to": hi}
        spec["col" if col is not None else "row"] = col if col is not None else row
        return json.loads(_take_string(load_library().pramana_parametric(self._h, json.dumps(spec).encode())))

    def family(self, *, col: Optional[str] = None, row: Optional[str] = None, kind: str = "cost",
               lo: float = 0.0, hi: float = 1.0, cases: int = 32, allow_gpu: bool = True) -> dict:
        spec = {"kind": kind, "from": lo, "to": hi, "cases": cases, "allow_gpu": allow_gpu}
        spec["col" if col is not None else "row"] = col if col is not None else row
        return json.loads(_take_string(load_library().pramana_family(self._h, json.dumps(spec).encode())))


def solve_file(path: str, **options) -> Result:
    return Model.read(path).solve(**options)


def parametric(path: str, **kw) -> dict:
    return Model.read(path).parametric(**kw)


def family(path: str, **kw) -> dict:
    return Model.read(path).family(**kw)


def gpu_info() -> dict:
    return json.loads(_take_string(load_library().pramana_gpu_info()))


def version() -> str:
    return load_library().pramana_version().decode()
