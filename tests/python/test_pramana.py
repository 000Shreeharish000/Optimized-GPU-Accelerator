"""End-to-end tests: Python bindings, CLI, independent exact verifier, generators."""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "gen"))
EXE = ROOT / "build" / ("pramana.exe" if os.name == "nt" else "pramana")

import pramana  # noqa: E402
from pramana import verify as V  # noqa: E402


def cli(*args, timeout=600):
    return subprocess.run([str(EXE), *map(str, args)], capture_output=True, text=True, timeout=timeout)


def test_bindings_lp_max_with_duals():
    m = pramana.Model("demo", maximize=True)
    x1 = m.add_var(0, pramana.INF, 3, name="x1")
    x2 = m.add_var(0, pramana.INF, 5, name="x2")
    x3 = m.add_var(0, 1.5, 4, name="x3")
    m.add_constraint({x1: 1, x2: 2, x3: 1}, ub=8, name="cap1")
    m.add_constraint({x1: 2, x2: 1, x3: 3}, lb=5, ub=9, name="cap2")
    m.add_constraint({x1: 1, x2: -1}, lb=0, ub=0, name="bal")
    r = m.solve(algorithm="dual")
    assert r.status == "OPTIMAL" and r.certified
    assert abs(r.objective - 22.0) < 1e-9
    assert len(r.row_duals) == 3 and r.certificate["accepted"]


def test_bindings_mip_and_qp():
    m = pramana.Model("knap", maximize=True)
    w = [3, 4, 5, 8, 9]
    p = [4, 5, 6, 10, 11]
    xs = [m.add_var(0, 1, p[i], integer=True) for i in range(5)]
    m.add_constraint({x: w[i] for i, x in enumerate(xs)}, ub=14)
    r = m.solve()
    assert r.status == "OPTIMAL" and r.certified
    best = max(sum(p[i] for i in range(5) if s >> i & 1) for s in range(32) if sum(w[i] for i in range(5) if s >> i & 1) <= 14)
    assert abs(r.objective - best) < 1e-9
    q = pramana.Model("qp")
    x = q.add_var(-pramana.INF, pramana.INF, -4)
    y = q.add_var(-pramana.INF, pramana.INF, -6)
    q.add_quadratic(x, x, 2)
    q.add_quadratic(y, y, 2)
    q.add_constraint({x: 1, y: 1}, ub=2)
    rq = q.solve()
    assert rq.status == "OPTIMAL" and abs(rq.objective + 8.5) < 1e-7


def test_nonconvex_qp_refused_with_evidence():
    q = pramana.Model("nc")
    x = q.add_var(0, 1)
    y = q.add_var(0, 1)
    q.add_quadratic(x, y, 2.0)  # Q = [[0,2],[2,0]] indefinite
    q.add_quadratic(x, x, 1.0)
    q.add_quadratic(y, y, 1.0)
    r = q.solve()
    assert r.status == "NONCONVEX"


def test_cli_json_and_exact_verifier(tmp_path):
    out = tmp_path / "afiro.json"
    p = cli(ROOT / "data/netlib/afiro.mps.gz", "--json", out, "--vectors", "--log", "0")
    assert p.returncode == 0, p.stdout + p.stderr
    rep = V.verify(str(ROOT / "data/netlib/afiro.mps.gz"), str(out), exact_basis=True)
    assert rep["verified"], rep
    assert any(c["check"] == "exact_basis" and c["verdict"] == "PASS" for c in rep["checks"])
    # Tamper with the claimed solution: the exact verifier must reject it.
    j = json.loads(out.read_text())
    k = next(iter(j["x"]))
    j["x"][k] += 1.0
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(j))
    assert not V.verify(str(ROOT / "data/netlib/afiro.mps.gz"), str(bad))["verified"]


def test_exact_farkas_on_infeasible():
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        out = Path(d) / "r.json"
        cli(ROOT / "data/netlib_infeas/klein1.mps", "--json", out, "--vectors", "--log", "0")
        rep = V.verify(str(ROOT / "data/netlib_infeas/klein1.mps"), str(out))
        assert rep["claimed_status"] == "INFEASIBLE" and rep["verified"], rep


def test_generators_valid_models(tmp_path):
    import refinery
    refinery.planning(4, 3, False, 1, str(tmp_path / "p.mps"))
    refinery.unloading(2, 2, 2, 6, "tight", 1, str(tmp_path / "u.mps"))
    refinery.dispatch(5, 4, 1, str(tmp_path / "d.qps"))
    for f in ("p.mps", "u.mps", "d.qps"):
        r = pramana.solve_file(str(tmp_path / f), time_limit=60)
        assert r.status == "OPTIMAL" and r.certified, (f, r)


def test_parametric_binding():
    res = pramana.parametric(str(ROOT / "data/netlib/afiro.mps.gz"), col="X02", kind="cost", lo=-2, hi=1)
    assert res["complete"] and len(res["segments"]) == 2
    assert all(s["certified"] for s in res["segments"])
    assert abs(res["breakpoints"][0] - (-0.0552286)) < 1e-6


def test_verify_command_roundtrip(tmp_path):
    out = tmp_path / "r.json"
    cli(ROOT / "tests/data/tiny_max.mps", "--json", out, "--vectors", "--log", "0")
    p = cli("verify", ROOT / "tests/data/tiny_max.mps", out)
    assert p.returncode == 0 and "ACCEPTED" in p.stdout


def test_terminal_console_scripted():
    """The terminal console (python -m pramana.tui) solves, analyses, debugs and verifies from scripted input."""
    env = dict(os.environ, PYTHONPATH=str(ROOT / "python"), NO_COLOR="1", COLUMNS="110")
    script = "afiro\n/debug\n/verify\n/analyze p0033\n/exit\n"
    p = subprocess.run([sys.executable, "-m", "pramana.tui"], input=script, capture_output=True, text=True,
                       encoding="utf-8", env=env, cwd=str(ROOT), timeout=300)
    assert p.returncode == 0, p.stderr
    o = p.stdout
    assert "PRAMANA v" in o and "OPTIMAL" in o and "certified" in o
    assert "Debug" in o and "Simplex" in o
    assert "✓ accepted" in o and "✓ verified" in o
    assert "Analysis" in o and "MILP" in o
