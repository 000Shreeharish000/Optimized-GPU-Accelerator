"""Download public benchmark instances used by the PRAMANA benchmark harness.

Sets (all public, none are MRPL data):
  netlib        feasible Netlib LPs (coin-or-tools/Data-Netlib mirror, .mps.gz)
  netlib-infeas Chinneck infeasible Netlib LPs (netlib.org, 'emps' compressed)
  miplib3       MIPLIB 3 instances (coin-or-tools/Data-miplib3 mirror)
  miplib2017    MIPLIB 2017 instances listed in bench/sets/miplib2017_subset.txt
  maros         Maros-Meszaros convex QPs (qpsolvers/maros_meszaros_qpbenchmark, .mat -> .qps)

Only the Python standard library is needed, except the Maros-Meszaros
conversion which uses scipy.io.loadmat to read MATLAB files (data conversion
only; no solver code is involved).

Usage: python bench/fetch_data.py [set ...]   (default: netlib netlib-infeas miplib3 maros)
"""
from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

NETLIB_RAW = "https://raw.githubusercontent.com/coin-or-tools/Data-Netlib/master/"
MIPLIB3_RAW = "https://raw.githubusercontent.com/coin-or-tools/Data-miplib3/master/"
MAROS_RAW = "https://raw.githubusercontent.com/qpsolvers/maros_meszaros_qpbenchmark/main/data/"
NETLIB_ORG = "https://www.netlib.org/lp/"
MIPLIB2017 = "https://miplib.zib.de/WebData/instances/"

NETLIB_FEASIBLE = """25fv47 80bau3b adlittle afiro agg agg2 agg3 bandm beaconfd blend bnl1 bnl2 boeing1
boeing2 bore3d brandy capri cycle czprob d2q06c d6cube degen2 degen3 dfl001 e226 etamacro fffff800
finnis fit1d fit1p fit2d fit2p forplan ganges gfrd-pnc greenbea greenbeb grow15 grow22 grow7 israel
kb2 lotfi maros-r7 maros modszk1 nesm perold pilot pilot4 pilot87 pilotnov recipe sc105 sc205 sc50a
sc50b scagr25 scagr7 scfxm1 scfxm2 scfxm3 scorpion scrs8 scsd1 scsd6 scsd8 sctap1 sctap2 sctap3 seba
share1b share2b shell ship04l ship04s ship08l ship08s ship12l ship12s sierra stair standata standgub
standmps stocfor1 stocfor2 tuff vtpbase wood1p woodw""".split()

NETLIB_INFEAS = """bgdbg1 bgetam bgindy bgprtr box1 ceria3d chemcom cplex1 cplex2 ex72a ex73a forest6
galenet gosh gran greenbea itest2 itest6 klein1 klein2 klein3 mondou2 pang pilot4i qual reactor
refinery vol1 woodinfe""".split()

MIPLIB3 = """p0033 p0201 p0282 p0548 p2756 lseu mod008 mod010 stein27 stein45 egout enigma flugpl
gt2 bell3a bell5 rgn vpm1 vpm2 dcmulti fixnet6 khb05250 misc03 misc06 misc07 pk1 pp08a pp08aCUTS
qnet1 qnet1_o blend2 gesa2 gesa2_o gesa3 gesa3_o set1ch mas74 mas76 noswot 10teams air03 air04
air05 cap6000 fiber gen harp2 l152lav mitre mod011 modglob qiu rentacar markshare1 markshare2
arki001 danoint dano3mip dsbmip fast0507 mkc nw04 rout seymour swath""".split()

MAROS_SMALLISH = """AUG2D AUG2DC AUG2DCQP AUG2DQP AUG3D AUG3DC AUG3DCQP AUG3DQP CONT-050 CVXQP1_S CVXQP1_M
CVXQP2_S CVXQP2_M CVXQP3_S CVXQP3_M DPKLO1 DTOC3 DUAL1 DUAL2 DUAL3 DUAL4 DUALC1 DUALC2 DUALC5 DUALC8
GENHS28 GOULDQP2 GOULDQP3 HS118 HS21 HS268 HS35 HS35MOD HS51 HS52 HS53 HS76 HUES-MOD HUESTIS KSIP
LASER LISWET1 LISWET2 LISWET3 LISWET4 LISWET5 LISWET6 LISWET7 LISWET8 LISWET9 LISWET10 LISWET11
LISWET12 LOTSCHD MOSARQP1 MOSARQP2 POWELL20 PRIMAL1 PRIMAL2 PRIMAL3 PRIMAL4 PRIMALC1 PRIMALC2
PRIMALC5 PRIMALC8 Q25FV47 QADLITTL QAFIRO QBANDM QBEACONF QBORE3D QBRANDY QCAPRI QE226 QETAMACR
QFFFFF80 QFORPLAN QGFRDXPN QGROW15 QGROW22 QGROW7 QISRAEL QPCBLEND QPCBOEI1 QPCBOEI2 QPCSTAIR
QPILOTNO QPTEST QRECIPE QSC205 QSCAGR25 QSCAGR7 QSCFXM1 QSCFXM2 QSCFXM3 QSCORPIO QSCRS8 QSCSD1
QSCSD6 QSCSD8 QSCTAP1 QSCTAP2 QSCTAP3 QSEBA QSHARE1B QSHARE2B QSHELL QSHIP04L QSHIP04S QSHIP08L
QSHIP08S QSHIP12L QSHIP12S QSIERRA QSTAIR QSTANDAT S268 STADAT1 STADAT2 STADAT3 TAME UBH1 YAO
ZECEVIC2""".split()


def fetch(url: str, dest: Path, retries: int = 3) -> bool:
    if dest.exists() and dest.stat().st_size > 0:
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r, open(str(dest) + ".part", "wb") as f:
                shutil.copyfileobj(r, f)
            os.replace(str(dest) + ".part", dest)
            return True
        except Exception as e:  # noqa: BLE001
            if attempt == retries - 1:
                print(f"  FAILED {url}: {e}", file=sys.stderr)
    return False


def gunzip(src: Path, dest: Path) -> None:
    if dest.exists():
        return
    with gzip.open(src, "rb") as fi, open(dest, "wb") as fo:
        shutil.copyfileobj(fi, fo)


def parallel(jobs, workers: int = 8):
    with ThreadPoolExecutor(workers) as ex:
        return list(ex.map(lambda j: j(), jobs))


def get_netlib() -> None:
    out = DATA / "netlib"
    print(f"netlib -> {out}")
    jobs = [lambda n=n: fetch(NETLIB_RAW + n + ".mps.gz", out / (n + ".mps.gz")) for n in NETLIB_FEASIBLE]
    ok = parallel(jobs)
    print(f"  {sum(ok)}/{len(ok)} files")
    fetch(NETLIB_ORG + "data/readme", out / "README.netlib")


def build_emps() -> Path | None:
    tools = DATA / "tools"
    exe = tools / ("emps.exe" if os.name == "nt" else "emps")
    if exe.exists():
        return exe
    if not fetch(NETLIB_ORG + "data/emps.c", tools / "emps.c"):
        return None
    # emps.c is Netlib's *data decompressor* (not a solver). Built only for data preparation.
    if os.name == "nt":
        vcvars = os.environ.get(
            "PRAMANA_VCVARS",
            r"C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat")
        bat = tools / "build_emps.bat"
        bat.write_text("@echo off\r\n"
                       f'call "{vcvars}" >nul 2>nul\r\n'
                       'cd /d "%~dp0"\r\n'
                       "cl /nologo /O2 /w /D_CRT_SECURE_NO_WARNINGS emps.c /Fe:emps.exe\r\n")
        subprocess.run(["cmd", "/c", str(bat)], check=False, capture_output=True)
    else:
        subprocess.run(["cc", "-O2", "-w", "-o", str(exe), str(tools / "emps.c")], check=False)
    return exe if exe.exists() else None


def get_netlib_infeas() -> None:
    out = DATA / "netlib_infeas"
    print(f"netlib-infeas -> {out}")
    emps = build_emps()
    if emps is None:
        print("  cannot build emps decompressor; skipping", file=sys.stderr)
        return
    raw = out / "raw"
    jobs = [lambda n=n: fetch(NETLIB_ORG + "infeas/" + n, raw / n) for n in NETLIB_INFEAS]
    parallel(jobs)
    count = 0
    for n in NETLIB_INFEAS:
        src = raw / n
        dst = out / (n + ".mps")
        if not src.exists():
            continue
        if not dst.exists():
            with open(src, "rb") as fi, open(dst, "wb") as fo:
                subprocess.run([str(emps)], stdin=fi, stdout=fo, check=False)
        count += dst.exists() and dst.stat().st_size > 0
    print(f"  {count}/{len(NETLIB_INFEAS)} decompressed")


def get_miplib3() -> None:
    out = DATA / "miplib3"
    print(f"miplib3 -> {out}")
    jobs = [lambda n=n: fetch(MIPLIB3_RAW + n + ".gz", out / (n + ".mps.gz")) for n in MIPLIB3]
    ok = parallel(jobs)
    fetch(MIPLIB3_RAW + "miplib3.cat", out / "miplib3.cat")
    print(f"  {sum(ok)}/{len(ok)} files")


def get_miplib2017() -> None:
    lst = ROOT / "bench" / "sets" / "miplib2017_subset.txt"
    out = DATA / "miplib2017"
    names = [l.split()[0] for l in lst.read_text().splitlines() if l.strip() and not l.startswith("#")]
    print(f"miplib2017 ({len(names)}) -> {out}")
    jobs = [lambda n=n: fetch(MIPLIB2017 + n + ".mps.gz", out / (n + ".mps.gz")) for n in names]
    ok = parallel(jobs)
    print(f"  {sum(ok)}/{len(ok)} files")


def fnum(v) -> str:
    return repr(float(v))


def mat_to_qps(mat_path: Path, qps_path: Path, name: str) -> None:
    """Convert qpbenchmark .mat (OSQP form: min 1/2 x'Px + q'x + r, l<=Ax<=u, lb<=x<=ub) to QPS."""
    import numpy as np
    import scipy.io
    import scipy.sparse as sp

    m = scipy.io.loadmat(str(mat_path))
    P = sp.csc_matrix(m["P"]).astype(float)
    q = np.asarray(m["q"], dtype=float).ravel()
    r = float(np.asarray(m.get("r", [[0.0]])).ravel()[0]) if "r" in m else 0.0
    A = sp.csc_matrix(m["A"]).astype(float)
    l = np.asarray(m["l"], dtype=float).ravel()
    u = np.asarray(m["u"], dtype=float).ravel()
    lb = np.asarray(m["lb"], dtype=float).ravel() if "lb" in m else np.full(q.size, -np.inf)
    ub = np.asarray(m["ub"], dtype=float).ravel() if "ub" in m else np.full(q.size, np.inf)
    n = q.size
    big = 1e20
    l = np.where(l <= -big, -np.inf, l)
    u = np.where(u >= big, np.inf, u)
    lb = np.where(lb <= -big, -np.inf, lb)
    ub = np.where(ub >= big, np.inf, ub)
    lines = [f"NAME          {name}", "ROWS", " N  OBJ"]
    rtype = []
    for i in range(A.shape[0]):
        if np.isfinite(l[i]) and np.isfinite(u[i]):
            t = "E" if l[i] == u[i] else "L"   # L with RANGES for two-sided
        elif np.isfinite(l[i]):
            t = "G"
        elif np.isfinite(u[i]):
            t = "L"
        else:
            t = "N"
        rtype.append(t)
        lines.append(f" {t}  R{i}")
    lines.append("COLUMNS")
    for j in range(n):
        if q[j] != 0:
            lines.append(f"    C{j}  OBJ  {fnum(q[j])}")
        for k in range(A.indptr[j], A.indptr[j + 1]):
            lines.append(f"    C{j}  R{A.indices[k]}  {fnum(A.data[k])}")
        if q[j] == 0 and A.indptr[j] == A.indptr[j + 1]:
            lines.append(f"    C{j}  OBJ  0.0")
    lines.append("RHS")
    if r != 0:
        lines.append(f"    RHS  OBJ  {fnum(-r)}")
    for i, t in enumerate(rtype):
        if t == "E" or t == "L":
            if np.isfinite(u[i]) and u[i] != 0:
                lines.append(f"    RHS  R{i}  {fnum(u[i])}")
        elif t == "G" and l[i] != 0:
            lines.append(f"    RHS  R{i}  {fnum(l[i])}")
    rng = [(i, u[i] - l[i]) for i, t in enumerate(rtype) if t == "L" and np.isfinite(l[i]) and np.isfinite(u[i])]
    if rng:
        lines.append("RANGES")
        for i, v in rng:
            lines.append(f"    RNG  R{i}  {fnum(v)}")
    lines.append("BOUNDS")
    for j in range(n):
        lo, hi = lb[j], ub[j]
        if lo == hi:
            lines.append(f" FX BND  C{j}  {fnum(lo)}")
            continue
        if not np.isfinite(lo) and not np.isfinite(hi):
            lines.append(f" FR BND  C{j}")
            continue
        if not np.isfinite(lo):
            lines.append(f" MI BND  C{j}")
        elif lo != 0:
            lines.append(f" LO BND  C{j}  {fnum(lo)}")
        if np.isfinite(hi):
            lines.append(f" UP BND  C{j}  {fnum(hi)}")
    Pu = sp.triu(P).tocsc()  # QUADOBJ lists upper triangle once
    if Pu.nnz:
        lines.append("QUADOBJ")
        for j in range(n):
            for k in range(Pu.indptr[j], Pu.indptr[j + 1]):
                i = Pu.indices[k]
                lines.append(f"    C{i}  C{j}  {fnum(Pu.data[k])}")
    lines.append("ENDATA")
    qps_path.write_text("\n".join(lines) + "\n")


def get_maros() -> None:
    out = DATA / "maros"
    raw = out / "mat"
    print(f"maros -> {out}")
    jobs = [lambda n=n: fetch(MAROS_RAW + n + ".mat", raw / (n + ".mat")) for n in MAROS_SMALLISH]
    parallel(jobs)
    fetch(MAROS_RAW + "README.md", out / "README.md")
    count = 0
    for n in MAROS_SMALLISH:
        src = raw / (n + ".mat")
        dst = out / (n + ".qps")
        if src.exists() and not dst.exists():
            try:
                mat_to_qps(src, dst, n)
            except Exception as e:  # noqa: BLE001
                print(f"  convert {n} failed: {e}", file=sys.stderr)
        count += dst.exists()
    print(f"  {count}/{len(MAROS_SMALLISH)} converted")


SETS = {
    "netlib": get_netlib,
    "netlib-infeas": get_netlib_infeas,
    "miplib3": get_miplib3,
    "miplib2017": get_miplib2017,
    "maros": get_maros,
}

if __name__ == "__main__":
    wanted = sys.argv[1:] or ["netlib", "netlib-infeas", "miplib3", "maros"]
    for s in wanted:
        SETS[s]()
