"""PRAMANA terminal console.

An interactive terminal front end for the PRAMANA solver: solve and certify models, analyse their
structure and numerics, debug a run (degeneracy, scaling, refactorizations, cuts, GPU time), re-verify
results independently, draw certified parametric value curves, compare engines and browse benchmarks.

    python -m pramana.tui            (or start.bat / start.sh from the repository root)

No third-party packages: ANSI truecolor + a small line editor (Windows console and POSIX terminals).
Everything is computed by the solver binary (build/pramana[.exe]); this module only drives and renders it.
"""
from __future__ import annotations

import getpass
import glob
import json
import math
import os
import queue
import re
import shlex
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXE = ROOT / "build" / ("pramana.exe" if os.name == "nt" else "pramana")
STATE_DIR = ROOT / "build" / "tui"
STATE_FILE = STATE_DIR / "state.json"
VERSION = "1.0"

# --------------------------------------------------------------------------- colours
# Palette: orange (frames, accents), yellow (proven / highlights), white (text), red (failures).
USE_COLOR = "NO_COLOR" not in os.environ


def _fg(r, g, b):
    return f"\x1b[38;2;{r};{g};{b}m" if USE_COLOR else ""


def _bg(r, g, b):
    return f"\x1b[48;2;{r};{g};{b}m" if USE_COLOR else ""


ORANGE = _fg(255, 140, 26)
AMBER = _fg(255, 176, 46)
YELLOW = _fg(255, 214, 10)
WHITE = _fg(242, 242, 242)
RED = _fg(239, 68, 68)
DIM = _fg(135, 135, 135)
BLACK = _fg(20, 20, 20)
BG_ORANGE = _bg(255, 140, 26)
BG_YELLOW = _bg(255, 214, 10)
BG_WHITE = _bg(235, 235, 235)
BG_RED = _bg(239, 68, 68)
BG_INPUT = _bg(48, 44, 40)
BOLD = "\x1b[1m" if USE_COLOR else ""
ITALIC = "\x1b[3m" if USE_COLOR else ""
RESET = "\x1b[0m" if USE_COLOR else ""
ANSI_RE = re.compile(r"\x1b\[[0-9;?]*[A-Za-z]")

PROVEN = {"OPTIMAL", "INFEASIBLE", "UNBOUNDED"}
LIMITS = {"TIME_LIMIT", "ITERATION_LIMIT", "NODE_LIMIT", "NOT_SOLVED", "INTERRUPTED"}


def c(text, *styles):
    return "".join(styles) + str(text) + RESET if USE_COLOR else str(text)


def status_color(status):
    if status in PROVEN:
        return YELLOW
    if status in LIMITS or status == "NONCONVEX":
        return ORANGE
    return RED


def vlen(s):
    return len(ANSI_RE.sub("", s))


def truncate(s, width):
    """Cut a string with ANSI codes to `width` visible characters."""
    if vlen(s) <= width:
        return s
    out, n, i = [], 0, 0
    while i < len(s) and n < max(0, width - 1):
        m = ANSI_RE.match(s, i)
        if m:
            out.append(m.group())
            i = m.end()
            continue
        out.append(s[i])
        n += 1
        i += 1
    return "".join(out) + "…" + RESET


def pad(s, width, align="left"):
    s = truncate(s, width)
    gap = width - vlen(s)
    if align == "center":
        return " " * (gap // 2) + s + " " * (gap - gap // 2)
    if align == "right":
        return " " * gap + s
    return s + " " * gap


def term_width():
    return max(40, shutil.get_terminal_size((100, 30)).columns)


def panel_width():
    return min(term_width() - 1, 104)


def out(s=""):
    sys.stdout.write(s + "\n")
    sys.stdout.flush()


def fmt(v, digits=10):
    if v is None:
        return "—"
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, (int,)) or (isinstance(v, float) and v.is_integer() and abs(v) < 1e15):
        return f"{int(v):,}"
    if isinstance(v, float):
        if not math.isfinite(v) or abs(v) >= 1e300:
            return "∞" if v > 0 else "-∞"
        return f"{v:.{digits}g}"
    return str(v)


def fsec(t):
    if t is None:
        return "—"
    if t < 1e-3:
        return "<1 ms"
    if t < 1:
        return f"{t * 1e3:.0f} ms"
    if t < 120:
        return f"{t:.2f} s"
    return f"{t / 60:.1f} min"


def sci(v):
    if v is None:
        return "—"
    return f"{v:.1e}"


# --------------------------------------------------------------------------- drawing
def box(rows, width=None, title=None, color=ORANGE, title_style=None):
    width = width or panel_width()
    inner = width - 4
    head = "─" * 3
    if title:
        t = f" {title} "
        top = color + "╭" + head + RESET + (title_style or (BOLD + ORANGE)) + t + RESET + color + \
            "─" * max(0, width - 5 - vlen(t)) + "╮" + RESET
    else:
        top = color + "╭" + "─" * (width - 2) + "╮" + RESET
    lines = [top]
    for r in rows:
        if r == "---":
            lines.append(color + "├" + "─" * (width - 2) + "┤" + RESET)
        else:
            lines.append(color + "│ " + RESET + pad(r, inner) + color + " │" + RESET)
    lines.append(color + "╰" + "─" * (width - 2) + "╯" + RESET)
    for line in lines:
        out(line)


def section(title, sub=""):
    return c(title, BOLD, ORANGE) + (c("  " + sub, DIM) if sub else "")


def bar(frac, width, color=YELLOW, rest=DIM):
    frac = 0.0 if frac is None or not math.isfinite(frac) else max(0.0, min(1.0, frac))
    parts = "▏▎▍▌▋▊▉"
    full = int(frac * width)
    rem = frac * width - full
    s = "█" * full
    if full < width and rem > 1 / 8:
        s += parts[min(6, int(rem * 8) - 1)]
    return c(s, color) + c("░" * (width - vlen(s)), rest)


def kv(label, value, lw=16):
    return c(pad(label, lw), DIM) + value


def verdict_mark(v):
    v = (v or "").upper()
    if v == "PASS":
        return c("✓", BOLD, YELLOW)
    if v == "WARN":
        return c("!", BOLD, ORANGE)
    if v in ("N/A", "NA"):
        return c("·", DIM)
    return c("✗", BOLD, RED)


def finding(level, text):
    mark = {"ok": c("✓", BOLD, YELLOW), "info": c("•", ORANGE), "warn": c("▲", BOLD, ORANGE),
            "bad": c("✗", BOLD, RED)}[level]
    color = {"ok": WHITE, "info": WHITE, "warn": WHITE, "bad": RED}[level]
    return f"{mark} {c(text, color)}"


def wrap(text, width, indent=""):
    words, lines, cur = text.split(), [], ""
    for w in words:
        if vlen(cur) + len(w) + 1 > width and cur:
            lines.append(cur)
            cur = indent + w
        else:
            cur = (cur + " " + w) if cur else w
    if cur:
        lines.append(cur)
    return lines


# --------------------------------------------------------------------------- terminal setup / keys
def enable_vt():
    if os.name == "nt":
        try:
            import ctypes
            k = ctypes.windll.kernel32
            for std in (-11, -12):
                h = k.GetStdHandle(std)
                mode = ctypes.c_uint32()
                if k.GetConsoleMode(h, ctypes.byref(mode)):
                    k.SetConsoleMode(h, mode.value | 0x0004)
            k.SetConsoleOutputCP(65001)
        except Exception:  # noqa: BLE001
            pass
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass


class Keys:
    """Blocking single-key reader returning names ('enter', 'up', ...) or ('char', c)."""

    def __init__(self):
        self.win = os.name == "nt"
        if not self.win:
            import termios
            import tty  # noqa: F401
            self.fd = sys.stdin.fileno()
            self.old = termios.tcgetattr(self.fd)

    def pending(self):
        if self.win:
            import msvcrt
            return msvcrt.kbhit()
        import select
        return bool(select.select([sys.stdin], [], [], 0)[0])

    def read(self):
        if self.win:
            import ctypes
            import msvcrt
            k32 = ctypes.windll.kernel32
            h = k32.GetStdHandle(-10)
            mode = ctypes.c_uint32()
            have = k32.GetConsoleMode(h, ctypes.byref(mode))
            if have:  # read Ctrl+C as a key (0x03) instead of a signal while waiting at the prompt
                k32.SetConsoleMode(h, mode.value & ~0x0001)
            try:
                ch = msvcrt.getwch()
                if ch in ("\x00", "\xe0"):
                    code = msvcrt.getwch()
                    return {"H": "up", "P": "down", "K": "left", "M": "right", "G": "home", "O": "end",
                            "S": "delete"}.get(code, "none")
                return self._map(ch)
            finally:
                if have:
                    k32.SetConsoleMode(h, mode.value)
        import termios
        import tty
        tty.setcbreak(self.fd)
        try:
            new = termios.tcgetattr(self.fd)
            new[3] &= ~(termios.ECHO | termios.ISIG)
            termios.tcsetattr(self.fd, termios.TCSANOW, new)
            ch = os.read(self.fd, 1).decode("utf-8", "replace")
            if ch == "\x1b":
                import select
                if select.select([self.fd], [], [], 0.03)[0]:
                    seq = os.read(self.fd, 2).decode("utf-8", "replace")
                    if seq.startswith("["):
                        code = seq[1:]
                        if code in "ABCDHF":
                            return {"A": "up", "B": "down", "C": "right", "D": "left", "H": "home", "F": "end"}[code]
                        if code == "3":
                            os.read(self.fd, 1)
                            return "delete"
                    return "none"
                return "esc"
            if ord(ch[0]) >= 0xC0:  # rest of a UTF-8 sequence
                pass
            return self._map(ch)
        finally:
            termios.tcsetattr(self.fd, termios.TCSADRAIN, self.old)

    @staticmethod
    def _map(ch):
        return {"\r": "enter", "\n": "enter", "\x08": "backspace", "\x7f": "backspace", "\t": "tab",
                "\x1b": "esc", "\x03": "ctrl-c", "\x04": "ctrl-d", "\x0c": "ctrl-l", "\x15": "ctrl-u",
                "\x17": "ctrl-w", "\x01": "home", "\x05": "end"}.get(ch, ("char", ch))


# --------------------------------------------------------------------------- models, state
def model_index():
    order = ["netlib", "miplib3", "maros", "gen", "adversarial", "netlib_infeas"]
    idx = {}
    pats = ["data/*/*.mps", "data/*/*.mps.gz", "data/*/*.qps", "data/*/*.QPS", "data/*/*.SIF", "tests/data/*.mps",
            "tests/data/*.qps"]
    files = []
    for p in pats:
        files += glob.glob(str(ROOT / p))
    files.sort(key=lambda f: (order.index(Path(f).parent.name) if Path(f).parent.name in order else 99, f))
    for f in files:
        name = Path(f).name
        for ext in (".mps.gz", ".mps", ".qps", ".QPS", ".SIF"):
            if name.endswith(ext):
                name = name[: -len(ext)]
                break
        idx.setdefault(name.lower(), f)
    return idx


def load_state():
    try:
        return json.loads(STATE_FILE.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001
        return {"runs": [], "history": []}


def save_state(st):
    try:
        STATE_DIR.mkdir(parents=True, exist_ok=True)
        st["runs"] = st["runs"][-40:]
        st["history"] = st["history"][-200:]
        STATE_FILE.write_text(json.dumps(st, indent=1), encoding="utf-8")
    except Exception:  # noqa: BLE001
        pass


# --------------------------------------------------------------------------- process runner with spinner
SPIN = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"


def run_process(args, label, cwd=None, show_log=True):
    """Run a command, stream its output into a live spinner line. Returns (rc, lines, seconds)."""
    t0 = time.perf_counter()
    try:
        proc = subprocess.Popen([str(a) for a in args], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                cwd=str(cwd or ROOT), text=True, encoding="utf-8", errors="replace", bufsize=1)
    except OSError as e:
        return 127, [str(e)], 0.0
    q = queue.Queue()

    def reader():
        for line in proc.stdout:
            q.put(line.rstrip("\n"))
        q.put(None)

    threading.Thread(target=reader, daemon=True).start()
    lines, last, k, done = [], "", 0, False
    try:
        while not done:
            try:
                while True:
                    item = q.get(timeout=0.08)
                    if item is None:
                        done = True
                        break
                    lines.append(item)
                    if item.strip():
                        last = item.strip()
            except queue.Empty:
                pass
            el = time.perf_counter() - t0
            w = term_width() - 1
            head = f"{c(SPIN[k % len(SPIN)], BOLD, ORANGE)} {c(label, WHITE)} {c(f'{el:5.1f}s', AMBER)}"
            tail = ("  " + c(last, DIM)) if show_log and last else ""
            sys.stdout.write("\r\x1b[2K" + truncate(head + tail, w))
            sys.stdout.flush()
            k += 1
        proc.wait()
    except KeyboardInterrupt:
        proc.kill()
        proc.wait()
        sys.stdout.write("\r\x1b[2K")
        out(c("  ✗ interrupted", RED))
        return 130, lines, time.perf_counter() - t0
    sys.stdout.write("\r\x1b[2K")
    sys.stdout.flush()
    return proc.returncode, lines, time.perf_counter() - t0


# --------------------------------------------------------------------------- analysis (findings)
def result_findings(j):
    f = []
    cert = j.get("certificate", {})
    tel = j.get("telemetry", {})
    st, est = j.get("status"), j.get("engine_status")
    if cert and not cert.get("accepted", False) and st not in LIMITS:
        f.append(("bad", f"Certificate REJECTED: the engine's claim ({est}) could not be proven, "
                         f"so the reported status is {st}."))
    if est and st and est != st:
        f.append(("warn", f"Engine claimed {est}; the independent certifier downgraded it to {st}."))
    for key in ("simplex", "postsolve_cleanup"):
        s = tel.get(key)
        if not s:
            continue
        name = "Simplex" if key == "simplex" else "Cleanup simplex"
        deg = s.get("degenerate_fraction") or 0
        if deg >= 0.3:
            f.append(("warn", f"{name}: highly degenerate, {deg:.0%} of {fmt(s.get('iterations'))} pivots made no "
                              f"progress (cost perturbations: {fmt(s.get('perturbations'))}, bound flips: "
                              f"{fmt(s.get('bound_flips'))})."))
        rb, ra = s.get("scaling_range_before"), s.get("scaling_range_after")
        if rb and ra and rb >= 1e6:
            f.append(("warn", f"{name}: badly scaled input, coefficient range {sci(rb)} reduced to {sci(ra)} by "
                              f"scaling."))
        if s.get("numerical_refactors"):
            f.append(("warn", f"{name}: {s['numerical_refactors']} refactorizations were triggered by residual "
                              f"checks (LU accuracy loss caught and repaired)."))
        if s.get("singular_basis_repairs"):
            f.append(("warn", f"{name}: {s['singular_basis_repairs']} singular bases repaired with slack columns."))
        if s.get("stall_recoveries"):
            f.append(("warn", f"{name}: stalling detected {s['stall_recoveries']}× and broken by re-perturbation."))
        if s.get("rejected_farkas_rays"):
            f.append(("info", f"{name}: {s['rejected_farkas_rays']} infeasibility ray(s) rejected by "
                              f"certify-before-claim."))
        if s.get("rejected_pivots"):
            f.append(("info", f"{name}: {s['rejected_pivots']} unstable pivots rejected (tiny pivot elements)."))
        if (s.get("max_cost_shift") or 0) > 1e-6:
            f.append(("info", f"{name}: cost shifting up to {sci(s['max_cost_shift'])} (removed before the "
                              f"final answer)."))
    pre = tel.get("presolve", {})
    if pre.get("status") == "reduced":
        f.append(("info", f"Presolve removed {fmt(pre.get('rows_removed'))} rows / {fmt(pre.get('cols_removed'))} "
                          f"columns, tightened {fmt(pre.get('bounds_tightened'))} bounds and "
                          f"{fmt(pre.get('coefficients_tightened'))} coefficients."))
    mip = tel.get("mip")
    if mip:
        r0, r1 = mip.get("root_lp_bound"), mip.get("root_bound_after_cuts")
        best = j.get("objective")
        if r0 is not None and r1 is not None and best is not None and abs(best - r0) > 1e-9:
            closed = (r1 - r0) / (best - r0) if best != r0 else 0
            if math.isfinite(closed):
                f.append(("info", f"Cutting planes closed {max(0.0, closed):.0%} of the root gap "
                                  f"(GMI {mip.get('cuts_gomory', 0)}, MIR {mip.get('cuts_mir', 0)}, "
                                  f"cover {mip.get('cuts_cover', 0)})."))
        gap = j.get("relative_gap")
        if st in LIMITS and gap is not None and math.isfinite(gap):
            f.append(("warn", f"Stopped at the limit with a {gap:.2%} gap; the bound is safe "
                              f"(Neumaier–Shcherbina, {mip.get('safe_bound_corrections', 0)} corrections)."))
        if mip.get("cuts_rejected_numerics"):
            f.append(("info", f"{mip['cuts_rejected_numerics']} cuts rejected as numerically unsafe."))
        if mip.get("incumbent_source"):
            f.append(("info", f"Best solution found by: {mip['incumbent_source']} "
                              f"(first incumbent after {fsec(mip.get('first_incumbent_seconds'))})."))
    pd = tel.get("pdhg")
    if pd:
        k, tr = pd.get("kernel_seconds") or 0, pd.get("transfer_seconds") or 0
        f.append(("info", f"PDHG on {pd.get('backend', '?').upper()}: {fmt(pd.get('iterations'))} iterations, "
                          f"{pd.get('restarts', 0)} restarts, relative KKT {sci(pd.get('relative_kkt'))}; kernels "
                          f"{fsec(k)}, transfers {fsec(tr)}."))
    rt = tel.get("router", {})
    pred = (rt.get("predicted_seconds") or {}).get(rt.get("selected", ""))
    actual = tel.get("timing", {}).get("engine")
    if pred and actual and actual > 0.05:
        ratio = max(pred / actual, actual / pred)
        if ratio > 5:
            f.append(("warn", f"Router predicted {fsec(pred)} for {rt.get('selected')} but it took {fsec(actual)} "
                              f"(×{ratio:.0f} off)."))
    if not any(l in ("warn", "bad") for l, _ in f):
        f.insert(0, ("ok", "No numerical issues detected."))
    return f


def analysis_findings(a):
    f = []
    num, st, rt = a.get("numerics", {}), a.get("structure", {}), a.get("router", {})
    rng = num.get("matrix_range") or 1
    if rng >= 1e10:
        f.append(("bad", f"Extreme coefficient range {sci(rng)}: tolerances become meaningless without scaling; "
                         f"expect the certifier to matter."))
    elif rng >= 1e6:
        f.append(("warn", f"Wide coefficient range {sci(rng)}: ill-conditioning likely; geometric scaling applied "
                          f"automatically."))
    else:
        f.append(("ok", f"Coefficient range {sci(rng)}: well scaled."))
    ints, bins = a.get("integers", 0), a.get("binaries", 0)
    if ints and (num.get("matrix_max_abs") or 0) >= 1e3 and bins:
        f.append(("warn", f"Large coefficients (up to {sci(num.get('matrix_max_abs'))}) next to {bins} binaries: "
                          f"likely big-M constraints and a weak LP relaxation; presolve tightens them."))
    eqf = a.get("row_types", {}).get("equality", 0) / max(1, a.get("rows", 1))
    if eqf > 0.8 and a.get("rows", 0) > 1000:
        f.append(("info", f"{eqf:.0%} equality rows (flow/balance structure): typical of planning models; "
                          f"prone to dual degeneracy."))
    if st.get("empty_rows") or st.get("empty_cols"):
        f.append(("info", f"{st.get('empty_rows', 0)} empty rows and {st.get('empty_cols', 0)} empty columns "
                          f"(removed by presolve)."))
    if st.get("singleton_rows"):
        f.append(("info", f"{st['singleton_rows']} singleton rows become bounds in presolve."))
    fr = a.get("col_bounds", {}).get("free", 0)
    if fr:
        f.append(("info", f"{fr} free columns (never at a bound; handled as free nonbasics)."))
    eng = rt.get("engine")
    reason = rt.get("reason", "")
    if eng:
        f.append(("info", f"Router: {eng} for this model — {reason}"))
    return f


# --------------------------------------------------------------------------- the app
COMMANDS = [
    ("/solve", "<model> [flags]", "solve and certify a model (default for a bare model name)", True),
    ("/analyze", "<model>", "structure, numerics and engine prediction, without solving", True),
    ("/debug", "[run#]", "numerical diagnostics of the last (or a given) run", False),
    ("/verify", "[exact]", "re-check the last result independently (C++ + exact rationals)", False),
    ("/compare", "<model>", "run every engine on an LP and compare certified times", True),
    ("/param", "<model> <col> <cost|upper|lower> <from> <to>", "certified parametric value curve", True),
    ("/bench", "", "benchmark results vs HiGHS", False),
    ("/gpu", "[calibrate]", "GPU report and measured CPU/GPU crossover", False),
    ("/models", "[filter]", "list available benchmark and industrial models", False),
    ("/set", "<time|algo|threads|gpu> <value>", "session settings", True),
    ("/check", "", "full self-check (build, tests, audit, certified solves)", False),
    ("/history", "", "recent runs", False),
    ("/clear", "", "clear the screen", False),
    ("/help", "", "commands and keys", False),
    ("/exit", "", "quit", False),
]
MODEL_ARG = {"/solve", "/analyze", "/compare", "/param"}
TIPS = ['Try "afiro"', 'Try "/analyze pilot87"', 'Try "/compare plan_20x12"', 'Try "p0201"',
        'Try "/param plan_10x12 BUY_CR03_0 upper 0 120"', 'Try "/bench"', 'Try "scaled_afiro" then "/debug"']


class App:
    def __init__(self):
        self.state = load_state()
        self.models = model_index()
        self.settings = {"time": 300.0, "algo": "auto", "threads": 0, "gpu": True}
        self.last = None  # last run record
        self.sysinfo = {"threads": "?", "gpu": None}
        self.tip = 0
        self.keys = None
        self.area = False  # input area currently drawn

    # ----- system info
    def probe(self):
        if not EXE.exists():
            return
        try:
            p = subprocess.run([str(EXE), "info"], capture_output=True, text=True, timeout=30,
                               encoding="utf-8", errors="replace")
            for line in p.stdout.splitlines():
                if line.strip().startswith("threads"):
                    self.sysinfo["threads"] = line.split(":", 1)[1].strip()
                if line.strip().startswith("gpu "):
                    g = line.split(":", 1)[1].strip()
                    self.sysinfo["gpu"] = None if g.startswith("unavailable") else g
        except Exception:  # noqa: BLE001
            pass

    # ----- welcome screen
    def welcome(self):
        out()
        w = panel_width()
        two = w >= 92
        lw = 58 if two else w - 4
        rw = w - 4 - lw - 3 if two else 0
        user = os.environ.get("USERNAME") or os.environ.get("USER") or getpass.getuser()
        O, Y, R = ORANGE, YELLOW, RED
        logo = [
            f"   {O}▁▁▁▁▁{RESET}{BOLD}{Y}◆{RESET} {R}↗{RESET}  ",
            f" {O}╱{AMBER}░░░░░░{O}╲{RESET}    ",
            f"{O}▕{AMBER}░░░░░░░░{O}▏{RESET}   ",
            f" {O}╲{AMBER}░░░░░░{O}╱{RESET}    ",
            f"   {O}▔▔▔▔▔{RESET}      ",
        ]
        word = [["█▀█", "█▀█", "▄▀█", "█▀▄▀█", "▄▀█", "█▄ █", "▄▀█"],
                ["█▀▀", "█▀▄", "█▀█", "█ ▀ █", "█▀█", "█ ▀█", "█▀█"]]
        shades = [ORANGE, ORANGE, AMBER, AMBER, YELLOW, YELLOW, YELLOW]
        wordmark = [" ".join(BOLD + shades[i] + row[i] + RESET for i in range(7)) for row in word]
        art = [
            logo[0] + "",
            logo[1] + wordmark[0],
            logo[2] + wordmark[1],
            logo[3] + c("certified optimization engine", ITALIC, WHITE),
            logo[4],
        ]
        gpu = (self.sysinfo["gpu"] or "no GPU (CPU only)").replace("NVIDIA GeForce ", "").replace("NVIDIA ", "")
        left = ["", pad(c(f"Welcome back, {user}!", BOLD, WHITE), lw, "center"), ""]
        left += [pad(a, lw, "center") for a in art]
        left += ["",
                 pad(c("LP · MILP · QP", YELLOW) + c("  ·  ", DIM) + c(f"{self.sysinfo['threads']} threads", WHITE)
                     + c("  ·  ", DIM) + c(gpu, WHITE), lw, "center"),
                 pad(c(str(ROOT), DIM), lw, "center"), ""]
        right = [c("Getting started", BOLD, ORANGE),
                 c("Type a model name, e.g. ", WHITE) + c("afiro", YELLOW),
                 c("/", YELLOW) + c(" for commands · ", WHITE) + c("Tab", YELLOW) + c(" completes", WHITE),
                 c("/analyze", YELLOW) + c(" · ", DIM) + c("/debug", YELLOW) + c(" · ", DIM) + c("/verify", YELLOW),
                 c("─" * rw, ORANGE),
                 c("Recent activity", BOLD, ORANGE)]
        runs = self.state.get("runs", [])[-6:][::-1]
        if not runs:
            right.append(c("No recent activity", DIM))
        for r in runs:
            st = r.get("status", "?")
            right.append(pad(c(r.get("model", "?"), WHITE), max(8, rw - 22)) + " " +
                         pad(c(st, status_color(st)), 13) + " " + c(fsec(r.get("seconds")), DIM))
        rows = []
        if two:
            n = max(len(left), len(right))
            left += [""] * (n - len(left))
            right += [""] * (n - len(right))
            for a, b in zip(left, right):
                rows.append(pad(a, lw) + c(" │ ", ORANGE) + pad(b, rw))
        else:
            rows = left + ["---"] + right
        box(rows, w, title=f"PRAMANA v{VERSION}", title_style=BOLD + ORANGE)
        out()

    # ----- input area
    def suggestions(self, buf):
        if buf.startswith("/") and " " not in buf:
            return [("cmd", n, a, d) for n, a, d, _ in COMMANDS if n.startswith(buf.lower())][:9]
        toks = buf.split(" ")
        if len(toks) == 2 and toks[0] in MODEL_ARG and toks[1]:
            pre = toks[1].lower()
            return [("model", n, "", Path(p).parent.name) for n, p in self.models.items() if n.startswith(pre)][:7]
        if not buf.startswith("/") and " " not in buf and len(buf) >= 2:
            pre = buf.lower()
            return [("model", n, "", Path(p).parent.name) for n, p in self.models.items() if n.startswith(pre)][:7]
        return []

    def status_bar(self, width):
        def seg(text, bgc, fgc=BLACK):
            return bgc + fgc + BOLD + f" {text} " + RESET

        left = seg("PRAMANA", BG_ORANGE) + seg(ROOT.name, BG_YELLOW) + \
            seg(f"engine {self.settings['algo']} · {self.settings['time']:g}s", BG_WHITE)
        right = ""
        if self.last:
            st = self.last["status"]
            bgc = BG_YELLOW if st in PROVEN else (BG_ORANGE if st in LIMITS else BG_RED)
            fgc = BLACK if st in PROVEN or st in LIMITS else WHITE
            right += seg(f"● {self.last['model']} {st} {fsec(self.last.get('seconds'))}", bgc, fgc) + " "
        g = self.sysinfo["gpu"]
        gpu_txt = ("GPU " + g.replace("NVIDIA GeForce ", "")) if (g and self.settings["gpu"]) else "CPU only"
        right += seg(f"{gpu_txt} · {self.sysinfo['threads']} thr", BG_ORANGE)
        gap = width - vlen(left) - vlen(right)
        if gap < 1:
            return truncate(left, width)
        return left + " " * gap + right

    def draw_area(self, buf, cur, sel):
        w = term_width() - 1
        rule = c("─" * w, ORANGE)
        prompt = c("> ", BOLD, ORANGE)
        room = w - 3
        if buf:
            start = max(0, cur - room + 1)
            view = buf[start:start + room]
            line = prompt + c(view, WHITE)
            col = 2 + (cur - start)
        else:
            line = prompt + c(TIPS[self.tip % len(TIPS)], DIM)
            col = 2
        lines = [rule, line, rule, self.status_bar(w)]
        for i, (kind, name, args, desc) in enumerate(self.suggestions(buf)):
            hl = i == sel
            mark = c("▸ ", BOLD, YELLOW) if hl else "  "
            nm = c(pad(name, 12), BOLD, YELLOW) if hl else c(pad(name, 12), ORANGE if kind == "cmd" else WHITE)
            extra = (c(args + "  ", WHITE) if args else "") + c(desc, DIM)
            lines.append(truncate(mark + nm + "  " + extra, w))
        s = "\r\x1b[J" + "\n".join(truncate(x, w) for x in lines)
        up = len(lines) - 2
        s += (f"\x1b[{up}A" if up > 0 else "") + "\r" + (f"\x1b[{col}C" if col > 0 else "")
        if self.area:
            s = "\x1b[1A" + s
        sys.stdout.write(s)
        sys.stdout.flush()
        self.area = True

    def clear_area(self):
        if self.area:
            sys.stdout.write("\x1b[1A\r\x1b[J")
            sys.stdout.flush()
            self.area = False

    def read_line(self):
        hist = self.state.setdefault("history", [])
        buf, cur, sel, hpos, armed_exit = "", 0, 0, len(hist), False
        self.draw_area(buf, cur, sel)
        while True:
            try:
                k = self.keys.read()
            except KeyboardInterrupt:  # Windows delivers Ctrl+C as a signal, not a key
                k = "ctrl-c"
            sugg = self.suggestions(buf)
            if k == "enter":
                if sugg:
                    kind, name = sugg[min(sel, len(sugg) - 1)][:2]
                    if kind == "cmd" and name != buf.strip().lower():
                        takes = next(t for n, _, _, t in COMMANDS if n == name)
                        if takes:  # complete the command and wait for its arguments
                            buf, cur, sel = name + " ", len(name) + 1, 0
                            self.draw_area(buf, cur, sel)
                            continue
                        buf = name
                    elif kind == "model":
                        toks = buf.split(" ")
                        if not self.resolve(toks[-1]):
                            toks[-1] = name
                            buf = " ".join(toks)
                self.clear_area()
                return buf.strip()
            if k == "ctrl-c" or k == "ctrl-d":
                if buf:
                    buf, cur, sel = "", 0, 0
                elif armed_exit or k == "ctrl-d":
                    self.clear_area()
                    return "/exit"
                else:
                    armed_exit = True
                    self.clear_area()
                    out(c("  Press Ctrl+C again to exit", DIM))
            elif k == "esc":
                buf, cur, sel = "", 0, 0
            elif k == "ctrl-l":
                self.clear_area()
                sys.stdout.write("\x1b[2J\x1b[H")
                self.welcome()
            elif k == "ctrl-u":
                buf, cur = buf[cur:], 0
            elif k == "ctrl-w":
                left = buf[:cur].rstrip()
                i = left.rfind(" ") + 1
                buf, cur = buf[:i] + buf[cur:], i
            elif k == "backspace":
                if cur > 0:
                    buf, cur = buf[:cur - 1] + buf[cur:], cur - 1
                sel = 0
            elif k == "delete":
                buf = buf[:cur] + buf[cur + 1:]
            elif k == "left":
                cur = max(0, cur - 1)
            elif k == "right":
                cur = min(len(buf), cur + 1)
            elif k == "home":
                cur = 0
            elif k == "end":
                cur = len(buf)
            elif k == "up":
                if sugg:
                    sel = (sel - 1) % len(sugg)
                elif hist and hpos > 0:
                    hpos -= 1
                    buf = hist[hpos]
                    cur = len(buf)
            elif k == "down":
                if sugg:
                    sel = (sel + 1) % len(sugg)
                elif hpos < len(hist):
                    hpos += 1
                    buf = hist[hpos] if hpos < len(hist) else ""
                    cur = len(buf)
            elif k == "tab":
                if sugg:
                    kind, name = sugg[min(sel, len(sugg) - 1)][:2]
                    if kind == "cmd":
                        buf = name + " "
                    else:
                        toks = buf.split(" ")
                        toks[-1] = name
                        buf = " ".join(toks) + " "
                    cur, sel = len(buf), 0
            elif isinstance(k, tuple):
                ch = k[1]
                if ch.isprintable():
                    buf, cur = buf[:cur] + ch + buf[cur:], cur + 1
                    sel = 0
                    armed_exit = False
            if not self.keys.pending():
                self.draw_area(buf, cur, sel)

    # ----- helpers
    def resolve(self, name):
        p = Path(name)
        if p.exists():
            return str(p.resolve())
        if (ROOT / name).exists():
            return str(ROOT / name)
        return self.models.get(name.lower()) or self.models.get(Path(name).stem.lower())

    def need_model(self, args, usage):
        if not args:
            out(c(f"  usage: {usage}", ORANGE))
            return None
        path = self.resolve(args[0])
        if not path:
            near = [n for n in self.models if n.startswith(args[0].lower()[:3])][:6]
            out(c(f"  ✗ model '{args[0]}' not found", RED) +
                (c("   did you mean: ", DIM) + c(", ".join(near), YELLOW) if near else ""))
            return None
        return path

    def engine_flags(self):
        f = ["--time", f"{self.settings['time']:g}"]
        if self.settings["algo"] != "auto":
            f += ["--algo", self.settings["algo"]]
        if self.settings["threads"]:
            f += ["--threads", str(self.settings["threads"])]
        if not self.settings["gpu"]:
            f += ["--no-gpu"]
        return f

    # ----- commands
    def cmd_solve(self, args, quiet=False):
        path = self.need_model(args, "/solve <model> [--algo dual|ipm|pdhg-gpu|race] [--time s] [--gap g] ...")
        if not path:
            return None
        STATE_DIR.mkdir(parents=True, exist_ok=True)
        n = len(self.state["runs"]) + 1
        name = Path(path).name.split(".")[0]
        jf = STATE_DIR / f"run{n:03d}_{name}.json"
        extra = args[1:]
        flags = self.engine_flags()
        if "--algo" in extra and "--algo" in flags:
            i = flags.index("--algo")
            del flags[i:i + 2]
        rc, lines, secs = run_process([EXE, path, "--json", jf, "--vectors", "--log", "1"] + flags + extra,
                                      f"Solving {name}")
        if rc == 130:
            return None
        if not jf.exists():
            out(c(f"  ✗ solver error (exit {rc})", BOLD, RED))
            for l in lines[-8:]:
                out("    " + c(l, RED))
            return None
        j = json.loads(jf.read_text(encoding="utf-8"))
        rec = {"n": n, "model": name, "path": path, "json": str(jf), "status": j.get("status"),
               "objective": j.get("objective"), "seconds": j.get("seconds"), "engine": j.get("engine"),
               "when": time.strftime("%Y-%m-%d %H:%M")}
        self.state["runs"].append(rec)
        self.last = rec
        save_state(self.state)
        if not quiet:
            self.render_result(j, rec)
        return j

    def render_result(self, j, rec):
        tel = j.get("telemetry", {})
        md = tel.get("model", {})
        st = j.get("status", "?")
        cert = j.get("certificate", {})
        col = status_color(st)
        cls = "MILP" if md.get("is_mip") else ("QP" if md.get("is_qp") else "LP")
        rows = []
        certified = cert.get("accepted")
        badge = c(" certified ✓ ", BOLD, BLACK, BG_YELLOW) if certified else c(" not certified ✗ ", BOLD, WHITE,
                                                                                  BG_RED)
        rows.append(kv("Status", c(f"● {st}", BOLD, col) + "   " + badge +
                       (c(f"   engine claimed {j.get('engine_status')}", ORANGE)
                        if j.get("engine_status") not in (None, st) else "")))
        if j.get("objective") is not None and st not in ("INFEASIBLE",):
            rows.append(kv("Objective", c(fmt(j.get("objective"), 12), BOLD, WHITE)))
        if j.get("best_bound") is not None and math.isfinite(j.get("best_bound") or float("nan")):
            rows.append(kv("Certified bound", c(fmt(j.get("best_bound"), 12), WHITE) +
                           c(f"   gap {sci(j.get('relative_gap'))}", DIM)))
        rt = tel.get("router", {})
        rows.append(kv("Engine", c(j.get("engine", "?"), YELLOW) +
                       (c(f"   router: {rt.get('reason', '')}", DIM) if rt.get("requested") == "auto" else "")))
        rows.append(kv("Model", c(f"{md.get('name', rec['model'])}", WHITE) + c(
            f"   {cls} · {fmt(md.get('rows'))} rows · {fmt(md.get('cols'))} cols · {fmt(md.get('nnz'))} nnz"
            + (f" · {fmt(md.get('integers'))} integer" if md.get("integers") else ""), DIM)))
        timing = tel.get("timing", {})
        total = timing.get("total") or j.get("seconds") or 0
        rows.append(kv("Time", c(fsec(total), BOLD, WHITE) + "   " + self.timing_bar(timing, 30)))
        mip = tel.get("mip")
        if mip:
            rows.append(kv("Search", c(f"{fmt(mip.get('nodes'))} nodes · {fmt(mip.get('lp_iterations'))} LP iterations"
                                       f" · depth {fmt(mip.get('max_depth'))} · cuts "
                                       f"{(mip.get('cuts_gomory') or 0) + (mip.get('cuts_mir') or 0) + (mip.get('cuts_cover') or 0)}",
                                       WHITE)))
        elif tel.get("simplex"):
            s = tel["simplex"]
            rows.append(kv("Simplex", c(f"{fmt(s.get('iterations'))} iterations · {fmt(s.get('refactorizations'))} "
                                        f"refactorizations · {s.get('degenerate_fraction', 0):.0%} degenerate", WHITE)))
        rows.append("---")
        rows.append(section("Certificate", "independent re-check on the original, unscaled model"))
        for ch in cert.get("checks", []):
            rows.append(f"  {verdict_mark(ch.get('verdict'))} " + pad(c(ch.get("check", ""), WHITE), 30) +
                        c(pad(sci(ch.get("value")), 10, "right"), YELLOW if ch.get("verdict") == "PASS" else ORANGE)
                        + c(f"  tol {sci(ch.get('tolerance'))}", DIM) + c(f"  {ch.get('detail', '')}", DIM))
        rows.append("---")
        rows.append(section("Analysis"))
        for lvl, text in result_findings(j)[:6]:
            for i, line in enumerate(wrap(text, panel_width() - 8)):
                rows.append("  " + (finding(lvl, line) if i == 0 else "  " + c(line, WHITE if lvl != "bad" else RED)))
        box(rows, title=f"Run #{rec['n']} · {rec['model']}")
        out(c("  /debug", YELLOW) + c(" numerical diagnostics   ", DIM) + c("/verify", YELLOW) +
            c(" independent re-check   ", DIM) + c(f"JSON: {Path(rec['json']).relative_to(ROOT)}", DIM))

    @staticmethod
    def timing_bar(t, width):
        parts = [("presolve", ORANGE), ("engine", YELLOW), ("postsolve", AMBER), ("certify", WHITE)]
        tot = sum(t.get(k) or 0 for k, _ in parts)
        if tot <= 0:
            return ""
        s, legend = "", []
        for k, colr in parts:
            v = t.get(k) or 0
            n = int(round(v / tot * width))
            s += c("█" * n, colr)
            if v / tot >= 0.02:
                legend.append(c("■", colr) + c(f" {k} {fsec(v)}", DIM))
        return s + "  " + "  ".join(legend)

    def get_run(self, args):
        runs = self.state.get("runs", [])
        if args:
            try:
                n = int(args[0].lstrip("#"))
                return next((r for r in runs if r.get("n") == n), None)
            except ValueError:
                return None
        return self.last or (runs[-1] if runs else None)

    def cmd_debug(self, args):
        rec = self.get_run(args)
        if not rec or not Path(rec["json"]).exists():
            out(c("  no run to debug yet — solve a model first (e.g. afiro)", ORANGE))
            return
        j = json.loads(Path(rec["json"]).read_text(encoding="utf-8"))
        tel = j.get("telemetry", {})
        rows = [section("Findings")]
        for lvl, text in result_findings(j):
            for i, line in enumerate(wrap(text, panel_width() - 8)):
                rows.append("  " + (finding(lvl, line) if i == 0 else "  " + c(line, WHITE)))
        w2 = (panel_width() - 8) // 2

        def grid(title, items):
            if not items:
                return
            rows.append("---")
            rows.append(section(title))
            items = [(k, v) for k, v in items if v is not None]
            for i in range(0, len(items), 2):
                cells = [pad(c(pad(k, 24), DIM) + c(fmt(v, 6), WHITE), w2) for k, v in items[i:i + 2]]
                rows.append("  " + "  ".join(cells))

        s = tel.get("simplex")
        if s:
            grid("Simplex", [("iterations", s.get("iterations")), ("phase 1 / 2", f"{fmt(s.get('dual_phase1_iterations'))} / {fmt(s.get('dual_phase2_iterations'))}"),
                             ("primal cleanup", s.get("primal_iterations")), ("degenerate pivots", s.get("degenerate_pivots")),
                             ("bound flips", s.get("bound_flips")), ("refactorizations", s.get("refactorizations")),
                             ("numerical refactors", s.get("numerical_refactors")), ("rejected pivots", s.get("rejected_pivots")),
                             ("singular repairs", s.get("singular_basis_repairs")), ("stall recoveries", s.get("stall_recoveries")),
                             ("perturbations", s.get("perturbations")), ("max perturbation", s.get("max_perturbation")),
                             ("cost shifts", s.get("cost_shifts")), ("max cost shift", s.get("max_cost_shift")),
                             ("scaling range before", s.get("scaling_range_before")), ("scaling range after", s.get("scaling_range_after"))])
            deg = s.get("degenerate_fraction") or 0
            rows.append("  " + c(pad("degeneracy", 24), DIM) + bar(deg, 30, ORANGE if deg >= 0.3 else YELLOW) +
                        c(f" {deg:.0%}", WHITE))
            prof = s.get("profile_seconds") or {}
            if prof:
                tot = sum(v for v in prof.values() if isinstance(v, (int, float))) or 1
                rows.append("  " + c("where simplex time went", DIM))
                for k, v in sorted(prof.items(), key=lambda kv2: -kv2[1])[:6]:
                    rows.append("    " + c(pad(k, 20), WHITE) + bar(v / tot, 30) + c(f" {fsec(v)}", DIM))
        m = tel.get("mip")
        if m:
            grid("Branch-and-cut", [("nodes", m.get("nodes")), ("LP iterations", m.get("lp_iterations")),
                                    ("max depth", m.get("max_depth")), ("strong-branch LPs", m.get("strong_branch_lps")),
                                    ("root LP bound", m.get("root_lp_bound")), ("root bound after cuts", m.get("root_bound_after_cuts")),
                                    ("cut rounds", m.get("cut_rounds")), ("GMI / MIR / cover", f"{m.get('cuts_gomory', 0)} / {m.get('cuts_mir', 0)} / {m.get('cuts_cover', 0)}"),
                                    ("bound prunes", m.get("bound_prunes")), ("propagation prunes", m.get("propagation_prunes")),
                                    ("Farkas-certified prunes", m.get("farkas_certified_prunes")), ("integral leaves", m.get("integral_leaves")),
                                    ("heuristic solutions", m.get("heuristic_solutions")), ("incumbent source", m.get("incumbent_source")),
                                    ("safe-bound corrections", m.get("safe_bound_corrections")), ("bound is safe", m.get("bound_is_safe"))])
            prog = m.get("progress_t_incumbent_bound_minform") or []
            if len(prog) >= 2:
                rows.append("  " + c("gap over time", DIM))
                step = max(1, len(prog) // 8)
                for t, inc, bd in prog[::step][-7:] + [prog[-1]]:
                    gap = abs(inc - bd) / max(1.0, abs(inc)) if (inc is not None and bd is not None and abs(inc) < 1e300) else 1.0
                    rows.append("    " + c(pad(fsec(t), 9, "right"), DIM) + "  " + bar(1 - min(1, gap), 30) +
                                c(f"  gap {gap:.2%}", WHITE))
        for key in ("pdhg", "ipm"):
            if tel.get(key):
                items = [(k2, v) for k2, v in tel[key].items() if isinstance(v, (int, float, str, bool))]
                grid(key.upper(), items[:18])
        pre = tel.get("presolve")
        if pre:
            grid("Presolve", [(k2.replace("_", " "), v) for k2, v in pre.items() if v not in (0, None)][:14])
        rt = tel.get("router")
        if rt and not tel.get("model", {}).get("is_mip") and not tel.get("model", {}).get("is_qp"):
            rows.append("---")
            rows.append(section("Router", "predicted engine times for this model (log scale)"))
            pred = rt.get("predicted_seconds", {})
            self.log_bars(rows, pred, highlight=rt.get("selected"),
                          actual=(rt.get("selected"), tel.get("timing", {}).get("engine")))
        box(rows, title=f"Debug · run #{rec['n']} · {rec['model']}")

    @staticmethod
    def log_bars(rows, values, highlight=None, actual=None, width=34):
        vals = {k: v for k, v in values.items() if isinstance(v, (int, float)) and 0 < v < 1e12}
        if not vals:
            return
        missing = [k for k in values if k not in vals]
        lo, hi = math.log10(min(vals.values())) - 0.3, math.log10(max(vals.values()))
        for k, v in sorted(vals.items(), key=lambda kv2: kv2[1]):
            frac = (math.log10(v) - lo) / max(1e-9, hi - lo)
            hl = k == highlight
            rows.append("  " + (c("▸ ", BOLD, YELLOW) if hl else "  ") + c(pad(k, 10), BOLD if hl else "", YELLOW if hl else WHITE)
                        + bar(frac, width, YELLOW if hl else ORANGE) + c(f"  {fsec(v)}", WHITE) +
                        (c(f"   measured {fsec(actual[1])}", DIM) if actual and actual[0] == k and actual[1] else ""))
        for k in missing:
            rows.append("    " + c(pad(k, 10), DIM) + c("not probed / not available for this run", DIM))

    def cmd_analyze(self, args):
        path = self.need_model(args, "/analyze <model>")
        if not path:
            return
        flags = [] if self.settings["gpu"] else ["--no-gpu"]
        rc, lines, _ = run_process([EXE, "analyze", path, "--log", "0"] + flags, f"Analyzing {Path(path).name}")
        try:
            a = json.loads("\n".join(l for l in lines if not l.startswith("read ")))
        except ValueError:
            out(c(f"  ✗ analyze failed (exit {rc})", RED))
            return
        rows = []
        cls_col = {"LP": YELLOW, "MILP": ORANGE, "QP": AMBER}.get(a["class"], WHITE)
        rows.append(kv("Class", c(a["class"], BOLD, cls_col) + c(f"   {a['sense']}imize", DIM)))
        rows.append(kv("Size", c(f"{fmt(a['rows'])} rows · {fmt(a['cols'])} columns · {fmt(a['nnz'])} nonzeros", WHITE) +
                       c(f"   density {a['structure']['density']:.2e}", DIM)))
        if a["integers"]:
            rows.append(kv("Integers", c(f"{fmt(a['integers'])} ({fmt(a['binaries'])} binary)", WHITE)))
        if a["q_nnz"]:
            rows.append(kv("Quadratic", c(f"{fmt(a['q_nnz'])} Q nonzeros", WHITE)))
        rows.append("---")
        rows.append(section("Constraints and bounds"))
        rt_ = a["row_types"]
        tot = max(1, a["rows"])
        for k2, lab in (("equality", "= equality"), ("less_equal", "≤ less-equal"), ("greater_equal", "≥ greater-equal"),
                        ("range", "range"), ("free", "free (objective-like)")):
            if rt_.get(k2):
                rows.append("  " + c(pad(lab, 22), WHITE) + bar(rt_[k2] / tot, 30) + c(f"  {fmt(rt_[k2])}", DIM))
        cb = a["col_bounds"]
        rows.append("  " + c("columns: ", DIM) + "  ".join(c(f"{k2.replace('_', ' ')} ", DIM) + c(fmt(v), WHITE)
                                                          for k2, v in cb.items() if v))
        rows.append("---")
        rows.append(section("Numerics", "orders of magnitude of |values|"))
        num = a["numerics"]

        def mag(lo, hi, label):
            if not hi:
                return
            span = math.log10(hi / lo) if lo else 0
            colr = RED if span >= 10 else (ORANGE if span >= 6 else YELLOW)
            rows.append("  " + c(pad(label, 12), WHITE) + c(pad(f"{lo:.1e} … {hi:.1e}", 22), WHITE) +
                        bar(min(1, span / 12), 24, colr) + c(f"  {span:.1f} decades", DIM))

        mag(num["matrix_min_abs"], num["matrix_max_abs"], "matrix")
        mag(num["cost"]["min_abs"], num["cost"]["max_abs"], "costs")
        mag(num["rhs"]["min_abs"], num["rhs"]["max_abs"], "rhs")
        mag(num["bounds"]["min_abs"], num["bounds"]["max_abs"], "bounds")
        s = a["structure"]
        rows.append("  " + c("max row / col nnz ", DIM) + c(f"{s['max_row_nnz']} / {s['max_col_nnz']}", WHITE) +
                    c("   singleton rows / cols ", DIM) + c(f"{s['singleton_rows']} / {s['singleton_cols']}", WHITE) +
                    c("   empty rows / cols ", DIM) + c(f"{s['empty_rows']} / {s['empty_cols']}", WHITE))
        r = a.get("router", {})
        if a["class"] == "LP":
            rows.append("---")
            rows.append(section("Engine prediction", "router cost model" + (" (calibrated)" if r.get("calibrated_model") else "")
                                + ("" if r.get("gpu_available") else " · no GPU")))
            self.log_bars(rows, r.get("predicted_seconds", {}), highlight=r.get("engine"))
        rows.append("---")
        rows.append(section("Findings"))
        for lvl, text in analysis_findings(a):
            for i, line in enumerate(wrap(text, panel_width() - 8)):
                rows.append("  " + (finding(lvl, line) if i == 0 else "  " + c(line, WHITE)))
        box(rows, title=f"Analysis · {a.get('name', Path(path).stem)}")

    def cmd_verify(self, args):
        rec = self.last or (self.state["runs"][-1] if self.state.get("runs") else None)
        if not rec or not Path(rec["json"]).exists():
            out(c("  no result to verify yet — solve a model first", ORANGE))
            return
        rows = [section("1 · C++ certifier", "interval arithmetic, outward rounding")]
        rc, lines, _ = run_process([EXE, "verify", rec["path"], rec["json"], "--log", "0"], "Re-certifying (C++)")
        for l in lines:
            if l.startswith("read "):
                continue
            t = l.strip()
            v = "PASS" if " PASS " in l else ("FAIL" if " FAIL " in l else ("WARN" if " WARN " in l else None))
            rows.append("  " + (verdict_mark(v) + " " if v else "") + c(t, WHITE if v != "FAIL" else RED))
        rows.append("  " + (c("✓ accepted", BOLD, YELLOW) if rc == 0 else c("✗ rejected", BOLD, RED)))
        rows.append("---")
        exact = bool(args and args[0].lower().startswith("exact"))
        rows.append(section("2 · Exact-rational verifier", "own MPS parser, Python fractions"
                            + (" + exact basis" if exact else "")))
        cmd = [sys.executable, "-m", "pramana.verify", rec["path"], rec["json"]] + (["--exact-basis"] if exact else [])
        rc2, lines2, secs = run_process(cmd, "Verifying in exact arithmetic", cwd=ROOT / "python", show_log=False)
        try:
            rep = json.loads("\n".join(lines2))
            for ch in rep.get("checks", []):
                val = ch.get("value")
                rows.append(f"  {verdict_mark(ch.get('verdict'))} " + pad(c(ch.get("check", ""), WHITE), 32) +
                            c(pad(sci(val) if isinstance(val, (int, float)) else "", 10, "right"), YELLOW) +
                            c(f"  {ch.get('detail', '')}", DIM))
            rows.append("  " + (c("✓ verified", BOLD, YELLOW) if rep.get("verified") else c("✗ NOT verified", BOLD, RED))
                        + c(f"   {fsec(secs)}", DIM))
        except ValueError:
            for l in lines2[-6:]:
                rows.append("  " + c(l, RED))
        if not exact:
            rows.append(c("  tip: /verify exact  also proves the final basis is an exact optimal vertex (small LPs)", DIM))
        box(rows, title=f"Verify · run #{rec['n']} · {rec['model']}")

    def cmd_compare(self, args):
        path = self.need_model(args, "/compare <model>")
        if not path:
            return
        engines = ["dual", "ipm", "pdhg-cpu"] + (["pdhg-gpu"] if self.sysinfo["gpu"] and self.settings["gpu"] else [])
        name = Path(path).name.split(".")[0]
        res = []
        tl = min(self.settings["time"], 120)
        for e in engines:
            jf = STATE_DIR / f"compare_{name}_{e}.json"
            STATE_DIR.mkdir(parents=True, exist_ok=True)
            rc, lines, secs = run_process([EXE, path, "--algo", e, "--time", f"{tl:g}", "--log", "1", "--json", jf],
                                          f"{name}: {e}")
            if rc == 130:
                return
            try:
                j = json.loads(jf.read_text(encoding="utf-8"))
            except Exception:  # noqa: BLE001
                res.append((e, "ERROR", None, secs, False))
                continue
            if j.get("telemetry", {}).get("model", {}).get("is_mip") or j.get("telemetry", {}).get("model", {}).get("is_qp"):
                out(c("  /compare is for LPs (MILP and QP have one engine each); showing a single solve", ORANGE))
                self.cmd_solve(args)
                return
            res.append((e, j.get("status"), j.get("objective"), j.get("seconds"), j.get("certificate", {}).get("accepted")))
        ok = [r for r in res if r[1] == "OPTIMAL" and r[4]]
        best = min((r[3] for r in ok), default=None)
        rows = [c(pad("engine", 12), DIM) + c(pad("status", 18), DIM) + c(pad("objective", 20), DIM) +
                c(pad("time", 10), DIM) + c("relative time (log)", DIM)]
        ts = [r[3] for r in res if r[3]]
        lo, hi = (math.log10(min(ts)) - 0.3, math.log10(max(ts))) if ts else (0, 1)
        for e, st, obj, secs, cert in res:
            win = best is not None and secs == best and st == "OPTIMAL"
            frac = (math.log10(secs) - lo) / max(1e-9, hi - lo) if secs else 0
            rows.append(c(pad(("★ " if win else "  ") + e, 12), BOLD if win else "", YELLOW if win else WHITE) +
                        c(pad(st + (" ✓" if cert else ""), 18), status_color(st)) + c(pad(fmt(obj, 11), 20), WHITE) +
                        c(pad(fsec(secs), 10), WHITE) + bar(frac, 26, YELLOW if win else ORANGE))
        rows.append("---")
        rows.append(c("  every time includes presolve, crossover to an exact vertex and certification; GPU times include "
                      "transfers", DIM))
        box(rows, title=f"Engine comparison · {name} · limit {tl:g}s")

    def cmd_param(self, args):
        if len(args) < 5:
            out(c("  usage: /param <model> <column> <cost|upper|lower|rhs> <from> <to>", ORANGE))
            out(c("  e.g.   /param plan_10x12 BUY_CR03_0 upper 0 120", DIM))
            return
        path = self.need_model(args, "/param ...")
        if not path:
            return
        col, kind, lo, hi = args[1], args[2], args[3], args[4]
        jf = STATE_DIR / "param.json"
        STATE_DIR.mkdir(parents=True, exist_ok=True)
        flag = "--row" if kind == "rhs" else "--col"
        rc, lines, secs = run_process([EXE, "parametric", path, flag, col, "--kind", kind, "--from", lo, "--to", hi,
                                       "--json", jf, "--log", "0"], f"Parametric sweep of {col}")
        if rc not in (0, 2) or not jf.exists():
            for l in lines[-4:]:
                out("  " + c(l, RED))
            return
        p = json.loads(jf.read_text(encoding="utf-8"))
        segs = p.get("segments", [])
        rows = [kv("Parameter", c(f"{col} ({kind})", YELLOW) + c(f"  from {lo} to {hi}", DIM)),
                kv("Result", c(f"{len(segs)} segments · {len(segs) - 1} breakpoints · {fmt(p.get('pivots'))} pivots · "
                               f"{fsec(p.get('seconds'))}", WHITE) +
                   ("  " + c(" all certified ✓ ", BOLD, BLACK, BG_YELLOW) if segs and all(s.get("certified") for s in segs)
                    else "")),
                "---", section("Value curve", "objective vs parameter; ┴ = breakpoint")]
        self.area_chart(rows, segs, float(lo), float(hi))
        rows.append("---")
        rows.append(c(pad("θ from", 14, "right") + pad("θ to", 14, "right") + pad("objective", 16, "right") +
                      pad("slope (dual)", 14, "right") + "   certificate", DIM))
        for s in segs[:14]:
            rows.append(c(pad(f"{s['theta_lo']:.5g}", 14, "right") + pad(f"{s['theta_hi']:.5g}", 14, "right") +
                          pad(f"{s['objective_lo']:.8g}", 16, "right") + pad(f"{s['slope']:.4g}", 14, "right"), WHITE) +
                        "   " + (c("✓ certified", YELLOW) + c(f" gap {sci(s.get('certified_gap'))}", DIM) if s.get("certified")
                                 else c(s.get("status", "?"), status_color(s.get("status")))))
        if len(segs) > 14:
            rows.append(c(f"  … {len(segs) - 14} more segments in {jf.relative_to(ROOT)}", DIM))
        box(rows, title=f"Parametric analysis · {Path(path).stem}")

    @staticmethod
    def area_chart(rows, segs, lo, hi, height=9):
        width = min(70, panel_width() - 20)
        if not segs or hi <= lo:
            return

        def val(th):
            for s in segs:
                if s["theta_lo"] - 1e-12 <= th <= s["theta_hi"] + 1e-12:
                    a, b = s["theta_lo"], s["theta_hi"]
                    return s["objective_lo"] + (0 if b == a else (th - a) / (b - a) * (s["objective_hi"] - s["objective_lo"]))
            return None

        xs = [lo + (hi - lo) * i / (width - 1) for i in range(width)]
        ys = [val(x) for x in xs]
        good = [y for y in ys if y is not None and math.isfinite(y)]
        if not good:
            return
        ymin, ymax = min(good), max(good)
        span = (ymax - ymin) or 1.0
        levels = " ▁▂▃▄▅▆▇█"
        grid = []
        for r in range(height, 0, -1):
            line = ""
            for y in ys:
                if y is None:
                    line += " "
                    continue
                h = (y - ymin) / span * (height - 1) + 1
                cell = h - (r - 1)
                line += "█" if cell >= 1 else (levels[max(0, int(cell * 8))] if cell > 0 else " ")
            grid.append(line)
        for i, line in enumerate(grid):
            label = f"{ymax:.6g}" if i == 0 else (f"{ymin:.6g}" if i == height - 1 else "")
            rows.append(c(pad(label, 12, "right"), DIM) + c(" ┤", ORANGE) + c(line, YELLOW))
        axis = ["─"] * width
        for s in segs[1:]:
            k = int(round((s["theta_lo"] - lo) / (hi - lo) * (width - 1)))
            if 0 <= k < width:
                axis[k] = "┴"
        rows.append(" " * 12 + c(" └", ORANGE) + "".join(c(a, BOLD, YELLOW) if a == "┴" else c(a, ORANGE) for a in axis))
        rows.append(" " * 14 + c(pad(f"{lo:.4g}", width // 2), DIM) + c(pad(f"{hi:.4g}", width - width // 2, "right"), DIM))

    def cmd_bench(self, args):
        sets = [("netlib", "Netlib LP (91 feasible)"), ("netlib-infeas", "Netlib infeasible (29)"),
                ("adversarial", "Adversarial suite (17)"), ("miplib3", "MIPLIB 3 subset (42)"),
                ("maros", "Maros–Mészáros QP (124)"), ("gen", "Industrial models (18)")]
        rows = [c(pad("benchmark", 26) + pad("solver", 17) + pad("solved", 26) + pad("SGM time", 11) + "wrong", DIM)]
        found = False
        # Ground truth for the adversarial suite (known by construction): "Score ...: PRAMANA [c, w, n], HiGHS [c, w, n]"
        truth = {}
        bm = ROOT / "docs" / "BENCHMARKS.md"
        if bm.exists():
            m2 = re.search(r"PRAMANA \[(\d+), (\d+), (\d+)\], HiGHS \[(\d+), (\d+), (\d+)\]", bm.read_text(encoding="utf-8"))
            if m2:
                truth = {"PRAMANA": m2.group(2), "HiGHS": m2.group(5)}
        for key, title in sets:
            f = ROOT / "results" / key / "summary.md"
            if not f.exists():
                continue
            found = True
            first = True
            for line in f.read_text(encoding="utf-8").splitlines():
                m = re.match(r"\|\s*(\w[\w-]*)\s*\|\s*(\d+)/(\d+)\s*\|\s*([\d-]+)\s*\|\s*([\d.]+)\s*\|\s*([\d-]+)\s*\|", line)
                if not m:
                    continue
                engine, a, b, _, sgm, wrong = m.groups()
                solver = "HiGHS" if engine == "highs" else "PRAMANA"
                label = solver if engine in ("highs", "auto") else f"PRAMANA {engine}"
                if key == "adversarial" and truth:
                    wrong = truth[solver]  # vs the known true answer, not vs HiGHS
                frac = int(a) / max(1, int(b))
                col = YELLOW if solver == "PRAMANA" else WHITE
                rows.append(c(pad(title if first else "", 26), WHITE) + c(pad(label, 17), BOLD if solver == "PRAMANA" else "", col) +
                            bar(frac, 14, YELLOW if solver == "PRAMANA" else WHITE) + c(pad(f" {a}/{b}", 12), WHITE) +
                            c(pad(f"{float(sgm):.2f}s", 11), DIM) + (c(wrong, BOLD, RED) if wrong not in ("0", "-") else c(wrong, DIM)))
                first = False
        if not found:
            out(c("  no benchmark results yet — run  bash bench/run_all.sh", ORANGE))
            return
        rows.append("---")
        rows.append(c("  'wrong': disagreements with HiGHS; for the adversarial suite, wrong answers vs the true answer known by "
                      "construction (HiGHS errs on the badly scaled models)", DIM))
        rows.append(c("  full tables: docs/BENCHMARKS.md · results/<set>/summary.md", DIM))
        box(rows, title="Benchmarks · same machine, same limits")

    def cmd_gpu(self, args):
        if args and args[0].startswith("cal"):
            rc, lines, _ = run_process([EXE, "calibrate", "--out", STATE_DIR / "calibrate.json"], "Measuring PDHG iteration cost")
            cal = STATE_DIR / "calibrate.json"
        else:
            cal = STATE_DIR / "calibrate.json"
            if not cal.exists():
                cal = ROOT / "results" / "gpu" / "calibrate.json"
        rc, lines, _ = run_process([EXE, "info"], "Querying device", show_log=False)
        info = {}
        try:
            info = json.loads("\n".join(lines[lines.index("{"):]) if "{" in lines else "{}")
        except ValueError:
            pass
        rows = []
        if info.get("available"):
            rows.append(kv("Device", c(info.get("name"), BOLD, YELLOW) + c(
                f"   compute {info.get('compute_capability')} · {info.get('multiprocessors')} SMs · "
                f"{info.get('memory_bytes', 0) / 1e9:.1f} GB · peak {info.get('peak_bandwidth_GBs', 0):.0f} GB/s", WHITE)))
            rows.append(kv("Kernels", c("own CUDA C, compiled at runtime by NVRTC", WHITE) +
                           c("  (PTX cached)" if info.get("ptx_from_cache") else "", DIM)))
        else:
            rows.append(kv("Device", c("no GPU available", ORANGE) + c(f"  {info.get('reason', '')}", DIM)))
        if cal.exists():
            pts = json.loads(cal.read_text(encoding="utf-8")).get("pdhg_iteration_timing", [])
            rows.append("---")
            rows.append(section("PDHG iteration cost", f"CPU all threads vs GPU · {cal.relative_to(ROOT)}"))
            rows.append(c(pad("nonzeros", 14, "right") + pad("CPU ms", 11, "right") + pad("GPU ms", 11, "right") +
                          "   GPU speed-up (1× = break-even)", DIM))
            for pt in pts:
                cpu = pt["cpu"]["seconds_per_iteration"] * 1e3
                g = pt.get("gpu", {}).get("seconds_per_iteration")
                sp = cpu / (g * 1e3) if g else None
                rows.append(c(pad(f"{pt['nnz']:,}", 14, "right") + pad(f"{cpu:.3f}", 11, "right") +
                              pad(f"{g * 1e3:.3f}" if g else "—", 11, "right"), WHITE) + "   " +
                            (bar(min(1, sp / 6), 24, YELLOW if sp >= 1 else RED) + c(f" {sp:.2f}×", YELLOW if sp >= 1 else RED)
                             if sp else ""))
            rows.append(c("  the GPU pays off only above the break-even size; the router uses this to pick engines", DIM))
        rows.append(c("  /gpu calibrate  re-measures on this machine", DIM))
        box(rows, title="GPU")

    def cmd_models(self, args):
        flt = args[0].lower() if args else ""
        groups = {}
        for n, p in self.models.items():
            if flt in n or flt in Path(p).parent.name:
                groups.setdefault(Path(p).parent.name, []).append(n)
        if not groups:
            out(c(f"  no models match '{flt}'", ORANGE))
            return
        rows = []
        w = panel_width() - 6
        for g, names in groups.items():
            rows.append(section(g, f"{len(names)} models"))
            line = "  "
            for n in names[:60]:
                if vlen(line) + len(n) + 2 > w:
                    rows.append(c(line, WHITE))
                    line = "  "
                line += n + "  "
            rows.append(c(line, WHITE) + (c(f"… +{len(names) - 60}", DIM) if len(names) > 60 else ""))
        box(rows, title="Models" + (f" · '{flt}'" if flt else ""))

    def cmd_set(self, args):
        if len(args) < 2:
            out(c("  settings: ", DIM) + "  ".join(c(f"{k}=", DIM) + c(str(v), YELLOW) for k, v in self.settings.items()))
            out(c("  usage: /set time 60 · /set algo auto|dual|ipm|pdhg-cpu|pdhg-gpu|race · /set threads 8 · /set gpu off", DIM))
            return
        k, v = args[0].lower(), args[1].lower()
        try:
            if k == "time":
                self.settings["time"] = float(v)
            elif k == "algo":
                assert v in ("auto", "dual", "ipm", "pdhg", "pdhg-cpu", "pdhg-gpu", "race")
                self.settings["algo"] = v
            elif k == "threads":
                self.settings["threads"] = int(v)
            elif k == "gpu":
                self.settings["gpu"] = v in ("on", "1", "true", "yes")
            else:
                raise ValueError
            out(c("  ✓ ", YELLOW) + c(f"{k} = {v}", WHITE))
        except (ValueError, AssertionError):
            out(c(f"  ✗ invalid setting {k} {v}", RED))

    def cmd_check(self, args):
        if os.name != "nt":
            out(c("  /check runs check.bat (Windows). On Linux/macOS: build/pramana_tests && pytest tests/python", ORANGE))
            return
        env = dict(os.environ, PYTHON=sys.executable)
        t0 = time.perf_counter()
        proc = subprocess.Popen(["cmd", "/c", str(ROOT / "check.bat")], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                cwd=str(ROOT), text=True, encoding="utf-8", errors="replace", env=env)
        rows_pass = rows_fail = 0
        try:
            for line in proc.stdout:
                line = line.rstrip()
                if "  PASS  " in line:
                    rows_pass += 1
                    out("  " + c("✓", BOLD, YELLOW) + " " + c(line.split("PASS", 1)[1].strip(), WHITE))
                elif "  FAIL  " in line:
                    rows_fail += 1
                    out("  " + c("✗", BOLD, RED) + " " + c(line.split("FAIL", 1)[1].strip(), RED))
                elif line.startswith("==="):
                    out(c("  " + line.strip("= "), BOLD, ORANGE))
                else:
                    sys.stdout.write("\r\x1b[2K" + truncate(c("  " + line.strip(), DIM), term_width() - 1))
                    sys.stdout.flush()
            proc.wait()
        except KeyboardInterrupt:
            proc.kill()
            out(c("\n  ✗ interrupted", RED))
            return
        sys.stdout.write("\r\x1b[2K")
        ok = proc.returncode == 0 and rows_fail == 0
        out()
        out("  " + (c(f" ALL {rows_pass} CHECKS PASSED ", BOLD, BLACK, BG_YELLOW) if ok else
                    c(f" {rows_fail} FAILED · {rows_pass} passed ", BOLD, WHITE, BG_RED)) +
            c(f"  {fsec(time.perf_counter() - t0)}", DIM))

    def cmd_history(self, args):
        runs = self.state.get("runs", [])[-20:][::-1]
        if not runs:
            out(c("  no runs yet", DIM))
            return
        rows = [c(pad("#", 5) + pad("model", 18) + pad("status", 18) + pad("objective", 20) + pad("time", 10) + "when", DIM)]
        for r in runs:
            st = r.get("status", "?")
            rows.append(c(pad(str(r.get("n")), 5), DIM) + c(pad(r.get("model", ""), 18), WHITE) +
                        c(pad(st, 18), status_color(st)) + c(pad(fmt(r.get("objective"), 10), 20), WHITE) +
                        c(pad(fsec(r.get("seconds")), 10), WHITE) + c(r.get("when", ""), DIM))
        rows.append(c("  /debug <#> shows the diagnostics of any run", DIM))
        box(rows, title="History")

    def cmd_help(self, args):
        rows = [section("Commands")]
        for n, a, d, _ in COMMANDS:
            rows.append("  " + c(pad(n, 10), BOLD, YELLOW) + c(pad(a, 34), WHITE) + c(d, DIM))
        rows.append("---")
        rows.append(section("Keys"))
        for k2, d in (("Enter", "run · a bare model name (e.g. afiro, p0201, pilot87) solves it"),
                      ("Tab", "complete command or model name"), ("↑ ↓", "history, or move through suggestions"),
                      ("Esc", "clear the line"), ("Ctrl+C", "cancel a running solve · twice on an empty line exits"),
                      ("Ctrl+L", "clear the screen")):
            rows.append("  " + c(pad(k2, 10), BOLD, ORANGE) + c(d, WHITE))
        rows.append("---")
        rows.append(c("  Solver flags pass straight through, e.g. ", DIM) + c("/solve pilot87 --algo ipm --no-presolve", YELLOW))
        box(rows, title="Help")

    # ----- main loop
    def dispatch(self, line):
        if not line:
            return True
        hist = self.state.setdefault("history", [])
        if not hist or hist[-1] != line:
            hist.append(line)
        out(BG_INPUT + c("> ", BOLD, ORANGE) + BG_INPUT + c(line + " ", WHITE) + RESET)
        try:
            toks = [t.strip('"') for t in shlex.split(line, posix=False)]
        except ValueError:
            toks = line.split()
        cmd, args = toks[0].lower(), toks[1:]
        table = {"/solve": self.cmd_solve, "/analyze": self.cmd_analyze, "/analyse": self.cmd_analyze,
                 "/debug": self.cmd_debug, "/verify": self.cmd_verify, "/compare": self.cmd_compare,
                 "/param": self.cmd_param, "/bench": self.cmd_bench, "/gpu": self.cmd_gpu, "/models": self.cmd_models,
                 "/set": self.cmd_set, "/check": self.cmd_check, "/history": self.cmd_history, "/help": self.cmd_help}
        try:
            if cmd in ("/exit", "/quit", "/q", "exit", "quit"):
                return False
            if cmd == "/clear":
                sys.stdout.write("\x1b[2J\x1b[3J\x1b[H")
                self.welcome()
            elif cmd in table:
                table[cmd](args)
            elif not cmd.startswith("/") and self.resolve(toks[0]):
                self.cmd_solve(toks)
            else:
                out(c(f"  ✗ unknown command '{toks[0]}'", RED) + c("   type / to see commands, or a model name", DIM))
        except KeyboardInterrupt:
            out(c("  ✗ interrupted", RED))
        except Exception as e:  # noqa: BLE001
            out(c(f"  ✗ {type(e).__name__}: {e}", RED))
        save_state(self.state)
        out()
        self.tip += 1
        return True

    def run(self):
        enable_vt()
        if not EXE.exists():
            out(c(f"✗ solver binary not found: {EXE}", RED))
            out(c("  build it first:  build.bat   (or start.bat, which builds and launches this console)", WHITE))
            return 1
        sys.stdout.write("\x1b[2J\x1b[3J\x1b[H")
        sys.stdout.write(c("  starting PRAMANA …", DIM))
        sys.stdout.flush()
        self.probe()
        sys.stdout.write("\r\x1b[2K")
        self.welcome()
        if not sys.stdin.isatty():
            for line in sys.stdin.read().splitlines():  # scripted use: commands on stdin
                if not self.dispatch(line.strip().lstrip("﻿")):
                    break
            return 0
        self.keys = Keys()
        try:
            while True:
                line = self.read_line()
                if not self.dispatch(line):
                    break
        except KeyboardInterrupt:
            self.clear_area()
        save_state(self.state)
        out(c("  goodbye — every answer was certified.", DIM))
        return 0


def main():
    if "--help" in sys.argv or "-h" in sys.argv:
        print(__doc__)
        return 0
    return App().run()


if __name__ == "__main__":
    sys.exit(main())
