"""Repeated-solve experiment (refinery crude valuation).

Question a planner asks: "what is crude CR03 worth / at which price does it enter
the slate?" Answered on refinery planning LPs by four strategies (pramana family):
  * exact certified parametric sweep (one run, every breakpoint, every segment certified)
  * K cold re-solves, a K-case warm-started simplex chain, batched PDHG (CPU / GPU)
Reports time, accuracy vs the exact value function, and breakpoints missed by sampling.

Writes results/family/*.json and a markdown table results/family/summary.md
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BUILD_EXE = ROOT / "build" / ("pramana.exe" if os.name == "nt" else "pramana")
OUT = ROOT / "results" / "family"

CASES = [
    # model, column, kind, from, to, K, description
    ("plan_10x12", "BUY_CR03_0", "cost", -95.0, -35.0, 64, "crude CR03 price (period 0), $/bbl, negated (max model)"),
    ("plan_10x12", "BUY_CR03_0", "upper", 0.0, 120.0, 64, "crude CR03 availability (period 0), kbbl"),
    ("plan_20x12", "BUY_CR07_3", "cost", -100.0, -30.0, 64, "crude CR07 price (period 3)"),
    ("plan_20x52", "BUY_CR05_10", "upper", 0.0, 200.0, 32, "crude CR05 availability (week 10)"),
]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    exe = ROOT / "results" / "bin" / ("pramana_family" + BUILD_EXE.suffix)
    exe.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(BUILD_EXE, exe)
    lines = ["# Repeated-solve (crude valuation) experiment", "",
             "Exact = certified parametric sweep. Errors are max relative errors vs the exact value function.", ""]
    for model, col, kind, lo, hi, k, desc in CASES:
        path = ROOT / "data" / "gen" / f"{model}.mps"
        if not path.exists():
            continue
        jf = OUT / f"{model}_{col}_{kind}.json"
        pf = OUT / f"{model}_{col}_{kind}_parametric.json"
        subprocess.run([str(exe), "family", str(path), "--col", col, "--kind", kind, "--from", str(lo), "--to", str(hi),
                        "--cases", str(k), "--json", str(jf), "--log", "0"], capture_output=True, text=True, timeout=7200)
        subprocess.run([str(exe), "parametric", str(path), "--col", col, "--kind", kind, "--from", str(lo), "--to", str(hi),
                        "--json", str(pf), "--log", "0"], capture_output=True, text=True, timeout=7200)
        if not jf.exists():
            continue
        j = json.loads(jf.read_text())
        p = json.loads(pf.read_text()) if pf.exists() else {}
        lines += [f"## {model}: {desc}", "",
                  f"- parameter `{col}` ({kind}) in [{lo}, {hi}], K = {k} sampled cases",
                  f"- exact breakpoints: {j['breakpoints']} ({j['breakpoints_missed_by_sampling']} segments contain no sample"
                  " -> invisible to any sampling strategy)",
                  f"- breakpoints: {', '.join(f'{b:.6g}' for b in p.get('breakpoints', [])[:12])}"
                  + (" ..." if len(p.get("breakpoints", [])) > 12 else ""), "",
                  "| strategy | seconds | solved | certified | max rel. error | note |", "|---|---|---|---|---|---|"]
        for s in j["strategies"]:
            lines.append(f"| {s['strategy']} | {s['seconds']:.3f} | {s['solved']} | {s['certified']} | "
                         f"{s['max_relative_error']:.2e} | {s['note']} |")
        lines.append("")
        print("\n".join(lines[-12:]), flush=True)
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
