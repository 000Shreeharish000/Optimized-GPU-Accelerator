"""CPU vs GPU experiment ("GPU where it provides MEASURABLE benefit").

1. Per-iteration PDHG time on CPU (16 threads) vs GPU for synthetic LPs of growing
   size (`pramana calibrate`): locates the nnz where the GPU starts to win per iteration.
2. End-to-end, certified time of every engine on the refinery planning LP family
   (6x4 ... 60x156 crudes x periods), including GPU transfers and crossover.
3. Accuracy trade-off: raw PDHG (no crossover) at 1e-4 / 1e-6 relative KKT.

Writes results/gpu/*.json|csv and docs/figures/gpu_*.png
"""
from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BUILD_EXE = ROOT / "build" / ("pramana.exe" if os.name == "nt" else "pramana")
OUT = ROOT / "results" / "gpu"
FIG = ROOT / "docs" / "figures"


def run(exe, args, timeout):
    return subprocess.run([str(exe)] + args, capture_output=True, text=True, timeout=timeout)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    exe = ROOT / "results" / "bin" / ("pramana_gpu" + BUILD_EXE.suffix)
    exe.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(BUILD_EXE, exe)
    # 1. per-iteration timing
    p = run(exe, ["calibrate", "--out", str(OUT / "calibrate.json")], 3600)
    cal = json.loads((OUT / "calibrate.json").read_text())
    # 2. + 3. planning family
    plans = sorted((ROOT / "data" / "gen").glob("plan_*x*.mps"), key=lambda p: p.stat().st_size)
    plans = [p for p in plans if "mip" not in p.name]
    rows = []
    configs = [("dual", []), ("ipm", []), ("pdhg-cpu", []), ("pdhg-gpu", []),
               ("pdhg-cpu", ["--no-crossover", "--pdhg-tol", "1e-4"]), ("pdhg-gpu", ["--no-crossover", "--pdhg-tol", "1e-4"]),
               ("pdhg-gpu", ["--no-crossover", "--pdhg-tol", "1e-6"])]
    for path in plans:
        for eng, extra in configs:
            tag = eng + ("" if not extra else "_raw" + extra[-1])
            jf = OUT / f"{path.stem}.{tag}.json"
            try:
                run(exe, [str(path), "--algo", eng, "--time", "600", "--log", "0", "--json", str(jf)] + extra, 1300)
                j = json.loads(jf.read_text())
            except Exception as e:  # noqa: BLE001
                rows.append({"model": path.stem, "config": tag, "status": f"ERROR {e}"})
                continue
            tel = j.get("telemetry", {})
            pd = tel.get("pdhg", {})
            rows.append({
                "model": path.stem, "config": tag, "status": j.get("status"), "objective": j.get("objective"),
                "seconds": j.get("seconds"), "nnz": tel.get("model", {}).get("nnz"),
                "rows": tel.get("model", {}).get("rows"), "cols": tel.get("model", {}).get("cols"),
                "pdhg_iterations": pd.get("iterations"), "pdhg_kernel_s": pd.get("kernel_seconds"),
                "pdhg_transfer_s": pd.get("transfer_seconds"), "pdhg_setup_s": pd.get("setup_seconds"),
                "pdhg_rel_kkt": pd.get("relative_kkt"),
                "crossover_s": tel.get("crossover", {}).get("seconds"),
                "cert_gap": j.get("certificate", {}).get("certified_relative_gap"),
            })
            print(rows[-1], flush=True)
    with open(OUT / "planning_engines.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}))
        w.writeheader()
        w.writerows(rows)
    # Figures
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        pts = cal["pdhg_iteration_timing"]
        nnz = [p["nnz"] for p in pts]
        fig, ax = plt.subplots(figsize=(6.4, 4))
        ax.loglog(nnz, [p["cpu"]["seconds_per_iteration"] for p in pts], "o-", label="PDHG iteration, CPU (all threads)")
        if all("gpu" in p for p in pts):
            ax.loglog(nnz, [p["gpu"]["seconds_per_iteration"] for p in pts], "s-", label="PDHG iteration, GPU (own kernels)")
        ax.set_xlabel("nonzeros of A")
        ax.set_ylabel("seconds / iteration")
        ax.set_title("Per-iteration cost: where the GPU starts to win")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(FIG / "gpu_iteration_crossover.png", dpi=130)
        fig, ax = plt.subplots(figsize=(6.4, 4))
        for tag, mk in [("dual", "o-"), ("ipm", "^-"), ("pdhg-cpu", "s--"), ("pdhg-gpu", "D-")]:
            xs = [r["nnz"] for r in rows if r["config"] == tag and r.get("status") == "OPTIMAL"]
            ys = [r["seconds"] for r in rows if r["config"] == tag and r.get("status") == "OPTIMAL"]
            if xs:
                ax.loglog(xs, ys, mk, label=f"{tag} (certified optimal)")
        ax.set_xlabel("nonzeros (refinery planning LP family)")
        ax.set_ylabel("end-to-end seconds incl. transfers + crossover + certificate")
        ax.set_title("Engine time vs size: refinery planning LPs")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(FIG / "gpu_planning_engines.png", dpi=130)
    except Exception as e:  # noqa: BLE001
        print("plot failed:", e)


if __name__ == "__main__":
    main()
