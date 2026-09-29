"""Fit and VALIDATE the CPU/GPU engine router from measured benchmark runs.

Data: every LP instance solved by every engine (results/netlib*, results/gen).
Model per engine e:  log t_e = w_e . phi(features)   (ridge least squares), with
phi identical to src/router/router.cpp. Unsolved runs are censored at 4x the limit.

Validation: 5-fold cross-validation grouped by instance family (the router never
sees the held-out instances), reporting
  * oracle agreement: how often the predicted-fastest engine is the truly fastest,
  * regret: time(chosen) / time(oracle), geometric mean and worst case,
compared with the fixed policies "always dual simplex", "always GPU PDHG", "always IPM".
Writes build/router_model.json (loaded by `pramana --router-model` or
PRAMANA_ROUTER_MODEL) and results/router/summary.md.
"""
from __future__ import annotations

import csv
import json
import math
import random
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
ENGINES = ["dual", "ipm", "pdhg-cpu", "pdhg-gpu"]
PROVEN = {"OPTIMAL", "INFEASIBLE", "UNBOUNDED"}


def phi(f):
    lg = lambda v: math.log(max(float(v or 0), 1.0))  # noqa: E731
    return [1.0, lg(f["rows"] + f["cols"]), lg(f["nnz"]), lg(f["avg_row_nnz"]), lg(f["max_row_nnz"]),
            lg(f["coefficient_range"]), f["equality_row_fraction"], f["free_column_fraction"]]


def load():
    data = {}  # instance -> {"family", "phi", times{engine: t}}
    for d in (ROOT / "results").glob("*"):
        if not (d / "runs.csv").exists() or d.name.startswith("miplib") or d.name.startswith("maros") \
                or d.name.startswith("adversarial") or "infeas" in d.name or "smoke" in d.name:
            continue
        fam = "gen" if d.name.startswith("gen") else "netlib"
        with open(d / "runs.csv") as fh:
            for r in csv.DictReader(fh):
                e = r["solver"]
                if e not in ENGINES:
                    continue
                js = d / "json" / f"{r['instance']}.{e}.json"
                if not js.exists():
                    continue
                j = json.loads(js.read_text())
                feats = j.get("telemetry", {}).get("router", {}).get("features")
                if not feats or feats.get("is_mip") or feats.get("is_qp"):
                    continue
                inst = data.setdefault(r["instance"], {"family": fam, "phi": phi(feats), "times": {}, "limit": {}})
                t = float(r["time"]) if r["time"] else None
                lim = 120.0 if fam == "gen" else (30.0 if e != "dual" else 300.0)
                ok = r["status"] in PROVEN and t is not None
                inst["times"][e] = t if ok else 4 * lim
    return data


def fit(X, y, lam=1e-2):
    X = np.asarray(X)
    y = np.asarray(y)
    A = X.T @ X + lam * np.eye(X.shape[1])
    return np.linalg.solve(A, X.T @ y)


def main():
    data = load()
    names = sorted(data)
    engines = [e for e in ENGINES if all(e in data[n]["times"] for n in names)]
    out = ROOT / "results" / "router"
    out.mkdir(parents=True, exist_ok=True)
    if len(names) < 10:
        print("not enough data yet:", len(names))
        return
    rng = random.Random(1)
    folds = {}
    for fam in ("netlib", "gen"):
        members = [n for n in names if data[n]["family"] == fam]
        rng.shuffle(members)
        for i, n in enumerate(members):
            folds[n] = i % 5
    choices, oracle, fixed = {}, {}, {e: {} for e in engines}
    for f in range(5):
        train = [n for n in names if folds[n] != f]
        test = [n for n in names if folds[n] == f]
        W = {e: fit([data[n]["phi"] for n in train], [math.log(data[n]["times"][e] + 1e-3) for n in train])
             for e in engines}
        for n in test:
            pred = {e: float(np.dot(W[e], data[n]["phi"])) for e in engines}
            choices[n] = min(pred, key=pred.get)
    for n in names:
        t = data[n]["times"]
        oracle[n] = min(t.values())
        for e in engines:
            fixed[e][n] = t[e]

    def regret(sel):
        r = [(data[n]["times"][sel[n]] + 1e-3) / (oracle[n] + 1e-3) for n in names]
        return math.exp(sum(math.log(x) for x in r) / len(r)), max(r), sum(1 for x in r if x <= 1.0001) / len(r)

    lines = ["# Router validation (5-fold cross-validation, held-out instances)", "",
             f"Instances: {len(names)} LPs ({sum(data[n]['family'] == 'netlib' for n in names)} Netlib, "
             f"{sum(data[n]['family'] == 'gen' for n in names)} refinery planning). Engines: {', '.join(engines)}.", "",
             "| policy | geo-mean regret (time / oracle time) | worst regret | picks the fastest engine |", "|---|---|---|---|"]
    g, w, a = regret(choices)
    lines.append(f"| **PRAMANA router (held-out)** | {g:.3f} | {w:.1f} | {100 * a:.0f}% |")
    for e in engines:
        g, w, a = regret({n: e for n in names})
        lines.append(f"| always {e} | {g:.3f} | {w:.1f} | {100 * a:.0f}% |")
    counts = {e: sum(1 for n in names if min(data[n]["times"], key=data[n]["times"].get) == e) for e in engines}
    lines += ["", "Oracle (fastest certified engine) counts: " + ", ".join(f"{e}: {c}" for e, c in counts.items()), ""]
    # Final model on all data.
    Wall = {e: fit([data[n]["phi"] for n in names], [math.log(data[n]["times"][e] + 1e-3) for n in names]).tolist()
            for e in engines}
    for e in ENGINES:
        Wall.setdefault(e, [50.0] + [0.0] * 7)  # engine without data: never chosen
    model = {"weights": Wall, "features": ["1", "log(m+n)", "log(nnz)", "log(nnz/row)", "log(max row nnz)",
                                           "log(coef range)", "eq fraction", "free fraction"],
             "instances": len(names)}
    (ROOT / "build" / "router_model.json").write_text(json.dumps(model, indent=1))
    (out / "router_model.json").write_text(json.dumps(model, indent=1))
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
