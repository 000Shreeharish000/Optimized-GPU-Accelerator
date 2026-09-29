#!/usr/bin/env bash
# One-command, reproducible benchmark campaign. Produces results/*/summary.md,
# results/*/runs.csv, docs/figures/*.png and docs/BENCHMARKS.md.
# Runs sequentially so timings never interfere with each other.
set -e
cd "$(dirname "$0")/.."
PY=${PYTHON:-python}
T_LP=${T_LP:-300}
T_MIP=${T_MIP:-60}

$PY bench/fetch_data.py netlib netlib-infeas miplib3 maros
$PY gen/refinery.py suite --out data/gen
$PY gen/adversarial.py --out data/adversarial

$PY bench/run_bench.py --set netlib        --engines dual --reference --time $T_LP
$PY bench/run_bench.py --set netlib-infeas --engines dual --reference --time $T_LP
$PY bench/run_bench.py --set adversarial   --engines auto --reference --time 60
$PY bench/run_bench.py --set maros         --engines auto --time 60
$PY bench/run_bench.py --set gen           --engines dual,ipm,pdhg-cpu,pdhg-gpu --reference --time 120
$PY bench/run_bench.py --set miplib3       --engines auto --reference --time $T_MIP
$PY bench/run_bench.py --set miplib3       --engines auto --time $T_MIP --extra "--no-cuts" --tag nocuts
$PY bench/run_bench.py --set netlib        --engines ipm,pdhg-cpu,pdhg-gpu --time 30 --tag engines
$PY bench/gpu_crossover.py
$PY bench/family_experiment.py
$PY bench/fit_router.py
$PY bench/report.py
