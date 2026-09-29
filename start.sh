#!/usr/bin/env bash
# PRAMANA - one command to start everything (Linux / macOS):
# builds the solver (incremental), prepares the benchmark data if missing, opens the terminal console.
set -e
cd "$(dirname "$0")"
mkdir -p build
PY=${PYTHON:-python3}
command -v "$PY" >/dev/null || { echo "Python 3.8+ not found (set PYTHON=...)"; exit 1; }
printf '\n \033[38;2;255;140;26mPRAMANA\033[0m \033[38;2;135;135;135m- certified LP / MILP / QP solver\033[0m\n\n'
if cmake -S . -B build -DCMAKE_BUILD_TYPE=Release >build/start_build.log 2>&1 && cmake --build build -j >>build/start_build.log 2>&1; then
  echo " ✓ build       build/pramana"
elif [ -x build/pramana ]; then
  echo " ! build       could not rebuild - using the existing build/pramana"
else
  echo " ✗ build failed - see build/start_build.log"; exit 1
fi
[ -f data/netlib/afiro.mps.gz ] || "$PY" bench/fetch_data.py netlib netlib-infeas miplib3 maros || true
[ -f data/gen/plan_10x12.mps ] || "$PY" gen/refinery.py suite --out data/gen >/dev/null
[ -f data/adversarial/scaled_afiro.mps ] || "$PY" gen/adversarial.py --out data/adversarial >/dev/null
echo " ✓ data        benchmark and industrial models ready"
PYTHONPATH="$PWD/python${PYTHONPATH:+:$PYTHONPATH}" exec "$PY" -m pramana.tui "$@"
