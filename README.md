# PRAMANA: a certified, from-scratch LP / MILP / QP solver with measured GPU acceleration

**Smart India Hackathon 2026 · Problem statement SIH26119 · Mangalore Refinery and Petrochemicals Ltd (MRPL)**
*Indigenous GPU-Accelerated Optimization Solver (sovereign alternative to Xpress / CPLEX)*

> *pramāṇa* (Sanskrit): "valid means of knowledge, proof".

PRAMANA is an optimization solver core written from mathematical foundations in C++20 (~11.6k lines),
with **no optimization or linear-algebra solver library linked or included**. It has its own:

- MPS/QPS reader
- sparse LU
- bounded dual and primal simplex
- sparse LDLᵀ with minimum-degree ordering
- interior-point method
- restarted PDHG with its own CUDA kernels
- presolve and postsolve
- branch-and-cut

Its defining property is that **it never asks to be trusted**:

| | |
|---|---|
| **Every answer is a proof** | OPTIMAL carries a rigorous, outward-rounded dual bound (Neumaier–Shcherbina). INFEASIBLE carries a verified Farkas certificate and UNBOUNDED a verified primal ray. A non-convex QP is rejected with a vector v where vᵀQv < 0. If a proof fails, the status is downgraded and never reported as proven. |
| **Two independent checkers** | A C++ interval-arithmetic certifier, and a Python verifier with its *own* MPS parser that works in exact rational arithmetic (including exact basis/vertex verification). |
| **GPU only where it measurably pays** | Own CUDA kernels compiled at runtime (no CUDA toolkit needed), batched PDHG, a measured CPU/GPU crossover, and an engine router validated on held-out instances. The same binary runs on CPU-only machines. |
| **Built for refinery planning** | Certified parametric analysis gives the exact crude-valuation curve in one run: every breakeven price, with every segment certified. |

---

## Contents

1. [Results at a glance](#1-results-at-a-glance)
2. [Requirements](#2-requirements)
3. [Build and run in 5 minutes](#3-build-and-run-in-5-minutes)
4. [How to check this submission](#4-how-to-check-this-submission)
5. [Using the solver](#5-using-the-solver)
6. [Reproducing every benchmark](#6-reproducing-every-benchmark)
7. [Problem statement coverage](#7-problem-statement-coverage)
8. [Repository layout](#8-repository-layout)
9. [What we claim, and what we do not](#9-what-we-claim-and-what-we-do-not)

---

## 1. Results at a glance

Machine: AMD Ryzen laptop (16 threads), NVIDIA RTX 4050 Laptop GPU (`results/machine.json`). Reference
solver: HiGHS (via scipy) on the same machine with the same limits. Full per-instance tables, including
every failure, are in [docs/BENCHMARKS.md](docs/BENCHMARKS.md) and `results/<set>/summary.md`.

| benchmark | PRAMANA | HiGHS |
|---|---|---|
| **Netlib LP**, all 91 feasible | **91/91 optimal, all certified**, 0 disagreements | 91/91 |
| **Netlib infeasible**, all 29 | **29/29 infeasible, all with verified Farkas proofs** | 28/29 |
| **Badly scaled Netlib** (rows/cols ×10^±6, true optimum known) | **7/8 true optimum + certificate, 0 wrong** (1 honestly unproven) | 0/8 correct: 4 wrong optima, 3 false INFEASIBLE, 1 false UNBOUNDED |
| **Adversarial suite** (17: Klee–Minty, degenerate, near-singular, thin-infeasible, unbounded, scaled) | 16/17 proven + certified, **0 wrong** | 17/17 answered, **8 wrong** |
| **MIPLIB 3** pre-registered subset (42, 60 s) | 29/42 optimal, all certified, 0 disagreements | 36/42 |
| **Maros–Mészáros** convex QP (124, 60 s) | 105/124 optimal + certified, 0 wrong | (no QP in scipy) |
| Cutting planes on/off (unit commitment MILP) | 0.09 s (1 node) vs 29.9 s (20,005 nodes) | |
| PDHG iteration, GPU vs 16-thread CPU | 0.5× at 8k nnz (GPU loses), 3.3× at 128k, **5.5× at 512k**, 4.8× at 2M | |
| Refinery planning LP, 5.6·10⁵ nnz, raw PDHG 1e-4 | GPU **11 s** vs CPU 27 s; rigorous bound within 1.2e-4 | exact optimum in 92 s |
| Crude valuation family (64 cases) | exact certified curve 0.16 s · warm simplex chain 0.19 s · batched GPU PDHG 4.8 s · cold solves 9.9 s | |
| Engine router (99 LPs, 5-fold held-out) | 1.06× the oracle's time; picks the fastest engine 89% of the time | |

**Honest reading.**
- PRAMANA is correct everywhere it claims anything, and more robust than HiGHS on badly scaled models.
- HiGHS is faster on large LPs and hard MIPs, often by an order of magnitude.
- The GPU wins per PDHG iteration from about 3·10⁴ nonzeros. For refinery-sized LPs at 1e-9 accuracy, warm-started simplex still wins end to end, which is exactly what the router learns.
- The MIPLIB rows were produced by the binary before the final big-M presolve improvement. The final binary is re-checked in [results/mip_regression_current_binary.md](results/mip_regression_current_binary.md): 0 wrong answers.

---

## 2. Requirements

| | required | notes |
|---|---|---|
| OS | Windows 10/11 x64 | Linux/macOS build with CMake (section 3) |
| Compiler | Visual Studio 2022 or newer, or Build Tools, with **"Desktop development with C++"** | ships MSVC, CMake and Ninja; `build.bat` finds it automatically |
| Python | 3.10+ with `numpy scipy matplotlib pytest` | only for tests, data download, generators, the exact verifier and benchmarks. **The solver itself needs no Python.** scipy is used only as the HiGHS *reference* in `bench/` |
| GPU (optional) | any NVIDIA GPU + driver | kernels are compiled at runtime by NVRTC: install the CUDA toolkit **or** `pip install nvidia-cuda-nvrtc-cu12`. Without a GPU everything runs on the CPU |
| Internet | once | to download Netlib / MIPLIB / Maros–Mészáros (`bench/fetch_data.py`); the repository already contains them under `data/` |

```bat
pip install numpy scipy matplotlib pytest
```

---

## 3. Build and run in 5 minutes

### One command: the PRAMANA console

```bat
start.bat            :: Windows      (Linux/macOS: ./start.sh)
```

`start.bat` does everything in one go:
1. finds Python;
2. builds the solver (incrementally; a no-op when nothing has changed);
3. downloads or generates the benchmark and industrial models if they are missing;
4. opens the **PRAMANA console**, an interactive terminal dashboard (no browser, no extra packages).

```
╭─── PRAMANA v1.0 ────────────────────────────────────────────────────────────────────────────────╮
│               Welcome back, <user>!                        │ Getting started                    │
│                          ▁▁▁▁▁◆ ↗                          │ Type a model name, e.g. afiro      │
│         ╱░░░░░░╲    █▀█ █▀█ ▄▀█ █▀▄▀█ ▄▀█ █▄ █ ▄▀█         │ / for commands · Tab completes     │
│        ▕░░░░░░░░▏   █▀▀ █▀▄ █▀█ █ ▀ █ █▀█ █ ▀█ █▀█         │ ────────────────────────────────── │
│          ╲░░░░░░╱    certified optimization engine         │ Recent activity                    │
│   LP · MILP · QP  ·  16 threads  ·  RTX 4050 Laptop GPU    │ p0201        OPTIMAL       2.59 s  │
╰─────────────────────────────────────────────────────────────────────────────────────────────────╯
────────────────────────────────────────────────────────────────────────────────────────────────────
> Try "/analyze pilot87"
────────────────────────────────────────────────────────────────────────────────────────────────────
 PRAMANA  GPU Acceltor  engine auto · 300s              ● afiro OPTIMAL <1 ms  GPU RTX 4050 · 16 thr
```

| type | what you get |
|---|---|
| `afiro`, `p0201`, `pilot87`, … (any model name or path) | solve + certify: status, objective, certified bound, engine and router reasoning, time breakdown, every certificate check, automatic **analysis** of what happened |
| `/analyze <model>` | structure and numerics **without solving**: row/bound types, coefficient ranges in decades, big-M detection, predicted time per engine |
| `/debug [run#]` | **numerical diagnostics** of a run: degeneracy, scaling, refactorizations, perturbations, stall recoveries, where simplex time went, cuts and root-gap closure, gap over time, PDHG kernel/transfer time, router prediction vs measured |
| `/verify [exact]` | independent re-check of the last result: C++ interval certifier + exact-rational verifier (optionally exact optimal basis) |
| `/compare <model>` | every LP engine (dual simplex, IPM, PDHG CPU, PDHG GPU) on the same model, certified times side by side |
| `/param <model> <col> <kind> <from> <to>` | certified parametric analysis with a value-curve chart and every breakpoint |
| `/bench` · `/gpu` · `/models` · `/history` · `/set` · `/check` · `/help` | benchmark results vs HiGHS, GPU crossover, model browser, run history, session settings, full self-check |

Keys: `Tab` completes commands and model names, `↑`/`↓` browse history or suggestions, `Esc` clears,
`Ctrl+C` cancels a running solve (twice on an empty line exits). Commands can also be piped in for
scripted use: `echo afiro | start.bat`.

### Manual build and first solve

```bat
build.bat                                   :: -> build\pramana.exe, build\libpramana.dll, build\pramana_tests.exe
build\pramana.exe info                      :: threads, GPU detection
build\pramana.exe data\netlib\afiro.mps.gz  :: first solve: status, objective, certificate
```

Expected output of the last command (abridged):

```
Engine           dual-simplex
Status           OPTIMAL
Objective        -464.753142857
Certificate      ACCEPTED
  primal_row_violation           PASS   ...
  safe_dual_bound_gap            PASS   ...
```

Linux/macOS: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j`
(the GPU is used via `libcuda.so` + `libnvrtc.so` when present).

If the benchmark data is missing (e.g. a trimmed checkout):

```bat
python bench\fetch_data.py                   :: Netlib, Netlib-infeasible, MIPLIB 3, Maros-Meszaros
python gen\refinery.py suite --out data\gen  :: refinery / scheduling / unit commitment / dispatch models
python gen\adversarial.py --out data\adversarial
```

---

## 4. How to check this submission

### One command

```bat
check.bat
```

`check.bat` builds, then runs everything below. It prints `PASS`/`FAIL` per item and exits with 0
only if all pass. On the reference machine: **21 passed, 0 failed, about 1 minute**.

| # | check | what it proves |
|---|---|---|
| 1 | C++ unit tests (24 groups, 32k assertions) | LU, LDLᵀ, interval arithmetic, certifier rejects false claims, simplex/IPM/PDHG on known optima, CPU = GPU, presolve round trip, MIP vs brute force, MIPLIB optima, parametric vs re-solves, determinism |
| 2 | Python end-to-end tests | bindings, exact verifier, generators, terminal console |
| 3 | from-scratch dependency audit | the binaries import only OS/C++ runtime DLLs, and no solver headers are included ([docs/FROM_SCRATCH_AUDIT.md](docs/FROM_SCRATCH_AUDIT.md)) |
| 4 | 14 certified solves | LP (AFIRO, ill-conditioned PILOT87, degenerate DEGEN3, Klee–Minty, badly scaled), infeasible (Farkas), unbounded (ray), MILP (MIPLIB P0033/P0282, unit commitment, crude unloading), QP (HS21, economic dispatch), refinery planning |
| 5 | independent re-verification | C++ `verify` of a saved result, exact-rational verification including the **exact optimal basis**, exact Farkas check |
| 6 | certified parametric analysis | refinery crude-valuation curve |

If Python is not on `PATH`: `set PYTHON=C:\path\to\python.exe` before `check.bat`.

### Check any result yourself

```bat
build\pramana.exe data\netlib\pilot87.mps.gz --json out.json --vectors   :: solve and save x, y, basis, certificate
build\pramana.exe verify data\netlib\pilot87.mps.gz out.json              :: re-check with the C++ certifier
cd python
python -m pramana.verify ..\data\netlib\pilot87.mps.gz ..\out.json --exact-basis   :: own parser, exact rationals
```

### Guided demo (10 steps, for evaluators)

`demo.bat` walks through the audit, an ill-conditioned LP, exact re-verification, a badly scaled model,
a Farkas proof, CPU-vs-GPU racing, cuts on/off, MILP racing and crude valuation. See
[docs/DEMO.md](docs/DEMO.md) for the script and talking points.

---

## 5. Using the solver

### Command line

```bat
:: LP / MILP / QP: the router picks the engine; certificate printed
build\pramana.exe data\netlib\pilot87.mps.gz
build\pramana.exe data\miplib3\bell5.mps.gz --time 60 --mip-threads 4
build\pramana.exe data\maros\CVXQP2_M.qps

:: choose the engine
build\pramana.exe data\gen\plan_20x52.mps --algo dual|ipm|pdhg-cpu|pdhg-gpu|race

:: raw GPU first-order solve (no crossover), still with a rigorous bound
build\pramana.exe data\gen\plan_40x104.mps --algo pdhg-gpu --no-crossover --pdhg-tol 1e-4

:: full JSON: certificate, telemetry, named primal/dual vectors, basis
build\pramana.exe model.mps --json result.json --vectors

:: refinery crude valuation: exact certified value curve / strategy comparison
build\pramana.exe parametric data\gen\plan_10x12.mps --col BUY_CR03_0 --kind cost  --from -95 --to -35 --json param.json
build\pramana.exe family     data\gen\plan_10x12.mps --col BUY_CR03_0 --kind upper --from 0 --to 120 --cases 64

:: structure / numerics report and engine prediction, without solving
build\pramana.exe analyze data\netlib\pilot87.mps.gz

:: measure where the GPU starts to win on this machine
build\pramana.exe calibrate

build\pramana.exe --help                    :: all options
```

Exit codes: `0` proven (OPTIMAL / INFEASIBLE / UNBOUNDED, certificate accepted), `2` not proven
(time/iteration limit), `3` numerical failure detected, `4` verification rejected, `10` error.

### Viewer

In the terminal, the PRAMANA console (`start.bat`, section 3) shows everything. For a browser view, open `web\index.html` in a browser and drop a `--json` result or a parametric JSON onto it. It shows
the status, the certificate checks, telemetry and the value curve.

### Python

```python
import sys; sys.path.insert(0, "python")
import pramana                                  # loads build\libpramana.dll
m = pramana.Model("demo", maximize=True)
x = m.add_var(0, 4, cost=3); y = m.add_var(0, cost=5, integer=True)
m.add_constraint({x: 1, y: 2}, ub=8)
r = m.solve()                                   # r.status, r.objective, r.certified, r.row_duals, r.certificate
r2 = pramana.solve_file("data/netlib/afiro.mps.gz", algorithm="dual")
```

### C

`include/pramana/pramana.h`: create/read a model, add columns, rows and Q entries, solve with a JSON
options string, and read status, objective, bound, certificate, x, duals, JSON, parametric and family
results. See [docs/API.md](docs/API.md).

---

## 6. Reproducing every benchmark

```bash
bash bench/run_all.sh          # several hours; sequential so timings do not interfere
```

This fetches the data, generates the industrial and adversarial models, runs PRAMANA and HiGHS on
every set, and runs the GPU crossover, the crude-valuation family experiment and router fitting. It
writes `results/*/runs.csv`, `results/*/summary.md`, `docs/figures/*.png` and `docs/BENCHMARKS.md`.

Individual pieces:

```bat
python bench\run_bench.py --set netlib  --engines dual --reference --time 300
python bench\run_bench.py --set miplib3 --engines auto --reference --time 60
python bench\run_bench.py --set miplib3 --engines auto --time 60 --extra=--no-cuts --tag nocuts
python bench\gpu_crossover.py
python bench\family_experiment.py
python bench\fit_router.py
python bench\report.py        :: rebuild docs\BENCHMARKS.md from results\
python bench\audit_deps.py    :: rebuild docs\FROM_SCRATCH_AUDIT.md
```

The instance lists are pre-registered in `bench/sets/`. Nothing is filtered: timeouts count at the
limit in shifted geometric means, and every failure is listed.

Timing note: the reference machine is a laptop, and timings vary by up to ~2× between runs (thermal and
power state). Compare solvers within the same run.

---

## 7. Problem statement coverage

The full mapping, clause by clause, is in [docs/PS_TRACEABILITY.md](docs/PS_TRACEABILITY.md).

| PS asks for | PRAMANA |
|---|---|
| LP, MILP, QP; extensible to MIQP / NLP / MINLP | all three, each certified; relaxation-agnostic B&B and shared KKT/LDLᵀ as extension points ([docs/ROADMAP.md](docs/ROADMAP.md)) |
| Revised simplex, interior point | bounded dual + primal simplex with own sparse LU; Mehrotra IPM with own LDLᵀ + ordering; plus restarted PDHG |
| B&B, B&C, cutting planes, presolve, heuristics, node selection | reliability branching, best-bound + plunging, GMI / MIR / cover cuts, bound propagation, rounding / diving / feasibility pump, presolve with primal + dual postsolve and big-M coefficient tightening |
| Sparse techniques, efficient numerical LA | CSC/CSR, hyper-sparse FTRAN/BTRAN, Markowitz LU with updates, minimum-degree LDLᵀ |
| Multi-core | parallel PDHG, CPU‖GPU LP racing, concurrent MIP racing (`--mip-threads`) |
| GPU where it gives *measurable* benefit | own CUDA kernels, measured per-iteration crossover, end-to-end timings incl. transfers, validated router |
| Numerical stability; degeneracy, ill-conditioning, weak relaxations | scaling, perturbation, Harris + bound-flipping ratio test, refactor on residual, certify-before-claim; adversarial suite, DEGEN/PILOT, big-M vs tight scheduling |
| MIPLIB / Netlib / Mittelmann / QPLIB benchmarks vs an established solver | Netlib (all), MIPLIB 3 subset, Maros–Mészáros QP, all vs HiGHS on the same machine |
| Built from scratch, no open-source solver library | dependency audit in CI style (`bench/audit_deps.py`), [docs/FROM_SCRATCH_AUDIT.md](docs/FROM_SCRATCH_AUDIT.md) |
| API or CLI | CLI, C API (DLL), Python package |
| Refinery / blending / scheduling / power / supply-chain cases | refinery planning LP/MILP, crude blending, crude unloading and tank scheduling, unit commitment, economic dispatch QP: literature structures, synthetic data (`gen/`) |

---

## 8. Repository layout

```
start.bat  start.sh             one command: build + data + PRAMANA console
build.bat  check.bat  demo.bat     build / self-check / guided demo (Windows)
CMakeLists.txt                     targets: pramana_core (static), libpramana (DLL), pramana (CLI), pramana_tests
include/pramana/pramana.h          C API
src/
  core/        sparse matrices, model
  io/          MPS/QPS reader + writer, own gzip decoder
  lp/          sparse LU, dual/primal simplex, scaling, crossover
  ipm/         minimum-degree ordering, LDLᵀ, Mehrotra IPM (LP + QP), convexity certificate
  pdhg/        restarted Halpern PDHG, CPU + CUDA backends, batched solves, CUDA kernel source
  gpu/         CUDA driver / NVRTC runtime loader (no CUDA toolkit needed at build time)
  mip/         branch-and-cut, cuts, propagation
  presolve/    reductions + primal/dual/basis postsolve
  cert/        interval-arithmetic certifier, safe dual bounds
  router/      engine cost model
  parametric/  certified parametric analysis, case-family strategies
  api/  cli/   solve() pipeline, C API, command line
  util/        JSON, logging, thread pool, timers
python/pramana/                    Python bindings, independent exact-rational verifier, terminal console (tui.py)
gen/                               refinery / scheduling / UC / dispatch / adversarial model generators
bench/                             data fetch, benchmark harness vs HiGHS, GPU + family experiments, report, audit
tests/                             C++ unit tests (own framework) + pytest end-to-end tests + small models
data/                              public benchmark sets (Netlib, MIPLIB 3, Maros–Mészáros) + generated models
results/                           raw benchmark results (runs.csv, summary.md per set; logs in results/logs)
docs/                              documentation (below)
web/index.html                     static result / certificate / value-curve viewer
```

| document | content |
|---|---|
| [ARCHITECTURE](docs/ARCHITECTURE.md) | modules, data flow, engine interfaces, extension points |
| [ALGORITHMS](docs/ALGORITHMS.md) | the mathematics of every engine and of the certificates |
| [PS_TRACEABILITY](docs/PS_TRACEABILITY.md) | problem statement → implementation → evidence |
| [VALIDATION](docs/VALIDATION.md) | validation matrix and how each claim is tested |
| [BENCHMARKS](docs/BENCHMARKS.md) | all benchmark tables and figures (generated) |
| [FROM_SCRATCH_AUDIT](docs/FROM_SCRATCH_AUDIT.md) | dependency audit (generated) |
| [API](docs/API.md) | CLI, C API, Python API, JSON formats |
| [DEMO](docs/DEMO.md) | 7-minute demo script |
| [ROADMAP](docs/ROADMAP.md) | honest limitations and next steps |

---

## 9. What we claim, and what we do not

**We claim**
- The solver core is our own code. The dependency audit is attached and reproducible.
- Every reported proof (optimal, infeasible, unbounded, MIP bound) is independently checkable, and exactly checkable for LPs.
- The measured results against HiGHS on public sets, with every failure listed.
- The GPU is faster than our multi-threaded CPU for PDHG above a measured size and slower below it, and the router chooses accordingly.

**We do not claim**
- Parity with Gurobi, CPLEX or HiGHS on speed. HiGHS is faster on large LPs and hard MIPs.
- A new GPU LP algorithm; PDHG follows PDLP / cuPDLPx.
- MRPL data: all industrial models are synthetic, generated from published model structures.
- That "millions of variables" are solved to certified optimality. The largest certified optimum is 86k columns / 62k rows / 252k nnz; at 192k columns / 562k nnz we give a rigorous 1e-4 bound. See [docs/ROADMAP.md](docs/ROADMAP.md).
