# PRAMANA — a certified, from-scratch LP / MILP / QP solver with measured GPU acceleration

*SIH26119 (MRPL): indigenous GPU-accelerated optimization solver.*
*pramāṇa (Sanskrit): "valid means of knowledge, proof".*

PRAMANA is an optimization solver written from mathematical foundations: its own MPS reader, sparse
LU, dual/primal simplex, sparse LDLᵀ with minimum-degree ordering, interior point method, restarted
PDHG with its own CUDA kernels, presolve, and branch-and-cut. Its defining property is that **it
never asks to be trusted**:

- **Every answer is a proof.** OPTIMAL comes with a rigorous (outward-rounded) Neumaier–Shcherbina
  dual bound, INFEASIBLE with a verified Farkas certificate, UNBOUNDED with a verified primal ray,
  NONCONVEX with a vector v where vᵀQv < 0. If a proof fails, the status is downgraded, never
  reported as proven.
- **Two independent checkers.** The C++ certifier (interval arithmetic) and a Python verifier with
  its *own* MPS parser working in exact rational arithmetic (including exact vertex verification).
- **GPU only where it measurably pays.** Own CUDA kernels, compiled at runtime (no CUDA toolkit
  needed), batched PDHG, a measured CPU/GPU crossover, and a router validated on held-out instances.
- **Built for refinery planning.** Certified parametric analysis gives the exact crude-valuation
  curve (every breakeven price, every segment certified) in one run instead of hundreds of case
  re-solves that miss breakpoints.

## Headline results (this machine: see `results/machine.json`; details in [docs/BENCHMARKS.md](docs/BENCHMARKS.md))

| benchmark | PRAMANA | reference (HiGHS via scipy) |
|---|---|---|
| Netlib LP, all 91 feasible | **91/91 optimal, all certified**, 0 disagreements | 91/91 |
| Netlib infeasible, all 29 | **29/29 infeasible, all with verified Farkas proofs** | 28/29 |
| Badly scaled Netlib (10^±6 rescaling, known true optimum) | **7/8 true optimum + certificate, 0 wrong answers** (1 honestly unproven) | 0/8 correct: 4 wrong optima, 3 false INFEASIBLE, 1 false UNBOUNDED (default tolerances) |
| Adversarial suite (17: Klee–Minty, degenerate, near-singular, thin-infeasible, unbounded, scaled) | 16/17 proven + certified, **0 wrong** | 17/17 answered, **8 wrong** (all on the scaled models) |
| MIPLIB 3 pre-registered subset (42, 60 s) | 29/42 optimal, all certified, 0 disagreements (SGM 14.1 s) | 36/42 (SGM 5.9 s) |
| Maros–Mészáros convex QP (124, 60 s) | 105/124 optimal + certified, 0 wrong (19 unsolved, listed) | (no QP in scipy) |
| Unit commitment MILP, cuts on vs off | 0.09 s (1 node) vs 29.9 s (20,005 nodes) | |
| PDHG iteration, GPU vs 16-thread CPU (own kernels, RTX 4050 Laptop) | 0.5× at 8k nnz (GPU loses), 3.3× at 128k, **5.5× at 512k**, 4.8× at 2M | |
| Refinery planning LP, 2.5·10⁵ nnz, certified optimum via PDHG + crossover | GPU 260 s vs CPU PDHG > 300 s; dual simplex 254 s | 26.6 s |
| Refinery planning LP, 5.6·10⁵ nnz, raw PDHG 1e-4 | GPU **11 s** vs CPU 27 s, rigorous bound within 1.2e-4 (no exact vertex in 300 s) | exact optimum in 91.8 s |
| Crude valuation family (plan 10x12, availability 0–120, 64 cases) | exact certified curve (19 segments) 0.16 s; warm simplex chain 0.19 s; batched GPU PDHG 4.8 s; cold solves 9.9 s | |
| Engine router (99 LPs, 5-fold held-out) | regret 1.06× vs oracle, picks the fastest engine 89% of the time (never picks GPU below the crossover) | |

**Honest reading.** PRAMANA is correct everywhere it claims anything, and more robust than HiGHS on
badly scaled models. HiGHS is faster on large LPs and hard MIPs, often by an order of magnitude. The GPU
wins per PDHG iteration from about 3·10⁴ nnz, but for refinery-sized LPs at 1e-9 accuracy, warm-started
simplex (or HiGHS) still wins end to end; that is exactly what the router learns.

Correctness is established independently of HiGHS: where HiGHS and PRAMANA disagree (badly scaled
suite), the ground truth is known by construction and both checkers confirm PRAMANA's answers.

## Quick start (Windows, Visual Studio 2022/2026 Build Tools)

```bat
build.bat                                   :: CMake + Ninja + MSVC -> build\pramana.exe, libpramana.dll, pramana_tests.exe
python bench\fetch_data.py                  :: Netlib, Netlib-infeasible, MIPLIB 3, Maros-Meszaros
python gen\refinery.py suite --out data\gen :: refinery / scheduling / UC / dispatch models
python gen\adversarial.py                   :: degenerate, ill-conditioned, infeasible, unbounded suite
build\pramana_tests.exe                     :: 23 unit-test groups (800+ checks)
python -m pytest tests\python -q            :: end-to-end tests (bindings, exact verifier, generators)
```

`python` is any Python ≥ 3.10 with `numpy`, `scipy`, `matplotlib` and `pytest`. The solver itself needs no Python;
scipy is only used as the HiGHS *reference* in `bench/`.

Linux/macOS: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build` (GPU via
`libcuda.so` + `libnvrtc.so` if present, otherwise CPU only).

**GPU:** any NVIDIA driver provides `nvcuda.dll`. For the kernel compiler either install the CUDA
toolkit or `pip install nvidia-cuda-nvrtc-cu12` (PRAMANA finds it automatically, or set
`PRAMANA_NVRTC`). The compiled PTX is cached, so NVRTC is only needed once. Without a GPU
everything runs on the CPU. Check with `build\pramana.exe info`.

## How to run

```bat
:: solve anything (LP / MILP / QP), router picks the engine, certificate printed
build\pramana.exe data\netlib\pilot87.mps.gz
build\pramana.exe data\miplib3\bell5.mps.gz --time 60
build\pramana.exe data\maros\CVXQP2_M.qps

:: choose the engine
build\pramana.exe data\gen\plan_20x52.mps --algo dual|ipm|pdhg-cpu|pdhg-gpu|race

:: full JSON (certificate, telemetry, named primal/dual vectors, basis)
build\pramana.exe model.mps --json result.json --vectors

:: re-check a result: C++ certifier, then the independent exact-rational verifier
build\pramana.exe verify model.mps result.json
cd python && python -m pramana.verify ..\model.mps ..\result.json --exact-basis

:: refinery crude valuation: exact certified value curve / strategy comparison
build\pramana.exe parametric data\gen\plan_10x12.mps --col BUY_CR03_0 --kind cost  --from -95 --to -35
build\pramana.exe family     data\gen\plan_10x12.mps --col BUY_CR03_0 --kind upper --from 0 --to 120 --cases 64

:: GPU per-iteration crossover measurement
build\pramana.exe calibrate

:: view any result / parametric JSON in the browser
start web\index.html
```

Python:

```python
import sys; sys.path.insert(0, "python")
import pramana
m = pramana.Model("demo", maximize=True)
x = m.add_var(0, 4, cost=3); y = m.add_var(0, cost=5, integer=True)
m.add_constraint({x: 1, y: 2}, ub=8)
r = m.solve()                      # r.status, r.objective, r.certified, r.row_duals, r.certificate
```

Full reproduction of every table and figure: `bash bench/run_all.sh` (then `python bench/report.py`).

## Repository map

| path | content |
|---|---|
| `src/` | solver (≈11.6k lines C++20): `lp` simplex+LU, `ipm` IPM+LDLᵀ, `pdhg` PDHG + CUDA kernels, `gpu` driver/NVRTC loader, `mip` branch-and-cut, `presolve`, `cert` certifier, `router`, `parametric`, `api`, `cli`, `io`, `core`, `util` |
| `include/pramana/pramana.h` | C API |
| `python/pramana/` | Python bindings + independent exact verifier |
| `gen/` | refinery planning / blending / unloading / UC / dispatch generators, adversarial suite |
| `bench/` | data fetch, benchmark harness (vs HiGHS), GPU crossover, family experiment, router fitting, report, dependency audit |
| `tests/` | C++ unit tests (own framework) and pytest end-to-end tests |
| `web/index.html` | static result / certificate / value-function viewer |
| `docs/` | [ARCHITECTURE](docs/ARCHITECTURE.md), [ALGORITHMS](docs/ALGORITHMS.md), [PS_TRACEABILITY](docs/PS_TRACEABILITY.md), [VALIDATION](docs/VALIDATION.md), [BENCHMARKS](docs/BENCHMARKS.md), [FROM_SCRATCH_AUDIT](docs/FROM_SCRATCH_AUDIT.md), [API](docs/API.md), [DEMO](docs/DEMO.md), [ROADMAP](docs/ROADMAP.md) |

## What we claim, and what we do not

- **We claim:** the solver core is our own code (audit attached); every reported proof is
  independently checkable; measured results against HiGHS on public sets with every failure listed;
  the GPU is faster than our multi-threaded CPU for PDHG above a measured size and slower below it.
- **We do not claim:** parity with Gurobi/CPLEX/HiGHS on speed, a new GPU LP algorithm (PDHG follows
  PDLP/cuPDLPx), or MRPL data (all industrial models are synthetic, from published structures).
