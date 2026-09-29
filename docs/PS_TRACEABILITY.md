# Problem statement → implementation → evidence (SIH26119, MRPL)

Numbers quoted here are produced by `bench/run_all.sh` and summarized in `docs/BENCHMARKS.md`
(per-instance tables, including every failure, in `results/<set>/summary.md`).

## Explicit requirements

| # | PS requirement | Implementation | Evidence |
|---|---|---|---|
| 1 | Solve **LP** | bounded dual simplex + primal cleanup (`src/lp/simplex.cpp`), IPM + crossover (`src/ipm`), PDHG CPU/GPU + crossover (`src/pdhg`) | Netlib 91/91 certified optimal, 29/29 infeasible certified (`results/netlib*`) |
| 2 | Solve **MILP** | branch-and-cut (`src/mip/bnb.cpp`, `cuts.cpp`, `propagate.cpp`) | MIPLIB 3 subset (`results/miplib3`), refinery MILPs (`results/gen`) |
| 3 | Solve **QP** (convex) | Mehrotra IPM on quasidefinite KKT (`src/ipm/ipm.cpp`), convexity check with certificate | Maros–Mészáros (`results/maros`), dispatch QP (`data/gen`) |
| 4 | Extensible to MIQP / NLP / MINLP | B&B is relaxation-agnostic (node = bounds + warm basis); IPM handles Q; parametric/sequential-LP driver structure; explicit error + hook for MIQP | `docs/ROADMAP.md` |
| 5 | Revised simplex / interior point | both implemented from scratch (own LU, own LDLᵀ + ordering) | unit tests `lu_*`, `ldl_*`, `simplex_*`, `ipm_*` |
| 6 | Branch-and-bound / branch-and-cut / cutting planes | reliability branching, GMI + c-MIR + cover cuts with safety filters | cut ablation in `docs/BENCHMARKS.md`; UC root solve 0.09 s with cuts vs 29.9 s without |
| 7 | Presolve | 10 reduction rules + full primal/dual/basis postsolve (QP-safe subset) | unit test `presolve_postsolve_roundtrip_certifies`; telemetry `presolve` in every result |
| 8 | Heuristics | rounding, fractional diving, feasibility pump, fix-and-resolve | telemetry `mip.heuristic_solutions`, `incumbent_source` |
| 9 | Node selection | best bound with depth-first plunging | `src/mip/bnb.cpp` |
| 10 | Sparse matrix techniques | CSC/CSR, hyper-sparse vectors, sparse LU with Markowitz, sparse LDLᵀ with minimum degree | `src/core/sparse.*`, `src/lp/lu.*`, `src/ipm/ldl.*` |
| 11 | Multi-core | thread-pool parallel PDHG (SpMV + fused vector ops), **race mode** (simplex ‖ GPU-PDHG/IPM, cooperative cancel) | `--algo race`, `--threads N`; telemetry `race` |
| 12 | GPU acceleration **where it provides measurable benefit** | own CUDA kernels (NVRTC at runtime), batched PDHG, measured per-iteration crossover, validated router | `results/gpu`, `results/router`, `docs/figures/gpu_*.png` |
| 13 | Numerical stability | scaling, perturbation, Harris/BFRT, per-variable unscaled tolerances, refactor on residual, rank repair, stall recovery, dual polish, IPM regularization + refinement | adversarial suite: badly-scaled models solved to the true optimum with certificates |
| 14 | Benchmarks on MIPLIB / Netlib / Mittelmann-type / QPLIB-type sets | Netlib (all), Netlib infeasible (all), MIPLIB 3 subset (pre-registered), Maros–Mészáros (124), refinery family up to 136k rows | `docs/BENCHMARKS.md` |
| 15 | Compared against ≥ 1 established solver | HiGHS (via scipy) on the same machine and limits, never linked into the solver | every `results/*/summary.md` |
| 16 | Degeneracy / weak relaxations / ill-conditioning demonstrated | DEGEN2/3, assignment & transportation (degenerate), Klee–Minty, big-M vs tight scheduling, scaled-badly suite, PILOT family | `results/adversarial`, `results/gen` |
| 17 | Not built on any open-source solver library | own code only; binaries import only OS + C++ runtime | `docs/FROM_SCRATCH_AUDIT.md` (`bench/audit_deps.py`) |
| 18 | API / CLI | `pramana` CLI, C API (`include/pramana/pramana.h`), Python package (`python/pramana`) | `docs/API.md`, `tests/python` |
| 19 | Industrial use cases (refining, energy, logistics) | refinery planning LP/MILP, crude blending, crude unloading & tank scheduling, unit commitment, economic dispatch QP (literature structures, synthetic data) | `gen/refinery.py`, `results/gen` |

## Implicit expectations

| expectation | how PRAMANA meets it |
|---|---|
| Correctness beats speed | certify-before-claim everywhere; two independent checkers (C++ interval certifier + Python exact-rational verifier with its own parser) |
| Honest GPU usage | GPU times include transfers and setup; the router picks CPU when the GPU loses; the family experiment shows where the GPU does **not** pay |
| Ownership = inspect & modify | 11.6k lines of readable C++ with structured telemetry (pivots, refactors, perturbations, cuts, kernel/transfer times) |
| Failure honesty | full failure tables; statuses downgraded when unproven |
| Modular, extensible | engine interfaces behind `solve()`, backend interface for PDHG (CPU/CUDA; HIP/SYCL future), relaxation-agnostic B&B |

## Differentiators (beyond the PS)

1. **Proof-carrying results**: every OPTIMAL / INFEASIBLE / UNBOUNDED claim carries a certificate
   (rigorous safe dual bound, Farkas ray, primal ray, QP linearization bound, MIP safe bound),
   re-checkable in exact rational arithmetic by `python -m pramana.verify`.
2. **Certified parametric analysis** for refinery crude valuation: the exact value function and
   every breakeven point in one run, each segment certified (vs hundreds of sampled case runs that
   miss breakpoints).
3. **Measured, validated CPU/GPU routing** with published regret on held-out instances.
4. **GPU without a CUDA toolkit**: kernels compiled at runtime by NVRTC through the driver API;
   the same binary runs CPU-only machines.
