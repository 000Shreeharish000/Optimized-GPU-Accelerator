# Roadmap and honest limitations

## Known limitations (measured, see `docs/BENCHMARKS.md`)

- **Speed vs mature solvers.** On small Netlib LPs PRAMANA's dual simplex is within ~2× of HiGHS
  (shifted geometric mean); on larger LPs (10⁴–10⁵ rows) HiGHS is roughly an order of magnitude
  faster. Implemented: DSE, BFRT, adaptive hyper-sparse FTRAN/BTRAN, incremental CHUZR.
  Missing: Forrest–Tomlin updates, partial/multiple pricing, parallel dual simplex.
- **QP IPM.** A subset of Maros–Mészáros (LISWET*, HUES*, KSIP, YAO, UBH1, ...) does
  not converge within the limit; these are reported as ITERATION_LIMIT, never as wrong answers.
  Planned: Gondzio multiple centrality correctors, homogeneous self-dual embedding, QP crossover.
- **MILP.** Hard MIPLIB instances time out (reported with certified gaps). Missing: conflict
  analysis, symmetry, restarts, more cut families (flow cover, clique), a shared-tree parallel search
  (today: concurrent diversified racing with `--mip-threads`).
- **Scope of MIP certificates.** The incumbent is checked independently (feasibility, integrality,
  objective), and the global bound is the minimum of Neumaier–Shcherbina safe node bounds. That bound
  is only as valid as the presolve reductions and cutting planes beneath it. A presolve bug found
  during development (coefficient tightening applied to the right-hand side but not the matrix) produced
  a wrong "certified" optimum on p2756 that only the reference comparison exposed. It is fixed and
  covered by unit tests; the planned remedy is VIPR-style tree certificates that re-derive every cut and
  reduction.
- **Weak big-M relaxations.** Presolve now shrinks big-M coefficients to implied bounds
  (Savelsbergh tightening for both coefficient signs). On the crude-unloading big-M model this lifts
  the root bound from 52 to 96 and finds the true optimum 109.49, but the final 6.9% gap is not closed
  in 60 s (HiGHS: 0.17 s). On MIPLIB it newly solves gesa2, gesa2_o and bell5. dcmulti and qnet1 now
  need more than 30 s: their bound reaches the optimum, but the heuristics miss the optimal incumbent
  (`results/mip_regression_current_binary.md`). Next: RINS / local branching in the tree, implied-bound
  cuts.
- **Largest planning LP (plan_60x156, 5.6·10⁵ nnz).** No PRAMANA engine proves optimality within
  300 s in the final run (HiGHS: 92 s). Raw GPU PDHG gives a rigorous bound within ~1e-4 in 11 s, but
  the crossover to an exact vertex does not finish in time. Next step: a PDHG-to-simplex handoff that
  starts crossover from a 1e-6 point with a partial basis instead of a full basis guess.
- **Timing variance and one unexplained GPU overrun.** The laptop used for the benchmarks varies by up
  to ~2× between runs (thermal/power state). One GPU run (plan_60x156, pdhg-gpu + crossover) took
  651 s under a 300 s limit, with 6× slower kernels than the neighbouring runs. The deadline is checked
  every 64 iterations, so we suspect a power-state event; it has not been reproduced. Results must be
  compared within one run, and the benchmark should be repeated on a desktop/datacenter GPU.
- **Tolerance semantics.** A model infeasible by less than the feasibility tolerance is proven
  infeasible only if a presolve activity bound yields a one-row Farkas certificate (e.g.
  `infeasthin_50`); otherwise it is reported as unproven, never as optimal.

## v2 (next 3 months)

- Forrest–Tomlin update; bound-flipping in primal; parallel strong branching.
- Deterministic epoch-parallel shared-tree B&B (beyond today's concurrent racing).
- Rational re-verification of the final LP basis inside the C++ binary (today: Python `--exact-basis`).
- VIPR-style MIP tree certificates for small instances.
- GPU: HIP/SYCL backends behind the existing `Backend` interface; FP32-storage SpMV experiment.

## v3

- **MIQP**: B&B over QP relaxations (the tree code is relaxation-agnostic; IPM warm starts needed).
- **NLP / MINLP**: sequential LP / distributive recursion driver for refinery pooling (warm-started
  simplex per pass), then a filter-line-search IPM reusing the LDLᵀ.
- Refinery structure detection (pooling-aware bounds) in generic MPS input.
