# Roadmap and honest limitations

## Known limitations (measured, see `docs/BENCHMARKS.md`)

- **Speed vs mature solvers.** On small Netlib LPs PRAMANA's dual simplex is within ~2× of HiGHS
  (shifted geometric mean); on larger LPs (10⁴–10⁵ rows) HiGHS is roughly an order of magnitude
  faster. Missing: Forrest–Tomlin updates, hyper-sparse FTRAN/BTRAN with DFS, partial pricing,
  parallel dual simplex.
- **QP IPM.** A subset of Maros–Mészáros (LISWET*, HUES*, POWELL20, QSCAGR25, YAO, UBH1, ...) does
  not converge within the limit; these are reported as ITERATION_LIMIT, never as wrong answers.
  Planned: Gondzio multiple centrality correctors, homogeneous self-dual embedding, QP crossover.
- **MILP.** Hard MIPLIB instances time out (reported with certified gaps). Missing: conflict
  analysis, symmetry, restarts, more cut families (flow cover, clique), parallel tree search.
- **Tolerance semantics.** Models infeasible by less than the feasibility tolerance (e.g.
  `infeasthin_50`, 1e-6 absolute) are reported as unproven rather than as infeasible.

## v2 (next 3 months)

- Forrest–Tomlin update and hyper-sparse solves; bound-flipping in primal; parallel strong branching.
- Deterministic epoch-parallel B&B (option `mipThreads` is reserved in `SolverOptions`).
- Rational re-verification of the final LP basis inside the C++ binary (today: Python `--exact-basis`).
- VIPR-style MIP tree certificates for small instances.
- GPU: HIP/SYCL backends behind the existing `Backend` interface; FP32-storage SpMV experiment.

## v3

- **MIQP**: B&B over QP relaxations (the tree code is relaxation-agnostic; IPM warm starts needed).
- **NLP / MINLP**: sequential LP / distributive recursion driver for refinery pooling (warm-started
  simplex per pass), then a filter-line-search IPM reusing the LDLᵀ.
- Refinery structure detection (pooling-aware bounds) in generic MPS input.
