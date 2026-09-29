# Validation matrix → where it is checked

| # | category | check | where |
|---|---|---|---|
| 1 | Parser correctness | RANGES on E/L/G rows (both signs), negative UP, MARKER, OBJSENSE, QUADOBJ, fixed-format names with blanks, gzip; writer round trip | `tests/cpp/test_core.cpp` (`mps_*`); FORPLAN/BOEING in Netlib run |
| 2 | LP correctness | all 91 Netlib vs published optima and HiGHS, rigorous safe-bound certificates | `results/netlib`, test `simplex_netlib_known_optima` |
| 3 | Infeasibility | all 29 Netlib infeasible with verified Farkas rays (C++ interval check + Python exact check) | `results/netlib-infeas`, test `test_exact_farkas_on_infeasible` |
| 4 | Unboundedness | constructed unbounded LPs with verified primal rays | `tests/data/unbounded.mps`, `data/adversarial/unbnd_30.mps` |
| 5 | Degeneracy | DEGEN2/3, assignment 80×80, degenerate transportation, Klee–Minty; stall detector + perturbation telemetry | `results/adversarial`, telemetry `simplex.degenerate_fraction`, `stall_recoveries` |
| 6 | Ill-conditioning | Netlib models rescaled by 10^U(-6,6): true optimum recovered and certified; PILOT family | `results/adversarial` ground-truth table |
| 7 | LU / update stability | refactor every 100 updates + on residual mismatch; rank repair; LU tests vs fresh factorization | `lu_*` tests, telemetry `numerical_refactors`, `singular_basis_repairs` |
| 8 | LP scalability | refinery planning family 480 → 136,656 rows, all engines vs HiGHS | `results/gen`, `results/gpu` |
| 9 | Presolve | postsolve round trip certified on original model; presolve bug found by the certifier is covered | test `presolve_postsolve_roundtrip_certifies`, Netlib run through `solve()` |
| 10 | GPU benefit | per-iteration crossover (calibrate), planning family, batched families | `results/gpu`, `results/family` |
| 11 | GPU non-benefit | small LPs: router picks CPU; GPU slower (reported) | `results/router`, `results/netlib_engines` |
| 12 | Accuracy trade-off | raw PDHG at 1e-4 / 1e-6 vs certified crossover | `results/gpu/planning_engines.csv` |
| 13 | Batched solve | batched PDHG vs warm simplex chain vs parametric | `results/family` |
| 14 | MILP correctness | MIPLIB 3 optima; brute-force knapsack comparison | `results/miplib3`, test `mip_matches_brute_force_knapsacks` |
| 15 | Weak relaxations | big-M vs tight unloading formulation; UC with/without cuts; markshare | `results/gen`, `results/miplib3_nocuts` |
| 16 | Branching ablation | `--branching mostfrac|pseudocost|reliability|strong` | CLI option; `bench/run_bench.py --extra` |
| 17 | Heuristics | time to first incumbent, incumbent source | telemetry `mip.first_incumbent_seconds`, `incumbent_source` |
| 18 | Safe bounds | every MIP node bound is a Neumaier–Shcherbina bound; random duals always give valid bounds | test `safe_bound_valid_for_any_dual`, telemetry `mip.safe_bound_corrections` |
| 19 | QP | Maros–Mészáros 124 instances, KKT + linearized safe bound | `results/maros` |
| 20 | Determinism | identical results and iteration counts on repeated runs | test `determinism_identical_runs` |
| 21 | Industrial realism | literature model structures; duals = marginal crude values (parametric slopes) | `gen/refinery.py`, `results/family` |
| 22 | From-scratch audit | import table + source scan | `bench/audit_deps.py` → `docs/FROM_SCRATCH_AUDIT.md` |
| 23 | Timeouts | incumbent + certified bound returned with TIME_LIMIT | `results/miplib3`, `uc_30x24` |
