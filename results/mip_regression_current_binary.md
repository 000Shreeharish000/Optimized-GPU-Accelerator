# MILP regression check of the final binary (30 s limit)

The MIPLIB tables in `docs/BENCHMARKS.md` were produced by the binary *before* the last presolve change
(big-M coefficient tightening for negative binary coefficients, `src/presolve/presolve.cpp`). The table
below re-runs every MIPLIB 3 subset instance and every industrial MILP with the final binary at a 30 s
limit and compares it with the HiGHS optimum. Command: `build\pramana.exe <model> --time 30`.

**Result: 0 wrong answers.** Every OPTIMAL agrees with HiGHS within the default 1e-4 MIP gap, and every
certificate is accepted.
- Newly solved: gesa2 (0.9 s), gesa2_o (8.2 s), bell5 (13.3 s).
- unload_bigm now finds the true optimum 109.49 (the old binary was stuck at 259–394) with a certified 6.9% gap.
- Lost within 30 s: dcmulti and qnet1. The bound reaches the optimum, but the primal heuristics find no better incumbent.

| instance | status | objective | HiGHS optimum | seconds |
|---|---|---|---|---|
| 10teams | TIME_LIMIT | – | 924 | 30.0 |
| bell3a | OPTIMAL | 878430.316 | 878430.316 | 4.13 |
| bell5 | OPTIMAL | 8966406.49 | 8966406.49 | 13.29 |
| blend2 | OPTIMAL | 7.598985 | 7.598985 | 1.95 |
| blend_12_mip | OPTIMAL | 67489.42086 | 67489.42086 | 0.00 |
| blend_30_mip | OPTIMAL | 66451.92534 | 66451.92534 | 0.01 |
| dcmulti | TIME_LIMIT | 188359.5 | 188182 | 30.0 |
| egout | OPTIMAL | 568.1007 | 568.1007 | 0.01 |
| enigma | OPTIMAL | 0 | 0 | 0.64 |
| fixnet6 | TIME_LIMIT | 3983.0 | 3983 | 30.1 |
| flugpl | OPTIMAL | 1201500 | 1201500 | 0.03 |
| gesa2 | OPTIMAL | 25780031.43 | 25779856.37 | 0.92 |
| gesa2_o | OPTIMAL | 25779856.37 | 25779856.37 | 8.23 |
| gesa3 | OPTIMAL | 27991042.65 | 27991042.65 | 1.18 |
| gesa3_o | OPTIMAL | 27991042.65 | 27991042.65 | 1.55 |
| gt2 | OPTIMAL | 21166 | 21166 | 0.05 |
| khb05250 | OPTIMAL | 106940226 | 106940226 | 0.18 |
| lseu | OPTIMAL | 1120 | 1120 | 0.87 |
| misc03 | OPTIMAL | 3360 | 3360 | 1.49 |
| misc06 | OPTIMAL | 12850.86 | 12851.08 | 0.13 |
| misc07 | TIME_LIMIT | 2810 | 2810 | 30.0 |
| mod008 | OPTIMAL | 307 | 307 | 0.54 |
| mod010 | OPTIMAL | 6548 | 6548 | 0.03 |
| p0033 | OPTIMAL | 3089 | 3089 | 0.01 |
| p0201 | OPTIMAL | 7615 | 7615 | 2.74 |
| p0282 | OPTIMAL | 258411 | 258411 | 0.97 |
| p0548 | OPTIMAL | 8691 | 8691 | 0.43 |
| p2756 | OPTIMAL | 3124 | 3124 | 1.99 |
| plan_10x12_mip | OPTIMAL | 40545.31158 | 40545.31158 | 24.46 |
| plan_10x6_mip | OPTIMAL | 21815.62546 | 21815.62545 | 4.49 |
| plan_6x4_mip | OPTIMAL | 13985.02229 | 13985.02229 | 0.59 |
| pp08a | TIME_LIMIT | 7370 | 7350 | 30.0 |
| pp08aCUTS | TIME_LIMIT | 7390 | 7350 | 30.0 |
| qnet1 | TIME_LIMIT | 16046.74 | 16029.69 | 30.0 |
| qnet1_o | OPTIMAL | 16030.99 | 16029.69 | 15.79 |
| rgn | OPTIMAL | 82.2 | 82.2 | 0.77 |
| set1ch | TIME_LIMIT | 96824.54 | 54542.25 | 30.0 |
| stein27 | OPTIMAL | 18 | 18 | 0.46 |
| stein45 | OPTIMAL | 30 | 30 | 19.80 |
| uc_10x24 | OPTIMAL | 644229.2463 | 644229.2463 | 0.05 |
| uc_30x24 | TIME_LIMIT | 2067219.14 | 1933392.05 | 30.0 |
| unload_bigm | TIME_LIMIT | 109.4911 | 109.4911 | 30.0 |
| unload_tight | OPTIMAL | 119.4911 | 119.4911 | 10.47 |
| vpm1 | OPTIMAL | 20 | 20 | 0.23 |
| vpm2 | TIME_LIMIT | 14.5 | 13.75 | 30.0 |

Differences such as gesa2 (6.8e-6 relative), misc06 (1.7e-5) and qnet1_o (8.1e-5) are within the 1e-4
relative MIP gap that both solvers use by default. Run with `--gap 0` for exact optima.
