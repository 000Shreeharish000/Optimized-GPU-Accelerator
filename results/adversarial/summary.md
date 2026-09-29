# Benchmark: adversarial

Instances: 17; time limit 60s; machine: see results/machine.json.
Status = certified status (PRAMANA) / solver status (HiGHS). Time excludes file reading.

| solver | solved (proven) | certified | SGM time (s, shift 10, unsolved at limit) | wrong vs reference |
|---|---|---|---|---|
| auto | 15/17 | 15 | 2.624 | 4 |
| highs | 17/17 | - | 0.007 | - |

## Full per-instance table (including every failure)

| instance | auto status | auto obj | auto time | highs status | highs obj | highs time |
|---|---|---|---|---|---|---|
| assign_30 | OPTIMAL | 30 | 0.001 | OPTIMAL | 30 | 0.008 |
| assign_80 | OPTIMAL | 80 | 0.674 | OPTIMAL | 80 | 0.039 |
| bigmknap_25 | OPTIMAL | 1893.396858 | 0.009 | OPTIMAL | 1893.396859 | 0.029 |
| infeasthin_50 | ITERATION_LIMIT | 50 | 0.000 | INFEASIBLE |  | 0.001 |
| kleeminty_12 | OPTIMAL | 244140625 | 0.000 | OPTIMAL | 244140625 | 0.002 |
| kleeminty_20 | OPTIMAL | 9.536743164e+13 | 0.000 | OPTIMAL | 9.536743164e+13 | 0.002 |
| nearsing_200 | OPTIMAL | -334.6634592 | 0.006 | OPTIMAL | -334.6634592 | 0.009 |
| scaled_adlittle | OPTIMAL | 225494.9632 | 0.001 | OPTIMAL | 187880.5044 | 0.003 |
| scaled_afiro | OPTIMAL | -464.7531429 | 0.000 | OPTIMAL | -458.9245714 | 0.002 |
| scaled_bandm | OPTIMAL | -158.6280185 | 0.008 | INFEASIBLE |  | 0.002 |
| scaled_blend | OPTIMAL | -0.00382226974 | 0.001 | UNBOUNDED |  | 0.003 |
| scaled_e226 | OPTIMAL | -11.63892907 | 0.005 | OPTIMAL | -18.51329414 | 0.007 |
| scaled_sc50a | OPTIMAL | -64.57507706 | 0.001 | OPTIMAL | -86.66666667 | 0.002 |
| scaled_scfxm1 | NUMERICAL_FAILURE | 18416.75903 | 0.005 | INFEASIBLE |  | 0.003 |
| scaled_share1b | OPTIMAL | -76589.31858 | 0.002 | INFEASIBLE |  | 0.002 |
| transpdeg_40 | OPTIMAL | 410 | 0.002 | OPTIMAL | 410 | 0.008 |
| unbnd_30 | UNBOUNDED | 0 | 0.000 | UNBOUNDED |  | 0.002 |

## Failures / disagreements

- infeasthin_50 / auto: ITERATION_LIMIT 
- scaled_scfxm1 / auto: NUMERICAL_FAILURE 
