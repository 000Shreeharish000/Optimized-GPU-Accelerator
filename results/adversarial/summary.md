# Benchmark: adversarial

Instances: 17; time limit 60s; machine: see results/machine.json.
Status = certified status (PRAMANA) / solver status (HiGHS). Time excludes file reading.

| solver | solved (proven) | certified | SGM time (s, shift 10, unsolved at limit) | wrong vs reference |
|---|---|---|---|---|
| auto | 16/17 | 16 | 1.257 | 4 |
| highs | 17/17 | - | 0.013 | - |

## Full per-instance table (including every failure)

| instance | auto status | auto obj | auto time | highs status | highs obj | highs time |
|---|---|---|---|---|---|---|
| assign_30 | OPTIMAL | 30 | 0.001 | OPTIMAL | 30 | 0.008 |
| assign_80 | OPTIMAL | 80 | 0.643 | OPTIMAL | 80 | 0.078 |
| bigmknap_25 | OPTIMAL | 1893.396858 | 0.010 | OPTIMAL | 1893.396859 | 0.029 |
| infeasthin_50 | INFEASIBLE |  | 0.000 | INFEASIBLE |  | 0.003 |
| kleeminty_12 | OPTIMAL | 244140625 | 0.000 | OPTIMAL | 244140625 | 0.002 |
| kleeminty_20 | OPTIMAL | 9.536743164e+13 | 0.001 | OPTIMAL | 9.536743164e+13 | 0.003 |
| nearsing_200 | OPTIMAL | -334.6634592 | 0.008 | OPTIMAL | -334.6634592 | 0.021 |
| scaled_adlittle | OPTIMAL | 225494.9632 | 0.001 | OPTIMAL | 187880.5044 | 0.006 |
| scaled_afiro | OPTIMAL | -464.7531429 | 0.000 | OPTIMAL | -458.9245714 | 0.004 |
| scaled_bandm | OPTIMAL | -158.6280185 | 0.010 | INFEASIBLE |  | 0.006 |
| scaled_blend | OPTIMAL | -30.81214985 | 0.001 | UNBOUNDED |  | 0.005 |
| scaled_e226 | OPTIMAL | -11.63892907 | 0.007 | OPTIMAL | -18.51329414 | 0.015 |
| scaled_sc50a | OPTIMAL | -64.57507706 | 0.001 | OPTIMAL | -86.66666667 | 0.005 |
| scaled_scfxm1 | NUMERICAL_FAILURE | 18416.75903 | 0.007 | INFEASIBLE |  | 0.006 |
| scaled_share1b | OPTIMAL | -76589.31858 | 0.002 | INFEASIBLE |  | 0.004 |
| transpdeg_40 | OPTIMAL | 410 | 0.002 | OPTIMAL | 410 | 0.019 |
| unbnd_30 | UNBOUNDED | 0 | 0.000 | UNBOUNDED |  | 0.003 |

## Failures / disagreements

- scaled_scfxm1 / auto: NUMERICAL_FAILURE 
