# Benchmark: netlib

Instances: 3; time limit 60s; machine: see results/machine.json.
Status = certified status (PRAMANA) / solver status (HiGHS). Time excludes file reading.

| solver | solved (proven) | certified | SGM time (s, shift 10, unsolved at limit) | wrong vs reference |
|---|---|---|---|---|
| dual | 3/3 | 3 | 0.368 | 0 |
| highs | 3/3 | - | 0.073 | - |

## Full per-instance table (including every failure)

| instance | dual status | dual obj | dual time | highs status | highs obj | highs time |
|---|---|---|---|---|---|---|
| 25fv47 | OPTIMAL | 5501.845888 | 0.201 | OPTIMAL | 5501.845888 | 0.113 |
| 80bau3b | OPTIMAL | 987224.1924 | 0.737 | OPTIMAL | 987224.1924 | 0.103 |
| adlittle | OPTIMAL | 225494.9632 | 0.176 | OPTIMAL | 225494.9632 | 0.003 |

## Failures / disagreements

- none
