# Benchmark: netlib-infeas

Instances: 29; time limit 300s; machine: see results/machine.json.
Status = certified status (PRAMANA) / solver status (HiGHS). Time excludes file reading.

| solver | solved (proven) | certified | SGM time (s, shift 10, unsolved at limit) | wrong vs reference |
|---|---|---|---|---|
| dual | 29/29 | 29 | 0.174 | 0 |
| highs | 28/29 | - | 1.270 | - |

## Full per-instance table (including every failure)

| instance | dual status | dual obj | dual time | highs status | highs obj | highs time |
|---|---|---|---|---|---|---|
| bgdbg1 | INFEASIBLE | 1.92 | 0.001 | INFEASIBLE |  | 0.003 |
| bgetam | INFEASIBLE | -762.9410132 | 0.001 | INFEASIBLE |  | 0.004 |
| bgindy | INFEASIBLE | 9534768.342 | 0.015 | INFEASIBLE |  | 0.031 |
| bgprtr | INFEASIBLE | 7472389.091 | 0.000 | INFEASIBLE |  | 0.002 |
| box1 | INFEASIBLE | 268 | 0.001 | INFEASIBLE |  | 0.002 |
| ceria3d | INFEASIBLE | -0.9583333333 | 0.214 | INFEASIBLE |  | 0.011 |
| chemcom | INFEASIBLE | 38955.15347 | 0.001 | INFEASIBLE |  | 0.006 |
| cplex1 | INFEASIBLE | -477351.3279 | 0.090 | INFEASIBLE |  | 0.103 |
| cplex2 | INFEASIBLE | 0.6570600371 | 0.004 | INFEASIBLE |  | 0.008 |
| ex72a | INFEASIBLE | 356 | 0.001 | INFEASIBLE |  | 0.002 |
| ex73a | INFEASIBLE | 305 | 0.001 | INFEASIBLE |  | 0.002 |
| forest6 | INFEASIBLE | 1388164.069 | 0.001 | INFEASIBLE |  | 0.002 |
| galenet | INFEASIBLE | 0 | 0.000 | INFEASIBLE |  | 0.002 |
| gosh | INFEASIBLE | 22.08935954 | 4.974 | INFEASIBLE |  | 0.068 |
| gran | INFEASIBLE | -6458.706985 | 0.069 | INFEASIBLE |  | 0.012 |
| greenbea | INFEASIBLE | -7240.661643 | 0.409 | INFEASIBLE |  | 0.018 |
| itest2 | INFEASIBLE | 0 | 0.000 | INFEASIBLE |  | 0.002 |
| itest6 | INFEASIBLE | 1747000 | 0.000 | INFEASIBLE |  | 0.002 |
| klein1 | INFEASIBLE | 0 | 0.001 | INFEASIBLE |  | 0.003 |
| klein2 | INFEASIBLE | 0 | 0.020 | INFEASIBLE |  | 0.011 |
| klein3 | INFEASIBLE | 0 | 0.096 | STATUS4 |  | 0.032 |
| mondou2 | INFEASIBLE | 353200553 | 0.003 | INFEASIBLE |  | 0.003 |
| pang | INFEASIBLE | 65151.37709 | 0.011 | INFEASIBLE |  | 0.008 |
| pilot4i | INFEASIBLE | -141.0246005 | 0.015 | INFEASIBLE |  | 0.005 |
| qual | INFEASIBLE | 389040570.5 | 0.011 | INFEASIBLE |  | 0.009 |
| reactor | INFEASIBLE | -365899.7398 | 0.001 | INFEASIBLE |  | 0.003 |
| refinery | INFEASIBLE | 35536237.18 | 0.008 | INFEASIBLE |  | 0.007 |
| vol1 | INFEASIBLE | 836482.7736 | 0.014 | INFEASIBLE |  | 0.009 |
| woodinfe | INFEASIBLE | 23720 | 0.000 | INFEASIBLE |  | 0.002 |

## Failures / disagreements

- klein3 / highs: STATUS4 
