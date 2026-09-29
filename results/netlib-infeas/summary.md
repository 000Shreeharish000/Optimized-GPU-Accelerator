# Benchmark: netlib-infeas

Instances: 29; time limit 300s; machine: see results/machine.json.
Status = certified status (PRAMANA) / solver status (HiGHS). Time excludes file reading.

| solver | solved (proven) | certified | SGM time (s, shift 10, unsolved at limit) | wrong vs reference |
|---|---|---|---|---|
| dual | 29/29 | 29 | 0.239 | 0 |
| highs | 28/29 | - | 1.279 | - |

## Full per-instance table (including every failure)

| instance | dual status | dual obj | dual time | highs status | highs obj | highs time |
|---|---|---|---|---|---|---|
| bgdbg1 | INFEASIBLE | 1.92 | 0.001 | INFEASIBLE |  | 0.004 |
| bgetam | INFEASIBLE | -762.9410132 | 0.002 | INFEASIBLE |  | 0.005 |
| bgindy | INFEASIBLE | 9534768.342 | 0.017 | INFEASIBLE |  | 0.036 |
| bgprtr | INFEASIBLE | 7472389.091 | 0.000 | INFEASIBLE |  | 0.002 |
| box1 | INFEASIBLE | 268 | 0.001 | INFEASIBLE |  | 0.002 |
| ceria3d | INFEASIBLE | -0.953125 | 0.508 | INFEASIBLE |  | 0.018 |
| chemcom | INFEASIBLE | 38955.15347 | 0.001 | INFEASIBLE |  | 0.010 |
| cplex1 | INFEASIBLE | -522642.6188 | 0.087 | INFEASIBLE |  | 0.158 |
| cplex2 | INFEASIBLE | 0.6570600371 | 0.005 | INFEASIBLE |  | 0.011 |
| ex72a | INFEASIBLE | 356 | 0.001 | INFEASIBLE |  | 0.004 |
| ex73a | INFEASIBLE | 305 | 0.001 | INFEASIBLE |  | 0.003 |
| forest6 | INFEASIBLE | 1388164.069 | 0.001 | INFEASIBLE |  | 0.004 |
| galenet | INFEASIBLE | 0 | 0.000 | INFEASIBLE |  | 0.003 |
| gosh | INFEASIBLE | 19.14195053 | 7.573 | INFEASIBLE |  | 0.109 |
| gran | INFEASIBLE | -6458.706985 | 0.072 | INFEASIBLE |  | 0.027 |
| greenbea | INFEASIBLE | 2531.558058 | 0.339 | INFEASIBLE |  | 0.040 |
| itest2 | INFEASIBLE | 0 | 0.000 | INFEASIBLE |  | 0.003 |
| itest6 | INFEASIBLE | 1747000 | 0.000 | INFEASIBLE |  | 0.003 |
| klein1 | INFEASIBLE | 0 | 0.001 | INFEASIBLE |  | 0.004 |
| klein2 | INFEASIBLE | 0 | 0.025 | INFEASIBLE |  | 0.015 |
| klein3 | INFEASIBLE | 0 | 0.109 | STATUS4 |  | 0.056 |
| mondou2 | INFEASIBLE | 353200553 | 0.003 | INFEASIBLE |  | 0.006 |
| pang | INFEASIBLE | 65151.37709 | 0.013 | INFEASIBLE |  | 0.019 |
| pilot4i | INFEASIBLE | -141.0246005 | 0.016 | INFEASIBLE |  | 0.010 |
| qual | INFEASIBLE | 389040570.5 | 0.012 | INFEASIBLE |  | 0.015 |
| reactor | INFEASIBLE | -365899.7398 | 0.001 | INFEASIBLE |  | 0.006 |
| refinery | INFEASIBLE | 35536237.18 | 0.009 | INFEASIBLE |  | 0.013 |
| vol1 | INFEASIBLE | 836482.7736 | 0.019 | INFEASIBLE |  | 0.021 |
| woodinfe | INFEASIBLE | 23720 | 0.000 | INFEASIBLE |  | 0.003 |

## Failures / disagreements

- klein3 / highs: STATUS4 
