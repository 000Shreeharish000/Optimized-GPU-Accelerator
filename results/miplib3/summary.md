# Benchmark: miplib3

Instances: 42; time limit 60s; machine: see results/machine.json.
Status = certified status (PRAMANA) / solver status (HiGHS). Time excludes file reading.

| solver | solved (proven) | certified | SGM time (s, shift 10, unsolved at limit) | wrong vs reference |
|---|---|---|---|---|
| auto | 29/42 | 29 | 14.107 | 0 |
| highs | 36/42 | - | 5.869 | - |

## Full per-instance table (including every failure)

| instance | auto status | auto obj | auto time | highs status | highs obj | highs time |
|---|---|---|---|---|---|---|
| 10teams | TIME_LIMIT |  | 60.057 | OPTIMAL | 924 | 25.749 |
| bell3a | OPTIMAL | 878430.316 | 7.398 | OPTIMAL | 878430.316 | 0.443 |
| bell5 | OPTIMAL | 8966406.492 | 21.693 | OPTIMAL | 8966406.492 | 0.709 |
| blend2 | OPTIMAL | 7.598985 | 3.765 | OPTIMAL | 7.598985 | 3.677 |
| dcmulti | TIME_LIMIT | 188359.5 | 60.014 | OPTIMAL | 188182 | 1.635 |
| egout | OPTIMAL | 568.1007 | 0.013 | OPTIMAL | 568.1007 | 0.010 |
| enigma | OPTIMAL | 0 | 1.202 | OPTIMAL | 0 | 0.545 |
| fixnet6 | OPTIMAL | 3982.999998 | 57.664 | OPTIMAL | 3983 | 4.727 |
| flugpl | OPTIMAL | 1201500 | 0.055 | OPTIMAL | 1201500 | 0.112 |
| gesa2 | OPTIMAL | 25780031.43 | 1.924 | OPTIMAL | 25779856.37 | 0.654 |
| gesa2_o | OPTIMAL | 25779856.37 | 13.700 | OPTIMAL | 25779856.37 | 1.072 |
| gesa3 | OPTIMAL | 27991042.65 | 2.402 | OPTIMAL | 27991042.65 | 2.802 |
| gesa3_o | OPTIMAL | 27991042.65 | 2.877 | OPTIMAL | 27991042.65 | 5.078 |
| gt2 | OPTIMAL | 21166 | 0.092 | OPTIMAL | 21166 | 0.079 |
| khb05250 | OPTIMAL | 106940226 | 0.339 | OPTIMAL | 106940226 | 0.395 |
| lseu | OPTIMAL | 1120 | 0.644 | OPTIMAL | 1120 | 0.407 |
| markshare1 | TIME_LIMIT | 19.99999954 | 60.371 | TIME_LIMIT | 11 | 60.009 |
| markshare2 | TIME_LIMIT | 16 | 60.327 | TIME_LIMIT | 43 | 60.008 |
| mas74 | TIME_LIMIT | 11857.37126 | 60.178 | TIME_LIMIT | 12052.1521 | 60.007 |
| mas76 | TIME_LIMIT | 40005.04626 | 60.302 | TIME_LIMIT | 40005.05414 | 60.013 |
| misc03 | OPTIMAL | 3360 | 2.723 | OPTIMAL | 3360 | 0.648 |
| misc06 | OPTIMAL | 12850.86074 | 0.240 | OPTIMAL | 12851.07629 | 0.143 |
| misc07 | TIME_LIMIT | 2810 | 60.015 | OPTIMAL | 2810 | 27.741 |
| mod008 | OPTIMAL | 307 | 1.027 | OPTIMAL | 307 | 1.487 |
| mod010 | OPTIMAL | 6548 | 0.049 | OPTIMAL | 6548 | 0.356 |
| noswot | TIME_LIMIT | -41 | 60.182 | TIME_LIMIT | -41 | 60.010 |
| p0033 | OPTIMAL | 3089 | 0.029 | OPTIMAL | 3089 | 0.042 |
| p0201 | OPTIMAL | 7615 | 5.392 | OPTIMAL | 7615 | 2.829 |
| p0282 | OPTIMAL | 258411 | 1.750 | OPTIMAL | 258411 | 0.352 |
| p0548 | OPTIMAL | 8691 | 0.372 | OPTIMAL | 8691 | 0.101 |
| p2756 | OPTIMAL | 3124 | 3.303 | OPTIMAL | 3124 | 0.480 |
| pk1 | TIME_LIMIT | 12 | 60.124 | TIME_LIMIT | 14 | 60.014 |
| pp08a | TIME_LIMIT | 7369.999992 | 60.016 | OPTIMAL | 7350 | 1.827 |
| pp08aCUTS | TIME_LIMIT | 7389.999996 | 60.018 | OPTIMAL | 7350 | 1.625 |
| qnet1 | OPTIMAL | 16029.69268 | 50.025 | OPTIMAL | 16029.69268 | 1.684 |
| qnet1_o | OPTIMAL | 16030.99268 | 24.510 | OPTIMAL | 16029.69268 | 1.214 |
| rgn | OPTIMAL | 82.19999913 | 1.531 | OPTIMAL | 82.19999924 | 0.383 |
| set1ch | TIME_LIMIT | 96824.54226 | 60.041 | OPTIMAL | 54542.25 | 0.700 |
| stein27 | OPTIMAL | 18 | 0.912 | OPTIMAL | 18 | 1.057 |
| stein45 | OPTIMAL | 30 | 30.831 | OPTIMAL | 30 | 36.614 |
| vpm1 | OPTIMAL | 20 | 0.381 | OPTIMAL | 20 | 0.023 |
| vpm2 | TIME_LIMIT | 14.5 | 60.043 | OPTIMAL | 13.75 | 2.808 |

## Failures / disagreements

- 10teams / auto: TIME_LIMIT 
- dcmulti / auto: TIME_LIMIT 
- markshare1 / auto: TIME_LIMIT 
- markshare1 / highs: TIME_LIMIT 
- markshare2 / auto: TIME_LIMIT 
- markshare2 / highs: TIME_LIMIT 
- mas74 / auto: TIME_LIMIT 
- mas74 / highs: TIME_LIMIT 
- mas76 / auto: TIME_LIMIT 
- mas76 / highs: TIME_LIMIT 
- misc07 / auto: TIME_LIMIT 
- noswot / auto: TIME_LIMIT 
- noswot / highs: TIME_LIMIT 
- pk1 / auto: TIME_LIMIT 
- pk1 / highs: TIME_LIMIT 
- pp08a / auto: TIME_LIMIT 
- pp08aCUTS / auto: TIME_LIMIT 
- set1ch / auto: TIME_LIMIT 
- vpm2 / auto: TIME_LIMIT 
