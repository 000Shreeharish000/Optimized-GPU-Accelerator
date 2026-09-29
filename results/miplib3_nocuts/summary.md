# Benchmark: miplib3

Instances: 42; time limit 60s; machine: see results/machine.json.
Status = certified status (PRAMANA) / solver status (HiGHS). Time excludes file reading.

| solver | solved (proven) | certified | SGM time (s, shift 10, unsolved at limit) | wrong vs reference |
|---|---|---|---|---|
| auto | 25/42 | 25 | 14.807 | 0 |

## Full per-instance table (including every failure)

| instance | auto status | auto obj | auto time |
|---|---|---|---|
| 10teams | TIME_LIMIT |  | 60.017 |
| bell3a | OPTIMAL | 878430.316 | 2.864 |
| bell5 | TIME_LIMIT | 8966768.209 | 60.826 |
| blend2 | OPTIMAL | 7.598985 | 1.066 |
| dcmulti | OPTIMAL | 188182 | 0.544 |
| egout | OPTIMAL | 568.1007 | 0.105 |
| enigma | OPTIMAL | 0 | 0.129 |
| fixnet6 | TIME_LIMIT | 6472.999981 | 60.168 |
| flugpl | OPTIMAL | 1201500 | 0.014 |
| gesa2 | TIME_LIMIT | 26504664.43 | 60.170 |
| gesa2_o | TIME_LIMIT |  | 60.171 |
| gesa3 | OPTIMAL | 27991042.65 | 6.631 |
| gesa3_o | OPTIMAL | 27991042.65 | 10.363 |
| gt2 | OPTIMAL | 21166 | 0.143 |
| khb05250 | OPTIMAL | 106940226 | 0.381 |
| lseu | OPTIMAL | 1120 | 1.061 |
| markshare1 | TIME_LIMIT | 17.99999911 | 60.866 |
| markshare2 | TIME_LIMIT | 28 | 60.831 |
| mas74 | TIME_LIMIT | 11801.17983 | 60.265 |
| mas76 | TIME_LIMIT | 40005.05414 | 60.302 |
| misc03 | OPTIMAL | 3360 | 0.451 |
| misc06 | OPTIMAL | 12850.86074 | 0.093 |
| misc07 | TIME_LIMIT | 2810 | 60.023 |
| mod008 | OPTIMAL | 307 | 1.069 |
| mod010 | OPTIMAL | 6548 | 0.779 |
| noswot | TIME_LIMIT | -41 | 60.225 |
| p0033 | OPTIMAL | 3089 | 0.009 |
| p0201 | OPTIMAL | 7615 | 0.327 |
| p0282 | OPTIMAL | 258411 | 0.130 |
| p0548 | OPTIMAL | 8691 | 9.289 |
| p2756 | TIME_LIMIT |  | 60.044 |
| pk1 | TIME_LIMIT | 11 | 60.288 |
| pp08a | TIME_LIMIT | 7879.999991 | 60.329 |
| pp08aCUTS | TIME_LIMIT | 7619.99999 | 60.175 |
| qnet1 | OPTIMAL | 16029.69268 | 1.260 |
| qnet1_o | OPTIMAL | 16029.69268 | 1.189 |
| rgn | OPTIMAL | 82.19999924 | 0.242 |
| set1ch | TIME_LIMIT | 83382.5 | 60.203 |
| stein27 | OPTIMAL | 18 | 0.291 |
| stein45 | OPTIMAL | 30 | 17.010 |
| vpm1 | OPTIMAL | 20 | 19.028 |
| vpm2 | TIME_LIMIT | 14.25 | 60.431 |

## Failures / disagreements

- 10teams / auto: TIME_LIMIT 
- bell5 / auto: TIME_LIMIT 
- fixnet6 / auto: TIME_LIMIT 
- gesa2 / auto: TIME_LIMIT 
- gesa2_o / auto: TIME_LIMIT 
- markshare1 / auto: TIME_LIMIT 
- markshare2 / auto: TIME_LIMIT 
- mas74 / auto: TIME_LIMIT 
- mas76 / auto: TIME_LIMIT 
- misc07 / auto: TIME_LIMIT 
- noswot / auto: TIME_LIMIT 
- p2756 / auto: TIME_LIMIT 
- pk1 / auto: TIME_LIMIT 
- pp08a / auto: TIME_LIMIT 
- pp08aCUTS / auto: TIME_LIMIT 
- set1ch / auto: TIME_LIMIT 
- vpm2 / auto: TIME_LIMIT 
