# Benchmark: gen

Instances: 18; time limit 120s; machine: see results/machine.json.
Status = certified status (PRAMANA) / solver status (HiGHS). Time excludes file reading.

| solver | solved (proven) | certified | SGM time (s, shift 10, unsolved at limit) | wrong vs reference |
|---|---|---|---|---|
| dual | 14/18 | 14 | 11.880 | 0 |
| ipm | 14/18 | 14 | 11.867 | 0 |
| pdhg-cpu | 14/18 | 14 | 16.728 | 0 |
| pdhg-gpu | 14/18 | 14 | 13.822 | 0 |
| highs | 17/18 | - | 6.321 | - |

## Full per-instance table (including every failure)

| instance | dual status | dual obj | dual time | ipm status | ipm obj | ipm time | pdhg-cpu status | pdhg-cpu obj | pdhg-cpu time | pdhg-gpu status | pdhg-gpu obj | pdhg-gpu time | highs status | highs obj | highs time |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| blend_12 | OPTIMAL | 66477.85424 | 0.000 | OPTIMAL | 66477.85424 | 0.000 | OPTIMAL | 66477.85424 | 0.055 | OPTIMAL | 66477.85424 | 0.118 | OPTIMAL | 66477.85424 | 0.003 |
| blend_12_mip | OPTIMAL | 67489.42086 | 0.001 | OPTIMAL | 67489.42086 | 0.001 | OPTIMAL | 67489.42086 | 0.001 | OPTIMAL | 67489.42086 | 0.001 | OPTIMAL | 67489.42086 | 0.008 |
| blend_30_mip | OPTIMAL | 66451.92534 | 0.012 | OPTIMAL | 66451.92534 | 0.012 | OPTIMAL | 66451.92534 | 0.012 | OPTIMAL | 66451.92534 | 0.012 | OPTIMAL | 66451.92534 | 0.030 |
| dispatch_50x24 | OPTIMAL | 3524877.432 | 0.010 | OPTIMAL | 3524877.432 | 0.012 | OPTIMAL | 3524877.432 | 0.014 | OPTIMAL | 3524877.432 | 0.011 | N/A |  |  |
| plan_10x12 | OPTIMAL | 42782.89797 | 0.140 | OPTIMAL | 42782.89797 | 0.087 | OPTIMAL | 42782.89797 | 2.000 | OPTIMAL | 42782.89797 | 1.381 | OPTIMAL | 42782.89797 | 0.053 |
| plan_10x12_mip | OPTIMAL | 40545.31158 | 21.022 | OPTIMAL | 40545.31158 | 21.118 | OPTIMAL | 40545.31158 | 21.068 | OPTIMAL | 40545.31158 | 21.141 | OPTIMAL | 40545.31158 | 24.019 |
| plan_10x6_mip | OPTIMAL | 21815.62546 | 3.804 | OPTIMAL | 21815.62546 | 3.960 | OPTIMAL | 21815.62546 | 3.797 | OPTIMAL | 21815.62546 | 3.815 | OPTIMAL | 21815.62545 | 2.635 |
| plan_20x12 | OPTIMAL | 149376.2855 | 0.464 | OPTIMAL | 149376.2855 | 0.290 | OPTIMAL | 149376.2855 | 5.668 | OPTIMAL | 149376.2855 | 2.869 | OPTIMAL | 149376.2855 | 0.134 |
| plan_20x52 | OPTIMAL | 524332.4409 | 9.613 | OPTIMAL | 524332.4409 | 9.316 | OPTIMAL | 524332.4409 | 85.252 | OPTIMAL | 524332.4409 | 29.514 | OPTIMAL | 524332.4409 | 1.352 |
| plan_30x52 | OPTIMAL | 933661.1098 | 19.783 | OPTIMAL | 933661.1098 | 20.291 | OPTIMAL | 933661.1098 | 97.684 | OPTIMAL | 933661.1098 | 33.982 | OPTIMAL | 933661.1098 | 3.885 |
| plan_40x104 | TIME_LIMIT | 2353429.389 | 120.013 | TIME_LIMIT | 1626719.665 | 121.753 | TIME_LIMIT | 2160802.613 | 120.086 | NUMERICAL_FAILURE | 2160639.312 | 120.084 | OPTIMAL | 2160633.512 | 26.646 |
| plan_60x156 | TIME_LIMIT | 24948143.41 | 120.044 | TIME_LIMIT | -3708873.131 | 137.319 | TIME_LIMIT | 5410048.211 | 120.146 | NUMERICAL_FAILURE | 5410116.994 | 120.196 | OPTIMAL | 5410126.647 | 91.764 |
| plan_6x4 | OPTIMAL | 14856.53514 | 0.007 | OPTIMAL | 14856.53514 | 0.010 | OPTIMAL | 14856.53514 | 1.781 | OPTIMAL | 14856.53514 | 1.127 | OPTIMAL | 14856.53514 | 0.008 |
| plan_6x4_mip | OPTIMAL | 13985.02229 | 0.572 | OPTIMAL | 13985.02229 | 0.554 | OPTIMAL | 13985.02229 | 0.556 | OPTIMAL | 13985.02229 | 0.550 | OPTIMAL | 13985.02229 | 0.429 |
| uc_10x24 | OPTIMAL | 644229.2463 | 0.046 | OPTIMAL | 644229.2463 | 0.046 | OPTIMAL | 644229.2463 | 0.046 | OPTIMAL | 644229.2463 | 0.046 | OPTIMAL | 644229.2463 | 0.712 |
| uc_30x24 | TIME_LIMIT | 2064403.136 | 120.010 | TIME_LIMIT | 2064403.136 | 120.269 | TIME_LIMIT | 2064403.136 | 120.039 | TIME_LIMIT | 2064403.136 | 120.075 | OPTIMAL | 1933392.047 | 7.358 |
| unload_bigm | TIME_LIMIT | 259.4911452 | 120.364 | TIME_LIMIT | 259.4911452 | 120.334 | TIME_LIMIT | 259.4911452 | 120.304 | TIME_LIMIT | 259.4911452 | 120.340 | OPTIMAL | 109.4911452 | 0.167 |
| unload_tight | OPTIMAL | 119.4911449 | 6.353 | OPTIMAL | 119.4911449 | 6.296 | OPTIMAL | 119.4911449 | 6.302 | OPTIMAL | 119.4911449 | 6.313 | OPTIMAL | 119.4911452 | 0.190 |

## Failures / disagreements

- dispatch_50x24 / highs: N/A QP not supported by scipy.optimize
- plan_40x104 / dual: TIME_LIMIT 
- plan_40x104 / ipm: TIME_LIMIT 
- plan_40x104 / pdhg-cpu: TIME_LIMIT 
- plan_40x104 / pdhg-gpu: NUMERICAL_FAILURE 
- plan_60x156 / dual: TIME_LIMIT 
- plan_60x156 / ipm: TIME_LIMIT 
- plan_60x156 / pdhg-cpu: TIME_LIMIT 
- plan_60x156 / pdhg-gpu: NUMERICAL_FAILURE 
- uc_30x24 / dual: TIME_LIMIT 
- uc_30x24 / ipm: TIME_LIMIT 
- uc_30x24 / pdhg-cpu: TIME_LIMIT 
- uc_30x24 / pdhg-gpu: TIME_LIMIT 
- unload_bigm / dual: TIME_LIMIT 
- unload_bigm / ipm: TIME_LIMIT 
- unload_bigm / pdhg-cpu: TIME_LIMIT 
- unload_bigm / pdhg-gpu: TIME_LIMIT 
