# Router validation (5-fold cross-validation, held-out instances)

Instances: 99 LPs (91 Netlib, 8 refinery planning). Engines: dual, ipm, pdhg-cpu, pdhg-gpu.

| policy | geo-mean regret (time / oracle time) | worst regret | picks the fastest engine |
|---|---|---|---|
| **PRAMANA router (held-out)** | 1.062 | 11.3 | 89% |
| always dual | 1.059 | 11.3 | 90% |
| always ipm | 1.796 | 11.7 | 12% |
| always pdhg-cpu | 47.685 | 4278.5 | 2% |
| always pdhg-gpu | 113.305 | 18498.7 | 2% |

Oracle (fastest certified engine) counts: dual: 89, ipm: 10, pdhg-cpu: 0, pdhg-gpu: 0

