# 7-minute demo script

Run from the repository root after `build.bat` and `python bench/fetch_data.py`.

| t | command | what it proves |
|---|---|---|
| 0:00 | `python bench/audit_deps.py` | binaries import only the OS + C++ runtime; no solver library; 11.6k lines of own code |
| 0:40 | `build\pramana.exe data\netlib\pilot87.mps.gz --algo dual --log 2` | a notoriously ill-conditioned Netlib LP solved from scratch, certificate table printed (gap ~1e-13) |
| 1:30 | `build\pramana.exe data\netlib\pilot87.mps.gz --json r.json --vectors --log 0` then `python -m pramana.verify data\netlib\pilot87.mps.gz r.json` (from `python/`) | an independent program with its own parser re-checks the claim in exact rational arithmetic |
| 2:10 | `build\pramana.exe data\adversarial\scaled_afiro.mps` | badly scaled data (10^±6): PRAMANA returns the true optimum −464.753 with a rigorous certificate (HiGHS via scipy returns −458.92) |
| 2:50 | `build\pramana.exe data\netlib_infeas\klein3.mps` | infeasibility with a verified Farkas certificate |
| 3:20 | `build\pramana.exe data\gen\plan_20x52.mps --algo race` then `--algo pdhg-gpu --log 3` | the GPU is real (own kernels, NVRTC), transfers are counted, the router explains its choice |
| 4:20 | `build\pramana.exe calibrate` and `docs/figures/gpu_iteration_crossover.png` | where the GPU starts to win per iteration, and where it loses |
| 5:00 | `build\pramana.exe family data\gen\plan_10x12.mps --col BUY_CR03_0 --kind upper --from 0 --to 120 --cases 32` | crude valuation: exact certified value curve in 0.16 s; uniform sampling misses 7 of 18 segments; warm simplex vs batched GPU PDHG measured |
| 5:50 | `build\pramana.exe data\gen\uc_10x24.mps` and `--no-cuts` | MILP: cuts close the root gap (0.09 s vs 29.9 s) |
| 6:20 | open `web/index.html`, drop `r.json` | certificate / telemetry / MIP progress in the browser |
| 6:40 | `docs/BENCHMARKS.md` | full tables vs HiGHS including every failure |
