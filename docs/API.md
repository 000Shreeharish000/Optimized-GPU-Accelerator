# Interfaces: CLI, Python, C

## CLI (`build/pramana.exe`)

```
pramana [solve] <model.mps[.gz]|model.qps> [options]
pramana verify <model> <result.json>              re-certify a result file (C++ certifier)
pramana parametric <model> (--col N | --row N) --kind cost|rhs|lower|upper --from A --to B [--json F]
pramana family <model> (--col N | --row N) --kind ... --from A --to B --cases K [--json F]
pramana info                                       build + GPU report
pramana calibrate [--out F]                        CPU vs GPU PDHG per-iteration timing
```

| option | meaning | default |
|---|---|---|
| `--algo auto\|dual\|ipm\|pdhg\|pdhg-cpu\|pdhg-gpu\|race` | engine (auto = router) | auto |
| `--time S` | time limit | none |
| `--threads N` | CPU threads (PDHG CPU backend) | all |
| `--no-presolve --no-scale --no-certify --no-crossover --no-gpu` | switches | on |
| `--tol T` | simplex primal/dual feasibility tolerance (scaled space) | 1e-7 |
| `--pdhg-tol T  --pdhg-iters N  --ipm-tol T` | first-order / IPM accuracy | 1e-6 / 1e6 / 1e-8 |
| `--gap G --nodes N --no-cuts --no-heuristics --branching reliability\|pseudocost\|mostfrac\|strong` | MIP | 1e-4 |
| `--mip-threads K` | concurrent diversified branch-and-cut racing on K threads | 1 |
| `--json F [--vectors]` | result + certificate + telemetry (+ named x, duals, basis, rays) | |
| `--cert F` | certificate only | |
| `--router-model F` | calibrated router weights (also auto-loaded from `router_model.json` next to the exe or `PRAMANA_ROUTER_MODEL`) | |
| `--log 0..3` | verbosity | 1 |

Exit codes: 0 proven (optimal/infeasible/unbounded), 2 unproven (limit), 3 numerical failure, 10 error.

Environment: `PRAMANA_DISABLE_GPU=1` hides the GPU; `PRAMANA_NVRTC=<path>` selects an NVRTC DLL;
`PRAMANA_GPU_DEVICE=k` selects the device.

## Python (`python/pramana`)

```python
import pramana
m = pramana.Model("plan", maximize=True)
x = m.add_var(lb=0, ub=10, cost=3, name="x")
y = m.add_var(lb=0, cost=5, integer=True, name="y")
m.add_constraint({x: 1, y: 2}, ub=8, name="cap")
r = m.solve(algorithm="auto", time_limit=60)        # options as in the C API
r.status, r.objective, r.certified, r.row_duals, r.certificate, r.telemetry
pramana.solve_file("model.mps", algorithm="dual")
pramana.parametric("plan.mps", col="BUY_CR03_0", kind="cost", lo=-95, hi=-35)
pramana.family("plan.mps", col="BUY_CR03_0", kind="upper", lo=0, hi=120, cases=64)
```

Independent verification: `python -m pramana.verify model.mps result.json [--exact-basis]`.

## C (`include/pramana/pramana.h`, `build/libpramana.dll`)

`pramana_model_create/read`, `pramana_add_col`, `pramana_add_row`, `pramana_add_q`,
`pramana_solve(model, options_json)`, `pramana_result_status/objective/bound/certified`,
`pramana_result_x/row_duals/reduced_costs/row_activity`, `pramana_result_json`,
`pramana_parametric`, `pramana_family`, `pramana_gpu_info`.
Options JSON keys: `algorithm, time_limit, presolve, certify, scale, threads, log_level, seed,
tolerance, pdhg_tolerance, pdhg_crossover, ipm_tolerance, mip_gap, node_limit, cuts, heuristics,
branching, mip_threads, allow_gpu, router_model`.
