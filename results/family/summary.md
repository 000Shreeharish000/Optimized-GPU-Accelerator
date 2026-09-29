# Repeated-solve (crude valuation) experiment

Exact = certified parametric sweep. Errors are max relative errors vs the exact value function.

## plan_10x12: crude CR03 price (period 0), $/bbl, negated (max model)

- parameter `BUY_CR03_0` (cost) in [-95.0, -35.0], K = 64 sampled cases
- exact breakpoints: 3 (2 segments contain no sample -> invisible to any sampling strategy)
- breakpoints: -68.2646, -68.1292, -68.1193

| strategy | seconds | solved | certified | max rel. error | note |
|---|---|---|---|---|---|
| parametric sweep (exact, certified) | 0.140 | 4 | 4 | 0.00e+00 | 4 segments, 3190 pivots; certificates per segment |
| cold dual simplex x K | 10.048 | 64 | 64 | 8.51e-16 | independent solves (incl. certification) |
| warm-started simplex chain | 0.190 | 64 | 64 | 3.39e-16 | 3169 pivots total |
| batched PDHG (CPU) 1e-6 + safe bounds | 18.974 | 64 | 64 | 2.73e-06 | 16 threads, 18560 iterations, kernel 17.313s, transfer 0.039s, setup 1.659s; certified = safe bound within 1e-4 |
| batched PDHG (GPU) 1e-6 + safe bounds | 3.787 | 64 | 64 | 2.73e-06 | NVIDIA GeForce RTX 4050 Laptop GPU, 18560 iterations, kernel 3.759s, transfer 0.076s, setup 0.024s; certified = safe bound within 1e-4 |

## plan_10x12: crude CR03 availability (period 0), kbbl

- parameter `BUY_CR03_0` (upper) in [0.0, 120.0], K = 64 sampled cases
- exact breakpoints: 18 (2 segments contain no sample -> invisible to any sampling strategy)
- breakpoints: 2.99514, 3.34439, 49.6512, 51.7112, 54.7239, 56.6473, 62.181, 73.1157, 86.4339, 87.2694, 87.9345, 90.9976 ...

| strategy | seconds | solved | certified | max rel. error | note |
|---|---|---|---|---|---|
| parametric sweep (exact, certified) | 0.163 | 19 | 19 | 0.00e+00 | 19 segments, 3291 pivots; certificates per segment |
| cold dual simplex x K | 9.907 | 64 | 64 | 1.36e-15 | independent solves (incl. certification) |
| warm-started simplex chain | 0.191 | 64 | 64 | 1.28e-14 | 3219 pivots total |
| batched PDHG (CPU) 1e-6 + safe bounds | 24.632 | 64 | 64 | 2.84e-06 | 16 threads, 23488 iterations, kernel 22.966s, transfer 0.042s, setup 1.664s; certified = safe bound within 1e-4 |
| batched PDHG (GPU) 1e-6 + safe bounds | 4.823 | 64 | 64 | 2.84e-06 | NVIDIA GeForce RTX 4050 Laptop GPU, 23488 iterations, kernel 4.792s, transfer 0.073s, setup 0.027s; certified = safe bound within 1e-4 |

## plan_20x12: crude CR07 price (period 3)

- parameter `BUY_CR07_3` (cost) in [-100.0, -30.0], K = 64 sampled cases
- exact breakpoints: 22 (21 segments contain no sample -> invisible to any sampling strategy)
- breakpoints: -69.7383, -69.7373, -69.7324, -69.7274, -69.7197, -69.7121, -69.7055, -69.7011, -69.695, -69.6817, -69.6348, -69.1328 ...

| strategy | seconds | solved | certified | max rel. error | note |
|---|---|---|---|---|---|
| parametric sweep (exact, certified) | 0.637 | 23 | 23 | 0.00e+00 | 23 segments, 7139 pivots; certificates per segment |
| cold dual simplex x K | 38.230 | 64 | 64 | 9.34e-16 | independent solves (incl. certification) |
| warm-started simplex chain | 0.688 | 64 | 64 | 5.84e-16 | 6375 pivots total |
| batched PDHG (CPU) 1e-6 + safe bounds | 80.931 | 64 | 64 | 2.58e-06 | 16 threads, 49920 iterations, kernel 80.310s, transfer 0.067s, setup 0.619s; certified = safe bound within 1e-4 |
| batched PDHG (GPU) 1e-6 + safe bounds | 19.419 | 64 | 64 | 2.58e-06 | NVIDIA GeForce RTX 4050 Laptop GPU, 49920 iterations, kernel 19.367s, transfer 0.149s, setup 0.049s; certified = safe bound within 1e-4 |

## plan_20x52: crude CR05 availability (week 10)

- parameter `BUY_CR05_10` (upper) in [0.0, 200.0], K = 32 sampled cases
- exact breakpoints: 3 (0 segments contain no sample -> invisible to any sampling strategy)
- breakpoints: 25.5207, 26.3146, 36.812

| strategy | seconds | solved | certified | max rel. error | note |
|---|---|---|---|---|---|
| parametric sweep (exact, certified) | 8.825 | 4 | 4 | 0.00e+00 | 4 segments, 28059 pivots; certificates per segment |
| cold dual simplex x K | 272.982 | 32 | 32 | 2.61e-09 | independent solves (incl. certification) |
| warm-started simplex chain | 9.434 | 32 | 32 | 2.61e-09 | 28026 pivots total |
| batched PDHG (CPU) 1e-6 + safe bounds | 803.309 | 0 | 0 | 0.00e+00 | 16 threads, 200000 iterations, kernel 801.504s, transfer 0.104s, setup 1.737s; certified = safe bound within 1e-4 |
| batched PDHG (GPU) 1e-6 + safe bounds | 167.071 | 0 | 0 | 0.00e+00 | NVIDIA GeForce RTX 4050 Laptop GPU, 200000 iterations, kernel 166.876s, transfer 0.158s, setup 0.097s; certified = safe bound within 1e-4 |

