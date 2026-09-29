# SIH26119 — Indigenous GPU-Accelerated Optimization Solver
## Adversarial Research & Solution-Design Study

*Prepared 29 Sep 2026. Evidence labels used throughout:*
**[F]** documented fact (source cited) · **[R]** published research finding · **[I]** my inference · **[Rec]** my recommendation.

> **Read this first — three things changed the strategy during research**
> 1. **You are not the only team, and at least one competitor is far along.** Public GitHub repos for SIH26119 already exist: SANKHYA (C++20; primal/dual simplex, Mehrotra IPM with own LDLᵀ, CUDA PDHG, B&B with Gomory/cover/MIR root cuts, presolve, independent verifier, claims Netlib 81/89 at 1e-6 and MIPLIB-2017 14/30 easy instances, GPU PDHG 4–10× on L4/A100), PIPEPYE (dual simplex + Markowitz LU, CUDA SpMV PDHG, 174 tests, "hardware router"), YUKTI (dense simplex, 3-layer verifier), INDRA (simplex/IPM/PDHG, Flask dashboard, no GPU, no benchmarks), NOMOS and several others ([search](https://github.com/thegoodengineers/SANKHYA), [PIPEPYE](https://github.com/Satyanshgaur/PIPEPYE), [YUKTI](https://github.com/Deepak-Kambala/YUKTI), [INDRA](https://github.com/ashrayr03/INDRA-SIH2026)). SANKHYA's issue tracker even has "safe dual bounds from any approximate dual" and "batched PDHG bounding of open nodes" ([#519](https://github.com/thegoodengineers/SANKHYA/issues/519), [PR #682](https://github.com/thegoodengineers/SANKHYA/pull/682)). **"Simplex + PDHG-on-GPU + B&B + verifier" is now the *baseline* for this PS, not the differentiator.** (Claims above are the repos' own; I did not re-run them.)
> 2. **The GPU-for-LP/MIP frontier moved a lot in 2025–26.** GPU PDHG is in Gurobi 13 (GPU fully supported from 13.0.2), Xpress 9.8 (beta), COPT, HiGHS and NVIDIA cuOpt; batched PDLP for strong branching shipped in cuOpt 26.02; batched GPU B&B for SOS2 problems (B³-PWL, Aug 2026); GPU presolve (cuPSLP, Sep 2026). Most "GPU + MIP" novelty ideas are now Category A/B.
> 3. **SIH 2026 internal idea submission closes ~30 Sep 2026; the 36-hour finale is in December 2026** ([timeline](https://techpathdaily.com/sih-2026-official-timeline-is-here/)). You have roughly **9–10 weeks** of real build time.

---

## 1. Executive Interpretation of the PS

MRPL is not asking for "an app that solves LPs." Stripped of framing, the PS asks: **can a small Indian team build — from mathematical foundations — a solver core whose *correctness* is trustworthy on the ugly models refineries actually produce, whose internals are inspectable and modifiable, and which uses modern hardware (GPU) honestly?**

The title name-checks *Xpress* and *CPLEX* (misspelt "Express/CEPLEX") rather than Gurobi. **[F]** Aspen PIMS — the dominant refinery planning LP system — runs on CPLEX or Xpress underneath, and Honeywell RPMS licenses Xpress ([patent text](https://patents.google.com/patent/US20030097243A1/en)). **[I]** MRPL's real pain point is very likely the solver embedded in its PIMS-type planning stack: medium-size (10⁴–10⁵ rows), highly degenerate, badly scaled LPs/MILPs solved **repeatedly** (successive-LP "distributive recursion" passes, crude-valuation case runs, monthly/rolling plans). That workload profile matters enormously for where a GPU helps (mostly: it doesn't, per-solve — see §12).

**The honest one-line thesis a strong team can defend:** *"We built the solver core ourselves; every answer it gives is independently certified; we measured exactly where the GPU helps and route work accordingly; and we targeted the repeated-solve structure of refinery planning."*

## 2. What MRPL Actually Wants

### A. Explicit asks
LP, MILP, QP; extensible to MIQP/NLP/MINLP; revised simplex / IPM; B&B / B&C / cuts / presolve / heuristics / node selection; sparse techniques; multi-core; GPU "where it provides measurable benefits"; numerical stability; MIPLIB/Netlib/Mittelmann/QPLIB benchmarks compared against ≥1 established solver; demonstrate degeneracy / weak relaxations / ill-conditioning; **not built on any existing open-source solver library**; API/CLI sufficient.

### B. Implicit expectations
| Implicit expectation | Evidence in wording |
|---|---|
| Correctness beats speed | "numerically robust", "reliable convergence", "consistently finds high-quality solutions" |
| Honest GPU usage — GPU must *earn* its place | "GPU acceleration **considered where it provides measurable benefits**" — this is a trap for "GPU everywhere" teams |
| Ownership = ability to inspect/modify internals | "cannot inspect, modify or tailor the solver internals" |
| Head-to-head comparison, not self-reported numbers | "compared against at least one established commercial or open-source solver" |
| Failure honesty | "where weaker implementations exhibit excessive computation times or fail to converge" |
| Industrial realism | named sectors; dataset note asks for literature case studies |
| Modular architecture | "modular architecture that can later be extended" |

### C. Mandatory vs aspirational
- **Mandatory:** LP + MILP + QP working on real benchmark files; from-scratch core; comparison vs a reference solver; robustness demonstration; API/CLI.
- **Aspirational/future:** MIQP/NLP/MINLP; "thousands to millions of variables" consistently; commercial-grade B&C; GPU wherever beneficial.

### D. Scoring-signal words
"numerically robust", "measurable benefits", "from scratch / mathematical foundation", "compared against", "degeneracy / ill-conditioned / weak LP relaxations", "transparent, extensible". **[I]** Judges will score *trustworthiness of claims* and *depth of the core*, and will discount UI.

### E–G. Difficulty triage
| Requirement | Difficulty | Hackathon realism |
|---|---|---|
| Sparse dual simplex robust on Netlib | Very hard (decades of engineering: LU updates, Harris ratio test, perturbation, DSE) | Realistic to reach ~80–90% of Netlib with 8–10 weeks of focused work |
| IPM with own sparse LDLᵀ + ordering | Hard | Realistic (medium-size) |
| GPU PDHG LP | Moderate (algorithm is simple; tuning isn't) | Very realistic — hence everyone does it |
| MILP B&B with pseudocost/reliability branching, root cuts, basic heuristics | Hard | Realistic for "easy" MIPLIB subset |
| Competitive with Gurobi/CPLEX on MIPLIB | Extremely hard (person-centuries) | **Impossible** — do not attempt, do not claim |
| "Millions of variables" MILP | Extremely hard | Impossible except LP via PDHG (to 1e-4 rel.) |
| MIQP/NLP/MINLP | Hard | Future work (design hooks only) |

### H–J. Weak vs serious vs exceptional team
- **Weak:** Python simplex on dense tableau (or secretly `scipy.optimize.linprog` = HiGHS), toy refinery LP, React dashboard, "100× GPU speedup" on a dense matrix multiply, no benchmark files, "AI-powered".
- **Serious:** C++ sparse simplex + B&B, MPS reader, Netlib/MIPLIB subset results vs HiGHS, CUDA PDHG with CPU-vs-GPU crossover plot, verifier. **(This is where SANKHYA/PIPEPYE already are.)**
- **Exceptional:** All of the above *plus* a sharp, measured contribution nobody else has: e.g., **certified** outputs (safe bounds / rational re-verification), a workload cost model that predicts CPU-vs-GPU choice and is validated, and a refinery repeated-solve story (warm-start + batch) with numbers, plus a published "failure table".

---

## 3. Complete Technical Problem Decomposition

Format per component: **Problem → why hard → conventional/commercial/open → naïve outcome → minimum credible → possible differentiator.**

| Component | Core problem & why hard | Conventional / commercial / OSS practice | Naïve approach outcome | Minimum credible | Differentiator potential |
|---|---|---|---|---|---|
| **Model representation** | Need bounds on every var/row (ranges, free, fixed), objective sense/offset, integrality, Q | All solvers: column-major matrix + lower/upper vectors for columns *and* rows (row activity as logical slack) | Converting everything to standard form `Ax=b, x≥0` → doubles variables, destroys bound structure, breaks bound-flipping | Internal form `min cᵀx + ½xᵀQx, L ≤ Ax ≤ U, l ≤ x ≤ u` | Low |
| **Parsing** | Free/fixed MPS dialects, RANGES semantics (sign-dependent!), negative upper bounds, MARKER INTORG, QMATRIX/QUADOBJ | Every solver has its own; CPLEX/Gurobi read MPS/LP | Wrong RANGES/BOUNDS semantics silently changes the problem → "wrong optimum" on Netlib | Full free-MPS + QPS reader; tested against published Netlib optima | Low |
| **Sparse storage** | Simplex needs column access (pricing rows need row access too); IPM needs AAᵀ/KKT; GPU needs CSR | CSC + row-wise copy (HiGHS, Clp); CSR on GPU | Dense matrices → Netlib ok, anything ≥10k rows dies | CSC master + CSR mirror, both 32-bit indices, 64-bit values | Low |
| **Presolve** | Reductions must be postsolvable exactly (dual values too); interactions; numerics of substitution | Removing empties/singletons, doubleton eqs, dominated cols, duplicate rows/cols, bound tightening, coefficient tightening, probing (MIP). Commercial presolve is a large fraction of MIP speed [R: Achterberg & Wunderling 2013] | Presolve without dual postsolve → no duals/reduced costs; wrong postsolve → infeasible "optimal" solutions | Singleton rows/cols, empty, fixed, forcing/redundant rows, bound tightening, dual fixing, with full primal+dual postsolve stack | Medium (lightweight GPU-oriented presolve [R: Cederberg & Boyd 2026]) |
| **Scaling** | Coefficients spanning 1e-6…1e6 cause tolerance meaninglessness | Geometric-mean + equilibration (Clp/HiGHS); Ruiz + Pock–Chambolle for PDHG (PDLP) | No scaling → simplex stalls, PDHG diverges/crawls | Geometric (simplex), Ruiz+PC (PDHG); tolerances checked on *unscaled* model | Low |
| **LP: dual simplex** | Degeneracy, cycling/stalling, numerical pivot errors, LU fill, update stability | Bounded dual simplex, DSE pricing, bound-flipping ratio test, Harris two-pass + cost perturbation, Forrest–Tomlin LU update, hyper-sparse solves [F: Koberstein 2005; Forrest & Goldfarb 1992; Huangfu & Hall 2018] | Textbook primal simplex with Dantzig pricing + Bland: 10–100× more iterations; stalls on degenerate models | Bounded dual (+ primal cleanup), dual Devex→DSE, Harris + perturbation, own Markowitz LU + FT or PFI update, refactor every ~100 pivots / on error | Medium (instrumented degeneracy telemetry) |
| **LP: IPM** | Ill-conditioned normal equations near optimum; dense columns; free vars; need crossover for a basis | Mehrotra predictor–corrector, homogeneous self-dual (MOSEK/Xpress), supernodal Cholesky, crossover | No ordering → catastrophic fill; no regularization → breakdown | Mehrotra PC on normal eqns or regularized KKT (quasidefinite LDLᵀ), AMD ordering, primal-dual regularization | Low–Medium |
| **LP: first-order (PDHG)** | Only low/medium accuracy (1e-4…1e-6 rel.), tail convergence slow on degenerate/ill-conditioned LPs, no basis | PDLP: restarts, adaptive step, primal weight, presolve, diagonal preconditioning; Halpern variants (cuPDLPx, HPR-LP) | Vanilla PDHG: thousands× slower, never converges on Netlib-hard | Restarted Halpern/averaged PDHG with adaptive restarts & primal weight on GPU | **Now Category A** |
| **Numerical LA / factorization** | Sparse LU with threshold pivoting; stability vs sparsity trade-off; update error growth | Markowitz with threshold partial pivoting (u≈0.1), singleton/triangular detection, FT update | Dense LU → O(m³); no threshold → blow-up | Own sparse LU (Markowitz + threshold), triangularization pass, FT update, residual-triggered refactor | Low |
| **Basis mgmt** | Singular bases after update; restoring feasibility; warm start from parent node | Basis repair by slack substitution; store basis status per node | Restart from slack basis each node → 10–100× more iterations | Basis status arrays, repair, warm start across nodes/solves | **Medium (repeated-solve industrial angle)** |
| **B&B** | Tree explodes; node LP cost; memory | Best-bound + depth-first plunging, dual-simplex warm starts | Pure DFS/BFS, cold node LPs → hopeless beyond toy | Best-estimate/best-bound hybrid with plunging, warm-start dual simplex | Low |
| **Branching** | Choice drives tree size by orders of magnitude | Reliability branching (pseudocosts initialized by strong branching) [F: Achterberg, Koch, Martin 2005] | Most-fractional ≈ random [R] | Pseudocost + limited strong branching (reliability) | Low (batched GPU strong branching = cuOpt 26.02, Blin et al. 2026 → A/B) |
| **Cuts** | Separation cost, numerics (dense, tiny coefficients), cut selection | GMI, MIR, knapsack cover, flow cover, clique, implied bound; efficacy/parallelism selection | Adding every GMI → numerically dangerous, LP bloat | Root GMI + MIR + cover with efficacy/orthogonality filter and coefficient safety checks | Low |
| **Heuristics** | Incumbent early prunes tree; weak LP relaxations need primal side | Rounding, diving, feasibility pump, RINS, local branching, feasibility jump (GPU in cuOpt) | No heuristics → no incumbent on big models | Simple/shift rounding, fractional + guided diving, feasibility pump (own LP), RINS if time | Low |
| **Node selection** | Balance bound improvement vs finding incumbents | Best-estimate with plunging | Pure BFS memory blow-up; pure DFS never proves | Best-bound with periodic plunging | Low |
| **Bound propagation** | Activity-based bound tightening at each node | Standard (SCIP, CPLEX), GPU version exists [R: Sofranac, Gleixner, Pokutta 2022] | None → weak nodes | Activity-based propagation with incremental activity updates | Low |
| **Parallelism** | Determinism, load balance | Deterministic parallel B&B (Gurobi), racing, UG/ParaSCIP | Nondeterministic threads → irreproducible demo | Parallel node LPs with deterministic sync points, or concurrent LP (simplex ‖ IPM ‖ PDHG) | Low–Medium |
| **GPU** | Irregular sparse memory, transfer costs, precision | PDHG on GPU (everyone), cuDSS barrier (cuOpt, COPT), GPU heuristics (cuOpt), batched PDLP | GPU dense simplex tableaus → impressive only on dense toys | Custom CSR SpMV + fused PDHG kernels, batched SpMM | **Only as a measured cost-model/routing contribution** |
| **Termination** | What does "optimal" mean? tolerances on scaled vs unscaled | Primal/dual feasibility 1e-6..1e-9 abs/rel, gap for MIP (1e-4 default in Gurobi/CPLEX) | Declaring optimal on scaled problem → unscaled violation | Two-level: solver tolerances + independent unscaled verification | Medium (certified termination) |
| **Verification** | Solver can "lie" (numerical false optimum/infeasible) | Rarely exposed; exact solvers (SoPlex iterative refinement, QSopt_ex, exact SCIP, VIPR certificates) | Trust log → hidden wrong answers | Independent checker: primal residuals, dual residuals, complementary slackness, safe dual bound, integrality | **High as engineering credibility** |
| **Benchmarking** | Cherry-picking, variance (MIP performance variability [R: Lodi & Tramontani 2013]) | Mittelmann shifted geometric means, fixed hardware, time limits | Selective showcase | Full-suite run with failure table, SGM(10s), multiple seeds for MIP | Medium |
| **Observability** | Understanding why a run failed | Solver logs; Gurobi "numerical focus" diagnostics | Opaque | Structured JSON telemetry (degenerate pivots, LU refactors, condition estimates, kernel timings) | Medium (fits "transparency") |

---

## 4. Mathematical Foundations (what must be implemented, precisely)

**LP (internal form):** `min cᵀx s.t. L ≤ Ax ≤ U, l ≤ x ≤ u`. Introduce row activities `r = Ax` as logical variables with bounds `[L,U]` → `[A −I][x;r] = 0`. Basis `B` of size m; nonbasic variables sit at a bound (or zero if free).

**Dual simplex iteration (bounded):** choose leaving row p with primal infeasibility (DSE: maximize `infeas_p²/‖e_pᵀB⁻¹‖²`); compute `ρ_p = e_pᵀB⁻¹` (BTRAN), pivot row `α_p = ρ_pᵀ A_N` (row-wise PRICE); bound-flipping ratio test over `d_j/α_pj` (Harris tolerance); FTRAN entering column; update primal, dual, DSE weights; update LU (Forrest–Tomlin). Three sparse triangular solves per iteration dominate.

**Optimality certificate (LP):** primal feasible x, dual (y, z) with `z = c − Aᵀy`, sign conditions consistent with active bounds, complementary slackness. **Infeasibility certificate:** Farkas ray y with `Aᵀy` bounds-consistent and positive infeasibility measure. **Unboundedness:** primal ray.

**Safe dual bound (Neumaier–Shcherbina 2004) [F]:** for *any* y (e.g., from PDHG or a slightly wrong simplex), with finite bounds `l ≤ x ≤ u`: `r = c − Aᵀy` computed with outward rounding; `LB = Σ_i (y_i>0 ? y_i L_i : y_i U_i) + Σ_j min(r_j l_j, r_j u_j)` evaluated with directed rounding is a **rigorous** lower bound ([paper](https://link.springer.com/article/10.1007/s10107-003-0433-3)). This is the single most useful piece of mathematics for an "honest GPU" story: low-accuracy GPU duals still give *valid* pruning bounds when variables are bounded.

**IPM (Mehrotra PC):** solve regularized KKT `[−(Q+X⁻¹Z+ρI) Aᵀ; A δI]` by quasidefinite LDLᵀ (static ordering, no pivoting needed thanks to regularization — Vanderbei's quasidefinite result), predictor/corrector with centering σ=(μ_aff/μ)³.

**PDHG for LP (PDLP family):** `x⁺ = proj_X(x − τ(c − Aᵀy))`, `y⁺ = proj_Y*(y + σ(b − A(2x⁺ − x)))`; restarts based on normalized duality gap / KKT error; adaptive steps; primal weight ω balancing τ/σ; Halpern anchoring (cuPDLPx, HPR-LP). Per iteration: one SpMV with A, one with Aᵀ, vector ops — memory-bandwidth-bound.

**QP:** convex `½xᵀQx + cᵀx`. Options: (i) IPM sharing LP KKT code (best accuracy), (ii) ADMM/OSQP-style with own quasidefinite LDLᵀ, (iii) PDHG-type (PDQP, rAPDHG) on GPU. Must check convexity (Q PSD) — LDLᵀ of Q+εI inertia or refusal with certificate.

**MILP:** LP-based B&B; node bound = LP value (or safe bound); prune if `LB ≥ incumbent − gap_abs`. Cuts: GMI from tableau row `x_B + Σ ā_j x_j = b̄`, MIR from aggregated rows, lifted cover from knapsack rows.

---
## 5. Existing Solver Ecosystem — Comparison Matrix

| System | Licence | Classes | LP algorithms | MIP | GPU (as of Sep 2026) | Known weaknesses | Comparable to us? |
|---|---|---|---|---|---|---|---|
| **Gurobi 13** | Commercial | LP, QP, QCP, MILP, MIQP, NLP (local barrier new in 13), global MINLP | Primal/dual simplex, barrier, **PDHG** (CPU+GPU) | Full B&C, deterministic parallel | GPU PDHG, fully supported since 13.0.2 [F: [Gurobi 13 notes](https://docs.gurobi.com/projects/optimizer/en/current/reference/releasenotes/changes.html)] | Cost, black box | Reference only |
| **CPLEX** | Commercial (IBM) | LP, QP, QCP, MILP, MIQP | Simplex, barrier, network | Full B&C, Benders | None documented | Slower development cadence [I] | Reference; likely MRPL incumbent via PIMS |
| **FICO Xpress 9.8** | Commercial | LP, QP, SOCP, MILP, MINLP, global | Simplex, barrier, **PDHG ("hybrid gradient")** | Full B&C | GPU PDHG beta, "up to 50×" on selected large LPs [F: [FICO blog](https://www.fico.com/blogs/gpu-acceleration-hybrid-gradient-algorithm-fico-xpress)] | Cost | Reference; embedded in PIMS/RPMS |
| **MOSEK** | Commercial | LP, conic (SOCP, SDP, exp), QP, MICP | Homogeneous self-dual IPM, simplex | B&B (weaker than MIP leaders) | None known | MIP not its strength | QP/conic reference |
| **COPT (Cardinal, China)** | Commercial | LP, QP, SOCP, SDP, MILP | Simplex, barrier, PDLP | Full B&C | GPU PDLP and **GPU barrier** (only commercial GPU barrier per FICO/GAMS comparison) [F] | — | **Key "sovereign solver" precedent** |
| **MindOpt (Alibaba), OptVerse (Huawei)** | Commercial (China) | LP, MILP, QP | Simplex/IPM | B&C | Varies | — | Sovereign precedent: national solvers reached top of Mittelmann LP tables within ~5 years [I from Mittelmann history] |
| **HiGHS** | MIT | LP, QP (convex, active set), MILP | Dual simplex (parallel PAMI/SIP), IPX IPM, new HiPO IPM, cuPDLP-C (GPU) | B&C (Gomory, MIR, cover…), symmetry | GPU PDLP | MIP well behind commercial | **Best reference solver for us** (open, same classes) |
| **SCIP 10 / SoPlex** | Apache 2.0 | MILP, MINLP, CIP, exact MIP | SoPlex simplex with iterative refinement / rational | Most extensible B&C framework; exact rational MIP; VIPR certificates | None native | Slower than commercial on MILP | Reference for exactness/certification |
| **CBC / Clp / Cgl** | EPL | LP, MILP | Clp primal/dual simplex, barrier | CBC B&C with Cgl cuts | None | Aging, slower | Weak baseline (don't compare only against this — reviewers see it as cherry-picking) |
| **GLPK** | GPL | LP, MILP | Simplex, IPM | Basic B&C | None | Very slow on hard MIP | **Too weak as sole baseline** |
| **NVIDIA cuOpt** | Apache 2.0 | LP, MILP, routing | PDLP, dual simplex, cuDSS barrier, concurrent | CPU B&B + GPU heuristics, root cuts (26.02+), batch-PDLP strong branching (26.02) [F: [release notes](https://archive.docs.nvidia.com/cuopt/user-guide/26.04.00/release-notes.html)] | Core design | Uses PaPILO + PSLP presolve (not from scratch) | **Closest architectural analogue** |
| **PDLP (OR-Tools)** | Apache 2.0 | LP, QP (diag.) | Restarted PDHG | — | CPU (GPU via cuPDLP family) | 1e-4..1e-8 accuracy, tail convergence | Algorithm source |
| **cuPDLP.jl / cuPDLP-C / cuPDLPx** | MIT/Apache | LP | GPU restarted (Halpern) PDHG | — | Yes | Same | Algorithm source |
| **HPR-LP** | Open | LP | GPU Halpern–Peaceman–Rachford | — | Yes | — | ~Parity w/ cuPDLPx on Mittelmann, faster on MIPLIB relaxations [R: [2509.23903](https://arxiv.org/pdf/2509.23903)] |
| **OSQP / cuOSQP** | Apache | Convex QP | ADMM + LDLᵀ / PCG on GPU | — | cuOSQP | Moderate accuracy | QP algorithm source |
| **qpOASES** | LGPL | Small dense QP | Parametric active set | — | — | Small problems only | MPC-style QP reference |
| **IPOPT / MadNLP** | EPL / MIT | NLP | Filter line-search IPM | — | MadNLP has GPU condensed-space IPM | — | NLP extension reference |

## 6. Commercial Solver Analysis (what they do better, honestly)

**[F/R]** Commercial MIP performance is an accumulation of hundreds of components; Achterberg & Wunderling (2013) show, for CPLEX, that cutting planes and presolve are each worth large factors on hard models, with branching, heuristics, node presolve and conflict analysis each contributing further. **Implication [I]:** a student solver will be **10²–10⁴× slower** than Gurobi on hard MIPLIB instances and will time out on most of the 240-instance benchmark set. Saying this out loud, with data, is a *credibility asset*.

What commercial systems do that you cannot match in 10 weeks: tuned DSE/LU/hypersparsity; decades of numerical safeguards; 20+ cut families with calibrated selection; deterministic parallel B&C; symmetry; conflict analysis; restarts; ML-tuned parameter defaults.

What they do **not** give MRPL: source visibility; the ability to add refinery-specific presolve/branching/cuts; control over tolerances semantics and certificates; hardware independence (NVIDIA-only GPU paths); freedom from per-core licensing; auditability of *why* a plan was chosen.

## 7. Open-Source Solver Analysis

- **HiGHS** is the realistic open reference: MIT licence, dual simplex by Huangfu & Hall with parallel variants [F: Huangfu & Hall, *Math. Prog. Comp.* 2018], IPX, active-set QP, B&C MIP, GPU PDLP via cuPDLP-C. **If MRPL just needs "not foreign-licensed", HiGHS already exists.** Your answer to "why not HiGHS" must therefore be about *capability you build and own*, not licence (see §7, §10).
- **SCIP/SoPlex** is the reference for **exactness**: iterative refinement for LP [R: Gleixner, Steffy, Wolter 2016], exact rational MIP [R: Eifler & Gleixner 2023], VIPR certificates [R: Cheung, Gleixner, Steffy 2017]. Borrow the *ideas* (safe bounds, certificate format), not code.
- **CBC/GLPK**: fine for correctness cross-checks; comparing *speed* only against them looks like cherry-picking.

### 7b. "Why not just use HiGHS / SCIP / CBC?" — the strongest honest answer
- **Today, for raw performance, they are better than anything a student team will build.** HiGHS's dual simplex, IPM and MIP embody a decade of specialist work; SCIP is the most complete open B&C framework. Say this first.
- **Licence is not the gap** — HiGHS is MIT, SCIP is Apache 2.0. MRPL *could* adopt them tomorrow.
- **The gap we target is capability, not code:** (1) a team in India that understands and can modify every line of an LP/MIP engine (HiGHS/SCIP are maintained abroad; forking ≠ owning the know-how); (2) **certified-by-default outputs** incl. GPU first-order results (neither ships safe-bound certification on its default path; exact SCIP does, at large cost); (3) a **measured CPU↔GPU routing model** — HiGHS integrates cuPDLP-C but does not claim a validated router; (4) **repeated-solve refinery workflows** as a first-class API (case families, warm-start chains, batched GPU).
- **What we will *not* claim:** better performance, broader features, or more robustness than HiGHS/SCIP. We claim *equal-or-honestly-worse* performance on a published benchmark, *with* certification and routing they don't expose.

### 7c. "Why not just use Gurobi / CPLEX / Xpress?"
- **They are better today** on essentially every performance metric, and Gurobi/Xpress already ship GPU PDHG.
- **What they don't give MRPL:** source access and the ability to add refinery-specific presolve, branching or cuts; independent verifiability (their "optimal" is a black-box claim); control of tolerance semantics; freedom from per-core/per-user licences and export/licensing policy changes; hardware choice (their GPU paths are NVIDIA-CUDA-only).
- **What our prototype makes controllable:** algorithms, tolerances, certificates, hardware routing, extension points. It does *not* make MRPL independent tomorrow — it's a foundation, and the deck should say so.

### 7d. Is "sovereign" a technical property? Only if you define it as one
| Property | Engineering definition | Evidence you can show |
|---|---|---|
| Source ownership | 100% of solver-core code authored by team; no solver-library deps | SBOM, link audit, git history |
| Dependency control | Allow-listed deps only (C++ std, CUDA runtime); every dep replaceable | `deps.allow` + CI |
| Auditability | Every result ships a checkable certificate; independent checker ~1–2k LoC anyone can read | certificate.json + checker |
| Reproducibility | Deterministic runs; pinned builds; one-command benchmarks | bit-identical logs |
| Deployment independence | Runs offline on commodity CPU; GPU optional | CPU fallback demo |
| Hardware portability | Backend interface; CUDA today, HIP/SYCL later (CUDA-only is itself lock-in) | `LinOpBackend` abstraction |
| Algorithmic transparency | Structured telemetry explaining pivots, cuts, branching decisions | JSON logs |
| Extensibility | Plugin API for branching rules, cuts, heuristics, presolve | 30-line custom rule demo |
| Maintainability | Tests, docs, design notes per module | coverage report |

**[F]** International precedent: China's COPT, MindOpt and OptVerse show a national solver can reach the top tier — but it took funded professional teams years, not a hackathon.

## 8. GPU Optimization Landscape (2025–26 state)

| Area | State | Key refs |
|---|---|---|
| GPU LP by first-order methods | **Mature, commercialized**: PDLP → cuPDLP.jl (2023) → cuPDLP-C → cuPDLPx (2025, 2.5–6.8× over cuPDLP) and HPR-LP; in Gurobi, Xpress, COPT, HiGHS, cuOpt | [PDLP](https://arxiv.org/abs/2106.04756), [cuPDLP.jl](https://arxiv.org/pdf/2311.12180), [cuPDLPx](https://arxiv.org/abs/2507.14051), [HPR-LP](https://link.springer.com/article/10.1007/s12532-025-00292-0), [survey](https://arxiv.org/pdf/2506.02174) |
| GPU barrier | Real but narrow: cuOpt barrier with cuDSS (25.10), COPT GPU barrier | cuOpt notes |
| GPU simplex | ~15 years of papers (Spampinato & Elster 2009; Lalami et al.; Ploskas & Samaras) — wins mostly on **dense** LPs; sparse revised simplex has little exploitable parallelism per iteration [R: Hall 2010; Huangfu & Hall 2018] | [1803.04378](https://arxiv.org/pdf/1803.04378), [2211.10979](https://arxiv.org/pdf/2211.10979) |
| Batched LP on GPU | Gurung & Ray (2018/19) small dense batches; **Blin, Gualandi, Maes, Lodi, Stellato (Jan 2026): batched PDHG via SpMM for strong branching and OBBT**; shipped as batch-PDLP strong branching in cuOpt 26.02 | [1802.08557](https://arxiv.org/pdf/1802.08557), [2601.21990](https://arxiv.org/abs/2601.21990) |
| GPU B&B | Problem-specific: B³-PWL (Aug 2026) GPU-batched B&B for SOS2/PWL, 9.25× GM over cuOpt; k-sparse GLM batched B&B (May 2026) | [2608.28988](https://arxiv.org/abs/2608.28988), [2605.22188](https://arxiv.org/abs/2605.22188) |
| GPU MIP heuristics | cuOpt: GPU feasibility pump/jump, fix-and-propagate, bound propagation with probing [R: Çördük et al. Oct 2025]; Kempke & Koch: low-precision PDLP + fix-and-propagate, 243M-nnz unit commitment within 2% where commercial solvers found no feasible solution in 2 days; CHAP (May 2026) hybrid GPU/CPU portfolio beats Gurobi & cuOpt on heuristic feasibility (47 vs 44 vs 43 of 50) | [2510.20499](https://arxiv.org/abs/2510.20499), [2503.10344](https://arxiv.org/abs/2503.10344), [2605.05086](https://arxiv.org/abs/2605.05086) |
| GPU propagation | Sofranac, Gleixner, Pokutta (2022) GPU domain propagation | *Parallel Computing* 2022 |
| GPU presolve | Cederberg & Boyd: PSLP (lightweight presolve for GPU FOMs, Apr 2026) and cuPSLP GPU presolve, 11× (Mittelmann) / 42× (GAMS large) over CPU PSLP (Sep 2026) | [2604.23951](https://arxiv.org/abs/2604.23951), [2609.16182](https://arxiv.org/abs/2609.16182) |
| GPU QP | cuOSQP (2020), PDQP / rAPDHG, PDHCG | survey above |
| GPU NLP | MadNLP condensed-space IPM on GPU (Shin, Pacaud, Anitescu) for ACOPF | — |
| Autotuning FOMs | "Parameter tuning with generalization guarantees for GPU LP" (2026) | [2606.08638](https://arxiv.org/pdf/2606.08638) |

## 9. Prior-Art / Literature Audit

| Idea | Category | Evidence |
|---|---|---|
| GPU-accelerated LP (PDHG) | **A** | Commercial in 5 solvers |
| GPU simplex | **A** (incl. negative results for sparse) | 15 yrs of papers |
| GPU interior point | **A/B** | cuOpt+cuDSS, COPT, MadNLP |
| GPU MILP heuristics | **A** | cuOpt, CHAP, Kempke & Koch |
| GPU branch-and-bound (general MILP) | **B/C** | B³-PWL (PWL only), GLM B&B; general-MILP node LPs still on CPU simplex in cuOpt |
| Batched LP for strong branching / OBBT | **A** (since Jan 2026) | Blin et al.; cuOpt 26.02 |
| Batched PDHG *node bounding with safe bounds* | **B/C** | Neumaier–Shcherbina (2004) + batched PDHG (2026) — combination is obvious; a competitor repo already has a PR for exactly this |
| Parallel / heterogeneous MIP | **A** | UG/ParaSCIP, Gurobi, cuOpt |
| Adaptive CPU/GPU solver selection | **B** | Concurrent LP (Gurobi, cuOpt) = "run all, take first"; ML algorithm selection studied; PIPEPYE claims a "router" |
| GPU presolve | **A/B** (Sep 2026) | cuPSLP |
| GPU cut generation | **C** | Little found — but cut separation is irregular and cheap relative to LP; low payoff [I] |
| Accelerator-aware node scheduling | **C** | B³-PWL batch scheduling; general MILP weakly explored |
| Numerical stabilization for GPU FOMs | **B** | Mixed/FP32 in cuOpt 26.04; restarts/precond. literature |
| Refinery scheduling/planning MILP | **A** | Lee, Pinto, Grossmann & Park 1996; Pinto, Joly & Moro 2000; Jia & Ierapetritou 2003–04; review in *Front. Eng. Manag.* 2020 |
| Crude blending / pooling | **A** | Haverly 1978; distributive recursion; Penalty-DR (Dai group, Nov 2024) [R: [2411.09554](https://arxiv.org/abs/2411.09554)] |
| Indigenous/sovereign solvers | **A** internationally (COPT, MindOpt, OptVerse); in India, only student/hackathon prototypes found | GitHub SIH26119 repos |
| Optimization-as-a-service / orchestration | **A** | NEOS, Gurobi cloud, cuOpt server |
| Benchmarking frameworks | **A** | Mittelmann, MIPLIB scripts, Performance profiles (Dolan & Moré 2002) |
| **Repeated-solve acceleration for refinery case studies (warm-start simplex vs batched GPU FOM, with certified bounds)** | **C** | Parts exist (warm starts: universal; batched FOM: Blin et al.; refinery SLP: PIMS); I found **no published study measuring the crossover for refinery-planning LP families** — but absence of evidence after a short search is not proof |
| **Measured CPU/GPU crossover cost model for LP algorithm routing, validated on held-out instances** | **C** | Many papers report speedups; few give a *predictive, validated* routing model; tuning papers exist (2606.08638) |

### D/E candidates
**Nothing I found survives as Category E.** The best candidates are **C**, i.e., *known pieces, under-explored combination on a specific workload*. For each:

**C1. Repeated-solve refinery LP families (crude valuation / SLP recursion / scenario cases).**
- Closest prior work: Blin et al. 2026 (batched PDHG for *related* LPs in MIP); MPAX batching; PIMS/PDR for SLP.
- Overlap: batching related LPs on GPU via SpMM.
- Difference: workload = many *coefficient- and bound-perturbed* planning LPs where (a) warm-started dual simplex is the incumbent method and very strong; (b) planners need **duals/marginal values** at high accuracy, which FOMs deliver poorly.
- Why still possibly interesting: the question "when does batched GPU FOM + CPU simplex polish beat sequential warm-started simplex, at *equal certified accuracy*?" is operationally real and not answered for this class.
- Proof experiment: family of 100–1000 crude-valuation LPs from a literature refinery model at 3 scales; compare (i) cold simplex, (ii) warm dual simplex chain, (iii) batched GPU PDHG → crossover/polish; report wall-clock to certified 1e-6 gap and dual accuracy.

**C2. Validated routing model.** Predict `t_GPU-PDHG / t_CPU-simplex` from cheap features (nnz, m, n, row/col nnz variance, estimated condition proxy, degeneracy proxy after short simplex probe). Proof: train on half of Netlib+Mittelmann+MIPLIB relaxations+generated refinery LPs, test on the other half; report regret vs oracle choice and vs "always concurrent".

## 10. Novelty Gap Analysis — blunt summary
- "Solver from scratch by students" = **engineering contribution**, not research.
- "GPU PDHG" = **engineering reproduction** of 2021–2025 research.
- "Sovereign/Indian" = **strategic framing**, not technical novelty.
- "Refinery case study" = **application contribution** (only valuable if the model is literature-grounded and results are validated).
- Real gap you *can* credibly occupy: **(a) certified-by-construction outputs across all algorithms including GPU FOMs; (b) a measured, predictive CPU↔GPU routing model; (c) repeated-solve refinery workload characterization.** (a) is engineering-grade, (b)(c) are systems/empirical research at "workshop paper" level if done rigorously.

## 11. Potential Differentiation Opportunities

| Candidate | Exists? | Closeness | Exact difference we'd implement | Matters? | Demonstrable? | Defensible? |
|---|---|---|---|---|---|---|
| **Certified results layer** (safe bounds for every LP/MIP bound, rational re-check of final basis, Farkas/ray checks, VIPR-style MIP log) | Yes in exact SCIP/SoPlex/QSopt_ex; partly in competitor repos | High | Apply uniformly to *every* engine incl. GPU PDHG and batched node bounds; report "certified gap" not "claimed gap" | **Yes** — directly addresses "numerically robust" | Yes: adversarial suite where naïve solver claims wrong optimum, ours flags/repairs | Strong |
| **Measured CPU/GPU routing model** | Concurrent racing exists; ML algorithm selection exists | Medium | Cheap-feature cost model + validated regret on held-out set | Yes (PS literally says "where measurable") | Yes | Strong if validation honest |
| **Refinery repeated-solve engine** (warm-start across cases/recursion passes, batched GPU for many cases) | Pieces exist | Medium | Case-family API; warm-start graph; batch GPU; certified duals for marginal crude values | **Yes for MRPL** | Yes | Medium–strong |
| Degeneracy/ill-conditioning telemetry ("why was this hard?") | Solver logs only | Low–medium | Per-run numerical health report: degenerate pivot %, LU refactor reasons, κ estimates, bound-shift magnitudes | Moderately (transparency) | Yes | Strong, cheap |
| Batched GPU strong branching | cuOpt 26.02, Blin 2026 | Very high | — | Only as reproduction | Yes | Weak as novelty |
| GPU heuristics (FJ/FP) | cuOpt, CHAP | Very high | — | — | — | Weak |
| Domain-aware refinery cuts/presolve (e.g., tank-state disjunctions, pooling-aware bounds) | Literature exists per model | Medium | Detect structures automatically in generic MPS | Potentially, but hard in 10 weeks | Partially | Medium |
| LLM/AI front end, dashboard | — | — | — | **No** | — | Counter-productive |

## 12. CPU vs GPU Workload Analysis

Order-of-magnitude hardware facts **[F, vendor specs]**: A100 HBM ~1.5–2 TB/s, L4 ~300 GB/s, RTX 30/40 consumer ~200–1000 GB/s; PCIe 4.0 ×16 ~25 GB/s practical; kernel launch ~3–10 µs; desktop CPU DRAM ~50–100 GB/s. FP64 throughput on consumer/L4 GPUs is **1/32–1/64 of FP32** — critical: PDHG in FP64 on an RTX card is bandwidth-bound anyway, but dense factorizations in FP64 are crippled on consumer GPUs.

| Operation | Arith. intensity | Parallelism | Irregularity / divergence | Sync | CPU fit | GPU fit | Transfer |
|---|---|---|---|---|---|---|---|
| SpMV `Ax`, `Aᵀy` (PDHG, IPM-CG, propagation) | ~0.15–0.25 flop/B | nnz-wide | Row-length variance → load imbalance | Per iter reductions | Good | **Best** when nnz ≳ 10⁶ | One-time upload |
| Batched SpMM (k related LPs) | ↑ with k (A reused) | nnz × k | Same | Per iter | OK | **Best** | One-time |
| Simplex FTRAN/BTRAN (sparse triangular, hypersparse) | Very low | Tiny (often <100 nnz touched) | Extreme | Each iter | **Best** | Poor | — |
| Simplex pricing (DSE over m rows) | Low | m-wide | Low | Each iter | Good | OK only for very large m, but latency kills it | per-iter |
| Row-wise PRICE `ρᵀA_N` | Low | nnz(ρ)-dependent | High | Each iter | Good (hypersparse) | Poor | — |
| Sparse LU factorization | Moderate | Limited (elimination tree levels) | High | Many | **Best** | Poor at simplex sizes | — |
| Supernodal Cholesky/LDLᵀ (IPM) | High in supernodes (dense BLAS3) | Tree-level + dense | Moderate | Tree | Good | Good for large, fill-heavy KKT (cuDSS) | Per-factorization if CPU-side |
| B&B tree management | ~0 | Node-level | Very high | Global incumbent | **Best** | Poor | — |
| Node LP (warm dual simplex, ~10–100 pivots) | Low | None | High | — | **Best** | Poor | — |
| Many node LPs in batch (FOM) | Moderate | Batch | Moderate | Batch | OK | Good *if* approximate bounds suffice (with safe-bound correction) | Bounds per batch |
| Bound propagation | Low | Constraint-wide | Moderate | Fixpoint iterations | Good | Good for large models [R: Sofranac 2022] | Per node — bad |
| Cut separation (GMI/MIR) | Low | Row-wise | High | — | **Best** | Poor | — |
| Heuristics (FJ, local search) | Low | Massive (many moves/candidates) | Moderate | Weak | OK | Good [R: cuOpt, CHAP] | Low |
| Presolve | Low | Rule-dependent | High | Many passes | Good | Good only at very large scale [R: cuPSLP] | — |
| QP ADMM (with factorization) | Mixed | — | — | — | Good | OK with PCG | — |

**Workload model [I, to be calibrated by you]:** For one PDHG iteration, `t_GPU ≈ max(2·bytes(A)/BW_gpu, L_launch·k_kernels)` and `t_CPU ≈ 2·bytes(A)/BW_cpu`. With bytes(A) ≈ 12·nnz (8 B value + 4 B index), GPU wins per-iteration only when `12·nnz·2/BW_gpu > ~k·5µs` → nnz ≳ 10⁵–10⁶ for ~10 fused kernels. **But the real comparison is PDHG-GPU vs simplex-CPU**, not PDHG-GPU vs PDHG-CPU: for typical refinery planning LPs (10⁴–10⁵ rows, 10⁵–10⁶ nnz), a warm-started dual simplex converges to 1e-9 in fewer, cheaper steps than PDHG reaches 1e-6. Competitor data points agree: SANKHYA reports GPU PDHG *slower below ~2,000 rows*, and only 1.37× on an RTX 5050 at 10⁴ rows.

## 12b. Challenging the GPU Assumption

**Where GPU makes sense:** (1) very large LPs (≥10⁶ nnz) needing moderate accuracy; (2) many related LPs solved together (batching amortizes launch overhead and reuses A in cache); (3) primal heuristics with massive independent moves; (4) large IPM factorizations with heavy fill (via cuDSS — but that's a library solver, see §15).

**Where GPU makes things worse:** simplex iterations; small/medium LPs; node-by-node B&B; anything needing per-iteration host round-trips (PCIe latency ~10 µs dominates a 20 µs kernel); FP64-heavy dense work on consumer GPUs; degenerate/ill-conditioned LPs where FOM tails blow up; anything requiring an optimal *basis* (crossover needed — CPU).

**Design rules [Rec]:** keep the whole PDHG loop on device (no per-iteration host sync; check convergence every ~64 iterations); fuse vector ops; CSR with row-length-aware kernels (vector-CSR/merge-path) for imbalanced rows; FP64 by default, FP32 SpMV + FP64 accumulation as an experiment; pinned memory + async copies for batch updates; **always report GPU results with the transfer and setup time included**.

## 13. Proposed Architecture Options

| Arch | Idea | Strengths | Weaknesses | Complexity | Novelty | Hackathon feasibility | MRPL relevance | Risk |
|---|---|---|---|---|---|---|---|---|
| **A** CPU-first, GPU for linear algebra | Simplex/IPM/B&B on CPU; GPU for SpMV/factorization | Mathematically sound | GPU adds little for simplex; cuDSS dependency questionable | Medium | Low | High | High | GPU looks decorative |
| **B** GPU-first | PDHG everywhere, GPU B&B | Flashy speedups on big LPs | Low accuracy, no bases, poor on refinery-size degenerate LPs; MIP proofs hard | High | Low (A-category) | Medium | **Low** | Judges catch accuracy/tolerance games |
| **C** Heterogeneous by role | CPU: simplex, IPM, B&B, cuts; GPU: PDHG (LP/QP), batched bounds, heuristics; certifier on CPU | Each op where it fits; honest | Engineering breadth | High | Low–Med | Medium | High | Breadth → shallow everywhere |
| **D** Adaptive runtime / router | C + cost model picks engine or races them | Directly answers "measurable benefit" | Needs calibration data | High | **Med (C-cat)** | Medium | High | Model may be unconvincing if small data |
| **E** Industrial repeated-solve specialization | Case-family API, warm-start chains, batched scenarios | Speaks MRPL's language | Narrower | Medium | **Med (C-cat)** | High | **Very high** | Requires a credible refinery model |
| **F** Certified solver ("trust-first") | Every answer carries checkable certificate; rational re-verification | Unique credibility, cheap to add | Not a speed story | Low–Med | Low (A/B) but rare in practice | **Very high** | High | None major |

## 14. Recommended Architecture and Why
**[Rec] C + F core, with D and E as the two measured contributions.**
1. **F** is cheap, defends against every "how do you know it's correct" attack, and turns robustness from adjective into evidence.
2. **C** is the evidence-aligned hardware split (§12).
3. **D** turns the PS's "measurable benefit" clause into an experiment you own.
4. **E** is the one axis where the competitor repos I saw are *not* focused, and it's MRPL's actual workload [I].

Do **not** try to out-feature SANKHYA on cut families or NLP. Out-*prove* and out-*measure* everyone.

---
## 15. From-Scratch Implementation Boundary

PS wording: *"shall not be built upon any existing open source solver library but shall be built from scratch from mathematical foundation."* The operative object is **"solver library"**. It says nothing about language runtimes or general numerical kernels. This is interpretation, not legal advice — **ask the MRPL mentor/SPOC in writing** and keep the reply.

| Tier | Items | Reasoning |
|---|---|---|
| **Definitely safe** | C++ standard library, CUDA runtime/driver API, own CUDA kernels, OpenMP/std::thread, pybind11 (binding only), CMake, GoogleTest, JSON libs, plotting in Python, CUB/Thrust primitives (sort, scan, reduce) | Language/runtime/plumbing; no optimization algorithm inside |
| **Probably safe, disclose** | BLAS level-1/2/3 (OpenBLAS/MKL, cuBLAS) for dense sub-blocks; cuSPARSE SpMV (but writing your own is easy and better for the story); Eigen for *tests* only | General LA primitives; still, the more your core math is yours, the stronger the "from scratch" claim |
| **Questionable** | LAPACK dense factorizations inside the core; cuSOLVER dense; **cuDSS / CHOLMOD / MUMPS / SuiteSparse (sparse direct solvers)**; AMD/METIS orderings | Sparse factorization is *the* numerical heart of simplex/IPM; judges may view outsourcing it as outsourcing the solver. Write your own LU/LDLᵀ and minimum-degree ordering |
| **Definitely risky / violating** | HiGHS, SCIP/SoPlex, CBC/Clp/Cgl, GLPK, OSQP, qpOASES, IPOPT, PDLP/cuPDLP code, PaPILO/PSLP presolve, cuOpt, OR-Tools, `scipy.optimize.linprog/milp` (=HiGHS), PuLP/Pyomo *as the solver* | These are solver libraries |

Use reference solvers **only** in a separate `bench/` harness, never linked into the core. **Proof mechanism:** CI job that runs `nm`/`ldd`/`otool -L` on the solver binary and a dependency allow-list check; SBOM (SPDX) in repo; `git log` authorship.

## 16. Numerical Robustness Architecture

### Tolerances (defaults, all on *unscaled* model at the end)
primal feas 1e-6 abs + 1e-9 rel·‖row‖; dual feas 1e-7; pivot tolerance 1e-7 (LU threshold u=0.1, raise to 0.5/0.9 on instability); Harris tolerance 1e-9 → 1e-7 adaptively; integrality 1e-5 (report 1e-6 separately); MIP rel. gap 1e-4 (report also 0 gap proofs).

### Mechanisms
1. **Scaling:** geometric-mean iterations (≤4) + equilibration for simplex; Ruiz (10 passes) + Pock–Chambolle (α=1) for PDHG; record scale factors, unscale before verification.
2. **Anti-degeneracy:** random cost perturbation (dual) and bound perturbation (primal) of order 1e-7·(1+|c_j|), removed at end with cleanup iterations; Harris two-pass ratio test; bound-flipping ratio test; stall detector (objective unchanged for k pivots) → increase perturbation.
3. **Anti-cycling:** perturbation + randomized tie-breaking; Bland fallback after stall threshold (last resort).
4. **LU stability:** threshold Markowitz; refactor every ~100 updates or when `‖B x − a_q‖/‖a_q‖ > 1e-9`; FT update growth monitoring; singular basis → replace offending columns with slacks (basis repair), log event.
5. **Residual checks:** after each refactor recompute x_B and y from scratch; compare with updated values; discrepancy > tol → refactor + possibly tighten pivot tolerance.
6. **Condition estimation:** Hager/Higham 1-norm estimator on B (cheap, uses FTRAN/BTRAN) at each refactor → reported in telemetry.
7. **Iterative refinement:** on final basis, solve `Bx=b` with residual in long double / double-double; optionally **rational re-verification** of the final basis via exact rational LU (GMP or own big-rational) for small/medium problems — the SoPlex/QSopt_ex idea, from scratch.
8. **Infeasibility/unboundedness:** Farkas ray from dual simplex (unbounded dual ray) / primal ray; verify ray independently before reporting status.
9. **IPM:** primal-dual regularization (1e-8…1e-10), Mehrotra step-length safeguard, detect non-convergence via residual stagnation → hand off to simplex.

### Independent Numerical Safety / Validation Layer ("Certifier")
Separate module, separate code path, re-reads the **original** MPS (own second parser or the same parser + checksum), receives only `(status, x, y, z, basis?, ray?)`.

| Detects | Check |
|---|---|
| Infeasible "optimal" solution | max row/bound violation, absolute and relative (unscaled) |
| Integrality violation | max |x_j − round(x_j)| over integers; also re-evaluate objective with rounded integers and re-check feasibility |
| False optimality (LP) | dual feasibility of (y,z) + complementary slackness + |primal obj − dual obj| gap; **safe dual bound** (N–S) with directed rounding |
| False MIP optimality | recompute objective; gap = incumbent − certified global bound (min safe bound over open leaves at termination); optionally VIPR-style certificate for small instances |
| False infeasibility | verify Farkas ray: `yᵀA` bounds-consistent and certificate value > tol |
| Numerically suspicious | large basis κ estimate (>1e10), large perturbation shifts, many refactors, large |x| (>1e9), big residual ratios, status flip-flops between runs with different seeds |
| Poor conditioning | row/col coefficient ranges, κ(B) estimate, scaled-vs-unscaled violation divergence |

Output a `certificate.json` with PASS/WARN/FAIL per check. **A run whose certificate FAILs is reported as `NUMERICAL_FAILURE`, never `OPTIMAL`.**

## 17. LP Design

- **Default engine:** bounded **dual simplex** (Phase 1 via artificial bounding "big box" or composite method) + primal simplex for cleanup after perturbation removal and for warm starts when primal feasible.
- **Pricing:** dual Devex first, DSE when stable (DSE is worth it on most Netlib-hard) [F: Forrest & Goldfarb 1992].
- **Ratio test:** bound-flipping (long-step) dual ratio test + Harris [F: Koberstein 2005 thesis — the single most useful reference for implementers].
- **Factorization:** own sparse LU: (1) singleton/triangular extraction ("bump" isolation), (2) Markowitz with threshold on the bump, (3) Forrest–Tomlin update (PFI acceptable in v1).
- **Hypersparsity:** sparse FTRAN/BTRAN using DFS on L/U graphs (Gilbert–Peierls) when rhs is very sparse [F: Hall & McKinnon 2005]. Big speed factor on large sparse LPs.
- **IPM (v2):** Mehrotra PC, regularized quasidefinite LDLᵀ, own approximate minimum degree ordering, crossover to basis (primal + dual push phases).
- **PDHG (GPU):** restarted Halpern PDHG with adaptive restart (KKT-error based), primal weight updates (PID à la cuPDLPx is optional), feasibility polishing optional; termination at 1e-4 / 1e-6 relative; **crossover handoff** to CPU simplex using PDHG x,y to pick a starting basis.
- **Concurrent / routed mode (D):** router picks from features; "race" mode runs simplex on CPU while PDHG runs on GPU; first certified finish wins.

## 18. MILP Design

**Necessary for credibility (v1):** LP relaxation via warm dual simplex; best-bound node selection with depth-first plunging; pseudocost branching initialized by limited strong branching (reliability, η_rel≈4–8); incumbent handling with objective cutoff; pruning by (safe) bound; root presolve + activity-based bound propagation at each node; simple rounding, shift-rounding, fractional diving, feasibility pump; root cuts: GMI + MIR + knapsack cover with efficacy/parallelism filtering and max 5–10 rounds; node limit/time limit returning best incumbent + certified gap.

**Biggest credibility-per-effort items (ranked) [I]:** (1) warm-started node LPs (else nothing works); (2) reliability branching; (3) bound propagation; (4) root GMI+MIR cuts; (5) feasibility pump/diving; (6) presolve with probing; (7) RINS. Symmetry, conflict analysis, restarts, local cuts → v3.

**Parallel:** v1 single-threaded deterministic; v2 parallel node LPs with deterministic epochs; GPU used for: batched safe-bound evaluation of candidate children (strong branching) and GPU heuristics *only if time* — and label them as reproductions of cuOpt/Blin et al.

**Node data structure:** `{parent_id, depth, bound_changes: small vector<(var, lb, ub)>, basis_delta or full basis status (2 bits/var), lp_bound, safe_bound, estimate}`; store diffs from parent to keep memory O(depth) per node.

## 19. QP Design
Convex QP only (declare and certify: attempt LDLᵀ of Q+εI, report inertia; refuse nonconvex with evidence). v1: **IPM** reusing the LP KKT machinery (`Q` block added to (1,1)) → accurate, handles Maros–Mészáros. GPU option: PDHG-for-QP (PDQP-style: linearize Qx term with momentum) — benchmark honestly, mark accuracy. MIQP → v3 (B&B over QP relaxations; same tree code).

## 20. Presolve Design
Rules (v1): remove empty rows/cols; fixed columns; singleton rows → bounds; singleton columns (free/implied free) with dual postsolve; forcing & redundant rows via activity bounds; doubleton equations (substitution, watch fill & coefficient size); dominated columns / dual fixing; duplicate rows/cols (hashing); bound tightening. MIP: coefficient tightening, probing on binaries (limited), clique detection (v2).
**Postsolve stack** records every reduction with the data needed to recover primal **and dual** values; unit test: presolve→solve→postsolve must reproduce the no-presolve objective and pass the certifier on the original model.
Lightweight GPU-oriented presolve is known to be sufficient for FOMs [R: Cederberg & Boyd 2026] — a good experiment: measure presolve effect separately for simplex and PDHG.

## 21. GPU Runtime Design
**Stack [Rec]:** CUDA C++ with own kernels; CUB for reductions/scans; no cuSPARSE in core (optional A/B comparison in bench). **Your dev machine is a Mac (darwin) — no CUDA.** Develop CPU code on Mac; GPU work on an NVIDIA box (college lab, Colab/Kaggle T4/L4, cloud A10/L4/A100). Put the GPU backend behind an interface (`LinOpBackend`) so the CPU reference implementation of every kernel exists and results can be diffed bit-for-bit (up to reduction order). Note the sovereignty angle: CUDA-only is itself vendor lock-in; mention HIP/SYCL portability as future work.

| Kernel | Parallel over | Memory pattern | Sync | Precision | Transfer | Placement |
|---|---|---|---|---|---|---|
| CSR SpMV `Ax` (row-per-warp for long rows, row-per-thread for short; binned by row length) | rows | gather on x | none | FP64 (FP32-store experiment) | none in loop | GPU-native |
| CSR SpMV of `Aᵀ` (store Aᵀ explicitly as CSR — avoids atomics) | cols | gather on y | none | FP64 | none | GPU-native |
| Fused PDHG primal step: `x⁺ = clamp(x − τ(c − Aᵀy), l, u)`, and `2x⁺ − x` | vars | streaming | none | FP64 | none | GPU |
| Fused dual step with projection onto row-bound dual cone | rows | streaming | none | FP64 | none | GPU |
| Norm/KKT-error reductions (every k iters) | block reduce | streaming | grid-level | FP64 (compensated) | 1 scalar every k iters | GPU |
| Batched SpMM for k LPs (A shared, X column-major n×k) | nnz×k | A reused, coalesced on X | none | FP64 | per-batch bounds only | GPU |
| Safe-bound evaluation (N–S) | vars | streaming | reduce | **directed rounding** (`__dadd_rd` / `__dmul_rd` intrinsics exist in CUDA) | 1 scalar | GPU or CPU |
| Ruiz scaling iterations | rows/cols | two passes | per pass | FP64 | once | GPU if large |
| Bound propagation (v3) | constraints | activity recompute | fixpoint | FP64 + safety margins | per node (bad) | GPU only for large root probing |

Adaptive selection: router decides `{CPU simplex, CPU IPM, GPU PDHG(+crossover), race}` from features; log decision + predicted vs actual time (this *is* the D experiment).

## 21b. Data Structures

| Object | CPU representation | GPU representation | Why |
|---|---|---|---|
| Constraint matrix | CSC master (`colptr, rowidx, val`), plus CSR copy | CSR for A and for Aᵀ (32-bit idx, 64-bit vals), row-length bins | Simplex needs columns (FTRAN entering col) and rows (PRICE); GPU wants row-parallel gathers without atomics |
| COO | Only in parser/presolve builders | — | Easy assembly, converts to CSC |
| Bounds/vars | SoA arrays: `lb, ub, cost, type, status(2-bit)` | same SoA | Coalescing, cache |
| Basis | `basic_index[m]`, `var_status[n+m]`, position map | — | O(1) lookups |
| LU factors | L as column etas; U row- and column-wise (for FT); permutations; update etas | — | FT needs row access to U |
| Dual/reduced costs | dense `y[m]`, `d[n+m]`; DSE weights `w[m]` | `y, x` dense | — |
| Hypersparse vectors | `HVector {index list, dense values, count}` | — | Hall & McKinnon hypersparsity |
| Nodes | bound-change diffs + compressed basis; priority queue keyed by (bound, estimate) | batch buffers of bound vectors for batched evaluation | Memory O(depth) |
| Cuts | pool: CSR rows + efficacy, age, origin; global vs local flag | — | Aging/removal |
| Incumbent | `x*`, obj, source heuristic, time found; certified flag | — | Reporting |
| KKT (IPM) | symmetric lower CSC + supernode partition | — | LDLᵀ |
| Working set (active-set QP, v3) | index set + incremental factorization | — | — |

## 22. Industrial Use Cases (literature-grounded; do not fabricate MRPL data)

State explicitly: *"synthetic data generated from published model structures and public assays/price ranges; not MRPL data."*

**A. Crude blending (LP, optional MILP)** — crudes i, blend properties k (sulfur, API, TAN, RVP…). `min Σ c_i x_i` s.t. `Σ x_i = D`, linear-by-volume/weight property limits `Σ (q_ik − Q_k^max) x_i ≤ 0`, availability, contract min-lifts via binaries `x_i ≤ U_i y_i, x_i ≥ m_i y_i`, cargo size integrality (`x_i = S·n_i`, n integer). Size: tens–hundreds of vars; **small** — use it for correctness & explainability (duals = marginal crude values), not performance.
**B. Multi-period refinery planning (LP → MILP)** — units (CDU, VDU, FCC, reformer, HDS, blending pools), modes (binary), yields per crude per mode (fixed-yield LP à la PIMS base), inventories, product specs, demand, T periods. Parametric generator: 20 crudes × 12 units × 15 products × T∈{12, 52, 365} → 10⁴–10⁶ vars, nnz/col ~5–15, strongly block-staircase. *This is the scaling + degeneracy workhorse.* Crude valuation family: 20–200 variants perturbing one crude's price/availability → **repeated-solve experiment (C1)**.
**C. Crude-oil unloading & tank scheduling (MILP)** — Lee, Pinto, Grossmann & Park (1996) benchmark structure: vessels, storage tanks, charging tanks, CDUs, time slots; binaries for transfers; big-M → **weak LP relaxation** (classic). Use for B&B/cuts demonstration; compare tightened vs big-M formulation.
**D. Unit commitment / economic dispatch (MILP / LP / QP)** — generators g, hours t; binaries u,v,w (on/start/stop), min up/down, ramping, reserve; quadratic cost → QP variant (convex). Use tight formulation (Morales-España et al. 2013) vs 3-bin standard to show relaxation strength. Sizes 10³–10⁶ with 100–1000 units × 24–168 h.
**E. Supply chain (MILP)** — plants → depots → customers, fixed-charge facility location + flows; huge LP relaxations, weak without cuts → cover/flow-cover effect.
**F. Adversarial suite** — (i) Klee–Minty (exponential Dantzig pricing), (ii) highly degenerate assignment/transportation LPs (massive dual degeneracy), (iii) Netlib degenerate/hard subset (`degen3, greenbea, greenbeb, pilot87, pilot-ja, d2q06c, perold, scfxm*`), (iv) scaled-badly models: multiply random rows/cols by 10^U(−6,6), (v) near-singular: duplicated rows with 1e-12 perturbations, (vi) weak relaxations: big-M scheduling and `markshare` from MIPLIB (tiny but nearly unsolvable), (vii) infeasible Netlib set (`infeas/`) for Farkas certificates, (viii) unbounded constructed models.

For each: publish generator + seed + MPS files, and evidence: certified objective, runtime vs reference, dual values sanity (e.g., marginal crude value = price change breakeven, verified by re-solve).

## 23. Benchmarking Strategy

**Sets:**
- **Netlib LP** (all ~94 feasible + infeasible set) — correctness vs published optima (Koch 2004 exact values).
- **Mittelmann LP feasible benchmark** (large) — scalability; also Mittelmann's first-order/PDLP benchmark subsets.
- **MIPLIB 2017** — "benchmark" set is 240 instances [F: Gleixner et al., *Math. Prog. Comp.* 2021]; pre-register an "easy & small" subset rule (e.g., status=easy, ≤10⁵ nnz, solved by HiGHS < 60 s) *before* running yours; plus LP relaxations of the full benchmark set for PDHG.
- **Maros–Mészáros** convex QP (138) + **QPLIB** convex continuous subset.
- Generated industrial families (§22) + adversarial suite.

**Reference solvers:** HiGHS (must), SCIP/SoPlex (exactness reference), optionally Gurobi academic licence (only if licence terms allow publishing numbers — check), cuOpt / cuPDLPx for GPU LP comparison (fairness: same GPU).

**Metrics & why:** wall time (incl. read/presolve/transfer) — honest end-to-end; objective & rel. diff vs reference — correctness; certified gap — trust; max primal violation, dual residual, integrality violation — robustness; iterations/nodes/cuts — algorithmic behavior; LU refactors & κ estimates — numerics; memory peak; CPU utilization; GPU kernel time vs transfer time vs idle (Nsight Systems) — proves where time goes.

**Protocol:** fixed hardware described precisely; time limits (LP 1 h large / 300 s medium; MIP 600 s; QP 300 s); 3 seeds for MIP (performance variability); **failure criteria**: wrong status, certifier FAIL, >1e-6 rel. objective error (LP), timeouts counted as time limit; **aggregate**: shifted geometric mean (shift 10 s) over *all* instances with timeouts at limit, plus solved-count and performance profiles (Dolan–Moré); **speedup**: SGM ratio, never best-case only; **reproducibility**: one script (`bench/run_all.sh`) produces all tables from pinned commit + container; raw logs committed.

**Inclusion rules published up front**, and a **full failure table** in the deck appendix.

## 24. Validation Matrix

| # | Category | Instances | Expected evidence | Pass criterion |
|---|---|---|---|---|
| 1 | Parser correctness | Netlib RANGES/BOUNDS-heavy (`boeing1/2`, `forplan`, `tuff`) | Objective match | 1e-8 rel |
| 2 | LP correctness | Netlib full | Match Koch exact optima | ≥85% at 1e-6 rel, 0 false OPTIMAL |
| 3 | Infeasibility | Netlib `infeas/` | Verified Farkas rays | 100% certificates verify or status UNKNOWN (never wrong) |
| 4 | Unboundedness | constructed | Verified rays | 100% |
| 5 | Degeneracy | `degen3`, assignment LPs, Klee–Minty | Iterations, degenerate-pivot %, no cycling | Terminates; telemetry shows perturbation effect (ablation on/off) |
| 6 | Ill-conditioning | scaled-badly variants, `pilot*` | Violations before/after scaling; κ estimates | Certifier PASS with scaling; documented FAIL/flag without |
| 7 | LU/update stability | long runs (`d2q06c`, `dfl001`) | Refactor counts, residual growth | No silent drift (residual checks catch) |
| 8 | LP scalability | Mittelmann + refinery T∈{12,52,365} | Runtime vs nnz log–log | Report slope; compare HiGHS |
| 9 | Presolve | all LP/MIP | Reduction %, time, postsolve correctness | 100% postsolve certify |
| 10 | GPU benefit | large LPs ≥10⁶ nnz | CPU-PDHG vs GPU-PDHG vs CPU-simplex, incl. transfers | GPU faster where predicted |
| 11 | GPU non-benefit | small/medium LPs | Router picks CPU; GPU slower | Router regret < X% vs oracle |
| 12 | Accuracy trade-off | same LPs at 1e-4/1e-6/1e-8 | PDHG time explodes at tight tol; crossover cost | Honest curve |
| 13 | Batched solve | crude-valuation families k∈{1,8,32,128} | Throughput vs k; vs warm simplex chain | Crossover point found |
| 14 | MILP correctness | MIPLIB easy subset | Optimal/gap vs published | 0 wrong optimum claims |
| 15 | Weak relaxation | crude unloading big-M vs tight; `markshare` | Root gap, nodes, effect of cuts | Cuts reduce root gap (ablation) |
| 16 | Branching ablation | same | most-fractional vs pseudocost vs reliability | Node counts |
| 17 | Heuristics ablation | same | Time to first incumbent | — |
| 18 | Safe bounds | node bounds from simplex vs PDHG | Fraction of prunes kept after N–S correction | 0 invalid prunes (verified vs exact rational on small) |
| 19 | QP | Maros–Mészáros | Objective vs reference, KKT residuals | ≥80% at 1e-6 |
| 20 | Determinism | repeated runs | Bit-identical logs | Same answer & path |
| 21 | Industrial realism | refinery models | Plausible yields/specs; duals interpretable | Domain review (faculty/mentor) |
| 22 | From-scratch audit | binary | No solver symbols linked | CI green |
| 23 | Timeouts | hard MIPLIB | Returns incumbent + certified gap | Graceful |

---
## 25. MVP Scope (v1 — the smallest thing that honestly satisfies the PS)

| Aspect | v1 |
|---|---|
| Classes | LP, MILP, convex QP |
| Input | Free/fixed MPS, QPS; C++ API; Python API (pybind11) building models in memory |
| LP | Bounded dual simplex (Devex→DSE, BFRT + Harris, perturbation), own sparse LU (Markowitz threshold) + PFI/FT updates |
| QP | Mehrotra IPM with own quasidefinite LDLᵀ + minimum-degree ordering (also used as LP IPM) |
| MILP | B&B, best-bound + plunging, reliability branching, bound propagation, root GMI+MIR+cover, rounding/diving/feasibility pump |
| GPU | CUDA restarted PDHG for LP (own CSR SpMV, fused kernels), crossover handoff to CPU simplex; batched PDHG (k LPs) for the repeated-solve experiment |
| Sparse | CSC + CSR mirror; GPU CSR(A), CSR(Aᵀ) |
| Presolve | Rules in §20 with primal+dual postsolve |
| Robustness | Scaling, tolerances, residual checks, κ estimate, **certifier with safe bounds**, Farkas/ray checks |
| CLI | `solver model.mps --algo {auto,dual,ipm,pdhg-gpu,race} --tol ... --time ... --json out.json --cert cert.json` |
| Outputs | status, objective, x, y, reduced costs, basis, certificate, telemetry JSON |
| Diagnostics | per-phase timing, degenerate-pivot %, refactors, κ, kernel/transfer times |
| Bench | `bench/` harness running ours + HiGHS on pinned instance lists, producing tables + performance profiles |

**v2 (credibility multipliers):** router cost model trained/validated; FT update + hypersparse FTRAN/BTRAN; crossover from IPM; deterministic parallel node LPs; presolve probing & cliques; rational re-verification of final LP basis; VIPR-style MIP certificate for small instances; batched GPU safe-bound strong branching (labelled as reproduction of Blin et al./cuOpt).
**v3 (research/product-grade):** MIQP; NLP via filter IPM; conflict analysis, symmetry, restarts; GPU heuristics portfolio; HIP backend; refinery structure detection (pooling-aware SLP driver).

## 26. Advanced Scope
See v2/v3. Also: rolling-horizon planning driver with warm starts across windows; scenario-stochastic planning (two-stage LP via batched GPU across scenarios); MRPL-facing "explainability" report (binding constraints, marginal values, sensitivity ranges from the basis — ranging is cheap once you own the basis).

## 27. Research / Paper Opportunity

**Title:** *When Should a Refinery Planning Solver Use the GPU? Certified Heterogeneous Solving of Repeated LP Families.*
- **Research question:** For families of closely related planning LPs (crude valuation, SLP recursion passes, scenarios), under which (size, batch size, perturbation magnitude, required accuracy) regimes does batched GPU first-order solving with safe-bound certification and CPU simplex polish beat sequential warm-started dual simplex?
- **Hypothesis:** Warm-started dual simplex dominates for single solves and small families of medium-size LPs; batched GPU PDHG dominates above a batch-size × nnz threshold *only* when the required output accuracy is ≤1e-6 rel. and duals are post-certified; a cheap feature model predicts the regime with low regret.
- **Prior art:** PDLP/cuPDLPx/HPR-LP; Blin et al. 2026 batched PDHG (MIP subroutines); Neumaier–Shcherbina 2004; PIMS-style SLP/distributive recursion; Penalty-DR 2024; tuning FOMs 2026.
- **Gap:** no measured crossover study on refinery-planning LP families with *certified* accuracy parity (to the extent of my search).
- **Contribution (systems/empirical):** open benchmark family + generator; certified comparison protocol; validated routing model.
- **Methodology:** 3 model scales × 4 family types × batch sizes {1..512} × tol {1e-4,1e-6,1e-8}; 2 GPUs (consumer + datacenter) + 1 CPU; report SGM, regret, Nsight breakdowns.
- **Metrics:** time-to-certified-accuracy, throughput (LPs/s), dual-value error vs exact, router regret.
- **Threats to validity:** synthetic models may not reflect MRPL's real PIMS models; implementation quality (our simplex slower than HiGHS inflates GPU's relative win → **must also compare against HiGHS warm-start chain**); hardware specificity.
- **Academic value:** modest but real (workshop / applied OR journal). **MRPL value:** tells them whether buying GPUs helps *their* workload.

## 28. 10× SIH Winner Strategy

**Do NOT build:** web dashboards beyond one validation page; LLM chat; NLP/MINLP; GPU simplex; your own modeling language; 15 cut families; Docker-microservice "platform".
**Build real:** parser, presolve/postsolve, dual simplex + own LU, IPM + own LDLᵀ, B&B with reliability branching and root cuts, CUDA PDHG + batching, certifier, bench harness, refinery generator.
**Mock/stub honestly:** nothing in the solver path. UI may be minimal static HTML reading JSON.
**Must be mathematically genuine:** every status, every certificate, every speedup.
**Visible in demo:** raw CLI logs, certificate PASS/FAIL, Nsight timeline, failure table, CPU-vs-GPU crossover plot, refinery duals.
**Defer to future work (without weakening):** MIQP/NLP/MINLP, parallel B&C, HIP backend, GPU heuristics.

### Roadmap (≈10 weeks; critical path in **bold**)
| Wk | Stage | Output | Test | Failure condition |
|---|---|---|---|---|
| 1 | **Foundation**: MPS/QPS parser, model struct, CSC/CSR, certifier v0, bench harness running HiGHS | reads all Netlib | objective of known x recomputed | RANGES bugs |
| 1–3 | **Dual simplex + LU** | solves ≥60% Netlib | vs Koch optima | stalls on degenerate — add perturbation early |
| 3–4 | **Robustness pass**: scaling, Harris/BFRT, DSE, residual checks, refactor policy | ≥85% Netlib | + adversarial | — |
| 2–4 | GPU PDHG (parallel track) on NVIDIA box | Mittelmann large LPs to 1e-4 | vs cuPDLPx objective | no NVIDIA access → secure by week 1 |
| 4–5 | Presolve + postsolve | reductions + certify | round-trip tests | dual postsolve bugs |
| 4–6 | **MILP B&B**, propagation, reliability branching, heuristics | MIPLIB easy subset | vs published optima | node LP warm start broken |
| 5–6 | IPM + LDLᵀ (LP & QP) | Maros–Mészáros | vs reference | ordering fill |
| 6–7 | Root cuts; safe bounds; crossover PDHG→simplex | ablations | root gap drop | numerically bad cuts |
| 7–8 | Router + batched PDHG + refinery families | crossover plots | held-out regret | too few instances |
| 8–9 | Full benchmark runs, failure table, performance profiles | frozen tables | reproducible script | last-minute regressions → freeze code wk 9 |
| 9–10 | Demo, deck, rehearsal, mentor review | 7-min demo | dry runs | — |

### Team (4–5)
| Role | Owns | Depends on |
|---|---|---|
| P1 LP/numerics lead | dual simplex, LU, scaling, tolerances | parser |
| P2 MILP lead | B&B, branching, cuts, heuristics, propagation | P1's warm-start API (week 3) |
| P3 GPU/HPC | PDHG kernels, batching, Nsight profiling, router features | parser + CSR; NVIDIA hardware |
| P4 IPM/QP + presolve | LDLᵀ, ordering, IPM, presolve/postsolve | model struct |
| P5 Bench/industrial/validation (may be shared) | harness, certifier, refinery generators, adversarial suite, deck | all |

## 29. 10× SIH Evaluator / Red-Team Review — hostile questions with honest answers

1. **"How do I know you didn't wrap HiGHS?"** — Dependency audit (`otool -L`/`ldd`, CI allow-list), SBOM, commit history, and live: step through `ratio_test.cpp` with a breakpoint on a Netlib instance; iteration logs differ from HiGHS's.
2. **"Is the GPU real?"** — Nsight Systems timeline live, kernel names are ours; CPU and GPU PDHG produce the same iterates (up to FP reduction order).
3. **"Where exactly is the speedup?"** — Only in large-LP PDHG and batched families; show the crossover plot including where GPU *loses*.
4. **"GPU PDHG vs *what* baseline?"** — Vs our CPU PDHG (same algorithm), vs our simplex, vs HiGHS simplex/IPM. Speedup vs our own weak CPU code alone would be misleading; we show all.
5. **"PDHG to 1e-4 isn't 'optimal'."** — Correct. We report accuracy with every time; "optimal" only after crossover/certification at 1e-6/1e-9.
6. **"How do you know the answer is correct?"** — Independent certifier on the original model; safe dual bound; rational re-check on small ones.
7. **"Show a failure."** — Failure table: e.g., instances where we time out or where the certifier flagged our own wrong answers during development (and what fixed them).
8. **"What happens on degenerate LPs?"** — Perturbation + BFRT + Harris; show `degen3` telemetry with/without perturbation.
9. **"Ill-conditioned?"** — Scaled-badly suite; κ estimates; residuals before/after; we flag rather than lie.
10. **"Infeasible model?"** — Farkas ray printed and verified.
11. **"Unbounded?"** — Ray verified.
12. **"Timeout?"** — Returns incumbent + certified gap, status TIME_LIMIT.
13. **"How far behind Gurobi are you?"** — Orders of magnitude on hard MIPs. Honest numbers on the easy subset; we don't claim parity.
14. **"Why not just use HiGHS?"** — (see §7b): licence isn't the point; control, certification, workload specialization, GPU routing, and building national capability.
15. **"Why not Gurobi/CPLEX?"** — (see §7c).
16. **"What part is original?"** — All solver code is ours (engineering); original *contributions*: certified-everywhere design including GPU FOMs, validated routing model, refinery repeated-solve study.
17. **"Isn't batched PDHG already published?"** — Yes (Blin et al. 2026; cuOpt 26.02). We reproduce it and apply it to refinery families with certification; we don't claim invention.
18. **"How does it scale?"** — log–log runtime vs nnz plots; memory; up to X rows (report the largest solved and the largest failed).
19. **"Is it deterministic?"** — Yes, single-thread and deterministic GPU reductions (fixed order); show two identical logs.
20. **"Can we reproduce your benchmarks?"** — One command, pinned commit, instance hashes.
21. **"Did you choose instances favorably?"** — Inclusion rules committed before runs (git timestamp); all Netlib included.
22. **"Industrial data is fake."** — Synthetic from published structures; labelled; we'd validate on MRPL models under NDA.
23. **"Why does MRPL care about 1e-9?"** — Planners use duals (marginal values); inaccurate duals → wrong crude purchase decisions. Certified duals matter.
24. **"Does your MILP actually use cuts?"** — Ablation: root gap closed with/without.
25. **"Your B&B is just DFS."** — Show node selection + reliability branching ablation.
26. **"Why FP64 on consumer GPUs?"** — PDHG is bandwidth-bound; FP64 costs ~2× bytes not 32× flops; we measured FP32-storage variant too.
27. **"What if transfer dominates?"** — We include it; for single small LPs it does, router picks CPU.
28. **"Does cuSPARSE count as a solver?"** — We don't use it in the core; own SpMV (and we show ours vs cuSPARSE in bench).
29. **"What if the GPU isn't present?"** — Same binary falls back to CPU; demo it by hiding the device (`CUDA_VISIBLE_DEVICES=`).
30. **"How is QP convexity checked?"** — LDLᵀ inertia of Q+εI; nonconvex refused with evidence.
31. **"Presolve correctness?"** — Postsolve round-trip certified on all instances.
32. **"Numerical tolerances: scaled or unscaled?"** — Unscaled final check, both reported.
33. **"How many Netlib problems fail?"** — Exact number, listed.
34. **"What's your MIPLIB solved count vs HiGHS?"** — Exact table, same limits, same machine.
35. **"Performance variability?"** — 3 seeds, report spread.
36. **"How is the router validated?"** — Train/test split by instance family, regret vs oracle.
37. **"Wasn't a similar project on GitHub?"** — Yes; we know them; our differentiation is X (certification depth, routing validation, repeated-solve), with numbers.
38. **"Is 'sovereign' just a label?"** — We define it as properties: source ownership, zero proprietary deps, reproducible builds, auditable certificates, hardware-portable design.
39. **"Can MRPL extend it?"** — Plugin interfaces for branching rules, cuts, heuristics, presolve rules; show a 30-line custom branching rule.
40. **"What would you do with 1 more year?"** — Parallel B&C, MIQP, HIP backend, real MRPL model validation.
41. **"Maintained after SIH?"** — Honest: needs institutional owner; propose academic/MRPL partnership.
42. **"What's the largest problem you solved to certified optimality?"** — Report it.

**Credibility-destroying question:** "Your GPU speedup is against your own slow CPU code, isn't it?" → pre-empt by always including HiGHS on the same plot.

## 30. Demo Strategy (7 minutes)

| t | Screen | Action | Proves |
|---|---|---|---|
| 0:00 | Terminal | `otool -L build/solver` / `ldd`; `cloc src/` by module | Own core, no solver deps |
| 0:40 | Terminal | `solver netlib/pilot87.mps --algo dual --json` — live log: presolve reductions, scaling ranges, perturbation, refactors | Real algorithm on real benchmark |
| 1:40 | Terminal | `certify pilot87.mps out.json` → PASS table (primal viol, dual viol, safe bound, gap) | Mathematical validity |
| 2:20 | Terminal | Adversarial: badly-scaled model — naïve mode claims OPTIMAL, certifier FAIL; robust mode PASS | Robustness is measured |
| 3:00 | Terminal | Netlib infeasible model → Farkas ray verified | Infeasibility handled |
| 3:30 | Terminal + Nsight | Large LP: `--algo pdhg-gpu` vs `--algo dual`, then `--algo auto` (router explains its choice); Nsight timeline | GPU real & honest; fallback |
| 4:30 | Plot | CPU/GPU crossover chart (incl. HiGHS), batched family throughput | Measured benefit & non-benefit |
| 5:15 | Terminal | Refinery MILP (crude unloading) — B&B log: nodes, gap, cuts; duals of planning LP as marginal crude values | Industrial credibility |
| 6:00 | Table | Benchmark summary + **failure table** | Honesty |
| 6:30 | Diagram | Architecture + extension points | Extensible foundation |

**Minimal frontend:** a single static HTML page reading run JSONs: Problem (m, n, nnz, class, coefficient ranges) · Solver (presolve reduction, iterations, nodes, gap, objective, status) · Hardware (CPU time, GPU kernel time, transfer time, idle) · Validation (max residuals, integrality, safe bound, reference objective, diff, PASS/WARN/FAIL).

## 31. Risk Register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Dual simplex numerically fragile late | High | High | Build certifier + adversarial tests in week 1; perturbation early |
| No NVIDIA GPU access (dev machine is a Mac) | Medium | High | Secure lab/cloud in week 1; keep CPU PDHG reference |
| MILP too slow for any MIPLIB | Medium | Medium | Pre-registered easy subset; emphasize certified gap on timeouts |
| Scope creep (NLP, dashboards) | High | High | This document's "do not build" list |
| Competitor has stronger numbers | High | Medium | Differentiate on proof + routing + repeated-solve; acknowledge them |
| "From scratch" interpretation dispute | Medium | High | Written clarification from mentor; no sparse direct solver libs |
| Benchmark accusations of cherry-picking | Medium | High | Pre-registered inclusion rules, full tables |
| Finale 36h environment lacks GPU | Medium | Medium | Recorded Nsight traces + pre-run logs as backup, plus live CPU runs |
| Refinery models judged unrealistic | Medium | Medium | Cite literature for every structure; get a chemical-engg faculty review |

## 32. Final Technical Thesis + Claim → Evidence → Experiment

**Thesis:** *A from-scratch LP/MILP/QP core whose every result is independently certifiable, that uses the GPU only where a validated cost model predicts benefit, and that is engineered for the repeated-solve structure of refinery planning.*

| Claim | Evidence | Experiment |
|---|---|---|
| From scratch | Dependency audit, SBOM, code walkthrough | CI link check; live `otool` |
| Correct | Certifier PASS; matches published optima | Netlib/MIPLIB/MM vs Koch/MIPLIB/reference values |
| Numerically robust | Adversarial suite residuals; no false OPTIMAL | Exp. 5–7, 18 |
| Handles infeasible/unbounded | Verified rays | Exp. 3–4 |
| GPU gives measurable benefit where predicted | Crossover plot; Nsight breakdown | Exp. 10–13 |
| Adaptive fallback | Router regret | Exp. 11, 36 |
| MILP engine works | Solved count, gaps, ablations | Exp. 14–17 |
| Industrial relevance | Literature-based models, interpretable duals | Exp. 21 |
| Competitive? | **Only** "solves X/Y of set Z vs HiGHS X'/Y'" | Exp. 8, 14 |
| Reproducible | One-command bench | Fresh-clone rerun |

## 33. Open Questions Requiring Further Research
1. MRPL's actual stack (PIMS? GRTMPS? in-house?) and typical model sizes — ask the mentor; it decides whether GPU matters at all.
2. MRPL/SIH interpretation of "from scratch" re: BLAS/cuSPARSE/cuDSS — get it in writing.
3. Exact solved counts of HiGHS on your pre-registered MIPLIB subset on your hardware.
4. Whether batched-PDHG + safe bounds beats warm simplex for *any* realistic refinery family (C1) — unknown until measured.
5. Accuracy needed for planning duals in practice (1e-6? 1e-9?).
6. Whether SANKHYA / others' published claims replicate (don't cite their numbers as fact in your deck).
7. Gurobi/CPLEX academic licence terms on publishing benchmark comparisons.
8. Whether Mittelmann's current PDLP benchmark set is appropriate for your GPU class.

---

# Final Synthesis

**A. The actual problem.** Build and *prove* a trustworthy, owned LP/MILP/QP core for degenerate, badly scaled, repeatedly solved industrial models; use the GPU only where it measurably pays.

**B. Hardest 10 technical problems.** (1) Stable sparse LU + updates; (2) degeneracy/stalling in simplex; (3) scaling & tolerance semantics; (4) presolve with correct dual postsolve; (5) warm-started node LPs in B&B; (6) branching quality (reliability); (7) numerically safe cuts; (8) IPM KKT factorization with ordering/regularization; (9) PDHG tail convergence & accuracy/crossover; (10) honest benchmarking & certification.

**C. 10 most important existing approaches.** Bounded dual simplex with DSE+BFRT (Koberstein; Forrest–Goldfarb); hypersparse LU (Hall–McKinnon); Mehrotra IPM; PDLP/cuPDLPx/HPR-LP; reliability branching (Achterberg–Koch–Martin); GMI/MIR cuts; feasibility pump / fix-and-propagate; safe bounds (Neumaier–Shcherbina); iterative refinement / exact MIP (SoPlex, exact SCIP, VIPR); GPU heuristics & batched PDHG (cuOpt, Blin et al., CHAP).

**D. 5 strongest differentiation directions.** (1) Certified outputs everywhere incl. GPU; (2) validated CPU↔GPU routing model; (3) refinery repeated-solve engine (warm-start chains + batched GPU + certified duals); (4) numerical-health telemetry/explainability; (5) extension points MRPL can use (custom branching/cuts/presolve).

**E. 3 most credible architectures.** C+F (heterogeneous by role + certifier) ← recommended; D (adaptive router on top); E (industrial repeated-solve specialization).

**F. Strongest hackathon direction.** CPU dual simplex/IPM/B&B core + GPU PDHG (single + batched) + certifier + router, demonstrated on Netlib/MIPLIB-easy/Maros–Mészáros/refinery families with full failure tables vs HiGHS.

**G. Strongest research candidate.** "When should a refinery planning solver use the GPU? Certified heterogeneous solving of repeated LP families" (Category C: systems/empirical; not novel algorithmically).

**H. 10 hardest objections.** Wrapper? GPU real? Speedup vs weak baseline? PDHG accuracy? Correctness proof? Degenerate/ill-conditioned behavior? Gap to Gurobi? Why not HiGHS? Similar GitHub projects? Fake industrial data?

**I. Must-run experiments.** Netlib full correctness; infeasible set certificates; adversarial scaling suite; MIPLIB pre-registered subset vs HiGHS; Maros–Mészáros; CPU/GPU crossover incl. HiGHS; batched family throughput; router held-out regret; cut/branching ablations; determinism & reproducibility rerun.

**J. Minimum defensible prototype.** MPS reader + presolve/postsolve + dual simplex (own LU) + B&B (reliability, propagation, root GMI) + IPM for QP + CUDA PDHG + certifier + bench harness vs HiGHS.

**K. Claims to make.** "Our solver core is implemented by us from mathematical foundations; dependency audit attached." "Every reported result is independently certified; here's the failure table." "Solves X/94 Netlib, Y/Z MIPLIB-easy, W/138 Maros–Mészáros within stated tolerances vs HiGHS on the same machine." "GPU is faster for LPs above ~N nnz / batches above k, slower below; our router chooses with R% regret."

**L. Claims NOT to make.** "Replacement for CPLEX/Gurobi." "Faster than commercial solvers." "Novel GPU LP algorithm." "First GPU MIP solver." "First Indian solver" (other SIH teams exist; you can't prove first). "Solves millions of variables" (unless certified, and say which class/tolerance). "AI-powered." Any speedup without baseline, accuracy and transfer time. "Optimal" for PDHG output without crossover/certificate.

**M. Name & thesis.** **PRAMĀṆA** (Sanskrit: *valid means of knowledge / proof*) — *"A sovereign optimization core that proves its answers: from-scratch LP/MILP/QP, certified results, and GPU acceleration exactly where measurement says it pays."* (Check name availability; "Pramana" is used by some unrelated products.)

---

## Key Sources
- PS text (MRPL, SIH26119) — provided by user.
- Competing SIH26119 repos: [SANKHYA](https://github.com/thegoodengineers/SANKHYA), [SANKHYA #519 safe bounds](https://github.com/thegoodengineers/SANKHYA/issues/519), [SANKHYA PR #682](https://github.com/thegoodengineers/SANKHYA/pull/682), [PIPEPYE](https://github.com/Satyanshgaur/PIPEPYE), [YUKTI](https://github.com/Deepak-Kambala/YUKTI), [INDRA](https://github.com/ashrayr03/INDRA-SIH2026), [NOMOS](https://github.com/Mage-100/nomos), [SIH26119 (prakyath108)](https://github.com/prakyath108/SIH26119)
- [Gurobi 13.0 release notes](https://docs.gurobi.com/projects/optimizer/en/current/reference/releasenotes/changes.html) · [Gurobi 13 FAQ](https://www.gurobi.com/resources/faq/gurobi-13-0)
- [FICO: GPU acceleration of hybrid gradient in Xpress](https://www.fico.com/blogs/gpu-acceleration-hybrid-gradient-algorithm-fico-xpress) · [GAMS + cuOpt](https://www.gams.com/blog/2025/09/gpu-accelerated-optimization-with-gams-and-nvidia-cuopt/)
- [NVIDIA cuOpt 26.04 release notes](https://archive.docs.nvidia.com/cuopt/user-guide/26.04.00/release-notes.html)
- Applegate et al., PDLP ([arXiv 2106.04756](https://arxiv.org/abs/2106.04756)); Lu & Yang, [cuPDLP.jl](https://arxiv.org/pdf/2311.12180); [cuPDLPx](https://arxiv.org/abs/2507.14051); [HPR-LP](https://link.springer.com/article/10.1007/s12532-025-00292-0); [GPU FOM survey](https://arxiv.org/pdf/2506.02174); [relationships among GPU FOMs](https://arxiv.org/pdf/2509.23903)
- Çördük et al., [GPU-Accelerated Primal Heuristics for MIP](https://arxiv.org/abs/2510.20499) · Kempke & Koch, [Fix-and-propagate with low-precision FO LP](https://arxiv.org/abs/2503.10344) · Tjusila et al., [CHAP](https://arxiv.org/abs/2605.05086)
- Blin, Gualandi, Maes, Lodi, Stellato, [Batched First-Order Methods for Parallel LP Solving in MIP](https://arxiv.org/abs/2601.21990)
- Guan et al., [B³-PWL](https://arxiv.org/abs/2608.28988) · [Batched B&B for k-sparse GLMs](https://arxiv.org/abs/2605.22188) · [Tuning GPU LP](https://arxiv.org/pdf/2606.08638)
- Cederberg & Boyd, [Presolving for GPU FO LP solvers (PSLP)](https://arxiv.org/abs/2604.23951), [GPU-accelerated presolving (cuPSLP)](https://arxiv.org/abs/2609.16182)
- Neumaier & Shcherbina, [Safe bounds in LP and MILP](https://link.springer.com/article/10.1007/s10107-003-0433-3), *Math. Prog.* 99 (2004)
- Gurung & Ray, [Batched LPs on GPU](https://arxiv.org/pdf/1802.08557); GPU simplex: [1803.04378](https://arxiv.org/pdf/1803.04378), [2211.10979](https://arxiv.org/pdf/2211.10979)
- Zhang, …, Dai, [Distributed Recursion Revisited](https://arxiv.org/abs/2411.09554) · [Aspen PIMS / RPMS use of CPLEX/Xpress (patent)](https://patents.google.com/patent/US20030097243A1/en) · [Refinery planning/scheduling review](https://journal.hep.com.cn/fem/EN/10.1007/s42524-020-0123-3)
- [COPT user guide](https://arxiv.org/pdf/2208.14314) · [SIH 2026 timeline](https://techpathdaily.com/sih-2026-official-timeline-is-here/)
- Cited from knowledge (verify page numbers before quoting in the deck): Koberstein (2005) PhD thesis on dual simplex; Forrest & Goldfarb (1992) steepest edge; Hall & McKinnon (2005) hypersparsity; Huangfu & Hall (2018) *Math. Prog. Comp.*; Achterberg, Koch & Martin (2005) *OR Letters*; Achterberg & Wunderling (2013); Gleixner et al. MIPLIB 2017 (*Math. Prog. Comp.* 2021); Gleixner, Steffy & Wolter (2016) iterative refinement; Eifler & Gleixner exact MIP; Cheung, Gleixner & Steffy (2017) VIPR; Sofranac, Gleixner & Pokutta (2022) GPU propagation; Lee, Pinto, Grossmann & Park (1996); Pinto, Joly & Moro (2000); Morales-España et al. (2013); Stellato et al. OSQP (2020); Maros & Mészáros (1999); Furini et al. QPLIB (2019); Klotz & Newman (2013); Dolan & Moré (2002); Lodi & Tramontani (2013).
