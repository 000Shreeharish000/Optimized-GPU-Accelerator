# Algorithms and references

Every algorithm is implemented in this repository; the references are the mathematics followed.

| component | what is implemented | references |
|---|---|---|
| Sparse LU (`lp/lu.cpp`) | Markowitz pivoting with threshold partial pivoting (u = 0.1) and count buckets; column/row singletons first; rank-deficiency report; product-form (PFI) updates; residual-triggered refactorization | Markowitz (1957); Suhl & Suhl (1990); Duff, Erisman & Reid, *Direct Methods for Sparse Matrices* |
| Dual simplex (`lp/simplex.cpp`) | bounded dual simplex on `[A I]`, dual steepest edge with exact reference weights, bound-flipping ratio test + Harris two-pass, cost perturbation, cost shifting, dual phase 1 by the artificial-bounding subproblem, stall detection with re-perturbation, 100× tighter dual polish, per-variable scaled/unscaled tolerances, rigorous Farkas check before any infeasibility claim | Koberstein (2005) PhD thesis; Forrest & Goldfarb (1992); Harris (1973); Maros (2003) |
| Primal simplex | Devex pricing, Harris ratio test with bound flips, primal ray for unboundedness | Harris (1973); Forrest & Goldfarb (1992) |
| Scaling (`lp/scaling.cpp`) | geometric mean + equilibration (simplex/IPM), Ruiz + Pock–Chambolle (PDHG); powers of two (exact) | Curtis & Reid (1972); Ruiz (2001); Pock & Chambolle (2011) |
| Presolve (`presolve/presolve.cpp`) | singleton/empty/redundant/forcing rows, fixed/empty/dominated columns, implied-free column singletons; postsolve of x, y, d and basis | Andersen & Andersen (1995); Achterberg et al. (2020) |
| Ordering (`ipm/ldl.cpp`) | quotient-graph approximate minimum degree (element absorption, external-degree bounds) | Amestoy, Davis & Duff (1996) |
| LDLᵀ | elimination tree + column counts, up-looking numeric factorization, dynamic regularization of quasidefinite systems, iterative refinement | Davis (2005) LDL; Vanderbei (1995) quasidefinite matrices |
| IPM (`ipm/ipm.cpp`) | Mehrotra predictor–corrector on the regularized augmented system, general bounds, fixed-column elimination, centered start, centrality safeguard, best-iterate recovery, crossover | Mehrotra (1992); Wright (1997); Gondzio (1996) |
| Convexity check | LDLᵀ inertia; certificate v = Pᵀ L⁻ᵀ eₖ with vᵀQv = dₖ < 0 verified by multiplication | — |
| PDHG (`pdhg/*`) | restarted reflected Halpern PDHG, constant step 0.998/‖A‖₂ (power iteration), primal weight updates at restarts, fixed-point-residual restarts, KKT termination on the unscaled problem; CPU and CUDA backends; batched SpMM | Applegate et al. (2021) PDLP; Lu & Yang (2023) cuPDLP; Lu, Peng & Yang (2025) cuPDLPx; Chambolle & Pock (2011) |
| Branch-and-cut (`mip/*`) | warm-started node LPs, reliability branching, best bound + plunging, activity propagation, root GMI/c-MIR/cover cuts with efficacy/parallelism selection, rounding, diving, feasibility pump, fix-and-resolve, safe node bounds | Achterberg, Koch & Martin (2005); Achterberg (2007); Marchand & Wolsey (2001); Fischetti, Glover & Lodi (2005) |
| Safe bounds & certificates (`cert/*`) | Neumaier–Shcherbina bound with outward rounding (TwoSum / FMA TwoProduct), sign-cone projection of y, backward-error reporting, implied-bound fallback for inexact duals, Farkas / ray / QP linearization checks | Neumaier & Shcherbina (2004); Ogita, Rump & Oishi (2005) |
| Parametric analysis (`parametric/*`) | tangent-intersection enumeration of the breakpoints of the optimal value function, every oracle point certified | Gal (1995) *Postoptimal Analyses, Parametric Programming*; Eichfelder (2008) sandwich methods |
| Exact verifier (`python/pramana/verify.py`) | own MPS parser into exact rationals, exact feasibility / safe bound / Farkas / ray checks, exact rational vertex (basis) verification | Applegate et al. QSopt_ex (2007); Gleixner, Steffy & Wolter (2016) |
