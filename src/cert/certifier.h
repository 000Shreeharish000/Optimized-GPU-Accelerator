// Independent result certifier ("PRAMANA layer").
//
// Receives only the ORIGINAL model and the claimed (status, x, y, rays) — never
// solver internals — and re-derives every claim:
//
//   OPTIMAL (LP)     primal feasibility on the unscaled model; rigorous
//                    Neumaier–Shcherbina dual bound  LB(y) <= min c^T x  computed
//                    with outward rounding (valid for ANY y, so it also certifies
//                    low-accuracy GPU/PDHG duals); gap = c^T x - LB(y).
//   INFEASIBLE       Farkas certificate: for multipliers y, the interval hull of
//                    (A^T y)^T x over the column box is disjoint from y^T r over
//                    the row box  =>  no x satisfies both.
//   UNBOUNDED        primal ray d in the recession cone with c^T d < 0, plus a
//                    feasible point.
//   OPTIMAL (MIP)    incumbent feasibility + integrality; the reported bound
//                    must come from safe node bounds (mip/bnb.cpp).
//   OPTIMAL (QP)     KKT residuals + linearization bound
//                    f(x) >= f(x^) + g^T (x - x^), g = c + Q x^, then the LP safe
//                    bound with cost g (valid because convexity was certified).
//
// A claim that fails a check is DOWNGRADED (Optimal -> NumericalFailure etc.).
#pragma once

#include <string>
#include <vector>

#include "core/model.h"
#include "util/json.h"

namespace pramana {

enum class Verdict { Pass, Warn, Fail, NotApplicable };
const char* verdictName(Verdict v);

struct CertCheck {
  std::string name;
  Verdict verdict = Verdict::NotApplicable;
  double value = 0;
  double tolerance = 0;
  std::string detail;
};

struct CertTolerances {
  double primalFeas = 1e-6;   // absolute + relative (scaled by 1+|bound|)
  double dualFeas = 1e-6;
  double relGap = 1e-6;       // LP optimality gap (relative)
  double integrality = 1e-5;
  double mipRelGap = 1e-4;    // accepted MIP gap for OPTIMAL
};

struct Certificate {
  Status claimed = Status::NotSolved;
  Status certified = Status::NotSolved;
  bool accepted = false;
  std::vector<CertCheck> checks;
  double primalObjective = 0;
  double safeDualBound = -kInf;   // rigorous lower bound on the optimum (min form, model sense)
  double certifiedGap = kInf;     // relative gap between primal objective and safe bound
  double maxPrimalViolation = 0;
  double maxDualViolation = 0;
  double maxIntegralityViolation = 0;
  double costPerturbation = 0;    // backward error of the safe bound (relative cost change)
  std::string summary() const;
  Json toJson() const;
};

struct SolutionClaim {
  Status status = Status::NotSolved;
  std::vector<double> x;          // column values
  std::vector<double> rowDual;    // y (min form): reduced cost d = c - A^T y
  std::vector<double> farkas;     // row multipliers (infeasibility)
  std::vector<double> primalRay;  // column direction (unboundedness)
  double bestBound = -kInf;       // MIP: claimed global lower bound (min form)
  bool boundIsSafe = false;       // MIP: bound computed from safe node bounds
};

// `model` is the original model in its own sense; claims are interpreted in the
// model's sense (duals of the minimization form are converted internally).
Certificate certify(const Model& model, const SolutionClaim& claim, const CertTolerances& tol = {});

// ---- Rigorous building blocks (also used by B&B for safe node bounds) ----
struct SafeBoundResult {
  double bound = -kInf;       // rigorous lower bound on min c^T x (+1/2 x'Qx not included)
  int infiniteTerms = 0;      // columns that forced -inf
  double maxDualInfeasibility = 0;
  double rowProjection = 0;   // max change from projecting y onto its sign cone (exact, free)
  // Wrong-sign reduced costs on infinite bounds are clipped; the bound is then
  // rigorous for a cost vector perturbed by at most this (relative) amount.
  double costPerturbation = 0;
  double maxAllowedPerturbation = 1e-9;
};
// Lower bound on  min { c^T x : L <= Ax <= U, l <= x <= u }  valid for any y.
SafeBoundResult safeDualBound(const SparseMatrix& A, const std::vector<double>& cost,
                              const std::vector<double>& colLower, const std::vector<double>& colUpper,
                              const std::vector<double>& rowLower, const std::vector<double>& rowUpper,
                              const std::vector<double>& y);
// True if y proves infeasibility of the row/column boxes.
bool verifyFarkas(const SparseMatrix& A, const std::vector<double>& colLower, const std::vector<double>& colUpper,
                  const std::vector<double>& rowLower, const std::vector<double>& rowUpper,
                  const std::vector<double>& y, double* margin = nullptr, double* perturbation = nullptr);

}  // namespace pramana
