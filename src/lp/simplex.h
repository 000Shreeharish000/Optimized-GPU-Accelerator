// Bounded revised simplex engine (dual simplex main algorithm + primal
// simplex cleanup), written from scratch.
//
// Internal form: variables 0..n-1 are structural, n..n+m-1 are row logicals
// with column e_i, so that   [A I] (x; s) = 0,   s_i = -a_i^T x in [-U_i, -L_i].
//
// Algorithmic components (see docs/ALGORITHMS.md for references):
//   * dual phase 1 by the artificial-bounding subproblem (all variables boxed
//     in {[0,0],[0,1],[-1,0],[-1000,1000]}), solved by dual phase 2 itself;
//   * dual phase 2 with dual steepest-edge pricing (exact reference weights),
//     bound-flipping ("long-step") ratio test combined with Harris' two-pass
//     tolerance, cost perturbation against dual degeneracy, cost shifting for
//     residual dual infeasibilities;
//   * primal phase 2 with Devex pricing and Harris ratio test (removes the
//     perturbation/shift residue, detects unboundedness with a primal ray);
//   * own sparse LU (lp/lu.h) with PFI updates, residual-triggered
//     refactorization and rank-deficiency repair with logicals;
//   * Farkas certificate for infeasibility, primal ray for unboundedness;
//   * warm start from any basis after bound / cost / row changes (B&B, cuts,
//     parametric analysis).
#pragma once

#include <memory>
#include <vector>

#include "core/model.h"
#include "lp/lu.h"
#include "util/common.h"
#include "util/rng.h"
#include "util/timer.h"

namespace pramana {

struct SimplexOptions {
  double primalFeasTol = 1e-7;
  double dualFeasTol = 1e-7;
  double pivotTol = 1e-9;          // min |alpha| accepted in ratio tests
  long long maxIterations = 1LL << 60;
  int refactorInterval = 100;
  bool perturb = true;
  bool scale = true;
  bool useBfrt = true;
  uint64_t seed = 12345;
  int logInterval = 0;             // 0 = no iteration log
  double objectiveCutoff = kInf;   // (min form) stop early when proven worse
};

struct SimplexStats {
  long long iterations = 0;
  long long dualPhase1Iterations = 0;
  long long dualPhase2Iterations = 0;
  long long primalIterations = 0;
  long long degeneratePivots = 0;
  long long boundFlips = 0;
  long long refactorizations = 0;
  long long numericalRefactors = 0;
  long long rejectedPivots = 0;
  long long rejectedRays = 0;
  long long polishIterations = 0;  // successful dual-polish passes
  long long stallRecoveries = 0;   // re-perturbations triggered by the stall detector      // dual rays that failed the rigorous Farkas check
  long long singularRepairs = 0;
  long long costShifts = 0;
  long long perturbations = 0;
  double maxPerturbation = 0;
  double maxCostShift = 0;
  double scalingRangeBefore = 0, scalingRangeAfter = 0;
  double lastConditionEstimate = 0;
  double timeSeconds = 0;
};

class Simplex {
 public:
  // `lp` is interpreted as a minimization LP (Q and integrality are ignored).
  Simplex(const Model& lp, const SimplexOptions& opts = {});

  Status solve(const Deadline* deadline = nullptr);

  // --- modification (all in unscaled model units) ---
  void setColBounds(int j, double lower, double upper);
  void setRowBounds(int i, double lower, double upper);
  void setColCost(int j, double cost);
  void setObjectiveOffset(double off) { offset_ = off; }
  // Appends rows; new logicals enter the basis (keeps dual feasibility).
  void addRows(const std::vector<double>& lower, const std::vector<double>& upper,
               const std::vector<int>& starts, const std::vector<int>& cols,
               const std::vector<double>& vals);
  // Removes rows whose logical is basic (e.g. slack cuts); returns count removed.
  int removeBasicRows(const std::vector<int>& rows);

  // --- basis ---
  void getBasis(std::vector<BasisStatus>& colStatus, std::vector<BasisStatus>& rowStatus) const;
  void setBasis(const std::vector<BasisStatus>& colStatus, const std::vector<BasisStatus>& rowStatus);
  bool hasBasis() const { return haveBasis_; }

  // --- results (unscaled, min form incl. offset) ---
  double objective() const;
  std::vector<double> colValues() const;
  std::vector<double> rowActivities() const;
  std::vector<double> rowDuals() const;      // y: d = c - A^T y
  std::vector<double> reducedCosts() const;
  const std::vector<double>& farkasRay() const { return farkas_; }   // row multipliers
  const std::vector<double>& primalRay() const { return primalRay_; }  // column direction

  // Row `pos` of B^{-1}[A I] in unscaled variables: x_B(pos) + sum_j t_j v_j = 0
  // over all n+m variables (structural then logical s = -Ax). Returns basic var.
  int tableauRow(int pos, std::vector<double>& coefs) const;
  int basicVariable(int pos) const { return basicIndex_[pos]; }
  // Sensitivity of the basic solution to a nonbasic variable v (unscaled):
  // x_B(pos) changes by out[pos] * (change of x_v). Logical v = n+i is s_i = -a_i^T x.
  void basicDirection(int v, std::vector<double>& out) const;
  // Unscaled variable scale (x_v = scale * internal value).
  double variableScale(int v) const { return v < n_ ? colScale_[v] : 1.0 / rowScale_[v - n_]; }
  int numRows() const { return m_; }
  int numCols() const { return n_; }
  double variableValue(int v) const;  // unscaled (logical = -row activity)

  const SimplexStats& stats() const { return stats_; }
  SimplexOptions& options() { return opt_; }

 private:
  // setup
  void initSlackBasis();
  void computeTolerances();
  void setNonbasicValue(int j, int preferredMove);
  bool factor();
  void computePrimal();
  void computeDual();
  double computeObjectiveWork() const;
  void flipToDualFeasibility(bool recomputePrimal);
  int countDualInfeasibilities(double* sum = nullptr) const;
  void shiftCostsForDualFeasibility();
  void perturbCosts(double factor = 1.0);
  void restoreCosts();

  // algorithms
  enum class PhaseResult { Optimal, Infeasible, Unbounded, Limit, Numerical, Cutoff };
  PhaseResult dualPhase2(bool phase1);
  PhaseResult primalPhase2();
  void dualPhase1();

  // helpers
  void columnOf(int j, HVector& out) const;  // [A I] column into row-indexed HVector
  void priceRow(const HVector& rho, std::vector<double>& rowAp, std::vector<int>& rowIdx);
  void updateBasis(int pos, int enter, const HVector& colAq);
  bool limitReached(const Deadline* dl);
  void logIteration(const char* phase);
  double primalInfeasibilitySum() const;
  bool farkasHolds(const std::vector<double>& yScaled) const;

  SimplexOptions opt_;
  SimplexStats stats_;
  Rng rng_;
  Timer timer_;

  int n_ = 0, m_ = 0;
  double offset_ = 0;
  SparseMatrix A_;    // scaled, CSC
  SparseMatrix AT_;   // scaled, row-wise (CSC of transpose)
  bool atValid_ = false;
  std::vector<double> colScale_, rowScale_;
  std::vector<double> cost_, costWork_;  // size n+m (scaled)
  std::vector<double> lower_, upper_;    // size n+m (scaled)
  std::vector<double> tolP_, tolD_;      // per-variable feasibility tolerances (scaled)

  std::vector<int> basicIndex_;
  std::vector<int> basicPos_;  // var -> position or -1
  std::vector<int8_t> move_;   // nonbasic: +1 at lower, -1 at upper, 0 free/fixed
  std::vector<double> value_;  // all n+m variables
  std::vector<double> dual_;   // reduced costs (nonbasic), 0 for basic
  std::vector<double> y_;      // row duals (scaled)
  std::vector<double> dseWeight_;
  std::vector<double> devex_;

  LuFactor lu_;
  bool haveBasis_ = false;
  bool factored_ = false;
  bool costsModified_ = false;  // perturbed or shifted
  const Deadline* deadline_ = nullptr;

  std::vector<double> farkas_, primalRay_;

  // work
  HVector col_, rho_, tau_, flipCol_;
  std::vector<double> rowAp_;
  std::vector<int> rowIdx_;
  std::vector<char> rowMark_;
};

// Convenience: solve an LP model with the simplex and fill a solution vector set.
struct LpResult {
  Status status = Status::NotSolved;
  double objective = 0;
  std::vector<double> x, rowActivity, rowDual, reducedCost;
  std::vector<BasisStatus> colStatus, rowStatus;
  std::vector<double> farkas, primalRay;
  SimplexStats stats;
};

LpResult solveLpSimplex(const Model& lp, const SimplexOptions& opts = {},
                        const Deadline* deadline = nullptr,
                        const std::vector<BasisStatus>* warmCol = nullptr,
                        const std::vector<BasisStatus>* warmRow = nullptr);

}  // namespace pramana
