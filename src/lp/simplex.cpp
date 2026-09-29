#include "lp/simplex.h"

#include <algorithm>
#include <cmath>

#include "cert/certifier.h"
#include "lp/scaling.h"
#include "util/log.h"

namespace pramana {

namespace {
inline bool finite(double v) { return std::fabs(v) < kInfBoundThreshold; }
}  // namespace

Simplex::Simplex(const Model& lp, const SimplexOptions& opts) : opt_(opts), rng_(opts.seed) {
  n_ = lp.numCols();
  m_ = lp.numRows();
  offset_ = lp.objOffset;
  A_ = lp.A;
  A_.numRows = m_;
  colScale_.assign(n_, 1.0);
  rowScale_.assign(m_, 1.0);
  if (opt_.scale && A_.nnz() > 0) {
    ScalingResult s = geometricScaling(A_);
    applyScaling(A_, s);
    colScale_ = s.colScale;
    rowScale_ = s.rowScale;
    stats_.scalingRangeBefore = s.rangeBefore;
    stats_.scalingRangeAfter = s.rangeAfter;
  }
  const int N = n_ + m_;
  cost_.assign(N, 0.0);
  lower_.assign(N, 0.0);
  upper_.assign(N, 0.0);
  for (int j = 0; j < n_; ++j) {
    cost_[j] = lp.colCost[j] * colScale_[j];
    lower_[j] = finite(lp.colLower[j]) ? lp.colLower[j] / colScale_[j] : -kInf;
    upper_[j] = finite(lp.colUpper[j]) ? lp.colUpper[j] / colScale_[j] : kInf;
  }
  for (int i = 0; i < m_; ++i) {
    lower_[n_ + i] = finite(lp.rowUpper[i]) ? -lp.rowUpper[i] * rowScale_[i] : -kInf;
    upper_[n_ + i] = finite(lp.rowLower[i]) ? -lp.rowLower[i] * rowScale_[i] : kInf;
  }
  costWork_ = cost_;
  computeTolerances();
  AT_ = A_.transpose();
  atValid_ = true;
  lu_.setup(m_, n_, &A_);
  col_.setup(m_);
  rho_.setup(m_);
  tau_.setup(m_);
  flipCol_.setup(m_);
  rowAp_.assign(N, 0.0);
  rowIdx_.reserve(N);
  rowMark_.assign(N, 0);
  value_.assign(N, 0.0);
  dual_.assign(N, 0.0);
  y_.assign(m_, 0.0);
  move_.assign(N, 0);
  basicIndex_.assign(m_, 0);
  basicPos_.assign(N, -1);
  dseWeight_.assign(m_, 1.0);
  devex_.assign(N, 1.0);
}

// ---------------------------------------------------------------------------
// Basis setup
// ---------------------------------------------------------------------------

// Per-variable tolerances in the scaled space chosen so that a variable within
// tolerance is ALSO within tolerance on the unscaled model (x = vs * x', d = d' / vs).
void Simplex::computeTolerances() {
  const int N = n_ + m_;
  tolP_.assign(N, opt_.primalFeasTol);
  tolD_.assign(N, opt_.dualFeasTol);
  for (int v = 0; v < N; ++v) {
    double vs = v < n_ ? colScale_[v] : 1.0 / rowScale_[v - n_];
    tolP_[v] = opt_.primalFeasTol * std::min(1.0, 1.0 / vs);
    tolD_[v] = opt_.dualFeasTol * std::min(1.0, vs);
  }
}

void Simplex::setNonbasicValue(int j, int preferredMove) {
  double lo = lower_[j], up = upper_[j];
  if (lo == up) {
    value_[j] = lo;
    move_[j] = 0;
  } else if (finite(lo) && finite(up)) {
    if (preferredMove < 0) {
      value_[j] = up;
      move_[j] = -1;
    } else {
      value_[j] = lo;
      move_[j] = 1;
    }
  } else if (finite(lo)) {
    value_[j] = lo;
    move_[j] = 1;
  } else if (finite(up)) {
    value_[j] = up;
    move_[j] = -1;
  } else {
    value_[j] = 0.0;
    move_[j] = 0;
  }
}

void Simplex::initSlackBasis() {
  const int N = n_ + m_;
  std::fill(basicPos_.begin(), basicPos_.end(), -1);
  for (int i = 0; i < m_; ++i) {
    basicIndex_[i] = n_ + i;
    basicPos_[n_ + i] = i;
    move_[n_ + i] = 0;
  }
  for (int j = 0; j < n_; ++j) setNonbasicValue(j, cost_[j] < 0 ? -1 : 1);
  std::fill(dseWeight_.begin(), dseWeight_.end(), 1.0);  // exact for B = I
  std::fill(devex_.begin(), devex_.end(), 1.0);
  (void)N;
  haveBasis_ = true;
  factored_ = false;
}

bool Simplex::factor() {
  for (int attempt = 0; attempt < 4; ++attempt) {
    int def = lu_.factorize(basicIndex_);
    stats_.refactorizations++;
    if (def == 0) {
      factored_ = true;
      return true;
    }
    // Basis repair: replace unpivoted columns by logicals of unpivoted rows.
    for (int k = 0; k < def; ++k) {
      int pos = lu_.deficientPositions[k];
      int row = lu_.deficientRows[k];
      int old = basicIndex_[pos];
      int nv = n_ + row;
      basicPos_[old] = -1;
      setNonbasicValue(old, dual_[old] < 0 ? -1 : 1);
      basicIndex_[pos] = nv;
      basicPos_[nv] = pos;
      move_[nv] = 0;
      dseWeight_[pos] = 1.0;
    }
    stats_.singularRepairs += def;
    PLOG_DETAIL("simplex: basis rank deficiency %d repaired with logicals", def);
    // A repaired basis generally invalidates DSE weights; reset to 1 (Devex-like).
    std::fill(dseWeight_.begin(), dseWeight_.end(), 1.0);
  }
  factored_ = false;
  return false;
}

void Simplex::columnOf(int j, HVector& out) const {
  out.clear();
  if (j >= n_) {
    out.array[j - n_] = 1.0;
    out.index[0] = j - n_;
    out.count = 1;
    return;
  }
  for (int k = A_.start[j]; k < A_.start[j + 1]; ++k) {
    out.array[A_.index[k]] = A_.value[k];
    out.index[out.count++] = A_.index[k];
  }
}

void Simplex::computePrimal() {
  col_.clear();
  double* rhs = col_.array.data();
  for (int j = 0; j < n_; ++j) {
    if (basicPos_[j] >= 0) continue;
    double v = value_[j];
    if (v == 0.0) continue;
    for (int k = A_.start[j]; k < A_.start[j + 1]; ++k) rhs[A_.index[k]] -= A_.value[k] * v;
  }
  for (int i = 0; i < m_; ++i)
    if (basicPos_[n_ + i] < 0) rhs[i] -= value_[n_ + i];
  lu_.ftran(col_);
  for (int p = 0; p < m_; ++p) value_[basicIndex_[p]] = col_.array[p];
  col_.clear();
}

void Simplex::computeDual() {
  rho_.clear();
  for (int p = 0; p < m_; ++p) rho_.array[p] = costWork_[basicIndex_[p]];
  lu_.btran(rho_);
  for (int i = 0; i < m_; ++i) y_[i] = rho_.array[i];
  rho_.clear();
  for (int j = 0; j < n_; ++j) {
    if (basicPos_[j] >= 0) {
      dual_[j] = 0.0;
      continue;
    }
    double s = 0;
    for (int k = A_.start[j]; k < A_.start[j + 1]; ++k) s += A_.value[k] * y_[A_.index[k]];
    dual_[j] = costWork_[j] - s;
  }
  for (int i = 0; i < m_; ++i) dual_[n_ + i] = basicPos_[n_ + i] >= 0 ? 0.0 : costWork_[n_ + i] - y_[i];
}

double Simplex::computeObjectiveWork() const {
  double f = 0;
  for (int j = 0; j < n_ + m_; ++j) f += costWork_[j] * value_[j];
  return f;
}

void Simplex::flipToDualFeasibility(bool recomputePrimal) {
  bool any = false;
  for (int j = 0; j < n_ + m_; ++j) {
    if (basicPos_[j] >= 0) continue;
    double lo = lower_[j], up = upper_[j];
    if (lo == up || !finite(lo) || !finite(up)) continue;
    const double tol = tolD_[j];
    if (move_[j] == 1 && dual_[j] < -tol) {
      value_[j] = up;
      move_[j] = -1;
      any = true;
    } else if (move_[j] == -1 && dual_[j] > tol) {
      value_[j] = lo;
      move_[j] = 1;
      any = true;
    }
  }
  if (any && recomputePrimal) computePrimal();
}

int Simplex::countDualInfeasibilities(double* sum) const {
  int cnt = 0;
  double s = 0;
  for (int j = 0; j < n_ + m_; ++j) {
    if (basicPos_[j] >= 0) continue;
    if (lower_[j] == upper_[j]) continue;
    const double tol = tolD_[j];
    double d = dual_[j];
    double inf = 0;
    if (move_[j] == 1) inf = -d;
    else if (move_[j] == -1) inf = d;
    else inf = std::fabs(d);  // free
    if (inf > tol) {
      ++cnt;
      s += inf;
    }
  }
  if (sum) *sum = s;
  return cnt;
}

void Simplex::shiftCostsForDualFeasibility() {
  for (int j = 0; j < n_ + m_; ++j) {
    if (basicPos_[j] >= 0 || lower_[j] == upper_[j]) continue;
    const double tol = tolD_[j];
    double d = dual_[j];
    double target;
    if (move_[j] == 1) {
      if (d >= -tol) continue;
      target = tol * (1.0 + rng_.uniform());
    } else if (move_[j] == -1) {
      if (d <= tol) continue;
      target = -tol * (1.0 + rng_.uniform());
    } else {
      if (std::fabs(d) <= tol) continue;
      target = 0.0;
    }
    double shift = target - d;
    costWork_[j] += shift;
    dual_[j] = target;
    costsModified_ = true;
    stats_.costShifts++;
    stats_.maxCostShift = std::max(stats_.maxCostShift, std::fabs(shift));
  }
}

void Simplex::perturbCosts(double factor) {
  // Perturb only nonbasic columns in the direction that *increases* their dual
  // feasibility margin: y is unchanged, so dual feasibility is preserved.
  double maxCost = 0;
  for (int j = 0; j < n_; ++j) maxCost = std::max(maxCost, std::fabs(cost_[j]));
  if (maxCost == 0) maxCost = 1;  // pure feasibility problem: still perturb (dual degenerate)
  const double base = factor * 5e-7 * std::min(std::max(maxCost, 1.0), 1e3);
  for (int j = 0; j < n_; ++j) {
    if (basicPos_[j] >= 0 || move_[j] == 0) continue;
    double delta = base * (1.0 + std::fabs(cost_[j])) * (1.0 + rng_.uniform());
    delta = std::min(delta, 1e-3 * (1.0 + std::fabs(cost_[j])));
    costWork_[j] += move_[j] * delta;
    dual_[j] += move_[j] * delta;
    stats_.maxPerturbation = std::max(stats_.maxPerturbation, delta);
  }
  costsModified_ = true;
  stats_.perturbations++;
}

void Simplex::restoreCosts() {
  costWork_ = cost_;
  costsModified_ = false;
  computeDual();
}

double Simplex::primalInfeasibilitySum() const {
  double s = 0;
  for (int p = 0; p < m_; ++p) {
    int v = basicIndex_[p];
    double x = value_[v];
    if (x < lower_[v] - tolP_[v]) s += lower_[v] - x;
    else if (x > upper_[v] + tolP_[v]) s += x - upper_[v];
  }
  return s;
}

bool Simplex::limitReached(const Deadline* dl) {
  if (stats_.iterations >= opt_.maxIterations) return true;
  if (dl && (stats_.iterations & 15) == 0 && dl->expired()) return true;
  return false;
}

void Simplex::logIteration(const char* phase) {
  if (opt_.logInterval <= 0 || stats_.iterations % opt_.logInterval != 0) return;
  double dsum = 0;
  countDualInfeasibilities(&dsum);
  PLOG_INFO("  %-7s it %8lld  obj % .10e  pinf %.2e  dinf %.2e  t %.2fs", phase, stats_.iterations,
            computeObjectiveWork() + offset_, primalInfeasibilitySum(), dsum, timer_.seconds());
}

// ---------------------------------------------------------------------------
// PRICE: pivot row alpha_p = rho^T [A I] for nonbasic variables
// ---------------------------------------------------------------------------
void Simplex::priceRow(const HVector& rho, std::vector<double>& rowAp, std::vector<int>& rowIdx) {
  for (int j : rowIdx) rowAp[j] = 0.0, rowMark_[j] = 0;
  rowIdx.clear();
  const double density = static_cast<double>(rho.count) / std::max(1, m_);
  if (density < 0.1) {
    for (int k = 0; k < rho.count; ++k) {
      int i = rho.index[k];
      double r = rho.array[i];
      for (int e = AT_.start[i]; e < AT_.start[i + 1]; ++e) {
        int j = AT_.index[e];
        if (!rowMark_[j]) {
          rowMark_[j] = 1;
          rowIdx.push_back(j);
        }
        rowAp[j] += r * AT_.value[e];
      }
      int lj = n_ + i;
      if (!rowMark_[lj]) {
        rowMark_[lj] = 1;
        rowIdx.push_back(lj);
      }
      rowAp[lj] += r;
    }
  } else {
    const double* r = rho.array.data();
    for (int j = 0; j < n_; ++j) {
      if (basicPos_[j] >= 0) continue;
      double s = 0;
      for (int k = A_.start[j]; k < A_.start[j + 1]; ++k) s += A_.value[k] * r[A_.index[k]];
      if (s != 0.0) {
        rowAp[j] = s;
        rowMark_[j] = 1;
        rowIdx.push_back(j);
      }
    }
    for (int i = 0; i < m_; ++i) {
      if (r[i] == 0.0 || basicPos_[n_ + i] >= 0) continue;
      rowAp[n_ + i] = r[i];
      rowMark_[n_ + i] = 1;
      rowIdx.push_back(n_ + i);
    }
  }
}

bool Simplex::farkasHolds(const std::vector<double>& y) const {
  std::vector<double> cl(lower_.begin(), lower_.begin() + n_), cu(upper_.begin(), upper_.begin() + n_);
  std::vector<double> rl(m_), ru(m_);
  for (int i = 0; i < m_; ++i) {
    rl[i] = -upper_[n_ + i];
    ru[i] = -lower_[n_ + i];
  }
  double margin = 0, pert = 0;
  bool ok = verifyFarkas(A_, cl, cu, rl, ru, y, &margin, &pert);
  if (!ok) PLOG_DEBUG("simplex: Farkas ray rejected (margin %.3e, backward error %.3e)", margin, pert);
  return ok;
}

void Simplex::updateBasis(int pos, int enter, const HVector& colAq) {
  int leave = basicIndex_[pos];
  lu_.update(colAq, pos);
  basicIndex_[pos] = enter;
  basicPos_[enter] = pos;
  basicPos_[leave] = -1;
  move_[enter] = 0;
  dual_[enter] = 0.0;
}

// ---------------------------------------------------------------------------
// Dual simplex phase 2 (also used for the phase-1 subproblem)
// ---------------------------------------------------------------------------
Simplex::PhaseResult Simplex::dualPhase2(bool phase1) {
  const double tolP = opt_.primalFeasTol;
  const double tolD = opt_.dualFeasTol;
  std::vector<char> taboo(m_, 0);
  int tabooCount = 0;
  int allTabooEvents = 0;
  double stallBest = -kInf, stallFactor = 1.0;
  long long stallIter = stats_.iterations;
  struct Cand {
    int j;
    double ratio, harris, absAlpha;
  };
  std::vector<Cand> cands;
  std::vector<int> flips, group;

  auto reinvert = [&]() {
    factor();
    computePrimal();
    computeDual();
    flipToDualFeasibility(true);
    if (countDualInfeasibilities() > 0) shiftCostsForDualFeasibility();
  };

  for (;;) {
    if (limitReached(deadline_)) return PhaseResult::Limit;
    if (lu_.numUpdates() >= opt_.refactorInterval) {
      reinvert();
      if (tabooCount) {
        std::fill(taboo.begin(), taboo.end(), 0);
        tabooCount = 0;
      }
    }

    // ---- CHUZR: dual steepest edge ----
    int p = -1;
    double best = 0;
    for (int i = 0; i < m_; ++i) {
      if (taboo[i]) continue;
      int v = basicIndex_[i];
      double x = value_[v], inf;
      if (x < lower_[v] - tolP_[v]) inf = lower_[v] - x;
      else if (x > upper_[v] + tolP_[v]) inf = x - upper_[v];
      else continue;
      double score = inf * inf / dseWeight_[i];
      if (score > best) {
        best = score;
        p = i;
      }
    }
    if (p < 0) {
      if (tabooCount > 0) {
        // Remaining infeasible rows had only unstable pivots / unprovable rays.
        std::fill(taboo.begin(), taboo.end(), 0);
        tabooCount = 0;
        if (++allTabooEvents > 3) {
          PLOG_DETAIL("simplex: all infeasible rows tabooed (rejected pivots %lld, rejected rays %lld, pinf %.3e)",
                      stats_.rejectedPivots, stats_.rejectedRays, primalInfeasibilitySum());
          return PhaseResult::Numerical;
        }
        reinvert();
        if (primalInfeasibilitySum() > 0) continue;
      }
      return PhaseResult::Optimal;
    }
    const int leave = basicIndex_[p];
    const double xl = value_[leave];
    const int moveOut = xl < lower_[leave] ? -1 : 1;
    double delta = moveOut < 0 ? xl - lower_[leave] : xl - upper_[leave];

    // ---- BTRAN ----
    rho_.clear();
    rho_.array[p] = 1.0;
    rho_.index[0] = p;
    rho_.count = 1;
    lu_.btran(rho_);
    const double wp = std::max(rho_.norm2sq(), 1e-12);
    dseWeight_[p] = wp;

    // ---- PRICE ----
    priceRow(rho_, rowAp_, rowIdx_);

    // ---- CHUZC: bound-flipping ratio test with Harris tolerance ----
    cands.clear();
    for (int j : rowIdx_) {
      if (basicPos_[j] >= 0) continue;
      double lo = lower_[j], up = upper_[j];
      if (lo == up) continue;
      double a = rowAp_[j];
      double absA = std::fabs(a);
      if (absA < opt_.pivotTol) continue;
      int mv = move_[j];
      if (mv == 0) {  // free nonbasic: blocks in either direction
        double dj = std::fabs(dual_[j]);
        cands.push_back({j, dj / absA, (dj + tolD) / absA, absA});
      } else if (mv * a * moveOut > 0) {
        double dj = dual_[j] * mv;
        cands.push_back({j, std::max(dj, 0.0) / absA, std::max(dj + tolD, 0.0) / absA, absA});
      }
    }
    double slope = std::fabs(delta);
    flips.clear();
    group.clear();
    int enter = -1;
    size_t remaining = cands.size();
    while (remaining > 0) {
      double thMax = kInf;
      for (size_t k = 0; k < remaining; ++k) thMax = std::min(thMax, cands[k].harris);
      int bestK = -1;
      double bestAlpha = 0, total = 0;
      bool unboxed = false;
      for (size_t k = 0; k < remaining; ++k) {
        const Cand& c = cands[k];
        if (c.ratio > thMax) continue;
        if (c.absAlpha > bestAlpha) {
          bestAlpha = c.absAlpha;
          bestK = static_cast<int>(k);
        }
        double range = upper_[c.j] - lower_[c.j];
        if (!finite(lower_[c.j]) || !finite(upper_[c.j])) unboxed = true;
        else total += c.absAlpha * range;
      }
      // Pass a group only if the primal infeasibility left afterwards is significant;
      // otherwise round-off can "prove" infeasibility of a barely feasible row.
      if (!opt_.useBfrt || unboxed || slope - total <= tolP_[leave]) {
        enter = cands[bestK].j;
        for (size_t k = 0; k < remaining; ++k)
          if (cands[k].ratio <= thMax && static_cast<int>(k) != bestK) group.push_back(cands[k].j);
        break;
      }
      // Pass the whole group: all of it flips to the opposite bound.
      slope -= total;
      size_t w = 0;
      for (size_t k = 0; k < remaining; ++k) {
        if (cands[k].ratio <= thMax) flips.push_back(cands[k].j);
        else cands[w++] = cands[k];
      }
      remaining = w;
    }

    if (enter < 0) {
      if (phase1) return PhaseResult::Optimal;  // cannot happen for the boxed subproblem
      // Dual unbounded ray => primal infeasible, but only if the ray is a *rigorous*
      // Farkas certificate (checked on the scaled data; validity is scale invariant).
      std::vector<double> ys(m_);
      for (int i = 0; i < m_; ++i) ys[i] = rho_.array[i] * moveOut;
      if (farkasHolds(ys)) {
        farkas_.assign(m_, 0.0);
        for (int i = 0; i < m_; ++i) farkas_[i] = ys[i] * rowScale_[i];
        return PhaseResult::Infeasible;
      }
      stats_.rejectedRays++;
      if (lu_.numUpdates() > 0) {
        reinvert();
      } else {
        taboo[p] = 1;
        ++tabooCount;
      }
      continue;
    }

    // ---- FTRAN entering column + stability check ----
    columnOf(enter, col_);
    lu_.ftran(col_);
    const double alphaCol = col_.array[p];
    const double alphaRow = rowAp_[enter];
    const double err = std::fabs(alphaCol - alphaRow);
    if (std::fabs(alphaCol) < opt_.pivotTol || err > 1e-6 * std::max(1.0, std::fabs(alphaCol))) {
      stats_.rejectedPivots++;
      if (lu_.numUpdates() > 0) {
        stats_.numericalRefactors++;
        reinvert();
      } else {
        taboo[p] = 1;
        ++tabooCount;
      }
      continue;
    }

    // ---- Dual step ----
    double thetaD = dual_[enter] / alphaRow;
    if (thetaD * moveOut < 0) {  // entering reduced cost slightly infeasible: shift it to zero
      double shift = -dual_[enter];
      costWork_[enter] += shift;
      costsModified_ = true;
      stats_.costShifts++;
      stats_.maxCostShift = std::max(stats_.maxCostShift, std::fabs(shift));
      dual_[enter] = 0.0;
      thetaD = 0.0;
    }
    // Group members passed by Harris that would become dual infeasible: flip or shift.
    for (int j : group) {
      double nd = (dual_[j] - thetaD * rowAp_[j]);
      int mv = move_[j];
      double inf = mv == 1 ? -nd : (mv == -1 ? nd : std::fabs(nd));
      if (inf <= tolD_[j]) continue;
      if (finite(lower_[j]) && finite(upper_[j]) && mv != 0) {
        flips.push_back(j);
      }
    }
    for (int j : rowIdx_) {
      if (basicPos_[j] >= 0) continue;
      dual_[j] -= thetaD * rowAp_[j];
    }
    dual_[enter] = 0.0;
    dual_[leave] = -thetaD;
    for (int j : group) {  // residual infeasibilities on unboxed group members: shift costs
      if (std::find(flips.begin(), flips.end(), j) != flips.end()) continue;
      int mv = move_[j];
      double d = dual_[j];
      double inf = mv == 1 ? -d : (mv == -1 ? d : std::fabs(d));
      if (inf > tolD_[j]) {
        // Shift onto the feasible side with a random margin (not to exactly 0,
        // which would recreate the dual degeneracy perturbation removed).
        double target = mv == 0 ? 0.0 : mv * tolD_[j] * (1.0 + rng_.uniform());
        costWork_[j] += target - d;
        dual_[j] = target;
        costsModified_ = true;
        stats_.costShifts++;
        stats_.maxCostShift = std::max(stats_.maxCostShift, std::fabs(d));
      }
    }

    // ---- Bound flips (primal update) ----
    if (!flips.empty()) {
      flipCol_.clear();
      double* fc = flipCol_.array.data();
      for (int j : flips) {
        double from = value_[j];
        double to = move_[j] == 1 ? upper_[j] : lower_[j];
        double dx = to - from;
        value_[j] = to;
        move_[j] = static_cast<int8_t>(-move_[j]);
        if (j >= n_) {
          fc[j - n_] += dx;
        } else {
          for (int k = A_.start[j]; k < A_.start[j + 1]; ++k) fc[A_.index[k]] += A_.value[k] * dx;
        }
      }
      lu_.ftran(flipCol_);
      for (int k = 0; k < flipCol_.count; ++k) {
        int pos = flipCol_.index[k];
        value_[basicIndex_[pos]] -= flipCol_.array[pos];
      }
      flipCol_.clear();
      stats_.boundFlips += static_cast<long long>(flips.size());
    }

    // ---- DSE: tau = B^{-1} rho ----
    tau_.clear();
    for (int k = 0; k < rho_.count; ++k) {
      int i = rho_.index[k];
      tau_.array[i] = rho_.array[i];
    }
    lu_.ftran(tau_);

    // ---- Primal step ----
    const double bound = moveOut < 0 ? lower_[leave] : upper_[leave];
    delta = value_[leave] - bound;
    const double thetaP = delta / alphaCol;
    if (std::fabs(thetaP) < 1e-12) stats_.degeneratePivots++;
    for (int k = 0; k < col_.count; ++k) {
      int pos = col_.index[k];
      value_[basicIndex_[pos]] -= thetaP * col_.array[pos];
    }
    value_[enter] += thetaP;
    value_[leave] = bound;

    // ---- DSE weight update ----
    for (int k = 0; k < col_.count; ++k) {
      int i = col_.index[k];
      if (i == p) continue;
      double r = col_.array[i] / alphaCol;
      double w = dseWeight_[i] + r * (r * wp - 2.0 * tau_.array[i]);
      dseWeight_[i] = std::max(w, 1e-4);
    }
    dseWeight_[p] = std::max(wp / (alphaCol * alphaCol), 1e-4);
    tau_.clear();

    // ---- Basis change ----
    updateBasis(p, enter, col_);
    if (lower_[leave] == upper_[leave]) move_[leave] = 0;
    else move_[leave] = static_cast<int8_t>(moveOut < 0 ? 1 : -1);
    col_.clear();
    rho_.clear();
    stats_.iterations++;
    if (phase1) stats_.dualPhase1Iterations++;
    else stats_.dualPhase2Iterations++;
    if (tabooCount) {
      std::fill(taboo.begin(), taboo.end(), 0);
      tabooCount = 0;
    }
    // ---- Stall detection: no dual objective progress => stronger perturbation ----
    {
      double obj = 0;
      if ((stats_.iterations & 63) == 0) {
        obj = computeObjectiveWork();
        if (obj > stallBest + 1e-12 * (1.0 + std::fabs(stallBest))) {
          stallBest = obj;
          stallIter = stats_.iterations;
        } else if (stats_.iterations - stallIter > std::max<long long>(1000, m_ / 2)) {
          stallFactor = std::min(stallFactor * 10.0, 1e4);
          perturbCosts(stallFactor);
          stats_.stallRecoveries++;
          stallIter = stats_.iterations;
          stallBest = -kInf;
          PLOG_DETAIL("simplex: stall detected at it %lld, re-perturbing (x%g)", stats_.iterations, stallFactor);
        }
      }
    }
    logIteration(phase1 ? "dual-1" : "dual");
  }
}

// ---------------------------------------------------------------------------
// Dual phase 1: artificial bounding subproblem, solved by dual phase 2.
// ---------------------------------------------------------------------------
void Simplex::dualPhase1() {
  const int N = n_ + m_;
  std::vector<double> saveLo = lower_, saveUp = upper_;
  for (int j = 0; j < N; ++j) {
    double lo = lower_[j], up = upper_[j];
    if (lo == up || (finite(lo) && finite(up))) {
      lower_[j] = upper_[j] = 0.0;
    } else if (finite(lo)) {
      lower_[j] = 0.0;
      upper_[j] = 1.0;
    } else if (finite(up)) {
      lower_[j] = -1.0;
      upper_[j] = 0.0;
    } else {
      lower_[j] = -1000.0;
      upper_[j] = 1000.0;
    }
  }
  for (int j = 0; j < N; ++j)
    if (basicPos_[j] < 0) setNonbasicValue(j, dual_[j] < 0 ? -1 : 1);
  // The subproblem is highly dual degenerate: perturb costs (removed afterwards)
  // and cap the effort; residual dual infeasibilities are handled by cost shifting.
  if (opt_.perturb) perturbCosts();
  computePrimal();
  const long long saveMax = opt_.maxIterations;
  opt_.maxIterations = std::min(saveMax, stats_.iterations + 20LL * (n_ + m_) + 1000);
  PhaseResult r = dualPhase2(true);
  opt_.maxIterations = saveMax;
  if (r == PhaseResult::Limit) PLOG_DETAIL("simplex: dual phase 1 capped; continuing with cost shifting");
  if (costsModified_) restoreCosts();
  lower_ = saveLo;
  upper_ = saveUp;
  for (int j = 0; j < N; ++j)
    if (basicPos_[j] < 0) setNonbasicValue(j, dual_[j] < 0 ? -1 : 1);
  computePrimal();
}

// ---------------------------------------------------------------------------
// Primal simplex phase 2 (Devex pricing, Harris ratio test with bound flips)
// ---------------------------------------------------------------------------
Simplex::PhaseResult Simplex::primalPhase2() {
  const double tolP = opt_.primalFeasTol;
  const double tolD = opt_.dualFeasTol;
  std::vector<char> taboo(n_ + m_, 0);
  int tabooCount = 0;
  auto reinvert = [&]() {
    factor();
    computePrimal();
    computeDual();
  };
  for (;;) {
    if (limitReached(deadline_)) return PhaseResult::Limit;
    if (lu_.numUpdates() >= opt_.refactorInterval) reinvert();

    // ---- CHUZC: Devex ----
    int q = -1;
    double best = 0;
    for (int j = 0; j < n_ + m_; ++j) {
      if (basicPos_[j] >= 0 || lower_[j] == upper_[j] || taboo[j]) continue;
      double d = dual_[j];
      double inf;
      if (move_[j] == 1) inf = -d;
      else if (move_[j] == -1) inf = d;
      else inf = std::fabs(d);
      if (inf <= tolD_[j]) continue;
      double score = inf * inf / devex_[j];
      if (score > best) {
        best = score;
        q = j;
      }
    }
    if (q < 0) {
      if (tabooCount > 0) {
        std::fill(taboo.begin(), taboo.end(), 0);
        tabooCount = 0;
        reinvert();
        if (countDualInfeasibilities() > 0 && stats_.rejectedPivots > 50 * (m_ + 10)) return PhaseResult::Numerical;
        continue;
      }
      return PhaseResult::Optimal;
    }
    const int dir = dual_[q] < 0 ? 1 : -1;

    // ---- FTRAN ----
    columnOf(q, col_);
    lu_.ftran(col_);

    // ---- Harris ratio test ----
    double thMax = kInf;
    for (int k = 0; k < col_.count; ++k) {
      int i = col_.index[k];
      double a = col_.array[i] * dir;
      if (std::fabs(a) < opt_.pivotTol) continue;
      int v = basicIndex_[i];
      if (a > 0) {
        if (finite(lower_[v])) thMax = std::min(thMax, (value_[v] - lower_[v] + tolP) / a);
      } else {
        if (finite(upper_[v])) thMax = std::min(thMax, (upper_[v] - value_[v] + tolP) / -a);
      }
    }
    const double range = (finite(lower_[q]) && finite(upper_[q])) ? upper_[q] - lower_[q] : kInf;
    if (thMax == kInf && range == kInf) {
      // Unbounded: primal ray (structural components, unscaled).
      primalRay_.assign(n_, 0.0);
      if (q < n_) primalRay_[q] = dir * colScale_[q];
      for (int k = 0; k < col_.count; ++k) {
        int i = col_.index[k];
        int v = basicIndex_[i];
        if (v < n_) primalRay_[v] = -dir * col_.array[i] * colScale_[v];
      }
      col_.clear();
      return PhaseResult::Unbounded;
    }
    int p = -1;
    double theta = kInf, bestA = 0;
    for (int k = 0; k < col_.count; ++k) {
      int i = col_.index[k];
      double a = col_.array[i] * dir;
      if (std::fabs(a) < opt_.pivotTol) continue;
      int v = basicIndex_[i];
      double r;
      if (a > 0) {
        if (!finite(lower_[v])) continue;
        r = std::max(value_[v] - lower_[v], 0.0) / a;
      } else {
        if (!finite(upper_[v])) continue;
        r = std::max(upper_[v] - value_[v], 0.0) / -a;
      }
      if (r <= thMax && std::fabs(a) > bestA) {
        bestA = std::fabs(a);
        p = i;
        theta = r;
      }
    }
    if (range <= theta) {
      // Bound flip of the entering variable, no basis change.
      for (int k = 0; k < col_.count; ++k) {
        int i = col_.index[k];
        value_[basicIndex_[i]] -= dir * range * col_.array[i];
      }
      value_[q] = move_[q] == 1 ? upper_[q] : lower_[q];
      move_[q] = static_cast<int8_t>(-move_[q]);
      col_.clear();
      stats_.boundFlips++;
      stats_.iterations++;
      stats_.primalIterations++;
      // Duals are unchanged by a bound flip.
      continue;
    }
    if (p < 0) {  // only tiny pivots available
      taboo[q] = 1;
      ++tabooCount;
      stats_.rejectedPivots++;
      col_.clear();
      continue;
    }

    // ---- BTRAN + PRICE for dual update ----
    rho_.clear();
    rho_.array[p] = 1.0;
    rho_.index[0] = p;
    rho_.count = 1;
    lu_.btran(rho_);
    priceRow(rho_, rowAp_, rowIdx_);
    const double alphaCol = col_.array[p];
    const double alphaRow = q >= n_ ? rho_.array[q - n_] : rowAp_[q];
    if (std::fabs(alphaCol - alphaRow) > 1e-6 * std::max(1.0, std::fabs(alphaCol))) {
      stats_.rejectedPivots++;
      if (lu_.numUpdates() > 0) {
        stats_.numericalRefactors++;
        reinvert();
      } else {
        taboo[q] = 1;
        ++tabooCount;
      }
      col_.clear();
      rho_.clear();
      continue;
    }

    // ---- Primal update ----
    const int leave = basicIndex_[p];
    if (theta < 1e-12) stats_.degeneratePivots++;
    for (int k = 0; k < col_.count; ++k) {
      int i = col_.index[k];
      value_[basicIndex_[i]] -= dir * theta * col_.array[i];
    }
    value_[q] += dir * theta;
    const double aLeave = col_.array[p] * dir;
    const int leaveMove = aLeave > 0 ? 1 : -1;  // decreasing -> leaves at lower
    value_[leave] = leaveMove == 1 ? lower_[leave] : upper_[leave];

    // ---- Dual update ----
    const double thetaD = dual_[q] / alphaCol;
    for (int j : rowIdx_) {
      if (basicPos_[j] >= 0) continue;
      dual_[j] -= thetaD * rowAp_[j];
    }
    // Devex reference weights.
    const double gq = devex_[q];
    for (int j : rowIdx_) {
      if (basicPos_[j] >= 0) continue;
      double r = rowAp_[j] / alphaCol;
      devex_[j] = std::max(devex_[j], r * r * gq);
    }
    devex_[leave] = std::max(gq / (alphaCol * alphaCol), 1.0);

    updateBasis(p, q, col_);
    dual_[leave] = -thetaD;
    move_[leave] = lower_[leave] == upper_[leave] ? 0 : static_cast<int8_t>(leaveMove);
    col_.clear();
    rho_.clear();
    stats_.iterations++;
    stats_.primalIterations++;
    if (tabooCount) {
      std::fill(taboo.begin(), taboo.end(), 0);
      tabooCount = 0;
    }
    logIteration("primal");
  }
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------
Status Simplex::solve(const Deadline* deadline) {
  deadline_ = deadline;
  timer_.reset();
  farkas_.clear();
  primalRay_.clear();
  if (!haveBasis_) initSlackBasis();
  if (!factored_ && !factor()) {
    initSlackBasis();
    if (!factor()) return Status::NumericalFailure;
  }
  costWork_ = cost_;
  costsModified_ = false;
  computePrimal();
  computeDual();
  bool triedPhase1 = false;
  auto limitStatus = [&]() {
    stats_.timeSeconds += timer_.seconds();
    return stats_.iterations >= opt_.maxIterations ? Status::IterationLimit : Status::TimeLimit;
  };
  // Dual polish: re-run the primal simplex with a 100x tighter dual tolerance so
  // the final duals admit a rigorous (not merely tolerance-based) certificate.
  auto finishOptimal = [&]() {
    std::vector<double> saveTol = tolD_;
    for (auto& t : tolD_) t *= 0.01;
    if (countDualInfeasibilities() > 0) {
      const long long saveMax = opt_.maxIterations;
      opt_.maxIterations = std::min(saveMax, stats_.iterations + n_ + m_ + 100);
      std::vector<double> saveVal = value_;
      std::vector<int> saveBasic = basicIndex_;
      PhaseResult pr = primalPhase2();
      opt_.maxIterations = saveMax;
      if (pr == PhaseResult::Optimal) {
        stats_.polishIterations++;
      } else {
        // Polish failed: keep the (tolerance-optimal) basis we had.
        basicIndex_ = saveBasic;
        std::fill(basicPos_.begin(), basicPos_.end(), -1);
        for (int p = 0; p < m_; ++p) basicPos_[basicIndex_[p]] = p;
        value_ = saveVal;
        for (int p = 0; p < m_; ++p) move_[basicIndex_[p]] = 0;
        for (int j = 0; j < n_ + m_; ++j)
          if (basicPos_[j] < 0 && move_[j] == 0 && lower_[j] != upper_[j] && (finite(lower_[j]) || finite(upper_[j])))
            setNonbasicValue(j, value_[j] == upper_[j] ? -1 : 1);
        factor();
        computePrimal();
      }
      computeDual();
    }
    tolD_ = saveTol;
    stats_.timeSeconds += timer_.seconds();
    return Status::Optimal;
  };
  for (int round = 0; round < 8; ++round) {
    flipToDualFeasibility(false);
    int nInf = countDualInfeasibilities();
    if (nInf > 0 && !triedPhase1) {
      triedPhase1 = true;
      dualPhase1();
      if (limitReached(deadline_)) return limitStatus();
      flipToDualFeasibility(false);
      nInf = countDualInfeasibilities();
    }
    if (nInf > 0) shiftCostsForDualFeasibility();
    if (opt_.perturb && round == 0) perturbCosts();
    computePrimal();
    PhaseResult r = dualPhase2(false);
    if (r == PhaseResult::Limit) return limitStatus();
    if (r == PhaseResult::Infeasible) {
      if (costsModified_) restoreCosts();
      stats_.timeSeconds += timer_.seconds();
      return Status::Infeasible;
    }
    if (r == PhaseResult::Numerical) {
      PLOG_DETAIL("simplex: numerical trouble in dual phase 2, refactoring (round %d)", round);
      factor();
      computePrimal();
      costWork_ = cost_;
      costsModified_ = false;
      computeDual();
      continue;
    }
    if (costsModified_) restoreCosts();
    else computeDual();
    if (countDualInfeasibilities() == 0) return finishOptimal();
    r = primalPhase2();
    if (r == PhaseResult::Limit) return limitStatus();
    if (r == PhaseResult::Unbounded) {
      if (primalInfeasibilitySum() == 0) {
        stats_.timeSeconds += timer_.seconds();
        return Status::Unbounded;
      }
      // Ray found from a primal-infeasible point: unbounded or infeasible; resolve feasibility.
      costWork_.assign(n_ + m_, 0.0);
      computeDual();
      PhaseResult f = dualPhase2(false);
      costWork_ = cost_;
      computeDual();
      if (f == PhaseResult::Infeasible) {
        stats_.timeSeconds += timer_.seconds();
        return Status::Infeasible;
      }
      continue;
    }
    factor();
    computePrimal();
    computeDual();
    if (primalInfeasibilitySum() == 0 && countDualInfeasibilities() == 0) return finishOptimal();
  }
  stats_.timeSeconds += timer_.seconds();
  return Status::NumericalFailure;
}

// ---------------------------------------------------------------------------
// Modification
// ---------------------------------------------------------------------------
void Simplex::setColBounds(int j, double lo, double up) {
  lower_[j] = finite(lo) ? lo / colScale_[j] : -kInf;
  upper_[j] = finite(up) ? up / colScale_[j] : kInf;
  if (haveBasis_ && basicPos_[j] < 0) setNonbasicValue(j, move_[j] == -1 ? -1 : 1);
}

void Simplex::setRowBounds(int i, double lo, double up) {
  int v = n_ + i;
  lower_[v] = finite(up) ? -up * rowScale_[i] : -kInf;
  upper_[v] = finite(lo) ? -lo * rowScale_[i] : kInf;
  if (haveBasis_ && basicPos_[v] < 0) setNonbasicValue(v, move_[v] == -1 ? -1 : 1);
}

void Simplex::setColCost(int j, double c) {
  cost_[j] = c * colScale_[j];
  costWork_[j] = cost_[j];
}

void Simplex::addRows(const std::vector<double>& lower, const std::vector<double>& upper,
                      const std::vector<int>& starts, const std::vector<int>& cols,
                      const std::vector<double>& vals) {
  const int k = static_cast<int>(lower.size());
  if (k == 0) return;
  const int oldM = m_;
  // Scale factors for the new rows.
  std::vector<double> newScale(k, 1.0);
  for (int r = 0; r < k; ++r) {
    double mx = 0;
    for (int e = starts[r]; e < starts[r + 1]; ++e) mx = std::max(mx, std::fabs(vals[e] * colScale_[cols[e]]));
    if (opt_.scale && mx > 0) newScale[r] = std::ldexp(1.0, -static_cast<int>(std::lround(std::log2(mx))));
  }
  // Rebuild scaled A with appended rows.
  std::vector<int> ri, ci;
  std::vector<double> v;
  ri.reserve(A_.nnz() + cols.size());
  for (int j = 0; j < n_; ++j)
    for (int e = A_.start[j]; e < A_.start[j + 1]; ++e) {
      ri.push_back(A_.index[e]);
      ci.push_back(j);
      v.push_back(A_.value[e]);
    }
  for (int r = 0; r < k; ++r)
    for (int e = starts[r]; e < starts[r + 1]; ++e) {
      ri.push_back(oldM + r);
      ci.push_back(cols[e]);
      v.push_back(vals[e] * colScale_[cols[e]] * newScale[r]);
    }
  m_ = oldM + k;
  A_ = SparseMatrix::fromTriplets(m_, n_, ri, ci, v);
  AT_ = A_.transpose();
  for (int r = 0; r < k; ++r) {
    rowScale_.push_back(newScale[r]);
    double lo = finite(upper[r]) ? -upper[r] * newScale[r] : -kInf;
    double up = finite(lower[r]) ? -lower[r] * newScale[r] : kInf;
    cost_.push_back(0.0);
    costWork_.push_back(0.0);
    lower_.push_back(lo);
    upper_.push_back(up);
    value_.push_back(0.0);
    dual_.push_back(0.0);
    move_.push_back(0);
    basicPos_.push_back(oldM + r);
    basicIndex_.push_back(n_ + oldM + r);
    dseWeight_.push_back(1.0);
    devex_.push_back(1.0);
    rowAp_.push_back(0.0);
    rowMark_.push_back(0);
    y_.push_back(0.0);
  }
  lu_.setup(m_, n_, &A_);
  col_.setup(m_);
  rho_.setup(m_);
  tau_.setup(m_);
  flipCol_.setup(m_);
  computeTolerances();
  factored_ = false;
}

int Simplex::removeBasicRows(const std::vector<int>& rows) {
  std::vector<char> drop(m_, 0);
  int cnt = 0;
  for (int i : rows)
    if (i >= 0 && i < m_ && basicPos_[n_ + i] >= 0 && !drop[i]) {
      drop[i] = 1;
      ++cnt;
    }
  if (cnt == 0) return 0;
  std::vector<int> newRow(m_, -1);
  int nm = 0;
  for (int i = 0; i < m_; ++i)
    if (!drop[i]) newRow[i] = nm++;
  std::vector<int> ri, ci;
  std::vector<double> v;
  for (int j = 0; j < n_; ++j)
    for (int e = A_.start[j]; e < A_.start[j + 1]; ++e) {
      int r = newRow[A_.index[e]];
      if (r < 0) continue;
      ri.push_back(r);
      ci.push_back(j);
      v.push_back(A_.value[e]);
    }
  A_ = SparseMatrix::fromTriplets(nm, n_, ri, ci, v);
  AT_ = A_.transpose();
  auto compactN = [&](auto& vec) {
    for (int i = 0; i < m_; ++i)
      if (newRow[i] >= 0) vec[n_ + newRow[i]] = vec[n_ + i];
    vec.resize(n_ + nm);
  };
  auto compactM = [&](auto& vec) {
    for (int i = 0; i < m_; ++i)
      if (newRow[i] >= 0) vec[newRow[i]] = vec[i];
    vec.resize(nm);
  };
  compactN(cost_);
  compactN(costWork_);
  compactN(lower_);
  compactN(upper_);
  compactN(value_);
  compactN(dual_);
  compactN(move_);
  compactN(devex_);
  compactM(rowScale_);
  compactM(y_);
  // Basis: drop positions holding the removed logicals, remap logical indices.
  std::vector<int> nb;
  std::vector<double> nw;
  for (int p = 0; p < m_; ++p) {
    int v2 = basicIndex_[p];
    if (v2 >= n_) {
      int r = v2 - n_;
      if (drop[r]) continue;
      v2 = n_ + newRow[r];
    }
    nb.push_back(v2);
    nw.push_back(dseWeight_[p]);
  }
  m_ = nm;
  basicIndex_ = nb;
  dseWeight_ = nw;
  basicPos_.assign(n_ + m_, -1);
  for (int p = 0; p < m_; ++p) basicPos_[basicIndex_[p]] = p;
  rowAp_.assign(n_ + m_, 0.0);
  rowMark_.assign(n_ + m_, 0);
  rowIdx_.clear();
  lu_.setup(m_, n_, &A_);
  col_.setup(m_);
  rho_.setup(m_);
  tau_.setup(m_);
  flipCol_.setup(m_);
  computeTolerances();
  factored_ = false;
  return cnt;
}

// ---------------------------------------------------------------------------
// Basis I/O
// ---------------------------------------------------------------------------
void Simplex::getBasis(std::vector<BasisStatus>& colStatus, std::vector<BasisStatus>& rowStatus) const {
  colStatus.assign(n_, BasisStatus::Lower);
  rowStatus.assign(m_, BasisStatus::Basic);
  for (int j = 0; j < n_; ++j) {
    if (basicPos_[j] >= 0) colStatus[j] = BasisStatus::Basic;
    else if (lower_[j] == upper_[j]) colStatus[j] = BasisStatus::Fixed;
    else if (move_[j] == 1) colStatus[j] = BasisStatus::Lower;
    else if (move_[j] == -1) colStatus[j] = BasisStatus::Upper;
    else colStatus[j] = BasisStatus::Zero;
  }
  for (int i = 0; i < m_; ++i) {
    int v = n_ + i;
    if (basicPos_[v] >= 0) rowStatus[i] = BasisStatus::Basic;
    else if (lower_[v] == upper_[v]) rowStatus[i] = BasisStatus::Fixed;
    else if (move_[v] == 1) rowStatus[i] = BasisStatus::Upper;  // logical at its lower = row at upper
    else if (move_[v] == -1) rowStatus[i] = BasisStatus::Lower;
    else rowStatus[i] = BasisStatus::Zero;
  }
}

void Simplex::setBasis(const std::vector<BasisStatus>& colStatus, const std::vector<BasisStatus>& rowStatus) {
  int basics = 0;
  for (auto s : colStatus) basics += s == BasisStatus::Basic;
  for (auto s : rowStatus) basics += s == BasisStatus::Basic;
  if (static_cast<int>(colStatus.size()) != n_ || static_cast<int>(rowStatus.size()) != m_ || basics != m_) {
    initSlackBasis();
    return;
  }
  std::fill(basicPos_.begin(), basicPos_.end(), -1);
  int p = 0;
  for (int j = 0; j < n_; ++j) {
    if (colStatus[j] == BasisStatus::Basic) {
      basicIndex_[p] = j;
      basicPos_[j] = p++;
      move_[j] = 0;
    } else {
      setNonbasicValue(j, colStatus[j] == BasisStatus::Upper ? -1 : 1);
    }
  }
  for (int i = 0; i < m_; ++i) {
    int v = n_ + i;
    if (rowStatus[i] == BasisStatus::Basic) {
      basicIndex_[p] = v;
      basicPos_[v] = p++;
      move_[v] = 0;
    } else {
      setNonbasicValue(v, rowStatus[i] == BasisStatus::Lower ? -1 : 1);
    }
  }
  std::fill(dseWeight_.begin(), dseWeight_.end(), 1.0);
  std::fill(devex_.begin(), devex_.end(), 1.0);
  haveBasis_ = true;
  factored_ = false;
}

// ---------------------------------------------------------------------------
// Results
// ---------------------------------------------------------------------------
double Simplex::variableValue(int v) const {
  if (v < n_) return value_[v] * colScale_[v];
  return value_[v] / rowScale_[v - n_];
}

double Simplex::objective() const {
  double f = offset_;
  for (int j = 0; j < n_; ++j) f += cost_[j] * value_[j];  // scaled c'x' == c x
  return f;
}

std::vector<double> Simplex::colValues() const {
  std::vector<double> x(n_);
  for (int j = 0; j < n_; ++j) x[j] = value_[j] * colScale_[j];
  return x;
}

std::vector<double> Simplex::rowActivities() const {
  std::vector<double> r(m_);
  for (int i = 0; i < m_; ++i) r[i] = -value_[n_ + i] / rowScale_[i];
  return r;
}

std::vector<double> Simplex::rowDuals() const {
  std::vector<double> y(m_);
  for (int i = 0; i < m_; ++i) y[i] = y_[i] * rowScale_[i];
  return y;
}

std::vector<double> Simplex::reducedCosts() const {
  std::vector<double> d(n_);
  for (int j = 0; j < n_; ++j) d[j] = basicPos_[j] >= 0 ? 0.0 : dual_[j] / colScale_[j];
  return d;
}

void Simplex::basicDirection(int v, std::vector<double>& out) const {
  HVector col;
  col.setup(m_);
  if (v >= n_) {
    col.array[v - n_] = 1.0;
  } else {
    for (int k = A_.start[v]; k < A_.start[v + 1]; ++k) col.array[A_.index[k]] = A_.value[k];
  }
  lu_.ftran(col);
  out.assign(m_, 0.0);
  const double sv = variableScale(v);
  for (int p = 0; p < m_; ++p) out[p] = -col.array[p] * variableScale(basicIndex_[p]) / sv;
}

int Simplex::tableauRow(int pos, std::vector<double>& coefs) const {
  HVector rho;
  rho.setup(m_);
  rho.array[pos] = 1.0;
  rho.index[0] = pos;
  rho.count = 1;
  lu_.btran(rho);
  const int bv = basicIndex_[pos];
  auto vs = [&](int v) { return v < n_ ? colScale_[v] : 1.0 / rowScale_[v - n_]; };
  const double sB = vs(bv);
  coefs.assign(n_ + m_, 0.0);
  for (int j = 0; j < n_; ++j) {
    if (basicPos_[j] >= 0) continue;
    double s = 0;
    for (int k = A_.start[j]; k < A_.start[j + 1]; ++k) s += A_.value[k] * rho.array[A_.index[k]];
    coefs[j] = s * sB / vs(j);
  }
  for (int i = 0; i < m_; ++i) {
    int v = n_ + i;
    if (basicPos_[v] >= 0) continue;
    coefs[v] = rho.array[i] * sB / vs(v);
  }
  return bv;
}

LpResult solveLpSimplex(const Model& lp, const SimplexOptions& opts, const Deadline* deadline,
                        const std::vector<BasisStatus>* warmCol, const std::vector<BasisStatus>* warmRow) {
  Simplex s(lp, opts);
  if (warmCol && warmRow) s.setBasis(*warmCol, *warmRow);
  LpResult res;
  res.status = s.solve(deadline);
  res.stats = s.stats();
  res.x = s.colValues();
  res.rowActivity = s.rowActivities();
  res.rowDual = s.rowDuals();
  res.reducedCost = s.reducedCosts();
  res.objective = s.objective();
  s.getBasis(res.colStatus, res.rowStatus);
  res.farkas = s.farkasRay();
  res.primalRay = s.primalRay();
  return res;
}

}  // namespace pramana
