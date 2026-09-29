#include "presolve/presolve.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "util/log.h"
#include "util/timer.h"

namespace pramana {

namespace {
inline bool fin(double v) { return std::fabs(v) < kInfBoundThreshold; }
constexpr double kFeasTol = 1e-9;
}  // namespace

double Presolver::activityMin(int i) const {
  double s = 0;
  for (int e = AT_.start[i]; e < AT_.start[i + 1]; ++e) {
    int j = AT_.index[e];
    if (!colActive_[j]) continue;
    double a = AT_.value[e];
    double b = a > 0 ? lo_[j] : up_[j];
    if (!fin(b)) return -kInf;
    s += a * b;
  }
  return s;
}

double Presolver::activityMax(int i) const {
  double s = 0;
  for (int e = AT_.start[i]; e < AT_.start[i + 1]; ++e) {
    int j = AT_.index[e];
    if (!colActive_[j]) continue;
    double a = AT_.value[e];
    double b = a > 0 ? up_[j] : lo_[j];
    if (!fin(b)) return kInf;
    s += a * b;
  }
  return s;
}

void Presolver::removeRow(int i) {
  rowActive_[i] = 0;
  stats_.rowsRemoved++;
  for (int e = AT_.start[i]; e < AT_.start[i + 1]; ++e) {
    int j = AT_.index[e];
    if (colActive_[j]) {
      colCount_[j]--;
      stats_.nnzRemoved++;
    }
  }
}

void Presolver::fixColumn(int j, double v, bool record) {
  const SparseMatrix& A = orig_.A;
  for (int k = A.start[j]; k < A.start[j + 1]; ++k) {
    int i = A.index[k];
    if (!rowActive_[i]) continue;
    double a = A.value[k];
    if (fin(rl_[i])) rl_[i] -= a * v;
    if (fin(ru_[i])) ru_[i] -= a * v;
    rowCount_[i]--;
    stats_.nnzRemoved++;
  }
  offset_ += cost_[j] * v;
  colActive_[j] = 0;
  stats_.colsRemoved++;
  if (record) {
    Record r;
    r.kind = Kind::FixedCol;
    r.col = j;
    r.value = v;
    r.cost = cost_[j];
    stack_.push_back(std::move(r));
  }
}

bool Presolver::pass(bool mip) {
  const int m = orig_.numRows(), n = orig_.numCols();
  const SparseMatrix& A = orig_.A;
  bool changed = false;

  // ---- Rows ----
  for (int i = 0; i < m && !infeasible_; ++i) {
    if (!rowActive_[i]) continue;
    if (rowCount_[i] == 0 || (!fin(rl_[i]) && !fin(ru_[i]))) {
      if (rowCount_[i] == 0 && ((fin(rl_[i]) && rl_[i] > 1e-7 * (1 + std::fabs(rl_[i]))) ||
                                (fin(ru_[i]) && ru_[i] < -1e-7 * (1 + std::fabs(ru_[i]))))) {
        infeasible_ = true;
        break;
      }
      Record r;
      r.kind = Kind::EmptyRow;
      r.row = i;
      stack_.push_back(r);
      removeRow(i);
      stats_.emptyRows++;
      changed = true;
      continue;
    }
    if (rowCount_[i] == 1) {
      int j = -1;
      double a = 0;
      for (int e = AT_.start[i]; e < AT_.start[i + 1]; ++e)
        if (colActive_[AT_.index[e]]) {
          j = AT_.index[e];
          a = AT_.value[e];
        }
      double ilo, iup;
      if (a > 0) {
        ilo = fin(rl_[i]) ? rl_[i] / a : -kInf;
        iup = fin(ru_[i]) ? ru_[i] / a : kInf;
      } else {
        ilo = fin(ru_[i]) ? ru_[i] / a : -kInf;
        iup = fin(rl_[i]) ? rl_[i] / a : kInf;
      }
      bool isInt = orig_.colType[j] == VarType::Integer;
      if (isInt) {
        if (fin(ilo)) ilo = std::ceil(ilo - 1e-9);
        if (fin(iup)) iup = std::floor(iup + 1e-9);
      }
      Record r;
      r.kind = Kind::SingletonRow;
      r.row = i;
      r.col = j;
      r.a = a;
      r.oldLower = lo_[j];
      r.oldUpper = up_[j];
      r.rowLower = rl_[i];
      r.rowUpper = ru_[i];
      r.cost = cost_[j];
      double nlo = std::max(lo_[j], ilo), nup = std::min(up_[j], iup);
      if (nlo > nup + kFeasTol * (1 + std::fabs(nlo))) {
        infeasible_ = true;
        break;
      }
      if (nlo > nup) nlo = nup;  // round-off
      if (ilo > lo_[j]) stats_.boundsTightened++;
      if (iup < up_[j]) stats_.boundsTightened++;
      lo_[j] = nlo;
      up_[j] = nup;
      stack_.push_back(r);
      removeRow(i);
      stats_.singletonRows++;
      changed = true;
      continue;
    }
    double amin = activityMin(i), amax = activityMax(i);
    double tol = 1e-9 * (1 + std::max(std::fabs(rl_[i]) * fin(rl_[i]), std::fabs(ru_[i]) * fin(ru_[i])));
    // Strict violation of a row's range by its activity bounds is a one-row Farkas
    // certificate candidate (checked rigorously later on the original model).
    if (certRow_ < 0 && ((fin(ru_[i]) && fin(amin) && amin > ru_[i]) || (fin(rl_[i]) && fin(amax) && amax < rl_[i]))) {
      certRow_ = i;
      certSign_ = (fin(rl_[i]) && fin(amax) && amax < rl_[i]) ? 1 : -1;
    }
    if ((fin(ru_[i]) && amin > ru_[i] + 1e3 * tol) || (fin(rl_[i]) && amax < rl_[i] - 1e3 * tol)) {
      infeasible_ = true;
      break;
    }
    if ((!fin(rl_[i]) || amin >= rl_[i] - tol) && (!fin(ru_[i]) || amax <= ru_[i] + tol)) {
      Record r;
      r.kind = Kind::RedundantRow;
      r.row = i;
      stack_.push_back(r);
      removeRow(i);
      stats_.redundantRows++;
      changed = true;
      continue;
    }
    int side = 0;
    if (fin(ru_[i]) && fin(amin) && amin >= ru_[i] - tol) side = 1;       // forced to min activity
    else if (fin(rl_[i]) && fin(amax) && amax <= rl_[i] + tol) side = -1;  // forced to max activity
    if (side != 0 && !qp_) {  // forcing rows fix columns: not valid with a quadratic objective
      Record r;
      r.kind = Kind::ForcingRow;
      r.row = i;
      r.side = side;
      r.rowLower = rl_[i];
      r.rowUpper = ru_[i];
      std::vector<std::pair<int, double>> fixes;
      for (int e = AT_.start[i]; e < AT_.start[i + 1]; ++e) {
        int j = AT_.index[e];
        if (!colActive_[j]) continue;
        double a = AT_.value[e];
        double v = side == 1 ? (a > 0 ? lo_[j] : up_[j]) : (a > 0 ? up_[j] : lo_[j]);
        r.cols.push_back(j);
        r.vals.push_back(a);
        r.colVals.push_back(v);
        fixes.emplace_back(j, v);
      }
      std::vector<double> stageCost;
      for (auto& f : fixes) stageCost.push_back(cost_[f.first]);
      removeRow(i);
      for (auto& f : fixes) fixColumn(f.first, f.second, false);
      // Stage costs of the fixed columns are kept in forcingCosts_[value].
      forcingCosts_.push_back(stageCost);
      r.value = static_cast<double>(forcingCosts_.size() - 1);
      stack_.push_back(std::move(r));
      stats_.forcingRows++;
      changed = true;
      continue;
    }
  }
  if (infeasible_) return false;

  // ---- Columns ----
  for (int j = 0; j < n && !infeasible_ && !unbounded_; ++j) {
    if (!colActive_[j]) continue;
    bool isInt = orig_.colType[j] == VarType::Integer;
    if (isInt) {
      double nl = fin(lo_[j]) ? std::ceil(lo_[j] - 1e-9) : lo_[j];
      double nu = fin(up_[j]) ? std::floor(up_[j] + 1e-9) : up_[j];
      if (nl != lo_[j] || nu != up_[j]) {
        lo_[j] = nl;
        up_[j] = nu;
        changed = true;
      }
    }
    if (lo_[j] > up_[j] + kFeasTol * (1 + std::fabs(lo_[j]))) {
      infeasible_ = true;
      break;
    }
    // QP: column reductions interact with Q (the IPM eliminates fixed columns itself).
    if (qp_) continue;
    if (fin(lo_[j]) && fin(up_[j]) && up_[j] - lo_[j] <= 1e-12 * (1 + std::fabs(lo_[j]))) {
      fixColumn(j, lo_[j]);
      stats_.fixedCols++;
      changed = true;
      continue;
    }
    const double c = cost_[j];
    if (colCount_[j] == 0) {
      double v;
      if (c > 0) v = lo_[j];
      else if (c < 0) v = up_[j];
      else v = fin(lo_[j]) ? lo_[j] : (fin(up_[j]) ? up_[j] : 0.0);
      if (!fin(v)) {
        unbounded_ = true;
        break;
      }
      fixColumn(j, v);
      stats_.emptyCols++;
      changed = true;
      continue;
    }
    // Dominated column / dual fixing.
    bool downSafe = true, upSafe = true;  // moving the column down (up) never hurts feasibility
    for (int k = A.start[j]; k < A.start[j + 1]; ++k) {
      int i = A.index[k];
      if (!rowActive_[i]) continue;
      double a = A.value[k];
      if (a > 0) {
        if (fin(rl_[i])) downSafe = false;
        if (fin(ru_[i])) upSafe = false;
      } else {
        if (fin(ru_[i])) downSafe = false;
        if (fin(rl_[i])) upSafe = false;
      }
    }
    if (c >= 0 && downSafe && fin(lo_[j])) {
      fixColumn(j, lo_[j]);
      stats_.dominatedCols++;
      changed = true;
      continue;
    }
    if (c <= 0 && upSafe && fin(up_[j])) {
      fixColumn(j, up_[j]);
      stats_.dominatedCols++;
      changed = true;
      continue;
    }
    // (Implied) free column singleton in an equality row.
    if (!isInt && colCount_[j] == 1) {
      int i = -1;
      double a = 0;
      for (int k = A.start[j]; k < A.start[j + 1]; ++k)
        if (rowActive_[A.index[k]]) {
          i = A.index[k];
          a = A.value[k];
        }
      if (i >= 0 && fin(rl_[i]) && rl_[i] == ru_[i] && std::fabs(a) > 1e-3) {
        // Implied bounds of x_j from the row with the other columns at their bounds.
        double restMin = 0, restMax = 0;
        bool rminInf = false, rmaxInf = false;
        for (int e = AT_.start[i]; e < AT_.start[i + 1]; ++e) {
          int k = AT_.index[e];
          if (!colActive_[k] || k == j) continue;
          double ak = AT_.value[e];
          double bmin = ak > 0 ? lo_[k] : up_[k], bmax = ak > 0 ? up_[k] : lo_[k];
          if (!fin(bmin)) rminInf = true; else restMin += ak * bmin;
          if (!fin(bmax)) rmaxInf = true; else restMax += ak * bmax;
        }
        double b = rl_[i];
        double impLo, impUp;
        if (a > 0) {
          impLo = rmaxInf ? -kInf : (b - restMax) / a;
          impUp = rminInf ? kInf : (b - restMin) / a;
        } else {
          impLo = rminInf ? -kInf : (b - restMin) / a;
          impUp = rmaxInf ? kInf : (b - restMax) / a;
        }
        const double t = 1e-9;
        bool impliedFree = (!fin(lo_[j]) || (fin(impLo) && impLo >= lo_[j] - t * (1 + std::fabs(lo_[j])))) &&
                           (!fin(up_[j]) || (fin(impUp) && impUp <= up_[j] + t * (1 + std::fabs(up_[j]))));
        if (impliedFree) {
          Record r;
          r.kind = Kind::FreeColSingleton;
          r.row = i;
          r.col = j;
          r.a = a;
          r.cost = c;
          r.rowLower = b;
          r.oldLower = lo_[j];
          r.oldUpper = up_[j];
          for (int e = AT_.start[i]; e < AT_.start[i + 1]; ++e) {
            int k = AT_.index[e];
            if (!colActive_[k] || k == j) continue;
            r.cols.push_back(k);
            r.vals.push_back(AT_.value[e]);
            cost_[k] -= c * AT_.value[e] / a;  // substitute x_j = (b - sum a_k x_k)/a
          }
          offset_ += c * b / a;
          colActive_[j] = 0;
          stats_.colsRemoved++;
          removeRow(i);
          stack_.push_back(std::move(r));
          stats_.freeColSingletons++;
          changed = true;
          continue;
        }
      }
    }
  }
  return changed && !infeasible_ && !unbounded_;
}

PresolveStatus Presolver::run(const Model& model, bool mip) {
  Timer t;
  orig_ = model;
  qp_ = model.isQp();
  const int m = model.numRows(), n = model.numCols();
  AT_ = model.A.transpose();
  lo_ = model.colLower;
  up_ = model.colUpper;
  rl_ = model.rowLower;
  ru_ = model.rowUpper;
  cost_ = model.colCost;
  offset_ = model.objOffset;
  rowActive_.assign(m, 1);
  colActive_.assign(n, 1);
  rowCount_.assign(m, 0);
  colCount_.assign(n, 0);
  for (int j = 0; j < n; ++j) colCount_[j] = model.A.start[j + 1] - model.A.start[j];
  for (int i = 0; i < m; ++i) rowCount_[i] = AT_.start[i + 1] - AT_.start[i];
  stack_.clear();
  forcingCosts_.clear();
  certRow_ = -1;
  certSign_ = 0;
  stats_ = PresolveStats();
  infeasible_ = unbounded_ = false;

  for (int p = 0; p < 20; ++p) {
    stats_.passes++;
    if (!pass(mip)) break;
  }
  stats_.seconds = t.seconds();
  if (infeasible_) return PresolveStatus::Infeasible;
  if (unbounded_) return PresolveStatus::Unbounded;

  // ---- Build the reduced model ----
  colOrig_.clear();
  rowOrig_.clear();
  std::vector<int> colNew(n, -1), rowNew(m, -1);
  for (int j = 0; j < n; ++j)
    if (colActive_[j]) {
      colNew[j] = static_cast<int>(colOrig_.size());
      colOrig_.push_back(j);
    }
  for (int i = 0; i < m; ++i)
    if (rowActive_[i]) {
      rowNew[i] = static_cast<int>(rowOrig_.size());
      rowOrig_.push_back(i);
    }
  reduced_ = Model();
  reduced_.name = model.name + "_presolved";
  reduced_.objOffset = offset_;
  for (int j : colOrig_) {
    reduced_.colCost.push_back(cost_[j]);
    reduced_.colLower.push_back(lo_[j]);
    reduced_.colUpper.push_back(up_[j]);
    reduced_.colType.push_back(model.colType[j]);
    reduced_.colNames.push_back(j < static_cast<int>(model.colNames.size()) ? model.colNames[j] : "");
  }
  for (int i : rowOrig_) {
    reduced_.rowLower.push_back(rl_[i]);
    reduced_.rowUpper.push_back(ru_[i]);
    reduced_.rowNames.push_back(i < static_cast<int>(model.rowNames.size()) ? model.rowNames[i] : "");
  }
  std::vector<int> ri, ci;
  std::vector<double> v;
  for (int j : colOrig_)
    for (int k = model.A.start[j]; k < model.A.start[j + 1]; ++k) {
      int i = model.A.index[k];
      if (!rowActive_[i]) continue;
      ri.push_back(rowNew[i]);
      ci.push_back(colNew[j]);
      v.push_back(model.A.value[k]);
    }
  reduced_.A = SparseMatrix::fromTriplets(static_cast<int>(rowOrig_.size()), static_cast<int>(colOrig_.size()), ri, ci, v);
  reduced_.Q = SparseMatrix(0, 0);
  if (qp_) {
    std::vector<int> qi, qj;
    std::vector<double> qv;
    for (int j : colOrig_)
      for (int k = model.Q.start[j]; k < model.Q.start[j + 1]; ++k) {
        int i = model.Q.index[k];
        if (colNew[i] < 0) continue;
        qi.push_back(colNew[i]);
        qj.push_back(colNew[j]);
        qv.push_back(model.Q.value[k]);
      }
    int nn = static_cast<int>(colOrig_.size());
    reduced_.Q = SparseMatrix::fromTriplets(nn, nn, qi, qj, qv);
  }
  if (colOrig_.empty() && rowOrig_.empty()) return PresolveStatus::Empty;
  return stack_.empty() ? PresolveStatus::Unchanged : PresolveStatus::Reduced;
}

Presolver::Solution Presolver::postsolve(const Solution& rs) const {
  const int n = orig_.numCols(), m = orig_.numRows();
  const SparseMatrix& A = orig_.A;
  Solution s;
  s.hasDual = rs.hasDual;
  s.hasBasis = rs.hasBasis;
  s.x.assign(n, 0.0);
  s.rowDual.assign(m, 0.0);
  s.reducedCost.assign(n, 0.0);
  s.colStatus.assign(n, BasisStatus::Lower);
  s.rowStatus.assign(m, BasisStatus::Basic);
  for (size_t k = 0; k < colOrig_.size(); ++k) {
    s.x[colOrig_[k]] = rs.x[k];
    if (rs.hasBasis) s.colStatus[colOrig_[k]] = rs.colStatus[k];
  }
  for (size_t r = 0; r < rowOrig_.size(); ++r) {
    if (rs.hasDual) s.rowDual[rowOrig_[r]] = rs.rowDual[r];
    if (rs.hasBasis) s.rowStatus[rowOrig_[r]] = rs.rowStatus[r];
  }
  // Stage reduced cost: stage cost minus contributions of rows restored so far.
  // QP: reduced cost d = c + Qx - A^T y (columns are never removed in QP mode, so x is known).
  std::vector<double> qx(n, 0.0);
  if (qp_) orig_.Q.multiply(s.x.data(), qx.data());
  auto stageD = [&](int j, double stageCost) {
    double d = stageCost + qx[j];
    for (int k = A.start[j]; k < A.start[j + 1]; ++k) d -= A.value[k] * s.rowDual[A.index[k]];
    return d;
  };
  auto boundStatus = [&](int j) {
    double l = orig_.colLower[j], u = orig_.colUpper[j], x = s.x[j];
    if (l == u) return BasisStatus::Fixed;
    if (fin(l) && std::fabs(x - l) <= 1e-9 * (1 + std::fabs(l))) return BasisStatus::Lower;
    if (fin(u) && std::fabs(x - u) <= 1e-9 * (1 + std::fabs(u))) return BasisStatus::Upper;
    if (!fin(l) && !fin(u)) return BasisStatus::Zero;
    return BasisStatus::Basic;  // strictly between bounds: must be basic
  };
  auto rowSideStatus = [&](int i, double act) {
    double l = orig_.rowLower[i], u = orig_.rowUpper[i];
    if (l == u) return BasisStatus::Fixed;
    if (fin(l) && std::fabs(act - l) <= std::fabs(act - u)) return BasisStatus::Lower;
    return BasisStatus::Upper;
  };
  auto rowAct = [&](int i) {
    double a = 0;
    for (int e = AT_.start[i]; e < AT_.start[i + 1]; ++e) a += AT_.value[e] * s.x[AT_.index[e]];
    return a;
  };

  for (auto it = stack_.rbegin(); it != stack_.rend(); ++it) {
    const Record& r = *it;
    switch (r.kind) {
      case Kind::FixedCol: {
        s.x[r.col] = r.value;
        s.colStatus[r.col] = boundStatus(r.col);
        if (s.colStatus[r.col] == BasisStatus::Basic) s.colStatus[r.col] = BasisStatus::Lower;
        break;
      }
      case Kind::EmptyRow:
      case Kind::RedundantRow: {
        s.rowDual[r.row] = 0.0;
        s.rowStatus[r.row] = BasisStatus::Basic;
        break;
      }
      case Kind::SingletonRow: {
        const int j = r.col, i = r.row;
        const double x = s.x[j];
        s.rowDual[i] = 0.0;
        s.rowStatus[i] = BasisStatus::Basic;
        // Bounds implied by this row (before intersection with the column's own bounds).
        double rowImpLo = r.a > 0 ? (fin(r.rowLower) ? r.rowLower / r.a : -kInf) : (fin(r.rowUpper) ? r.rowUpper / r.a : -kInf);
        double rowImpUp = r.a > 0 ? (fin(r.rowUpper) ? r.rowUpper / r.a : kInf) : (fin(r.rowLower) ? r.rowLower / r.a : kInf);
        const double tx = 1e-9 * (1 + std::fabs(x));
        bool atRowLo = fin(rowImpLo) && std::fabs(x - rowImpLo) <= tx &&
                       (!fin(r.oldLower) || rowImpLo > r.oldLower + tx);
        bool atRowUp = fin(rowImpUp) && std::fabs(x - rowImpUp) <= tx &&
                       (!fin(r.oldUpper) || rowImpUp < r.oldUpper - tx);
        const bool colNonbasic = !s.hasBasis || s.colStatus[j] != BasisStatus::Basic;
        if (!s.hasBasis) {
          // Interior / first-order solution (no basis): the sign of the reduced cost says
          // which bound's multiplier it carries; if that bound came from this row,
          // transfer it to the row dual (x need not sit exactly on the bound).
          double d = stageD(j, r.cost);
          bool lowerFromRow = fin(rowImpLo) && (!fin(r.oldLower) || rowImpLo > r.oldLower);
          bool upperFromRow = fin(rowImpUp) && (!fin(r.oldUpper) || rowImpUp < r.oldUpper);
          if ((d > 0 && lowerFromRow) || (d < 0 && upperFromRow)) s.rowDual[i] = d / r.a;
          break;
        }
        if (colNonbasic && (atRowLo || atRowUp)) {
          // x_j sits on a bound owned by the row: the row becomes nonbasic and takes
          // over the reduced cost (y_i = d_j / a), the column becomes basic.
          double d = stageD(j, r.cost);
          s.rowDual[i] = d / r.a;
          s.colStatus[j] = BasisStatus::Basic;
          s.rowStatus[i] = rowSideStatus(i, r.a * x);
        } else if (colNonbasic && s.hasBasis) {
          BasisStatus own = boundStatus(j);
          s.colStatus[j] = own == BasisStatus::Basic ? BasisStatus::Lower : own;
        }
        break;
      }
      case Kind::ForcingRow: {
        const int i = r.row;
        const std::vector<double>& stageCost = forcingCosts_[static_cast<size_t>(r.value)];
        for (size_t k = 0; k < r.cols.size(); ++k) s.x[r.cols[k]] = r.colVals[k];
        double y = 0;
        int arg = -1;
        for (size_t k = 0; k < r.cols.size(); ++k) {
          double d = stageD(r.cols[k], stageCost[k]);
          double ratio = d / r.vals[k];
          if (r.side == 1) {
            if (ratio < y) { y = ratio; arg = static_cast<int>(k); }
          } else {
            if (ratio > y) { y = ratio; arg = static_cast<int>(k); }
          }
        }
        s.rowDual[i] = y;
        for (size_t k = 0; k < r.cols.size(); ++k) {
          BasisStatus st = boundStatus(r.cols[k]);
          if (st == BasisStatus::Basic) st = BasisStatus::Lower;
          s.colStatus[r.cols[k]] = st;
        }
        if (arg >= 0) {
          s.colStatus[r.cols[arg]] = BasisStatus::Basic;
          s.rowStatus[i] = r.side == 1 ? (orig_.rowLower[i] == orig_.rowUpper[i] ? BasisStatus::Fixed : BasisStatus::Upper)
                                       : (orig_.rowLower[i] == orig_.rowUpper[i] ? BasisStatus::Fixed : BasisStatus::Lower);
        } else {
          s.rowStatus[i] = BasisStatus::Basic;
        }
        break;
      }
      case Kind::FreeColSingleton: {
        const int j = r.col, i = r.row;
        double rest = 0;
        for (size_t k = 0; k < r.cols.size(); ++k) rest += r.vals[k] * s.x[r.cols[k]];
        s.x[j] = (r.rowLower - rest) / r.a;
        s.rowDual[i] = r.cost / r.a;
        s.colStatus[j] = BasisStatus::Basic;
        s.rowStatus[i] = BasisStatus::Fixed;
        break;
      }
    }
  }
  if (std::getenv("PRESOLVE_DEBUG"))
    for (int j = 0; j < n; ++j)
      if (!std::isfinite(s.x[j]) || std::fabs(s.x[j]) > 1e20) PLOG_INFO("postsolve: x[%d] = %g", j, s.x[j]);
  // Final reduced costs on the original model.
  for (int j = 0; j < n; ++j) {
    double d = orig_.colCost[j] + qx[j];
    for (int k = A.start[j]; k < A.start[j + 1]; ++k) d -= A.value[k] * s.rowDual[A.index[k]];
    s.reducedCost[j] = d;
  }
  (void)rowAct;
  return s;
}

}  // namespace pramana
