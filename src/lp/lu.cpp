#include "lp/lu.h"

#include <algorithm>
#include <cmath>
#include <utility>

#include "util/common.h"

namespace pramana {

namespace {

// Doubly linked lists of items bucketed by their current count.
struct CountLists {
  std::vector<int> head, next, prev, count;
  void init(int items, int maxCount) {
    head.assign(maxCount + 2, -1);
    next.assign(items, -1);
    prev.assign(items, -1);
    count.assign(items, 0);
  }
  void insert(int j, int c) {
    count[j] = c;
    next[j] = head[c];
    prev[j] = -1;
    if (head[c] >= 0) prev[head[c]] = j;
    head[c] = j;
  }
  void remove(int j) {
    int c = count[j];
    if (prev[j] >= 0) next[prev[j]] = next[j];
    else head[c] = next[j];
    if (next[j] >= 0) prev[next[j]] = prev[j];
    next[j] = prev[j] = -1;
  }
  void change(int j, int c) {
    remove(j);
    insert(j, c);
  }
};

}  // namespace

void LuFactor::setup(int numRows, int numStructural, const SparseMatrix* A) {
  m_ = numRows;
  n_ = numStructural;
  A_ = A;
  work_.assign(m_, 0.0);
  mark_.assign(m_, 0);
  clearFactors();
}

void LuFactor::clearFactors() {
  pivotRow_.clear();
  pivotPos_.clear();
  pivotValue_.clear();
  lStart_.assign(1, 0);
  lIndex_.clear();
  lValue_.clear();
  uStart_.assign(1, 0);
  uIndex_.clear();
  uValue_.clear();
  etaPos_.clear();
  etaStart_.assign(1, 0);
  etaIndex_.clear();
  etaPivot_.clear();
  etaValue_.clear();
}

int LuFactor::factorize(const std::vector<int>& basicIndex) {
  const int m = m_;
  clearFactors();
  deficientPositions.clear();
  deficientRows.clear();
  stats.factorizations++;
  if (static_cast<int>(work_.size()) != m) {
    work_.assign(m, 0.0);
    mark_.assign(m, 0);
  }

  // ---- Build the active matrix (column entries + row patterns). ----
  std::vector<std::vector<std::pair<int, double>>> col(m);
  std::vector<std::vector<int>> rowPat(m);
  double maxOrig = 0;
  for (int p = 0; p < m; ++p) {
    int v = basicIndex[p];
    if (v >= n_) {
      col[p].emplace_back(v - n_, 1.0);
    } else {
      for (int k = A_->start[v]; k < A_->start[v + 1]; ++k)
        if (A_->value[k] != 0.0) col[p].emplace_back(A_->index[k], A_->value[k]);
    }
    for (auto& e : col[p]) {
      rowPat[e.first].push_back(p);
      maxOrig = std::max(maxOrig, std::fabs(e.second));
    }
  }

  CountLists colList, rowList;
  colList.init(m, m);
  rowList.init(m, m);
  for (int p = 0; p < m; ++p) colList.insert(p, static_cast<int>(col[p].size()));
  for (int r = 0; r < m; ++r) rowList.insert(r, static_cast<int>(rowPat[r].size()));
  std::vector<char> colDone(m, 0), rowDone(m, 0);

  auto colMax = [&](int p) {
    double mx = 0;
    for (auto& e : col[p]) mx = std::max(mx, std::fabs(e.second));
    return mx;
  };
  auto valueAt = [&](int p, int r) -> double {
    for (auto& e : col[p])
      if (e.first == r) return e.second;
    return 0.0;
  };

  std::vector<int> markPos(m, -1);
  std::vector<std::pair<int, double>> uRow, lCol;
  double maxU = 0;
  int bump = 0;

  // Empty columns (all-zero basis columns) are immediately deficient.
  for (int step = 0; step < m; ++step) {
    // ---- Markowitz pivot search. ----
    int bestRow = -1, bestCol = -1;
    double bestCost = 1e300, bestAbs = 0;
    int searched = 0;
    bool done = false;
    for (int cnt = 1; cnt <= m && !done; ++cnt) {
      const double colBound = static_cast<double>(cnt - 1) * (cnt - 1);
      const double rowBound = static_cast<double>(cnt - 1) * cnt;
      for (int p = colList.head[cnt]; p >= 0 && !done; p = colList.next[p]) {
        double mx = colMax(p);
        double thr = std::max(pivotThreshold * mx, absPivotTol);
        for (auto& e : col[p]) {
          double a = std::fabs(e.second);
          if (a < thr) continue;
          double cost = static_cast<double>(rowList.count[e.first] - 1) * (cnt - 1);
          if (cost < bestCost || (cost == bestCost && a > bestAbs)) {
            bestCost = cost;
            bestRow = e.first;
            bestCol = p;
            bestAbs = a;
          }
        }
        ++searched;
        if (bestRow >= 0 && (searched >= 4 || bestCost <= colBound)) done = true;
      }
      for (int r = rowList.head[cnt]; r >= 0 && !done; r = rowList.next[r]) {
        for (int p : rowPat[r]) {
          double a = std::fabs(valueAt(p, r));
          if (a < std::max(pivotThreshold * colMax(p), absPivotTol)) continue;
          double cost = static_cast<double>(cnt - 1) * (colList.count[p] - 1);
          if (cost < bestCost || (cost == bestCost && a > bestAbs)) {
            bestCost = cost;
            bestRow = r;
            bestCol = p;
            bestAbs = a;
          }
        }
        ++searched;
        if (bestRow >= 0 && (searched >= 4 || bestCost <= rowBound)) done = true;
      }
    }
    if (bestRow < 0) break;  // remaining active submatrix is (numerically) singular

    const int r = bestRow, c = bestCol;
    const double piv = valueAt(c, r);
    if (bestCost > 0) ++bump;

    // ---- Extract U row (row r without column c) and remove row r from columns. ----
    uRow.clear();
    for (int p : rowPat[r]) {
      if (p == c) continue;
      auto& ce = col[p];
      for (size_t k = 0; k < ce.size(); ++k) {
        if (ce[k].first == r) {
          uRow.emplace_back(p, ce[k].second);
          maxU = std::max(maxU, std::fabs(ce[k].second));
          ce[k] = ce.back();
          ce.pop_back();
          break;
        }
      }
    }
    // ---- Extract L column (column c without row r) and remove column c from rows. ----
    lCol.clear();
    for (auto& e : col[c]) {
      if (e.first == r) continue;
      lCol.emplace_back(e.first, e.second / piv);
      auto& rp = rowPat[e.first];
      for (size_t k = 0; k < rp.size(); ++k)
        if (rp[k] == c) {
          rp[k] = rp.back();
          rp.pop_back();
          break;
        }
    }
    colDone[c] = 1;
    rowDone[r] = 1;
    colList.remove(c);
    rowList.remove(r);
    col[c].clear();
    rowPat[r].clear();

    // ---- Schur complement update: col_j -= u_rj * l for each j in U row. ----
    if (!lCol.empty()) {
      for (auto& ue : uRow) {
        int p = ue.first;
        double u = ue.second;
        auto& ce = col[p];
        for (size_t k = 0; k < ce.size(); ++k) markPos[ce[k].first] = static_cast<int>(k);
        for (auto& le : lCol) {
          double delta = -le.second * u;
          int i = le.first;
          if (markPos[i] >= 0) {
            ce[markPos[i]].second += delta;
          } else {
            ce.emplace_back(i, delta);
            markPos[i] = static_cast<int>(ce.size()) - 1;
            rowPat[i].push_back(p);
          }
        }
        for (auto& e : ce) markPos[e.first] = -1;
      }
    }
    for (auto& ue : uRow) colList.change(ue.first, static_cast<int>(col[ue.first].size()));
    for (auto& le : lCol) rowList.change(le.first, static_cast<int>(rowPat[le.first].size()));

    // ---- Record factors. ----
    pivotRow_.push_back(r);
    pivotPos_.push_back(c);
    pivotValue_.push_back(piv);
    for (auto& le : lCol) {
      if (le.second == 0.0) continue;
      lIndex_.push_back(le.first);
      lValue_.push_back(le.second);
    }
    lStart_.push_back(static_cast<int>(lIndex_.size()));
    for (auto& ue : uRow) {
      if (ue.second == 0.0) continue;
      uIndex_.push_back(ue.first);
      uValue_.push_back(ue.second);
    }
    uStart_.push_back(static_cast<int>(uIndex_.size()));
  }

  int K = static_cast<int>(pivotRow_.size());
  if (K < m) {
    for (int p = 0; p < m; ++p)
      if (!colDone[p]) deficientPositions.push_back(p);
    for (int r = 0; r < m; ++r)
      if (!rowDone[r]) deficientRows.push_back(r);
    stats.rankDeficiencies++;
  }
  posStep_.assign(m, -1);
  for (int k = 0; k < K; ++k) posStep_[pivotPos_[k]] = k;
  buildUColumns();
  stats.lNnz = static_cast<long long>(lIndex_.size());
  stats.uNnz = static_cast<long long>(uIndex_.size()) + K;
  stats.bumpSize = bump;
  stats.lastGrowth = maxOrig > 0 ? maxU / maxOrig : 0;
  return m - K;
}

void LuFactor::buildUColumns() {
  int K = static_cast<int>(pivotRow_.size());
  ucStart_.assign(K + 1, 0);
  for (int k = 0; k < K; ++k)
    for (int e = uStart_[k]; e < uStart_[k + 1]; ++e) {
      int s = posStep_[uIndex_[e]];
      if (s >= 0) ucStart_[s + 1]++;
    }
  for (int k = 0; k < K; ++k) ucStart_[k + 1] += ucStart_[k];
  ucIndex_.assign(ucStart_[K], 0);
  ucValue_.assign(ucStart_[K], 0.0);
  std::vector<int> pos(ucStart_.begin(), ucStart_.end() - 1);
  for (int k = 0; k < K; ++k)
    for (int e = uStart_[k]; e < uStart_[k + 1]; ++e) {
      int s = posStep_[uIndex_[e]];
      if (s < 0) continue;
      int q = pos[s]++;
      ucIndex_[q] = k;
      ucValue_[q] = uValue_[e];
    }
}

void LuFactor::ftran(HVector& v) const {
  double* b = v.array.data();
  const int K = static_cast<int>(pivotRow_.size());
  // Forward elimination with L.
  for (int k = 0; k < K; ++k) {
    double br = b[pivotRow_[k]];
    if (br == 0.0) continue;
    for (int e = lStart_[k]; e < lStart_[k + 1]; ++e) b[lIndex_[e]] -= lValue_[e] * br;
  }
  // Back substitution with U (column oriented), result indexed by basic position.
  double* x = work_.data();
  for (int k = K - 1; k >= 0; --k) {
    int r = pivotRow_[k];
    double br = b[r];
    if (br == 0.0) continue;
    b[r] = 0.0;
    double xv = br / pivotValue_[k];
    x[pivotPos_[k]] = xv;
    for (int e = ucStart_[k]; e < ucStart_[k + 1]; ++e) b[pivotRow_[ucIndex_[e]]] -= ucValue_[e] * xv;
  }
  // Rows never pivoted (deficient factor) may hold residue: clear them.
  if (K < m_)
    for (int r : deficientRows) b[r] = 0.0;
  v.array.swap(work_);
  // Product-form etas.
  double* y = v.array.data();
  const int T = static_cast<int>(etaPos_.size());
  for (int t = 0; t < T; ++t) {
    int p = etaPos_[t];
    double xp = y[p];
    if (xp == 0.0) continue;
    xp /= etaPivot_[t];
    y[p] = xp;
    for (int e = etaStart_[t]; e < etaStart_[t + 1]; ++e) y[etaIndex_[e]] -= etaValue_[e] * xp;
  }
  v.rebuildIndex(1e-14);
}

void LuFactor::btran(HVector& v) const {
  double* e = v.array.data();
  const int T = static_cast<int>(etaPos_.size());
  for (int t = T - 1; t >= 0; --t) {
    int p = etaPos_[t];
    double s = e[p];
    for (int q = etaStart_[t]; q < etaStart_[t + 1]; ++q) s -= etaValue_[q] * e[etaIndex_[q]];
    e[p] = s / etaPivot_[t];
  }
  const int K = static_cast<int>(pivotRow_.size());
  double* w = work_.data();
  for (int k = 0; k < K; ++k) {
    int c = pivotPos_[k];
    double val = e[c];
    if (val == 0.0) continue;
    e[c] = 0.0;
    double wv = val / pivotValue_[k];
    w[pivotRow_[k]] = wv;
    for (int q = uStart_[k]; q < uStart_[k + 1]; ++q) e[uIndex_[q]] -= uValue_[q] * wv;
  }
  if (K < m_)
    for (int p : deficientPositions) e[p] = 0.0;
  for (int k = K - 1; k >= 0; --k) {
    int r = pivotRow_[k];
    double s = w[r];
    for (int q = lStart_[k]; q < lStart_[k + 1]; ++q) s -= lValue_[q] * w[lIndex_[q]];
    w[r] = s;
  }
  v.array.swap(work_);
  v.rebuildIndex(1e-14);
}

void LuFactor::update(const HVector& column, int pos) {
  etaPos_.push_back(pos);
  etaPivot_.push_back(column.array[pos]);
  for (int k = 0; k < column.count; ++k) {
    int i = column.index[k];
    if (i == pos) continue;
    double a = column.array[i];
    if (std::fabs(a) <= 1e-14) continue;
    etaIndex_.push_back(i);
    etaValue_.push_back(a);
  }
  etaStart_.push_back(static_cast<int>(etaIndex_.size()));
}

}  // namespace pramana
