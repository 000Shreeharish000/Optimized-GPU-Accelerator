#include "ipm/ldl.h"

#include <algorithm>
#include <cmath>

#include "util/common.h"

namespace pramana {

// ---------------------------------------------------------------------------
// Approximate minimum degree on the quotient graph.
// ---------------------------------------------------------------------------
std::vector<int> minimumDegreeOrdering(const SparseMatrix& S) {
  const int n = S.numCols;
  std::vector<std::vector<int>> adjV(n), adjE(n), elem(n);
  {
    // Symmetrize the pattern and drop the diagonal.
    std::vector<std::vector<int>> tmp(n);
    for (int j = 0; j < n; ++j)
      for (int k = S.start[j]; k < S.start[j + 1]; ++k) {
        int i = S.index[k];
        if (i == j) continue;
        tmp[i].push_back(j);
        tmp[j].push_back(i);
      }
    for (int i = 0; i < n; ++i) {
      auto& v = tmp[i];
      std::sort(v.begin(), v.end());
      v.erase(std::unique(v.begin(), v.end()), v.end());
      adjV[i] = std::move(v);
    }
  }
  std::vector<int> degree(n), head(n + 1, -1), next(n, -1), prev(n, -1);
  std::vector<char> eliminated(n, 0), deadElem(n, 0);
  std::vector<int> mark(n, -1);
  auto bucketInsert = [&](int i, int d) {
    d = std::min(d, n);
    degree[i] = d;
    next[i] = head[d];
    prev[i] = -1;
    if (head[d] >= 0) prev[head[d]] = i;
    head[d] = i;
  };
  auto bucketRemove = [&](int i) {
    int d = degree[i];
    if (prev[i] >= 0) next[prev[i]] = next[i];
    else head[d] = next[i];
    if (next[i] >= 0) prev[next[i]] = prev[i];
  };
  for (int i = 0; i < n; ++i) bucketInsert(i, static_cast<int>(adjV[i].size()));

  std::vector<int> perm;
  perm.reserve(n);
  int minDeg = 0;
  std::vector<int> Lp;
  for (int k = 0; k < n; ++k) {
    while (minDeg <= n && head[minDeg] < 0) ++minDeg;
    int p = head[minDeg];
    bucketRemove(p);
    eliminated[p] = 1;
    perm.push_back(p);

    // Lp = adjV[p] U (union of elements of p) minus eliminated nodes.
    Lp.clear();
    mark[p] = k;
    for (int v : adjV[p])
      if (!eliminated[v] && mark[v] != k) {
        mark[v] = k;
        Lp.push_back(v);
      }
    for (int e : adjE[p]) {
      if (deadElem[e]) continue;
      for (int v : elem[e])
        if (!eliminated[v] && mark[v] != k) {
          mark[v] = k;
          Lp.push_back(v);
        }
      deadElem[e] = 1;  // absorbed into p
      elem[e].clear();
      elem[e].shrink_to_fit();
    }
    elem[p] = Lp;
    adjV[p].clear();
    adjV[p].shrink_to_fit();
    adjE[p].clear();

    const int remaining = n - k - 1;
    for (int i : Lp) {
      // Elements: drop dead ones, add p.
      auto& E = adjE[i];
      size_t w = 0;
      for (int e : E)
        if (!deadElem[e]) E[w++] = e;
      E.resize(w);
      E.push_back(p);
      // Variables: drop eliminated and those now covered by element p (in Lp).
      auto& V = adjV[i];
      w = 0;
      for (int v : V)
        if (!eliminated[v] && mark[v] != k) V[w++] = v;
      V.resize(w);
      // Approximate external degree.
      long long d = static_cast<long long>(V.size());
      for (int e : E) d += static_cast<long long>(elem[e].size()) - 1;
      d = std::min<long long>(d, remaining);
      bucketRemove(i);
      bucketInsert(i, static_cast<int>(std::max<long long>(d, 0)));
      if (degree[i] < minDeg) minDeg = degree[i];
    }
  }
  return perm;
}

// ---------------------------------------------------------------------------
// Symbolic analysis (elimination tree + column counts), Davis' up-looking LDL.
// ---------------------------------------------------------------------------
void LdlFactor::analyze(const SparseMatrix& K, const std::vector<int>& perm) {
  n_ = K.numCols;
  perm_ = perm;
  pinv_.assign(n_, 0);
  for (int k = 0; k < n_; ++k) pinv_[perm_[k]] = k;
  parent_.assign(n_, -1);
  lnz_.assign(n_, 0);
  std::vector<int> flag(n_, -1);
  for (int k = 0; k < n_; ++k) {
    parent_[k] = -1;
    flag[k] = k;
    int kk = perm_[k];
    for (int p = K.start[kk]; p < K.start[kk + 1]; ++p) {
      int i = pinv_[K.index[p]];
      if (i >= k) continue;
      for (; flag[i] != k; i = parent_[i]) {
        if (parent_[i] == -1) parent_[i] = k;
        lnz_[i]++;
        flag[i] = k;
      }
    }
  }
  lp_.assign(n_ + 1, 0);
  for (int k = 0; k < n_; ++k) lp_[k + 1] = lp_[k] + lnz_[k];
  li_.assign(lp_[n_], 0);
  lx_.assign(lp_[n_], 0.0);
  d_.assign(n_, 0.0);
  work_.assign(n_, 0.0);
}

int LdlFactor::factorize(const SparseMatrix& K, const std::vector<int8_t>& signs, double reg) {
  std::vector<double> Y(n_, 0.0);
  std::vector<int> pattern(n_), flag(n_, -1);
  std::fill(lnz_.begin(), lnz_.end(), 0);
  int regularized = 0;
  for (int k = 0; k < n_; ++k) {
    Y[k] = 0.0;
    int top = n_;
    flag[k] = k;
    int kk = perm_[k];
    for (int p = K.start[kk]; p < K.start[kk + 1]; ++p) {
      int i = pinv_[K.index[p]];
      if (i > k) continue;
      Y[i] += K.value[p];
      int len = 0;
      for (; flag[i] != k; i = parent_[i]) {
        pattern[len++] = i;
        flag[i] = k;
      }
      while (len > 0) pattern[--top] = pattern[--len];
    }
    double dk = Y[k];
    Y[k] = 0.0;
    for (; top < n_; ++top) {
      int i = pattern[top];
      double yi = Y[i];
      Y[i] = 0.0;
      int p2 = lp_[i] + lnz_[i];
      for (int p = lp_[i]; p < p2; ++p) Y[li_[p]] -= lx_[p] * yi;
      double lki = yi / d_[i];
      dk -= lki * yi;
      li_[p2] = k;
      lx_[p2] = lki;
      lnz_[i]++;
    }
    int s = signs[kk];
    if (reg > 0) {
      if (s < 0 ? dk > -reg : dk < reg) {
        dk = s < 0 ? -reg : reg;
        ++regularized;
      }
    } else if (dk == 0.0) {  // inertia mode: keep true pivots, avoid division by zero
      dk = 1e-300;
      ++regularized;
    }
    d_[k] = dk;
  }
  return regularized;
}

void LdlFactor::solve(std::vector<double>& b) const {
  std::vector<double>& x = work_;
  for (int k = 0; k < n_; ++k) x[k] = b[perm_[k]];
  for (int j = 0; j < n_; ++j) {
    double xj = x[j];
    if (xj == 0.0) continue;
    for (int p = lp_[j]; p < lp_[j + 1]; ++p) x[li_[p]] -= lx_[p] * xj;
  }
  for (int j = 0; j < n_; ++j) x[j] /= d_[j];
  for (int j = n_ - 1; j >= 0; --j) {
    double s = x[j];
    for (int p = lp_[j]; p < lp_[j + 1]; ++p) s -= lx_[p] * x[li_[p]];
    x[j] = s;
  }
  for (int k = 0; k < n_; ++k) b[perm_[k]] = x[k];
}

void LdlFactor::solveUnitLt(int k, std::vector<double>& out) const {
  std::vector<double> z(n_, 0.0);
  z[k] = 1.0;
  for (int j = n_ - 1; j >= 0; --j) {
    double sum = z[j];
    for (int p = lp_[j]; p < lp_[j + 1]; ++p) sum -= lx_[p] * z[li_[p]];
    z[j] = sum;
  }
  out.assign(n_, 0.0);
  for (int q = 0; q < n_; ++q) out[perm_[q]] = z[q];
}

int LdlFactor::negativePivots() const {
  int c = 0;
  for (double d : d_) c += d < 0;
  return c;
}

}  // namespace pramana
