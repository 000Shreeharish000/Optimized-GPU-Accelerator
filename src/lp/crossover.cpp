#include "lp/crossover.h"

#include <algorithm>
#include <cmath>

namespace pramana {

namespace {
inline bool fin(double v) { return std::fabs(v) < kInfBoundThreshold; }
}

void guessBasis(const Model& model, const std::vector<double>& x, const std::vector<double>& rowDual,
                std::vector<BasisStatus>& cs, std::vector<BasisStatus>& rs) {
  const int n = model.numCols(), m = model.numRows();
  std::vector<double> act = model.rowActivity(x);
  std::vector<double> d = model.colCost;
  if (static_cast<int>(rowDual.size()) == m) model.A.multiplyTranspose(rowDual.data(), d.data(), -1.0);
  cs.assign(n, BasisStatus::Lower);
  rs.assign(m, BasisStatus::Basic);
  struct Cand {
    double score;
    int idx;
  };
  std::vector<Cand> cands;
  cands.reserve(n + m);
  for (int j = 0; j < n; ++j) {
    double l = model.colLower[j], u = model.colUpper[j], v = x[j];
    if (l == u) {
      cs[j] = BasisStatus::Fixed;
      continue;
    }
    double dl = fin(l) ? v - l : kInf, du = fin(u) ? u - v : kInf;
    if (!fin(l) && !fin(u)) cs[j] = BasisStatus::Zero;
    else cs[j] = dl <= du ? BasisStatus::Lower : BasisStatus::Upper;
    double dist = std::min(dl, du) / (1.0 + std::fabs(v));
    double score = (!fin(l) && !fin(u)) ? kInf : dist / (std::fabs(d[j]) + 1e-12);
    cands.push_back({score, j});
  }
  for (int i = 0; i < m; ++i) {
    double l = model.rowLower[i], u = model.rowUpper[i], v = act[i];
    double dl = fin(l) ? v - l : kInf, du = fin(u) ? u - v : kInf;
    double yi = static_cast<int>(rowDual.size()) == m ? std::fabs(rowDual[i]) : 0.0;
    double dist = std::min(dl, du) / (1.0 + std::fabs(v));
    double score = (!fin(l) && !fin(u)) ? kInf : dist / (yi + 1e-12);
    rs[i] = l == u ? BasisStatus::Fixed : (dl <= du ? BasisStatus::Lower : BasisStatus::Upper);
    cands.push_back({score, n + i});
  }
  std::stable_sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) {
    if (a.score != b.score) return a.score > b.score;
    return a.idx > b.idx;  // prefer logicals on ties
  });
  for (int k = 0; k < m && k < static_cast<int>(cands.size()); ++k) {
    int idx = cands[k].idx;
    if (idx < n) cs[idx] = BasisStatus::Basic;
    else rs[idx - n] = BasisStatus::Basic;
  }
}

LpResult crossover(const Model& model, const std::vector<double>& x, const std::vector<double>& rowDual,
                   const SimplexOptions& opts, const Deadline* deadline) {
  std::vector<BasisStatus> cs, rs;
  guessBasis(model, x, rowDual, cs, rs);
  return solveLpSimplex(model, opts, deadline, &cs, &rs);
}

}  // namespace pramana
