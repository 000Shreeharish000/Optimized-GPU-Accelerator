#include "mip/cuts.h"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "lp/simplex.h"

namespace pramana {

namespace {
inline bool fin(double v) { return std::fabs(v) < kInfBoundThreshold; }
inline double frac(double v) { return v - std::floor(v); }
}  // namespace

CutGenerator::CutGenerator(const Model& model) : model_(model) {
  AT_ = model.A.transpose();
  isInt_.assign(model.numCols(), 0);
  for (int j = 0; j < model.numCols(); ++j) isInt_[j] = model.colType[j] == VarType::Integer;
}

bool CutGenerator::finalize(std::vector<double>& dense, double rhs, const std::vector<double>& x,
                            const std::vector<double>& lo, const std::vector<double>& up, const std::string& type,
                            Cut& out) {
  const int n = model_.numCols();
  double maxAbs = 0;
  for (int j = 0; j < n; ++j) maxAbs = std::max(maxAbs, std::fabs(dense[j]));
  if (maxAbs < 1e-9 || !std::isfinite(rhs)) return false;
  // Remove tiny coefficients conservatively: a_j x_j <= max over the box, so
  // dropping the term needs rhs -= max(a_j l_j, a_j u_j).
  for (int j = 0; j < n; ++j) {
    double a = dense[j];
    if (a == 0.0 || std::fabs(a) >= 1e-9 * maxAbs) continue;
    double m = std::max(fin(lo[j]) ? a * lo[j] : (a < 0 ? kInf : -kInf), fin(up[j]) ? a * up[j] : (a > 0 ? kInf : -kInf));
    if (!std::isfinite(m)) {
      stats.rejectedNumerics++;
      return false;
    }
    rhs -= m;
    dense[j] = 0.0;
  }
  double minAbs = kInf, norm2 = 0, act = 0;
  out.idx.clear();
  out.val.clear();
  for (int j = 0; j < n; ++j) {
    double a = dense[j];
    if (a == 0.0) continue;
    minAbs = std::min(minAbs, std::fabs(a));
    norm2 += a * a;
    act += a * x[j];
    out.idx.push_back(j);
    out.val.push_back(a);
  }
  if (out.idx.empty() || maxAbs / minAbs > 1e6) {
    stats.rejectedNumerics++;
    return false;
  }
  rhs -= 1e-9 * (1.0 + std::fabs(rhs));  // safety relaxation
  out.rhs = rhs;
  out.type = type;
  out.efficacy = (rhs - act) / std::sqrt(norm2);
  return out.efficacy > 0;
}

// ---------------------------------------------------------------------------
// Gomory mixed-integer cuts
// ---------------------------------------------------------------------------
std::vector<Cut> CutGenerator::gomory(const Simplex& lp, const std::vector<double>& x,
                                      const std::vector<double>& lo, const std::vector<double>& up, int maxCuts) {
  const int n = model_.numCols(), m = model_.numRows();
  std::vector<Cut> cuts;
  std::vector<BasisStatus> cs, rs;
  lp.getBasis(cs, rs);
  // Candidate rows: fractional basic integer structurals, most fractional first.
  std::vector<std::pair<double, int>> cand;
  for (int p = 0; p < lp.numRows(); ++p) {
    int v = lp.basicVariable(p);
    if (v >= n || !isInt_[v]) continue;
    double f = frac(x[v]);
    if (f < 0.005 || f > 0.995) continue;
    cand.emplace_back(-std::min(f, 1 - f), p);
  }
  std::sort(cand.begin(), cand.end());
  if (static_cast<int>(cand.size()) > 4 * maxCuts) cand.resize(4 * maxCuts);
  std::vector<double> t, dense(n);
  const int mTab = lp.numRows();  // may include earlier cut rows (logicals beyond m)
  for (auto& c : cand) {
    int p = c.second;
    int bv = lp.tableauRow(p, t);
    double beta = x[bv];
    double f0 = frac(beta);
    std::fill(dense.begin(), dense.end(), 0.0);
    double cst = 0;
    bool ok = true;
    for (int v = 0; v < n + mTab && ok; ++v) {
      double tv = t[v];
      if (tv == 0.0 || std::fabs(tv) < 1e-12) continue;
      // Nonbasic status and bounds of v (logical s_i = -row_i(x) has bounds [-U,-L]).
      BasisStatus st;
      double lv, uv;
      if (v < n) {
        st = cs[v];
        lv = lo[v];
        uv = up[v];
      } else {
        int i = v - n;
        st = rs[i];
        if (i >= m) {  // logical of a previously added cut row: skip row (not tracked here)
          ok = false;
          break;
        }
        lv = fin(model_.rowUpper[i]) ? -model_.rowUpper[i] : -kInf;
        uv = fin(model_.rowLower[i]) ? -model_.rowLower[i] : kInf;
        // Row status "Lower" means activity at L => logical at its upper (-L).
        if (st == BasisStatus::Lower) st = BasisStatus::Upper;
        else if (st == BasisStatus::Upper) st = BasisStatus::Lower;
      }
      if (st == BasisStatus::Basic) continue;
      if (st == BasisStatus::Fixed || (fin(lv) && lv == uv)) continue;  // constant, x' = 0
      if (st == BasisStatus::Zero) {
        ok = false;  // free nonbasic: no sign on x'
        break;
      }
      bool atLower = st == BasisStatus::Lower;
      double abar = atLower ? tv : -tv;
      bool intVar = v < n && isInt_[v] && fin(lv) && fin(uv) && lv == std::floor(lv) && uv == std::floor(uv);
      double g;
      if (intVar) {
        double fj = frac(abar);
        g = fj <= f0 ? fj / f0 : (1 - fj) / (1 - f0);
      } else {
        g = abar >= 0 ? abar / f0 : -abar / (1 - f0);
      }
      if (g == 0.0) continue;
      // Substitute x' back:  sum g x' >= 1.
      if (v < n) {
        if (atLower) {
          dense[v] += g;
          cst -= g * lv;
        } else {
          dense[v] -= g;
          cst += g * uv;
        }
      } else {
        int i = v - n;
        // s_i = -a_i^T x
        double sign = atLower ? -1.0 : 1.0;  // x' = s - l  or  u - s
        for (int e = AT_.start[i]; e < AT_.start[i + 1]; ++e) dense[AT_.index[e]] += sign * g * AT_.value[e];
        cst += atLower ? -g * lv : g * uv;
      }
    }
    if (!ok) continue;
    Cut cut;
    if (finalize(dense, 1.0 - cst, x, lo, up, "gmi", cut)) {
      cuts.push_back(std::move(cut));
      stats.gomory++;
    }
    if (static_cast<int>(cuts.size()) >= 2 * maxCuts) break;
  }
  return cuts;
}

// ---------------------------------------------------------------------------
// Complemented MIR on single rows
// ---------------------------------------------------------------------------
std::vector<Cut> CutGenerator::mir(const std::vector<double>& x, const std::vector<double>& lo,
                                   const std::vector<double>& up, int maxCuts) {
  const int n = model_.numCols(), m = model_.numRows();
  std::vector<Cut> cuts;
  std::vector<double> dense(n);
  for (int i = 0; i < m; ++i) {
    for (int side = 0; side < 2; ++side) {
      double b;
      double sgn;
      if (side == 0) {
        if (!fin(model_.rowUpper[i])) continue;
        b = model_.rowUpper[i];
        sgn = 1.0;
      } else {
        if (!fin(model_.rowLower[i])) continue;
        b = -model_.rowLower[i];
        sgn = -1.0;
      }
      // Bound substitution: x_j = l + x' (lower) or u - x' (upper, complemented).
      struct Term {
        int j;
        double a;       // coefficient of x' (>= 0 variable)
        bool comp;      // complemented (x = u - x')
        bool integer;
        double xp;      // value of x'
      };
      std::vector<Term> terms;
      double bb = b;
      bool ok = true;
      int nInt = 0;
      for (int e = AT_.start[i]; e < AT_.start[i + 1]; ++e) {
        int j = AT_.index[e];
        double a = sgn * AT_.value[e];
        bool integer = isInt_[j] && fin(lo[j]) && fin(up[j]);
        bool useLower;
        if (fin(lo[j]) && fin(up[j])) useLower = (x[j] - lo[j]) <= (up[j] - x[j]);
        else if (fin(lo[j])) useLower = true;
        else if (fin(up[j])) useLower = false;
        else {
          ok = false;
          break;
        }
        if (isInt_[j] && !integer) {
          ok = false;
          break;
        }
        if (useLower) {
          bb -= a * lo[j];
          terms.push_back({j, a, false, integer, x[j] - lo[j]});
        } else {
          bb -= a * up[j];
          terms.push_back({j, -a, true, integer, up[j] - x[j]});
        }
        nInt += integer;
      }
      if (!ok || nInt == 0) continue;
      // Continuous part: s = sum_{a<0} (-a) x' >= 0 ; positive continuous terms dropped.
      double sval = 0;
      for (auto& t : terms)
        if (!t.integer && t.a < 0) sval += -t.a * t.xp;
      // Candidate divisors.
      std::vector<double> deltas{1.0};
      for (auto& t : terms)
        if (t.integer && t.xp > 1e-6 && std::fabs(t.a) > 1e-6) deltas.push_back(std::fabs(t.a));
      std::sort(deltas.begin(), deltas.end());
      deltas.erase(std::unique(deltas.begin(), deltas.end()), deltas.end());
      if (deltas.size() > 8) deltas.resize(8);
      double bestViol = 1e-6;
      std::vector<double> bestDense;
      double bestRhs = 0;
      for (double delta : deltas) {
        double bt = bb / delta;
        double f = frac(bt);
        if (f < 0.05 || f > 0.95) continue;
        // MIR in x' space:  sum F(a/delta) x'_int - s/(delta (1-f)) <= floor(bt)
        double lhs = 0;
        for (auto& t : terms) {
          if (t.integer) {
            double at = t.a / delta;
            double fj = frac(at);
            double F = std::floor(at) + std::max(0.0, fj - f) / (1 - f);
            lhs += F * t.xp;
          }
        }
        lhs -= sval / (delta * (1 - f));
        double viol = lhs - std::floor(bt);
        if (viol <= bestViol) continue;
        bestViol = viol;
        // Build in original x space as a ">=" cut: -(MIR lhs) >= -floor(bt)
        std::fill(dense.begin(), dense.end(), 0.0);
        double rhsX = std::floor(bt);  // MIR: sum c'_j x'_j <= rhsX
        for (auto& t : terms) {
          double coef;
          if (t.integer) {
            double at = t.a / delta;
            double fj = frac(at);
            coef = std::floor(at) + std::max(0.0, fj - f) / (1 - f);
          } else {
            if (t.a >= 0) continue;
            coef = t.a / (delta * (1 - f));  // -(-a)/(delta(1-f))
          }
          if (coef == 0.0) continue;
          // x' = x - l  or  u - x
          if (!t.comp) {
            dense[t.j] += coef;
            rhsX += coef * lo[t.j];
          } else {
            dense[t.j] -= coef;
            rhsX -= coef * up[t.j];
          }
        }
        for (int j = 0; j < n; ++j) dense[j] = -dense[j];
        bestDense = dense;
        bestRhs = -rhsX;
      }
      if (bestDense.empty()) continue;
      Cut cut;
      if (finalize(bestDense, bestRhs, x, lo, up, "mir", cut)) {
        cuts.push_back(std::move(cut));
        stats.mir++;
      }
    }
    if (static_cast<int>(cuts.size()) >= 4 * maxCuts) break;
  }
  return cuts;
}

// ---------------------------------------------------------------------------
// Knapsack cover cuts
// ---------------------------------------------------------------------------
std::vector<Cut> CutGenerator::covers(const std::vector<double>& x, const std::vector<double>& lo,
                                      const std::vector<double>& up, int maxCuts) {
  const int n = model_.numCols(), m = model_.numRows();
  std::vector<Cut> cuts;
  std::vector<double> dense(n);
  for (int i = 0; i < m; ++i) {
    for (int side = 0; side < 2; ++side) {
      double b, sgn;
      if (side == 0) {
        if (!fin(model_.rowUpper[i])) continue;
        b = model_.rowUpper[i];
        sgn = 1;
      } else {
        if (!fin(model_.rowLower[i])) continue;
        b = -model_.rowLower[i];
        sgn = -1;
      }
      struct Item {
        int j;
        double a;
        bool comp;
        double xv;
      };
      std::vector<Item> items;
      bool ok = true;
      for (int e = AT_.start[i]; e < AT_.start[i + 1] && ok; ++e) {
        int j = AT_.index[e];
        double a = sgn * AT_.value[e];
        if (!isInt_[j] || lo[j] != 0.0 || up[j] != 1.0) {
          ok = false;
          break;
        }
        if (a > 0) items.push_back({j, a, false, x[j]});
        else if (a < 0) {
          b -= a;  // x = 1 - x̄
          items.push_back({j, -a, true, 1 - x[j]});
        }
      }
      if (!ok || items.size() < 2 || b < 0) continue;
      double total = 0;
      for (auto& it : items) total += it.a;
      if (total <= b + 1e-9) continue;  // no cover exists
      std::sort(items.begin(), items.end(), [](const Item& p, const Item& q) {
        double sp = (1 - p.xv) / p.a, sq = (1 - q.xv) / q.a;
        return sp < sq;
      });
      std::vector<Item> cover;
      double w = 0;
      for (auto& it : items) {
        cover.push_back(it);
        w += it.a;
        if (w > b + 1e-9) break;
      }
      if (w <= b + 1e-9) continue;
      // Extended cover: add items with a_j >= max cover coefficient.
      double amax = 0;
      for (auto& c : cover) amax = std::max(amax, c.a);
      std::vector<Item> ext = cover;
      for (auto& it : items) {
        bool in = false;
        for (auto& c : cover)
          if (c.j == it.j) in = true;
        if (!in && it.a >= amax) ext.push_back(it);
      }
      double lhs = 0;
      for (auto& c : ext) lhs += c.xv;
      double rhsC = static_cast<double>(cover.size()) - 1;
      if (lhs <= rhsC + 1e-6) continue;
      // sum x̃ <= |C|-1  ->  as ">=" cut on x:  -sum x̃ >= -(|C|-1)
      std::fill(dense.begin(), dense.end(), 0.0);
      double rhs = -rhsC;
      for (auto& c : ext) {
        if (!c.comp) {
          dense[c.j] -= 1.0;
        } else {  // x̃ = 1 - x
          dense[c.j] += 1.0;
          rhs += 1.0;
        }
      }
      Cut cut;
      if (finalize(dense, rhs, x, lo, up, "cover", cut)) {
        cuts.push_back(std::move(cut));
        stats.cover++;
      }
    }
    if (static_cast<int>(cuts.size()) >= 4 * maxCuts) break;
  }
  return cuts;
}

std::vector<Cut> CutGenerator::select(std::vector<Cut> cuts, int maxCuts, double minEfficacy, double maxParallel) {
  std::sort(cuts.begin(), cuts.end(), [](const Cut& a, const Cut& b) { return a.efficacy > b.efficacy; });
  std::vector<Cut> out;
  std::vector<double> norms;
  for (auto& c : cuts) {
    if (c.efficacy < minEfficacy) break;
    double nc = 0;
    for (double v : c.val) nc += v * v;
    nc = std::sqrt(nc);
    bool parallel = false;
    for (size_t k = 0; k < out.size() && !parallel; ++k) {
      // sparse dot product (idx sorted ascending)
      const Cut& o = out[k];
      double dot = 0;
      size_t p = 0, q = 0;
      while (p < c.idx.size() && q < o.idx.size()) {
        if (c.idx[p] == o.idx[q]) dot += c.val[p++] * o.val[q++];
        else if (c.idx[p] < o.idx[q]) ++p;
        else ++q;
      }
      if (std::fabs(dot) / (nc * norms[k]) > maxParallel) parallel = true;
    }
    if (parallel) continue;
    out.push_back(c);
    norms.push_back(nc);
    if (static_cast<int>(out.size()) >= maxCuts) break;
  }
  return out;
}

}  // namespace pramana
