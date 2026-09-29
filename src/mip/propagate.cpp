#include "mip/propagate.h"

#include <algorithm>
#include <cmath>

namespace pramana {

namespace {
inline bool fin(double v) { return std::fabs(v) < kInfBoundThreshold; }
}

void Propagator::setup(const Model& model) {
  m_ = model.numRows();
  n_ = model.numCols();
  A_ = model.A;
  AT_ = model.A.transpose();
  rl_ = model.rowLower;
  ru_ = model.rowUpper;
  isInt_.assign(n_, 0);
  for (int j = 0; j < n_; ++j) isInt_[j] = model.colType[j] == VarType::Integer;
}

PropagationResult Propagator::propagate(std::vector<double>& lo, std::vector<double>& up, int maxRounds) const {
  PropagationResult res;
  std::vector<char> inQueue(m_, 1);
  std::vector<int> queue(m_);
  for (int i = 0; i < m_; ++i) queue[i] = i;
  const double intTol = 1e-6;
  for (int round = 0; round < maxRounds && !queue.empty(); ++round) {
    res.rounds++;
    std::vector<int> next;
    for (int i : queue) {
      inQueue[i] = 0;
      double minFin = 0, maxFin = 0;
      int minInf = 0, maxInf = 0;
      for (int e = AT_.start[i]; e < AT_.start[i + 1]; ++e) {
        int j = AT_.index[e];
        double a = AT_.value[e];
        double bmin = a > 0 ? lo[j] : up[j];
        double bmax = a > 0 ? up[j] : lo[j];
        if (fin(bmin)) minFin += a * bmin; else ++minInf;
        if (fin(bmax)) maxFin += a * bmax; else ++maxInf;
      }
      const double L = rl_[i], U = ru_[i];
      const double scale = 1.0 + std::max(fin(L) ? std::fabs(L) : 0.0, fin(U) ? std::fabs(U) : 0.0);
      if (minInf == 0 && fin(U) && minFin > U + 1e-6 * scale) {
        res.infeasible = true;
        return res;
      }
      if (maxInf == 0 && fin(L) && maxFin < L - 1e-6 * scale) {
        res.infeasible = true;
        return res;
      }
      for (int e = AT_.start[i]; e < AT_.start[i + 1]; ++e) {
        int j = AT_.index[e];
        double a = AT_.value[e];
        // Residual activities excluding j.
        double bmin = a > 0 ? lo[j] : up[j];
        double bmax = a > 0 ? up[j] : lo[j];
        double resMin, resMax;
        bool resMinOk, resMaxOk;
        if (fin(bmin)) {
          resMinOk = minInf == 0;
          resMin = minFin - a * bmin;
        } else {
          resMinOk = minInf == 1;
          resMin = minFin;
        }
        if (fin(bmax)) {
          resMaxOk = maxInf == 0;
          resMax = maxFin - a * bmax;
        } else {
          resMaxOk = maxInf == 1;
          resMax = maxFin;
        }
        double newLo = -kInf, newUp = kInf;
        if (fin(U) && resMinOk) {
          double lim = (U - resMin) / a;  // a x_j <= U - resMin
          if (a > 0) newUp = lim;
          else newLo = lim;
        }
        if (fin(L) && resMaxOk) {
          double lim = (L - resMax) / a;  // a x_j >= L - resMax
          if (a > 0) newLo = std::max(newLo, lim);
          else newUp = std::min(newUp, lim);
        }
        bool changed = false;
        if (isInt_[j]) {
          if (fin(newLo)) newLo = std::ceil(newLo - intTol);
          if (fin(newUp)) newUp = std::floor(newUp + intTol);
          if (newLo > lo[j]) {
            lo[j] = newLo;
            changed = true;
          }
          if (newUp < up[j]) {
            up[j] = newUp;
            changed = true;
          }
        } else {
          // Continuous: only significant tightenings, relaxed by a safety margin.
          double range = (fin(lo[j]) && fin(up[j])) ? up[j] - lo[j] : kInf;
          if (fin(newLo)) {
            double margin = 1e-7 * (1 + std::fabs(newLo));
            double cand = newLo - margin;
            if (cand > lo[j] + 1e-3 * std::min(range, 1.0 + std::fabs(cand))) {
              lo[j] = cand;
              changed = true;
            }
          }
          if (fin(newUp)) {
            double margin = 1e-7 * (1 + std::fabs(newUp));
            double cand = newUp + margin;
            double r2 = (fin(lo[j]) && fin(up[j])) ? up[j] - lo[j] : kInf;
            if (cand < up[j] - 1e-3 * std::min(r2, 1.0 + std::fabs(cand))) {
              up[j] = cand;
              changed = true;
            }
          }
        }
        if (lo[j] > up[j] + 1e-6 * (1 + std::fabs(lo[j]))) {
          res.infeasible = true;
          return res;
        }
        if (lo[j] > up[j]) lo[j] = up[j];
        if (changed) {
          res.tightened++;
          for (int k = A_.start[j]; k < A_.start[j + 1]; ++k) {
            int r = A_.index[k];
            if (r != i && !inQueue[r]) {
              inQueue[r] = 1;
              next.push_back(r);
            }
          }
          // Refresh activities of the current row for subsequent variables.
          minFin = 0;
          maxFin = 0;
          minInf = maxInf = 0;
          for (int e2 = AT_.start[i]; e2 < AT_.start[i + 1]; ++e2) {
            int j2 = AT_.index[e2];
            double a2 = AT_.value[e2];
            double b1 = a2 > 0 ? lo[j2] : up[j2];
            double b2 = a2 > 0 ? up[j2] : lo[j2];
            if (fin(b1)) minFin += a2 * b1; else ++minInf;
            if (fin(b2)) maxFin += a2 * b2; else ++maxInf;
          }
        }
      }
    }
    queue.swap(next);
  }
  return res;
}

}  // namespace pramana
