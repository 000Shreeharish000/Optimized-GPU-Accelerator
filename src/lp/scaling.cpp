#include "lp/scaling.h"

#include <algorithm>
#include <cmath>

namespace pramana {

namespace {
double pow2Round(double v) {
  if (!(v > 0) || !std::isfinite(v)) return 1.0;
  int e = static_cast<int>(std::lround(std::log2(v)));
  e = std::max(-60, std::min(60, e));
  return std::ldexp(1.0, e);
}
}  // namespace

double coefficientRange(const SparseMatrix& A) {
  double mx = 0, mn = 1e300;
  for (double v : A.value) {
    double a = std::fabs(v);
    if (a == 0) continue;
    mx = std::max(mx, a);
    mn = std::min(mn, a);
  }
  return mx > 0 ? mx / mn : 1.0;
}

void applyScaling(SparseMatrix& A, const ScalingResult& s) {
  for (int j = 0; j < A.numCols; ++j)
    for (int k = A.start[j]; k < A.start[j + 1]; ++k)
      A.value[k] *= s.rowScale[A.index[k]] * s.colScale[j];
}

ScalingResult geometricScaling(const SparseMatrix& A, int passes) {
  const int m = A.numRows, n = A.numCols;
  ScalingResult res;
  res.rowScale.assign(m, 1.0);
  res.colScale.assign(n, 1.0);
  res.rangeBefore = coefficientRange(A);
  if (A.nnz() == 0) return res;
  std::vector<double> rmin(m), rmax(m);
  double prevRange = res.rangeBefore;
  for (int pass = 0; pass < passes; ++pass) {
    // Row pass.
    std::fill(rmin.begin(), rmin.end(), 1e300);
    std::fill(rmax.begin(), rmax.end(), 0.0);
    for (int j = 0; j < n; ++j)
      for (int k = A.start[j]; k < A.start[j + 1]; ++k) {
        double a = std::fabs(A.value[k]) * res.colScale[j] * res.rowScale[A.index[k]];
        if (a == 0) continue;
        rmin[A.index[k]] = std::min(rmin[A.index[k]], a);
        rmax[A.index[k]] = std::max(rmax[A.index[k]], a);
      }
    for (int i = 0; i < m; ++i)
      if (rmax[i] > 0) res.rowScale[i] /= std::sqrt(rmin[i] * rmax[i]);
    // Column pass.
    double gmax = 0, gmin = 1e300;
    for (int j = 0; j < n; ++j) {
      double cmin = 1e300, cmax = 0;
      for (int k = A.start[j]; k < A.start[j + 1]; ++k) {
        double a = std::fabs(A.value[k]) * res.colScale[j] * res.rowScale[A.index[k]];
        if (a == 0) continue;
        cmin = std::min(cmin, a);
        cmax = std::max(cmax, a);
      }
      if (cmax > 0) {
        double f = 1.0 / std::sqrt(cmin * cmax);
        res.colScale[j] *= f;
        gmax = std::max(gmax, cmax * f);
        gmin = std::min(gmin, cmin * f);
      }
    }
    double range = gmax > 0 ? gmax / gmin : 1.0;
    if (range > 0.9 * prevRange && pass >= 1) break;
    prevRange = range;
  }
  // Column equilibration: max |a_ij| in each column becomes ~1.
  for (int j = 0; j < n; ++j) {
    double cmax = 0;
    for (int k = A.start[j]; k < A.start[j + 1]; ++k)
      cmax = std::max(cmax, std::fabs(A.value[k]) * res.colScale[j] * res.rowScale[A.index[k]]);
    if (cmax > 0) res.colScale[j] /= cmax;
  }
  for (auto& v : res.rowScale) v = pow2Round(v);
  for (auto& v : res.colScale) v = pow2Round(v);
  SparseMatrix tmp = A;
  applyScaling(tmp, res);
  res.rangeAfter = coefficientRange(tmp);
  if (res.rangeAfter > res.rangeBefore) {  // never make things worse
    std::fill(res.rowScale.begin(), res.rowScale.end(), 1.0);
    std::fill(res.colScale.begin(), res.colScale.end(), 1.0);
    res.rangeAfter = res.rangeBefore;
  }
  return res;
}

ScalingResult ruizScaling(const SparseMatrix& A, int iterations, bool pockChambolle) {
  const int m = A.numRows, n = A.numCols;
  ScalingResult res;
  res.rowScale.assign(m, 1.0);
  res.colScale.assign(n, 1.0);
  res.rangeBefore = coefficientRange(A);
  std::vector<double> rmax(m), cmax(n);
  for (int it = 0; it < iterations; ++it) {
    std::fill(rmax.begin(), rmax.end(), 0.0);
    for (int j = 0; j < n; ++j) {
      double c = 0;
      for (int k = A.start[j]; k < A.start[j + 1]; ++k) {
        double a = std::fabs(A.value[k]) * res.colScale[j] * res.rowScale[A.index[k]];
        c = std::max(c, a);
        rmax[A.index[k]] = std::max(rmax[A.index[k]], a);
      }
      cmax[j] = c;
    }
    for (int i = 0; i < m; ++i)
      if (rmax[i] > 0) res.rowScale[i] /= std::sqrt(rmax[i]);
    for (int j = 0; j < n; ++j)
      if (cmax[j] > 0) res.colScale[j] /= std::sqrt(cmax[j]);
  }
  if (pockChambolle) {  // alpha = 1: row/col 1-norms
    std::vector<double> rsum(m, 0.0), csum(n, 0.0);
    for (int j = 0; j < n; ++j)
      for (int k = A.start[j]; k < A.start[j + 1]; ++k) {
        double a = std::fabs(A.value[k]) * res.colScale[j] * res.rowScale[A.index[k]];
        csum[j] += a;
        rsum[A.index[k]] += a;
      }
    for (int i = 0; i < m; ++i)
      if (rsum[i] > 0) res.rowScale[i] /= std::sqrt(rsum[i]);
    for (int j = 0; j < n; ++j)
      if (csum[j] > 0) res.colScale[j] /= std::sqrt(csum[j]);
  }
  for (auto& v : res.rowScale) v = pow2Round(v);
  for (auto& v : res.colScale) v = pow2Round(v);
  SparseMatrix tmp = A;
  applyScaling(tmp, res);
  res.rangeAfter = coefficientRange(tmp);
  return res;
}

}  // namespace pramana
