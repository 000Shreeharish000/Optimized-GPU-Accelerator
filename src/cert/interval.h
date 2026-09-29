// Rigorous outward-rounded interval arithmetic that does not depend on the FPU
// rounding mode. Each round-to-nearest result is tested for exactness with an
// error-free transformation (TwoSum for +/-, FMA-based TwoProduct for *); only
// inexact results are widened by one ulp in the safe direction. This is exact
// as long as the compiler does not contract a*b+c into FMA on its own (MSVC
// /fp:precise default, GCC/Clang -ffp-contract=off; see CMakeLists.txt).
#pragma once

#include <cmath>
#include <limits>

namespace pramana {

inline double roundDown(double v) {
  if (std::isinf(v) || std::isnan(v)) return v;
  return std::nextafter(v, -std::numeric_limits<double>::infinity());
}
inline double roundUp(double v) {
  if (std::isinf(v) || std::isnan(v)) return v;
  return std::nextafter(v, std::numeric_limits<double>::infinity());
}

// Error of s = fl(a+b):  a + b = s + err exactly (Knuth's TwoSum).
inline double twoSumError(double a, double b, double s) {
  double bb = s - a;
  return (a - (s - bb)) + (b - bb);
}

// Lower / upper bounds of the exact a + b.
inline double addDown(double a, double b) {
  double s = a + b;
  if (std::isinf(s) || std::isnan(s)) return s;
  double e = twoSumError(a, b, s);
  return e < 0 ? roundDown(s) : s;
}
inline double addUp(double a, double b) {
  double s = a + b;
  if (std::isinf(s) || std::isnan(s)) return s;
  double e = twoSumError(a, b, s);
  return e > 0 ? roundUp(s) : s;
}
inline double mulDown(double a, double b) {
  if (a == 0 || b == 0) return 0.0;
  double p = a * b;
  if (std::isinf(p)) return p;
  double e = std::fma(a, b, -p);  // exact remainder a*b - p
  return e < 0 ? roundDown(p) : p;
}
inline double mulUp(double a, double b) {
  if (a == 0 || b == 0) return 0.0;
  double p = a * b;
  if (std::isinf(p)) return p;
  double e = std::fma(a, b, -p);
  return e > 0 ? roundUp(p) : p;
}

struct Interval {
  double lo = 0, hi = 0;
  static Interval point(double v) { return {v, v}; }
  bool containsZero() const { return lo <= 0 && hi >= 0; }
};

// [a,a] * [b,b] for exact doubles a, b.
inline Interval mulExact(double a, double b) { return {mulDown(a, b), mulUp(a, b)}; }
inline Interval add(Interval x, Interval y) { return {addDown(x.lo, y.lo), addUp(x.hi, y.hi)}; }
inline Interval sub(Interval x, Interval y) { return {addDown(x.lo, -y.hi), addUp(x.hi, -y.lo)}; }

// Rigorous lower bound of min{ z*x : z in [z.lo,z.hi], x in [l,u] } (bounds may be infinite).
inline double minProduct(Interval z, double l, double u) {
  const double inf = std::numeric_limits<double>::infinity();
  if (z.lo == 0 && z.hi == 0) return 0.0;
  if (std::isinf(l) && z.hi > 0) return -inf;  // x -> -inf with z > 0
  if (std::isinf(u) && z.lo < 0) return -inf;  // x -> +inf with z < 0
  double best = inf;
  const double zs[2] = {z.lo, z.hi};
  const double xs[2] = {l, u};
  for (double zz : zs)
    for (double xx : xs) {
      double p;
      if (zz == 0 || xx == 0) p = 0.0;
      else if (std::isinf(xx)) p = zz * xx;  // sign is safe here (checked above)
      else p = mulDown(zz, xx);
      if (p < best) best = p;
    }
  return best;
}

inline double maxProduct(Interval z, double l, double u) {
  Interval neg{-z.hi, -z.lo};
  return -minProduct(neg, l, u);
}

// Neumaier/Kahan–Babuska compensated summation (for accurate residuals).
struct CompensatedSum {
  double sum = 0, comp = 0;
  void add(double v) {
    double t = sum + v;
    if (std::fabs(sum) >= std::fabs(v)) comp += (sum - t) + v;
    else comp += (v - t) + sum;
    sum = t;
  }
  double value() const { return sum + comp; }
};

}  // namespace pramana
