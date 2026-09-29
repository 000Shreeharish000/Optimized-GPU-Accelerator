// Cutting planes for the root node, written from scratch:
//   * Gomory mixed-integer (GMI) cuts from simplex tableau rows of fractional
//     basic integer variables (nonbasic columns shifted/complemented to x' >= 0);
//   * complemented mixed-integer rounding (c-MIR) on single rows, trying
//     several divisors delta (Marchand-Wolsey);
//   * lifted-by-extension knapsack cover cuts on pure-binary rows.
// Every cut goes through a numerical safety filter: bounded dynamism, removal of
// tiny coefficients by *conservative* rhs relaxation using variable bounds, and a
// small rhs relaxation so round-off cannot cut off feasible integer points.
#pragma once

#include <string>
#include <vector>

#include "core/model.h"

namespace pramana {

class Simplex;

// A cut in the form  sum_k val[k] * x[idx[k]] >= rhs.
struct Cut {
  std::vector<int> idx;
  std::vector<double> val;
  double rhs = 0;
  std::string type;
  double efficacy = 0;  // violation / ||val||_2 at the separated point
};

struct CutStats {
  int gomory = 0, mir = 0, cover = 0, rejectedNumerics = 0;
};

class CutGenerator {
 public:
  explicit CutGenerator(const Model& model);
  std::vector<Cut> gomory(const Simplex& lp, const std::vector<double>& x, const std::vector<double>& lo,
                          const std::vector<double>& up, int maxCuts);
  std::vector<Cut> mir(const std::vector<double>& x, const std::vector<double>& lo, const std::vector<double>& up,
                       int maxCuts);
  std::vector<Cut> covers(const std::vector<double>& x, const std::vector<double>& lo, const std::vector<double>& up,
                          int maxCuts);
  // Efficacy threshold + parallelism filter (cosine < maxParallel), best first.
  static std::vector<Cut> select(std::vector<Cut> cuts, int maxCuts, double minEfficacy, double maxParallel);
  CutStats stats;

 private:
  // Finalizes a dense cut (coefficients over structural columns): safety filter,
  // efficacy at x. Returns false if rejected.
  bool finalize(std::vector<double>& dense, double rhs, const std::vector<double>& x, const std::vector<double>& lo,
                const std::vector<double>& up, const std::string& type, Cut& out);
  const Model& model_;
  SparseMatrix AT_;
  std::vector<char> isInt_;
};

}  // namespace pramana
