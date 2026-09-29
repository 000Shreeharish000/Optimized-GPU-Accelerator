// Sparse LU factorization of a simplex basis matrix B, written from scratch.
//
//   * Basis columns come from [A | I]: variable v < n is column v of A, variable
//     v >= n is the unit column e_(v-n) (row logical).
//   * Pass 1 extracts column singletons and row singletons (the triangular part
//     of B, typically most of it: no fill, no numerical choice).
//   * Pass 2 factors the remaining "bump" with Markowitz pivoting under a
//     threshold test |a_ij| >= u * max_k |a_kj|  (u = 0.1 by default, raised on
//     instability), using count buckets to find low-Markowitz-cost pivots.
//   * Rank deficiency is reported (positions + rows) so the caller can repair
//     the basis by substituting row logicals.
//   * Updates use the product form (PFI eta file); the caller refactors after
//     `maxUpdates` or when a residual check fails.
//
// Factor storage (pivot order k = 0..K-1): pivot row r_k, basic position c_k,
// pivot value; L_k = multipliers (row i, l_ik); U row-wise (positions pivoted
// later) and U column-wise (rows pivoted earlier) for sparse FTRAN.
#pragma once

#include <vector>

#include "core/sparse.h"

namespace pramana {

struct LuStats {
  int factorizations = 0;
  int rankDeficiencies = 0;
  long long lNnz = 0, uNnz = 0;
  int bumpSize = 0;
  double lastGrowth = 0;
};

class LuFactor {
 public:
  void setup(int numRows, int numStructural, const SparseMatrix* A);

  // Returns the number of basis columns that could not be pivoted (0 = full rank).
  int factorize(const std::vector<int>& basicIndex);

  // Solves B x = b. In: `rhs` indexed by row. Out: indexed by basic position.
  void ftran(HVector& rhs) const;
  // Solves B^T y = e. In: indexed by basic position. Out: indexed by row.
  void btran(HVector& rhs) const;
  // Records a basis change at `pos`; `column` must be B^{-1} a_q (FTRAN result).
  void update(const HVector& column, int pos);
  int numUpdates() const { return static_cast<int>(etaPos_.size()); }

  double pivotThreshold = 0.1;
  double absPivotTol = 1e-11;
  std::vector<int> deficientPositions;  // basic positions not pivoted
  std::vector<int> deficientRows;       // rows not pivoted
  LuStats stats;

 private:
  void clearFactors();
  void buildUColumns();

  int m_ = 0, n_ = 0;
  const SparseMatrix* A_ = nullptr;

  std::vector<int> pivotRow_, pivotPos_;
  std::vector<double> pivotValue_;
  std::vector<int> lStart_, lIndex_;
  std::vector<double> lValue_;
  std::vector<int> uStart_, uIndex_;  // row-wise: positions
  std::vector<double> uValue_;
  std::vector<int> ucStart_, ucIndex_;  // column-wise: steps (row = pivotRow_[step])
  std::vector<double> ucValue_;
  std::vector<int> posStep_;  // basic position -> step

  // PFI etas
  std::vector<int> etaPos_, etaStart_, etaIndex_;
  std::vector<double> etaPivot_, etaValue_;

  // work arrays
  mutable std::vector<double> work_;
  mutable std::vector<char> mark_;
};

}  // namespace pramana
