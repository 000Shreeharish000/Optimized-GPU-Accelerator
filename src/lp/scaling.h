#pragma once

#include <vector>

#include "core/sparse.h"

namespace pramana {

struct ScalingResult {
  std::vector<double> rowScale;  // A' = diag(rowScale) * A * diag(colScale)
  std::vector<double> colScale;
  double rangeBefore = 0;  // max|a| / min|a| (nonzeros)
  double rangeAfter = 0;
};

// Geometric-mean scaling (alternating row/column passes) followed by column
// equilibration. All factors are rounded to powers of two so that scaling and
// unscaling introduce *no* rounding error in the matrix, bounds or costs.
ScalingResult geometricScaling(const SparseMatrix& A, int passes = 6);

// Ruiz equilibration (infinity norm) + optional Pock-Chambolle (alpha=1) step,
// as used by PDLP-type first-order methods. Factors are powers of two.
ScalingResult ruizScaling(const SparseMatrix& A, int iterations = 10, bool pockChambolle = true);

// Applies scaling in place: A <- R A C.
void applyScaling(SparseMatrix& A, const ScalingResult& s);
double coefficientRange(const SparseMatrix& A);

}  // namespace pramana
