// Sparse matrix storage. The master copy of A is column-major (CSC) because the
// simplex method needs column access (FTRAN of an entering column); a row-wise
// copy (CSR == CSC of the transpose) is built on demand for PRICE, activity
// computations and the GPU kernels.
#pragma once

#include <vector>

namespace pramana {

struct SparseMatrix {
  int numRows = 0;
  int numCols = 0;
  std::vector<int> start;    // size numCols+1 (column pointers)
  std::vector<int> index;    // row indices
  std::vector<double> value;

  SparseMatrix() { start.assign(1, 0); }
  SparseMatrix(int rows, int cols) : numRows(rows), numCols(cols) { start.assign(cols + 1, 0); }

  int nnz() const { return start.empty() ? 0 : start.back(); }
  int colBegin(int j) const { return start[j]; }
  int colEnd(int j) const { return start[j + 1]; }

  // Build from (row, col, value) triplets; duplicates are summed, zeros dropped.
  static SparseMatrix fromTriplets(int rows, int cols, const std::vector<int>& ri,
                                   const std::vector<int>& ci, const std::vector<double>& v);

  SparseMatrix transpose() const;
  // y += alpha * A x
  void multiply(const double* x, double* y, double alpha = 1.0) const;
  // y += alpha * A^T x
  void multiplyTranspose(const double* x, double* y, double alpha = 1.0) const;
  double maxAbs() const;
  // Append a column (indices must be < numRows).
  void appendColumn(const std::vector<int>& idx, const std::vector<double>& val);
  // Keep only the listed columns / rows (in the given order for columns).
  SparseMatrix selectColumns(const std::vector<int>& cols) const;
  void sortIndices();
};

// Dense-index sparse vector with explicit nonzero pattern (Hall & McKinnon
// "hyper-sparse" vector): values stored densely, nonzero positions listed.
struct HVector {
  int size = 0;
  int count = 0;                // number of listed nonzeros (valid if !dense)
  std::vector<int> index;       // nonzero positions
  std::vector<double> array;    // dense values
  bool packedValid = false;

  void setup(int n) {
    size = n;
    count = 0;
    index.assign(n, 0);
    array.assign(n, 0.0);
  }
  void clear();
  void rebuildIndex(double dropTol = 0.0);  // recompute index from array
  double norm2sq() const;
};

}  // namespace pramana
