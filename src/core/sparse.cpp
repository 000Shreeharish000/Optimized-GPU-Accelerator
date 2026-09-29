#include "core/sparse.h"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "util/common.h"

namespace pramana {

SparseMatrix SparseMatrix::fromTriplets(int rows, int cols, const std::vector<int>& ri,
                                        const std::vector<int>& ci, const std::vector<double>& v) {
  SparseMatrix m(rows, cols);
  std::vector<int> count(cols + 1, 0);
  for (int c : ci) count[c + 1]++;
  for (int j = 0; j < cols; ++j) count[j + 1] += count[j];
  std::vector<int> idx(v.size());
  std::vector<double> val(v.size());
  std::vector<int> pos(count.begin(), count.end() - 1);
  for (size_t k = 0; k < v.size(); ++k) {
    int p = pos[ci[k]]++;
    idx[p] = ri[k];
    val[p] = v[k];
  }
  // Sort each column by row index and merge duplicates.
  std::vector<int> order;
  m.start[0] = 0;
  for (int j = 0; j < cols; ++j) {
    int b = count[j], e = count[j + 1];
    order.resize(e - b);
    std::iota(order.begin(), order.end(), b);
    std::sort(order.begin(), order.end(), [&](int a, int c) { return idx[a] < idx[c]; });
    int lastRow = -1;
    for (int k : order) {
      if (idx[k] == lastRow) {
        m.value.back() += val[k];
      } else {
        m.index.push_back(idx[k]);
        m.value.push_back(val[k]);
        lastRow = idx[k];
      }
    }
    // Drop explicit zeros created by merging (or present in input).
    int w = m.start[j];
    for (int k = m.start[j]; k < static_cast<int>(m.index.size()); ++k) {
      if (m.value[k] != 0.0) {
        m.index[w] = m.index[k];
        m.value[w] = m.value[k];
        ++w;
      }
    }
    m.index.resize(w);
    m.value.resize(w);
    m.start[j + 1] = w;
  }
  return m;
}

SparseMatrix SparseMatrix::transpose() const {
  SparseMatrix t(numCols, numRows);
  std::vector<int> count(numRows + 1, 0);
  for (int k = 0; k < nnz(); ++k) count[index[k] + 1]++;
  for (int i = 0; i < numRows; ++i) count[i + 1] += count[i];
  t.start = count;
  t.index.resize(nnz());
  t.value.resize(nnz());
  std::vector<int> pos(count.begin(), count.end() - 1);
  for (int j = 0; j < numCols; ++j)
    for (int k = start[j]; k < start[j + 1]; ++k) {
      int p = pos[index[k]]++;
      t.index[p] = j;
      t.value[p] = value[k];
    }
  return t;
}

void SparseMatrix::multiply(const double* x, double* y, double alpha) const {
  for (int j = 0; j < numCols; ++j) {
    double xj = x[j];
    if (xj == 0.0) continue;
    xj *= alpha;
    for (int k = start[j]; k < start[j + 1]; ++k) y[index[k]] += value[k] * xj;
  }
}

void SparseMatrix::multiplyTranspose(const double* x, double* y, double alpha) const {
  for (int j = 0; j < numCols; ++j) {
    double s = 0;
    for (int k = start[j]; k < start[j + 1]; ++k) s += value[k] * x[index[k]];
    y[j] += alpha * s;
  }
}

double SparseMatrix::maxAbs() const {
  double m = 0;
  for (double v : value) m = std::max(m, std::fabs(v));
  return m;
}

void SparseMatrix::appendColumn(const std::vector<int>& idx, const std::vector<double>& val) {
  for (size_t k = 0; k < idx.size(); ++k) {
    if (val[k] == 0.0) continue;
    index.push_back(idx[k]);
    value.push_back(val[k]);
  }
  start.push_back(static_cast<int>(index.size()));
  ++numCols;
}

SparseMatrix SparseMatrix::selectColumns(const std::vector<int>& cols) const {
  SparseMatrix m(numRows, 0);
  for (int j : cols) {
    for (int k = start[j]; k < start[j + 1]; ++k) {
      m.index.push_back(index[k]);
      m.value.push_back(value[k]);
    }
    m.start.push_back(static_cast<int>(m.index.size()));
    ++m.numCols;
  }
  return m;
}

void SparseMatrix::sortIndices() {
  std::vector<std::pair<int, double>> tmp;
  for (int j = 0; j < numCols; ++j) {
    tmp.clear();
    for (int k = start[j]; k < start[j + 1]; ++k) tmp.emplace_back(index[k], value[k]);
    std::sort(tmp.begin(), tmp.end());
    for (int k = start[j], t = 0; k < start[j + 1]; ++k, ++t) {
      index[k] = tmp[t].first;
      value[k] = tmp[t].second;
    }
  }
}

void HVector::clear() {
  if (count >= 0 && count < size / 4) {
    for (int k = 0; k < count; ++k) array[index[k]] = 0.0;
  } else {
    std::fill(array.begin(), array.end(), 0.0);
  }
  count = 0;
}

void HVector::rebuildIndex(double dropTol) {
  count = 0;
  for (int i = 0; i < size; ++i) {
    if (array[i] != 0.0) {
      if (std::fabs(array[i]) <= dropTol) {
        array[i] = 0.0;
      } else {
        index[count++] = i;
      }
    }
  }
}

double HVector::norm2sq() const {
  double s = 0;
  for (int k = 0; k < count; ++k) s += array[index[k]] * array[index[k]];
  return s;
}

}  // namespace pramana
