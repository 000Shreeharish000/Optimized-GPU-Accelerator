#include "core/model.h"

#include <algorithm>
#include <cmath>

#include "util/log.h"

namespace pramana {

bool Model::isMip() const {
  for (auto t : colType)
    if (t == VarType::Integer) return true;
  return false;
}

int Model::numIntegers() const {
  int k = 0;
  for (auto t : colType) k += t == VarType::Integer;
  return k;
}

int Model::addColumn(double cost, double lower, double upper, const std::vector<int>& rows,
                     const std::vector<double>& vals, VarType type, const std::string& nm) {
  A.numRows = numRows();
  colCost.push_back(cost);
  colLower.push_back(lower);
  colUpper.push_back(upper);
  colType.push_back(type);
  colNames.push_back(nm.empty() ? "C" + std::to_string(colCost.size() - 1) : nm);
  A.appendColumn(rows, vals);
  if (Q.numCols > 0) {
    Q.numRows++;
    Q.start.push_back(Q.start.back());
    Q.numCols++;
  }
  return numCols() - 1;
}

void Model::addRows(const std::vector<double>& lower, const std::vector<double>& upper,
                    const std::vector<int>& starts, const std::vector<int>& cols,
                    const std::vector<double>& vals) {
  int first = numRows();
  int added = static_cast<int>(lower.size());
  std::vector<int> ri, ci;
  std::vector<double> v;
  ri.reserve(A.nnz() + cols.size());
  for (int j = 0; j < A.numCols; ++j)
    for (int k = A.start[j]; k < A.start[j + 1]; ++k) {
      ri.push_back(A.index[k]);
      ci.push_back(j);
      v.push_back(A.value[k]);
    }
  for (int r = 0; r < added; ++r) {
    for (int k = starts[r]; k < starts[r + 1]; ++k) {
      ri.push_back(first + r);
      ci.push_back(cols[k]);
      v.push_back(vals[k]);
    }
    rowLower.push_back(lower[r]);
    rowUpper.push_back(upper[r]);
    rowNames.push_back("R" + std::to_string(first + r));
  }
  A = SparseMatrix::fromTriplets(numRows(), numCols(), ri, ci, v);
}

int Model::addRow(double lower, double upper, const std::vector<int>& cols,
                  const std::vector<double>& vals, const std::string& nm) {
  std::vector<int> starts{0, static_cast<int>(cols.size())};
  addRows({lower}, {upper}, starts, cols, vals);
  if (!nm.empty()) rowNames.back() = nm;
  return numRows() - 1;
}

double Model::objective(const std::vector<double>& x) const {
  double f = objOffset;
  for (int j = 0; j < numCols(); ++j) f += colCost[j] * x[j];
  if (isQp()) {
    double q = 0;
    for (int j = 0; j < Q.numCols; ++j)
      for (int k = Q.start[j]; k < Q.start[j + 1]; ++k) q += x[Q.index[k]] * Q.value[k] * x[j];
    f += 0.5 * q;
  }
  return f;
}

std::vector<double> Model::rowActivity(const std::vector<double>& x) const {
  std::vector<double> r(numRows(), 0.0);
  A.multiply(x.data(), r.data());
  return r;
}

Model Model::toMinimization() const {
  Model m = *this;
  if (sense == ObjSense::Maximize) {
    for (double& c : m.colCost) c = -c;
    for (double& q : m.Q.value) q = -q;
    m.objOffset = -m.objOffset;
    m.sense = ObjSense::Minimize;
  }
  return m;
}

void Model::validate() const {
  int n = numCols(), m = numRows();
  PRAMANA_CHECK(static_cast<int>(colLower.size()) == n && static_cast<int>(colUpper.size()) == n,
                "column bound size");
  PRAMANA_CHECK(static_cast<int>(colType.size()) == n, "column type size");
  PRAMANA_CHECK(static_cast<int>(rowUpper.size()) == m, "row bound size");
  PRAMANA_CHECK(A.numCols == n && (A.numRows == m || (n == 0)), "matrix dimensions");
  PRAMANA_CHECK(Q.nnz() == 0 || (Q.numCols == n && Q.numRows == n), "Q dimensions");
  for (int j = 0; j < n; ++j) {
    PRAMANA_CHECK(!std::isnan(colCost[j]) && !std::isnan(colLower[j]) && !std::isnan(colUpper[j]),
                  "NaN in column data");
  }
  for (double v : A.value) PRAMANA_CHECK(std::isfinite(v), "non-finite matrix coefficient");
}

std::string Model::summary() const {
  return formatString("%s: %d rows, %d cols (%d integer), %d nnz%s", name.c_str(), numRows(),
                      numCols(), numIntegers(), nnz(),
                      isQp() ? formatString(", Q nnz %d", Q.nnz()).c_str() : "");
}

void Model::ensureNames() {
  colNames.resize(numCols());
  rowNames.resize(numRows());
  for (int j = 0; j < numCols(); ++j)
    if (colNames[j].empty()) colNames[j] = "C" + std::to_string(j);
  for (int i = 0; i < numRows(); ++i)
    if (rowNames[i].empty()) rowNames[i] = "R" + std::to_string(i);
}

}  // namespace pramana
