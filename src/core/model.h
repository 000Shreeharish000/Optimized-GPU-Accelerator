// Internal model form used by every engine:
//
//     min  c^T x + 1/2 x^T Q x + offset
//     s.t. rowLower <= A x <= rowUpper
//          colLower <= x   <= colUpper
//          x_j integer for colType[j] == Integer
//
// Row activities are *not* converted to standard form; bounds are kept on both
// columns and rows (every production simplex code does this, and converting to
// Ax=b, x>=0 would double the variables and destroy the bound-flipping ratio test).
#pragma once

#include <string>
#include <vector>

#include "core/sparse.h"
#include "util/common.h"

namespace pramana {

enum class VarType : uint8_t { Continuous = 0, Integer = 1 };
enum class ObjSense : int8_t { Minimize = 1, Maximize = -1 };

struct Model {
  std::string name = "model";
  std::string objName = "OBJ";
  ObjSense sense = ObjSense::Minimize;
  double objOffset = 0.0;

  std::vector<double> colCost, colLower, colUpper;
  std::vector<double> rowLower, rowUpper;
  std::vector<VarType> colType;
  std::vector<std::string> colNames, rowNames;
  SparseMatrix A;  // numRows x numCols, CSC
  SparseMatrix Q;  // numCols x numCols, CSC, full symmetric (both triangles). Empty => LP.

  int numCols() const { return static_cast<int>(colCost.size()); }
  int numRows() const { return static_cast<int>(rowLower.size()); }
  int nnz() const { return A.nnz(); }
  bool isMip() const;
  bool isQp() const { return Q.nnz() > 0; }
  int numIntegers() const;

  // Adds a column; returns its index.
  int addColumn(double cost, double lower, double upper, const std::vector<int>& rows,
                const std::vector<double>& vals, VarType type = VarType::Continuous,
                const std::string& name = "");
  // Adds a row over existing columns; returns its index. (Rebuilds CSC: O(nnz).)
  int addRow(double lower, double upper, const std::vector<int>& cols,
             const std::vector<double>& vals, const std::string& name = "");
  void addRows(const std::vector<double>& lower, const std::vector<double>& upper,
               const std::vector<int>& starts, const std::vector<int>& cols,
               const std::vector<double>& vals);

  // Objective in the model's own sense (includes offset).
  double objective(const std::vector<double>& x) const;
  std::vector<double> rowActivity(const std::vector<double>& x) const;
  // Converts to internal minimization (negates c, Q, offset if maximizing).
  Model toMinimization() const;

  void validate() const;  // throws on inconsistent dimensions / NaN / lower>upper
  std::string summary() const;
  void ensureNames();
};

}  // namespace pramana
