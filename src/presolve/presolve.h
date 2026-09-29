// Presolve with a full postsolve stack (primal values, row duals, reduced
// costs AND basis statuses), written from scratch.
//
// Reductions (each records exactly what postsolve needs):
//   empty rows / empty columns, fixed columns, singleton rows -> bounds,
//   redundant rows (activity bounds inside [L,U]), forcing rows (activity
//   bound meets the rhs: all columns fixed, dual recovered by a ratio test),
//   dominated columns (dual fixing), free/implied-free column singletons in
//   equality rows (column substituted out, dual y_i = c_j / a_ij),
//   MIP: integer bound rounding.
//
// Infeasibility / unboundedness detected in presolve is NOT reported directly:
// the orchestrator re-solves the original model so the claim carries a checkable
// certificate (Farkas ray / primal ray) from the simplex.
#pragma once

#include <vector>

#include "core/model.h"

namespace pramana {

struct PresolveStats {
  int rowsRemoved = 0, colsRemoved = 0, nnzRemoved = 0;
  int emptyRows = 0, emptyCols = 0, fixedCols = 0, singletonRows = 0, redundantRows = 0;
  int forcingRows = 0, dominatedCols = 0, freeColSingletons = 0, boundsTightened = 0;
  int passes = 0;
  double seconds = 0;
};

enum class PresolveStatus { Reduced, Unchanged, Infeasible, Unbounded, Empty };

class Presolver {
 public:
  // `model` must be in minimization form.
  PresolveStatus run(const Model& model, bool mip);
  const Model& reduced() const { return reduced_; }
  const PresolveStats& stats() const { return stats_; }

  // Map a reduced-space solution back to the original model.
  struct Solution {
    std::vector<double> x, rowDual, reducedCost;
    std::vector<BasisStatus> colStatus, rowStatus;
    bool hasDual = false, hasBasis = false;
  };
  Solution postsolve(const Solution& reducedSol) const;
  // Reduced -> original index maps.
  const std::vector<int>& colMap() const { return colOrig_; }
  const std::vector<int>& rowMap() const { return rowOrig_; }
  // A row whose activity bounds strictly miss its range (candidate Farkas ray +-e_i), or -1.
  int certificateRow() const { return certRow_; }
  int certificateSign() const { return certSign_; }

 private:
  enum class Kind { FixedCol, EmptyRow, SingletonRow, RedundantRow, ForcingRow, FreeColSingleton };
  struct Record {
    Kind kind;
    int row = -1, col = -1;
    double value = 0, a = 0;
    double oldLower = 0, oldUpper = 0, rowLower = 0, rowUpper = 0;
    int side = 0;  // singleton/forcing: which side became active (-1 lower, +1 upper)
    std::vector<int> cols;           // forcing row: fixed columns; free singleton: other cols in row
    std::vector<double> vals, colVals;
    double cost = 0;
  };

  double activityMin(int i) const;
  double activityMax(int i) const;
  void removeRow(int i);
  void fixColumn(int j, double v, bool record = true);
  bool pass(bool mip);

  Model orig_, reduced_;
  PresolveStats stats_;
  SparseMatrix AT_;  // row-wise copy of original A
  std::vector<double> lo_, up_, rl_, ru_, cost_;
  double offset_ = 0;
  std::vector<char> rowActive_, colActive_;
  std::vector<int> rowCount_, colCount_;
  std::vector<Record> stack_;
  std::vector<std::vector<double>> forcingCosts_;
  std::vector<int> colOrig_, rowOrig_;
  bool infeasible_ = false, unbounded_ = false;
  bool qp_ = false;
  int certRow_ = -1, certSign_ = 0;  // quadratic objective: only row reductions (singleton / empty / redundant)
};

}  // namespace pramana
