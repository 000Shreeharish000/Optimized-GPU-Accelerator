// Activity-based bound propagation (domain propagation) for MIP nodes.
//
// For a row  L <= sum_j a_j x_j <= U  with minimum / maximum activities
// minAct / maxAct over the current box, each variable satisfies
//     a_j > 0:  x_j <= l_j + (U - minAct)/a_j,   x_j >= u_j + (L - maxAct)/a_j
//     a_j < 0:  x_j >= u_j + (U - minAct)/a_j,   x_j <= l_j + (L - maxAct)/a_j
// Integer columns are rounded (with tolerance); continuous tightenings are only
// applied when significant (and relaxed by a safety margin) so floating-point
// noise can never cut off a feasible point. Infinite activity contributions are
// tracked by count so a single unbounded variable can still be tightened.
#pragma once

#include <vector>

#include "core/model.h"

namespace pramana {

struct PropagationResult {
  bool infeasible = false;
  int tightened = 0;
  int rounds = 0;
};

class Propagator {
 public:
  void setup(const Model& model);  // model: rows/cols/integrality (min form)
  // Tightens lower/upper in place. Rows given by `extraRows` (cuts) are ignored.
  PropagationResult propagate(std::vector<double>& lower, std::vector<double>& upper, int maxRounds = 5) const;
  int numRows() const { return m_; }

 private:
  int m_ = 0, n_ = 0;
  SparseMatrix AT_;  // row-wise
  SparseMatrix A_;   // column-wise
  std::vector<double> rl_, ru_;
  std::vector<char> isInt_;
};

}  // namespace pramana
