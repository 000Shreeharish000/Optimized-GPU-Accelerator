// Sparse symmetric LDL^T factorization for quasidefinite KKT systems, written
// from scratch:
//
//   * fill-reducing ordering: approximate minimum degree on the quotient graph
//     (elements absorb eliminated neighbourhoods; degrees are the usual
//     external-degree upper bound, as in AMD without supervariables);
//   * symbolic analysis: elimination tree + exact column counts (up-looking);
//   * numeric factorization: up-looking row-by-row sparse triangular solves;
//   * dynamic regularization: each pivot must carry its expected sign
//     (quasidefinite matrix [-H A^T; A dI]: first block negative, second
//     positive). Tiny or wrong-sign pivots are replaced by +-reg and counted,
//     which is the standard static-pivoting remedy (Vanderbei's quasidefinite
//     theory guarantees existence of the factorization for any ordering);
//   * iterative refinement is done by the caller against the unregularized system.
#pragma once

#include <vector>

#include "core/sparse.h"

namespace pramana {

// Minimum-degree ordering of the pattern of a symmetric matrix (both triangles
// or either one may be given). Returns perm: perm[k] = original index eliminated k-th.
std::vector<int> minimumDegreeOrdering(const SparseMatrix& symPattern);

class LdlFactor {
 public:
  // K: full symmetric matrix (both triangles) in CSC. signs[i] = -1 / +1 expected pivot sign.
  void analyze(const SparseMatrix& K, const std::vector<int>& perm);
  // Returns number of regularized pivots. regularization <= 0 selects inertia mode
  // (pivots kept as computed; used by the convexity check).
  int factorize(const SparseMatrix& K, const std::vector<int8_t>& signs, double regularization);
  void solve(std::vector<double>& b) const;  // in place, original ordering
  long long nnzL() const { return static_cast<long long>(li_.size()); }
  int dimension() const { return n_; }
  std::vector<double> pivots() const { return d_; }  // permuted order
  const std::vector<int>& permutation() const { return perm_; }
  int negativePivots() const;
  // out = P^T L^{-T} e_k (original ordering): out^T K out = D_k exactly.
  void solveUnitLt(int k, std::vector<double>& out) const;

 private:
  int n_ = 0;
  std::vector<int> perm_, pinv_;
  std::vector<int> parent_, lp_, lnz_;
  std::vector<int> li_;
  std::vector<double> lx_, d_;
  mutable std::vector<double> work_;
};

}  // namespace pramana
