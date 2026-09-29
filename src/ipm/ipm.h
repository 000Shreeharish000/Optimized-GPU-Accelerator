// Primal-dual interior point method (Mehrotra predictor-corrector) for LP and
// convex QP, written from scratch on top of the own sparse LDL^T (ipm/ldl.h).
//
//   min c^T x + 1/2 x^T Q x   s.t.  L <= Ax <= U,  l <= x <= u
//
// Inequality rows get slacks (a_i^T x - w_i = 0, L_i <= w_i <= U_i), fixed
// columns are eliminated, bounds are handled with explicit complementarity
// pairs (s_l z_l = mu, s_u z_u = mu). Each iteration factors the regularized
// quasidefinite augmented system
//
//        [ -(Q + Theta^{-1} + rho I)   A^T  ] [dx]   [r1]
//        [          A                 delta I] [dy] = [r2]
//
// once and solves it twice (predictor, corrector) with iterative refinement
// against the unregularized matrix. For LP a crossover step identifies a basis
// from the interior solution and hands it to the simplex (lp/simplex.h).
#pragma once

#include <vector>

#include "core/model.h"
#include "util/timer.h"

namespace pramana {

struct IpmOptions {
  double tolerance = 1e-8;      // relative primal/dual residual and gap
  int maxIterations = 200;
  double primalReg = 1e-8;
  double dualReg = 1e-8;
  bool scale = true;
  bool crossover = true;        // LP only
  int logLevel = 0;
};

struct IpmStats {
  int iterations = 0;
  int regularizedPivots = 0;
  long long nnzL = 0;
  double orderingSeconds = 0, factorSeconds = 0, solveSeconds = 0;
  double finalPrimalResidual = 0, finalDualResidual = 0, finalGap = 0, finalMu = 0;
  int crossoverIterations = 0;
  bool crossoverDone = false;
};

struct IpmResult {
  Status status = Status::NotSolved;
  double objective = 0;
  std::vector<double> x, rowActivity, rowDual, reducedCost;
  std::vector<BasisStatus> colStatus, rowStatus;  // filled after crossover
  IpmStats stats;
};

// `model` must be in minimization form. Q must be convex (see checkConvexity).
IpmResult solveIpm(const Model& model, const IpmOptions& opts = {}, const Deadline* deadline = nullptr);

struct ConvexityResult {
  bool convex = true;
  double curvature = 0;              // v^T Q v for the certificate (negative => nonconvex)
  std::vector<double> certificate;   // direction v with v^T Q v < 0 when nonconvex
  int negativePivots = 0;
};
// LDL^T of Q (+ tiny shift); a clearly negative pivot yields a direction v with
// v^T Q v < 0, which is verified by direct multiplication before refusing.
ConvexityResult checkConvexity(const SparseMatrix& Q);

}  // namespace pramana
