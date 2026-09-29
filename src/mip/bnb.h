// LP-based branch-and-bound / branch-and-cut for MILP, written from scratch.
//
//   * node LPs: dual simplex warm-started from the parent basis (lp/simplex.h);
//     pruning of infeasible nodes uses the simplex's verified Farkas rays;
//   * node bounds are Neumaier–Shcherbina SAFE bounds (cert/certifier.h), so the
//     reported global lower bound is rigorous, not merely a floating-point LP value;
//   * branching: reliability branching (pseudocosts initialized by limited strong
//     branching), alternatives: pseudocost, most-fractional, full strong;
//   * node selection: best bound with depth-first plunging;
//   * activity-based bound propagation at every node (mip/propagate.h);
//   * root cutting-plane loop: GMI + c-MIR + knapsack covers (mip/cuts.h) with
//     efficacy / parallelism selection;
//   * primal heuristics: simple rounding (every node), fractional diving,
//     feasibility pump (root);
//   * deterministic parallel mode: nodes are processed in epochs by worker
//     threads that each own an LP; results are merged in node-id order, so the
//     search is reproducible for a fixed thread count.
#pragma once

#include <string>
#include <vector>

#include "core/model.h"
#include "lp/simplex.h"
#include "util/json.h"
#include "util/timer.h"

namespace pramana {

struct MipOptions {
  double relGap = 1e-4;
  double absGap = 1e-9;
  long long nodeLimit = -1;
  bool cuts = true;
  int cutRounds = 8;
  int maxCutsPerRound = 100;
  bool heuristics = true;
  std::string branching = "reliability";  // reliability | pseudocost | mostfrac | strong
  int reliability = 4;                     // strong-branching evaluations per direction before trusting pseudocosts
  int strongCandidates = 8;
  int strongIterations = 60;
  int threads = 1;                         // > 1: deterministic epoch-parallel search
  int logLevel = 1;
  uint64_t seed = 12345;
  SimplexOptions lp;
};

struct MipStats {
  long long nodes = 0, lpIterations = 0, strongBranchLps = 0;
  long long farkasPrunes = 0, boundPrunes = 0, propagationPrunes = 0, integralLeaves = 0;
  long long propagationTightenings = 0;
  int maxDepth = 0, cutRoundsDone = 0;
  int cutsGomory = 0, cutsMir = 0, cutsCover = 0, cutsRejected = 0;
  int heuristicSolutions = 0;
  double rootLpBound = -kInf, rootBoundAfterCuts = -kInf;
  double rootSeconds = 0, firstIncumbentSeconds = -1;
  int incumbentUpdates = 0;
  long long safeBoundCorrections = 0;  // nodes where the safe bound was below the LP value
  double maxSafeBoundCorrection = 0;
};

struct MipResult {
  Status status = Status::NotSolved;
  double objective = kInf;    // min form, incl. offset
  double bestBound = -kInf;   // min form, rigorous
  bool boundSafe = true;
  std::vector<double> x;
  std::string incumbentSource;
  MipStats stats;
  Json progress = Json::array();  // [t, incumbent, bound] samples
};

// `model` must be a minimization MILP (integrality in colType).
MipResult solveMip(const Model& model, const MipOptions& opts = {}, const Deadline* deadline = nullptr);

}  // namespace pramana
