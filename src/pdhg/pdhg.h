// Restarted Halpern PDHG (primal-dual hybrid gradient) for LP, written from
// scratch, following the PDLP -> cuPDLP -> cuPDLPx line of work:
//
//   x+ = proj_[l,u]( x - tau (c - A^T y) )
//   y+ = v + sigma proj_[L,U](-v/sigma),   v = y - sigma A(2x+ - x)
//   z  = lam((1+rho) z+ - rho z) + (1-lam) z0      (reflected Halpern, lam=(k+1)/(k+2))
//
// with Ruiz + Pock-Chambolle diagonal preconditioning, objective / rhs
// rescaling, constant step eta = 0.998/||A||_2 split by a primal weight omega
// (updated at restarts), adaptive restarts on the fixed-point residual, and
// termination on the relative KKT error measured on the UNSCALED problem.
//
// Two backends execute the *same* operation sequence: a multi-threaded CPU
// backend and a CUDA backend (own kernels, pdhg/kernels.cpp). The batched API
// solves k LPs that share A (repeated-solve families: price / availability
// cases) as one SpMM per iteration.
//
// PDHG output is approximate (1e-4..1e-8 relative); the orchestrator turns it
// into certified answers by crossover (lp/crossover.h) and/or the safe dual
// bound of the certifier, which is valid for ANY dual vector.
#pragma once

#include <string>
#include <vector>

#include "core/model.h"
#include "util/json.h"
#include "util/timer.h"

namespace pramana {

struct PdhgOptions {
  double tolerance = 1e-6;       // relative KKT tolerance
  long long maxIterations = 1000000;
  int checkEvery = 64;           // iterations between (device->host) convergence checks
  bool useGpu = false;
  int ruizIterations = 10;
  bool pockChambolle = true;
  bool rescaleObjectiveAndBounds = true;
  double reflection = 1.0;       // rho in [0,1]
  int logInterval = 0;           // log every N checks (0 = off)
};

struct PdhgStats {
  long long iterations = 0;
  int restarts = 0;
  double primalResidual = 0, dualResidual = 0, gap = 0;
  double primalWeight = 1;
  double setupSeconds = 0, solveSeconds = 0;
  double kernelSeconds = 0;        // GPU: device time of the iteration loop (events)
  double transferSeconds = 0;      // host<->device copies (setup + checks + results)
  long long bytesUploaded = 0, bytesDownloaded = 0, kernelLaunches = 0;
  double spectralNorm = 0;
  std::string backend = "cpu", device;
  Json toJson() const;
};

struct PdhgResult {
  Status status = Status::NotSolved;  // Optimal (to tolerance) | IterationLimit | TimeLimit
  double objective = 0, dualObjective = 0;
  std::vector<double> x, rowDual;     // original space, min form
  long long iterations = 0;
  double relativeKkt = 0;
  PdhgStats stats;                    // single solve: full stats
};

// A case in a family: empty vectors mean "same as the base model".
struct LpCase {
  std::vector<double> cost, colLower, colUpper, rowLower, rowUpper;
};

// `lp` must be a minimization LP.
PdhgResult solvePdhg(const Model& lp, const PdhgOptions& opts = {}, const Deadline* deadline = nullptr);

// Solves all cases in one batched run (shared A). `batchStats` gets totals.
std::vector<PdhgResult> solvePdhgBatch(const Model& base, const std::vector<LpCase>& cases, const PdhgOptions& opts,
                                       const Deadline* deadline = nullptr, PdhgStats* batchStats = nullptr);

bool gpuAvailable();
Json gpuInfoJson();

}  // namespace pramana
