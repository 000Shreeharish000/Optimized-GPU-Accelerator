// PRAMANA top-level solve API: routes a model to an engine, runs presolve /
// postsolve, cleans up with a warm-started simplex, and certifies the result
// on the ORIGINAL model. Every public status passes through the certifier:
// a claim that fails certification is downgraded (never reported as proven).
#pragma once

#include <string>
#include <vector>

#include "cert/certifier.h"
#include "core/model.h"
#include "util/json.h"

namespace pramana {

struct SolverOptions {
  // Engine: auto | dual | primal | ipm | pdhg | pdhg-cpu | pdhg-gpu | race
  std::string algorithm = "auto";
  double timeLimit = 0;          // seconds, 0 = unlimited
  bool presolve = true;
  bool certify = true;
  bool scale = true;
  int threads = 0;               // 0 = all hardware threads
  int logLevel = 1;              // 0 quiet, 1 summary, 2 progress, 3 detail
  uint64_t seed = 12345;
  // LP tolerances (engine side; certification uses CertTolerances)
  double primalTol = 1e-7, dualTol = 1e-7;
  double ipmTol = 1e-8;
  // First-order (PDHG)
  double pdhgTol = 1e-6;         // relative KKT tolerance
  bool pdhgCrossover = true;     // polish PDHG output with the simplex
  long long pdhgMaxIter = 1000000;
  // MIP
  double mipRelGap = 1e-4;
  double mipAbsGap = 1e-9;
  long long nodeLimit = -1;
  bool cuts = true;
  bool heuristics = true;
  std::string branching = "reliability";  // reliability | pseudocost | mostfrac | strong
  std::string nodeSelection = "bestbound-plunge";
  int mipThreads = 0;            // parallel node workers (deterministic epochs)
  // GPU
  bool allowGpu = true;
  std::string routerModelPath;   // optional calibrated router coefficients (JSON)
  CertTolerances certTol;
};

struct SolveResult {
  Status status = Status::NotSolved;      // certified status
  Status engineStatus = Status::NotSolved;
  std::string engine;                     // engine that produced the answer
  double objective = 0;                   // model sense
  double bestBound = -kInf;               // model sense (MIP / safe LP bound)
  double relativeGap = kInf;
  std::vector<double> x, rowActivity;
  std::vector<double> rowDual, reducedCost;  // model sense
  std::vector<BasisStatus> colStatus, rowStatus;
  std::vector<double> farkas, primalRay;
  Certificate certificate;
  Json telemetry = Json::object();
  double seconds = 0;
  bool hasBasis() const { return !colStatus.empty(); }
};

SolveResult solve(const Model& model, const SolverOptions& options = {});

// Result as JSON (vectors optional: large models produce large files).
Json resultToJson(const Model& model, const SolveResult& result, bool includeVectors);
// Structural features of a model (used by the router and written to telemetry).
Json modelFeatures(const Model& model);

}  // namespace pramana
