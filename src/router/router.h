// Measured CPU <-> GPU engine router.
//
// For an LP it predicts the wall time of each engine from cheap structural
// features with a log-linear cost model
//     log t_e = w_e . phi(model),   phi = [1, log(m+n), log(nnz), log(nnz/row), log(range), eqFrac, ...]
// whose weights are fitted offline on benchmark runs (bench/fit_router.py) and
// validated on held-out instance families (regret vs. the oracle engine).
// Built-in defaults are used when no calibrated model file is present.
// Every decision is logged with predicted times so it can be audited against
// the measured time of the engine that actually ran.
#pragma once

#include <string>
#include <vector>

#include "core/model.h"
#include "util/json.h"

namespace pramana {

struct RouterFeatures {
  double m = 0, n = 0, nnz = 0, avgRowNnz = 0, maxRowNnz = 0, coeffRange = 1, eqFraction = 0, freeFraction = 0;
  bool isMip = false, isQp = false;
  std::vector<double> phi() const;  // regression feature vector
  Json toJson() const;
};

struct RouterDecision {
  std::string engine;  // dual | ipm | pdhg-gpu | pdhg-cpu | bnb | ipm-qp
  std::string reason;
  double predictedDual = 0, predictedIpm = 0, predictedPdhgCpu = 0, predictedPdhgGpu = 0;
  bool calibrated = false;
  RouterFeatures features;
  Json toJson() const;
};

RouterFeatures computeFeatures(const Model& model);

class Router {
 public:
  Router();
  bool load(const std::string& path);  // calibrated weights (JSON), returns success
  RouterDecision decide(const Model& model, bool gpuAvailable, double tolerance) const;
  static const std::vector<std::string>& engines();
  std::vector<double> weights(const std::string& engine) const;

 private:
  std::vector<std::vector<double>> w_;  // per engine
  bool calibrated_ = false;
};

}  // namespace pramana
