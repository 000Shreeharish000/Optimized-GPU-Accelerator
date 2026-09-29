#include "router/router.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "util/log.h"

namespace pramana {

namespace {
const std::vector<std::string> kEngines = {"dual", "ipm", "pdhg-cpu", "pdhg-gpu"};
inline double lg(double v) { return std::log(std::max(v, 1.0)); }
}  // namespace

const std::vector<std::string>& Router::engines() { return kEngines; }

std::vector<double> RouterFeatures::phi() const {
  return {1.0,
          lg(m + n),
          lg(nnz),
          lg(avgRowNnz),
          lg(maxRowNnz),
          lg(coeffRange),
          eqFraction,
          freeFraction};
}

Json RouterFeatures::toJson() const {
  Json j = Json::object();
  j["rows"] = m;
  j["cols"] = n;
  j["nnz"] = nnz;
  j["avg_row_nnz"] = avgRowNnz;
  j["max_row_nnz"] = maxRowNnz;
  j["coefficient_range"] = coeffRange;
  j["equality_row_fraction"] = eqFraction;
  j["free_column_fraction"] = freeFraction;
  j["is_mip"] = isMip;
  j["is_qp"] = isQp;
  return j;
}

Json RouterDecision::toJson() const {
  Json j = Json::object();
  j["engine"] = engine;
  j["reason"] = reason;
  j["calibrated_model"] = calibrated;
  Json p = Json::object();
  p["dual"] = predictedDual;
  p["ipm"] = predictedIpm;
  p["pdhg-cpu"] = predictedPdhgCpu;
  p["pdhg-gpu"] = predictedPdhgGpu;
  j["predicted_seconds"] = p;
  j["features"] = features.toJson();
  return j;
}

RouterFeatures computeFeatures(const Model& model) {
  RouterFeatures f;
  f.m = model.numRows();
  f.n = model.numCols();
  f.nnz = model.nnz();
  f.isMip = model.isMip();
  f.isQp = model.isQp();
  std::vector<int> rowCount(model.numRows(), 0);
  double mx = 0, mn = kInf;
  for (int j = 0; j < model.numCols(); ++j)
    for (int k = model.A.start[j]; k < model.A.start[j + 1]; ++k) {
      rowCount[model.A.index[k]]++;
      double a = std::fabs(model.A.value[k]);
      if (a > 0) {
        mx = std::max(mx, a);
        mn = std::min(mn, a);
      }
    }
  for (int c : rowCount) f.maxRowNnz = std::max(f.maxRowNnz, static_cast<double>(c));
  f.avgRowNnz = f.m > 0 ? f.nnz / f.m : 0;
  f.coeffRange = mx > 0 ? mx / mn : 1;
  int eq = 0;
  for (int i = 0; i < model.numRows(); ++i) eq += model.rowLower[i] == model.rowUpper[i];
  f.eqFraction = f.m > 0 ? eq / f.m : 0;
  int fr = 0;
  for (int j = 0; j < model.numCols(); ++j)
    fr += !(std::fabs(model.colLower[j]) < kInfBoundThreshold) && !(std::fabs(model.colUpper[j]) < kInfBoundThreshold);
  f.freeFraction = f.n > 0 ? fr / f.n : 0;
  return f;
}

Router::Router() {
  // Built-in defaults (log seconds). Rough priors from our own Netlib/Mittelmann-size
  // runs on a laptop CPU + RTX-class GPU; replaced by bench/fit_router.py output.
  //                 1      log(m+n) log(nnz) log(nnz/row) log(maxrow) log(range) eq    free
  w_ = {
      {-13.0, 1.20, 0.35, 0.30, 0.00, 0.05, 0.20, 0.30},   // dual simplex
      {-12.0, 0.40, 1.00, 0.40, 0.10, 0.03, 0.00, 0.00},   // ipm
      {-11.5, 0.60, 0.80, 0.00, 0.10, 0.15, 0.30, 0.00},   // pdhg-cpu
      {-6.0, 0.20, 0.35, 0.00, 0.05, 0.15, 0.30, 0.00},    // pdhg-gpu (launch-latency floor)
  };
  const char* env = std::getenv("PRAMANA_ROUTER_MODEL");
  if (env) load(env);
}

bool Router::load(const std::string& path) {
  try {
    Json j = Json::parse(readTextFile(path));
    std::vector<std::vector<double>> w;
    for (auto& e : kEngines) {
      if (!j.at("weights").has(e)) return false;
      w.push_back(j.at("weights").at(e).numVector());
    }
    for (auto& v : w)
      if (v.size() != w_[0].size()) return false;
    w_ = w;
    calibrated_ = true;
    return true;
  } catch (std::exception& e) {
    PLOG_DETAIL("router: cannot load %s (%s)", path.c_str(), e.what());
    return false;
  }
}

std::vector<double> Router::weights(const std::string& engine) const {
  for (size_t e = 0; e < kEngines.size(); ++e)
    if (kEngines[e] == engine) return w_[e];
  return {};
}

RouterDecision Router::decide(const Model& model, bool gpu, double tolerance) const {
  RouterDecision d;
  d.features = computeFeatures(model);
  d.calibrated = calibrated_;
  if (d.features.isMip) {
    d.engine = "bnb";
    d.reason = "integer variables: LP-based branch-and-cut (node LPs by warm-started dual simplex on CPU)";
    return d;
  }
  if (d.features.isQp) {
    d.engine = "ipm-qp";
    d.reason = "quadratic objective: interior point on the quasidefinite KKT system";
    return d;
  }
  std::vector<double> phi = d.features.phi();
  auto pred = [&](size_t e) {
    double s = 0;
    for (size_t q = 0; q < phi.size(); ++q) s += w_[e][q] * phi[q];
    return std::exp(s);
  };
  d.predictedDual = pred(0);
  d.predictedIpm = pred(1);
  d.predictedPdhgCpu = pred(2);
  d.predictedPdhgGpu = gpu ? pred(3) : kInf;
  // First-order output must be certified: add the (predicted) crossover polish,
  // a fraction of a simplex solve from a near-optimal warm start.
  const double polish = 0.15 * d.predictedDual;
  double bestFirstOrder = std::min(d.predictedPdhgCpu, d.predictedPdhgGpu) + polish;
  // Tight tolerances favour the simplex (PDHG tails are slow).
  double tightPenalty = tolerance < 1e-7 ? 4.0 : 1.0;
  bestFirstOrder *= tightPenalty;
  struct Opt {
    const char* name;
    double t;
  };
  std::vector<Opt> opts = {{"dual", d.predictedDual},
                           {"ipm", d.predictedIpm * 1.1},  // + crossover
                           {"pdhg-cpu", d.predictedPdhgCpu * tightPenalty + polish},
                           {"pdhg-gpu", d.predictedPdhgGpu * tightPenalty + polish}};
  auto best = std::min_element(opts.begin(), opts.end(), [](const Opt& a, const Opt& b) { return a.t < b.t; });
  d.engine = best->name;
  d.reason = formatString("min predicted time %.3gs (dual %.3g, ipm %.3g, pdhg-cpu %.3g, pdhg-gpu %s)%s", best->t,
                          d.predictedDual, d.predictedIpm, d.predictedPdhgCpu,
                          gpu ? formatString("%.3g", d.predictedPdhgGpu).c_str() : "n/a",
                          calibrated_ ? " [calibrated]" : " [default priors]");
  (void)bestFirstOrder;
  return d;
}

}  // namespace pramana
