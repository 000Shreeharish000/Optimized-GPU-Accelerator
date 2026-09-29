// Certified parametric LP analysis (PRAMANA's repeated-solve engine).
//
// Refinery planners answer questions such as "at which crude price does crude
// K enter / leave the slate?" or "what is crude K worth as its availability
// grows?" with HUNDREDS of case re-solves (PIMS-style crude valuation). Those
// questions are exactly the optimal value function v(theta) of a one-parameter
// LP family, which is piecewise linear (concave in a cost, convex in a bound).
//
// parametricAnalysis() computes v(theta) EXACTLY over [from, to] by tangent
// intersection ("sandwich"): each oracle call is a warm-started simplex solve
// whose optimal value AND subgradient (marginal value: x_j for a price, the
// reduced cost / row dual for an availability or rhs) are certified; two
// supporting lines intersect at the only candidate breakpoint between them, and
// convexity/concavity of v proves a segment linear when the value there lies on
// both lines. No basis ranging is used, so primal/dual degeneracy - endemic in
// refinery LPs - cannot stall the sweep. Infeasible or unbounded parameter
// regions are delimited by bisection.
//
// compareFamilyStrategies() answers the same question four ways on K sampled
// parameter values - cold re-solves, a warm-start chain, one batched PDHG run
// (GPU if available) - and against the exact parametric sweep, reporting time,
// accuracy, and the breakpoints that sampling misses.
#pragma once

#include <string>
#include <vector>

#include "core/model.h"
#include "util/json.h"

namespace pramana {

struct ParametricSpec {
  std::string kind = "cost";  // cost | rhs | lower | upper
  bool isRow = false;
  int index = -1;             // column (cost/lower/upper) or row (rhs/lower/upper)
  double from = 0, to = 1;
  int maxSegments = 1000;
};

struct ParametricSegment {
  double thetaLo = 0, thetaHi = 0;
  double objLo = 0, objHi = 0;  // model sense
  double slope = 0;             // d objective / d theta (marginal value)
  Status status = Status::NotSolved;
  bool certified = false;
  double certifiedGap = 0;
  int pivots = 0;               // simplex pivots to reach this segment's basis
  double paramVarValue = 0;     // value of the parameterized column / row activity in the segment (at midpoint)
};

struct ParametricResult {
  bool ok = false;
  std::string message;
  std::vector<ParametricSegment> segments;
  long long totalPivots = 0;
  double seconds = 0;
  double valueAt(double theta) const;  // exact piecewise-linear interpolation
  std::vector<double> breakpoints() const;
  std::string summary(const Model& model) const;
  Json toJson(const Model& model) const;
};

ParametricResult parametricAnalysis(const Model& model, const ParametricSpec& spec);

// Applies the parameter value to a copy of the model.
Model withParameter(const Model& model, const ParametricSpec& spec, double theta);

struct FamilyStrategy {
  std::string name;
  double seconds = 0;
  double maxRelError = 0;     // vs exact parametric value function
  int certified = 0;          // cases whose answer carries an accepted certificate
  int solved = 0;
  std::string note;
};

struct FamilyComparison {
  int cases = 0;
  std::vector<double> thetas;
  int breakpoints = 0;
  int breakpointsMissedBySampling = 0;
  std::vector<FamilyStrategy> strategies;
  std::string device;
  std::string summary() const;
  Json toJson() const;
};

FamilyComparison compareFamilyStrategies(const Model& model, const ParametricSpec& spec, int cases, bool allowGpu);

// Measures per-iteration PDHG time on CPU and GPU for growing synthetic LPs
// (the empirical CPU/GPU crossover used to calibrate the router).
Json calibrateBandwidth();

}  // namespace pramana
