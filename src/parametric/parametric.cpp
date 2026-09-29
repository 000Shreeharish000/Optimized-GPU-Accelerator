#include "parametric/parametric.h"

#include <algorithm>
#include <cmath>

#include "cert/certifier.h"
#include "lp/simplex.h"
#include "pdhg/pdhg.h"
#include "util/log.h"
#include "util/parallel.h"
#include "util/rng.h"
#include "util/timer.h"

namespace pramana {

namespace {

inline bool fin(double v) { return std::fabs(v) < kInfBoundThreshold; }

// Which variable bound (in the simplex's [A I] space) the parameter drives, and
// d(bound)/d(theta).  Logical s_i = -a_i^T x has bounds [-U_i, -L_i].
struct BoundParam {
  int var = -1;        // n + i for rows
  bool upperSide = true;
  double dBound = 1.0;
  bool both = false;   // equality row: both bounds move
};

BoundParam boundParam(const Model& m, const ParametricSpec& s) {
  BoundParam b;
  const int n = m.numCols();
  if (!s.isRow) {
    b.var = s.index;
    b.upperSide = s.kind == "upper";
    b.dBound = 1.0;
    return b;
  }
  const int i = s.index;
  b.var = n + i;
  b.dBound = -1.0;
  double L = m.rowLower[i], U = m.rowUpper[i];
  if (s.kind == "rhs") {
    if (L == U) {
      b.both = true;
      b.upperSide = true;
    } else if (fin(U)) {
      b.upperSide = false;  // row upper U -> logical lower -U
    } else {
      b.upperSide = true;   // row lower L -> logical upper -L
    }
  } else if (s.kind == "upper") {
    b.upperSide = false;
  } else {
    b.upperSide = true;
  }
  return b;
}

}  // namespace

Model withParameter(const Model& model, const ParametricSpec& s, double theta) {
  Model m = model;
  if (!s.isRow) {
    int j = s.index;
    if (s.kind == "cost") m.colCost[j] = theta;
    else if (s.kind == "upper") m.colUpper[j] = theta;
    else if (s.kind == "lower") m.colLower[j] = theta;
  } else {
    int i = s.index;
    double L = model.rowLower[i], U = model.rowUpper[i];
    if (s.kind == "rhs") {
      if (L == U) m.rowLower[i] = m.rowUpper[i] = theta;
      else if (fin(L) && fin(U)) {
        m.rowUpper[i] = theta;
        m.rowLower[i] = theta - (U - L);
      } else if (fin(U)) m.rowUpper[i] = theta;
      else m.rowLower[i] = theta;
    } else if (s.kind == "upper") m.rowUpper[i] = theta;
    else if (s.kind == "lower") m.rowLower[i] = theta;
  }
  return m;
}

double ParametricResult::valueAt(double theta) const {
  for (auto& s : segments) {
    if (s.status != Status::Optimal) continue;
    if (theta >= s.thetaLo - 1e-12 && theta <= s.thetaHi + 1e-12) return s.objLo + s.slope * (theta - s.thetaLo);
  }
  return std::nan("");
}

std::vector<double> ParametricResult::breakpoints() const {
  std::vector<double> b;
  for (size_t k = 1; k < segments.size(); ++k) b.push_back(segments[k].thetaLo);
  return b;
}

std::string ParametricResult::summary(const Model& model) const {
  std::string s = formatString("parametric analysis of %s: %zu segments, %lld pivots, %.3fs  [%s]\n", model.name.c_str(),
                               segments.size(), totalPivots, seconds, ok ? "complete" : message.c_str());
  s += "      theta_lo        theta_hi          obj(lo)          obj(hi)       slope   param_var  certificate\n";
  for (auto& g : segments) {
    if (g.status != Status::Optimal) {
      s += formatString("  %14.6g  %14.6g   %s\n", g.thetaLo, g.thetaHi, statusName(g.status));
      continue;
    }
    s += formatString("  %14.6g  %14.6g  %15.8g  %15.8g  %10.4g  %10.4g  %s (gap %.1e)\n", g.thetaLo, g.thetaHi, g.objLo,
                      g.objHi, g.slope, g.paramVarValue, g.certified ? "CERTIFIED" : "unverified", g.certifiedGap);
  }
  return s;
}

Json ParametricResult::toJson(const Model& model) const {
  Json j = Json::object();
  j["model"] = model.name;
  j["complete"] = ok;
  j["message"] = message;
  j["pivots"] = static_cast<double>(totalPivots);
  j["seconds"] = seconds;
  Json arr = Json::array();
  for (auto& g : segments) {
    Json o = Json::object();
    o["theta_lo"] = g.thetaLo;
    o["theta_hi"] = g.thetaHi;
    o["status"] = statusName(g.status);
    o["objective_lo"] = g.objLo;
    o["objective_hi"] = g.objHi;
    o["slope"] = g.slope;
    o["param_variable_value"] = g.paramVarValue;
    o["certified"] = g.certified;
    o["certified_gap"] = g.certifiedGap;
    o["pivots"] = g.pivots;
    arr.push(o);
  }
  j["segments"] = arr;
  j["breakpoints"] = Json(breakpoints());
  return j;
}

ParametricResult parametricAnalysis(const Model& model, const ParametricSpec& spec) {
  // Tangent-intersection ("sandwich") enumeration of the breakpoints of the
  // piecewise-linear value function v(theta). v is convex (bound / rhs
  // parameter, min form) or concave (cost parameter); an oracle call returns a
  // CERTIFIED optimal value v and a subgradient s (the marginal value) at theta.
  // Two supporting lines at a < b intersect at theta*: if v(theta*) lies on
  // them, theta* is the only breakpoint in (a, b) (by convexity); otherwise the
  // interval is split. Robust to primal/dual degeneracy (no basis ranging), and
  // needs ~2 warm-started solves per breakpoint.
  Timer timer;
  ParametricResult res;
  const Model minModel = model.toMinimization();
  const double sense = model.sense == ObjSense::Maximize ? -1.0 : 1.0;
  SimplexOptions so;
  so.perturb = true;
  Simplex lp(minModel, so);
  const bool costKind = !spec.isRow && spec.kind == "cost";
  const BoundParam bp = boundParam(model, spec);
  const int n = model.numCols();

  struct Eval {
    double theta = 0, v = 0, s = 0, paramVar = 0, gap = 0;
    Status status = Status::NotSolved;
    bool certified = false;
    int pivots = 0;
  };
  int evaluations = 0;
  auto oracle = [&](double theta) {
    Eval e;
    e.theta = theta;
    Model pm = withParameter(model, spec, theta);
    if (costKind) lp.setColCost(spec.index, sense * theta);
    else if (!spec.isRow) lp.setColBounds(spec.index, pm.colLower[spec.index], pm.colUpper[spec.index]);
    else lp.setRowBounds(spec.index, pm.rowLower[spec.index], pm.rowUpper[spec.index]);
    long long i0 = lp.stats().iterations;
    e.status = lp.solve();
    e.pivots = static_cast<int>(lp.stats().iterations - i0);
    res.totalPivots += e.pivots;
    ++evaluations;
    if (e.status != Status::Optimal) return e;
    std::vector<double> x = lp.colValues(), act = lp.rowActivities();
    e.v = sense * lp.objective();
    e.paramVar = spec.isRow ? act[spec.index] : x[spec.index];
    if (costKind) {
      e.s = x[spec.index];  // d v / d c_j = x_j (model sense)
    } else if (!spec.isRow) {
      std::vector<double> d = lp.reducedCosts();
      // Reduced cost is the marginal value of the column's active bound (min form).
      const int j = spec.index;
      bool active = bp.upperSide ? std::fabs(x[j] - pm.colUpper[j]) <= 1e-9 * (1 + std::fabs(x[j]))
                                 : std::fabs(x[j] - pm.colLower[j]) <= 1e-9 * (1 + std::fabs(x[j]));
      e.s = active ? sense * d[j] : 0.0;
    } else {
      std::vector<double> y = lp.rowDuals();
      e.s = sense * y[spec.index];  // d v / d rhs_i = y_i
    }
    SolutionClaim c;
    c.status = Status::Optimal;
    c.x = x;
    c.rowDual = lp.rowDuals();
    Certificate cert = certify(pm, c);
    e.certified = cert.accepted;
    e.gap = cert.certifiedGap;
    return e;
  };
  (void)n;

  const double span = std::max(spec.to - spec.from, 1e-300);
  const double thetaTol = 1e-10 * std::max(1.0, std::max(std::fabs(spec.from), std::fabs(spec.to)));
  auto valTol = [](double v) { return 1e-8 * (1.0 + std::fabs(v)); };

  // ---- Feasible sub-interval (bisection at the ends) ----
  std::vector<Eval> samples;
  for (int k = 0; k <= 8; ++k) samples.push_back(oracle(spec.from + span * k / 8.0));
  int firstOk = -1, lastOk = -1;
  for (int k = 0; k <= 8; ++k)
    if (samples[k].status == Status::Optimal) {
      if (firstOk < 0) firstOk = k;
      lastOk = k;
    }
  if (firstOk < 0) {
    ParametricSegment g;
    g.thetaLo = spec.from;
    g.thetaHi = spec.to;
    g.status = samples[4].status;
    res.segments.push_back(g);
    res.message = "no feasible/bounded parameter value found in the range";
    res.seconds = timer.seconds();
    return res;
  }
  auto boundary = [&](Eval good, Eval bad) {  // bisection to the feasibility boundary
    for (int it = 0; it < 60 && std::fabs(good.theta - bad.theta) > thetaTol; ++it) {
      Eval mid = oracle(0.5 * (good.theta + bad.theta));
      if (mid.status == Status::Optimal) good = mid;
      else bad = mid;
    }
    return std::make_pair(good, bad);
  };
  Eval left = samples[firstOk], right = samples[lastOk];
  Status leftStatus = Status::Optimal, rightStatus = Status::Optimal;
  if (firstOk > 0) {
    auto pr = boundary(left, samples[firstOk - 1]);
    left = pr.first;
    leftStatus = pr.second.status;
  }
  if (lastOk < 8) {
    auto pr = boundary(right, samples[lastOk + 1]);
    right = pr.first;
    rightStatus = pr.second.status;
  }

  // ---- Sandwich recursion ----
  std::vector<Eval> pts{left};  // evaluated points in increasing theta (breakpoint candidates)
  struct Job {
    Eval a, b;
  };
  std::vector<Eval> found;  // confirmed breakpoints (points where the slope changes)
  std::vector<Job> stack{{left, right}};
  std::vector<Eval> allPoints{left, right};
  int guard = 0;
  while (!stack.empty() && guard++ < 4 * spec.maxSegments) {
    Job jb = stack.back();
    stack.pop_back();
    const Eval& a = jb.a;
    const Eval& b = jb.b;
    if (b.theta - a.theta <= thetaTol) continue;
    double ds = a.s - b.s;
    // Linear on [a,b]?  chord slope equals both subgradients.
    double chord = (b.v - a.v) / (b.theta - a.theta);
    if (std::fabs(ds) <= 1e-9 * (1 + std::fabs(a.s) + std::fabs(b.s)) ||
        (std::fabs(chord - a.s) <= 1e-9 * (1 + std::fabs(chord)) && std::fabs(chord - b.s) <= 1e-9 * (1 + std::fabs(chord))))
      continue;
    double ts = (b.v - a.v + a.s * a.theta - b.s * b.theta) / (a.s - b.s);
    if (!(ts > a.theta + thetaTol && ts < b.theta - thetaTol)) {
      // Numerical corner: split in the middle.
      ts = 0.5 * (a.theta + b.theta);
    }
    Eval e = oracle(ts);
    allPoints.push_back(e);
    if (e.status != Status::Optimal) {  // should not happen inside a feasible interval (convexity)
      stack.push_back({a, e});
      continue;
    }
    double tangentA = a.v + a.s * (ts - a.theta);
    double tangentB = b.v + b.s * (ts - b.theta);
    if (std::fabs(e.v - tangentA) <= valTol(e.v) && std::fabs(e.v - tangentB) <= valTol(e.v)) {
      found.push_back(e);  // the unique breakpoint between a and b
    } else {
      stack.push_back({e, b});
      stack.push_back({a, e});
    }
  }
  // ---- Assemble segments between consecutive breakpoints (chord slopes) ----
  std::vector<Eval> nodes{left};
  std::sort(found.begin(), found.end(), [](const Eval& p, const Eval& q) { return p.theta < q.theta; });
  for (auto& f : found)
    if (f.theta > nodes.back().theta + thetaTol) nodes.push_back(f);
  if (right.theta > nodes.back().theta + thetaTol) nodes.push_back(right);
  std::sort(allPoints.begin(), allPoints.end(), [](const Eval& p, const Eval& q) { return p.theta < q.theta; });
  if (spec.from < left.theta - thetaTol) {
    ParametricSegment g;
    g.thetaLo = spec.from;
    g.thetaHi = left.theta;
    g.status = leftStatus;
    res.segments.push_back(g);
  }
  for (size_t k = 0; k + 1 < nodes.size() || (nodes.size() == 1 && k == 0); ++k) {
    ParametricSegment g;
    const Eval& p = nodes[k];
    const Eval& q = nodes.size() == 1 ? nodes[0] : nodes[k + 1];
    g.thetaLo = p.theta;
    g.thetaHi = q.theta;
    g.objLo = p.v;
    g.objHi = q.v;
    g.slope = q.theta > p.theta ? (q.v - p.v) / (q.theta - p.theta) : p.s;
    g.status = Status::Optimal;
    // Certificates: every evaluated point inside the segment must be certified.
    g.certified = true;
    g.certifiedGap = 0;
    double mid = 0.5 * (p.theta + q.theta);
    double bestDist = kInf;
    for (auto& e : allPoints) {
      if (e.theta < p.theta - thetaTol || e.theta > q.theta + thetaTol) continue;
      if (e.status == Status::Optimal) {
        g.certified = g.certified && e.certified;
        g.certifiedGap = std::max(g.certifiedGap, e.gap);
        if (std::fabs(e.theta - mid) < bestDist) {
          bestDist = std::fabs(e.theta - mid);
          g.paramVarValue = e.paramVar;
        }
        g.pivots += e.pivots;
      }
    }
    res.segments.push_back(g);
    if (nodes.size() == 1) break;
  }
  if (spec.to > right.theta + thetaTol) {
    ParametricSegment g;
    g.thetaLo = right.theta;
    g.thetaHi = spec.to;
    g.status = rightStatus;
    res.segments.push_back(g);
  }
  res.ok = true;
  res.message = formatString("%d certified solves", evaluations);
  res.seconds = timer.seconds();
  return res;
}

// ---------------------------------------------------------------------------
// Family comparison (the repeated-solve experiment)
// ---------------------------------------------------------------------------
std::string FamilyComparison::summary() const {
  std::string s = formatString("case family: %d cases, %d exact breakpoints (%d missed by uniform sampling)\n", cases,
                               breakpoints, breakpointsMissedBySampling);
  s += "  strategy                               time(s)   solved  certified   max rel. error   note\n";
  for (auto& st : strategies)
    s += formatString("  %-36s %9.4f  %6d  %9d   %14.3e   %s\n", st.name.c_str(), st.seconds, st.solved, st.certified,
                      st.maxRelError, st.note.c_str());
  return s;
}

Json FamilyComparison::toJson() const {
  Json j = Json::object();
  j["cases"] = cases;
  j["breakpoints"] = breakpoints;
  j["breakpoints_missed_by_sampling"] = breakpointsMissedBySampling;
  j["device"] = device;
  j["thetas"] = Json(thetas);
  Json arr = Json::array();
  for (auto& st : strategies) {
    Json o = Json::object();
    o["strategy"] = st.name;
    o["seconds"] = st.seconds;
    o["solved"] = st.solved;
    o["certified"] = st.certified;
    o["max_relative_error"] = st.maxRelError;
    o["note"] = st.note;
    arr.push(o);
  }
  j["strategies"] = arr;
  return j;
}

FamilyComparison compareFamilyStrategies(const Model& model, const ParametricSpec& spec, int K, bool allowGpu) {
  FamilyComparison fc;
  fc.cases = K;
  for (int k = 0; k < K; ++k) fc.thetas.push_back(K == 1 ? spec.from : spec.from + (spec.to - spec.from) * k / (K - 1));
  const Model minModel = model.toMinimization();
  const double sense = model.sense == ObjSense::Maximize ? -1.0 : 1.0;

  // Exact reference: parametric sweep.
  ParametricResult pr = parametricAnalysis(model, spec);
  fc.breakpoints = static_cast<int>(pr.breakpoints().size());
  for (auto& g : pr.segments) {
    bool sampled = false;
    for (double t : fc.thetas)
      if (t >= g.thetaLo && t <= g.thetaHi) sampled = true;
    if (!sampled) fc.breakpointsMissedBySampling++;
  }
  auto relErr = [&](double v, double theta) {
    double ref = pr.valueAt(theta);
    if (std::isnan(ref)) return 0.0;
    return std::fabs(v - ref) / std::max(1.0, std::fabs(ref));
  };
  {
    FamilyStrategy s;
    s.name = "parametric sweep (exact, certified)";
    s.seconds = pr.seconds;
    for (auto& g : pr.segments) {
      s.solved += g.status == Status::Optimal;
      s.certified += g.certified;
    }
    s.note = formatString("%zu segments, %lld pivots; certificates per segment", pr.segments.size(), pr.totalPivots);
    fc.strategies.push_back(s);
  }
  // Cold re-solves.
  {
    FamilyStrategy s;
    s.name = "cold dual simplex x K";
    Timer t;
    for (double th : fc.thetas) {
      Model pm = withParameter(model, spec, th);
      LpResult r = solveLpSimplex(pm.toMinimization());
      if (r.status == Status::Optimal) {
        s.solved++;
        s.maxRelError = std::max(s.maxRelError, relErr(sense * r.objective, th));
        SolutionClaim c;
        c.status = r.status;
        c.x = r.x;
        c.rowDual = r.rowDual;
        s.certified += certify(pm, c).accepted;
      }
    }
    s.seconds = t.seconds();
    s.note = "independent solves (incl. certification)";
    fc.strategies.push_back(s);
  }
  // Warm-start chain.
  {
    FamilyStrategy s;
    s.name = "warm-started simplex chain";
    Timer t;
    Simplex lp(minModel);
    long long piv = 0;
    for (double th : fc.thetas) {
      Model pm = withParameter(model, spec, th);
      if (!spec.isRow && spec.kind == "cost") lp.setColCost(spec.index, sense * th);
      else if (!spec.isRow) lp.setColBounds(spec.index, pm.colLower[spec.index], pm.colUpper[spec.index]);
      else lp.setRowBounds(spec.index, pm.rowLower[spec.index], pm.rowUpper[spec.index]);
      long long i0 = lp.stats().iterations;
      Status st = lp.solve();
      piv += lp.stats().iterations - i0;
      if (st == Status::Optimal) {
        s.solved++;
        s.maxRelError = std::max(s.maxRelError, relErr(sense * lp.objective(), th));
        SolutionClaim c;
        c.status = st;
        c.x = lp.colValues();
        c.rowDual = lp.rowDuals();
        s.certified += certify(pm, c).accepted;
      }
    }
    s.seconds = t.seconds();
    s.note = formatString("%lld pivots total", piv);
    fc.strategies.push_back(s);
  }
  // Batched PDHG (one run for all cases).
  for (int g = 0; g < 2; ++g) {
    bool useGpu = g == 1;
    if (useGpu && !(allowGpu && gpuAvailable())) continue;
    FamilyStrategy s;
    s.name = useGpu ? "batched PDHG (GPU) 1e-6 + safe bounds" : "batched PDHG (CPU) 1e-6 + safe bounds";
    std::vector<LpCase> cases;
    for (double th : fc.thetas) {
      Model pm = withParameter(model, spec, th).toMinimization();
      LpCase c;
      c.cost = pm.colCost;
      c.colLower = pm.colLower;
      c.colUpper = pm.colUpper;
      c.rowLower = pm.rowLower;
      c.rowUpper = pm.rowUpper;
      cases.push_back(std::move(c));
    }
    PdhgOptions po;
    po.useGpu = useGpu;
    po.tolerance = 1e-6;
    po.maxIterations = 200000;
    PdhgStats bst;
    Timer t;
    std::vector<PdhgResult> rr = solvePdhgBatch(minModel, cases, po, nullptr, &bst);
    s.seconds = t.seconds();
    for (size_t k = 0; k < rr.size(); ++k) {
      if (rr[k].status != Status::Optimal) continue;
      s.solved++;
      double th = fc.thetas[k];
      Model pm = withParameter(model, spec, th);
      s.maxRelError = std::max(s.maxRelError, relErr(sense * rr[k].objective, th));
      // A first-order answer is only "certified" to the accuracy its safe dual bound proves.
      Model pmin = pm.toMinimization();
      SafeBoundResult sb = safeDualBound(pmin.A, pmin.colCost, pmin.colLower, pmin.colUpper, pmin.rowLower,
                                         pmin.rowUpper, rr[k].rowDual);
      double lb = sb.bound + pmin.objOffset;
      if (std::isfinite(sb.bound) && std::fabs(rr[k].objective - lb) <= 1e-4 * (1 + std::fabs(rr[k].objective)))
        s.certified++;
    }
    s.note = formatString("%s, %lld iterations, kernel %.3fs, transfer %.3fs, setup %.3fs; certified = safe bound within 1e-4",
                          bst.device.c_str(), bst.iterations, bst.kernelSeconds, bst.transferSeconds, bst.setupSeconds);
    if (useGpu) fc.device = bst.device;
    fc.strategies.push_back(s);
  }
  return fc;
}

// ---------------------------------------------------------------------------
// CPU/GPU PDHG per-iteration crossover micro-benchmark
// ---------------------------------------------------------------------------
Json calibrateBandwidth() {
  Json out = Json::object();
  out["gpu"] = gpuInfoJson();
  out["cpu_threads"] = ThreadPool::instance().threads();
  Json rows = Json::array();
  Rng rng(7);
  for (int scale : {1000, 4000, 16000, 64000, 256000, 1000000}) {
    // Random sparse LP (transportation-like): ~8 nnz per column.
    int m = scale, n = 2 * scale;
    std::vector<int> ri, ci;
    std::vector<double> v;
    for (int j = 0; j < n; ++j)
      for (int e = 0; e < 4; ++e) {
        ri.push_back(rng.below(m));
        ci.push_back(j);
        v.push_back(0.5 + rng.uniform());
      }
    Model lp;
    lp.name = "calib";
    lp.A = SparseMatrix::fromTriplets(m, n, ri, ci, v);
    lp.colCost.resize(n);
    lp.colLower.assign(n, 0.0);
    lp.colUpper.assign(n, 10.0);
    lp.colType.assign(n, VarType::Continuous);
    for (int j = 0; j < n; ++j) lp.colCost[j] = rng.uniform() - 0.3;
    lp.rowLower.assign(m, 1.0);
    lp.rowUpper.assign(m, kInf);
    PdhgOptions po;
    po.tolerance = 1e-30;  // run a fixed number of iterations
    po.maxIterations = 512;
    Json r = Json::object();
    r["rows"] = m;
    r["cols"] = n;
    r["nnz"] = lp.nnz();
    for (int g = 0; g < 2; ++g) {
      po.useGpu = g == 1;
      if (po.useGpu && !gpuAvailable()) continue;
      Timer t;
      PdhgResult pr = solvePdhg(lp, po);
      double per = pr.stats.solveSeconds / std::max<long long>(1, pr.stats.iterations);
      double bytes = 2.0 * 12.0 * lp.nnz() + 10.0 * 8.0 * (n + m);
      Json e = Json::object();
      e["seconds_per_iteration"] = per;
      e["effective_GBs"] = bytes / per / 1e9;
      e["setup_seconds"] = pr.stats.setupSeconds;
      e["total_seconds"] = t.seconds();
      r[po.useGpu ? "gpu" : "cpu"] = e;
    }
    rows.push(r);
    PLOG_INFO("calibrate: nnz %d done", lp.nnz());
  }
  out["pdhg_iteration_timing"] = rows;
  return out;
}

}  // namespace pramana
