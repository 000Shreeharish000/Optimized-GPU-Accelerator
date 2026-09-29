#include "api/solver.h"

#include <atomic>
#include <cmath>
#include <mutex>
#include <thread>

#include "ipm/ipm.h"
#include "lp/crossover.h"
#include "lp/scaling.h"
#include "lp/simplex.h"
#include "mip/bnb.h"
#include "pdhg/pdhg.h"
#include "presolve/presolve.h"
#include "router/router.h"
#include "util/log.h"
#include "util/parallel.h"
#include "util/timer.h"

namespace pramana {

namespace {

inline bool fin(double v) { return std::fabs(v) < kInfBoundThreshold; }

// Engine output in the space of the model it was run on (minimization form).
struct EngineOut {
  Status status = Status::NotSolved;
  std::string engine;
  double objective = 0, bound = -kInf;
  bool boundSafe = false;
  std::vector<double> x, y, d;
  std::vector<BasisStatus> cs, rs;
  std::vector<double> farkas, ray;
  Json stats = Json::object();
};

Json simplexStatsJson(const SimplexStats& s) {
  Json j = Json::object();
  j["iterations"] = static_cast<double>(s.iterations);
  j["dual_phase1_iterations"] = static_cast<double>(s.dualPhase1Iterations);
  j["dual_phase2_iterations"] = static_cast<double>(s.dualPhase2Iterations);
  j["primal_iterations"] = static_cast<double>(s.primalIterations);
  j["degenerate_pivots"] = static_cast<double>(s.degeneratePivots);
  j["degenerate_fraction"] = s.iterations ? static_cast<double>(s.degeneratePivots) / s.iterations : 0.0;
  j["bound_flips"] = static_cast<double>(s.boundFlips);
  j["refactorizations"] = static_cast<double>(s.refactorizations);
  j["numerical_refactors"] = static_cast<double>(s.numericalRefactors);
  j["rejected_pivots"] = static_cast<double>(s.rejectedPivots);
  j["rejected_farkas_rays"] = static_cast<double>(s.rejectedRays);
  j["singular_basis_repairs"] = static_cast<double>(s.singularRepairs);
  j["cost_shifts"] = static_cast<double>(s.costShifts);
  j["max_cost_shift"] = s.maxCostShift;
  j["perturbations"] = static_cast<double>(s.perturbations);
  j["max_perturbation"] = s.maxPerturbation;
  j["stall_recoveries"] = static_cast<double>(s.stallRecoveries);
  j["dual_polish_passes"] = static_cast<double>(s.polishIterations);
  j["scaling_range_before"] = s.scalingRangeBefore;
  j["scaling_range_after"] = s.scalingRangeAfter;
  j["seconds"] = s.timeSeconds;
  Json prof = Json::object();
  prof["chuzr"] = s.tChuzr;
  prof["btran"] = s.tBtran;
  prof["price"] = s.tPrice;
  prof["chuzc"] = s.tChuzc;
  prof["ftran_dse"] = s.tFtran;
  prof["updates"] = s.tUpdate;
  prof["factorize"] = s.tFactor;
  j["profile_seconds"] = prof;
  return j;
}

EngineOut fromLp(const LpResult& r, const char* name) {
  EngineOut o;
  o.status = r.status;
  o.engine = name;
  o.objective = r.objective;
  o.x = r.x;
  o.y = r.rowDual;
  o.d = r.reducedCost;
  o.cs = r.colStatus;
  o.rs = r.rowStatus;
  o.farkas = r.farkas;
  o.ray = r.primalRay;
  o.stats["simplex"] = simplexStatsJson(r.stats);
  return o;
}

SimplexOptions simplexOptionsFrom(const SolverOptions& opt) {
  SimplexOptions so;
  so.primalFeasTol = opt.primalTol;
  so.dualFeasTol = opt.dualTol;
  so.scale = opt.scale;
  so.seed = opt.seed;
  so.logInterval = opt.logLevel >= 3 ? 1000 : 0;
  return so;
}

// Trivial LP without rows: every column at the bound its cost prefers.
EngineOut solveNoRows(const Model& m) {
  EngineOut o;
  o.engine = "trivial";
  o.status = Status::Optimal;
  o.x.assign(m.numCols(), 0.0);
  o.d = m.colCost;
  o.cs.assign(m.numCols(), BasisStatus::Lower);
  for (int j = 0; j < m.numCols(); ++j) {
    double c = m.colCost[j], l = m.colLower[j], u = m.colUpper[j];
    double v = c > 0 ? l : (c < 0 ? u : (fin(l) ? l : (fin(u) ? u : 0.0)));
    if (!fin(v)) {
      o.status = Status::Unbounded;
      o.ray.assign(m.numCols(), 0.0);
      o.ray[j] = c > 0 ? -1.0 : 1.0;
      v = fin(l) ? l : (fin(u) ? u : 0.0);
    }
    o.x[j] = v;
    o.cs[j] = v == l ? BasisStatus::Lower : (v == u ? BasisStatus::Upper : BasisStatus::Zero);
  }
  o.objective = m.objective(o.x);
  return o;
}

EngineOut runLpEngine(const std::string& engine, const Model& work, const SolverOptions& opt, const Deadline* dl) {
  if (work.numRows() == 0) return solveNoRows(work);
  SimplexOptions so = simplexOptionsFrom(opt);
  if (engine == "dual" || engine == "primal" || engine == "simplex") {
    return fromLp(solveLpSimplex(work, so, dl), "dual-simplex");
  }
  if (engine == "ipm") {
    IpmOptions io;
    io.tolerance = opt.ipmTol;
    io.scale = opt.scale;
    io.logLevel = opt.logLevel >= 3 ? 1 : 0;
    IpmResult r = solveIpm(work, io, dl);
    EngineOut o;
    o.engine = "ipm+crossover";
    o.status = r.status;
    o.objective = r.objective;
    o.x = r.x;
    o.y = r.rowDual;
    o.d = r.reducedCost;
    o.cs = r.colStatus;
    o.rs = r.rowStatus;
    Json s = Json::object();
    s["iterations"] = r.stats.iterations;
    s["nnz_L"] = static_cast<double>(r.stats.nnzL);
    s["regularized_pivots"] = r.stats.regularizedPivots;
    s["ordering_seconds"] = r.stats.orderingSeconds;
    s["factor_seconds"] = r.stats.factorSeconds;
    s["solve_seconds"] = r.stats.solveSeconds;
    s["final_relative_primal_residual"] = r.stats.finalPrimalResidual;
    s["final_relative_dual_residual"] = r.stats.finalDualResidual;
    s["final_relative_gap"] = r.stats.finalGap;
    s["crossover_iterations"] = r.stats.crossoverIterations;
    o.stats["ipm"] = s;
    return o;
  }
  if (engine == "pdhg" || engine == "pdhg-cpu" || engine == "pdhg-gpu") {
    PdhgOptions po;
    po.tolerance = opt.pdhgTol;
    po.maxIterations = opt.pdhgMaxIter;
    po.useGpu = engine != "pdhg-cpu" && opt.allowGpu;
    po.logInterval = opt.logLevel >= 3 ? 20 : 0;
    PdhgResult r = solvePdhg(work, po, dl);
    EngineOut o;
    o.engine = r.stats.backend == "gpu" ? "pdhg-gpu" : "pdhg-cpu";
    o.status = r.status;
    o.objective = r.objective;
    o.x = r.x;
    o.y = r.rowDual;
    o.stats["pdhg"] = r.stats.toJson();
    o.stats["pdhg"]["relative_kkt"] = r.relativeKkt;
    o.stats["pdhg"]["objective"] = r.objective;
    o.stats["pdhg"]["dual_objective"] = r.dualObjective;
    if (opt.pdhgCrossover && !(dl && dl->cancelled())) {
      Timer ct;
      LpResult c = crossover(work, r.x, r.rowDual, so, dl);
      Json cj = simplexStatsJson(c.stats);
      cj["seconds"] = ct.seconds();
      o.stats["crossover"] = cj;
      if (c.status == Status::Optimal || c.status == Status::Infeasible || c.status == Status::Unbounded) {
        EngineOut f = fromLp(c, "");
        f.engine = o.engine + "+crossover";
        f.stats = o.stats;
        return f;
      }
    }
    return o;
  }
  throw PramanaError("unknown LP engine '" + engine + "'");
}

bool isProof(Status s) { return s == Status::Optimal || s == Status::Infeasible || s == Status::Unbounded; }

EngineOut raceLp(const Model& work, const SolverOptions& opt, const Deadline* outer, Json& raceInfo) {
  std::string other = (opt.allowGpu && gpuAvailable()) ? "pdhg-gpu" : "ipm";
  std::atomic<bool> stop{false};
  std::mutex mu;
  EngineOut winner;
  bool have = false;
  std::string first;
  Timer t;
  auto runner = [&](std::string eng) {
    Deadline dl(outer ? outer->remaining() : 0);
    if (outer && outer->limit() <= 0) dl = Deadline(0);
    dl.setCancelFlag(&stop);
    EngineOut o;
    try {
      o = runLpEngine(eng, work, opt, &dl);
    } catch (std::exception& e) {
      o.status = Status::Error;
      o.engine = eng;
    }
    std::lock_guard<std::mutex> lock(mu);
    raceInfo[eng] = formatString("%s after %.3fs", statusName(o.status), t.seconds());
    if (!have && isProof(o.status)) {
      have = true;
      winner = std::move(o);
      first = eng;
      stop = true;
    } else if (!have && winner.engine.empty()) {
      winner = std::move(o);  // keep something in case nobody proves anything
    }
  };
  std::thread a(runner, std::string("dual")), b(runner, other);
  a.join();
  b.join();
  raceInfo["winner"] = first.empty() ? "none" : first;
  return winner;
}

}  // namespace

Json modelFeatures(const Model& model) {
  RouterFeatures f = computeFeatures(model);
  Json j = f.toJson();
  j["name"] = model.name;
  j["integers"] = model.numIntegers();
  j["q_nnz"] = model.Q.nnz();
  j["sense"] = model.sense == ObjSense::Maximize ? "max" : "min";
  return j;
}

SolveResult solve(const Model& input, const SolverOptions& opt) {
  Timer total;
  SolveResult res;
  Deadline dl(opt.timeLimit);
  if (opt.threads > 0) ThreadPool::setThreads(opt.threads);
  input.validate();
  const double sense = input.sense == ObjSense::Maximize ? -1.0 : 1.0;
  const Model minModel = input.toMinimization();
  const bool isMip = input.isMip(), isQp = input.isQp();
  Json& tel = res.telemetry;
  tel["model"] = modelFeatures(input);
  Json timing = Json::object();

  if (isMip && isQp) {
    res.status = res.engineStatus = Status::Error;
    tel["error"] = "MIQP is not supported in v1 (B&B over QP relaxations is a design hook; see docs/ROADMAP.md)";
    return res;
  }

  // ---- Convexity (QP) ----
  if (isQp) {
    Timer t;
    ConvexityResult cv = checkConvexity(minModel.Q);
    Json c = Json::object();
    c["convex"] = cv.convex;
    c["negative_pivots"] = cv.negativePivots;
    c["certificate_curvature"] = cv.curvature;
    tel["convexity"] = c;
    timing["convexity"] = t.seconds();
    if (!cv.convex) {
      res.status = res.engineStatus = Status::Nonconvex;
      res.primalRay = cv.certificate;  // direction v with v^T Q v < 0 (verifiable)
      tel["convexity"]["certificate_direction_norm"] = static_cast<double>(cv.certificate.size());
      res.seconds = total.seconds();
      tel["timing"] = timing;
      return res;
    }
  }

  // ---- Routing ----
  Router router;
  if (!opt.routerModelPath.empty()) router.load(opt.routerModelPath);
  // Only LPs can be routed to the GPU; avoid initializing CUDA/NVRTC otherwise.
  // In auto mode small LPs skip the probe (~0.5-1s of cuInit + JIT): below ~5e4
  // nonzeros kernel-launch latency alone makes a GPU first-order method lose.
  const bool probe = opt.algorithm == "auto" ? input.nnz() >= 50000
                                             : (opt.algorithm.rfind("pdhg", 0) == 0 || opt.algorithm == "race");
  const bool gpu = opt.allowGpu && !isMip && !isQp && probe ? gpuAvailable() : false;
  RouterDecision dec = router.decide(minModel, gpu, std::min(opt.primalTol, opt.dualTol));
  std::string engine = opt.algorithm == "auto" ? dec.engine : opt.algorithm;
  if (isMip) engine = "bnb";
  if (isQp) engine = "ipm-qp";
  if (engine == "pdhg") engine = gpu ? "pdhg-gpu" : "pdhg-cpu";
  tel["router"] = dec.toJson();
  tel["router"]["requested"] = opt.algorithm;
  tel["router"]["selected"] = engine;
  if (engine == "pdhg-gpu" || engine == "race") tel["gpu"] = gpuInfoJson();  // (initializes CUDA)

  // ---- Presolve ----
  Presolver pre;
  bool usePre = opt.presolve;  // QP: row reductions only (see presolve.h)
  const Model* work = &minModel;
  if (usePre) {
    Timer t;
    PresolveStatus ps = pre.run(minModel, isMip);
    timing["presolve"] = t.seconds();
    const PresolveStats& s = pre.stats();
    Json pj = Json::object();
    pj["status"] = ps == PresolveStatus::Reduced ? "reduced"
                   : ps == PresolveStatus::Unchanged ? "unchanged"
                   : ps == PresolveStatus::Infeasible ? "infeasible (re-solving original for a certificate)"
                   : ps == PresolveStatus::Unbounded ? "unbounded (re-solving original for a certificate)"
                                                     : "empty";
    pj["rows_removed"] = s.rowsRemoved;
    pj["cols_removed"] = s.colsRemoved;
    pj["nnz_removed"] = s.nnzRemoved;
    pj["singleton_rows"] = s.singletonRows;
    pj["forcing_rows"] = s.forcingRows;
    pj["redundant_rows"] = s.redundantRows;
    pj["empty_rows"] = s.emptyRows;
    pj["fixed_cols"] = s.fixedCols;
    pj["empty_cols"] = s.emptyCols;
    pj["dominated_cols"] = s.dominatedCols;
    pj["free_col_singletons"] = s.freeColSingletons;
    pj["bounds_tightened"] = s.boundsTightened;
    pj["passes"] = s.passes;
    tel["presolve"] = pj;
    if (ps == PresolveStatus::Reduced || ps == PresolveStatus::Empty) {
      work = &pre.reduced();
      if (opt.logLevel >= 1)
        PLOG_INFO("presolve: %d rows, %d cols removed (%.3fs)", s.rowsRemoved, s.colsRemoved, t.seconds());
    } else {
      usePre = false;
    }
  }

  // ---- Engine ----
  Timer et;
  EngineOut out;
  if (opt.logLevel >= 1)
    PLOG_INFO("engine: %s  (%s)", engine.c_str(), opt.algorithm == "auto" ? dec.reason.c_str() : "user selected");
  if (isMip) {
    MipOptions mo;
    mo.relGap = opt.mipRelGap;
    mo.absGap = opt.mipAbsGap;
    mo.nodeLimit = opt.nodeLimit;
    mo.cuts = opt.cuts;
    mo.heuristics = opt.heuristics;
    mo.branching = opt.branching;
    mo.threads = opt.mipThreads;
    mo.logLevel = opt.logLevel;
    mo.seed = opt.seed;
    mo.lp = simplexOptionsFrom(opt);
    MipResult r;
    const int racers = opt.mipThreads > 1 ? opt.mipThreads : 1;
    Json raceJson = Json::object();
    if (racers == 1) {
      r = solveMip(*work, mo, &dl);
    } else {
      // Concurrent MIP racing: diversified branch-and-cut runs on separate threads.
      // First proof wins and cancels the others; otherwise combine the best
      // incumbent with the MAXIMUM of the rigorous lower bounds (all valid).
      static const char* rules[] = {"reliability", "pseudocost", "reliability", "mostfrac", "strong", "pseudocost"};
      std::vector<MipResult> res(racers);
      std::atomic<bool> stop{false};
      std::atomic<int> winner{-1};
      std::vector<std::thread> th;
      for (int k = 0; k < racers; ++k)
        th.emplace_back([&, k]() {
          MipOptions mk = mo;
          mk.branching = rules[k % 6];
          mk.seed = opt.seed + 7919ULL * k;
          mk.lp.seed = mk.seed;
          mk.logLevel = k == 0 ? mo.logLevel : 0;
          if (k % 2 == 1) mk.cutRounds = std::max(2, mo.cutRounds / 2);
          Deadline d(opt.timeLimit);
          d.setCancelFlag(&stop);
          res[k] = solveMip(*work, mk, &d);
          if (res[k].status == Status::Optimal || res[k].status == Status::Infeasible) {
            int expected = -1;
            if (winner.compare_exchange_strong(expected, k)) stop = true;
          }
        });
      for (auto& t : th) t.join();
      int w = winner.load();
      if (w >= 0) {
        r = res[w];
      } else {
        int best = 0;
        double bound = -kInf;
        bool safe = true;
        for (int k = 0; k < racers; ++k) {
          if (res[k].objective < res[best].objective) best = k;
          if (res[k].bestBound > bound) bound = res[k].bestBound, safe = res[k].boundSafe;
        }
        r = res[best];
        r.bestBound = std::min(bound, r.objective);
        r.boundSafe = safe;
        r.status = dl.expired() ? Status::TimeLimit : r.status;
      }
      raceJson["racers"] = racers;
      raceJson["winner"] = w;
      raceJson["winner_branching"] = w >= 0 ? rules[w % 6] : "none";
    }
    out.engine = "branch-and-cut";
    out.status = r.status;
    out.objective = r.objective;
    out.bound = r.bestBound;
    out.boundSafe = r.boundSafe;
    out.x = r.x;
    const MipStats& s = r.stats;
    Json mj = Json::object();
    mj["nodes"] = static_cast<double>(s.nodes);
    mj["lp_iterations"] = static_cast<double>(s.lpIterations);
    mj["strong_branch_lps"] = static_cast<double>(s.strongBranchLps);
    mj["farkas_certified_prunes"] = static_cast<double>(s.farkasPrunes);
    mj["bound_prunes"] = static_cast<double>(s.boundPrunes);
    mj["propagation_prunes"] = static_cast<double>(s.propagationPrunes);
    mj["propagation_tightenings"] = static_cast<double>(s.propagationTightenings);
    mj["integral_leaves"] = static_cast<double>(s.integralLeaves);
    mj["max_depth"] = s.maxDepth;
    mj["cut_rounds"] = s.cutRoundsDone;
    mj["cuts_gomory"] = s.cutsGomory;
    mj["cuts_mir"] = s.cutsMir;
    mj["cuts_cover"] = s.cutsCover;
    mj["cuts_rejected_numerics"] = s.cutsRejected;
    mj["root_lp_bound"] = sense * s.rootLpBound;
    mj["root_bound_after_cuts"] = sense * s.rootBoundAfterCuts;
    mj["root_seconds"] = s.rootSeconds;
    mj["first_incumbent_seconds"] = s.firstIncumbentSeconds;
    mj["heuristic_solutions"] = s.heuristicSolutions;
    mj["incumbent_updates"] = s.incumbentUpdates;
    mj["incumbent_source"] = r.incumbentSource;
    mj["safe_bound_corrections"] = static_cast<double>(s.safeBoundCorrections);
    mj["max_safe_bound_correction"] = s.maxSafeBoundCorrection;
    mj["bound_is_safe"] = r.boundSafe;
    mj["progress_t_incumbent_bound_minform"] = r.progress;
    if (racers > 1) mj["concurrent_race"] = raceJson;
    out.stats["mip"] = mj;
  } else if (isQp) {
    IpmOptions io;
    io.tolerance = opt.ipmTol;
    io.scale = opt.scale;
    io.logLevel = opt.logLevel >= 3 ? 1 : 0;
    IpmResult r = solveIpm(*work, io, &dl);
    out.engine = "ipm-qp";
    out.status = r.status;
    out.objective = r.objective;
    out.x = r.x;
    out.y = r.rowDual;
    out.d = r.reducedCost;
    Json s = Json::object();
    s["iterations"] = r.stats.iterations;
    s["nnz_L"] = static_cast<double>(r.stats.nnzL);
    s["regularized_pivots"] = r.stats.regularizedPivots;
    s["factor_seconds"] = r.stats.factorSeconds;
    s["final_relative_primal_residual"] = r.stats.finalPrimalResidual;
    s["final_relative_dual_residual"] = r.stats.finalDualResidual;
    s["final_relative_gap"] = r.stats.finalGap;
    out.stats["ipm"] = s;
  } else if (engine == "race") {
    Json raceInfo = Json::object();
    out = raceLp(*work, opt, &dl, raceInfo);
    out.stats["race"] = raceInfo;
  } else {
    out = runLpEngine(engine, *work, opt, &dl);
  }
  timing["engine"] = et.seconds();
  res.engineStatus = out.status;
  res.engine = out.engine;
  for (auto& kv : out.stats.items()) tel[kv.first] = kv.second;

  // ---- Postsolve + cleanup to the original model ----
  Timer pt;
  EngineOut fin_ = out;
  if (usePre) {
    if (out.status == Status::Infeasible || out.status == Status::Unbounded ||
        out.status == Status::InfeasibleOrUnbounded) {
      // Rays live in the reduced space; re-derive the proof on the original model.
      if (isQp) {
        // No LP-based proof for a QP (it would ignore Q): keep the engine status, unproven.
        fin_.status = out.status == Status::InfeasibleOrUnbounded ? Status::InfeasibleOrUnbounded : out.status;
        fin_.x.clear();
      } else if (isMip) {
        MipOptions mo;
        mo.logLevel = 0;
        mo.lp = simplexOptionsFrom(opt);
        mo.nodeLimit = 0;
        // Infeasible MIP: the LP relaxation of the original often proves it; if not, keep engine status.
        LpResult lr = solveLpSimplex(minModel, mo.lp, &dl);
        if (lr.status == Status::Infeasible) fin_ = fromLp(lr, "dual-simplex (original)");
        else fin_.status = out.status, fin_.x.clear();
        fin_.engine = out.engine;
      } else {
        fin_ = fromLp(solveLpSimplex(minModel, simplexOptionsFrom(opt), &dl), "dual-simplex (original, proof)");
        fin_.engine = out.engine + " -> proof by dual simplex on original";
      }
    } else if (!out.x.empty() || (work->numCols() == 0 && out.status == Status::Optimal)) {
      Presolver::Solution rs;
      rs.x = out.x;
      rs.rowDual = out.y;
      rs.reducedCost = out.d;
      rs.colStatus = out.cs;
      rs.rowStatus = out.rs;
      rs.hasDual = !out.y.empty();
      rs.hasBasis = !out.cs.empty();
      Presolver::Solution full = pre.postsolve(rs);
      fin_.x = full.x;
      fin_.y = full.hasDual ? full.rowDual : std::vector<double>();
      fin_.d = full.hasDual ? full.reducedCost : std::vector<double>();
      fin_.cs.clear();
      fin_.rs.clear();
      if (!isMip && !isQp && out.status == Status::Optimal) {
        // Warm-started cleanup on the original model (usually 0 pivots): gives a
        // basis, duals and a certificate consistent with the ORIGINAL model.
        SimplexOptions so = simplexOptionsFrom(opt);
        so.perturb = false;
        LpResult c = full.hasBasis ? solveLpSimplex(minModel, so, &dl, &full.colStatus, &full.rowStatus)
                                   : crossover(minModel, full.x, full.rowDual, so, &dl);
        Json cj = simplexStatsJson(c.stats);
        tel["postsolve_cleanup"] = cj;
        if (c.status == Status::Optimal) {
          std::string eng = fin_.engine;
          fin_ = fromLp(c, "");
          fin_.engine = eng;
        }
      }
      if (isMip) fin_.objective = minModel.objective(fin_.x);
    }
  }
  timing["postsolve"] = pt.seconds();

  // ---- Certification on the ORIGINAL model ----
  Timer ct;
  SolutionClaim claim;
  claim.status = fin_.status;
  claim.x = fin_.x;
  claim.rowDual = fin_.y;
  claim.farkas = fin_.farkas;
  claim.primalRay = fin_.ray;
  if (isMip) {
    claim.bestBound = sense * out.bound;
    claim.boundIsSafe = out.boundSafe;
    if (claim.x.empty() && claim.status == Status::Optimal) claim.status = Status::NumericalFailure;
  }
  if (opt.certify) {
    res.certificate = certify(input, claim, opt.certTol);
    res.status = res.certificate.certified;
    // A feasible point whose optimality could not be proven (e.g. a first-order
    // answer at 1e-6 that crossover did not finish): report it as approximate,
    // with the rigorous bound, rather than as a numerical failure.
    if (!res.certificate.accepted && fin_.status == Status::Optimal && !fin_.x.empty() &&
        res.certificate.maxPrimalViolation <= opt.certTol.primalFeas)
      res.status = Status::IterationLimit;
  } else {
    res.status = fin_.status;
  }
  // Unproven LP but presolve saw a row whose activity bounds exclude its range:
  // try that one-row Farkas certificate (rigorously checked on the original model).
  if (opt.certify && !isQp && pre.certificateRow() >= 0 &&
      !(res.status == Status::Optimal || res.status == Status::Infeasible || res.status == Status::Unbounded)) {
    SolutionClaim fc;
    fc.status = Status::Infeasible;
    fc.farkas.assign(input.numRows(), 0.0);
    fc.farkas[pre.certificateRow()] = pre.certificateSign();
    Certificate c2 = certify(input, fc, opt.certTol);
    if (c2.accepted) {
      res.certificate = c2;
      res.status = Status::Infeasible;
      fin_.farkas = fc.farkas;
      fin_.x.clear();
      tel["presolve"]["farkas_certificate_row"] = pre.certificateRow();
    }
  }
  timing["certify"] = ct.seconds();

  // ---- Results in model sense ----
  res.x = fin_.x;
  if (!res.x.empty()) {
    res.objective = input.objective(res.x);
    res.rowActivity = input.rowActivity(res.x);
  }
  res.rowDual = fin_.y;
  for (double& v : res.rowDual) v *= sense;
  res.reducedCost = fin_.d;
  for (double& v : res.reducedCost) v *= sense;
  res.colStatus = fin_.cs;
  res.rowStatus = fin_.rs;
  res.farkas = fin_.farkas;
  res.primalRay = fin_.ray;
  if (isMip) {
    res.bestBound = sense * out.bound;
  } else if (opt.certify && std::isfinite(res.certificate.safeDualBound)) {
    res.bestBound = res.certificate.safeDualBound;
  }
  if (!res.x.empty() && std::isfinite(res.bestBound))
    res.relativeGap = std::fabs(res.objective - res.bestBound) / std::max(1.0, std::fabs(res.objective));
  res.seconds = total.seconds();
  timing["total"] = res.seconds;
  tel["timing"] = timing;
  return res;
}

Json resultToJson(const Model& model, const SolveResult& r, bool includeVectors) {
  Json j = Json::object();
  j["solver"] = "PRAMANA 1.0";
  j["model"] = model.name;
  j["status"] = statusName(r.status);
  j["engine_status"] = statusName(r.engineStatus);
  j["engine"] = r.engine;
  if (!r.x.empty()) j["objective"] = r.objective;
  if (std::isfinite(r.bestBound)) j["best_bound"] = r.bestBound;
  if (std::isfinite(r.relativeGap)) j["relative_gap"] = r.relativeGap;
  j["seconds"] = r.seconds;
  j["certificate"] = r.certificate.toJson();
  j["telemetry"] = r.telemetry;
  if (includeVectors) {
    Model m = model;
    m.ensureNames();
    auto named = [&](const std::vector<double>& v, const std::vector<std::string>& names) {
      Json o = Json::object();
      for (size_t k = 0; k < v.size() && k < names.size(); ++k) o[names[k]] = v[k];
      return o;
    };
    if (!r.x.empty()) j["x"] = named(r.x, m.colNames);
    if (!r.rowDual.empty()) j["row_duals"] = named(r.rowDual, m.rowNames);
    if (!r.reducedCost.empty()) j["reduced_costs"] = named(r.reducedCost, m.colNames);
    if (!r.farkas.empty()) j["farkas_ray"] = named(r.farkas, m.rowNames);
    if (!r.primalRay.empty()) j["primal_ray"] = named(r.primalRay, m.colNames);
    if (!r.colStatus.empty()) {
      Json b = Json::object();
      for (size_t k = 0; k < r.colStatus.size(); ++k) b[m.colNames[k]] = basisStatusName(r.colStatus[k]);
      j["column_basis"] = b;
      Json rb = Json::object();
      for (size_t k = 0; k < r.rowStatus.size(); ++k) rb[m.rowNames[k]] = basisStatusName(r.rowStatus[k]);
      j["row_basis"] = rb;
    }
  }
  return j;
}

}  // namespace pramana
