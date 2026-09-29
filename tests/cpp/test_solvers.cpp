// Tests: simplex, presolve round trip, IPM, PDHG (CPU & GPU agreement), MIP vs
// brute force, parametric analysis vs re-solves, determinism, router.
#include <cmath>
#include <string>

#include "api/solver.h"
#include "cert/certifier.h"
#include "io/mps.h"
#include "ipm/ipm.h"
#include "lp/simplex.h"
#include "mip/bnb.h"
#include "parametric/parametric.h"
#include "pdhg/pdhg.h"
#include "presolve/presolve.h"
#include "router/router.h"
#include "testing.h"
#include "util/rng.h"

using namespace pramana;

static Model netlib(const std::string& n) { return readMps(dataPath("data/netlib/" + n + ".mps.gz")); }

PTEST(simplex_netlib_known_optima) {
  struct Ref {
    const char* name;
    double obj;
  } refs[] = {{"afiro", -464.753142857}, {"sc50a", -64.5750770585}, {"adlittle", 225494.963162},
              {"blend", -30.8121498458}, {"share2b", -415.732240741}, {"bore3d", 1373.08039421}};
  for (auto& r : refs) {
    Model m = netlib(r.name).toMinimization();
    LpResult res = solveLpSimplex(m);
    EXPECT(res.status == Status::Optimal);
    EXPECT_NEAR(res.objective, r.obj, 1e-9);
    SolutionClaim c{Status::Optimal, res.x, res.rowDual};
    EXPECT(certify(m, c).accepted);
  }
}

PTEST(simplex_infeasible_and_unbounded_certificates) {
  Model inf = readMps(dataPath("data/netlib_infeas/itest6.mps"));
  LpResult r = solveLpSimplex(inf.toMinimization());
  EXPECT(r.status == Status::Infeasible);
  SolutionClaim c;
  c.status = Status::Infeasible;
  c.farkas = r.farkas;
  EXPECT(certify(inf, c).accepted);
  Model unb = readMps(dataPath("tests/data/unbounded.mps"));
  LpResult u = solveLpSimplex(unb.toMinimization());
  EXPECT(u.status == Status::Unbounded);
  SolutionClaim cu;
  cu.status = Status::Unbounded;
  cu.primalRay = u.primalRay;
  cu.x = u.x;
  EXPECT(certify(unb, cu).accepted);
}

PTEST(simplex_warm_start_after_bound_change) {
  Model m = netlib("sc105").toMinimization();
  Simplex s(m);
  EXPECT(s.solve() == Status::Optimal);
  long long cold = s.stats().iterations;
  std::vector<double> x = s.colValues();
  int j = 0;
  while (j < m.numCols() && x[j] < 1e-3) ++j;
  s.setColBounds(j, 0, x[j] * 0.5);  // cut off the current optimum
  long long before = s.stats().iterations;
  EXPECT(s.solve() == Status::Optimal);
  long long warm = s.stats().iterations - before;
  EXPECT(warm < cold);  // reoptimization is cheaper than a cold solve
  Model m2 = m;
  m2.colUpper[j] = x[j] * 0.5;
  LpResult ref = solveLpSimplex(m2);
  EXPECT_NEAR(s.objective(), ref.objective, 1e-9);
}

PTEST(presolve_postsolve_roundtrip_certifies) {
  for (const char* n : {"afiro", "adlittle", "sc205", "scagr7", "share1b", "ship04s", "stocfor1", "bore3d"}) {
    Model orig = netlib(n);
    SolverOptions o;
    o.logLevel = 0;
    o.presolve = true;
    SolveResult a = solve(orig, o);
    o.presolve = false;
    SolveResult b = solve(orig, o);
    EXPECT(a.status == Status::Optimal && b.status == Status::Optimal);
    EXPECT_NEAR(a.objective, b.objective, 1e-9);
    EXPECT(a.certificate.accepted);
  }
}

PTEST(ipm_lp_and_qp) {
  IpmOptions io;
  io.crossover = false;
  IpmResult r = solveIpm(netlib("afiro").toMinimization(), io);
  EXPECT(r.status == Status::Optimal);
  EXPECT_NEAR(r.objective, -464.753142857, 1e-7);
  // min x^2 + y^2 - 4x - 6y  s.t. x + y <= 2  -> x=0.5, y=1.5, obj = -8.5
  Model q;
  q.addColumn(-4, -kInf, kInf, {}, {});
  q.addColumn(-6, -kInf, kInf, {}, {});
  q.addRow(-kInf, 2, {0, 1}, {1, 1});
  q.Q = SparseMatrix::fromTriplets(2, 2, {0, 1}, {0, 1}, {2, 2});
  IpmResult qr = solveIpm(q, io);
  EXPECT(qr.status == Status::Optimal);
  EXPECT_NEAR(qr.x[0], 0.5, 1e-6);
  EXPECT_NEAR(qr.x[1], 1.5, 1e-6);
  EXPECT_NEAR(qr.objective, -8.5, 1e-8);
  SolutionClaim c{Status::Optimal, qr.x, qr.rowDual};
  EXPECT(certify(q, c).accepted);
}

PTEST(pdhg_cpu_matches_optimum_and_gpu_agrees) {
  Model m = netlib("sc50a").toMinimization();
  PdhgOptions po;
  po.tolerance = 1e-7;
  PdhgResult cpu = solvePdhg(m, po);
  EXPECT(cpu.status == Status::Optimal);
  EXPECT_NEAR(cpu.objective, -64.5750770585, 1e-5);
  // The safe dual bound from the approximate PDHG duals is rigorous.
  SafeBoundResult sb = safeDualBound(m.A, m.colCost, m.colLower, m.colUpper, m.rowLower, m.rowUpper, cpu.rowDual);
  EXPECT(sb.bound <= -64.5750770585 + 1e-9);
  if (gpuAvailable()) {
    po.useGpu = true;
    PdhgResult gpu = solvePdhg(m, po);
    EXPECT(gpu.status == Status::Optimal);
    EXPECT_NEAR(gpu.objective, cpu.objective, 1e-5);
    // Same algorithm, same iteration count up to reduction-order effects.
    EXPECT(std::llabs(gpu.iterations - cpu.iterations) <= 4 * po.checkEvery);
  }
}

PTEST(pdhg_batch_equals_individual_solves) {
  Model m = netlib("afiro").toMinimization();
  std::vector<LpCase> cases(4);
  for (int k = 0; k < 4; ++k) {
    cases[k].cost = m.colCost;
    cases[k].cost[1] = -0.4 + 0.1 * k;
  }
  PdhgOptions po;
  po.tolerance = 1e-7;
  std::vector<PdhgResult> rb = solvePdhgBatch(m, cases, po);
  for (int k = 0; k < 4; ++k) {
    Model mk = m;
    mk.colCost = cases[k].cost;
    LpResult ref = solveLpSimplex(mk);
    EXPECT(rb[k].status == Status::Optimal);
    EXPECT_NEAR(rb[k].objective, ref.objective, 1e-5);
  }
}

PTEST(mip_matches_brute_force_knapsacks) {
  Rng rng(11);
  for (int inst = 0; inst < 8; ++inst) {
    const int n = 12;
    Model m;
    m.sense = ObjSense::Maximize;
    std::vector<double> w(n), p(n);
    std::vector<int> idx(n);
    for (int j = 0; j < n; ++j) {
      w[j] = 1 + rng.below(20);
      p[j] = 1 + rng.below(30);
      idx[j] = j;
      m.addColumn(p[j], 0, 1, {}, {}, VarType::Integer);
    }
    double cap = 0;
    for (double t : w) cap += t;
    cap = std::floor(cap / 3);
    m.addRow(-kInf, cap, idx, w);
    m.addRow(-kInf, 5, idx, std::vector<double>(n, 1.0));  // cardinality
    double best = 0;
    for (int mask = 0; mask < (1 << n); ++mask) {
      double ww = 0, pp = 0;
      int card = 0;
      for (int j = 0; j < n; ++j)
        if (mask >> j & 1) {
          ww += w[j];
          pp += p[j];
          ++card;
        }
      if (ww <= cap && card <= 5) best = std::max(best, pp);
    }
    SolverOptions o;
    o.logLevel = 0;
    o.mipRelGap = 0;
    SolveResult r = solve(m, o);
    EXPECT(r.status == Status::Optimal);
    EXPECT_NEAR(r.objective, best, 1e-9);
  }
}

PTEST(mip_miplib_small_optima) {
  struct Ref {
    const char* name;
    double obj;
  } refs[] = {{"p0033", 3089}, {"lseu", 1120}, {"egout", 568.1007}, {"gt2", 21166},
              {"p2756", 3124}, {"p0548", 8691}, {"p0282", 258411}};  // p2756/p0548: coefficient tightening
  for (auto& r : refs) {
    Model m = readMps(dataPath(std::string("data/miplib3/") + r.name + ".mps.gz"));
    SolverOptions o;
    o.logLevel = 0;
    o.mipRelGap = 1e-6;
    SolveResult s = solve(m, o);
    EXPECT(s.status == Status::Optimal);
    EXPECT_NEAR(s.objective, r.obj, 1e-6);
  }
}

PTEST(parametric_matches_resolves) {
  Model m = netlib("afiro");
  m.ensureNames();
  ParametricSpec spec;
  spec.kind = "rhs";
  spec.isRow = true;
  for (int i = 0; i < m.numRows(); ++i)
    if (m.rowNames[i] == "X05") spec.index = i;
  spec.from = 1;
  spec.to = 150;
  ParametricResult pr = parametricAnalysis(m, spec);
  EXPECT(pr.ok);
  for (auto& s : pr.segments) EXPECT(s.status != Status::Optimal || s.certified);
  for (double th : {1.0, 20.0, 54.5, 60.0, 89.0, 100.0, 149.0}) {
    Model mm = withParameter(m, spec, th);
    LpResult r = solveLpSimplex(mm.toMinimization());
    EXPECT_NEAR(pr.valueAt(th), r.objective, 1e-7);
  }
}

PTEST(determinism_identical_runs) {
  Model m = netlib("scfxm1");
  SolverOptions o;
  o.logLevel = 0;
  o.algorithm = "dual";
  SolveResult a = solve(m, o), b = solve(m, o);
  EXPECT(a.objective == b.objective);
  EXPECT(a.x == b.x);
  EXPECT(a.telemetry.at("simplex").at("iterations").num() == b.telemetry.at("simplex").at("iterations").num());
}

PTEST(router_features_and_rules) {
  Router r;
  Model mip = readMps(dataPath("data/miplib3/p0033.mps.gz"));
  EXPECT(r.decide(mip, false, 1e-7).engine == "bnb");
  Model lp = netlib("afiro");
  RouterDecision d = r.decide(lp, false, 1e-7);
  EXPECT(d.engine != "pdhg-gpu");  // no GPU offered
  EXPECT(d.features.nnz == 83);
}
