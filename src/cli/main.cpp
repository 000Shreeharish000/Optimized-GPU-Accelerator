// PRAMANA command line interface.
//
//   pramana solve <model.mps[.gz]|.qps> [options]     (default command)
//   pramana verify <model> <result.json>               re-certify a result file
//   pramana parametric <model> --col NAME|--row NAME --kind cost|rhs|lower|upper --from A --to B
//   pramana family <model> --col NAME --kind cost --from A --to B --cases K   (batched GPU vs warm simplex)
//   pramana analyze <model>                            structure report + router prediction (no solve)
//   pramana info                                       build, dependency and GPU report
//   pramana calibrate [--out router_model.json]        SpMV bandwidth micro-benchmark
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "api/solver.h"
#include "cert/certifier.h"
#include "io/mps.h"
#include "parametric/parametric.h"
#include "router/router.h"
#include "pdhg/pdhg.h"
#include "util/json.h"
#include "util/log.h"
#include "util/parallel.h"
#include "util/timer.h"

using namespace pramana;

namespace {

void usage() {
  std::printf(
      "PRAMANA 1.0 - certified LP/MILP/QP solver (from-scratch core, CPU+GPU)\n\n"
      "usage: pramana [solve] <model.mps[.gz]|model.qps> [options]\n"
      "       pramana verify <model> <result.json>\n"
      "       pramana parametric <model> (--col NAME | --row NAME) --kind cost|rhs|lower|upper --from A --to B [--json F]\n"
      "       pramana family <model> (--col NAME | --row NAME) --kind cost|rhs|lower|upper --from A --to B --cases K [--json F]\n"
      "       pramana info\n"
      "       pramana calibrate\n\n"
      "solve options:\n"
      "  --algo auto|dual|ipm|pdhg|pdhg-cpu|pdhg-gpu|race   engine (auto = measured router)\n"
      "  --time SEC            time limit            --threads N    CPU threads\n"
      "  --no-presolve --no-scale --no-certify --no-crossover\n"
      "  --tol T               simplex primal/dual feasibility tolerance (scaled)\n"
      "  --pdhg-tol T  --ipm-tol T  --pdhg-iters N\n"
      "  --gap G  --nodes N  --no-cuts  --no-heuristics  --branching reliability|pseudocost|mostfrac|strong\n"
      "  --mip-threads N       MILP: N concurrent diversified branch-and-cut searches (first proof wins)\n"
      "  --json FILE           write result/certificate/telemetry JSON (--vectors to include x, y, basis)\n"
      "  --cert FILE           write the certificate JSON only\n"
      "  --router-model FILE   calibrated router weights\n"
      "  --log 0..3            verbosity (default 1)\n"
      "  --seed N\n");
}

struct Args {
  std::vector<std::string> pos;
  std::map<std::string, std::string> kv;
  bool has(const std::string& k) const { return kv.count(k) > 0; }
  std::string get(const std::string& k, const std::string& d = "") const {
    auto it = kv.find(k);
    return it == kv.end() ? d : it->second;
  }
  double num(const std::string& k, double d) const { return has(k) ? std::stod(get(k)) : d; }
};

Args parse(int argc, char** argv, int from) {
  static const char* flags[] = {"--no-presolve", "--no-scale", "--no-certify", "--no-crossover", "--no-cuts",
                                "--no-heuristics", "--vectors", "--quiet", "--no-gpu"};
  Args a;
  for (int i = from; i < argc; ++i) {
    std::string s = argv[i];
    if (s.rfind("--", 0) == 0) {
      bool isFlag = false;
      for (auto* f : flags)
        if (s == f) isFlag = true;
      if (isFlag || i + 1 >= argc) a.kv[s] = "1";
      else a.kv[s] = argv[++i];
    } else {
      a.pos.push_back(s);
    }
  }
  return a;
}

Model loadModel(const std::string& path) {
  Timer t;
  Model m = readMps(path);
  PLOG_INFO("read %s: %s (%.2fs)", path.c_str(), m.summary().c_str(), t.seconds());
  return m;
}

int cmdSolve(const Args& a) {
  if (a.pos.empty()) {
    usage();
    return 1;
  }
  Model model = loadModel(a.pos[0]);
  SolverOptions o;
  o.algorithm = a.get("--algo", "auto");
  o.timeLimit = a.num("--time", 0);
  o.threads = static_cast<int>(a.num("--threads", 0));
  o.presolve = !a.has("--no-presolve");
  o.scale = !a.has("--no-scale");
  o.certify = !a.has("--no-certify");
  o.pdhgCrossover = !a.has("--no-crossover");
  o.allowGpu = !a.has("--no-gpu");
  if (a.has("--tol")) o.primalTol = o.dualTol = a.num("--tol", 1e-7);
  o.pdhgTol = a.num("--pdhg-tol", o.pdhgTol);
  o.ipmTol = a.num("--ipm-tol", o.ipmTol);
  o.pdhgMaxIter = static_cast<long long>(a.num("--pdhg-iters", static_cast<double>(o.pdhgMaxIter)));
  o.mipRelGap = a.num("--gap", o.mipRelGap);
  o.nodeLimit = static_cast<long long>(a.num("--nodes", -1));
  o.cuts = !a.has("--no-cuts");
  o.heuristics = !a.has("--no-heuristics");
  o.branching = a.get("--branching", o.branching);
  o.mipThreads = static_cast<int>(a.num("--mip-threads", 1));
  o.routerModelPath = a.get("--router-model", "");
  o.logLevel = static_cast<int>(a.num("--log", 1));
  o.seed = static_cast<uint64_t>(a.num("--seed", 12345));
  SolveResult r = solve(model, o);
  const Certificate& c = r.certificate;
  std::printf("\n");
  std::printf("Model            %s\n", model.summary().c_str());
  std::printf("Engine           %s\n", r.engine.c_str());
  std::printf("Status           %s%s\n", statusName(r.status),
              r.status != r.engineStatus ? formatString("   (engine claimed %s)", statusName(r.engineStatus)).c_str() : "");
  if (!r.x.empty()) std::printf("Objective        %.12g\n", r.objective);
  if (std::isfinite(r.bestBound)) std::printf("Certified bound  %.12g\n", r.bestBound);
  if (std::isfinite(r.relativeGap)) std::printf("Relative gap     %.3e\n", r.relativeGap);
  std::printf("Time             %.3fs\n", r.seconds);
  if (o.certify) {
    std::printf("Certificate      %s\n", c.accepted ? "ACCEPTED" : "REJECTED");
    for (auto& ch : c.checks)
      std::printf("  %-30s %-4s %11.3e  (tol %.1e)\n", ch.name.c_str(), verdictName(ch.verdict), ch.value, ch.tolerance);
  }
  if (a.has("--json")) {
    writeTextFile(a.get("--json"), resultToJson(model, r, a.has("--vectors")).dump(1));
    std::printf("wrote %s\n", a.get("--json").c_str());
  }
  if (a.has("--cert")) writeTextFile(a.get("--cert"), c.toJson().dump(1));
  bool proven = r.status == Status::Optimal || r.status == Status::Infeasible || r.status == Status::Unbounded;
  return proven ? 0 : (r.status == Status::NumericalFailure ? 3 : 2);
}

int cmdVerify(const Args& a) {
  if (a.pos.size() < 2) {
    usage();
    return 1;
  }
  Model model = loadModel(a.pos[0]);
  Json j = Json::parse(readTextFile(a.pos[1]));
  model.ensureNames();
  SolutionClaim claim;
  std::string st = j.at("status").str();
  if (st == "OPTIMAL") claim.status = Status::Optimal;
  else if (st == "INFEASIBLE") claim.status = Status::Infeasible;
  else if (st == "UNBOUNDED") claim.status = Status::Unbounded;
  else if (st == "TIME_LIMIT") claim.status = Status::TimeLimit;
  else claim.status = Status::NotSolved;
  auto vec = [&](const char* key, const std::vector<std::string>& names) {
    std::vector<double> v;
    if (!j.has(key)) return v;
    const Json& o = j.at(key);
    v.assign(names.size(), 0.0);
    for (size_t k = 0; k < names.size(); ++k) v[k] = o.at(names[k]).num();
    return v;
  };
  claim.x = vec("x", model.colNames);
  const double sense = model.sense == ObjSense::Maximize ? -1.0 : 1.0;
  claim.rowDual = vec("row_duals", model.rowNames);
  for (double& v : claim.rowDual) v *= sense;  // stored in model sense
  claim.farkas = vec("farkas_ray", model.rowNames);
  claim.primalRay = vec("primal_ray", model.colNames);
  if (j.has("best_bound")) {
    claim.bestBound = j.at("best_bound").num();
    claim.boundIsSafe = j.at("telemetry").at("mip").at("bound_is_safe").boolean(false);
  }
  Certificate c = certify(model, claim);
  std::printf("%s\n", c.summary().c_str());
  return c.accepted ? 0 : 4;
}

// Structure / numerics report of a model without solving it (used by the terminal UI).
int cmdAnalyze(const Args& a) {
  if (a.pos.empty()) {
    usage();
    return 1;
  }
  Model model = loadModel(a.pos[0]);
  const int m = model.numRows(), n = model.numCols();
  auto fin = [](double v) { return std::isfinite(v); };
  Json j = Json::object();
  j["name"] = model.name;
  j["sense"] = model.sense == ObjSense::Maximize ? "max" : "min";
  j["rows"] = m;
  j["cols"] = n;
  j["nnz"] = model.A.nnz();
  j["integers"] = model.numIntegers();
  j["q_nnz"] = model.Q.nnz();
  j["class"] = model.Q.nnz() > 0 ? (model.numIntegers() > 0 ? "MIQP" : "QP") : (model.numIntegers() > 0 ? "MILP" : "LP");
  int binaries = 0;
  for (int k = 0; k < n; ++k)
    if (model.colType[k] == VarType::Integer && model.colLower[k] == 0 && model.colUpper[k] == 1) ++binaries;
  j["binaries"] = binaries;
  int eq = 0, le = 0, ge = 0, rng = 0, freeRows = 0;
  for (int i = 0; i < m; ++i) {
    bool l = fin(model.rowLower[i]), u = fin(model.rowUpper[i]);
    if (l && u) (model.rowLower[i] == model.rowUpper[i] ? eq : rng)++;
    else if (u) ++le;
    else if (l) ++ge;
    else ++freeRows;
  }
  Json rt = Json::object();
  rt["equality"] = eq; rt["less_equal"] = le; rt["greater_equal"] = ge; rt["range"] = rng; rt["free"] = freeRows;
  j["row_types"] = rt;
  int fr = 0, lo = 0, up = 0, box = 0, fx = 0;
  for (int k = 0; k < n; ++k) {
    bool l = fin(model.colLower[k]), u = fin(model.colUpper[k]);
    if (l && u) (model.colLower[k] == model.colUpper[k] ? fx : box)++;
    else if (l) ++lo;
    else if (u) ++up;
    else ++fr;
  }
  Json cb = Json::object();
  cb["free"] = fr; cb["lower_only"] = lo; cb["upper_only"] = up; cb["boxed"] = box; cb["fixed"] = fx;
  j["col_bounds"] = cb;
  double amin = 1e300, amax = 0;
  std::vector<int> rowCount(m, 0);
  int emptyCols = 0, singletonCols = 0, maxCol = 0;
  for (int k = 0; k < n; ++k) {
    int c = 0;
    for (int e = model.A.start[k]; e < model.A.start[k + 1]; ++e) {
      double v = std::fabs(model.A.value[e]);
      if (v == 0) continue;
      amin = std::min(amin, v);
      amax = std::max(amax, v);
      ++rowCount[model.A.index[e]];
      ++c;
    }
    emptyCols += c == 0;
    singletonCols += c == 1;
    maxCol = std::max(maxCol, c);
  }
  int emptyRows = 0, singletonRows = 0, maxRow = 0;
  for (int c : rowCount) {
    emptyRows += c == 0;
    singletonRows += c == 1;
    maxRow = std::max(maxRow, c);
  }
  auto range = [&](const std::vector<double>& v) {
    double lo2 = 1e300, hi2 = 0;
    for (double x : v)
      if (fin(x) && x != 0) {
        lo2 = std::min(lo2, std::fabs(x));
        hi2 = std::max(hi2, std::fabs(x));
      }
    Json r = Json::object();
    r["min_abs"] = hi2 > 0 ? lo2 : 0.0;
    r["max_abs"] = hi2;
    return r;
  };
  Json num = Json::object();
  num["matrix_min_abs"] = amax > 0 ? amin : 0.0;
  num["matrix_max_abs"] = amax;
  num["matrix_range"] = amax > 0 ? amax / amin : 1.0;
  num["cost"] = range(model.colCost);
  std::vector<double> rhs(model.rowLower);
  rhs.insert(rhs.end(), model.rowUpper.begin(), model.rowUpper.end());
  num["rhs"] = range(rhs);
  std::vector<double> bnds(model.colLower);
  bnds.insert(bnds.end(), model.colUpper.begin(), model.colUpper.end());
  num["bounds"] = range(bnds);
  j["numerics"] = num;
  Json st = Json::object();
  st["density"] = (m > 0 && n > 0) ? static_cast<double>(model.A.nnz()) / (static_cast<double>(m) * n) : 0.0;
  st["max_row_nnz"] = maxRow;
  st["max_col_nnz"] = maxCol;
  st["empty_rows"] = emptyRows;
  st["empty_cols"] = emptyCols;
  st["singleton_rows"] = singletonRows;
  st["singleton_cols"] = singletonCols;
  j["structure"] = st;
  Router router;
  if (a.has("--router-model")) router.load(a.get("--router-model"));
  const bool lp = model.numIntegers() == 0 && model.Q.nnz() == 0;
  const bool gpu = lp && !a.has("--no-gpu") && gpuAvailable();
  RouterDecision d = router.decide(model.toMinimization(), gpu, 1e-7);
  Json r = d.toJson();
  r["gpu_available"] = gpu;
  if (!lp) r["engine"] = model.Q.nnz() > 0 ? "ipm-qp" : "bnb";
  j["router"] = r;
  std::printf("%s\n", j.dump(1).c_str());
  return 0;
}

int cmdInfo() {
  std::printf("PRAMANA 1.0\n");
  std::printf("  solver core      : own code (C++20), no solver libraries linked\n");
  std::printf("  threads          : %d\n", ThreadPool::instance().threads());
  Json g = gpuInfoJson();
  std::printf("  gpu              : %s\n", g.at("available").boolean() ? g.at("name").str().c_str() : ("unavailable - " + g.at("reason").str()).c_str());
  if (g.at("available").boolean()) {
    std::printf("  compute cap.     : %s, %d SMs, %.1f GB, peak %.0f GB/s\n", g.at("compute_capability").str().c_str(),
                static_cast<int>(g.at("multiprocessors").num()), g.at("memory_bytes").num() / 1e9,
                g.at("peak_bandwidth_GBs").num());
    std::printf("  kernels          : own CUDA C, %s\n",
                g.at("ptx_from_cache").boolean() ? "PTX from cache" : ("compiled by NVRTC " + g.at("nvrtc").str()).c_str());
  }
  std::printf("%s\n", g.dump(1).c_str());
  return 0;
}

int cmdParametric(const Args& a, bool family) {
  if (a.pos.empty()) {
    usage();
    return 1;
  }
  Model model = loadModel(a.pos[0]);
  model.ensureNames();
  ParametricSpec spec;
  spec.kind = a.get("--kind", "cost");
  std::string col = a.get("--col", ""), row = a.get("--row", "");
  spec.index = -1;
  if (!col.empty()) {
    for (int j = 0; j < model.numCols(); ++j)
      if (model.colNames[j] == col) spec.index = j;
    spec.isRow = false;
  } else {
    for (int i = 0; i < model.numRows(); ++i)
      if (model.rowNames[i] == row) spec.index = i;
    spec.isRow = true;
  }
  if (spec.index < 0) {
    std::printf("unknown %s '%s'\n", col.empty() ? "row" : "column", col.empty() ? row.c_str() : col.c_str());
    return 1;
  }
  spec.from = a.num("--from", 0);
  spec.to = a.num("--to", 1);
  if (!family) {
    ParametricResult r = parametricAnalysis(model, spec);
    std::printf("%s\n", r.summary(model).c_str());
    if (a.has("--json")) writeTextFile(a.get("--json"), r.toJson(model).dump(1));
    return r.ok ? 0 : 2;
  }
  int cases = static_cast<int>(a.num("--cases", 64));
  FamilyComparison f = compareFamilyStrategies(model, spec, cases, !a.has("--no-gpu"));
  std::printf("%s\n", f.summary().c_str());
  if (a.has("--json")) writeTextFile(a.get("--json"), f.toJson().dump(1));
  return 0;
}

int cmdCalibrate(const Args& a) {
  Json j = calibrateBandwidth();
  std::printf("%s\n", j.dump(1).c_str());
  if (a.has("--out")) writeTextFile(a.get("--out"), j.dump(1));
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    usage();
    return 1;
  }
  std::string cmd = argv[1];
  try {
    Args a;
    int rc;
    if (cmd == "-h" || cmd == "--help" || cmd == "help") {
      usage();
      return 0;
    }
    if (cmd == "info") return cmdInfo();
    if (cmd == "verify" || cmd == "analyze" || cmd == "parametric" || cmd == "family" || cmd == "calibrate" || cmd == "solve") {
      a = parse(argc, argv, 2);
    } else {
      a = parse(argc, argv, 1);
      cmd = "solve";
    }
    Log::setLevel(a.has("--quiet") ? LogLevel::Error : static_cast<LogLevel>(std::min(4.0, 1.0 + a.num("--log", 1))));
    if (cmd == "verify") rc = cmdVerify(a);
    else if (cmd == "analyze") rc = cmdAnalyze(a);
    else if (cmd == "parametric") rc = cmdParametric(a, false);
    else if (cmd == "family") rc = cmdParametric(a, true);
    else if (cmd == "calibrate") rc = cmdCalibrate(a);
    else rc = cmdSolve(a);
    return rc;
  } catch (std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 10;
  }
}
