// C API implementation (thin layer over api/solver.h).
#include "pramana/pramana.h"

#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>

#include "api/solver.h"
#include "io/mps.h"
#include "parametric/parametric.h"
#include "pdhg/pdhg.h"
#include "util/json.h"
#include "util/log.h"

using namespace pramana;

namespace {
thread_local std::string gLastError;

char* dupString(const std::string& s) {
  char* p = static_cast<char*>(std::malloc(s.size() + 1));
  if (p) std::memcpy(p, s.c_str(), s.size() + 1);
  return p;
}
}  // namespace

struct pramana_model {
  Model model;
  // Buffered rows (appended in bulk before solve to keep row insertion O(nnz)).
  std::vector<double> rl, ru;
  std::vector<int> starts{0}, cols;
  std::vector<double> vals;
  std::vector<std::string> rowNames;
  std::vector<int> qi, qj;
  std::vector<double> qv;
  void flush() {
    if (!rl.empty()) {
      int first = model.numRows();
      model.addRows(rl, ru, starts, cols, vals);
      for (size_t k = 0; k < rowNames.size(); ++k)
        if (!rowNames[k].empty()) model.rowNames[first + k] = rowNames[k];
      rl.clear();
      ru.clear();
      starts.assign(1, 0);
      cols.clear();
      vals.clear();
      rowNames.clear();
    }
    if (!qv.empty()) {
      int n = model.numCols();
      std::vector<int> ri(qi), ci(qj);
      std::vector<double> v(qv);
      for (int j = 0; j < model.Q.numCols; ++j)
        for (int k = model.Q.start[j]; k < model.Q.start[j + 1]; ++k) {
          ri.push_back(model.Q.index[k]);
          ci.push_back(j);
          v.push_back(model.Q.value[k]);
        }
      model.Q = SparseMatrix::fromTriplets(n, n, ri, ci, v);
      qi.clear();
      qj.clear();
      qv.clear();
    }
  }
};

struct pramana_result {
  SolveResult result;
  std::string json;
  const Model* model = nullptr;
  Model modelCopy;
};

extern "C" {

PRAMANA_API const char* pramana_version(void) { return "PRAMANA 1.0.0"; }
PRAMANA_API const char* pramana_last_error(void) { return gLastError.c_str(); }

PRAMANA_API pramana_model* pramana_model_create(const char* name, int maximize) {
  auto* m = new pramana_model();
  m->model.name = name ? name : "model";
  m->model.sense = maximize ? ObjSense::Maximize : ObjSense::Minimize;
  m->model.A = SparseMatrix(0, 0);
  return m;
}

PRAMANA_API pramana_model* pramana_model_read(const char* path) {
  try {
    auto* m = new pramana_model();
    m->model = readMps(path);
    return m;
  } catch (std::exception& e) {
    gLastError = e.what();
    return nullptr;
  }
}

PRAMANA_API void pramana_model_free(pramana_model* m) { delete m; }

PRAMANA_API int pramana_add_col(pramana_model* m, double cost, double lower, double upper, int integer,
                                const char* name) {
  m->flush();
  return m->model.addColumn(cost, lower <= -1e20 ? -kInf : lower, upper >= 1e20 ? kInf : upper, {}, {},
                            integer ? VarType::Integer : VarType::Continuous, name ? name : "");
}

PRAMANA_API int pramana_add_row(pramana_model* m, double lower, double upper, int nnz, const int* cols,
                                const double* vals, const char* name) {
  for (int k = 0; k < nnz; ++k) {
    if (cols[k] < 0 || cols[k] >= m->model.numCols()) {
      gLastError = "pramana_add_row: column index out of range";
      return -1;
    }
    m->cols.push_back(cols[k]);
    m->vals.push_back(vals[k]);
  }
  m->starts.push_back(static_cast<int>(m->cols.size()));
  m->rl.push_back(lower <= -1e20 ? -kInf : lower);
  m->ru.push_back(upper >= 1e20 ? kInf : upper);
  m->rowNames.push_back(name ? name : "");
  return m->model.numRows() + static_cast<int>(m->rl.size()) - 1;
}

PRAMANA_API int pramana_add_q(pramana_model* m, int i, int j, double v) {
  m->qi.push_back(i);
  m->qj.push_back(j);
  m->qv.push_back(v);
  if (i != j) {
    m->qi.push_back(j);
    m->qj.push_back(i);
    m->qv.push_back(v);
  }
  return 0;
}

PRAMANA_API void pramana_set_offset(pramana_model* m, double offset) { m->model.objOffset = offset; }
PRAMANA_API int pramana_num_cols(pramana_model* m) { return m->model.numCols(); }
PRAMANA_API int pramana_num_rows(pramana_model* m) {
  return m->model.numRows() + static_cast<int>(m->rl.size());
}

PRAMANA_API int pramana_write_mps(pramana_model* m, const char* path) {
  try {
    m->flush();
    writeMps(m->model, path);
    return 0;
  } catch (std::exception& e) {
    gLastError = e.what();
    return -1;
  }
}

}  // extern "C"

// C++ helpers (C++ return types / may throw): kept outside the extern "C" block.
static SolverOptions optionsFromJson(const char* text) {
  SolverOptions o;
  o.logLevel = 0;
  if (!text || !*text) return o;
  Json j = Json::parse(text);
  if (j.has("algorithm")) o.algorithm = j.at("algorithm").str();
  if (j.has("time_limit")) o.timeLimit = j.at("time_limit").num();
  if (j.has("presolve")) o.presolve = j.at("presolve").boolean(true);
  if (j.has("certify")) o.certify = j.at("certify").boolean(true);
  if (j.has("scale")) o.scale = j.at("scale").boolean(true);
  if (j.has("threads")) o.threads = static_cast<int>(j.at("threads").num());
  if (j.has("log_level")) o.logLevel = static_cast<int>(j.at("log_level").num());
  if (j.has("seed")) o.seed = static_cast<uint64_t>(j.at("seed").num());
  if (j.has("tolerance")) o.primalTol = o.dualTol = j.at("tolerance").num();
  if (j.has("pdhg_tolerance")) o.pdhgTol = j.at("pdhg_tolerance").num();
  if (j.has("pdhg_crossover")) o.pdhgCrossover = j.at("pdhg_crossover").boolean(true);
  if (j.has("ipm_tolerance")) o.ipmTol = j.at("ipm_tolerance").num();
  if (j.has("mip_gap")) o.mipRelGap = j.at("mip_gap").num();
  if (j.has("node_limit")) o.nodeLimit = static_cast<long long>(j.at("node_limit").num());
  if (j.has("cuts")) o.cuts = j.at("cuts").boolean(true);
  if (j.has("heuristics")) o.heuristics = j.at("heuristics").boolean(true);
  if (j.has("branching")) o.branching = j.at("branching").str();
  if (j.has("mip_threads")) o.mipThreads = static_cast<int>(j.at("mip_threads").num());
  if (j.has("allow_gpu")) o.allowGpu = j.at("allow_gpu").boolean(true);
  if (j.has("router_model")) o.routerModelPath = j.at("router_model").str();
  return o;
}

static ParametricSpec specFromJson(pramana_model* m, const char* text);

extern "C" {

PRAMANA_API pramana_result* pramana_solve(pramana_model* m, const char* options_json) {
  try {
    m->flush();
    SolverOptions o = optionsFromJson(options_json);
    Log::setLevel(o.logLevel <= 0 ? LogLevel::Error : static_cast<LogLevel>(std::min(4, 1 + o.logLevel)));
    auto* r = new pramana_result();
    r->modelCopy = m->model;
    r->model = &r->modelCopy;
    r->result = solve(m->model, o);
    return r;
  } catch (std::exception& e) {
    gLastError = e.what();
    return nullptr;
  }
}

PRAMANA_API int pramana_result_status(const pramana_result* r) { return static_cast<int>(r->result.status); }
PRAMANA_API const char* pramana_result_status_name(const pramana_result* r) { return statusName(r->result.status); }
PRAMANA_API double pramana_result_objective(const pramana_result* r) { return r->result.objective; }
PRAMANA_API double pramana_result_bound(const pramana_result* r) { return r->result.bestBound; }
PRAMANA_API int pramana_result_certified(const pramana_result* r) { return r->result.certificate.accepted ? 1 : 0; }

static int copyOut(const std::vector<double>& v, double* out, int len) {
  if (out)
    for (int k = 0; k < len && k < static_cast<int>(v.size()); ++k) out[k] = v[k];
  return static_cast<int>(v.size());
}
PRAMANA_API int pramana_result_x(const pramana_result* r, double* out, int len) { return copyOut(r->result.x, out, len); }
PRAMANA_API int pramana_result_row_duals(const pramana_result* r, double* out, int len) {
  return copyOut(r->result.rowDual, out, len);
}
PRAMANA_API int pramana_result_reduced_costs(const pramana_result* r, double* out, int len) {
  return copyOut(r->result.reducedCost, out, len);
}
PRAMANA_API int pramana_result_row_activity(const pramana_result* r, double* out, int len) {
  return copyOut(r->result.rowActivity, out, len);
}
PRAMANA_API const char* pramana_result_json(pramana_result* r, int include_vectors) {
  r->json = resultToJson(*r->model, r->result, include_vectors != 0).dump(1);
  return r->json.c_str();
}
PRAMANA_API void pramana_result_free(pramana_result* r) { delete r; }

}  // extern "C"

static ParametricSpec specFromJson(pramana_model* m, const char* text) {
  Json j = Json::parse(text ? text : "{}");
  ParametricSpec s;
  s.kind = j.has("kind") ? j.at("kind").str() : "cost";
  s.from = j.at("from").num(0);
  s.to = j.at("to").num(1);
  m->model.ensureNames();
  if (j.has("col")) {
    s.isRow = false;
    for (int k = 0; k < m->model.numCols(); ++k)
      if (m->model.colNames[k] == j.at("col").str()) s.index = k;
  } else if (j.has("row")) {
    s.isRow = true;
    for (int k = 0; k < m->model.numRows(); ++k)
      if (m->model.rowNames[k] == j.at("row").str()) s.index = k;
  }
  if (s.index < 0) throw PramanaError("parametric: unknown column/row");
  return s;
}

extern "C" {

PRAMANA_API char* pramana_parametric(pramana_model* m, const char* spec_json) {
  try {
    m->flush();
    ParametricSpec s = specFromJson(m, spec_json);
    return dupString(parametricAnalysis(m->model, s).toJson(m->model).dump(1));
  } catch (std::exception& e) {
    gLastError = e.what();
    return nullptr;
  }
}

PRAMANA_API char* pramana_family(pramana_model* m, const char* spec_json) {
  try {
    m->flush();
    ParametricSpec s = specFromJson(m, spec_json);
    Json j = Json::parse(spec_json ? spec_json : "{}");
    int cases = static_cast<int>(j.at("cases").num(32));
    bool gpu = j.has("allow_gpu") ? j.at("allow_gpu").boolean(true) : true;
    return dupString(compareFamilyStrategies(m->model, s, cases, gpu).toJson().dump(1));
  } catch (std::exception& e) {
    gLastError = e.what();
    return nullptr;
  }
}

PRAMANA_API char* pramana_gpu_info(void) { return dupString(gpuInfoJson().dump(1)); }
PRAMANA_API void pramana_free_string(char* s) { std::free(s); }

}  // extern "C"
