#include "mip/bnb.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <queue>
#include <thread>

#include "cert/certifier.h"
#include "cert/interval.h"
#include "mip/cuts.h"
#include "mip/propagate.h"
#include "util/log.h"
#include "util/rng.h"

namespace pramana {

namespace {

inline bool fin(double v) { return std::fabs(v) < kInfBoundThreshold; }
constexpr double kIntTol = 1e-6;

struct BoundChange {
  int var;
  double lo, up;
};

using BasisVec = std::vector<int8_t>;

struct Node {
  long long id = 0;
  int depth = 0;
  double bound = -kInf;     // rigorous lower bound (min form)
  double estimate = 0;
  std::vector<BoundChange> changes;
  std::shared_ptr<BasisVec> basis;
  int branchVar = -1;
  int branchDir = 0;        // -1 down, +1 up
  double parentObj = 0;
  double branchFrac = 0;    // distance moved by the branching bound
};

struct HeapCmp {
  bool operator()(const std::shared_ptr<Node>& a, const std::shared_ptr<Node>& b) const {
    if (a->bound != b->bound) return a->bound > b->bound;
    return a->id > b->id;
  }
};

class BranchAndBound {
 public:
  BranchAndBound(const Model& model, const MipOptions& opt, const Deadline* dl)
      : model_(model), lpModel_(model), opt_(opt), dl_(dl), rng_(opt.seed) {
    n_ = model.numCols();
    isInt_.assign(n_, 0);
    for (int j = 0; j < n_; ++j) isInt_[j] = model.colType[j] == VarType::Integer;
    rootLo_ = model.colLower;
    rootUp_ = model.colUpper;
    pcDown_.assign(n_, 0.0);
    pcUp_.assign(n_, 0.0);
    pcDownN_.assign(n_, 0);
    pcUpN_.assign(n_, 0);
    // Objective integrality: integer costs on integer columns, zero on continuous.
    objIntegral_ = true;
    for (int j = 0; j < n_; ++j) {
      double c = model.colCost[j];
      if (c == 0) continue;
      if (!isInt_[j] || c != std::floor(c)) objIntegral_ = false;
    }
    if (model.objOffset != std::floor(model.objOffset)) objIntegral_ = false;
  }

  MipResult run();

 private:
  // ---- LP helpers ----
  void applyBounds(const std::vector<double>& lo, const std::vector<double>& up) {
    for (int j = 0; j < n_; ++j)
      if (lo[j] != curLo_[j] || up[j] != curUp_[j]) {
        lp_->setColBounds(j, lo[j], up[j]);
        curLo_[j] = lo[j];
        curUp_[j] = up[j];
      }
  }
  Status solveLp(long long iterLimit = -1) {
    long long before = lp_->stats().iterations;
    lp_->options().maxIterations = iterLimit > 0 ? before + iterLimit : (1LL << 60);
    Status st = lp_->solve(dl_);
    stats_.lpIterations += lp_->stats().iterations - before;
    lp_->options().maxIterations = 1LL << 60;
    return st;
  }
  std::shared_ptr<BasisVec> saveBasis() {
    std::vector<BasisStatus> cs, rs;
    lp_->getBasis(cs, rs);
    auto b = std::make_shared<BasisVec>(cs.size() + rs.size());
    for (size_t k = 0; k < cs.size(); ++k) (*b)[k] = static_cast<int8_t>(cs[k]);
    for (size_t k = 0; k < rs.size(); ++k) (*b)[cs.size() + k] = static_cast<int8_t>(rs[k]);
    return b;
  }
  void loadBasis(const BasisVec& b) {
    int m = lp_->numRows();
    if (static_cast<int>(b.size()) != n_ + m) return;
    std::vector<BasisStatus> cs(n_), rs(m);
    for (int k = 0; k < n_; ++k) cs[k] = static_cast<BasisStatus>(b[k]);
    for (int k = 0; k < m; ++k) rs[k] = static_cast<BasisStatus>(b[n_ + k]);
    lp_->setBasis(cs, rs);
  }
  // Rigorous lower bound for the current LP (node bounds lo/up), min form incl. offset.
  double safeBound(const std::vector<double>& lo, const std::vector<double>& up, double lpObj, bool* safe) {
    std::vector<double> y = lp_->rowDuals();
    SafeBoundResult sb =
        safeDualBound(lpModel_.A, lpModel_.colCost, lo, up, lpModel_.rowLower, lpModel_.rowUpper, y);
    if (std::isinf(sb.bound)) {
      *safe = false;
      return lpObj - 1e-9 * (1 + std::fabs(lpObj));
    }
    *safe = true;
    double b = addDown(sb.bound, lpModel_.objOffset);
    if (b < lpObj - 1e-9 * (1 + std::fabs(lpObj))) {
      stats_.safeBoundCorrections++;
      stats_.maxSafeBoundCorrection = std::max(stats_.maxSafeBoundCorrection, lpObj - b);
    }
    return std::min(b, lpObj);
  }
  double roundBound(double b) const {
    if (objIntegral_ && fin(b)) return std::ceil(b - 1e-6);
    return b;
  }
  double cutoff() const {
    if (!fin(incumbentObj_)) return kInf;
    double c = incumbentObj_ - std::max(opt_.absGap, 1e-9 * std::fabs(incumbentObj_));
    if (objIntegral_) c = incumbentObj_ - 1.0 + 1e-6;
    return c;
  }

  // ---- incumbent handling ----
  bool rowsFeasible(const std::vector<double>& x, double tol) const {
    std::vector<double> act(model_.numRows(), 0.0);
    model_.A.multiply(x.data(), act.data());
    for (int i = 0; i < model_.numRows(); ++i) {
      double l = model_.rowLower[i], u = model_.rowUpper[i];
      if (fin(l) && act[i] < l - tol * (1 + std::fabs(l))) return false;
      if (fin(u) && act[i] > u + tol * (1 + std::fabs(u))) return false;
    }
    for (int j = 0; j < n_; ++j) {
      if (fin(model_.colLower[j]) && x[j] < model_.colLower[j] - tol * (1 + std::fabs(model_.colLower[j]))) return false;
      if (fin(model_.colUpper[j]) && x[j] > model_.colUpper[j] + tol * (1 + std::fabs(model_.colUpper[j]))) return false;
    }
    return true;
  }
  bool offerIncumbent(std::vector<double> x, const char* source) {
    for (int j = 0; j < n_; ++j)
      if (isInt_[j]) {
        if (std::fabs(x[j] - std::round(x[j])) > kIntTol) return false;
        x[j] = std::round(x[j]);
      }
    if (!rowsFeasible(x, 5e-7)) return false;
    double obj = model_.objective(x);
    if (obj >= incumbentObj_ - 1e-12 * (1 + std::fabs(obj))) return false;
    incumbentObj_ = obj;
    incumbent_ = std::move(x);
    incumbentSource_ = source;
    stats_.incumbentUpdates++;
    if (stats_.firstIncumbentSeconds < 0) stats_.firstIncumbentSeconds = timer_.seconds();
    if (source != std::string("lp")) stats_.heuristicSolutions++;
    recordProgress();
    if (opt_.logLevel >= 2)
      PLOG_INFO("  mip: new incumbent %.10g (%s) at %.2fs", obj, source, timer_.seconds());
    return true;
  }
  void recordProgress() {
    Json p = Json::array();
    p.push(timer_.seconds());
    p.push(fin(incumbentObj_) ? incumbentObj_ : 1e308);
    p.push(fin(globalBound_) ? globalBound_ : -1e308);
    progress_.push(p);
  }

  // ---- heuristics ----
  void roundingHeuristic(const std::vector<double>& x) {
    std::vector<double> r = x;
    for (int j = 0; j < n_; ++j)
      if (isInt_[j]) r[j] = std::min(std::max(std::round(x[j]), curLo_[j]), curUp_[j]);
    if (rowsFeasible(r, 5e-7)) offerIncumbent(r, "rounding");
  }
  // Integers of an (almost) integral LP point are fixed at their rounded values
  // and the continuous part is re-optimized, so round-off in the LP solution can
  // never cost us a feasible incumbent. LP state is restored afterwards.
  bool fixAndResolve(const std::vector<double>& x, const std::vector<double>& lo, const std::vector<double>& up,
                     const char* source) {
    auto basis = saveBasis();
    std::vector<double> l2 = lo, u2 = up;
    for (int j = 0; j < n_; ++j)
      if (isInt_[j]) {
        double v = std::min(std::max(std::round(x[j]), lo[j]), up[j]);
        l2[j] = u2[j] = v;
      }
    applyBounds(l2, u2);
    Status st = solveLp(20000);
    bool ok = st == Status::Optimal && offerIncumbent(lp_->colValues(), source);
    restoreLp(lo, up, *basis);
    return ok;
  }
  // Integer columns with any deviation from integrality (tight tolerance).
  std::vector<int> deviations(const std::vector<double>& x) const {
    std::vector<int> f;
    for (int j = 0; j < n_; ++j)
      if (isInt_[j] && std::fabs(x[j] - std::round(x[j])) > 1e-9) f.push_back(j);
    return f;
  }
  std::vector<int> fractionals(const std::vector<double>& x) const {
    std::vector<int> f;
    for (int j = 0; j < n_; ++j)
      if (isInt_[j] && std::fabs(x[j] - std::round(x[j])) > kIntTol) f.push_back(j);
    return f;
  }
  // Restore LP to (lo, up, basis) and re-solve to the node optimum.
  void restoreLp(const std::vector<double>& lo, const std::vector<double>& up, const BasisVec& basis) {
    applyBounds(lo, up);
    loadBasis(basis);
    solveLp();
  }
  void diving(const std::vector<double>& lo0, const std::vector<double>& up0, const char* name) {
    auto basis = saveBasis();
    std::vector<double> lo = lo0, up = up0;
    for (int depth = 0; depth < 60; ++depth) {
      if (dl_ && dl_->expired()) break;
      std::vector<double> x = lp_->colValues();
      std::vector<int> fr = fractionals(x);
      if (fr.empty()) {
        offerIncumbent(x, name);
        break;
      }
      // Fractional diving: fix the least fractional variable to its nearest integer.
      int best = fr[0];
      double bestF = 1;
      for (int j : fr) {
        double f = std::fabs(x[j] - std::round(x[j]));
        if (f < bestF) {
          bestF = f;
          best = j;
        }
      }
      double v = std::min(std::max(std::round(x[best]), lo[best]), up[best]);
      lo[best] = up[best] = v;
      applyBounds(lo, up);
      Status st = solveLp(2000);
      if (st != Status::Optimal) break;
      if (lp_->objective() >= cutoff()) break;
    }
    restoreLp(lo0, up0, *basis);
  }
  void feasibilityPump(const std::vector<double>& lo, const std::vector<double>& up) {
    // Binary part only; general integers are rounded.
    std::vector<int> bins;
    for (int j = 0; j < n_; ++j)
      if (isInt_[j] && lo[j] == 0 && up[j] == 1) bins.push_back(j);
    if (bins.empty()) return;
    auto basis = saveBasis();
    std::vector<double> x = lp_->colValues();
    std::vector<double> xt(n_, 0.0), prev;
    for (int it = 0; it < 30; ++it) {
      if (dl_ && dl_->expired()) break;
      for (int j = 0; j < n_; ++j) xt[j] = isInt_[j] ? std::round(x[j]) : x[j];
      if (!prev.empty() && prev == xt) {  // cycle: flip the most uncertain binaries
        std::vector<std::pair<double, int>> sc;
        for (int j : bins) sc.emplace_back(-std::fabs(x[j] - xt[j]), j);
        std::sort(sc.begin(), sc.end());
        int T = 1 + rng_.below(std::max(1, static_cast<int>(bins.size()) / 10 + 1));
        for (int k = 0; k < T && k < static_cast<int>(sc.size()); ++k) xt[sc[k].second] = 1 - xt[sc[k].second];
      }
      prev = xt;
      if (rowsFeasible(xt, 5e-7) && offerIncumbent(xt, "feaspump")) break;
      for (int j = 0; j < n_; ++j) lp_->setColCost(j, 0.0);
      for (int j : bins) lp_->setColCost(j, xt[j] < 0.5 ? 1.0 : -1.0);
      Status st = solveLp(3000);
      if (st != Status::Optimal) break;
      x = lp_->colValues();
      if (fractionals(x).empty() && offerIncumbent(x, "feaspump")) break;
    }
    for (int j = 0; j < n_; ++j) lp_->setColCost(j, lpModel_.colCost[j]);
    restoreLp(lo, up, *basis);
  }

  // ---- branching ----
  double pcScore(int j, double f) const {
    double avgD = pcAvg(pcDown_, pcDownN_), avgU = pcAvg(pcUp_, pcUpN_);
    double d = (pcDownN_[j] ? pcDown_[j] / pcDownN_[j] : avgD) * f;
    double u = (pcUpN_[j] ? pcUp_[j] / pcUpN_[j] : avgU) * (1 - f);
    return std::max(d, 1e-6) * std::max(u, 1e-6);
  }
  static double pcAvg(const std::vector<double>& s, const std::vector<int>& c) {
    double sum = 0;
    long long cnt = 0;
    for (size_t j = 0; j < s.size(); ++j)
      if (c[j]) {
        sum += s[j] / c[j];
        ++cnt;
      }
    return cnt ? sum / cnt : 1.0;
  }
  void updatePc(int j, int dir, double gain, double dist) {
    if (dist <= 1e-9 || !std::isfinite(gain)) return;
    double g = std::max(gain, 0.0) / dist;
    if (dir < 0) {
      pcDown_[j] += g;
      pcDownN_[j]++;
    } else {
      pcUp_[j] += g;
      pcUpN_[j]++;
    }
  }
  int selectBranch(const std::vector<double>& x, const std::vector<int>& fr, const std::vector<double>& lo,
                   const std::vector<double>& up, double parentObj) {
    const std::string& rule = opt_.branching;
    if (rule == "mostfrac") {
      int best = fr[0];
      double bf = -1;
      for (int j : fr) {
        double f = std::fabs(x[j] - std::floor(x[j]) - 0.5);
        if (-f > bf) {
          bf = -f;
          best = j;
        }
      }
      return best;
    }
    // Rank by pseudocost score.
    std::vector<std::pair<double, int>> cand;
    for (int j : fr) cand.emplace_back(-pcScore(j, x[j] - std::floor(x[j])), j);
    std::sort(cand.begin(), cand.end());
    if (rule == "pseudocost") return cand[0].second;
    // Strong branching on unreliable (or all, for "strong") top candidates.
    const bool full = rule == "strong";
    auto basis = saveBasis();
    int best = cand[0].second;
    double bestScore = -1;
    int evaluated = 0;
    for (auto& c : cand) {
      int j = c.second;
      bool reliable = std::min(pcDownN_[j], pcUpN_[j]) >= opt_.reliability;
      double f = x[j] - std::floor(x[j]);
      double score;
      if ((full || !reliable) && evaluated < opt_.strongCandidates) {
        ++evaluated;
        double gains[2];
        for (int d = 0; d < 2; ++d) {
          std::vector<double> l2 = lo, u2 = up;
          if (d == 0) u2[j] = std::floor(x[j]);
          else l2[j] = std::ceil(x[j]);
          applyBounds(l2, u2);
          Status st = solveLp(opt_.strongIterations);
          stats_.strongBranchLps++;
          if (st == Status::Infeasible) gains[d] = 1e20;
          else gains[d] = std::max(0.0, lp_->objective() - parentObj);
          if (st == Status::Optimal || st == Status::IterationLimit)
            updatePc(j, d == 0 ? -1 : 1, gains[d], d == 0 ? f : 1 - f);
          applyBounds(lo, up);
          loadBasis(*basis);
        }
        score = std::max(gains[0], 1e-6) * std::max(gains[1], 1e-6);
      } else {
        score = -c.first;
      }
      if (score > bestScore) {
        bestScore = score;
        best = j;
      }
    }
    if (evaluated > 0) solveLp();  // back to the node optimum
    return best;
  }

  // ---- root cuts ----
  void rootCuts(std::vector<double>& x, double& obj) {
    for (int round = 0; round < opt_.cutRounds; ++round) {
      if (dl_ && dl_->expired()) break;
      if (fractionals(x).empty()) break;
      CutGenerator gen(lpModel_);
      std::vector<Cut> all = gen.gomory(*lp_, x, rootLo_, rootUp_, opt_.maxCutsPerRound);
      std::vector<Cut> mir = gen.mir(x, rootLo_, rootUp_, opt_.maxCutsPerRound);
      std::vector<Cut> cov = gen.covers(x, rootLo_, rootUp_, opt_.maxCutsPerRound);
      all.insert(all.end(), mir.begin(), mir.end());
      all.insert(all.end(), cov.begin(), cov.end());
      stats_.cutsRejected += gen.stats.rejectedNumerics;
      std::vector<Cut> sel = CutGenerator::select(all, opt_.maxCutsPerRound, 1e-4, 0.99);
      if (sel.empty()) break;
      std::vector<double> lower, upper, vals;
      std::vector<int> starts{0}, cols;
      for (auto& c : sel) {
        lower.push_back(c.rhs);
        upper.push_back(kInf);
        cols.insert(cols.end(), c.idx.begin(), c.idx.end());
        vals.insert(vals.end(), c.val.begin(), c.val.end());
        starts.push_back(static_cast<int>(cols.size()));
        if (c.type == "gmi") stats_.cutsGomory++;
        else if (c.type == "mir") stats_.cutsMir++;
        else stats_.cutsCover++;
      }
      lp_->addRows(lower, upper, starts, cols, vals);
      lpModel_.addRows(lower, upper, starts, cols, vals);
      Status st = solveLp();
      stats_.cutRoundsDone++;
      if (st != Status::Optimal) {
        PLOG_DETAIL("mip: LP not optimal after cut round (%s); stopping cuts", statusName(st));
        break;
      }
      double newObj = lp_->objective();
      x = lp_->colValues();
      double improve = newObj - obj;
      obj = newObj;
      if (opt_.logLevel >= 2)
        PLOG_INFO("  mip: cut round %d: +%zu cuts, bound %.10g", round + 1, sel.size(), obj);
      if (improve < 1e-4 * (1 + std::fabs(obj)) && round >= 1) break;
    }
  }

  const Model& model_;
  Model lpModel_;  // model + cuts
  MipOptions opt_;
  const Deadline* dl_;
  Rng rng_;
  int n_ = 0;
  std::vector<char> isInt_;
  std::vector<double> rootLo_, rootUp_, curLo_, curUp_;
  std::unique_ptr<Simplex> lp_;
  Propagator prop_;
  double incumbentObj_ = kInf;
  std::vector<double> incumbent_;
  std::string incumbentSource_;
  std::vector<double> pcDown_, pcUp_;
  std::vector<int> pcDownN_, pcUpN_;
  bool objIntegral_ = false;
  double globalBound_ = -kInf;
  bool boundSafe_ = true;
  MipStats stats_;
  Timer timer_;
  Json progress_ = Json::array();
};

MipResult BranchAndBound::run() {
  MipResult res;
  timer_.reset();
  lp_ = std::make_unique<Simplex>(lpModel_, opt_.lp);
  curLo_ = rootLo_;
  curUp_ = rootUp_;
  prop_.setup(model_);
  {
    PropagationResult pr = prop_.propagate(rootLo_, rootUp_, 10);
    stats_.propagationTightenings += pr.tightened;
    if (!pr.infeasible) applyBounds(rootLo_, rootUp_);
  }

  // ---- Root LP ----
  Status st = solveLp();
  if (st == Status::Infeasible) {
    res.status = Status::Infeasible;
    res.stats = stats_;
    return res;
  }
  if (st == Status::Unbounded) {
    res.status = Status::InfeasibleOrUnbounded;
    res.stats = stats_;
    return res;
  }
  if (st != Status::Optimal) {
    res.status = st == Status::TimeLimit ? Status::TimeLimit : Status::NumericalFailure;
    res.stats = stats_;
    return res;
  }
  std::vector<double> x = lp_->colValues();
  double obj = lp_->objective();
  stats_.rootLpBound = obj;
  if (opt_.logLevel >= 1)
    PLOG_INFO("mip: root LP %.10g (%lld it, %.2fs)", obj, stats_.lpIterations, timer_.seconds());
  if (opt_.heuristics) roundingHeuristic(x);
  if (opt_.cuts) rootCuts(x, obj);
  stats_.rootBoundAfterCuts = obj;
  if (fractionals(x).empty() && !offerIncumbent(x, "lp")) fixAndResolve(x, rootLo_, rootUp_, "lp-fixed");
  if (opt_.heuristics && !fractionals(x).empty()) {
    diving(rootLo_, rootUp_, "diving");
    if (!fin(incumbentObj_)) feasibilityPump(rootLo_, rootUp_);
    x = lp_->colValues();
    obj = lp_->objective();
  }
  stats_.rootSeconds = timer_.seconds();
  bool safe = true;
  double rootBound = roundBound(safeBound(rootLo_, rootUp_, obj, &safe));
  boundSafe_ = safe;
  globalBound_ = rootBound;
  recordProgress();
  if (opt_.logLevel >= 1)
    PLOG_INFO("mip: root bound %.10g after %d cut rounds (gmi %d, mir %d, cover %d), incumbent %s", rootBound,
              stats_.cutRoundsDone, stats_.cutsGomory, stats_.cutsMir, stats_.cutsCover,
              fin(incumbentObj_) ? formatString("%.10g", incumbentObj_).c_str() : "none");

  // ---- Tree search ----
  std::priority_queue<std::shared_ptr<Node>, std::vector<std::shared_ptr<Node>>, HeapCmp> heap;
  long long nextId = 1;
  auto root = std::make_shared<Node>();
  root->id = 0;
  root->bound = rootBound;
  root->basis = saveBasis();
  std::shared_ptr<Node> plunge = root;
  std::vector<double> unsolvedBounds;  // nodes we could not solve: their bounds stay in the global bound
  bool limitHit = false;
  double lastLog = 0;
  auto openBound = [&]() {
    double b = kInf;
    if (plunge) b = std::min(b, plunge->bound);
    if (!heap.empty()) b = std::min(b, heap.top()->bound);
    for (double u : unsolvedBounds) b = std::min(b, u);
    return b;
  };

  while (plunge || !heap.empty()) {
    if (dl_ && dl_->expired()) {
      limitHit = true;
      res.status = Status::TimeLimit;
      break;
    }
    if (opt_.nodeLimit >= 0 && stats_.nodes >= opt_.nodeLimit) {
      limitHit = true;
      res.status = Status::NodeLimit;
      break;
    }
    // Global bound & termination by gap.
    double gb = std::min(openBound(), incumbentObj_);
    if (gb > globalBound_) {
      globalBound_ = gb;
      recordProgress();
    }
    if (fin(incumbentObj_)) {
      double gap = (incumbentObj_ - globalBound_) / std::max(1.0, std::fabs(incumbentObj_));
      if (gap <= opt_.relGap) break;
    }
    std::shared_ptr<Node> node;
    bool warm = false;
    if (plunge) {
      node = plunge;
      plunge.reset();
      warm = true;  // LP state is the parent's (or root's) optimum
    } else {
      node = heap.top();
      heap.pop();
    }
    if (node->bound >= cutoff()) {
      stats_.boundPrunes++;
      continue;
    }
    // Node bounds.
    std::vector<double> lo = rootLo_, up = rootUp_;
    for (auto& c : node->changes) {
      lo[c.var] = std::max(lo[c.var], c.lo);
      up[c.var] = std::min(up[c.var], c.up);
    }
    PropagationResult pr = prop_.propagate(lo, up, 3);
    stats_.propagationTightenings += pr.tightened;
    if (pr.infeasible) {
      stats_.propagationPrunes++;
      continue;
    }
    applyBounds(lo, up);
    if (!warm && node->basis) loadBasis(*node->basis);
    Status s = solveLp();
    stats_.nodes++;
    stats_.maxDepth = std::max(stats_.maxDepth, node->depth);
    if (s == Status::Infeasible) {
      stats_.farkasPrunes++;
      if (node->branchVar >= 0) updatePc(node->branchVar, node->branchDir, 1e3 * (1 + std::fabs(node->parentObj)), node->branchFrac);
      continue;
    }
    if (s != Status::Optimal) {
      if (s == Status::TimeLimit) {
        unsolvedBounds.push_back(node->bound);
        limitHit = true;
        res.status = Status::TimeLimit;
        break;
      }
      // Could not solve this node reliably: keep its (safe) parent bound in the global bound.
      unsolvedBounds.push_back(node->bound);
      PLOG_DETAIL("mip: node %lld LP %s; keeping parent bound", node->id, statusName(s));
      continue;
    }
    x = lp_->colValues();
    obj = lp_->objective();
    if (node->branchVar >= 0) updatePc(node->branchVar, node->branchDir, obj - node->parentObj, node->branchFrac);
    bool nodeSafe = true;
    double nb = roundBound(std::max(node->bound, safeBound(lo, up, obj, &nodeSafe)));
    if (!nodeSafe) boundSafe_ = false;
    if (nb >= cutoff()) {
      stats_.boundPrunes++;
      continue;
    }
    std::vector<int> fr = fractionals(x);
    if (fr.empty()) {
      stats_.integralLeaves++;
      if (offerIncumbent(x, "lp") || fixAndResolve(x, lo, up, "lp-fixed") || nb >= cutoff()) continue;
      // Integral within tolerance but not acceptable even after re-optimizing the
      // continuous part: never discard a node without proof - keep branching.
      fr = deviations(x);
      if (fr.empty()) {
        unsolvedBounds.push_back(nb);  // cannot branch further: its bound stays in the global bound
        continue;
      }
    }
    if (opt_.heuristics) {
      roundingHeuristic(x);
      if (stats_.nodes % 200 == 0) diving(lo, up, "diving");
      if (nb >= cutoff()) {
        stats_.boundPrunes++;
        continue;
      }
    }
    int j = selectBranch(x, fr, lo, up, obj);
    x = lp_->colValues();  // strong branching re-solved the node
    double xj = x[j];
    if (std::fabs(xj - std::round(xj)) <= kIntTol) {  // strong branching changed the solution
      fr = fractionals(x);
      if (fr.empty()) fr = deviations(x);
      if (fr.empty()) {
        if (!offerIncumbent(x, "lp") && !fixAndResolve(x, lo, up, "lp-fixed") && nb < cutoff())
          unsolvedBounds.push_back(nb);
        continue;
      }
      j = fr[0];
      xj = x[j];
    }
    auto basis = saveBasis();
    auto down = std::make_shared<Node>();
    auto upn = std::make_shared<Node>();
    for (auto* c : {down.get(), upn.get()}) {
      c->depth = node->depth + 1;
      c->bound = nb;
      c->changes = node->changes;
      c->basis = basis;
      c->branchVar = j;
      c->parentObj = obj;
    }
    down->id = nextId++;
    upn->id = nextId++;
    down->changes.push_back({j, -kInf, std::floor(xj)});
    upn->changes.push_back({j, std::ceil(xj), kInf});
    down->branchDir = -1;
    upn->branchDir = 1;
    double f = xj - std::floor(xj);
    down->branchFrac = f;
    upn->branchFrac = 1 - f;
    down->estimate = obj + pcScore(j, f);
    upn->estimate = down->estimate;
    // Plunge into the child in the rounding direction, queue the other.
    if (f >= 0.5) {
      plunge = upn;
      heap.push(down);
    } else {
      plunge = down;
      heap.push(upn);
    }
    if (node->depth > 200) {  // stop plunging very deep
      heap.push(plunge);
      plunge.reset();
    }
    if (opt_.logLevel >= 1 && timer_.seconds() - lastLog > 2.0) {
      lastLog = timer_.seconds();
      double ob = std::min(openBound(), incumbentObj_);
      PLOG_INFO("  mip: %8lld nodes  %7zu open  inc %-16s bound %.10g  gap %s  %.1fs", stats_.nodes, heap.size(),
                fin(incumbentObj_) ? formatString("%.10g", incumbentObj_).c_str() : "-", ob,
                fin(incumbentObj_) ? formatString("%.3g%%", 100 * (incumbentObj_ - ob) / std::max(1.0, std::fabs(incumbentObj_))).c_str() : "-",
                timer_.seconds());
    }
  }

  // ---- Final status ----
  double finalBound = std::min(openBound(), incumbentObj_);
  if (!limitHit && heap.empty() && !plunge && unsolvedBounds.empty()) finalBound = incumbentObj_;
  if (!limitHit && fin(incumbentObj_)) {
    double gap = (incumbentObj_ - std::min(finalBound, incumbentObj_)) / std::max(1.0, std::fabs(incumbentObj_));
    if (gap <= opt_.relGap || (heap.empty() && !plunge)) {
      res.status = unsolvedBounds.empty() || gap <= opt_.relGap ? Status::Optimal : Status::NumericalFailure;
      if (heap.empty() && !plunge && unsolvedBounds.empty()) finalBound = incumbentObj_;
    }
  } else if (!limitHit && !fin(incumbentObj_)) {
    res.status = unsolvedBounds.empty() ? Status::Infeasible : Status::NumericalFailure;
  }
  if (res.status == Status::NotSolved) res.status = Status::TimeLimit;
  globalBound_ = std::max(globalBound_, std::min(finalBound, incumbentObj_));
  if (res.status == Status::Optimal && globalBound_ > incumbentObj_) globalBound_ = incumbentObj_;
  recordProgress();
  res.objective = incumbentObj_;
  res.bestBound = globalBound_;
  res.boundSafe = boundSafe_;
  res.x = incumbent_;
  res.incumbentSource = incumbentSource_;
  res.stats = stats_;
  res.progress = progress_;
  if (opt_.logLevel >= 1)
    PLOG_INFO("mip: %s  incumbent %s  bound %.10g  nodes %lld  lp it %lld  %.2fs", statusName(res.status),
              fin(incumbentObj_) ? formatString("%.10g", incumbentObj_).c_str() : "none", globalBound_, stats_.nodes,
              stats_.lpIterations, timer_.seconds());
  return res;
}

}  // namespace

MipResult solveMip(const Model& model, const MipOptions& opts, const Deadline* deadline) {
  BranchAndBound bnb(model, opts, deadline);
  return bnb.run();
}

}  // namespace pramana
