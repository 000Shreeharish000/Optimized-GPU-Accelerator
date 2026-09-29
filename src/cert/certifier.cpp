#include "cert/certifier.h"

#include <algorithm>
#include <cmath>

#include "cert/interval.h"
#include "mip/propagate.h"
#include "util/log.h"

namespace pramana {

const char* verdictName(Verdict v) {
  switch (v) {
    case Verdict::Pass: return "PASS";
    case Verdict::Warn: return "WARN";
    case Verdict::Fail: return "FAIL";
    case Verdict::NotApplicable: return "N/A";
  }
  return "?";
}

namespace {

inline bool fin(double v) { return std::fabs(v) < kInfBoundThreshold; }
inline double clean(double v) { return fin(v) ? v : (v > 0 ? kInf : -kInf); }

// z_j = cost_j - (A^T y)_j as rigorous intervals, with cost given as an interval.
std::vector<Interval> reducedCostIntervals(const SparseMatrix& A, const std::vector<double>& costLo,
                                           const std::vector<double>& costHi, const std::vector<double>& y) {
  std::vector<Interval> z(A.numCols);
  for (int j = 0; j < A.numCols; ++j) {
    Interval acc{costLo[j], costHi[j]};
    for (int k = A.start[j]; k < A.start[j + 1]; ++k) acc = sub(acc, mulExact(A.value[k], y[A.index[k]]));
    z[j] = acc;
  }
  return z;
}

// Row multipliers can be projected onto their sign cone for free (the bound is
// valid for ANY y); this removes round-off residue on one-sided rows exactly.
std::vector<double> projectRowDuals(const std::vector<double>& y, const std::vector<double>& rowLower,
                                    const std::vector<double>& rowUpper, double* moved) {
  std::vector<double> p = y;
  double mv = 0;
  for (size_t i = 0; i < y.size(); ++i) {
    bool lo = fin(rowLower[i]), up = fin(rowUpper[i]);
    double v = y[i];
    if (!lo && !up) v = 0;
    else if (!lo) v = std::min(v, 0.0);   // only upper side: y <= 0
    else if (!up) v = std::max(v, 0.0);   // only lower side: y >= 0
    mv = std::max(mv, std::fabs(v - y[i]));
    p[i] = v;
  }
  if (moved) *moved = mv;
  return p;
}

SafeBoundResult safeBoundInterval(const SparseMatrix& A, const std::vector<double>& costLo,
                                  const std::vector<double>& costHi, const std::vector<double>& colLower,
                                  const std::vector<double>& colUpper, const std::vector<double>& rowLower,
                                  const std::vector<double>& rowUpper, const std::vector<double>& yIn) {
  SafeBoundResult res;
  std::vector<double> y = projectRowDuals(yIn, rowLower, rowUpper, &res.rowProjection);
  std::vector<Interval> z = reducedCostIntervals(A, costLo, costHi, y);
  double lb = 0;
  bool minusInf = false;
  for (int i = 0; i < A.numRows; ++i)
    lb = addDown(lb, minProduct(Interval::point(y[i]), clean(rowLower[i]), clean(rowUpper[i])));
  for (int j = 0; j < A.numCols; ++j) {
    Interval zj = z[j];
    // Wrong-sign reduced cost on an infinite bound: clip it and record the cost
    // perturbation needed (the bound is then rigorous for costs c + delta).
    double shift = 0;
    if (!fin(colLower[j]) && zj.hi > 0) {
      shift = std::max(shift, zj.hi);
      zj.hi = 0;
      zj.lo = std::min(zj.lo, 0.0);
    }
    if (!fin(colUpper[j]) && zj.lo < 0) {
      shift = std::max(shift, -zj.lo);
      zj.lo = 0;
      zj.hi = std::max(zj.hi, 0.0);
    }
    if (shift > 0) {
      double rel = shift / (1.0 + std::fabs(costHi[j]));
      res.costPerturbation = std::max(res.costPerturbation, rel);
      res.maxDualInfeasibility = std::max(res.maxDualInfeasibility, shift);
      if (rel > res.maxAllowedPerturbation) {
        res.infiniteTerms++;
        minusInf = true;
        continue;
      }
    }
    double t = minProduct(zj, clean(colLower[j]), clean(colUpper[j]));
    if (std::isinf(t)) {
      res.infiniteTerms++;
      minusInf = true;
      continue;
    }
    lb = addDown(lb, t);
  }
  res.bound = minusInf ? -kInf : lb;
  return res;
}

void addCheck(Certificate& c, const std::string& name, Verdict v, double value, double tol,
              const std::string& detail = "") {
  c.checks.push_back({name, v, value, tol, detail});
}

}  // namespace

SafeBoundResult safeDualBound(const SparseMatrix& A, const std::vector<double>& cost,
                              const std::vector<double>& colLower, const std::vector<double>& colUpper,
                              const std::vector<double>& rowLower, const std::vector<double>& rowUpper,
                              const std::vector<double>& y) {
  return safeBoundInterval(A, cost, cost, colLower, colUpper, rowLower, rowUpper, y);
}

bool verifyFarkas(const SparseMatrix& A, const std::vector<double>& colLower, const std::vector<double>& colUpper,
                  const std::vector<double>& rowLower, const std::vector<double>& rowUpper,
                  const std::vector<double>& yIn, double* margin, double* perturbation) {
  // For every feasible x: (A^T y)^T x = y^T r with r = Ax in [L,U]. Infeasible if
  //   sup_x (A^T y)^T x < inf_r y^T r   (tested for y and -y).
  // Each sign first projects y onto the cone making inf_r y^T r finite (any y is allowed).
  double best = -kInf, bestPert = 0;
  for (int sgn = -1; sgn <= 1; sgn += 2) {
    std::vector<double> y(yIn.size());
    for (size_t i = 0; i < y.size(); ++i) y[i] = sgn * yIn[i];
    y = projectRowDuals(y, rowLower, rowUpper, nullptr);
    std::vector<double> zero(A.numCols, 0.0);
    std::vector<Interval> z = reducedCostIntervals(A, zero, zero, y);  // z = -A^T y
    double supX = 0, infR = 0, relPert = 0;
    bool bad = false;
    for (int j = 0; j < A.numCols && !bad; ++j) {
      Interval w{-z[j].hi, -z[j].lo};  // A^T y
      // sup_x w x is finite only if w <= 0 when u = inf and w >= 0 when l = -inf. Round-off
      // residue is clipped and reported as backward error: relative to the magnitude of the
      // contributions sum_i |a_ij y_i| (i.e. a relative perturbation of column j of A).
      double mag = 0;
      for (int k = A.start[j]; k < A.start[j + 1]; ++k) mag += std::fabs(A.value[k] * y[A.index[k]]);
      double shift = 0;
      if (!fin(colUpper[j]) && w.hi > 0) {
        shift = std::max(shift, w.hi);
        w.hi = 0;
        w.lo = std::min(w.lo, 0.0);
      }
      if (!fin(colLower[j]) && w.lo < 0) {
        shift = std::max(shift, -w.lo);
        w.lo = 0;
        w.hi = std::max(w.hi, 0.0);
      }
      if (shift > 0) relPert = std::max(relPert, mag > 0 ? shift / mag : kInf);
      double mx = maxProduct(w, clean(colLower[j]), clean(colUpper[j]));
      if (std::isinf(mx)) bad = true;
      else supX = addUp(supX, mx);
    }
    if (bad) continue;
    for (size_t i = 0; i < y.size(); ++i)
      infR = addDown(infR, minProduct(Interval::point(y[i]), clean(rowLower[i]), clean(rowUpper[i])));
    if (relPert > 1e-9) continue;
    double m = infR - supX;
    if (m > best) {
      best = m;
      bestPert = relPert;
    }
  }
  if (margin) *margin = best;
  if (perturbation) *perturbation = bestPert;
  return best > 0;
}

std::string Certificate::summary() const {
  std::string s = formatString("certificate: claimed %s -> %s (%s)", statusName(claimed), statusName(certified),
                               accepted ? "ACCEPTED" : "REJECTED");
  for (auto& c : checks)
    s += formatString("\n  %-28s %-4s value %.3e  tol %.1e  %s", c.name.c_str(), verdictName(c.verdict), c.value,
                      c.tolerance, c.detail.c_str());
  return s;
}

Json Certificate::toJson() const {
  Json j = Json::object();
  j["claimed_status"] = statusName(claimed);
  j["certified_status"] = statusName(certified);
  j["accepted"] = accepted;
  j["primal_objective"] = primalObjective;
  j["safe_dual_bound"] = safeDualBound;
  j["certified_relative_gap"] = certifiedGap;
  j["max_primal_violation"] = maxPrimalViolation;
  j["max_dual_violation"] = maxDualViolation;
  j["max_integrality_violation"] = maxIntegralityViolation;
  j["cost_perturbation_backward_error"] = costPerturbation;
  Json arr = Json::array();
  for (auto& c : checks) {
    Json o = Json::object();
    o["check"] = c.name;
    o["verdict"] = verdictName(c.verdict);
    o["value"] = c.value;
    o["tolerance"] = c.tolerance;
    if (!c.detail.empty()) o["detail"] = c.detail;
    arr.push(o);
  }
  j["checks"] = arr;
  return j;
}

Certificate certify(const Model& model, const SolutionClaim& claim, const CertTolerances& tol) {
  Certificate cert;
  cert.claimed = claim.status;
  cert.certified = claim.status;
  const int n = model.numCols(), m = model.numRows();
  const double sense = model.sense == ObjSense::Maximize ? -1.0 : 1.0;
  // Minimization-form cost.
  std::vector<double> cmin(n);
  for (int j = 0; j < n; ++j) cmin[j] = sense * model.colCost[j];

  auto checkPrimal = [&](const std::vector<double>& x, bool integrality) -> bool {
    if (static_cast<int>(x.size()) != n) {
      addCheck(cert, "primal_dimension", Verdict::Fail, static_cast<double>(x.size()), n);
      return false;
    }
    double colViol = 0, rowViol = 0, intViol = 0;
    int worstRow = -1;
    for (int j = 0; j < n; ++j) {
      double v = 0;
      if (!std::isfinite(x[j])) {  // NaN / inf can never be a feasible point
        colViol = kInf;
        continue;
      }
      if (fin(model.colLower[j]) && x[j] < model.colLower[j])
        v = (model.colLower[j] - x[j]) / (1 + std::fabs(model.colLower[j]));
      if (fin(model.colUpper[j]) && x[j] > model.colUpper[j])
        v = std::max(v, (x[j] - model.colUpper[j]) / (1 + std::fabs(model.colUpper[j])));
      colViol = std::max(colViol, v);
      if (integrality && model.colType[j] == VarType::Integer)
        intViol = std::max(intViol, std::fabs(x[j] - std::round(x[j])));
    }
    std::vector<CompensatedSum> act(m);
    for (int j = 0; j < n; ++j)
      for (int k = model.A.start[j]; k < model.A.start[j + 1]; ++k) act[model.A.index[k]].add(model.A.value[k] * x[j]);
    for (int i = 0; i < m; ++i) {
      double r = act[i].value();
      double v = std::isfinite(r) ? 0.0 : kInf;
      if (fin(model.rowLower[i]) && r < model.rowLower[i])
        v = (model.rowLower[i] - r) / (1 + std::fabs(model.rowLower[i]));
      if (fin(model.rowUpper[i]) && r > model.rowUpper[i])
        v = std::max(v, (r - model.rowUpper[i]) / (1 + std::fabs(model.rowUpper[i])));
      if (v > rowViol) {
        rowViol = v;
        worstRow = i;
      }
    }
    cert.maxPrimalViolation = std::max(colViol, rowViol);
    addCheck(cert, "primal_bound_violation", colViol <= tol.primalFeas ? Verdict::Pass : Verdict::Fail, colViol,
             tol.primalFeas, "max relative column-bound violation (unscaled)");
    addCheck(cert, "primal_row_violation", rowViol <= tol.primalFeas ? Verdict::Pass : Verdict::Fail, rowViol,
             tol.primalFeas,
             worstRow >= 0 ? "worst row " + (worstRow < static_cast<int>(model.rowNames.size())
                                                 ? model.rowNames[worstRow] : std::to_string(worstRow))
                           : "compensated row activities");
    bool ok = colViol <= tol.primalFeas && rowViol <= tol.primalFeas;
    if (integrality) {
      cert.maxIntegralityViolation = intViol;
      addCheck(cert, "integrality_violation", intViol <= tol.integrality ? Verdict::Pass : Verdict::Fail, intViol,
               tol.integrality);
      ok = ok && intViol <= tol.integrality;
    }
    cert.primalObjective = model.objective(x);
    return ok;
  };

  const bool isMip = model.isMip();
  switch (claim.status) {
    case Status::Optimal:
    case Status::TimeLimit:
    case Status::NodeLimit:
    case Status::IterationLimit: {
      if (claim.x.empty()) {
        addCheck(cert, "solution_present", claim.status == Status::Optimal ? Verdict::Fail : Verdict::NotApplicable, 0, 0,
                 "no primal solution");
        cert.accepted = claim.status != Status::Optimal;
        if (claim.status == Status::Optimal) cert.certified = Status::NumericalFailure;
        break;
      }
      bool primalOk = checkPrimal(claim.x, isMip);
      if (isMip) {
        double bound = claim.bestBound;
        double gap = std::fabs(cert.primalObjective - bound) / std::max(1.0, std::fabs(cert.primalObjective));
        if (!std::isfinite(bound)) gap = kInf;
        cert.safeDualBound = bound;
        cert.certifiedGap = gap;
        addCheck(cert, "mip_bound_is_safe", claim.boundIsSafe ? Verdict::Pass : Verdict::Warn,
                 claim.boundIsSafe ? 1 : 0, 1,
                 claim.boundIsSafe ? "global bound = min over open nodes of Neumaier-Shcherbina safe LP bounds"
                                   : "bound from floating-point LP values");
        Verdict gv = gap <= tol.mipRelGap ? Verdict::Pass : (claim.status == Status::Optimal ? Verdict::Fail : Verdict::Warn);
        addCheck(cert, "mip_relative_gap", gv, gap, tol.mipRelGap, "|incumbent - bound| / max(1,|incumbent|)");
        if (claim.status == Status::Optimal) {
          cert.accepted = primalOk && gap <= tol.mipRelGap;
          if (!cert.accepted) cert.certified = primalOk ? Status::TimeLimit : Status::NumericalFailure;
        } else {
          cert.accepted = primalOk;
          if (!primalOk) cert.certified = Status::NumericalFailure;
        }
        break;
      }
      if (claim.status != Status::Optimal) {
        cert.accepted = true;
        break;
      }
      // LP / QP optimality: rigorous dual bound.
      if (static_cast<int>(claim.rowDual.size()) != m) {
        addCheck(cert, "dual_present", Verdict::Fail, 0, 0, "no dual solution supplied");
        cert.accepted = false;
        cert.certified = Status::NumericalFailure;
        break;
      }
      std::vector<double> clo = cmin, chi = cmin;
      double constLo = 0;  // lower bound on constant part (min form)
      if (model.isQp()) {
        // g = c + Q x^ (interval), constant = -1/2 x^T Q x^ (lower bound).
        Interval quad{0, 0};
        for (int j = 0; j < n; ++j) {
          Interval gj{clo[j], chi[j]};
          for (int k = model.Q.start[j]; k < model.Q.start[j + 1]; ++k) {
            Interval t = mulExact(sense * model.Q.value[k], claim.x[model.Q.index[k]]);
            gj = add(gj, t);
            // x_j * Q_ij * x_i contributes to x^T Q x
            Interval xx = mulExact(t.lo, claim.x[j]);
            Interval xx2 = mulExact(t.hi, claim.x[j]);
            quad = add(quad, Interval{std::min(xx.lo, xx2.lo), std::max(xx.hi, xx2.hi)});
          }
          clo[j] = gj.lo;
          chi[j] = gj.hi;
        }
        constLo = mulDown(-0.5, quad.hi);
      }
      SafeBoundResult sb = safeBoundInterval(model.A, clo, chi, model.colLower, model.colUpper, model.rowLower,
                                             model.rowUpper, claim.rowDual);
      bool usedImplied = false;
      if (std::isinf(sb.bound)) {
        // Inexact duals (e.g. first-order methods) leave wrong-sign reduced costs on
        // unbounded columns. Bound those columns by constraint-implied bounds
        // (activity propagation, relaxed by a margin): valid for every feasible x.
        Model relaxed = model;
        for (auto& t : relaxed.colType) t = VarType::Continuous;
        Propagator prop;
        prop.setup(relaxed);
        std::vector<double> lo = model.colLower, up = model.colUpper;
        prop.propagate(lo, up, 30);
        for (int j = 0; j < n; ++j) {
          if (fin(lo[j]) && (!fin(model.colLower[j]) || lo[j] > model.colLower[j])) lo[j] -= 1e-6 * (1 + std::fabs(lo[j]));
          else lo[j] = model.colLower[j];
          if (fin(up[j]) && (!fin(model.colUpper[j]) || up[j] < model.colUpper[j])) up[j] += 1e-6 * (1 + std::fabs(up[j]));
          else up[j] = model.colUpper[j];
        }
        SafeBoundResult sb2 = safeBoundInterval(model.A, clo, chi, lo, up, model.rowLower, model.rowUpper, claim.rowDual);
        if (!std::isinf(sb2.bound)) {
          sb = sb2;
          usedImplied = true;
        }
      }
      double pobjMin = sense * cert.primalObjective;
      // min-form objective = c_min^T x + 1/2 x^T Q_min x + sense*offset
      double lbMin = std::isinf(sb.bound) ? -kInf : addDown(addDown(sb.bound, constLo), sense * model.objOffset);
      cert.maxDualViolation = sb.maxDualInfeasibility;
      double gap = std::isinf(lbMin) ? kInf : (pobjMin - lbMin) / std::max(1.0, std::fabs(pobjMin));
      cert.safeDualBound = std::isinf(lbMin) ? (sense > 0 ? -kInf : kInf) : sense * lbMin;
      cert.certifiedGap = gap;
      bool rigorous = !std::isinf(lbMin);
      cert.costPerturbation = sb.costPerturbation;
      Verdict sv = rigorous ? (gap <= tol.relGap ? Verdict::Pass : Verdict::Fail) : Verdict::Warn;
      addCheck(cert, model.isQp() ? "safe_dual_bound_qp_linearized" : "safe_dual_bound_gap", sv, gap, tol.relGap,
               rigorous ? formatString("rigorous bound %.12g (outward-rounded; exact for costs perturbed by <= %.1e rel.)%s",
                                       cert.safeDualBound, sb.costPerturbation,
                                       usedImplied ? " using constraint-implied column bounds" : "")
                        : formatString("%d unbounded directions with nonzero reduced cost", sb.infiniteTerms));
      bool dualOk = rigorous ? gap <= tol.relGap : sb.maxDualInfeasibility <= tol.dualFeas;
      if (!rigorous)
        addCheck(cert, "dual_feasibility", sb.maxDualInfeasibility <= tol.dualFeas ? Verdict::Pass : Verdict::Fail,
                 sb.maxDualInfeasibility, tol.dualFeas, "max wrong-sign reduced cost on infinite bounds");
      cert.accepted = primalOk && dualOk;
      if (!cert.accepted) cert.certified = Status::NumericalFailure;
      break;
    }
    case Status::Infeasible: {
      if (static_cast<int>(claim.farkas.size()) != m) {
        addCheck(cert, "farkas_present", Verdict::Fail, 0, 0, "no Farkas certificate");
        cert.accepted = false;
        cert.certified = Status::NumericalFailure;
        break;
      }
      double margin = 0, pert = 0;
      bool ok = verifyFarkas(model.A, model.colLower, model.colUpper, model.rowLower, model.rowUpper, claim.farkas,
                             &margin, &pert);
      addCheck(cert, "farkas_certificate", ok ? Verdict::Pass : Verdict::Fail, margin, 0,
               formatString("sup (A^T y)^T x < inf y^T r over the boxes (outward-rounded); backward error %.1e", pert));
      cert.accepted = ok;
      if (!ok) cert.certified = Status::NumericalFailure;
      break;
    }
    case Status::Unbounded: {
      bool ok = static_cast<int>(claim.primalRay.size()) == n;
      double dnorm = 0, cd = 0, viol = 0;
      if (ok) {
        for (int j = 0; j < n; ++j) dnorm = std::max(dnorm, std::fabs(claim.primalRay[j]));
        if (dnorm == 0) ok = false;
      }
      if (ok) {
        for (int j = 0; j < n; ++j) {
          double d = claim.primalRay[j] / dnorm;
          cd += cmin[j] * d;
          if (fin(model.colLower[j])) viol = std::max(viol, -d);
          if (fin(model.colUpper[j])) viol = std::max(viol, d);
        }
        std::vector<double> ad(m, 0.0), dn(n);
        for (int j = 0; j < n; ++j) dn[j] = claim.primalRay[j] / dnorm;
        model.A.multiply(dn.data(), ad.data());
        for (int i = 0; i < m; ++i) {
          if (fin(model.rowLower[i])) viol = std::max(viol, -ad[i]);
          if (fin(model.rowUpper[i])) viol = std::max(viol, ad[i]);
        }
      }
      if (ok && model.isQp()) {
        // A QP is unbounded along d only if d^T Q d = 0 (convex Q); otherwise the quadratic term dominates.
        std::vector<double> dn(n), qd(n, 0.0);
        for (int j = 0; j < n; ++j) dn[j] = claim.primalRay[j] / dnorm;
        model.Q.multiply(dn.data(), qd.data());
        double dqd = 0;
        for (int j = 0; j < n; ++j) dqd += dn[j] * qd[j];
        addCheck(cert, "primal_ray_zero_curvature", std::fabs(dqd) <= tol.primalFeas ? Verdict::Pass : Verdict::Fail,
                 dqd, tol.primalFeas, "d^T Q d must vanish along an unbounded direction of a convex QP");
        if (std::fabs(dqd) > tol.primalFeas) ok = false;
      }
      bool rayOk = ok && viol <= tol.primalFeas && cd < -tol.dualFeas;
      addCheck(cert, "primal_ray_recession", ok && viol <= tol.primalFeas ? Verdict::Pass : Verdict::Fail, viol,
               tol.primalFeas, "ray stays in the recession cone of all bounds");
      addCheck(cert, "primal_ray_descent", ok && cd < -tol.dualFeas ? Verdict::Pass : Verdict::Fail, cd, -tol.dualFeas,
               "c^T d < 0 (min form, normalized ray)");
      bool feasOk = !claim.x.empty() && checkPrimal(claim.x, false);
      if (claim.x.empty()) addCheck(cert, "feasible_point", Verdict::Warn, 0, 0, "no feasible point supplied");
      cert.accepted = rayOk && (claim.x.empty() || feasOk);
      if (!cert.accepted) cert.certified = rayOk ? Status::InfeasibleOrUnbounded : Status::NumericalFailure;
      break;
    }
    default:
      cert.accepted = true;  // no claim to certify
      break;
  }
  return cert;
}

}  // namespace pramana
