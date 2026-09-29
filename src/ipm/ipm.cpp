#include "ipm/ipm.h"

#include <algorithm>
#include <cmath>

#include "ipm/ldl.h"
#include "lp/scaling.h"
#include "lp/simplex.h"
#include "util/log.h"

namespace pramana {

namespace {
inline bool fin(double v) { return std::fabs(v) < kInfBoundThreshold; }

bool allFinite(const std::vector<double>& v) {
  for (double t : v)
    if (!std::isfinite(t)) return false;
  return true;
}

double infNorm(const std::vector<double>& v) {
  double m = 0;
  for (double x : v) m = std::max(m, std::fabs(x));
  return m;
}
}  // namespace

// ---------------------------------------------------------------------------
// Convexity check
// ---------------------------------------------------------------------------
ConvexityResult checkConvexity(const SparseMatrix& Q) {
  ConvexityResult res;
  const int n = Q.numCols;
  if (Q.nnz() == 0) return res;
  double maxDiag = 0;
  for (int j = 0; j < n; ++j)
    for (int k = Q.start[j]; k < Q.start[j + 1]; ++k)
      if (Q.index[k] == j) maxDiag = std::max(maxDiag, std::fabs(Q.value[k]));
  const double shift = 1e-12 * std::max(1.0, maxDiag);
  // Q + shift I with explicit diagonal.
  std::vector<int> ri, ci;
  std::vector<double> v;
  for (int j = 0; j < n; ++j) {
    ri.push_back(j);
    ci.push_back(j);
    v.push_back(shift);
    for (int k = Q.start[j]; k < Q.start[j + 1]; ++k) {
      ri.push_back(Q.index[k]);
      ci.push_back(j);
      v.push_back(Q.value[k]);
    }
  }
  SparseMatrix K = SparseMatrix::fromTriplets(n, n, ri, ci, v);
  LdlFactor f;
  f.analyze(K, minimumDegreeOrdering(K));
  std::vector<int8_t> signs(n, 1);
  // reg = -inf style: accept any positive pivot; negative pivots get replaced but counted.
  f.factorize(K, signs, 0.0);
  std::vector<double> d = f.pivots();
  const double tol = 1e-9 * std::max(1.0, maxDiag);
  int worst = -1;
  double worstVal = -tol;
  for (int k = 0; k < n; ++k) {
    if (d[k] < worstVal) {
      worstVal = d[k];
      worst = k;
    }
    if (d[k] < -tol) res.negativePivots++;
  }
  if (worst < 0) return res;
  // Certificate direction v = P^T L^{-T} e_k  =>  v^T (Q + shift I) v = d_k < 0.
  std::vector<double> rhs;
  f.solveUnitLt(worst, rhs);
  double vQv = 0, vv = 0;
  std::vector<double> Qv(n, 0.0);
  Q.multiply(rhs.data(), Qv.data());
  for (int j = 0; j < n; ++j) {
    vQv += rhs[j] * Qv[j];
    vv += rhs[j] * rhs[j];
  }
  if (vv > 0 && vQv < -1e-9 * std::max(1.0, maxDiag) * vv) {
    res.convex = false;
    res.curvature = vQv / vv;
    res.certificate = rhs;
  }
  return res;
}

// ---------------------------------------------------------------------------
// IPM
// ---------------------------------------------------------------------------
IpmResult solveIpm(const Model& model, const IpmOptions& opts, const Deadline* deadline) {
  Timer timer;
  IpmResult out;
  const int n0 = model.numCols(), m0 = model.numRows();
  const bool isQp = model.isQp();

  // ---- Build reduced problem: free (non-fixed) columns + slacks for inequality rows ----
  std::vector<int> colMap(n0, -1);  // original col -> reduced var
  std::vector<int> freeCols;
  std::vector<double> fixedVal(n0, 0.0);
  for (int j = 0; j < n0; ++j) {
    if (model.colLower[j] == model.colUpper[j]) {
      fixedVal[j] = model.colLower[j];
    } else {
      colMap[j] = static_cast<int>(freeCols.size());
      freeCols.push_back(j);
    }
  }
  std::vector<int> rowMap(m0, -1), keptRows, slackOf;
  for (int i = 0; i < m0; ++i) {
    if (!fin(model.rowLower[i]) && !fin(model.rowUpper[i])) continue;
    rowMap[i] = static_cast<int>(keptRows.size());
    keptRows.push_back(i);
  }
  const int nx = static_cast<int>(freeCols.size());
  const int m = static_cast<int>(keptRows.size());
  std::vector<int> slackRow;  // reduced var index nx+s -> kept row
  std::vector<int> rowSlack(m, -1);
  for (int r = 0; r < m; ++r) {
    int i = keptRows[r];
    if (model.rowLower[i] != model.rowUpper[i]) {
      rowSlack[r] = nx + static_cast<int>(slackRow.size());
      slackRow.push_back(r);
    }
  }
  const int n = nx + static_cast<int>(slackRow.size());

  std::vector<double> c(n, 0.0), lo(n), up(n), b(m, 0.0);
  double constant = model.objOffset;
  for (int j = 0; j < n0; ++j)
    if (colMap[j] < 0) constant += model.colCost[j] * fixedVal[j];
  for (int k = 0; k < nx; ++k) {
    int j = freeCols[k];
    c[k] = model.colCost[j];
    lo[k] = model.colLower[j];
    up[k] = model.colUpper[j];
  }
  for (size_t s = 0; s < slackRow.size(); ++s) {
    int i = keptRows[slackRow[s]];
    lo[nx + s] = model.rowLower[i];
    up[nx + s] = model.rowUpper[i];
  }
  // Q: reduced block + linear term from fixed columns + constant.
  std::vector<int> qi, qj;
  std::vector<double> qv;
  if (isQp) {
    for (int j = 0; j < n0; ++j)
      for (int k = model.Q.start[j]; k < model.Q.start[j + 1]; ++k) {
        int i = model.Q.index[k];
        double v = model.Q.value[k];
        if (colMap[i] >= 0 && colMap[j] >= 0) {
          qi.push_back(colMap[i]);
          qj.push_back(colMap[j]);
          qv.push_back(v);
        } else if (colMap[i] >= 0 && colMap[j] < 0) {
          c[colMap[i]] += v * fixedVal[j];
        } else if (colMap[i] < 0 && colMap[j] < 0) {
          constant += 0.5 * v * fixedVal[i] * fixedVal[j];
        }
      }
  }
  // A-hat with slack columns (-1).
  std::vector<int> ai, aj;
  std::vector<double> av;
  for (int j = 0; j < n0; ++j)
    for (int k = model.A.start[j]; k < model.A.start[j + 1]; ++k) {
      int r = rowMap[model.A.index[k]];
      if (r < 0) continue;
      if (colMap[j] >= 0) {
        ai.push_back(r);
        aj.push_back(colMap[j]);
        av.push_back(model.A.value[k]);
      } else {
        b[r] -= model.A.value[k] * fixedVal[j];
      }
    }
  for (int r = 0; r < m; ++r) {
    int i = keptRows[r];
    if (rowSlack[r] >= 0) {
      ai.push_back(r);
      aj.push_back(rowSlack[r]);
      av.push_back(-1.0);
    } else {
      b[r] += model.rowLower[i];
    }
  }
  SparseMatrix A = SparseMatrix::fromTriplets(m, n, ai, aj, av);
  SparseMatrix Q = SparseMatrix::fromTriplets(n, n, qi, qj, qv);

  // ---- Scaling (powers of two: exact) ----
  std::vector<double> R(m, 1.0), C(n, 1.0);
  if (opts.scale && A.nnz() > 0) {
    ScalingResult s = geometricScaling(A);
    R = s.rowScale;
    C = s.colScale;
    applyScaling(A, s);
  }
  for (int j = 0; j < n; ++j) {
    c[j] *= C[j];
    if (fin(lo[j])) lo[j] /= C[j]; else lo[j] = -kInf;
    if (fin(up[j])) up[j] /= C[j]; else up[j] = kInf;
  }
  for (int r = 0; r < m; ++r) b[r] *= R[r];
  for (int j = 0; j < n; ++j)
    for (int k = Q.start[j]; k < Q.start[j + 1]; ++k) Q.value[k] *= C[Q.index[k]] * C[j];
  SparseMatrix AT = A.transpose();

  // ---- KKT pattern: [-(Q+D) A^T; A dI], full symmetric CSC with explicit diagonal ----
  const int N = n + m;
  std::vector<int> ki, kj;
  std::vector<double> kv;
  for (int j = 0; j < n; ++j) {
    ki.push_back(j);
    kj.push_back(j);
    kv.push_back(-1.0);
    for (int k = Q.start[j]; k < Q.start[j + 1]; ++k)
      if (Q.index[k] != j) {
        ki.push_back(Q.index[k]);
        kj.push_back(j);
        kv.push_back(-Q.value[k]);
      }
    for (int k = A.start[j]; k < A.start[j + 1]; ++k) {
      ki.push_back(n + A.index[k]);
      kj.push_back(j);
      kv.push_back(A.value[k]);
      ki.push_back(j);
      kj.push_back(n + A.index[k]);
      kv.push_back(A.value[k]);
    }
  }
  for (int r = 0; r < m; ++r) {
    ki.push_back(n + r);
    kj.push_back(n + r);
    kv.push_back(1.0);
  }
  SparseMatrix K = SparseMatrix::fromTriplets(N, N, ki, kj, kv);
  std::vector<int> diagPos(N, -1);
  for (int j = 0; j < N; ++j)
    for (int k = K.start[j]; k < K.start[j + 1]; ++k)
      if (K.index[k] == j) diagPos[j] = k;
  std::vector<double> qDiag(n, 0.0);
  for (int j = 0; j < n; ++j)
    for (int k = Q.start[j]; k < Q.start[j + 1]; ++k)
      if (Q.index[k] == j) qDiag[j] = Q.value[k];
  std::vector<int8_t> signs(N, 1);
  for (int j = 0; j < n; ++j) signs[j] = -1;

  Timer ordT;
  LdlFactor ldl;
  ldl.analyze(K, minimumDegreeOrdering(K));
  out.stats.orderingSeconds = ordT.seconds();
  out.stats.nnzL = ldl.nnzL();

  // Unregularized copy for iterative refinement.
  SparseMatrix Kexact = K;
  auto setDiag = [&](const std::vector<double>& theta, double rho, double delta) {
    for (int j = 0; j < n; ++j) {
      K.value[diagPos[j]] = -(qDiag[j] + theta[j] + rho);
      Kexact.value[diagPos[j]] = -(qDiag[j] + theta[j]);
    }
    for (int r = 0; r < m; ++r) {
      K.value[diagPos[n + r]] = delta;
      Kexact.value[diagPos[n + r]] = 0.0;
    }
  };
  auto kktSolve = [&](std::vector<double>& rhs) {
    std::vector<double> sol = rhs;
    ldl.solve(sol);
    for (int it = 0; it < 3; ++it) {
      std::vector<double> res = rhs;
      // res -= Kexact * sol
      for (int j = 0; j < N; ++j) {
        double sj = sol[j];
        if (sj == 0.0) continue;
        for (int k = Kexact.start[j]; k < Kexact.start[j + 1]; ++k) res[Kexact.index[k]] -= Kexact.value[k] * sj;
      }
      double rn = infNorm(res), sn = infNorm(rhs);
      if (rn <= 1e-14 * (1.0 + sn)) break;
      ldl.solve(res);
      for (int j = 0; j < N; ++j) sol[j] += res[j];
    }
    rhs = sol;
  };

  // ---- Starting point ----
  std::vector<double> x(n, 0.0), y(m, 0.0), zl(n, 0.0), zu(n, 0.0);
  std::vector<char> hasL(n), hasU(n);
  for (int j = 0; j < n; ++j) {
    hasL[j] = fin(lo[j]);
    hasU[j] = fin(up[j]);
  }
  {
    std::vector<double> theta(n, 1.0);
    setDiag(theta, opts.primalReg, opts.dualReg);
    out.stats.regularizedPivots += ldl.factorize(K, signs, 1e-10);
    // min 1/2||x - x0||^2 + c^T x  s.t. Ax = b   (a regularized least-squares start)
    std::vector<double> rhs(N, 0.0);
    for (int j = 0; j < n; ++j) {
      double x0 = 0;
      if (hasL[j] && hasU[j]) x0 = 0.5 * (lo[j] + up[j]);
      else if (hasL[j]) x0 = lo[j];
      else if (hasU[j]) x0 = up[j];
      rhs[j] = -x0;
    }
    for (int r = 0; r < m; ++r) rhs[n + r] = b[r];
    std::vector<double> rhsX0(rhs.begin(), rhs.begin() + n);
    kktSolve(rhs);
    double scaleRef = 1.0 + infNorm(b);
    for (int j = 0; j < n; ++j) {
      if (hasL[j]) scaleRef = std::max(scaleRef, std::fabs(lo[j]));
      if (hasU[j]) scaleRef = std::max(scaleRef, std::fabs(up[j]));
    }
    bool lsOk = allFinite(rhs) && infNorm(rhs) < 1e8 * scaleRef;
    for (int j = 0; j < n; ++j) x[j] = lsOk ? rhs[j] : -rhsX0[j];
    for (int r = 0; r < m; ++r) y[r] = lsOk ? rhs[n + r] : 0.0;
    // Push into the interior.
    double shift = 0;
    for (int j = 0; j < n; ++j) {
      if (hasL[j] && hasU[j]) continue;
      if (hasL[j]) shift = std::max(shift, lo[j] - x[j]);
      if (hasU[j]) shift = std::max(shift, x[j] - up[j]);
    }
    shift = std::max(1.0, 1.5 * shift);
    for (int j = 0; j < n; ++j) {
      if (hasL[j] && hasU[j]) {
        double w = up[j] - lo[j];
        double mg = std::min(0.5 * w, std::max(0.1 * w, std::min(1.0, 0.5 * w)));
        x[j] = std::min(std::max(x[j], lo[j] + mg), up[j] - mg);
      } else if (hasL[j]) {
        x[j] = std::max(x[j], lo[j] + shift);
      } else if (hasU[j]) {
        x[j] = std::min(x[j], up[j] - shift);
      }
    }
    // Duals: z from the dual residual, kept positive.
    std::vector<double> rd = c;
    Q.multiply(x.data(), rd.data());
    AT.multiply(y.data(), rd.data(), -1.0);  // rd = c + Qx - A^T y
    double zshift = 1.0;
    for (int j = 0; j < n; ++j) zshift = std::max(zshift, std::fabs(rd[j]));
    zshift = std::min(zshift, 1e4);
    for (int j = 0; j < n; ++j) {
      if (hasL[j]) zl[j] = std::max(1.0, rd[j] > 0 ? rd[j] : 0.0) + 0.0 * zshift;
      if (hasU[j]) zu[j] = std::max(1.0, rd[j] < 0 ? -rd[j] : 0.0);
    }
  }

  // ---- Main loop ----
  const double bnorm = infNorm(b), cnorm = infNorm(c);
  int nComp = 0;
  for (int j = 0; j < n; ++j) nComp += hasL[j] + hasU[j];
  std::vector<double> rp(m), rd(n), theta(n), dx(N), dzl(n), dzu(n), rhs(N);
  std::vector<double> dxa(n), dzla(n), dzua(n);
  Status status = Status::IterationLimit;
  double pobj = 0, dobj = 0;
  double lastStep = 1.0;
  // Best iterate by KKT merit (restored on breakdown / stagnation).
  double bestMerit = kInf;
  int bestIter = 0;
  std::vector<double> bx, by, bzl, bzu;

  for (int iter = 0; iter < opts.maxIterations; ++iter) {
    if (deadline && deadline->expired()) {
      status = Status::TimeLimit;
      break;
    }
    // Residuals.
    std::fill(rp.begin(), rp.end(), 0.0);
    A.multiply(x.data(), rp.data(), -1.0);
    for (int r = 0; r < m; ++r) rp[r] += b[r];
    std::vector<double> Qx(n, 0.0);
    Q.multiply(x.data(), Qx.data());
    for (int j = 0; j < n; ++j) rd[j] = c[j] + Qx[j] - zl[j] + zu[j];
    AT.multiply(y.data(), rd.data(), -1.0);
    double mu = 0;
    for (int j = 0; j < n; ++j) {
      if (hasL[j]) mu += (x[j] - lo[j]) * zl[j];
      if (hasU[j]) mu += (up[j] - x[j]) * zu[j];
    }
    mu = nComp > 0 ? mu / nComp : 0.0;
    // Objectives (scaled space == unscaled values).
    double xQx = 0, cx = 0;
    for (int j = 0; j < n; ++j) {
      xQx += x[j] * Qx[j];
      cx += c[j] * x[j];
    }
    pobj = cx + 0.5 * xQx + constant;
    dobj = -0.5 * xQx + constant;
    for (int r = 0; r < m; ++r) dobj += b[r] * y[r];
    for (int j = 0; j < n; ++j) {
      if (hasL[j]) dobj += lo[j] * zl[j];
      if (hasU[j]) dobj -= up[j] * zu[j];
    }
    // Unscaled relative residuals.
    double pres = 0, dres = 0;
    for (int r = 0; r < m; ++r) pres = std::max(pres, std::fabs(rp[r]) / R[r]);
    for (int j = 0; j < n; ++j) dres = std::max(dres, std::fabs(rd[j]) / C[j]);
    double bn = 0, cn = 0;
    for (int r = 0; r < m; ++r) bn = std::max(bn, std::fabs(b[r]) / R[r]);
    for (int j = 0; j < n; ++j) cn = std::max(cn, std::fabs(c[j]) / C[j]);
    double relP = pres / (1.0 + bn), relD = dres / (1.0 + cn);
    double relGap = std::fabs(pobj - dobj) / (1.0 + std::fabs(pobj));
    out.stats.iterations = iter;
    out.stats.finalPrimalResidual = relP;
    out.stats.finalDualResidual = relD;
    out.stats.finalGap = relGap;
    out.stats.finalMu = mu;
    if (opts.logLevel > 0)
      PLOG_INFO("  ipm %3d  pobj % .10e  dobj % .10e  pres %.2e  dres %.2e  gap %.2e  mu %.2e", iter, pobj, dobj,
                relP, relD, relGap, mu);
    const double merit = std::max(relP, std::max(relD, relGap));
    if (!std::isfinite(merit)) {
      status = Status::NumericalFailure;
      break;
    }
    if (merit < bestMerit) {
      bestMerit = merit;
      bestIter = iter;
      bx = x;
      by = y;
      bzl = zl;
      bzu = zu;
    }
    if (relP <= opts.tolerance && relD <= opts.tolerance && relGap <= opts.tolerance) {
      status = Status::Optimal;
      break;
    }
    if (iter - bestIter > 30) {  // stagnation
      status = Status::IterationLimit;
      break;
    }
    // Divergence => infeasible or unbounded (no certificate from the IPM itself).
    if (iter > 5 && (infNorm(x) > 1e14 || infNorm(y) > 1e14)) {
      status = Status::InfeasibleOrUnbounded;
      break;
    }
    // Theta^{-1}.
    for (int j = 0; j < n; ++j) {
      double t = 0;
      if (hasL[j]) t += zl[j] / (x[j] - lo[j]);
      if (hasU[j]) t += zu[j] / (up[j] - x[j]);
      theta[j] = std::min(t, 1e30);
    }
    double rho = std::max(opts.primalReg, 1e-12), delta = std::max(opts.dualReg, 1e-12);
    setDiag(theta, rho, delta);
    Timer ft;
    out.stats.regularizedPivots += ldl.factorize(K, signs, 1e-11);
    out.stats.factorSeconds += ft.seconds();

    auto solveDir = [&](const std::vector<double>& rcl, const std::vector<double>& rcu) {
      for (int j = 0; j < n; ++j) {
        double v = rd[j];
        if (hasL[j]) v -= rcl[j] / (x[j] - lo[j]);
        if (hasU[j]) v += rcu[j] / (up[j] - x[j]);
        rhs[j] = v;
      }
      for (int r = 0; r < m; ++r) rhs[n + r] = rp[r];
      Timer st;
      kktSolve(rhs);
      out.stats.solveSeconds += st.seconds();
      for (int r = 0; r < m; ++r) dx[n + r] = rhs[n + r];  // dy
      for (int j = 0; j < n; ++j) {
        dx[j] = rhs[j];
        dzl[j] = hasL[j] ? (rcl[j] - zl[j] * dx[j]) / (x[j] - lo[j]) : 0.0;
        dzu[j] = hasU[j] ? (rcu[j] + zu[j] * dx[j]) / (up[j] - x[j]) : 0.0;
      }
    };
    auto stepLengths = [&](double& ap, double& ad) {
      ap = 1.0;
      ad = 1.0;
      for (int j = 0; j < n; ++j) {
        if (hasL[j] && dx[j] < 0) ap = std::min(ap, -(x[j] - lo[j]) / dx[j]);
        if (hasU[j] && dx[j] > 0) ap = std::min(ap, (up[j] - x[j]) / dx[j]);
        if (hasL[j] && dzl[j] < 0) ad = std::min(ad, -zl[j] / dzl[j]);
        if (hasU[j] && dzu[j] < 0) ad = std::min(ad, -zu[j] / dzu[j]);
      }
    };
    // Predictor.
    std::vector<double> rcl(n, 0.0), rcu(n, 0.0);
    for (int j = 0; j < n; ++j) {
      if (hasL[j]) rcl[j] = -(x[j] - lo[j]) * zl[j];
      if (hasU[j]) rcu[j] = -(up[j] - x[j]) * zu[j];
    }
    solveDir(rcl, rcu);
    double apA, adA;
    stepLengths(apA, adA);
    if (isQp) apA = adA = std::min(apA, adA);
    double muAff = 0;
    for (int j = 0; j < n; ++j) {
      if (hasL[j]) muAff += (x[j] - lo[j] + apA * dx[j]) * (zl[j] + adA * dzl[j]);
      if (hasU[j]) muAff += (up[j] - x[j] - apA * dx[j]) * (zu[j] + adA * dzu[j]);
    }
    muAff = nComp > 0 ? muAff / nComp : 0.0;
    double sigma = mu > 0 ? std::pow(std::max(muAff, 0.0) / mu, 3) : 0.0;
    sigma = std::min(sigma, 1.0);
    // Centrality safeguard: after short steps, recenter instead of pushing mu down.
    if (lastStep < 0.1) sigma = std::max(sigma, 0.3);
    if (lastStep < 0.01) sigma = std::max(sigma, 0.7);
    for (int j = 0; j < n; ++j) {
      dxa[j] = dx[j];
      dzla[j] = dzl[j];
      dzua[j] = dzu[j];
    }
    // Corrector.
    for (int j = 0; j < n; ++j) {
      if (hasL[j]) rcl[j] = sigma * mu - (x[j] - lo[j]) * zl[j] - dxa[j] * dzla[j];
      if (hasU[j]) rcu[j] = sigma * mu - (up[j] - x[j]) * zu[j] + dxa[j] * dzua[j];
    }
    solveDir(rcl, rcu);
    if (!allFinite(dx) || !allFinite(dzl) || !allFinite(dzu)) {
      status = Status::NumericalFailure;
      break;
    }
    double ap, ad;
    stepLengths(ap, ad);
    const double eta = std::max(0.9, 1.0 - 10.0 * mu);
    ap = std::min(1.0, eta * ap);
    ad = std::min(1.0, eta * ad);
    if (isQp) ap = ad = std::min(ap, ad);
    lastStep = std::min(ap, ad);
    for (int j = 0; j < n; ++j) {
      x[j] += ap * dx[j];
      zl[j] += ad * dzl[j];
      zu[j] += ad * dzu[j];
    }
    for (int r = 0; r < m; ++r) y[r] += ad * dx[n + r];
    // Keep strict interior against round-off.
    for (int j = 0; j < n; ++j) {
      if (hasL[j]) {
        if (x[j] - lo[j] < 1e-14 * (1 + std::fabs(lo[j]))) x[j] = lo[j] + 1e-14 * (1 + std::fabs(lo[j]));
        zl[j] = std::max(zl[j], 1e-30);
      }
      if (hasU[j]) {
        if (up[j] - x[j] < 1e-14 * (1 + std::fabs(up[j]))) x[j] = up[j] - 1e-14 * (1 + std::fabs(up[j]));
        zu[j] = std::max(zu[j], 1e-30);
      }
    }
    (void)bnorm;
    (void)cnorm;
  }

  // Restore the best iterate unless we terminated at the optimum.
  if (status != Status::Optimal && !bx.empty()) {
    x = bx;
    y = by;
    zl = bzl;
    zu = bzu;
    if (bestMerit <= std::max(1e-6, opts.tolerance)) status = Status::Optimal;  // tolerance-level optimum
    else if (status == Status::NumericalFailure) status = Status::IterationLimit;
    out.stats.finalGap = bestMerit;
  }
  // ---- Unscale and map back ----
  out.x.assign(n0, 0.0);
  for (int j = 0; j < n0; ++j) out.x[j] = colMap[j] >= 0 ? x[colMap[j]] * C[colMap[j]] : fixedVal[j];
  out.rowDual.assign(m0, 0.0);
  for (int r = 0; r < m; ++r) out.rowDual[keptRows[r]] = y[r] * R[r];
  out.rowActivity = model.rowActivity(out.x);
  out.reducedCost.assign(n0, 0.0);
  {
    // d = c + Qx - A^T y on the original model (exact recomputation).
    std::vector<double> d = model.colCost;
    if (isQp) model.Q.multiply(out.x.data(), d.data());
    model.A.multiplyTranspose(out.rowDual.data(), d.data(), -1.0);
    out.reducedCost = d;
  }
  out.objective = model.objective(out.x);
  out.status = status;

  // ---- Crossover (LP): basis identification + simplex cleanup ----
  if (!isQp && opts.crossover && (status == Status::Optimal || status == Status::IterationLimit)) {
    std::vector<BasisStatus> cs(n0, BasisStatus::Lower), rs(m0, BasisStatus::Basic);
    struct Cand {
      double score;
      int idx;  // < n0 column, >= n0 row
    };
    std::vector<Cand> cands;
    for (int j = 0; j < n0; ++j) {
      double l = model.colLower[j], u = model.colUpper[j], v = out.x[j];
      if (l == u) {
        cs[j] = BasisStatus::Fixed;
        continue;
      }
      double dl = fin(l) ? v - l : kInf, du = fin(u) ? u - v : kInf;
      double zj = std::fabs(out.reducedCost[j]);
      double dist = std::min(dl, du) / (1.0 + std::fabs(v));
      // Interior-ness: distance to nearest bound relative to the dual value.
      double score = dist / (zj + 1e-12);
      if (!fin(l) && !fin(u)) score = kInf;
      cs[j] = dl <= du ? BasisStatus::Lower : BasisStatus::Upper;
      if (!fin(l) && !fin(u)) cs[j] = BasisStatus::Zero;
      cands.push_back({score, j});
    }
    for (int i = 0; i < m0; ++i) {
      double l = model.rowLower[i], u = model.rowUpper[i], v = out.rowActivity[i];
      double dl = fin(l) ? v - l : kInf, du = fin(u) ? u - v : kInf;
      double yi = std::fabs(out.rowDual[i]);
      double dist = std::min(dl, du) / (1.0 + std::fabs(v));
      double score = (!fin(l) && !fin(u)) ? kInf : dist / (yi + 1e-12);
      rs[i] = l == u ? BasisStatus::Fixed : (dl <= du ? BasisStatus::Lower : BasisStatus::Upper);
      cands.push_back({score, n0 + i});
    }
    std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b2) {
      if (a.score != b2.score) return a.score > b2.score;
      return a.idx > b2.idx;  // prefer logicals on ties (well-conditioned)
    });
    for (int k = 0; k < static_cast<int>(cands.size()); ++k) {
      int idx = cands[k].idx;
      if (k < m0) {
        if (idx < n0) cs[idx] = BasisStatus::Basic;
        else rs[idx - n0] = BasisStatus::Basic;
      }
    }
    SimplexOptions so;
    so.perturb = true;
    LpResult lr = solveLpSimplex(model, so, deadline, &cs, &rs);
    out.stats.crossoverDone = true;
    out.stats.crossoverIterations = static_cast<int>(lr.stats.iterations);
    if (lr.status == Status::Optimal || lr.status == Status::Infeasible || lr.status == Status::Unbounded) {
      out.status = lr.status;
      out.x = lr.x;
      out.rowActivity = lr.rowActivity;
      out.rowDual = lr.rowDual;
      out.reducedCost = lr.reducedCost;
      out.colStatus = lr.colStatus;
      out.rowStatus = lr.rowStatus;
      out.objective = lr.objective;
    }
  }
  (void)timer;
  return out;
}

}  // namespace pramana
