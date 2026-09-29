#include "pdhg/pdhg.h"

#include <algorithm>
#include <cmath>
#include <memory>

#include "gpu/cuda_driver.h"
#include "lp/scaling.h"
#include "util/log.h"
#include "util/parallel.h"

namespace pramana {

namespace {

inline bool fin(double v) { return std::fabs(v) < kInfBoundThreshold; }
inline double clampInf(double v) { return fin(v) ? v : (v > 0 ? HUGE_VAL : -HUGE_VAL); }

struct Csr {
  int rows = 0, cols = 0;
  std::vector<int> ptr, idx;
  std::vector<double> val;
};

Csr csrOf(const SparseMatrix& csc) {  // CSR of the matrix whose CSC is given
  SparseMatrix t = csc.transpose();
  Csr c;
  c.rows = csc.numRows;
  c.cols = csc.numCols;
  c.ptr = t.start;
  c.idx = t.index;
  c.val = t.value;
  return c;
}
Csr csrOfTranspose(const SparseMatrix& csc) {  // CSR of A^T == CSC of A
  Csr c;
  c.rows = csc.numCols;
  c.cols = csc.numRows;
  c.ptr = csc.start;
  c.idx = csc.index;
  c.val = csc.value;
  return c;
}

enum V : int { X, XP, X0, ATY, ATYP, C, LO, UP, Y, YP, Y0, AX, AXP, AX0, RL, RU, NV };
inline bool rowVec(int v) { return v >= Y; }

// ---------------------------------------------------------------------------
// Backend interface: the driver issues the same sequence on CPU or GPU.
// ---------------------------------------------------------------------------
class Backend {
 public:
  virtual ~Backend() = default;
  virtual void init(const Csr& A, const Csr& AT, const std::vector<double>& rfac, const std::vector<double>& cfac,
                    int k) = 0;
  virtual void set(int v, const std::vector<double>& data) = 0;
  virtual void get(int v, std::vector<double>& out) = 0;
  virtual void spmv(int in, int out) = 0;   // out(m) = A in(n)
  virtual void spmvT(int in, int out) = 0;  // out(n) = A^T in(m)
  virtual void setSteps(const std::vector<double>& tau, const std::vector<double>& sigma) = 0;
  virtual void primalStep() = 0;
  virtual void dualStep() = 0;
  virtual void halpern(double rho) = 0;
  virtual void restart(const std::vector<int>& mask) = 0;
  virtual void kkt(std::vector<double>& rows3k, std::vector<double>& cols3k) = 0;
  virtual void diff2(int a, int b, std::vector<double>& out) = 0;
  virtual void startTiming() {}
  virtual double stopTiming() { return 0; }
  int n = 0, m = 0, k = 1;
};

// ---------------------------------------------------------------------------
// CPU backend (multi-threaded, deterministic static chunking)
// ---------------------------------------------------------------------------
class CpuBackend : public Backend {
 public:
  void init(const Csr& A, const Csr& AT, const std::vector<double>& rfac, const std::vector<double>& cfac,
            int kk) override {
    A_ = &A;
    AT_ = &AT;
    rfac_ = rfac;
    cfac_ = cfac;
    n = A.cols;
    m = A.rows;
    k = kk;
    for (int v = 0; v < NV; ++v) vec_[v].assign(static_cast<size_t>(rowVec(v) ? m : n) * k, 0.0);
    iter_.assign(k, 0);
    tau_.assign(k, 0);
    sigma_.assign(k, 0);
  }
  void set(int v, const std::vector<double>& d) override { vec_[v] = d; }
  void get(int v, std::vector<double>& out) override { out = vec_[v]; }
  void mult(const Csr& M, const std::vector<double>& in, std::vector<double>& out) {
    const int kk = k;
    ThreadPool::instance().parallelFor(M.rows, [&](int, int b, int e) {
      for (int i = b; i < e; ++i)
        for (int q = 0; q < kk; ++q) {
          double s = 0;
          for (int p = M.ptr[i]; p < M.ptr[i + 1]; ++p) s += M.val[p] * in[static_cast<size_t>(M.idx[p]) * kk + q];
          out[static_cast<size_t>(i) * kk + q] = s;
        }
    }, 512);
  }
  void spmv(int in, int out) override { mult(*A_, vec_[in], vec_[out]); }
  void spmvT(int in, int out) override { mult(*AT_, vec_[in], vec_[out]); }
  void setSteps(const std::vector<double>& t, const std::vector<double>& s) override {
    tau_ = t;
    sigma_ = s;
  }
  template <class F>
  void elementwise(size_t len, F f) {
    ThreadPool::instance().parallelFor(static_cast<int>(len), [&](int, int b, int e) {
      for (int t = b; t < e; ++t) f(static_cast<size_t>(t));
    }, 8192);
  }
  void primalStep() override {
    auto &x = vec_[X], &aty = vec_[ATY], &c = vec_[C], &lo = vec_[LO], &up = vec_[UP], &xp = vec_[XP];
    const int kk = k;
    elementwise(x.size(), [&](size_t t) {
      double v = x[t] - tau_[t % kk] * (c[t] - aty[t]);
      v = std::max(v, lo[t]);
      v = std::min(v, up[t]);
      xp[t] = v;
    });
  }
  void dualStep() override {
    auto &y = vec_[Y], &axp = vec_[AXP], &ax = vec_[AX], &rl = vec_[RL], &ru = vec_[RU], &yp = vec_[YP];
    const int kk = k;
    elementwise(y.size(), [&](size_t t) {
      double s = sigma_[t % kk];
      double v = y[t] - s * (2.0 * axp[t] - ax[t]);
      double w = -v / s;
      w = std::max(w, rl[t]);
      w = std::min(w, ru[t]);
      yp[t] = v + s * w;
    });
  }
  void halpernOne(int z, int zp, int z0, double rho) {
    auto &a = vec_[z], &p = vec_[zp], &o = vec_[z0];
    const int kk = k;
    elementwise(a.size(), [&](size_t t) {
      double kk2 = iter_[t % kk];
      double lam = (kk2 + 1.0) / (kk2 + 2.0);
      a[t] = lam * ((1.0 + rho) * p[t] - rho * a[t]) + (1.0 - lam) * o[t];
    });
  }
  void halpern(double rho) override {
    halpernOne(X, XP, X0, rho);
    halpernOne(Y, YP, Y0, rho);
    halpernOne(AX, AXP, AX0, rho);
    for (auto& it : iter_) ++it;
  }
  void restart(const std::vector<int>& mask) override {
    const int kk = k;
    auto cp = [&](int src, int d1, int d2) {
      auto &s = vec_[src], &a = vec_[d1], &b = vec_[d2];
      for (size_t t = 0; t < s.size(); ++t)
        if (mask[t % kk]) a[t] = b[t] = s[t];
    };
    cp(XP, X, X0);
    cp(YP, Y, Y0);
    cp(AXP, AX, AX0);
    for (int b = 0; b < k; ++b)
      if (mask[b]) iter_[b] = 0;
  }
  void kkt(std::vector<double>& r3, std::vector<double>& c3) override {
    r3.assign(3 * k, 0.0);
    c3.assign(3 * k, 0.0);
    auto &axp = vec_[AXP], &rl = vec_[RL], &ru = vec_[RU], &yp = vec_[YP];
    for (int i = 0; i < m; ++i)
      for (int b = 0; b < k; ++b) {
        size_t t = static_cast<size_t>(i) * k + b;
        double a = axp[t];
        double res = (std::max(rl[t] - a, 0.0) + std::max(a - ru[t], 0.0)) * rfac_[i];
        r3[3 * b] += res * res;
        double y = yp[t];
        if (y > 0 && rl[t] > -1e300) r3[3 * b + 1] += y * rl[t];
        else if (y < 0 && ru[t] < 1e300) r3[3 * b + 1] += y * ru[t];
        r3[3 * b + 2] += y * y;
      }
    auto &xp = vec_[XP], &c = vec_[C], &atyp = vec_[ATYP], &lo = vec_[LO], &up = vec_[UP];
    for (int j = 0; j < n; ++j)
      for (int b = 0; b < k; ++b) {
        size_t t = static_cast<size_t>(j) * k + b;
        double d = c[t] - atyp[t];
        bool hl = lo[t] > -1e300, hu = up[t] < 1e300;
        double res = 0;
        if (!hl && !hu) res = std::fabs(d);
        else if (!hl) res = std::max(d, 0.0);
        else if (!hu) res = std::max(-d, 0.0);
        res *= cfac_[j];
        c3[3 * b] += res * res;
        c3[3 * b + 1] += c[t] * xp[t];
        if (d > 0 && hl) c3[3 * b + 2] += d * lo[t];
        else if (d < 0 && hu) c3[3 * b + 2] += d * up[t];
      }
  }
  void diff2(int a, int b, std::vector<double>& out) override {
    out.assign(k, 0.0);
    auto &va = vec_[a], &vb = vec_[b];
    for (size_t t = 0; t < va.size(); ++t) {
      double d = va[t] - vb[t];
      out[t % k] += d * d;
    }
  }

 private:
  const Csr* A_ = nullptr;
  const Csr* AT_ = nullptr;
  std::vector<double> rfac_, cfac_;
  std::vector<double> vec_[NV];
  std::vector<int> iter_;
  std::vector<double> tau_, sigma_;
};

// ---------------------------------------------------------------------------
// CUDA backend (own kernels, driver API)
// ---------------------------------------------------------------------------
class GpuBackend : public Backend {
 public:
  GpuBackend() : ctx_(gpu::CudaContext::instance()) {}
  ~GpuBackend() override {
    for (auto p : bufs_) ctx_.free(p);
    if (evA_) ctx_.destroyEvent(evA_);
    if (evB_) ctx_.destroyEvent(evB_);
  }
  gpu::CUdeviceptr buf(size_t bytes) {
    gpu::CUdeviceptr p = ctx_.alloc(bytes);
    bufs_.push_back(p);
    return p;
  }
  template <class T>
  gpu::CUdeviceptr upload(const std::vector<T>& v) {
    gpu::CUdeviceptr p = buf(v.size() * sizeof(T));
    ctx_.upload(p, v.data(), v.size() * sizeof(T));
    return p;
  }
  void init(const Csr& A, const Csr& AT, const std::vector<double>& rfac, const std::vector<double>& cfac,
            int kk) override {
    n = A.cols;
    m = A.rows;
    k = kk;
    aPtr_ = upload(A.ptr);
    aIdx_ = upload(A.idx);
    aVal_ = upload(A.val);
    tPtr_ = upload(AT.ptr);
    tIdx_ = upload(AT.idx);
    tVal_ = upload(AT.val);
    avgA_ = m > 0 ? static_cast<double>(A.val.size()) / m : 0;
    avgT_ = n > 0 ? static_cast<double>(AT.val.size()) / n : 0;
    rfac_ = upload(rfac.empty() ? std::vector<double>(1, 0.0) : rfac);
    cfac_ = upload(cfac.empty() ? std::vector<double>(1, 0.0) : cfac);
    for (int v = 0; v < NV; ++v) {
      size_t len = static_cast<size_t>(rowVec(v) ? m : n) * k;
      d_[v] = buf(std::max<size_t>(len, 1) * sizeof(double));
      ctx_.memsetZero(d_[v], len * sizeof(double));
    }
    iter_ = buf(k * sizeof(int));
    ctx_.memsetZero(iter_, k * sizeof(int));
    mask_ = buf(k * sizeof(int));
    tau_ = buf(k * sizeof(double));
    sigma_ = buf(k * sizeof(double));
    partial_ = buf(static_cast<size_t>(kRedBlocks) * k * 3 * sizeof(double));
    fSpmm_ = ctx_.function("spmm_thread");
    fSpmvWarp_ = ctx_.function("spmv_warp");
    fPrimal_ = ctx_.function("pdhg_primal");
    fDual_ = ctx_.function("pdhg_dual");
    fHalpern_ = ctx_.function("halpern");
    fAdvance_ = ctx_.function("advance_iter");
    fRestart_ = ctx_.function("restart_copy");
    fResetIter_ = ctx_.function("reset_iter");
    fDiff2_ = ctx_.function("reduce_diff2");
    fKktRows_ = ctx_.function("kkt_rows");
    fKktCols_ = ctx_.function("kkt_cols");
    evA_ = ctx_.createEvent();
    evB_ = ctx_.createEvent();
  }
  static unsigned blocks(long long len) { return static_cast<unsigned>((len + 255) / 256); }
  void set(int v, const std::vector<double>& data) override {
    ctx_.upload(d_[v], data.data(), data.size() * sizeof(double));
  }
  void get(int v, std::vector<double>& out) override {
    out.resize(static_cast<size_t>(rowVec(v) ? m : n) * k);
    ctx_.download(out.data(), d_[v], out.size() * sizeof(double));
  }
  void mult(int rows, gpu::CUdeviceptr ptr, gpu::CUdeviceptr idx, gpu::CUdeviceptr val, double avg, int in, int out) {
    if (rows == 0) return;
    if (k == 1 && avg >= 16) {
      void* args[] = {&rows, &ptr, &idx, &val, &d_[in], &d_[out]};
      ctx_.launch(fSpmvWarp_, blocks(static_cast<long long>(rows) * 32), 1, 256, args);
    } else {
      void* args[] = {&rows, &k, &ptr, &idx, &val, &d_[in], &d_[out]};
      ctx_.launch(fSpmm_, blocks(static_cast<long long>(rows) * k), 1, 256, args);
    }
  }
  void spmv(int in, int out) override { mult(m, aPtr_, aIdx_, aVal_, avgA_, in, out); }
  void spmvT(int in, int out) override { mult(n, tPtr_, tIdx_, tVal_, avgT_, in, out); }
  void setSteps(const std::vector<double>& t, const std::vector<double>& s) override {
    ctx_.upload(tau_, t.data(), k * sizeof(double));
    ctx_.upload(sigma_, s.data(), k * sizeof(double));
  }
  void primalStep() override {
    void* args[] = {&n, &k, &d_[X], &d_[ATY], &d_[C], &d_[LO], &d_[UP], &tau_, &d_[XP]};
    ctx_.launch(fPrimal_, blocks(static_cast<long long>(n) * k), 1, 256, args);
  }
  void dualStep() override {
    void* args[] = {&m, &k, &d_[Y], &d_[AXP], &d_[AX], &d_[RL], &d_[RU], &sigma_, &d_[YP]};
    ctx_.launch(fDual_, blocks(static_cast<long long>(m) * k), 1, 256, args);
  }
  void halpern(double rho) override {
    auto one = [&](int len, int z, int zp, int z0) {
      void* args[] = {&len, &k, &iter_, &rho, &d_[z], &d_[zp], &d_[z0]};
      ctx_.launch(fHalpern_, blocks(static_cast<long long>(len) * k), 1, 256, args);
    };
    one(n, X, XP, X0);
    one(m, Y, YP, Y0);
    one(m, AX, AXP, AX0);
    void* args[] = {&k, &iter_};
    ctx_.launch(fAdvance_, blocks(k), 1, 256, args);
  }
  void restart(const std::vector<int>& mask) override {
    ctx_.upload(mask_, mask.data(), k * sizeof(int));
    auto cp = [&](int len, int src, int d1, int d2) {
      void* args[] = {&len, &k, &mask_, &d_[src], &d_[d1], &d_[d2]};
      ctx_.launch(fRestart_, blocks(static_cast<long long>(len) * k), 1, 256, args);
    };
    cp(n, XP, X, X0);
    cp(m, YP, Y, Y0);
    cp(m, AXP, AX, AX0);
    void* args[] = {&k, &mask_, &iter_};
    ctx_.launch(fResetIter_, blocks(k), 1, 256, args);
  }
  unsigned redBlocks(long long len) const {
    return static_cast<unsigned>(std::max<long long>(1, std::min<long long>(kRedBlocks, (len + 255) / 256)));
  }
  void kkt(std::vector<double>& r3, std::vector<double>& c3) override {
    r3.assign(3 * k, 0.0);
    c3.assign(3 * k, 0.0);
    std::vector<double> host;
    if (m > 0) {
      unsigned g = redBlocks(m);
      void* args[] = {&m, &k, &d_[AXP], &d_[RL], &d_[RU], &d_[YP], &rfac_, &partial_};
      ctx_.launch(fKktRows_, g, k, 256, args);
      host.resize(static_cast<size_t>(g) * k * 3);
      ctx_.download(host.data(), partial_, host.size() * sizeof(double));
      for (int b = 0; b < k; ++b)
        for (unsigned q = 0; q < g; ++q)
          for (int s = 0; s < 3; ++s) r3[3 * b + s] += host[(static_cast<size_t>(b) * g + q) * 3 + s];
    }
    unsigned g = redBlocks(n);
    void* args[] = {&n, &k, &d_[XP], &d_[C], &d_[ATYP], &d_[LO], &d_[UP], &cfac_, &partial_};
    ctx_.launch(fKktCols_, g, k, 256, args);
    host.resize(static_cast<size_t>(g) * k * 3);
    ctx_.download(host.data(), partial_, host.size() * sizeof(double));
    for (int b = 0; b < k; ++b)
      for (unsigned q = 0; q < g; ++q)
        for (int s = 0; s < 3; ++s) c3[3 * b + s] += host[(static_cast<size_t>(b) * g + q) * 3 + s];
  }
  void diff2(int a, int b, std::vector<double>& out) override {
    out.assign(k, 0.0);
    int len = rowVec(a) ? m : n;
    if (len == 0) return;
    unsigned g = redBlocks(len);
    void* args[] = {&len, &k, &d_[a], &d_[b], &partial_};
    ctx_.launch(fDiff2_, g, k, 256, args);
    std::vector<double> host(static_cast<size_t>(g) * k);
    ctx_.download(host.data(), partial_, host.size() * sizeof(double));
    for (int q = 0; q < k; ++q)
      for (unsigned p = 0; p < g; ++p) out[q] += host[static_cast<size_t>(q) * g + p];
  }
  void startTiming() override { ctx_.record(evA_); }
  double stopTiming() override {
    ctx_.record(evB_);
    return ctx_.elapsedMs(evA_, evB_) / 1000.0;
  }

 private:
  static constexpr int kRedBlocks = 256;
  gpu::CudaContext& ctx_;
  std::vector<gpu::CUdeviceptr> bufs_;
  gpu::CUdeviceptr d_[NV] = {};
  gpu::CUdeviceptr aPtr_ = 0, aIdx_ = 0, aVal_ = 0, tPtr_ = 0, tIdx_ = 0, tVal_ = 0;
  gpu::CUdeviceptr rfac_ = 0, cfac_ = 0, iter_ = 0, mask_ = 0, tau_ = 0, sigma_ = 0, partial_ = 0;
  double avgA_ = 0, avgT_ = 0;
  gpu::CUfunction fSpmm_ = nullptr, fSpmvWarp_ = nullptr, fPrimal_ = nullptr, fDual_ = nullptr, fHalpern_ = nullptr,
                  fAdvance_ = nullptr, fRestart_ = nullptr, fResetIter_ = nullptr, fDiff2_ = nullptr,
                  fKktRows_ = nullptr, fKktCols_ = nullptr;
  gpu::CUevent evA_ = nullptr, evB_ = nullptr;
};

// ---------------------------------------------------------------------------
// Problem preparation (scaling) shared by single and batched solves.
// ---------------------------------------------------------------------------
struct Prepared {
  int n = 0, m = 0;
  Csr A, AT;
  std::vector<double> R, Cs;   // Ruiz factors
  double beta = 1, gamma = 1;  // bound / objective rescaling
  std::vector<double> rfac, cfac;
  double normA = 1;
  double offset = 0;
};

double powerIteration(const SparseMatrix& A, int iters) {
  const int n = A.numCols, m = A.numRows;
  if (A.nnz() == 0) return 1.0;
  std::vector<double> x(n, 1.0 / std::sqrt(std::max(1, n))), y(m), z(n);
  double s = 1;
  for (int it = 0; it < iters; ++it) {
    std::fill(y.begin(), y.end(), 0.0);
    A.multiply(x.data(), y.data());
    std::fill(z.begin(), z.end(), 0.0);
    A.multiplyTranspose(y.data(), z.data());
    double nz = 0;
    for (double v : z) nz += v * v;
    nz = std::sqrt(nz);
    if (nz == 0) return 1.0;
    s = std::sqrt(nz);  // ||A^T A x|| -> sigma_max^2
    for (int j = 0; j < n; ++j) x[j] = z[j] / nz;
  }
  return s;
}

Prepared prepare(const Model& lp, const PdhgOptions& opt) {
  Prepared P;
  P.n = lp.numCols();
  P.m = lp.numRows();
  P.offset = lp.objOffset;
  SparseMatrix A = lp.A;
  A.numRows = P.m;
  ScalingResult s;
  if (opt.ruizIterations > 0 && A.nnz() > 0) {
    s = ruizScaling(A, opt.ruizIterations, opt.pockChambolle);
  } else {
    s.rowScale.assign(P.m, 1.0);
    s.colScale.assign(P.n, 1.0);
  }
  applyScaling(A, s);
  P.R = s.rowScale;
  P.Cs = s.colScale;
  if (opt.rescaleObjectiveAndBounds) {
    double cn = 0, bn = 0;
    for (int j = 0; j < P.n; ++j) cn += sqr(lp.colCost[j] * P.Cs[j]);
    for (int i = 0; i < P.m; ++i) {
      if (fin(lp.rowLower[i])) bn += sqr(lp.rowLower[i] * P.R[i]);
      if (fin(lp.rowUpper[i]) && lp.rowUpper[i] != lp.rowLower[i]) bn += sqr(lp.rowUpper[i] * P.R[i]);
    }
    P.gamma = 1.0 + std::sqrt(cn);
    P.beta = 1.0 + std::sqrt(bn);
  }
  P.A = csrOf(A);
  P.AT = csrOfTranspose(A);
  P.rfac.resize(P.m);
  P.cfac.resize(P.n);
  for (int i = 0; i < P.m; ++i) P.rfac[i] = P.beta / P.R[i];
  for (int j = 0; j < P.n; ++j) P.cfac[j] = P.gamma / P.Cs[j];
  P.normA = powerIteration(A, 40) * 1.01;
  return P;
}

struct CaseVectors {
  std::vector<double> c, lo, up, rl, ru;  // original space
  double cNorm = 0, rhsNorm = 0;
};

CaseVectors caseVectors(const Model& lp, const LpCase* cs) {
  CaseVectors v;
  v.c = (cs && !cs->cost.empty()) ? cs->cost : lp.colCost;
  v.lo = (cs && !cs->colLower.empty()) ? cs->colLower : lp.colLower;
  v.up = (cs && !cs->colUpper.empty()) ? cs->colUpper : lp.colUpper;
  v.rl = (cs && !cs->rowLower.empty()) ? cs->rowLower : lp.rowLower;
  v.ru = (cs && !cs->rowUpper.empty()) ? cs->rowUpper : lp.rowUpper;
  for (double c : v.c) v.cNorm += c * c;
  v.cNorm = std::sqrt(v.cNorm);
  for (size_t i = 0; i < v.rl.size(); ++i) {
    if (fin(v.rl[i])) v.rhsNorm += v.rl[i] * v.rl[i];
    if (fin(v.ru[i]) && v.ru[i] != v.rl[i]) v.rhsNorm += v.ru[i] * v.ru[i];
  }
  v.rhsNorm = std::sqrt(v.rhsNorm);
  return v;
}

std::vector<PdhgResult> runPdhg(const Model& lp, const std::vector<const LpCase*>& cases, const PdhgOptions& opt,
                                const Deadline* dl, PdhgStats& st) {
  Timer setupT;
  const int k = static_cast<int>(cases.size());
  Prepared P = prepare(lp, opt);
  const int n = P.n, m = P.m;
  std::unique_ptr<Backend> be;
  auto& gctx = gpu::CudaContext::instance();
  bool useGpu = opt.useGpu && gctx.available();
  if (useGpu) {
    be = std::make_unique<GpuBackend>();
    st.backend = "gpu";
    st.device = gctx.info().name;
  } else {
    be = std::make_unique<CpuBackend>();
    st.backend = "cpu";
    st.device = formatString("%d threads", ThreadPool::instance().threads());
  }
  size_t up0 = gctx.bytesUploaded, down0 = gctx.bytesDownloaded;
  long long launches0 = gctx.launches;
  Timer xferT;
  be->init(P.A, P.AT, P.rfac, P.cfac, k);

  // Case data in scaled space, k-interleaved.
  std::vector<CaseVectors> cv(k);
  std::vector<double> c(static_cast<size_t>(n) * k), lo(c.size()), upv(c.size()), rl(static_cast<size_t>(m) * k),
      ru(rl.size()), x(c.size()), y(rl.size(), 0.0);
  std::vector<double> omega(k, 1.0);
  for (int b = 0; b < k; ++b) {
    cv[b] = caseVectors(lp, cases[b]);
    double cn = 0, bn = 0;
    for (int j = 0; j < n; ++j) {
      size_t t = static_cast<size_t>(j) * k + b;
      c[t] = cv[b].c[j] * P.Cs[j] / P.gamma;
      lo[t] = fin(cv[b].lo[j]) ? cv[b].lo[j] / P.Cs[j] / P.beta : -HUGE_VAL;
      upv[t] = fin(cv[b].up[j]) ? cv[b].up[j] / P.Cs[j] / P.beta : HUGE_VAL;
      x[t] = std::min(std::max(0.0, lo[t]), upv[t]);
      cn += c[t] * c[t];
    }
    for (int i = 0; i < m; ++i) {
      size_t t = static_cast<size_t>(i) * k + b;
      rl[t] = fin(cv[b].rl[i]) ? cv[b].rl[i] * P.R[i] / P.beta : -HUGE_VAL;
      ru[t] = fin(cv[b].ru[i]) ? cv[b].ru[i] * P.R[i] / P.beta : HUGE_VAL;
      if (fin(rl[t])) bn += rl[t] * rl[t];
      if (fin(ru[t]) && ru[t] != rl[t]) bn += ru[t] * ru[t];
    }
    if (cn > 1e-20 && bn > 1e-20) omega[b] = std::sqrt(cn) / std::sqrt(bn);
  }
  be->set(C, c);
  be->set(LO, lo);
  be->set(UP, upv);
  be->set(RL, rl);
  be->set(RU, ru);
  be->set(X, x);
  be->set(X0, x);
  be->set(Y, y);
  be->set(Y0, y);
  be->spmv(X, AX);
  {
    std::vector<double> ax;
    be->get(AX, ax);
    be->set(AX0, ax);
  }
  const double eta = 0.998 / P.normA;
  st.spectralNorm = P.normA;
  std::vector<double> tau(k), sigma(k);
  auto updSteps = [&]() {
    for (int b = 0; b < k; ++b) {
      tau[b] = eta / omega[b];
      sigma[b] = eta * omega[b];
    }
    be->setSteps(tau, sigma);
  };
  updSteps();
  st.setupSeconds = setupT.seconds();
  st.transferSeconds += xferT.seconds();

  // ---- Iterations ----
  std::vector<PdhgResult> res(k);
  std::vector<char> done(k, 0);
  std::vector<double> fpRestart(k, kInf), fpPrev(k, kInf);
  std::vector<long long> sinceRestart(k, 0);
  int remaining = k;
  Timer solveT;
  be->startTiming();
  long long it = 0;
  std::vector<double> r3, c3, dx2, dy2, ax2, ay2;
  int checks = 0;
  auto extract = [&](int b, Status status, double relKkt, double pobj, double dobj) {
    std::vector<double> xs, ys;
    Timer tt;
    be->get(XP, xs);
    be->get(YP, ys);
    st.transferSeconds += tt.seconds();
    PdhgResult& r = res[b];
    r.status = status;
    r.x.resize(n);
    r.rowDual.resize(m);
    for (int j = 0; j < n; ++j) r.x[j] = xs[static_cast<size_t>(j) * k + b] * P.Cs[j] * P.beta;
    for (int i = 0; i < m; ++i) r.rowDual[i] = ys[static_cast<size_t>(i) * k + b] * P.R[i] * P.gamma;
    r.iterations = it;
    r.relativeKkt = relKkt;
    r.objective = pobj + P.offset;
    r.dualObjective = dobj + P.offset;
  };
  for (; it < opt.maxIterations; ++it) {
    be->spmvT(Y, ATY);
    be->primalStep();
    be->spmv(XP, AXP);
    be->dualStep();
    const bool check = ((it + 1) % opt.checkEvery) == 0 || it + 1 == opt.maxIterations;
    std::vector<int> mask;
    if (check) {
      ++checks;
      Timer tt;
      be->spmvT(YP, ATYP);
      be->kkt(r3, c3);
      be->diff2(X, XP, dx2);
      be->diff2(Y, YP, dy2);
      st.transferSeconds += 0;  // reductions are device-side; only O(k) partial sums move
      (void)tt;
      mask.assign(k, 0);
      bool anyRestart = false;
      for (int b = 0; b < k; ++b) {
        if (done[b]) continue;
        double scaleObj = P.beta * P.gamma;
        double pobj = c3[3 * b + 1] * scaleObj;
        double dobj = (r3[3 * b + 1] + c3[3 * b + 2]) * scaleObj;
        double pres = std::sqrt(r3[3 * b]), dres = std::sqrt(c3[3 * b]);
        double relP = pres / (1.0 + cv[b].rhsNorm);
        double relD = dres / (1.0 + cv[b].cNorm);
        double gap = std::fabs(pobj - dobj) / (1.0 + std::fabs(pobj) + std::fabs(dobj));
        double kktErr = std::max(relP, std::max(relD, gap));
        if (k == 1) {
          st.primalResidual = relP;
          st.dualResidual = relD;
          st.gap = gap;
        }
        if (opt.logInterval > 0 && k == 1 && checks % opt.logInterval == 0)
          PLOG_INFO("  pdhg %8lld  pobj % .8e  dobj % .8e  relP %.2e  relD %.2e  gap %.2e  w %.2e  restarts %d  %.2fs",
                    it + 1, pobj + P.offset, dobj + P.offset, relP, relD, gap, omega[b], st.restarts,
                    solveT.seconds());
        if (relP <= opt.tolerance && relD <= opt.tolerance && gap <= opt.tolerance) {
          extract(b, Status::Optimal, kktErr, pobj, dobj);
          done[b] = 1;
          --remaining;
          continue;
        }
        sinceRestart[b] += opt.checkEvery;
        double fp = std::sqrt(omega[b] * dx2[b] + dy2[b] / omega[b]);
        if (fpRestart[b] == kInf) fpRestart[b] = fp;
        bool restart = fp <= 0.2 * fpRestart[b] || (fp <= 0.8 * fpRestart[b] && fp > fpPrev[b]) ||
                       sinceRestart[b] >= static_cast<long long>(0.36 * (it + 1));
        fpPrev[b] = fp;
        if (restart) {
          mask[b] = 1;
          anyRestart = true;
        }
      }
      if (anyRestart) {
        // Primal weight update from the movement since the last restart.
        be->diff2(XP, X0, ax2);
        be->diff2(YP, Y0, ay2);
        for (int b = 0; b < k; ++b) {
          if (!mask[b]) continue;
          double dxn = std::sqrt(ax2[b]), dyn = std::sqrt(ay2[b]);
          if (dxn > 1e-10 && dyn > 1e-10) omega[b] = std::exp(0.5 * std::log(dyn / dxn) + 0.5 * std::log(omega[b]));
          fpRestart[b] = std::sqrt(omega[b] * dx2[b] + dy2[b] / omega[b]);
          fpPrev[b] = kInf;
          sinceRestart[b] = 0;
          st.restarts++;
        }
        updSteps();
      }
      if (remaining == 0) break;
      if (dl && dl->expired()) break;
    }
    be->halpern(opt.reflection);
    if (!mask.empty() && std::find(mask.begin(), mask.end(), 1) != mask.end()) be->restart(mask);
  }
  st.kernelSeconds = be->stopTiming();
  st.solveSeconds = solveT.seconds();
  st.iterations = it + (remaining == 0 ? 1 : 0);
  st.primalWeight = omega.empty() ? 1 : omega[0];
  // Unconverged cases: return the last PDHG output.
  if (remaining > 0) {
    be->spmvT(YP, ATYP);
    be->kkt(r3, c3);
    for (int b = 0; b < k; ++b) {
      if (done[b]) continue;
      double scaleObj = P.beta * P.gamma;
      double pobj = c3[3 * b + 1] * scaleObj, dobj = (r3[3 * b + 1] + c3[3 * b + 2]) * scaleObj;
      double relP = std::sqrt(r3[3 * b]) / (1 + cv[b].rhsNorm), relD = std::sqrt(c3[3 * b]) / (1 + cv[b].cNorm);
      double gap = std::fabs(pobj - dobj) / (1 + std::fabs(pobj) + std::fabs(dobj));
      extract(b, (dl && dl->expired()) ? Status::TimeLimit : Status::IterationLimit,
              std::max(relP, std::max(relD, gap)), pobj, dobj);
    }
  }
  st.bytesUploaded = static_cast<long long>(gctx.bytesUploaded - up0);
  st.bytesDownloaded = static_cast<long long>(gctx.bytesDownloaded - down0);
  st.kernelLaunches = gctx.launches - launches0;
  if (!useGpu) st.kernelSeconds = st.solveSeconds;
  for (auto& r : res) r.stats = st;
  return res;
}

}  // namespace

Json PdhgStats::toJson() const {
  Json j = Json::object();
  j["backend"] = backend;
  j["device"] = device;
  j["iterations"] = static_cast<double>(iterations);
  j["restarts"] = restarts;
  j["relative_primal_residual"] = primalResidual;
  j["relative_dual_residual"] = dualResidual;
  j["relative_gap"] = gap;
  j["primal_weight"] = primalWeight;
  j["spectral_norm_estimate"] = spectralNorm;
  j["setup_seconds"] = setupSeconds;
  j["solve_seconds"] = solveSeconds;
  j["kernel_seconds"] = kernelSeconds;
  j["transfer_seconds"] = transferSeconds;
  j["bytes_uploaded"] = static_cast<double>(bytesUploaded);
  j["bytes_downloaded"] = static_cast<double>(bytesDownloaded);
  j["kernel_launches"] = static_cast<double>(kernelLaunches);
  return j;
}

PdhgResult solvePdhg(const Model& lp, const PdhgOptions& opts, const Deadline* deadline) {
  PdhgStats st;
  std::vector<const LpCase*> cases{nullptr};
  std::vector<PdhgResult> r = runPdhg(lp, cases, opts, deadline, st);
  return r[0];
}

std::vector<PdhgResult> solvePdhgBatch(const Model& base, const std::vector<LpCase>& cases, const PdhgOptions& opts,
                                       const Deadline* deadline, PdhgStats* batchStats) {
  PdhgStats st;
  std::vector<const LpCase*> ptrs;
  for (auto& c : cases) ptrs.push_back(&c);
  if (ptrs.empty()) return {};
  std::vector<PdhgResult> r = runPdhg(base, ptrs, opts, deadline, st);
  if (batchStats) *batchStats = st;
  return r;
}

bool gpuAvailable() { return gpu::CudaContext::instance().available(); }

Json gpuInfoJson() {
  const auto& info = gpu::CudaContext::instance().info();
  Json j = Json::object();
  j["available"] = info.available;
  j["name"] = info.name;
  j["reason"] = info.reason;
  j["compute_capability"] = formatString("%d.%d", info.computeMajor, info.computeMinor);
  j["multiprocessors"] = info.multiprocessors;
  j["memory_bytes"] = static_cast<double>(info.totalMemory);
  j["peak_bandwidth_GBs"] = info.peakBandwidthGBs;
  j["driver_version"] = info.driverVersion;
  j["nvrtc"] = info.nvrtcPath;
  j["ptx_from_cache"] = info.ptxFromCache;
  return j;
}

}  // namespace pramana
