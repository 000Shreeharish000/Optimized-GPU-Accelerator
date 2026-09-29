// CUDA C source for PRAMANA's GPU kernels, compiled at runtime by NVRTC
// (see gpu/cuda_driver.cpp). All kernels are batched over k related LPs that
// share A: vectors are stored k-interleaved (element (i, b) at i*k + b), so a
// batch of k SpMVs becomes one coalesced SpMM that reads A once.
#include <string>

#include "gpu/cuda_driver.h"

namespace pramana {
namespace gpu {

const std::string& kernelSource() {
  static const std::string src = R"CUDA(
extern "C" {

#define TID ((long long)blockIdx.x * blockDim.x + threadIdx.x)

// y(i,b) = sum_e val[e] * x(col[e], b)        thread per (row, batch)
__global__ void spmm_thread(int rows, int k, const int* __restrict__ ptr, const int* __restrict__ col,
                            const double* __restrict__ val, const double* __restrict__ x, double* __restrict__ y) {
  long long t = TID;
  if (t >= (long long)rows * k) return;
  int i = (int)(t / k), b = (int)(t % k);
  double s = 0.0;
  for (int e = ptr[i]; e < ptr[i + 1]; ++e) s += val[e] * x[(long long)col[e] * k + b];
  y[t] = s;
}

// k == 1, long rows: one warp per row, shuffle reduction.
__global__ void spmv_warp(int rows, const int* __restrict__ ptr, const int* __restrict__ col,
                          const double* __restrict__ val, const double* __restrict__ x, double* __restrict__ y) {
  long long gt = TID;
  int warp = (int)(gt >> 5);
  int lane = threadIdx.x & 31;
  if (warp >= rows) return;
  double s = 0.0;
  for (int e = ptr[warp] + lane; e < ptr[warp + 1]; e += 32) s += val[e] * x[col[e]];
  for (int off = 16; off > 0; off >>= 1) s += __shfl_down_sync(0xffffffffu, s, off);
  if (lane == 0) y[warp] = s;
}

// x+ = proj_[lo,up]( x - tau_b (c - A^T y) )
__global__ void pdhg_primal(int n, int k, const double* __restrict__ x, const double* __restrict__ aty,
                            const double* __restrict__ c, const double* __restrict__ lo, const double* __restrict__ up,
                            const double* __restrict__ tau, double* __restrict__ xp) {
  long long t = TID;
  if (t >= (long long)n * k) return;
  int b = (int)(t % k);
  double v = x[t] - tau[b] * (c[t] - aty[t]);
  v = fmax(v, lo[t]);
  v = fmin(v, up[t]);
  xp[t] = v;
}

// y+ = v + sigma_b * proj_[L,U](-v / sigma_b),  v = y - sigma_b (2 A x+ - A x)
__global__ void pdhg_dual(int m, int k, const double* __restrict__ y, const double* __restrict__ axp,
                          const double* __restrict__ ax, const double* __restrict__ rl, const double* __restrict__ ru,
                          const double* __restrict__ sigma, double* __restrict__ yp) {
  long long t = TID;
  if (t >= (long long)m * k) return;
  int b = (int)(t % k);
  double s = sigma[b];
  double v = y[t] - s * (2.0 * axp[t] - ax[t]);
  double w = -v / s;
  w = fmax(w, rl[t]);
  w = fmin(w, ru[t]);
  yp[t] = v + s * w;
}

// Reflected Halpern step: z = lam((1+rho) z+ - rho z) + (1-lam) z0,  lam = (k+1)/(k+2)
__global__ void halpern(int len, int k, const int* __restrict__ iter, double rho, double* __restrict__ z,
                        const double* __restrict__ zp, const double* __restrict__ z0) {
  long long t = TID;
  if (t >= (long long)len * k) return;
  int b = (int)(t % k);
  double kk = (double)iter[b];
  double lam = (kk + 1.0) / (kk + 2.0);
  z[t] = lam * ((1.0 + rho) * zp[t] - rho * z[t]) + (1.0 - lam) * z0[t];
}

__global__ void advance_iter(int k, int* iter) {
  int b = blockIdx.x * blockDim.x + threadIdx.x;
  if (b < k) iter[b] += 1;
}

// Restart (masked per LP): dst1 = dst2 = src
__global__ void restart_copy(int len, int k, const int* __restrict__ mask, const double* __restrict__ src,
                             double* __restrict__ dst1, double* __restrict__ dst2) {
  long long t = TID;
  if (t >= (long long)len * k) return;
  int b = (int)(t % k);
  if (mask[b]) {
    dst1[t] = src[t];
    dst2[t] = src[t];
  }
}

__global__ void reset_iter(int k, const int* mask, int* iter) {
  int b = blockIdx.x * blockDim.x + threadIdx.x;
  if (b < k && mask[b]) iter[b] = 0;
}

// Block tree reduction helper (blockDim.x must be a power of two <= 256).
__device__ void block_sum(double* sh, int nq) {
  for (int s = blockDim.x / 2; s > 0; s >>= 1) {
    if (threadIdx.x < s)
      for (int q = 0; q < nq; ++q) sh[q * blockDim.x + threadIdx.x] += sh[q * blockDim.x + threadIdx.x + s];
    __syncthreads();
  }
}

// partial[(b * gridDim.x + blk)] = sum_i (a(i,b) - c(i,b))^2      grid = (blocks, k)
__global__ void reduce_diff2(int len, int k, const double* __restrict__ a, const double* __restrict__ c,
                             double* __restrict__ partial) {
  __shared__ double sh[256];
  int b = blockIdx.y;
  double s = 0.0;
  for (long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x; i < len; i += (long long)gridDim.x * blockDim.x) {
    double d = a[i * k + b] - c[i * k + b];
    s += d * d;
  }
  sh[threadIdx.x] = s;
  __syncthreads();
  block_sum(sh, 1);
  if (threadIdx.x == 0) partial[(long long)b * gridDim.x + blockIdx.x] = sh[0];
}

// Row part of the KKT error at (x+, y+): 3 partial sums per block:
//   [0] sum (unscaled primal residual)^2, [1] dual objective row part, [2] sum y^2
__global__ void kkt_rows(int m, int k, const double* __restrict__ axp, const double* __restrict__ rl,
                         const double* __restrict__ ru, const double* __restrict__ yp,
                         const double* __restrict__ rfac, double* __restrict__ partial) {
  __shared__ double sh[3 * 256];
  int b = blockIdx.y;
  double r2 = 0.0, dob = 0.0, y2 = 0.0;
  for (long long i = (long long)blockIdx.x * blockDim.x + threadIdx.x; i < m; i += (long long)gridDim.x * blockDim.x) {
    long long t = i * k + b;
    double a = axp[t];
    double res = fmax(rl[t] - a, 0.0) + fmax(a - ru[t], 0.0);
    res *= rfac[i];
    r2 += res * res;
    double y = yp[t];
    if (y > 0.0 && rl[t] > -1e300) dob += y * rl[t];
    else if (y < 0.0 && ru[t] < 1e300) dob += y * ru[t];
    y2 += y * y;
  }
  sh[threadIdx.x] = r2;
  sh[blockDim.x + threadIdx.x] = dob;
  sh[2 * blockDim.x + threadIdx.x] = y2;
  __syncthreads();
  block_sum(sh, 3);
  if (threadIdx.x == 0) {
    long long o = ((long long)b * gridDim.x + blockIdx.x) * 3;
    partial[o] = sh[0];
    partial[o + 1] = sh[blockDim.x];
    partial[o + 2] = sh[2 * blockDim.x];
  }
}

// Column part: [0] sum (unscaled dual residual)^2, [1] primal objective, [2] dual objective column part
__global__ void kkt_cols(int n, int k, const double* __restrict__ xp, const double* __restrict__ c,
                         const double* __restrict__ atyp, const double* __restrict__ lo, const double* __restrict__ up,
                         const double* __restrict__ cfac, double* __restrict__ partial) {
  __shared__ double sh[3 * 256];
  int b = blockIdx.y;
  double r2 = 0.0, pob = 0.0, dob = 0.0;
  for (long long j = (long long)blockIdx.x * blockDim.x + threadIdx.x; j < n; j += (long long)gridDim.x * blockDim.x) {
    long long t = j * k + b;
    double d = c[t] - atyp[t];
    bool hl = lo[t] > -1e300, hu = up[t] < 1e300;
    double res = 0.0;
    if (!hl && !hu) res = fabs(d);
    else if (!hl) res = fmax(d, 0.0);
    else if (!hu) res = fmax(-d, 0.0);
    res *= cfac[j];
    r2 += res * res;
    pob += c[t] * xp[t];
    if (d > 0.0 && hl) dob += d * lo[t];
    else if (d < 0.0 && hu) dob += d * up[t];
  }
  sh[threadIdx.x] = r2;
  sh[blockDim.x + threadIdx.x] = pob;
  sh[2 * blockDim.x + threadIdx.x] = dob;
  __syncthreads();
  block_sum(sh, 3);
  if (threadIdx.x == 0) {
    long long o = ((long long)b * gridDim.x + blockIdx.x) * 3;
    partial[o] = sh[0];
    partial[o + 1] = sh[blockDim.x];
    partial[o + 2] = sh[2 * blockDim.x];
  }
}

}  // extern "C"
)CUDA";
  return src;
}

}  // namespace gpu
}  // namespace pramana
