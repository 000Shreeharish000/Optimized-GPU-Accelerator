// Runtime loader for the NVIDIA CUDA *driver* API (nvcuda.dll / libcuda.so)
// and NVRTC (runtime CUDA C++ compiler). Nothing CUDA-related is linked at
// build time: the solver binary runs unchanged on machines without a GPU and
// falls back to the CPU backend. Kernels are PRAMANA's own CUDA C source,
// compiled once by NVRTC and cached as PTX (the driver JITs PTX to SASS).
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace pramana {
namespace gpu {

using CUdeviceptr = unsigned long long;
using CUresult = int;
using CUdevice = int;
struct CUctx_st;
struct CUmod_st;
struct CUfunc_st;
struct CUstream_st;
struct CUevent_st;
using CUcontext = CUctx_st*;
using CUmodule = CUmod_st*;
using CUfunction = CUfunc_st*;
using CUstream = CUstream_st*;
using CUevent = CUevent_st*;

struct DeviceInfo {
  bool available = false;
  std::string name;
  std::string reason;          // why unavailable
  int computeMajor = 0, computeMinor = 0;
  int multiprocessors = 0;
  size_t totalMemory = 0;
  double peakBandwidthGBs = 0;  // from memory clock * bus width
  int driverVersion = 0;
  std::string nvrtcPath;
  bool ptxFromCache = false;
};

class CudaContext {
 public:
  static CudaContext& instance();
  bool available() const { return info_.available; }
  const DeviceInfo& info() const { return info_; }

  CUdeviceptr alloc(size_t bytes);
  void free(CUdeviceptr p);
  void upload(CUdeviceptr dst, const void* src, size_t bytes);
  void download(void* dst, CUdeviceptr src, size_t bytes);
  void copyDevice(CUdeviceptr dst, CUdeviceptr src, size_t bytes);
  void memsetZero(CUdeviceptr p, size_t bytes);
  CUfunction function(const std::string& name);
  void launch(CUfunction f, unsigned gridX, unsigned gridY, unsigned blockX, void** args, unsigned sharedBytes = 0);
  void synchronize();
  // Event timing (milliseconds between two recorded events).
  CUevent createEvent();
  void record(CUevent e);
  float elapsedMs(CUevent a, CUevent b);
  void destroyEvent(CUevent e);

  size_t bytesUploaded = 0, bytesDownloaded = 0;
  long long launches = 0;

 private:
  CudaContext();
  bool loadDriver();
  bool loadNvrtcAndBuild(const std::string& source);
  void check(CUresult r, const char* what);

  DeviceInfo info_;
  void* driverLib_ = nullptr;
  void* nvrtcLib_ = nullptr;
  CUcontext ctx_ = nullptr;
  CUmodule module_ = nullptr;
  CUdevice device_ = 0;
  struct Fn;
  Fn* fn_ = nullptr;
};

// CUDA C source of all PRAMANA kernels (pdhg/kernels.cpp).
const std::string& kernelSource();

}  // namespace gpu
}  // namespace pramana
