#include "gpu/cuda_driver.h"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <sstream>

#include "util/common.h"
#include "util/log.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace pramana {
namespace gpu {

namespace {

void* openLib(const std::string& name) {
#ifdef _WIN32
  const size_t slash = name.find_last_of("\\/");
  if (slash != std::string::npos) {
    // NVRTC loads its companion nvrtc-builtins DLL by name at compile time:
    // make its own directory part of the DLL search path.
    SetDllDirectoryA(name.substr(0, slash).c_str());
  }
  HMODULE h = LoadLibraryExA(name.c_str(), nullptr, slash != std::string::npos ? LOAD_WITH_ALTERED_SEARCH_PATH : 0);
  return reinterpret_cast<void*>(h);
#else
  return dlopen(name.c_str(), RTLD_NOW | RTLD_LOCAL);
#endif
}

void* sym(void* lib, const char* name) {
#ifdef _WIN32
  return reinterpret_cast<void*>(GetProcAddress(reinterpret_cast<HMODULE>(lib), name));
#else
  return dlsym(lib, name);
#endif
}

bool fileExists(const std::string& p) {
  std::ifstream f(p, std::ios::binary);
  return f.good();
}

std::string envOr(const char* name, const std::string& dflt) {
  const char* v = std::getenv(name);
  return v ? std::string(v) : dflt;
}

std::vector<std::string> nvrtcCandidates() {
  std::vector<std::string> c;
  std::string env = envOr("PRAMANA_NVRTC", "");
  if (!env.empty()) c.push_back(env);
#ifdef _WIN32
  const char* names[] = {"nvrtc64_130_0.dll", "nvrtc64_120_0.dll", "nvrtc64_112_0.dll"};
  std::string cudaPath = envOr("CUDA_PATH", "");
  if (!cudaPath.empty())
    for (auto* n : names) {
      c.push_back(cudaPath + "\\bin\\" + n);
      c.push_back(cudaPath + "\\bin\\x64\\" + n);
    }
  // pip wheel "nvidia-cuda-nvrtc-cu12" under Python installations.
  std::vector<std::string> patterns;  // directories matching Python3*
  std::string lad = envOr("LOCALAPPDATA", "");
  std::string ad = envOr("APPDATA", "");
  if (!lad.empty()) patterns.push_back(lad + "\\Programs\\Python\\");
  if (!ad.empty()) patterns.push_back(ad + "\\Python\\");
  patterns.push_back("C:\\");
  for (auto& dir : patterns) {
    WIN32_FIND_DATAA fd;
    HANDLE h = FindFirstFileA((dir + "Python3*").c_str(), &fd);
    if (h == INVALID_HANDLE_VALUE) continue;
    do {
      std::string base = dir + fd.cFileName;
      for (auto* n : names) {
        c.push_back(base + "\\Lib\\site-packages\\nvidia\\cuda_nvrtc\\bin\\" + n);
        c.push_back(base + "\\site-packages\\nvidia\\cuda_nvrtc\\bin\\" + n);
      }
    } while (FindNextFileA(h, &fd));
    FindClose(h);
  }
  for (auto* n : names) c.push_back(n);  // finally: whatever is on PATH
#else
  const char* names[] = {"libnvrtc.so", "libnvrtc.so.13", "libnvrtc.so.12"};
  for (auto* n : names) c.push_back(n);
  c.push_back("/usr/local/cuda/lib64/libnvrtc.so");
#endif
  return c;
}

std::string cacheDir() {
#ifdef _WIN32
  std::string base = envOr("LOCALAPPDATA", envOr("TEMP", "."));
  std::string dir = base + "\\pramana";
  CreateDirectoryA(dir.c_str(), nullptr);
  return dir + "\\";
#else
  std::string base = envOr("XDG_CACHE_HOME", envOr("HOME", ".") + "/.cache");
  std::string dir = base + "/pramana";
  std::system(("mkdir -p '" + dir + "'").c_str());
  return dir + "/";
#endif
}

}  // namespace

// Function pointer table (driver API ABI; CUDAAPI == __stdcall on 32-bit only).
struct CudaContext::Fn {
  CUresult (*cuInit)(unsigned);
  CUresult (*cuDriverGetVersion)(int*);
  CUresult (*cuDeviceGetCount)(int*);
  CUresult (*cuDeviceGet)(CUdevice*, int);
  CUresult (*cuDeviceGetName)(char*, int, CUdevice);
  CUresult (*cuDeviceGetAttribute)(int*, int, CUdevice);
  CUresult (*cuDeviceTotalMem)(size_t*, CUdevice);
  CUresult (*cuDevicePrimaryCtxRetain)(CUcontext*, CUdevice);
  CUresult (*cuCtxSetCurrent)(CUcontext);
  CUresult (*cuCtxSynchronize)();
  CUresult (*cuModuleLoadDataEx)(CUmodule*, const void*, unsigned, int*, void**);
  CUresult (*cuModuleGetFunction)(CUfunction*, CUmodule, const char*);
  CUresult (*cuMemAlloc)(CUdeviceptr*, size_t);
  CUresult (*cuMemFree)(CUdeviceptr);
  CUresult (*cuMemcpyHtoD)(CUdeviceptr, const void*, size_t);
  CUresult (*cuMemcpyDtoH)(void*, CUdeviceptr, size_t);
  CUresult (*cuMemcpyDtoD)(CUdeviceptr, CUdeviceptr, size_t);
  CUresult (*cuMemsetD8)(CUdeviceptr, unsigned char, size_t);
  CUresult (*cuLaunchKernel)(CUfunction, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned,
                             CUstream, void**, void**);
  CUresult (*cuGetErrorString)(CUresult, const char**);
  CUresult (*cuEventCreate)(CUevent*, unsigned);
  CUresult (*cuEventRecord)(CUevent, CUstream);
  CUresult (*cuEventSynchronize)(CUevent);
  CUresult (*cuEventElapsedTime)(float*, CUevent, CUevent);
  CUresult (*cuEventDestroy)(CUevent);
};

CudaContext& CudaContext::instance() {
  static CudaContext ctx;
  return ctx;
}

CudaContext::CudaContext() {
  if (std::getenv("PRAMANA_DISABLE_GPU")) {
    info_.reason = "disabled by PRAMANA_DISABLE_GPU";
    return;
  }
  try {
    if (!loadDriver()) return;
    if (!loadNvrtcAndBuild(kernelSource())) return;
    info_.available = true;
  } catch (std::exception& e) {
    info_.available = false;
    info_.reason = e.what();
  }
  if (!info_.available) PLOG_DETAIL("gpu: unavailable (%s)", info_.reason.c_str());
}

void CudaContext::check(CUresult r, const char* what) {
  if (r == 0) return;
  const char* msg = "unknown";
  if (fn_ && fn_->cuGetErrorString) fn_->cuGetErrorString(r, &msg);
  throw PramanaError(std::string("CUDA driver error in ") + what + ": " + msg + " (" + std::to_string(r) + ")");
}

bool CudaContext::loadDriver() {
#ifdef _WIN32
  driverLib_ = openLib("nvcuda.dll");
#else
  driverLib_ = openLib("libcuda.so.1");
#endif
  if (!driverLib_) {
    info_.reason = "CUDA driver library not found";
    return false;
  }
  fn_ = new Fn();
  auto load = [&](auto& f, const char* n) {
    f = reinterpret_cast<std::remove_reference_t<decltype(f)>>(sym(driverLib_, n));
    if (!f) throw PramanaError(std::string("missing driver symbol ") + n);
  };
  load(fn_->cuInit, "cuInit");
  load(fn_->cuDriverGetVersion, "cuDriverGetVersion");
  load(fn_->cuDeviceGetCount, "cuDeviceGetCount");
  load(fn_->cuDeviceGet, "cuDeviceGet");
  load(fn_->cuDeviceGetName, "cuDeviceGetName");
  load(fn_->cuDeviceGetAttribute, "cuDeviceGetAttribute");
  load(fn_->cuDeviceTotalMem, "cuDeviceTotalMem_v2");
  load(fn_->cuDevicePrimaryCtxRetain, "cuDevicePrimaryCtxRetain");
  load(fn_->cuCtxSetCurrent, "cuCtxSetCurrent");
  load(fn_->cuCtxSynchronize, "cuCtxSynchronize");
  load(fn_->cuModuleLoadDataEx, "cuModuleLoadDataEx");
  load(fn_->cuModuleGetFunction, "cuModuleGetFunction");
  load(fn_->cuMemAlloc, "cuMemAlloc_v2");
  load(fn_->cuMemFree, "cuMemFree_v2");
  load(fn_->cuMemcpyHtoD, "cuMemcpyHtoD_v2");
  load(fn_->cuMemcpyDtoH, "cuMemcpyDtoH_v2");
  load(fn_->cuMemcpyDtoD, "cuMemcpyDtoD_v2");
  load(fn_->cuMemsetD8, "cuMemsetD8_v2");
  load(fn_->cuLaunchKernel, "cuLaunchKernel");
  load(fn_->cuGetErrorString, "cuGetErrorString");
  load(fn_->cuEventCreate, "cuEventCreate");
  load(fn_->cuEventRecord, "cuEventRecord");
  load(fn_->cuEventSynchronize, "cuEventSynchronize");
  load(fn_->cuEventElapsedTime, "cuEventElapsedTime");
  load(fn_->cuEventDestroy, "cuEventDestroy_v2");

  if (fn_->cuInit(0) != 0) {
    info_.reason = "cuInit failed (no NVIDIA device/driver)";
    return false;
  }
  fn_->cuDriverGetVersion(&info_.driverVersion);
  int count = 0;
  fn_->cuDeviceGetCount(&count);
  if (count <= 0) {
    info_.reason = "no CUDA device";
    return false;
  }
  int ordinal = std::atoi(envOr("PRAMANA_GPU_DEVICE", "0").c_str());
  check(fn_->cuDeviceGet(&device_, ordinal), "cuDeviceGet");
  char name[256] = {};
  fn_->cuDeviceGetName(name, 255, device_);
  info_.name = name;
  fn_->cuDeviceGetAttribute(&info_.computeMajor, 75, device_);  // CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR
  fn_->cuDeviceGetAttribute(&info_.computeMinor, 76, device_);
  fn_->cuDeviceGetAttribute(&info_.multiprocessors, 16, device_);
  int memClockKHz = 0, busWidth = 0;
  fn_->cuDeviceGetAttribute(&memClockKHz, 36, device_);  // MEMORY_CLOCK_RATE (kHz)
  fn_->cuDeviceGetAttribute(&busWidth, 37, device_);     // GLOBAL_MEMORY_BUS_WIDTH (bits)
  info_.peakBandwidthGBs = 2.0 * memClockKHz * 1e3 * (busWidth / 8.0) / 1e9;
  fn_->cuDeviceTotalMem(&info_.totalMemory, device_);
  check(fn_->cuDevicePrimaryCtxRetain(&ctx_, device_), "cuDevicePrimaryCtxRetain");
  check(fn_->cuCtxSetCurrent(ctx_), "cuCtxSetCurrent");
  return true;
}

bool CudaContext::loadNvrtcAndBuild(const std::string& source) {
  const int major = info_.computeMajor, minor = info_.computeMinor;
  int arch = major * 10 + minor;
  if (arch > 90) arch = 90;  // PTX is forward compatible through the driver JIT
  std::hash<std::string> h;
  std::ostringstream key;
  key << std::hex << h(source) << "_sm" << std::dec << arch;
  const std::string cachePath = cacheDir() + "kernels_" + key.str() + ".ptx";
  std::string ptx;
  if (fileExists(cachePath) && !std::getenv("PRAMANA_NO_PTX_CACHE")) {
    std::ifstream f(cachePath, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    ptx = ss.str();
    info_.ptxFromCache = true;
  } else {
    // Try every NVRTC we can find; for each, the newest PTX target it accepts
    // (PTX is forward compatible: the driver JIT-compiles it for this GPU).
    using nvrtcProgram = void*;
    std::string lastErr = "NVRTC not found (set PRAMANA_NVRTC or `pip install nvidia-cuda-nvrtc-cu12`) and no cached PTX";
    bool built = false;
    for (auto& cand : nvrtcCandidates()) {
      void* lib = openLib(cand);
      if (!lib) continue;
      int (*nvrtcCreateProgram)(nvrtcProgram*, const char*, const char*, int, const char**, const char**) = nullptr;
      int (*nvrtcCompileProgram)(nvrtcProgram, int, const char**) = nullptr;
      int (*nvrtcGetPTXSize)(nvrtcProgram, size_t*) = nullptr;
      int (*nvrtcGetPTX)(nvrtcProgram, char*) = nullptr;
      int (*nvrtcGetProgramLogSize)(nvrtcProgram, size_t*) = nullptr;
      int (*nvrtcGetProgramLog)(nvrtcProgram, char*) = nullptr;
      int (*nvrtcDestroyProgram)(nvrtcProgram*) = nullptr;
      bool symsOk = true;
      auto load = [&](auto& f, const char* n) {
        f = reinterpret_cast<std::remove_reference_t<decltype(f)>>(sym(lib, n));
        if (!f) symsOk = false;
      };
      load(nvrtcCreateProgram, "nvrtcCreateProgram");
      load(nvrtcCompileProgram, "nvrtcCompileProgram");
      load(nvrtcGetPTXSize, "nvrtcGetPTXSize");
      load(nvrtcGetPTX, "nvrtcGetPTX");
      load(nvrtcGetProgramLogSize, "nvrtcGetProgramLogSize");
      load(nvrtcGetProgramLog, "nvrtcGetProgramLog");
      load(nvrtcDestroyProgram, "nvrtcDestroyProgram");
      if (!symsOk) continue;
      std::vector<int> archs{arch};
      for (int a : {90, 89, 87, 86, 80, 75, 70, 60, 52})
        if (a < arch) archs.push_back(a);
      for (int a : archs) {
        nvrtcProgram prog = nullptr;
        if (nvrtcCreateProgram(&prog, source.c_str(), "pramana_kernels.cu", 0, nullptr, nullptr) != 0) break;
        std::string archOpt = "--gpu-architecture=compute_" + std::to_string(a);
        // No fast math, no FMA contraction: IEEE semantics for residuals that feed certificates.
        const char* opts[] = {archOpt.c_str(), "--std=c++14", "--fmad=false"};
        int rc = nvrtcCompileProgram(prog, 3, opts);
        if (rc != 0) {
          size_t ls = 0;
          nvrtcGetProgramLogSize(prog, &ls);
          std::string log(ls, '\0');
          nvrtcGetProgramLog(prog, log.data());
          nvrtcDestroyProgram(&prog);
          lastErr = "NVRTC (" + cand + ") compile failed for compute_" + std::to_string(a) + ": " + log;
          if (log.find("gpu-architecture") != std::string::npos) continue;  // try an older target
          break;
        }
        size_t ps = 0;
        nvrtcGetPTXSize(prog, &ps);
        ptx.assign(ps, '\0');
        nvrtcGetPTX(prog, ptx.data());
        nvrtcDestroyProgram(&prog);
        info_.nvrtcPath = cand + " (compute_" + std::to_string(a) + ")";
        nvrtcLib_ = lib;
        built = true;
        break;
      }
      if (built) break;
    }
    if (!built) {
      info_.reason = lastErr;
      return false;
    }
    std::ofstream f(cachePath, std::ios::binary);
    f << ptx;
  }
  CUresult r = fn_->cuModuleLoadDataEx(&module_, ptx.c_str(), 0, nullptr, nullptr);
  if (r != 0) {
    info_.reason = "cuModuleLoadDataEx failed (" + std::to_string(r) + ")";
    return false;
  }
  return true;
}

CUdeviceptr CudaContext::alloc(size_t bytes) {
  CUdeviceptr p = 0;
  check(fn_->cuMemAlloc(&p, bytes == 0 ? 8 : bytes), "cuMemAlloc");
  return p;
}
void CudaContext::free(CUdeviceptr p) {
  if (p) fn_->cuMemFree(p);
}
void CudaContext::upload(CUdeviceptr dst, const void* src, size_t bytes) {
  if (!bytes) return;
  check(fn_->cuMemcpyHtoD(dst, src, bytes), "cuMemcpyHtoD");
  bytesUploaded += bytes;
}
void CudaContext::download(void* dst, CUdeviceptr src, size_t bytes) {
  if (!bytes) return;
  check(fn_->cuMemcpyDtoH(dst, src, bytes), "cuMemcpyDtoH");
  bytesDownloaded += bytes;
}
void CudaContext::copyDevice(CUdeviceptr dst, CUdeviceptr src, size_t bytes) {
  if (bytes) check(fn_->cuMemcpyDtoD(dst, src, bytes), "cuMemcpyDtoD");
}
void CudaContext::memsetZero(CUdeviceptr p, size_t bytes) {
  if (bytes) check(fn_->cuMemsetD8(p, 0, bytes), "cuMemsetD8");
}
CUfunction CudaContext::function(const std::string& name) {
  CUfunction f = nullptr;
  check(fn_->cuModuleGetFunction(&f, module_, name.c_str()), name.c_str());
  return f;
}
void CudaContext::launch(CUfunction f, unsigned gridX, unsigned gridY, unsigned blockX, void** args, unsigned shared) {
  if (gridX == 0 || gridY == 0) return;
  check(fn_->cuLaunchKernel(f, gridX, gridY, 1, blockX, 1, 1, shared, nullptr, args, nullptr), "cuLaunchKernel");
  ++launches;
}
void CudaContext::synchronize() { check(fn_->cuCtxSynchronize(), "cuCtxSynchronize"); }
CUevent CudaContext::createEvent() {
  CUevent e = nullptr;
  check(fn_->cuEventCreate(&e, 0), "cuEventCreate");
  return e;
}
void CudaContext::record(CUevent e) { check(fn_->cuEventRecord(e, nullptr), "cuEventRecord"); }
float CudaContext::elapsedMs(CUevent a, CUevent b) {
  check(fn_->cuEventSynchronize(b), "cuEventSynchronize");
  float ms = 0;
  check(fn_->cuEventElapsedTime(&ms, a, b), "cuEventElapsedTime");
  return ms;
}
void CudaContext::destroyEvent(CUevent e) {
  if (e) fn_->cuEventDestroy(e);
}

}  // namespace gpu
}  // namespace pramana
