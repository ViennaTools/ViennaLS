// C++ (g++) interface to the GPU linear solves of the oxidation solvers.
//
// This header is safe to include from any .cpp file compiled by g++.  It
// declares free functions that call into ViennaLS_GPU (gpu/lsOxidationGpu.cu,
// compiled by nvcc) through the C ABI of lsOxidationGpuAbi.hpp.
//
// ── Why dlopen instead of linking ──────────────────────────────────────────
// ViennaLS must not carry a hard (DT_NEEDED) dependency on libcudart /
// libcusparse.  If it did, the dynamic loader would fail to load ViennaLS at
// all on a machine without a CUDA runtime — before a single line of our code
// runs, so no amount of runtime checking could recover.  That is exactly what
// broke the 5.8.3 PyPI wheel.
//
// Instead, the CUDA-linked code lives in a separate ViennaLS_GPU shared
// library which is opened lazily the first time the GPU path is requested.
// If it (or the CUDA runtime it needs) cannot be loaded, the solves report
// the GPU as unavailable and callers fall back to the CPU solver, with the
// reason available from gpuGetLastErrorMessage().

#pragma once

#ifdef VIENNALS_GPU_BICGSTAB

#include "lsOxidationGpuAbi.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>

#if defined(_WIN32)
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace viennals {
namespace gpu {
namespace detail {

#if defined(_WIN32)
using LibraryHandle = HMODULE;
inline LibraryHandle openLibrary(const char *path) {
  return LoadLibraryA(path);
}
inline void *findSymbol(LibraryHandle lib, const char *name) {
  return reinterpret_cast<void *>(GetProcAddress(lib, name));
}
inline std::string lastLoadError() {
  return "LoadLibrary failed (error " + std::to_string(GetLastError()) + ")";
}
inline const char *gpuLibraryName() { return "ViennaLS_GPU.dll"; }
inline char pathSeparator() { return '\\'; }
#else
using LibraryHandle = void *;
inline LibraryHandle openLibrary(const char *path) {
  return dlopen(path, RTLD_NOW | RTLD_LOCAL);
}
inline void *findSymbol(LibraryHandle lib, const char *name) {
  return dlsym(lib, name);
}
inline std::string lastLoadError() {
  const char *err = dlerror();
  return err ? std::string(err) : std::string("unknown dynamic loader error");
}
inline const char *gpuLibraryName() { return "libViennaLS_GPU.so"; }
inline char pathSeparator() { return '/'; }
#endif

/// Directory containing the binary this code is compiled into (e.g. the
/// Python extension module).  The GPU library is shipped alongside it.
inline std::string currentModuleDirectory() {
  std::string path;
#if defined(_WIN32)
  HMODULE module = nullptr;
  if (GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                             GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                         reinterpret_cast<LPCSTR>(&currentModuleDirectory),
                         &module) &&
      module) {
    char buffer[MAX_PATH] = {};
    if (GetModuleFileNameA(module, buffer, MAX_PATH))
      path = buffer;
  }
#else
  Dl_info info{};
  if (dladdr(reinterpret_cast<const void *>(&currentModuleDirectory), &info) &&
      info.dli_fname)
    path = info.dli_fname;
#endif
  const auto pos = path.find_last_of(pathSeparator());
  return pos == std::string::npos ? std::string() : path.substr(0, pos);
}

/// Lazily-loaded handle to ViennaLS_GPU plus the resolved entry points.
/// Constructed exactly once, on first use, by the function-local static in
/// runtime() below.
struct GpuRuntime {
  bool available = false;
  std::string status;

  int (*abiVersion)(void) = nullptr;
  const char *(*lastErrorMessage)(void) = nullptr;
  int (*diffusionSolve)(int64_t, const int64_t *, const int64_t *,
                        const double *, const double *, double *, unsigned,
                        double, unsigned *) = nullptr;
  void *(*macSetup)(int64_t, const int64_t *, const int64_t *, const double *,
                    int64_t, int, const double *) = nullptr;
  int (*macSolve)(void *, int64_t, const int64_t *, const int64_t *,
                  const double *, const double *, double *, unsigned, double,
                  unsigned *, double *) = nullptr;
  void (*macFree)(void *) = nullptr;

  GpuRuntime() { load(); }

  void load() {
    // Search order: explicit override, next to the calling binary, the
    // wheel's bundled library folder, then the default loader search path.
    std::string candidates[5];
    int count = 0;
    if (const char *override = std::getenv("VIENNALS_GPU_LIBRARY"))
      candidates[count++] = override;
    const std::string moduleDir = currentModuleDirectory();
    if (!moduleDir.empty()) {
      // Next to the caller (Python package layout), then the wheel's bundled
      // library folder, then ../lib for a bin/ + lib/ install prefix.
      candidates[count++] = moduleDir + pathSeparator() + gpuLibraryName();
      candidates[count++] = moduleDir + pathSeparator() + ".." +
                            pathSeparator() + "viennals.libs" +
                            pathSeparator() + gpuLibraryName();
      candidates[count++] = moduleDir + pathSeparator() + ".." +
                            pathSeparator() + "lib" + pathSeparator() +
                            gpuLibraryName();
    }
    candidates[count++] = gpuLibraryName();

    LibraryHandle lib = nullptr;
    std::string attempts;
    for (int i = 0; i < count && !lib; ++i) {
      lib = openLibrary(candidates[i].c_str());
      if (!lib) {
        // Report every attempt. The interesting failure is usually not the
        // last one: when the library is found but the CUDA runtime is not,
        // that candidate reports "libcudart.so.12: cannot open shared object
        // file", while the final bare-soname attempt only reports that the
        // library itself is missing. Keeping all of them means the message
        // names the real cause instead of the last symptom.
        attempts += "\n  " + candidates[i] + ": " + lastLoadError();
      }
    }

    if (!lib) {
      status = "ViennaLS_GPU could not be loaded, so the GPU solver is "
               "unavailable; using the CPU solver instead. Attempts:" +
               attempts;
      return;
    }

    if (!resolveAll(lib)) {
      status = "ViennaLS_GPU is missing expected entry points; it is probably "
               "from a different ViennaLS version. Using the CPU solver "
               "instead.";
      return;
    }

    if (abiVersion() != VIENNALS_GPU_ABI_VERSION) {
      status = "ViennaLS_GPU reports ABI version " +
               std::to_string(abiVersion()) + " but ViennaLS expects " +
               std::to_string(VIENNALS_GPU_ABI_VERSION) +
               ". Using the CPU solver instead.";
      return;
    }

    available = true;
    status = "ViennaLS_GPU loaded.";
  }

  /// Resolves every entry point; all-or-nothing so we never call through a
  /// null pointer once `available` is set.
  bool resolveAll(LibraryHandle lib) {
    return resolve(lib, "viennalsGpuAbiVersion", abiVersion) &&
           resolve(lib, "viennalsGpuGetLastErrorMessage", lastErrorMessage) &&
           resolve(lib, "viennalsGpuDiffusionSolve", diffusionSolve) &&
           resolve(lib, "viennalsGpuMacSetup", macSetup) &&
           resolve(lib, "viennalsGpuMacSolve", macSolve) &&
           resolve(lib, "viennalsGpuMacFree", macFree);
  }

  template <class Fn>
  static bool resolve(LibraryHandle lib, const char *name, Fn &target) {
    target = reinterpret_cast<Fn>(findSymbol(lib, name));
    return target != nullptr;
  }
};

/// The single process-wide runtime instance.  Initialization of a
/// function-local static is thread-safe since C++11, so the dlopen happens
/// exactly once even under concurrent first use.
inline const GpuRuntime &runtime() {
  static GpuRuntime instance;
  return instance;
}

} // namespace detail

/// True when the GPU library and its CUDA runtime were loaded successfully.
/// Callers that want to warn before doing GPU-specific setup can check this.
inline bool gpuRuntimeAvailable() { return detail::runtime().available; }

/// Human-readable reason describing why the GPU runtime is (un)available.
inline const char *gpuRuntimeStatusMessage() {
  return detail::runtime().status.c_str();
}

// Human-readable detail for the last GPU failure on this thread.  Falls back
// to the loader status when the library never loaded at all.
inline const char *gpuGetLastErrorMessage() {
  const auto &rt = detail::runtime();
  return rt.available ? rt.lastErrorMessage() : rt.status.c_str();
}

// The diffusion solve of OxidationDiffusion (see viennalsGpuDiffusionSolve):
// 1 when it ran and gave finite values, 0 when the GPU is unavailable (x
// untouched), -1 when it failed.
inline int gpuDiffusionSolve(std::int64_t n, const std::int64_t *ptr,
                             const std::int64_t *col, const double *val,
                             const double *rhs, double *x,
                             unsigned maxIterations, double absTolerance,
                             unsigned &outIterations) {
  const auto &rt = detail::runtime();
  if (!rt.available)
    return 0;
  return rt.diffusionSolve(n, ptr, col, val, rhs, x, maxIterations,
                           absTolerance, &outIterations);
}

// The 3-D velocity-pressure solve of OxidationDeformation: preconditioner
// setup for the equilibrated CSR system (see viennalsGpuMacSetup), nullptr
// when the GPU is unavailable or the setup fails.
inline void *gpuMacSetup(std::int64_t n, const std::int64_t *ptr,
                         const std::int64_t *col, const double *val,
                         std::int64_t nvel, int dim, const double *nullspace) {
  const auto &rt = detail::runtime();
  return rt.available ? rt.macSetup(n, ptr, col, val, nvel, dim, nullspace)
                      : nullptr;
}

// FGMRES with a setup's preconditioner on A x = rhs (x: initial guess on
// entry, solution on exit). False when the solve did not run or gave
// non-finite values.
inline bool gpuMacSolve(void *handle, std::int64_t n, const std::int64_t *ptr,
                        const std::int64_t *col, const double *val,
                        const double *rhs, double *x, unsigned maxIterations,
                        double tolerance, unsigned &outIterations,
                        double &outResidual) {
  const auto &rt = detail::runtime();
  return rt.available &&
         rt.macSolve(handle, n, ptr, col, val, rhs, x, maxIterations,
                     tolerance, &outIterations, &outResidual) != 0;
}

// Free a setup. Safe to call with nullptr.
inline void gpuMacFree(void *handle) {
  const auto &rt = detail::runtime();
  if (rt.available)
    rt.macFree(handle);
}

} // namespace gpu
} // namespace viennals

#endif // VIENNALS_GPU_BICGSTAB
