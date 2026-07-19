// C++ (g++) interface to the GPU BiCGSTAB solver.
//
// This header is safe to include from any .cpp file compiled by g++.
// It only forward-declares GpuBiCGSTABBuffers (opaque handle) and
// declares free functions that are implemented in
// ViennaLS_GPU (lsOxidationBiCGSTABKernels.cu, compiled by nvcc).
//
// The actual CUDA kernels live in lsOxidationBiCGSTAB.cuh.  That file
// must never be included from a .cpp compiled by g++ — only from .cu files.
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
// If it (or the CUDA runtime it needs) cannot be loaded, allocGpuBuffers()
// returns nullptr and callers transparently fall back to the CPU solver, with
// the reason available from gpuGetLastErrorMessage().

#pragma once

#ifdef VIENNALS_GPU_BICGSTAB

#include "lsOxidationBiCGSTABAbi.hpp"

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

// Sentinel value used in the neighbor-ID array to mark boundary / out-of-bounds
// faces (must match the constant in lsOxidationBiCGSTAB.cuh).
static constexpr uint32_t kNoNode = 0xFFFFFFFFu;

// Opaque handle — complete definition is in lsOxidationBiCGSTAB.cuh /
// lsOxidationBiCGSTABKernels.cu.  Consumers hold a raw pointer only.
struct GpuBiCGSTABBuffers;

namespace detail {

#if defined(_WIN32)
using LibraryHandle = HMODULE;
inline LibraryHandle openLibrary(const char *path) { return LoadLibraryA(path); }
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
  void *(*allocBuffers)(uint32_t, int, int) = nullptr;
  void (*freeBuffers)(void *) = nullptr;
  int (*isValid)(const void *) = nullptr;
  const char *(*lastErrorMessage)(void) = nullptr;
  int (*uploadNeighborIds)(void *, const uint32_t *, std::size_t) = nullptr;
  int (*setupCSR)(void *, const uint32_t *, uint32_t, int) = nullptr;
  int (*uploadSolverArrays)(void *, const double *, const double *,
                            const double *, uint32_t, std::size_t) = nullptr;
  int (*uploadRhs)(void *, const double *, uint32_t) = nullptr;
  int (*solveBiCGSTAB)(void *, double *, double, unsigned, double, unsigned *,
                       double *) = nullptr;

  GpuRuntime() { load(); }

  void load() {
    // Search order: explicit override, next to the calling binary, the
    // wheel's bundled library folder, then the default loader search path.
    std::string candidates[4];
    int count = 0;
    if (const char *override = std::getenv("VIENNALS_GPU_LIBRARY"))
      candidates[count++] = override;
    const std::string moduleDir = currentModuleDirectory();
    if (!moduleDir.empty()) {
      candidates[count++] = moduleDir + pathSeparator() + gpuLibraryName();
      candidates[count++] = moduleDir + pathSeparator() + ".." +
                            pathSeparator() + "viennals.libs" +
                            pathSeparator() + gpuLibraryName();
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
           resolve(lib, "viennalsGpuAllocBuffers", allocBuffers) &&
           resolve(lib, "viennalsGpuFreeBuffers", freeBuffers) &&
           resolve(lib, "viennalsGpuIsValid", isValid) &&
           resolve(lib, "viennalsGpuGetLastErrorMessage", lastErrorMessage) &&
           resolve(lib, "viennalsGpuUploadNeighborIds", uploadNeighborIds) &&
           resolve(lib, "viennalsGpuSetupCSR", setupCSR) &&
           resolve(lib, "viennalsGpuUploadSolverArrays", uploadSolverArrays) &&
           resolve(lib, "viennalsGpuUploadRhs", uploadRhs) &&
           resolve(lib, "viennalsGpuSolveBiCGSTAB", solveBiCGSTAB);
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

inline void *toHandle(GpuBiCGSTABBuffers *gpu) {
  return reinterpret_cast<void *>(gpu);
}
inline const void *toHandle(const GpuBiCGSTABBuffers *gpu) {
  return reinterpret_cast<const void *>(gpu);
}

} // namespace detail

/// True when the GPU library and its CUDA runtime were loaded successfully.
/// Callers that want to warn before doing GPU-specific setup can check this.
inline bool gpuRuntimeAvailable() { return detail::runtime().available; }

/// Human-readable reason describing why the GPU runtime is (un)available.
inline const char *gpuRuntimeStatusMessage() {
  return detail::runtime().status.c_str();
}

// Allocate GPU buffers for a solver with `n` nodes and `nFaces` (2*D) faces.
// Returns nullptr if CUDA is unavailable.
inline GpuBiCGSTABBuffers *allocGpuBuffers(uint32_t n, int nFaces,
                                           bool useIlu0Preconditioner) {
  const auto &rt = detail::runtime();
  if (!rt.available)
    return nullptr;
  return reinterpret_cast<GpuBiCGSTABBuffers *>(
      rt.allocBuffers(n, nFaces, useIlu0Preconditioner ? 1 : 0));
}

// Free previously allocated GPU buffers.  Safe to call with nullptr.
inline void freeGpuBuffers(GpuBiCGSTABBuffers *gpu) {
  const auto &rt = detail::runtime();
  if (rt.available)
    rt.freeBuffers(detail::toHandle(gpu));
}

// Human-readable detail for the last GPU wrapper failure on this thread.
// Falls back to the loader status when the library never loaded at all.
inline const char *gpuGetLastErrorMessage() {
  const auto &rt = detail::runtime();
  return rt.available ? rt.lastErrorMessage() : rt.status.c_str();
}

// Is the buffer handle valid (non-null and successfully allocated)?
inline bool gpuIsValid(const GpuBiCGSTABBuffers *gpu) {
  const auto &rt = detail::runtime();
  return rt.available && rt.isValid(detail::toHandle(gpu)) != 0;
}

// Upload geometry-fixed neighbor-ID array (face-major, kNoNode = 0xFFFFFFFF).
// `count` must equal nFaces * n.
inline bool gpuUploadNeighborIds(GpuBiCGSTABBuffers *gpu, const uint32_t *nb,
                                 std::size_t count) {
  const auto &rt = detail::runtime();
  return rt.available &&
         rt.uploadNeighborIds(detail::toHandle(gpu), nb, count) != 0;
}

// Build the CSR sparsity pattern from h_nb (face-major, length nFaces*n),
// upload to the device, and run CUSPARSE symbolic analysis for ILU(0) and
// the two triangular solves.  Must be called after gpuUploadNeighborIds and
// before the first gpuUploadSolverArrays / gpuSolveBiCGSTAB call.
inline bool gpuSetupCSR(GpuBiCGSTABBuffers *gpu, const uint32_t *h_nb,
                        uint32_t n, int nFaces) {
  const auto &rt = detail::runtime();
  return rt.available &&
         rt.setupCSR(detail::toHandle(gpu), h_nb, n, nFaces) != 0;
}

// Upload per-solve arrays (diag, b, faceCoeffs) and re-factorize ILU(0).
// `diagLen` == n, `coeffLen` == nFaces * n.
inline bool gpuUploadSolverArrays(GpuBiCGSTABBuffers *gpu, const double *diag,
                                  const double *b, const double *coeff,
                                  uint32_t diagLen, std::size_t coeffLen) {
  const auto &rt = detail::runtime();
  return rt.available && rt.uploadSolverArrays(detail::toHandle(gpu), diag, b,
                                               coeff, diagLen, coeffLen) != 0;
}

// Upload only the RHS vector (d_b).  Use when the matrix geometry is already
// uploaded and only the right-hand side changes (e.g. successive Stokes
// component solves that share the same stiffness matrix).
inline bool gpuUploadRhs(GpuBiCGSTABBuffers *gpu, const double *b, uint32_t n) {
  const auto &rt = detail::runtime();
  return rt.available && rt.uploadRhs(detail::toHandle(gpu), b, n) != 0;
}

// Run GPU BiCGSTAB.
//   x (length n, host): initial guess on entry, solution on exit.
//   outResidual is the raw (unnormalized) max-abs residual on exit.
// Returns true only when the GPU solve converged and produced finite values.
inline bool gpuSolveBiCGSTAB(GpuBiCGSTABBuffers *gpu, double *x, double diagEps,
                             unsigned maxIter, double tolerance,
                             unsigned &outIterations, double &outResidual) {
  const auto &rt = detail::runtime();
  return rt.available &&
         rt.solveBiCGSTAB(detail::toHandle(gpu), x, diagEps, maxIter, tolerance,
                          &outIterations, &outResidual) != 0;
}

} // namespace gpu
} // namespace viennals

#endif // VIENNALS_GPU_BICGSTAB
