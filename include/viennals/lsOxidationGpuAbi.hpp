// Stable C ABI exported by the ViennaLS_GPU shared library: the linear solves
// of the oxidation solvers on the GPU (gpu/lsOxidationGpu.cu).
//
// The library links against the CUDA runtime.  ViennaLS itself must NOT, so
// that a build with GPU support still loads on machines with no CUDA runtime
// installed.  The library is therefore opened lazily with dlopen() and these
// symbols are resolved with dlsym() -- see lsOxidationGpuInterface.hpp.
//
// Everything here is deliberately plain C: opaque void* handles, no
// references, int instead of bool.  That keeps the boundary independent of
// C++ name mangling and of the standard library, so the two sides stay
// compatible even when built by different compilers (nvcc vs g++).
//
// Bump VIENNALS_GPU_ABI_VERSION on any incompatible change to the signatures
// below; the loader refuses a library reporting a different version.

#pragma once

#include <cstddef>
#include <cstdint>

#define VIENNALS_GPU_ABI_VERSION 3

// The exporting side (lsOxidationGpu.cu) must keep these visible even under
// -fvisibility=hidden, or dlsym() cannot find them.
#if defined(_WIN32)
#define VIENNALS_GPU_ABI __declspec(dllexport)
#else
#define VIENNALS_GPU_ABI __attribute__((visibility("default")))
#endif

extern "C" {

/// Returns the ABI version this library was built against.
VIENNALS_GPU_ABI int viennalsGpuAbiVersion(void);

/// Human-readable detail for the last failure on this thread.
VIENNALS_GPU_ABI const char *viennalsGpuGetLastErrorMessage(void);

/// Diffusion solve of OxidationDiffusion, as on the CPU: BiCGSTAB
/// preconditioned by smoothed-aggregation AMG with an ILU0 smoother on the
/// CSR system (ptr, col, val) of `n` unknowns, stopping at
/// ||r||_2 <= absTolerance or after maxIterations.  `x` is the initial guess
/// on entry, the solution on exit.  Returns 1 when the solve ran and gave
/// finite values (whatever the residual), 0 when no CUDA device is usable
/// (`x` untouched), -1 when it failed.
VIENNALS_GPU_ABI int viennalsGpuDiffusionSolve(int64_t n, const int64_t *ptr,
                                               const int64_t *col,
                                               const double *val,
                                               const double *rhs, double *x,
                                               unsigned maxIterations,
                                               double absTolerance,
                                               unsigned *outIterations);

/// 3-D velocity-pressure solve of OxidationDeformation: set up the
/// preconditioner for the CSR system (ptr, col, val) of `n` unknowns, the
/// first `nvel` of them velocities; `nullspace` holds `dim` near-nullspace
/// values per velocity unknown.  Returns nullptr on failure.
VIENNALS_GPU_ABI void *viennalsGpuMacSetup(int64_t n, const int64_t *ptr,
                                           const int64_t *col,
                                           const double *val, int64_t nvel,
                                           int dim, const double *nullspace);

/// FGMRES on A x = rhs with the preconditioner of `handle`; A may differ from
/// the matrix of the setup (same unknowns).  `x` is the initial guess on
/// entry, the solution on exit.  Non-zero when the solve ran and gave finite
/// values, whatever the residual; iterations and relative residual returned.
VIENNALS_GPU_ABI int viennalsGpuMacSolve(void *handle, int64_t n,
                                         const int64_t *ptr, const int64_t *col,
                                         const double *val, const double *rhs,
                                         double *x, unsigned maxIterations,
                                         double tolerance,
                                         unsigned *outIterations,
                                         double *outResidual);

/// Free a setup.  Safe to call with nullptr.
VIENNALS_GPU_ABI void viennalsGpuMacFree(void *handle);

} // extern "C"
