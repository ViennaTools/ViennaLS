// Stable C ABI exported by the ViennaLS_GPU shared library.
//
// The GPU BiCGSTAB solver lives in a separate shared library that links
// against the CUDA runtime.  ViennaLS itself must NOT link against CUDA, so
// that a build with GPU support still loads on machines with no CUDA runtime
// installed.  The library is therefore opened lazily with dlopen() and these
// symbols are resolved with dlsym() -- see lsOxidationBiCGSTABInterface.hpp.
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

#define VIENNALS_GPU_ABI_VERSION 1

// The exporting side (lsOxidationBiCGSTABKernels.cu) must keep these visible
// even under -fvisibility=hidden, or dlsym() cannot find them.
#if defined(_WIN32)
#define VIENNALS_GPU_ABI __declspec(dllexport)
#else
#define VIENNALS_GPU_ABI __attribute__((visibility("default")))
#endif

extern "C" {

/// Returns the ABI version this library was built against.
VIENNALS_GPU_ABI int viennalsGpuAbiVersion(void);

/// Allocate GPU buffers for `n` nodes and `nFaces` (2*D) faces.
/// Returns nullptr if no usable CUDA device or context is available.
VIENNALS_GPU_ABI void *viennalsGpuAllocBuffers(uint32_t n, int nFaces,
                                               int useIlu0Preconditioner);

/// Free previously allocated buffers. Safe to call with nullptr.
VIENNALS_GPU_ABI void viennalsGpuFreeBuffers(void *handle);

/// Non-zero if the handle is valid (non-null and successfully allocated).
VIENNALS_GPU_ABI int viennalsGpuIsValid(const void *handle);

/// Human-readable detail for the last failure on this thread.
VIENNALS_GPU_ABI const char *viennalsGpuGetLastErrorMessage(void);

/// Upload the geometry-fixed neighbor-ID array; `count` must equal nFaces*n.
VIENNALS_GPU_ABI int viennalsGpuUploadNeighborIds(void *handle,
                                                  const uint32_t *nb,
                                                  std::size_t count);

/// Build the CSR sparsity pattern and run the cuSPARSE symbolic analysis.
VIENNALS_GPU_ABI int viennalsGpuSetupCSR(void *handle, const uint32_t *hNb,
                                         uint32_t n, int nFaces);

/// Upload per-solve arrays (diag, b, faceCoeffs) and re-factorize ILU(0).
VIENNALS_GPU_ABI int
viennalsGpuUploadSolverArrays(void *handle, const double *diag, const double *b,
                              const double *coeff, uint32_t diagLen,
                              std::size_t coeffLen);

/// Upload only the RHS vector, when the matrix geometry is already resident.
VIENNALS_GPU_ABI int viennalsGpuUploadRhs(void *handle, const double *b,
                                          uint32_t n);

/// Run GPU BiCGSTAB. `x` is the initial guess on entry, the solution on exit.
/// Non-zero only when the solve converged and produced finite values.
VIENNALS_GPU_ABI int viennalsGpuSolveBiCGSTAB(void *handle, double *x,
                                              double diagEps, unsigned maxIter,
                                              double tolerance,
                                              unsigned *outIterations,
                                              double *outResidual);

} // extern "C"
