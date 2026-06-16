// C++ (g++) interface to the GPU BiCGSTAB solver.
//
// This header is safe to include from any .cpp file compiled by g++.
// It only forward-declares GpuBiCGSTABBuffers (opaque handle) and
// declares free functions that are implemented in
// ViennaLS_GPU (lsOxidationBiCGSTABKernels.cu, compiled by nvcc).
//
// The actual CUDA kernels live in lsOxidationBiCGSTAB.cuh.  That file
// must never be included from a .cpp compiled by g++ — only from .cu files.

#pragma once

#ifdef VIENNALS_GPU_BICGSTAB

#include <cstddef>
#include <cstdint>

namespace viennals {
namespace gpu {

// Sentinel value used in the neighbor-ID array to mark boundary / out-of-bounds
// faces (must match the constant in lsOxidationBiCGSTAB.cuh).
static constexpr uint32_t kNoNode = 0xFFFFFFFFu;

// Opaque handle — complete definition is in lsOxidationBiCGSTAB.cuh /
// lsOxidationBiCGSTABKernels.cu.  Consumers hold a raw pointer only.
struct GpuBiCGSTABBuffers;

// Allocate GPU buffers for a solver with `n` nodes and `nFaces` (2*D) faces.
// Returns nullptr if CUDA is unavailable.
GpuBiCGSTABBuffers *allocGpuBuffers(uint32_t n, int nFaces,
                                    bool useIlu0Preconditioner);

// Free previously allocated GPU buffers.  Safe to call with nullptr.
void freeGpuBuffers(GpuBiCGSTABBuffers *gpu);

// Human-readable detail for the last GPU wrapper failure on this thread.
const char *gpuGetLastErrorMessage();

// Is the buffer handle valid (non-null and successfully allocated)?
bool gpuIsValid(const GpuBiCGSTABBuffers *gpu);

// Upload geometry-fixed neighbor-ID array (face-major, kNoNode = 0xFFFFFFFF).
// `count` must equal nFaces * n.
bool gpuUploadNeighborIds(GpuBiCGSTABBuffers *gpu, const uint32_t *nb,
                          std::size_t count);

// Build the CSR sparsity pattern from h_nb (face-major, length nFaces*n),
// upload to the device, and run CUSPARSE symbolic analysis for ILU(0) and
// the two triangular solves.  Must be called after gpuUploadNeighborIds and
// before the first gpuUploadSolverArrays / gpuSolveBiCGSTAB call.
bool gpuSetupCSR(GpuBiCGSTABBuffers *gpu, const uint32_t *h_nb, uint32_t n,
                 int nFaces);

// Upload per-solve arrays (diag, b, faceCoeffs) and re-factorize ILU(0).
// `diagLen` == n, `coeffLen` == nFaces * n.
bool gpuUploadSolverArrays(GpuBiCGSTABBuffers *gpu, const double *diag,
                           const double *b, const double *coeff,
                           uint32_t diagLen, std::size_t coeffLen);

// Upload only the RHS vector (d_b).  Use when the matrix geometry is already
// uploaded and only the right-hand side changes (e.g. successive Stokes
// component solves that share the same stiffness matrix).
bool gpuUploadRhs(GpuBiCGSTABBuffers *gpu, const double *b, uint32_t n);

// Run GPU BiCGSTAB.
//   x (length n, host): initial guess on entry, solution on exit.
//   outResidual is the raw (unnormalized) max-abs residual on exit.
// Returns true only when the GPU solve converged and produced finite values.
bool gpuSolveBiCGSTAB(GpuBiCGSTABBuffers *gpu, double *x, double diagEps,
                      unsigned maxIter, double tolerance,
                      unsigned &outIterations, double &outResidual);

} // namespace gpu
} // namespace viennals

#endif // VIENNALS_GPU_BICGSTAB
