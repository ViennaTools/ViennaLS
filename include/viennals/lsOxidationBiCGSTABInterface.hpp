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
GpuBiCGSTABBuffers* allocGpuBuffers(uint32_t n, int nFaces);

// Free previously allocated GPU buffers.  Safe to call with nullptr.
void freeGpuBuffers(GpuBiCGSTABBuffers* gpu);

// Is the buffer handle valid (non-null and successfully allocated)?
bool gpuIsValid(const GpuBiCGSTABBuffers* gpu);

// Upload geometry-fixed neighbor-ID array (face-major, kNoNode = 0xFFFFFFFF).
// `count` must equal nFaces * n.
void gpuUploadNeighborIds(GpuBiCGSTABBuffers* gpu,
                          const uint32_t* nb,
                          std::size_t count);

// Upload per-solve arrays (diag, b, faceCoeffs).
// `diagLen` == n, `coeffLen` == nFaces * n.
void gpuUploadSolverArrays(GpuBiCGSTABBuffers* gpu,
                           const float* diag,
                           const float* b,
                           const float* coeff,
                           uint32_t diagLen,
                           std::size_t coeffLen);

// Run GPU BiCGSTAB.
//   x (length n, host): initial guess on entry, solution on exit.
//   outResidual is the raw (unnormalized) max-abs residual on exit.
void gpuSolveBiCGSTAB(GpuBiCGSTABBuffers* gpu,
                      float*   x,
                      float    diagEps,
                      unsigned maxIter,
                      float    tolerance,
                      unsigned& outIterations,
                      double&   outResidual);

} // namespace gpu
} // namespace viennals

#endif // VIENNALS_GPU_BICGSTAB
