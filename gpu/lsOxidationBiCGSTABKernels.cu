// Compiled by nvcc.  Provides C++ wrapper functions that delegate to the
// CUDA kernel implementations in lsOxidationBiCGSTAB.cuh.
//
// This translation unit is the ONLY place where lsOxidationBiCGSTAB.cuh
// is included; no .hpp/.cpp file compiled by g++ ever sees the __global__
// kernel syntax.

#include <lsOxidationBiCGSTAB.cuh>
#include <lsOxidationBiCGSTABInterface.hpp>

#include <cstdlib>
#include <new>

namespace viennals {
namespace gpu {

GpuBiCGSTABBuffers* allocGpuBuffers(uint32_t n, int nFaces) {
    auto* b = new (std::nothrow) GpuBiCGSTABBuffers();
    if (!b) return nullptr;
    b->allocate(n, nFaces);
    if (!b->valid) {
        delete b;
        return nullptr;
    }
    return b;
}

void freeGpuBuffers(GpuBiCGSTABBuffers* gpu) {
    delete gpu; // calls ~GpuBiCGSTABBuffers() which calls gpu->free()
}

bool gpuIsValid(const GpuBiCGSTABBuffers* gpu) {
    return gpu != nullptr && gpu->valid;
}

void gpuUploadNeighborIds(GpuBiCGSTABBuffers* gpu,
                          const uint32_t* nb,
                          std::size_t count) {
    LS_CUDA_CHECK(cudaMemcpy(gpu->d_nb, nb, count * sizeof(uint32_t),
                             cudaMemcpyHostToDevice));
}

void gpuUploadSolverArrays(GpuBiCGSTABBuffers* gpu,
                           const float* diag,
                           const float* b,
                           const float* coeff,
                           uint32_t diagLen,
                           std::size_t coeffLen) {
    LS_CUDA_CHECK(cudaMemcpy(gpu->d_diag, diag, diagLen  * sizeof(float),
                             cudaMemcpyHostToDevice));
    LS_CUDA_CHECK(cudaMemcpy(gpu->d_b,    b,    diagLen  * sizeof(float),
                             cudaMemcpyHostToDevice));
    LS_CUDA_CHECK(cudaMemcpy(gpu->d_coeff, coeff, coeffLen * sizeof(float),
                             cudaMemcpyHostToDevice));
}

void gpuSolveBiCGSTAB(GpuBiCGSTABBuffers* gpu,
                      float*   x,
                      float    diagEps,
                      unsigned maxIter,
                      float    tolerance,
                      unsigned& outIterations,
                      double&   outResidual) {
    std::vector<float> xVec(x, x + gpu->n);
    solveBiCGSTAB(*gpu, xVec, diagEps, maxIter, tolerance,
                  outIterations, outResidual);
    std::copy(xVec.begin(), xVec.end(), x);
}

} // namespace gpu
} // namespace viennals
