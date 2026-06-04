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
                           const double* diag,
                           const double* b,
                           const double* coeff,
                           uint32_t diagLen,
                           std::size_t coeffLen) {
    LS_CUDA_CHECK(cudaMemcpy(gpu->d_diag, diag, diagLen  * sizeof(double),
                             cudaMemcpyHostToDevice));
    LS_CUDA_CHECK(cudaMemcpy(gpu->d_b,    b,    diagLen  * sizeof(double),
                             cudaMemcpyHostToDevice));
    LS_CUDA_CHECK(cudaMemcpy(gpu->d_coeff, coeff, coeffLen * sizeof(double),
                             cudaMemcpyHostToDevice));
}

bool gpuSolveBiCGSTAB(GpuBiCGSTABBuffers* gpu,
                      double*  x,
                      double   diagEps,
                      unsigned maxIter,
                      double   tolerance,
                      unsigned& outIterations,
                      double&   outResidual) {
    std::vector<double> xVec(x, x + gpu->n);
    const bool ok = solveBiCGSTAB(*gpu, xVec, diagEps, maxIter, tolerance,
                                  outIterations, outResidual);
    if (ok)
        std::copy(xVec.begin(), xVec.end(), x);
    return ok;
}

} // namespace gpu
} // namespace viennals
