// Compiled by nvcc.  Provides C++ wrapper functions that delegate to the
// CUDA kernel implementations in lsOxidationBiCGSTAB.cuh.
//
// This translation unit is the ONLY place where lsOxidationBiCGSTAB.cuh
// is included; no .hpp/.cpp file compiled by g++ ever sees the __global__
// kernel syntax.

#include <lsOxidationBiCGSTAB.cuh>
#include <lsOxidationBiCGSTABInterface.hpp>

#include <cstdlib>
#include <cstdio>
#include <new>
#include <string>

namespace viennals {
namespace gpu {

namespace {
thread_local std::string lastGpuError;

void setLastGpuError(const std::string& message) {
    lastGpuError = message;
    fprintf(stderr, "GpuBiCGSTABBuffers: %s\n", message.c_str());
}

void clearLastGpuError() {
    lastGpuError.clear();
}
} // namespace

GpuBiCGSTABBuffers* allocGpuBuffers(uint32_t n, int nFaces,
                                    bool useIlu0Preconditioner) {
    clearLastGpuError();
    // Bail out if no device exists or if a CUDA context cannot actually be
    // created (driver mismatch, permissions, etc.).  cudaGetDeviceCount only
    // enumerates devices; cudaFree(nullptr) forces real context initialization.
    int deviceCount = 0;
    const cudaError_t countErr = cudaGetDeviceCount(&deviceCount);
    if (countErr != cudaSuccess) {
        setLastGpuError(std::string("cudaGetDeviceCount failed: ") +
                        cudaGetErrorString(countErr));
        cudaGetLastError();
        return nullptr;
    }
    if (deviceCount == 0) {
        setLastGpuError("cudaGetDeviceCount reported zero CUDA devices");
        return nullptr;
    }
    int deviceId = 0;
    if (cudaGetDevice(&deviceId) != cudaSuccess || deviceId < 0 ||
        deviceId >= deviceCount) {
        cudaGetLastError();
        deviceId = 0;
    }
    const cudaError_t setErr = cudaSetDevice(deviceId);
    if (setErr != cudaSuccess) {
        setLastGpuError(std::string("cudaSetDevice(") +
                        std::to_string(deviceId) + ") failed: " +
                        cudaGetErrorString(setErr));
        cudaGetLastError();
        return nullptr;
    }
    const cudaError_t ctxErr = cudaFree(nullptr);
    if (ctxErr != cudaSuccess) {
        setLastGpuError(std::string("CUDA context initialization failed on "
                                    "device ") +
                        std::to_string(deviceId) + ": " +
                        cudaGetErrorString(ctxErr));
        cudaGetLastError();  // consume so later CUDA calls start clean
        return nullptr;
    }

    auto* b = new (std::nothrow) GpuBiCGSTABBuffers();
    if (!b) return nullptr;
    b->deviceId = deviceId;
    b->allocate(n, nFaces);
    b->useIlu0Preconditioner = useIlu0Preconditioner;
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

const char* gpuGetLastErrorMessage() {
    return lastGpuError.c_str();
}

bool gpuUploadNeighborIds(GpuBiCGSTABBuffers* gpu,
                          const uint32_t* nb,
                          std::size_t count) {
    const std::size_t expected =
        static_cast<std::size_t>(gpu->n) * static_cast<std::size_t>(gpu->nFaces);
    if (count != expected) {
        setLastGpuError(std::string("neighbor ID upload length mismatch (got ") +
                        std::to_string(count) + ", expected " +
                        std::to_string(expected) + ")");
        return false;
    }
    if (!gpu->activateDevice("gpuUploadNeighborIds"))
        return false;
    return gpu->checkCuda(
        cudaMemcpy(gpu->d_nb, nb, count * sizeof(uint32_t),
                   cudaMemcpyHostToDevice),
        "cudaMemcpy(d_nb)");
}

bool gpuUploadSolverArrays(GpuBiCGSTABBuffers* gpu,
                           const double* diag,
                           const double* b,
                           const double* coeff,
                           uint32_t diagLen,
                           std::size_t coeffLen) {
    const std::size_t expectedCoeff =
        static_cast<std::size_t>(gpu->n) * static_cast<std::size_t>(gpu->nFaces);
    if (diagLen != gpu->n || coeffLen != expectedCoeff) {
        setLastGpuError(std::string("solver array upload length mismatch "
                                    "(diag got ") +
                        std::to_string(diagLen) + ", expected " +
                        std::to_string(gpu->n) + "; coeff got " +
                        std::to_string(coeffLen) + ", expected " +
                        std::to_string(expectedCoeff) + ")");
        return false;
    }
    if (!gpu->activateDevice("gpuUploadSolverArrays"))
        return false;
    if (!gpu->checkCuda(cudaMemcpy(gpu->d_diag, diag, diagLen * sizeof(double),
                                   cudaMemcpyHostToDevice),
                        "cudaMemcpy(d_diag)"))
        return false;
    if (!gpu->checkCuda(cudaMemcpy(gpu->d_b, b, diagLen * sizeof(double),
                                   cudaMemcpyHostToDevice),
                        "cudaMemcpy(d_b)"))
        return false;
    if (!gpu->checkCuda(cudaMemcpy(gpu->d_coeff, coeff,
                                   coeffLen * sizeof(double),
                                   cudaMemcpyHostToDevice),
                        "cudaMemcpy(d_coeff)"))
        return false;
    if (!gpu->useIlu0Preconditioner)
        return true;
    // Refresh ILU(0) values and re-factorize whenever the coefficients change.
    if (gpu->iluReady)
        return gpu->fillAndFactorizeILU(diag, coeff, gpu->nFaces);
    return false;
}

bool gpuSetupCSR(GpuBiCGSTABBuffers* gpu,
                 const uint32_t* h_nb,
                 uint32_t n,
                 int nFaces) {
    if (!gpu->activateDevice("gpuSetupCSR"))
        return false;
    return gpu->setupCSR(h_nb, n, nFaces);
}

bool gpuSolveBiCGSTAB(GpuBiCGSTABBuffers* gpu,
                      double*  x,
                      double   diagEps,
                      unsigned maxIter,
                      double   tolerance,
                      unsigned& outIterations,
                      double&   outResidual) {
    if (!gpu->activateDevice("gpuSolveBiCGSTAB"))
        return false;
    std::vector<double> xVec(x, x + gpu->n);
    const bool ok = solveBiCGSTAB(*gpu, xVec, diagEps, maxIter, tolerance,
                                  outIterations, outResidual);
    if (ok)
        std::copy(xVec.begin(), xVec.end(), x);
    return ok;
}

} // namespace gpu
} // namespace viennals
