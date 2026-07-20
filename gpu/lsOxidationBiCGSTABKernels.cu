// Compiled by nvcc into the ViennaLS_GPU *shared* library.  Provides the
// stable C ABI declared in lsOxidationBiCGSTABAbi.hpp, delegating to the CUDA
// kernel implementations in lsOxidationBiCGSTAB.cuh.
//
// This translation unit is the ONLY place where lsOxidationBiCGSTAB.cuh
// is included; no .hpp/.cpp file compiled by g++ ever sees the __global__
// kernel syntax.
//
// This library is the only binary that links against the CUDA runtime.
// ViennaLS opens it lazily with dlopen() (see lsOxidationBiCGSTABInterface.hpp)
// so that a GPU-enabled build still loads where no CUDA runtime is installed.
// Do NOT include lsOxidationBiCGSTABInterface.hpp here: that header defines
// inline dispatch functions for the *caller* side.

#include <lsOxidationBiCGSTAB.cuh>
#include <lsOxidationBiCGSTABAbi.hpp>

#include <cstdio>
#include <cstdlib>
#include <new>
#include <string>

namespace {

using viennals::gpu::GpuBiCGSTABBuffers;

thread_local std::string lastGpuError;

void setLastGpuError(const std::string &message) {
  lastGpuError = message;
  fprintf(stderr, "GpuBiCGSTABBuffers: %s\n", message.c_str());
}

void clearLastGpuError() { lastGpuError.clear(); }

GpuBiCGSTABBuffers *asBuffers(void *handle) {
  return static_cast<GpuBiCGSTABBuffers *>(handle);
}

const GpuBiCGSTABBuffers *asBuffers(const void *handle) {
  return static_cast<const GpuBiCGSTABBuffers *>(handle);
}

} // namespace

extern "C" {

int viennalsGpuAbiVersion(void) { return VIENNALS_GPU_ABI_VERSION; }

void *viennalsGpuAllocBuffers(uint32_t n, int nFaces,
                              int useIlu0Preconditioner) {
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
    setLastGpuError(std::string("cudaSetDevice(") + std::to_string(deviceId) +
                    ") failed: " + cudaGetErrorString(setErr));
    cudaGetLastError();
    return nullptr;
  }
  const cudaError_t ctxErr = cudaFree(nullptr);
  if (ctxErr != cudaSuccess) {
    setLastGpuError(std::string("CUDA context initialization failed on "
                                "device ") +
                    std::to_string(deviceId) + ": " +
                    cudaGetErrorString(ctxErr));
    cudaGetLastError(); // consume so later CUDA calls start clean
    return nullptr;
  }

  auto *b = new (std::nothrow) GpuBiCGSTABBuffers();
  if (!b)
    return nullptr;
  b->deviceId = deviceId;
  b->allocate(n, nFaces);
  b->useIlu0Preconditioner = useIlu0Preconditioner != 0;
  if (!b->valid) {
    delete b;
    return nullptr;
  }
  return b;
}

void viennalsGpuFreeBuffers(void *handle) {
  delete asBuffers(handle); // calls ~GpuBiCGSTABBuffers() -> gpu->free()
}

int viennalsGpuIsValid(const void *handle) {
  const auto *gpu = asBuffers(handle);
  return (gpu != nullptr && gpu->valid) ? 1 : 0;
}

const char *viennalsGpuGetLastErrorMessage(void) { return lastGpuError.c_str(); }

int viennalsGpuUploadNeighborIds(void *handle, const uint32_t *nb,
                                 std::size_t count) {
  auto *gpu = asBuffers(handle);
  const std::size_t expected =
      static_cast<std::size_t>(gpu->n) * static_cast<std::size_t>(gpu->nFaces);
  if (count != expected) {
    setLastGpuError(std::string("neighbor ID upload length mismatch (got ") +
                    std::to_string(count) + ", expected " +
                    std::to_string(expected) + ")");
    return 0;
  }
  if (!gpu->activateDevice("gpuUploadNeighborIds"))
    return 0;
  return gpu->checkCuda(cudaMemcpy(gpu->d_nb, nb, count * sizeof(uint32_t),
                                   cudaMemcpyHostToDevice),
                        "cudaMemcpy(d_nb)")
             ? 1
             : 0;
}

int viennalsGpuUploadSolverArrays(void *handle, const double *diag,
                                  const double *b, const double *coeff,
                                  uint32_t diagLen, std::size_t coeffLen) {
  auto *gpu = asBuffers(handle);
  const std::size_t expectedCoeff =
      static_cast<std::size_t>(gpu->n) * static_cast<std::size_t>(gpu->nFaces);
  if (diagLen != gpu->n || coeffLen != expectedCoeff) {
    setLastGpuError(std::string("solver array upload length mismatch "
                                "(diag got ") +
                    std::to_string(diagLen) + ", expected " +
                    std::to_string(gpu->n) + "; coeff got " +
                    std::to_string(coeffLen) + ", expected " +
                    std::to_string(expectedCoeff) + ")");
    return 0;
  }
  if (!gpu->activateDevice("gpuUploadSolverArrays"))
    return 0;
  if (!gpu->checkCuda(cudaMemcpy(gpu->d_diag, diag, diagLen * sizeof(double),
                                 cudaMemcpyHostToDevice),
                      "cudaMemcpy(d_diag)"))
    return 0;
  if (!gpu->checkCuda(cudaMemcpy(gpu->d_b, b, diagLen * sizeof(double),
                                 cudaMemcpyHostToDevice),
                      "cudaMemcpy(d_b)"))
    return 0;
  if (!gpu->checkCuda(cudaMemcpy(gpu->d_coeff, coeff,
                                 coeffLen * sizeof(double),
                                 cudaMemcpyHostToDevice),
                      "cudaMemcpy(d_coeff)"))
    return 0;
  if (!gpu->useIlu0Preconditioner)
    return 1;
  // Refresh ILU(0) values and re-factorize whenever the coefficients change.
  if (gpu->iluReady)
    return gpu->fillAndFactorizeILU(diag, coeff, gpu->nFaces) ? 1 : 0;
  return 0;
}

int viennalsGpuUploadRhs(void *handle, const double *b, uint32_t n) {
  auto *gpu = asBuffers(handle);
  if (n != gpu->n) {
    setLastGpuError(std::string("RHS upload length mismatch (got ") +
                    std::to_string(n) + ", expected " +
                    std::to_string(gpu->n) + ")");
    return 0;
  }
  if (!gpu->activateDevice("gpuUploadRhs"))
    return 0;
  return gpu->checkCuda(cudaMemcpy(gpu->d_b, b,
                                   static_cast<std::size_t>(n) * sizeof(double),
                                   cudaMemcpyHostToDevice),
                        "cudaMemcpy(d_b)")
             ? 1
             : 0;
}

int viennalsGpuSetupCSR(void *handle, const uint32_t *hNb, uint32_t n,
                        int nFaces) {
  auto *gpu = asBuffers(handle);
  if (!gpu->activateDevice("gpuSetupCSR"))
    return 0;
  return gpu->setupCSR(hNb, n, nFaces) ? 1 : 0;
}

int viennalsGpuSolveBiCGSTAB(void *handle, double *x, double diagEps,
                             unsigned maxIter, double tolerance,
                             unsigned *outIterations, double *outResidual) {
  auto *gpu = asBuffers(handle);
  if (!gpu->activateDevice("gpuSolveBiCGSTAB"))
    return 0;
  std::vector<double> xVec(x, x + gpu->n);
  const bool ok = viennals::gpu::solveBiCGSTAB(
      *gpu, xVec, diagEps, maxIter, tolerance, *outIterations, *outResidual);
  if (ok)
    std::copy(xVec.begin(), xVec.end(), x);
  return ok ? 1 : 0;
}

} // extern "C"
