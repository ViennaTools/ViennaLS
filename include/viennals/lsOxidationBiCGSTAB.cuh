// GPU-accelerated BiCGSTAB solver for the oxidation diffusion equation.
//
// Activated by defining VIENNALS_GPU_BICGSTAB (set when VIENNALS_USE_GPU is ON in CMake).
// Requires CUDA Toolkit; does NOT require OptiX.
//
// All kernels use the CUDA Runtime API (cuda_runtime.h) and are templated so
// they compile cleanly in multiple translation units.
//
// Memory layout (face-major): index = faceId * n + nodeId
// This gives coalesced reads when all warp threads process the same face
// of consecutive nodes.

#pragma once

#ifdef VIENNALS_GPU_BICGSTAB

#include <cuda_runtime.h>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace viennals {
namespace gpu {

// kNoNode (0xFFFFFFFFu) is declared in lsOxidationBiCGSTABInterface.hpp and
// used by callers.  Kernels compare directly against the literal to avoid a
// duplicate-definition conflict when both headers are included in a .cu file.

static constexpr int kThreadsPerBlock = 256;

inline int blocksFor(int n) {
    return (n + kThreadsPerBlock - 1) / kThreadsPerBlock;
}

// ─────────────────────────── CUDA error check ────────────────────────────

#ifdef NDEBUG
#define LS_CUDA_CHECK(call) (call)
#else
#define LS_CUDA_CHECK(call)                                                    \
    do {                                                                       \
        cudaError_t _e = (call);                                               \
        if (_e != cudaSuccess)                                                 \
            fprintf(stderr, "CUDA error %s:%d – %s\n", __FILE__, __LINE__,    \
                    cudaGetErrorString(_e));                                    \
    } while (0)
#endif

// ─────────────────────────── kernels ─────────────────────────────────────

// Sparse matrix-vector product:
//   Av[i] = diag[i]*v[i] - sum_{fi: interior face}(coeff[fi*n+i] * v[nb[fi*n+i]])
//
// nFaces is a compile-time constant so the inner loop can be unrolled.
// Boundary faces are represented by nb[fi*n+i] == kNoNode and coeff == 0,
// so the branch is always predicted correctly for regular grids.
template <int nFaces>
__global__ void spmvKernel(
    const float* __restrict__ v,
    const float* __restrict__ diag,
    const uint32_t* __restrict__ nb,
    const float* __restrict__ coeff,
    float* __restrict__ Av,
    uint32_t n)
{
    const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    float res = diag[i] * v[i];
#pragma unroll
    for (int fi = 0; fi < nFaces; ++fi) {
        const uint32_t j = nb[static_cast<uint32_t>(fi) * n + i];
        if (j != 0xFFFFFFFFu)
            res -= coeff[static_cast<uint32_t>(fi) * n + i] * v[j];
    }
    Av[i] = res;
}

// r[i] = b[i] - Av[i]  (initial residual setup)
__global__ void residualKernel(
    const float* __restrict__ b,
    const float* __restrict__ Av,
    float* __restrict__ r,
    uint32_t n)
{
    const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    r[i] = b[i] - Av[i];
}

// p[i] = r[i] + beta*(p[i] - omega*v[i])
__global__ void updatePKernel(
    const float* __restrict__ r,
    float beta, float beta_omega,
    const float* __restrict__ v,
    float* __restrict__ p,
    uint32_t n)
{
    const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    p[i] = r[i] + beta * p[i] - beta_omega * v[i];
}

// out[i] = in[i] - alpha*v[i]  (used for s = r - alpha*v)
__global__ void subtractScaledKernel(
    const float* __restrict__ in,
    float alpha,
    const float* __restrict__ v,
    float* __restrict__ out,
    uint32_t n)
{
    const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = in[i] - alpha * v[i];
}

// Jacobi preconditioner: out[i] = in[i] / diag[i]  (guarded against zero diagonal)
__global__ void jacobiKernel(
    const float* __restrict__ in,
    const float* __restrict__ diag,
    float* __restrict__ out,
    float eps,
    uint32_t n)
{
    const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = (diag[i] > eps) ? in[i] / diag[i] : in[i];
}

// x[i] += alpha*y[i] + omega*z[i]
__global__ void updateXKernel(
    float* __restrict__ x,
    float alpha, const float* __restrict__ y,
    float omega, const float* __restrict__ z,
    uint32_t n)
{
    const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    x[i] += alpha * y[i] + omega * z[i];
}

// r[i] = s[i] - omega*t[i]  (residual update at end of BiCGSTAB iteration)
__global__ void updateRKernel(
    const float* __restrict__ s,
    float omega,
    const float* __restrict__ t,
    float* __restrict__ r,
    uint32_t n)
{
    const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    r[i] = s[i] - omega * t[i];
}

// Dot product: atomicAdd accumulates partial sums into result (must be zeroed before launch).
// Uses stride-loop so a fixed grid size (≤ 1024 blocks) handles any n.
__global__ void dotKernel(
    const float* __restrict__ a,
    const float* __restrict__ b,
    float* __restrict__ result,
    uint32_t n)
{
    extern __shared__ unsigned char smem[];
    float* sdata = reinterpret_cast<float*>(smem);

    const uint32_t tid = threadIdx.x;
    float sum = 0.f;
    for (uint32_t i = blockIdx.x * blockDim.x + tid;
         i < n;
         i += gridDim.x * blockDim.x)
        sum += a[i] * b[i];

    sdata[tid] = sum;
    __syncthreads();

    for (unsigned s = blockDim.x >> 1; s > 0; s >>= 1) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        __syncthreads();
    }

    if (tid == 0) atomicAdd(result, sdata[0]);
}

// Max absolute value using IEEE-754 bit trick:
// for non-negative floats, bit-pattern ordering equals float ordering.
// result (unsigned int*) must be zeroed before launch.
__global__ void maxAbsKernel(
    const float* __restrict__ v,
    unsigned int* __restrict__ result,
    uint32_t n)
{
    extern __shared__ unsigned char smem[];
    unsigned int* sdata = reinterpret_cast<unsigned int*>(smem);

    const uint32_t tid = threadIdx.x;
    unsigned int maxVal = 0u;
    for (uint32_t i = blockIdx.x * blockDim.x + tid;
         i < n;
         i += gridDim.x * blockDim.x) {
        unsigned int bits = __float_as_uint(fabsf(v[i]));
        if (bits > maxVal) maxVal = bits;
    }

    sdata[tid] = maxVal;
    __syncthreads();

    for (unsigned s = blockDim.x >> 1; s > 0; s >>= 1) {
        if (tid < s && sdata[tid + s] > sdata[tid])
            sdata[tid] = sdata[tid + s];
        __syncthreads();
    }

    if (tid == 0) atomicMax(result, sdata[0]);
}

// ─────────────────────────── GPU buffer manager ───────────────────────────

// Owns all device-side memory for one solver instance.
// Geometry arrays (d_nb, d_diag, d_b, d_coeff) are uploaded once per
// buildNodes() call.  Work vectors are reused across BiCGSTAB iterations
// within the same coupling sub-step.
struct GpuBiCGSTABBuffers {
    // Geometry (uploaded once when the node grid is rebuilt)
    float*        d_diag  = nullptr;
    float*        d_b     = nullptr;
    uint32_t*     d_nb    = nullptr;
    float*        d_coeff = nullptr;

    // BiCGSTAB work vectors (one float per node)
    float* d_x    = nullptr;
    float* d_r    = nullptr;
    float* d_rhat = nullptr;
    float* d_p    = nullptr;
    float* d_v    = nullptr;
    float* d_y    = nullptr;
    float* d_z    = nullptr;
    float* d_s    = nullptr;
    float* d_t    = nullptr;
    float* d_Ax   = nullptr;

    // Scalar accumulators for reductions (one element each)
    float*        d_dot  = nullptr;
    unsigned int* d_maxu = nullptr;

    uint32_t n      = 0;
    int      nFaces = 0;
    bool     valid  = false;

    // Non-copyable
    GpuBiCGSTABBuffers() = default;
    GpuBiCGSTABBuffers(const GpuBiCGSTABBuffers&) = delete;
    GpuBiCGSTABBuffers& operator=(const GpuBiCGSTABBuffers&) = delete;

    ~GpuBiCGSTABBuffers() { free(); }

    void free() {
        auto f = [](auto*& p) { if (p) { cudaFree(p); p = nullptr; } };
        f(d_diag); f(d_b); f(d_nb); f(d_coeff);
        f(d_x); f(d_r); f(d_rhat); f(d_p); f(d_v);
        f(d_y); f(d_z); f(d_s); f(d_t); f(d_Ax);
        f(d_dot); f(d_maxu);
        n = 0; nFaces = 0; valid = false;
    }

    void allocate(uint32_t nodeCount, int faceCount) {
        free();
        n      = nodeCount;
        nFaces = faceCount;
        const size_t N  = n;
        const size_t NF = N * static_cast<size_t>(nFaces);
        LS_CUDA_CHECK(cudaMalloc(&d_diag,  N  * sizeof(float)));
        LS_CUDA_CHECK(cudaMalloc(&d_b,     N  * sizeof(float)));
        LS_CUDA_CHECK(cudaMalloc(&d_nb,    NF * sizeof(uint32_t)));
        LS_CUDA_CHECK(cudaMalloc(&d_coeff, NF * sizeof(float)));
        LS_CUDA_CHECK(cudaMalloc(&d_x,     N  * sizeof(float)));
        LS_CUDA_CHECK(cudaMalloc(&d_r,     N  * sizeof(float)));
        LS_CUDA_CHECK(cudaMalloc(&d_rhat,  N  * sizeof(float)));
        LS_CUDA_CHECK(cudaMalloc(&d_p,     N  * sizeof(float)));
        LS_CUDA_CHECK(cudaMalloc(&d_v,     N  * sizeof(float)));
        LS_CUDA_CHECK(cudaMalloc(&d_y,     N  * sizeof(float)));
        LS_CUDA_CHECK(cudaMalloc(&d_z,     N  * sizeof(float)));
        LS_CUDA_CHECK(cudaMalloc(&d_s,     N  * sizeof(float)));
        LS_CUDA_CHECK(cudaMalloc(&d_t,     N  * sizeof(float)));
        LS_CUDA_CHECK(cudaMalloc(&d_Ax,    N  * sizeof(float)));
        LS_CUDA_CHECK(cudaMalloc(&d_dot,   sizeof(float)));
        LS_CUDA_CHECK(cudaMalloc(&d_maxu,  sizeof(unsigned int)));
        valid = true;
    }

    // ── primitive reductions ───────────────────────────────────────────

    // Dot product: result = sum(a[i] * b[i])
    float dotProduct(const float* a, const float* b) const {
        LS_CUDA_CHECK(cudaMemset(d_dot, 0, sizeof(float)));
        const int blk = std::min(blocksFor(static_cast<int>(n)), 1024);
        dotKernel<<<blk, kThreadsPerBlock,
                    kThreadsPerBlock * sizeof(float)>>>(a, b, d_dot, n);
        float result = 0.f;
        LS_CUDA_CHECK(cudaMemcpy(&result, d_dot, sizeof(float),
                                 cudaMemcpyDeviceToHost));
        return result;
    }

    // Max absolute value: result = max(|v[i]|)
    float maxAbs(const float* v) const {
        LS_CUDA_CHECK(cudaMemset(d_maxu, 0, sizeof(unsigned int)));
        const int blk = std::min(blocksFor(static_cast<int>(n)), 1024);
        maxAbsKernel<<<blk, kThreadsPerBlock,
                       kThreadsPerBlock * sizeof(unsigned int)>>>(
            v, d_maxu, n);
        unsigned int bits = 0u;
        LS_CUDA_CHECK(cudaMemcpy(&bits, d_maxu, sizeof(unsigned int),
                                 cudaMemcpyDeviceToHost));
        // Reverse the IEEE-754 bit-cast trick used in the kernel
        float result;
        __builtin_memcpy(&result, &bits, sizeof(float));
        return result;
    }

    // ── SpMV dispatch ──────────────────────────────────────────────────

    void spmv(const float* v_in, float* Av_out) const {
        const int blk = std::min(blocksFor(static_cast<int>(n)), 65535);
        if (nFaces == 4)
            spmvKernel<4><<<blk, kThreadsPerBlock>>>(
                v_in, d_diag, d_nb, d_coeff, Av_out, n);
        else
            spmvKernel<6><<<blk, kThreadsPerBlock>>>(
                v_in, d_diag, d_nb, d_coeff, Av_out, n);
    }
};

// ────────────────────── BiCGSTAB orchestration ───────────────────────────

// Run GPU-accelerated preconditioned BiCGSTAB.
//
//   x (host, in/out) : initial guess on entry, solution on exit.
//   diag_gpu must already hold d_diag and d_b, and d_nb / d_coeff must
//   have been uploaded to gpu before this call.
//
// Scalars (rho, alpha, omega) stay on the CPU for numerical stability;
// only dot-product scalars and the max-abs residual are downloaded.
inline void solveBiCGSTAB(
    GpuBiCGSTABBuffers&    gpu,
    std::vector<float>&    x,           // size n, host
    float                  diagEps,
    unsigned               maxIter,
    float                  tolerance,
    unsigned&              outIterations,
    double&                outResidual)
{
    assert(gpu.valid);
    const uint32_t n   = gpu.n;
    const int      blk = std::min(blocksFor(static_cast<int>(n)), 65535);

    // Upload initial guess to device
    LS_CUDA_CHECK(cudaMemcpy(gpu.d_x, x.data(), n * sizeof(float),
                             cudaMemcpyHostToDevice));
    // Zero all work vectors.  cudaMalloc does not zero GPU memory, and d_z in
    // particular is read (multiplied by omega=0) on the early-convergence path
    // before it has been written by jacobiKernel — giving 0*NaN = NaN if the
    // GPU heap happened to contain NaN from a prior allocation.
    LS_CUDA_CHECK(cudaMemset(gpu.d_p,   0, n * sizeof(float)));
    LS_CUDA_CHECK(cudaMemset(gpu.d_v,   0, n * sizeof(float)));
    LS_CUDA_CHECK(cudaMemset(gpu.d_y,   0, n * sizeof(float)));
    LS_CUDA_CHECK(cudaMemset(gpu.d_z,   0, n * sizeof(float)));
    LS_CUDA_CHECK(cudaMemset(gpu.d_s,   0, n * sizeof(float)));
    LS_CUDA_CHECK(cudaMemset(gpu.d_t,   0, n * sizeof(float)));
    LS_CUDA_CHECK(cudaMemset(gpu.d_Ax,  0, n * sizeof(float)));

    // r0 = b - A*x0
    gpu.spmv(gpu.d_x, gpu.d_Ax);
    residualKernel<<<blk, kThreadsPerBlock>>>(gpu.d_b, gpu.d_Ax, gpu.d_r, n);
    // r_hat = r0 (shadow residual; never changes)
    LS_CUDA_CHECK(cudaMemcpy(gpu.d_rhat, gpu.d_r, n * sizeof(float),
                             cudaMemcpyDeviceToDevice));

    // Normalise both convergence checks by ||b||_inf, matching the CPU path which
    // uses `residual < tolerance * b_norm`.  Without this the GPU checks the
    // raw (unnormalised) residual, which is 10–1000× tighter than the CPU
    // threshold.  BiCGSTAB then keeps iterating on an already-converged system
    // until breakdown (ρ → 0 → α = NaN) poisons the solution.
    // Both the s-check (early exit) and r-check (end-of-loop) use eff_tol.
    const float b_norm = gpu.maxAbs(gpu.d_b);
    const float eff_tol = (b_norm > 1e-37f) ? tolerance * b_norm : tolerance;

    double rho   = 1.0;
    double alpha = 1.0;
    double omega = 1.0;
    outIterations = 0;
    outResidual   = 0.0;

    for (unsigned iter = 0; iter < maxIter; ++iter) {

        // rho_new = <r_hat, r>
        const double rho_new =
            static_cast<double>(gpu.dotProduct(gpu.d_rhat, gpu.d_r));
        if (std::abs(rho_new) < 1e-100) break;

        // p = r + beta*(p - omega*v),  beta = (rho_new/rho)*(alpha/omega)
        const float beta       = static_cast<float>((rho_new / rho) * (alpha / omega));
        const float beta_omega = static_cast<float>(beta * omega);
        rho = rho_new;
        updatePKernel<<<blk, kThreadsPerBlock>>>(
            gpu.d_r, beta, beta_omega, gpu.d_v, gpu.d_p, n);

        // y = M^{-1} p  (Jacobi preconditioner)
        jacobiKernel<<<blk, kThreadsPerBlock>>>(
            gpu.d_p, gpu.d_diag, gpu.d_y, diagEps, n);

        // v = A y
        gpu.spmv(gpu.d_y, gpu.d_v);

        // alpha = rho_new / <r_hat, v>
        const double r_hat_v =
            static_cast<double>(gpu.dotProduct(gpu.d_rhat, gpu.d_v));
        if (std::abs(r_hat_v) < 1e-100) break;
        alpha = rho_new / r_hat_v;

        // s = r - alpha*v
        subtractScaledKernel<<<blk, kThreadsPerBlock>>>(
            gpu.d_r, static_cast<float>(alpha), gpu.d_v, gpu.d_s, n);

        // Check convergence on s
        outResidual = static_cast<double>(gpu.maxAbs(gpu.d_s));
        if (outResidual < eff_tol) {
            // x += alpha*y  (z term is zero since we converged early)
            updateXKernel<<<blk, kThreadsPerBlock>>>(
                gpu.d_x, static_cast<float>(alpha), gpu.d_y,
                0.f, gpu.d_z, n);
            ++outIterations;
            break;
        }

        // z = M^{-1} s  (Jacobi)
        jacobiKernel<<<blk, kThreadsPerBlock>>>(
            gpu.d_s, gpu.d_diag, gpu.d_z, diagEps, n);

        // t = A z
        gpu.spmv(gpu.d_z, gpu.d_t);

        // omega = <t,s> / <t,t>
        const double t_s = static_cast<double>(gpu.dotProduct(gpu.d_t, gpu.d_s));
        const double t_t = static_cast<double>(gpu.dotProduct(gpu.d_t, gpu.d_t));
        omega = (t_t > 1e-100) ? t_s / t_t : 0.0;

        // x += alpha*y + omega*z
        updateXKernel<<<blk, kThreadsPerBlock>>>(
            gpu.d_x, static_cast<float>(alpha), gpu.d_y,
            static_cast<float>(omega), gpu.d_z, n);

        // r = s - omega*t
        updateRKernel<<<blk, kThreadsPerBlock>>>(
            gpu.d_s, static_cast<float>(omega), gpu.d_t, gpu.d_r, n);

        // Check convergence on r
        outResidual = static_cast<double>(gpu.maxAbs(gpu.d_r));
        ++outIterations;
        if (outResidual < eff_tol) break;
    }

    // Download solution to host
    LS_CUDA_CHECK(cudaMemcpy(x.data(), gpu.d_x, n * sizeof(float),
                             cudaMemcpyDeviceToHost));
}

} // namespace gpu
} // namespace viennals

#endif // VIENNALS_GPU_BICGSTAB
