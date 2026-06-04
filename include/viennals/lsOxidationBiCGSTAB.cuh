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
#include <cusparse.h>
#include <algorithm>
#include <cassert>
#include <cmath>
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

// ─────────────────────── CUSPARSE error check ────────────────────────────

#ifdef NDEBUG
#define LS_CUSPARSE_CHECK(call) (call)
#else
#define LS_CUSPARSE_CHECK(call)                                                \
    do {                                                                       \
        cusparseStatus_t _s = (call);                                          \
        if (_s != CUSPARSE_STATUS_SUCCESS)                                     \
            fprintf(stderr, "CUSPARSE error %s:%d – %d\n", __FILE__,          \
                    __LINE__, static_cast<int>(_s));                           \
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
    const double* __restrict__ v,
    const double* __restrict__ diag,
    const uint32_t* __restrict__ nb,
    const double* __restrict__ coeff,
    double* __restrict__ Av,
    uint32_t n)
{
    const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    double res = diag[i] * v[i];
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
    const double* __restrict__ b,
    const double* __restrict__ Av,
    double* __restrict__ r,
    uint32_t n)
{
    const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    r[i] = b[i] - Av[i];
}

// p[i] = r[i] + beta*(p[i] - omega*v[i])
__global__ void updatePKernel(
    const double* __restrict__ r,
    double beta, double beta_omega,
    const double* __restrict__ v,
    double* __restrict__ p,
    uint32_t n)
{
    const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    p[i] = r[i] + beta * p[i] - beta_omega * v[i];
}

// out[i] = in[i] - alpha*v[i]  (used for s = r - alpha*v)
__global__ void subtractScaledKernel(
    const double* __restrict__ in,
    double alpha,
    const double* __restrict__ v,
    double* __restrict__ out,
    uint32_t n)
{
    const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = in[i] - alpha * v[i];
}

// Jacobi preconditioner: out[i] = in[i] / diag[i]  (guarded against zero diagonal)
__global__ void jacobiKernel(
    const double* __restrict__ in,
    const double* __restrict__ diag,
    double* __restrict__ out,
    double eps,
    uint32_t n)
{
    const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    out[i] = (diag[i] > eps) ? in[i] / diag[i] : in[i];
}

// x[i] += alpha*y[i] + omega*z[i]
__global__ void updateXKernel(
    double* __restrict__ x,
    double alpha, const double* __restrict__ y,
    double omega, const double* __restrict__ z,
    uint32_t n)
{
    const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    x[i] += alpha * y[i] + omega * z[i];
}

// r[i] = s[i] - omega*t[i]  (residual update at end of BiCGSTAB iteration)
__global__ void updateRKernel(
    const double* __restrict__ s,
    double omega,
    const double* __restrict__ t,
    double* __restrict__ r,
    uint32_t n)
{
    const uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    r[i] = s[i] - omega * t[i];
}

// Atomic add for double. Native double atomicAdd requires sm_60; the project
// may be compiled for older compatibility architectures, so provide the
// standard CAS fallback.
__device__ inline double atomicAddDouble(double* address, double val) {
#if __CUDA_ARCH__ >= 600
    return atomicAdd(address, val);
#else
    unsigned long long int* address_as_ull =
        reinterpret_cast<unsigned long long int*>(address);
    unsigned long long int old = *address_as_ull;
    unsigned long long int assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(
                            val + __longlong_as_double(assumed)));
    } while (assumed != old);

    return __longlong_as_double(old);
#endif
}

// Dot product: atomicAdd accumulates partial sums into result (must be zeroed before launch).
// Uses stride-loop so a fixed grid size (≤ 1024 blocks) handles any n.
__global__ void dotKernel(
    const double* __restrict__ a,
    const double* __restrict__ b,
    double* __restrict__ result,
    uint32_t n)
{
    extern __shared__ double sdata[];

    const uint32_t tid = threadIdx.x;
    double sum = 0.0;
    for (uint32_t i = blockIdx.x * blockDim.x + tid;
         i < n;
         i += gridDim.x * blockDim.x) {
        sum += a[i] * b[i];
    }

    sdata[tid] = sum;
    __syncthreads();

    for (unsigned s = blockDim.x >> 1; s > 0; s >>= 1) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        __syncthreads();
    }

    if (tid == 0) atomicAddDouble(result, sdata[0]);
}

// Max absolute value using IEEE-754 bit trick:
// for non-negative doubles, bit-pattern ordering equals double ordering.
// result (unsigned long long*) must be zeroed before launch.
__global__ void maxAbsKernel(
    const double* __restrict__ v,
    unsigned long long int* __restrict__ result,
    uint32_t n)
{
    extern __shared__ unsigned char smem[];
    unsigned long long int* sdata =
        reinterpret_cast<unsigned long long int*>(smem);

    const uint32_t tid = threadIdx.x;
    unsigned long long int maxVal = 0ull;
    for (uint32_t i = blockIdx.x * blockDim.x + tid;
         i < n;
         i += gridDim.x * blockDim.x) {
        unsigned long long int bits =
            static_cast<unsigned long long int>(__double_as_longlong(fabs(v[i])));
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
    double*       d_diag  = nullptr;
    double*       d_b     = nullptr;
    uint32_t*     d_nb    = nullptr;
    double*       d_coeff = nullptr;

    // BiCGSTAB work vectors (one double per node)
    double* d_x    = nullptr;
    double* d_r    = nullptr;
    double* d_rhat = nullptr;
    double* d_p    = nullptr;
    double* d_v    = nullptr;
    double* d_y    = nullptr;
    double* d_z    = nullptr;
    double* d_s    = nullptr;
    double* d_t    = nullptr;
    double* d_Ax   = nullptr;

    // Scalar accumulators for reductions (one element each)
    double*                 d_dot  = nullptr;
    unsigned long long int* d_maxu = nullptr;

    // ── ILU(0) preconditioner state ───────────────────────────────────────
    // CUSPARSE handles/descriptors (created once per buffer lifetime).
    cusparseHandle_t  csHandle = nullptr;
    // Matrix descriptors: A (general, for ILU factorization), L, U (for SpSV).
    cusparseMatDescr_t   descrA   = nullptr;  // legacy descriptor; still used by csrilu02
    cusparseSpMatDescr_t spMatL   = nullptr;  // generic CSR descriptor for L factor
    cusparseSpMatDescr_t spMatU   = nullptr;  // generic CSR descriptor for U factor
    // ILU(0) factorization info and triangular-solve descriptors.
    // csrilu02Info_t is deprecated in CUDA 12 but is the only supported
    // in-place ILU(0) factorization interface; use it with the warning suppressed.
    csrilu02Info_t      iluInfo  = nullptr;
    cusparseSpSVDescr_t svDescrL = nullptr;  // SpSV solve descriptor for L
    cusparseSpSVDescr_t svDescrU = nullptr;  // SpSV solve descriptor for U
    // CSR arrays for the ILU(0) factor (pattern set in setupCSR; values refreshed each solve).
    int*    d_csrRowPtr = nullptr;  // n+1 entries
    int*    d_csrColInd = nullptr;  // nnz entries
    double* d_iluVal    = nullptr;  // nnz entries — overwritten in-place by ILU(0)
    // Intermediate vector for the two-step triangular solve: L*svTmp = rhs, U*out = svTmp.
    double* d_svTmp     = nullptr;  // n entries
    // CUSPARSE internal work buffers (sizes determined at analysis time).
    void*   d_iluBuf    = nullptr;
    void*   d_svBufL    = nullptr;
    void*   d_svBufU    = nullptr;
    // Number of non-zeros in the ILU pattern.
    int     nnz         = 0;
    // True when setupCSR and analysis have completed successfully.
    bool    iluReady    = false;
    // Host-side CSR index arrays (kept for fillAndFactorizeILU value assembly).
    std::vector<int> h_csrRowPtr;
    std::vector<int> h_csrColInd;
    // h_csrFaceIdx[k] = face index fi for off-diagonal entry k, or -1 for diagonal.
    std::vector<int> h_csrFaceIdx;

    uint32_t n      = 0;
    int      nFaces = 0;
    bool     valid  = false;

    // Non-copyable
    GpuBiCGSTABBuffers() = default;
    GpuBiCGSTABBuffers(const GpuBiCGSTABBuffers&) = delete;
    GpuBiCGSTABBuffers& operator=(const GpuBiCGSTABBuffers&) = delete;

    ~GpuBiCGSTABBuffers() { free(); }

    void free() {
        // Destroy CUSPARSE objects before freeing device memory.
        iluReady = false;
        if (svDescrU)  { cusparseSpSV_destroyDescr(svDescrU);   svDescrU = nullptr; }
        if (svDescrL)  { cusparseSpSV_destroyDescr(svDescrL);   svDescrL = nullptr; }
        if (spMatU)    { cusparseDestroySpMat(spMatU);           spMatU   = nullptr; }
        if (spMatL)    { cusparseDestroySpMat(spMatL);           spMatL   = nullptr; }
        if (iluInfo)   { cusparseDestroyCsrilu02Info(iluInfo);   iluInfo  = nullptr; }
        if (descrA)    { cusparseDestroyMatDescr(descrA);        descrA   = nullptr; }
        if (csHandle)  { cusparseDestroy(csHandle);              csHandle = nullptr; }

        auto f = [](auto*& p) { if (p) { cudaFree(p); p = nullptr; } };
        f(d_iluBuf); f(d_svBufL); f(d_svBufU);
        f(d_svTmp); f(d_iluVal); f(d_csrColInd); f(d_csrRowPtr);
        f(d_diag); f(d_b); f(d_nb); f(d_coeff);
        f(d_x); f(d_r); f(d_rhat); f(d_p); f(d_v);
        f(d_y); f(d_z); f(d_s); f(d_t); f(d_Ax);
        f(d_dot); f(d_maxu);
        nnz = 0;
        h_csrRowPtr.clear(); h_csrColInd.clear(); h_csrFaceIdx.clear();
        n = 0; nFaces = 0; valid = false;
    }

    void allocate(uint32_t nodeCount, int faceCount) {
        free();
        n      = nodeCount;
        nFaces = faceCount;
        const size_t N  = n;
        const size_t NF = N * static_cast<size_t>(nFaces);
        // Every cudaMalloc is checked; any failure marks the buffers invalid.
        // This prevents the GPU path from being entered with null device pointers
        // when the CUDA context cannot be initialized (e.g. driver mismatch).
#define TRY_MALLOC(ptr, sz) \
        if (cudaMalloc(&(ptr), (sz)) != cudaSuccess) { cudaGetLastError(); return; }
        TRY_MALLOC(d_diag,  N  * sizeof(double))
        TRY_MALLOC(d_b,     N  * sizeof(double))
        TRY_MALLOC(d_nb,    NF * sizeof(uint32_t))
        TRY_MALLOC(d_coeff, NF * sizeof(double))
        TRY_MALLOC(d_x,     N  * sizeof(double))
        TRY_MALLOC(d_r,     N  * sizeof(double))
        TRY_MALLOC(d_rhat,  N  * sizeof(double))
        TRY_MALLOC(d_p,     N  * sizeof(double))
        TRY_MALLOC(d_v,     N  * sizeof(double))
        TRY_MALLOC(d_y,     N  * sizeof(double))
        TRY_MALLOC(d_z,     N  * sizeof(double))
        TRY_MALLOC(d_s,     N  * sizeof(double))
        TRY_MALLOC(d_t,     N  * sizeof(double))
        TRY_MALLOC(d_Ax,    N  * sizeof(double))
        TRY_MALLOC(d_dot,   sizeof(double))
        TRY_MALLOC(d_maxu,  sizeof(unsigned long long int))
#undef TRY_MALLOC
        valid = true;
    }

    // ── primitive reductions ───────────────────────────────────────────

    // Dot product: result = sum(a[i] * b[i]).
    double dotProduct(const double* a, const double* b) const {
        LS_CUDA_CHECK(cudaMemset(d_dot, 0, sizeof(double)));
        const int blk = std::min(blocksFor(static_cast<int>(n)), 1024);
        dotKernel<<<blk, kThreadsPerBlock,
                    kThreadsPerBlock * sizeof(double)>>>(a, b, d_dot, n);
        double result = 0.0;
        LS_CUDA_CHECK(cudaMemcpy(&result, d_dot, sizeof(double),
                                 cudaMemcpyDeviceToHost));
        return result;
    }

    // Max absolute value: result = max(|v[i]|)
    double maxAbs(const double* v) const {
        LS_CUDA_CHECK(cudaMemset(d_maxu, 0, sizeof(unsigned long long int)));
        const int blk = std::min(blocksFor(static_cast<int>(n)), 1024);
        maxAbsKernel<<<blk, kThreadsPerBlock,
                       kThreadsPerBlock * sizeof(unsigned long long int)>>>(
            v, d_maxu, n);
        unsigned long long int bits = 0ull;
        LS_CUDA_CHECK(cudaMemcpy(&bits, d_maxu, sizeof(unsigned long long int),
                                 cudaMemcpyDeviceToHost));
        // Reverse the IEEE-754 bit-cast trick used in the kernel
        double result;
        __builtin_memcpy(&result, &bits, sizeof(double));
        return result;
    }

    // ── SpMV dispatch ──────────────────────────────────────────────────

    void spmv(const double* v_in, double* Av_out) const {
        const int blk = std::min(blocksFor(static_cast<int>(n)), 65535);
        if (nFaces == 4)
            spmvKernel<4><<<blk, kThreadsPerBlock>>>(
                v_in, d_diag, d_nb, d_coeff, Av_out, n);
        else
            spmvKernel<6><<<blk, kThreadsPerBlock>>>(
                v_in, d_diag, d_nb, d_coeff, Av_out, n);
    }

    // ── ILU(0) setup ───────────────────────────────────────────────────

    // Build the CSR sparsity pattern from the face-major neighbor array h_nb,
    // upload index arrays to the device, create CUSPARSE handles/descriptors,
    // and run symbolic analysis for both the ILU(0) factorization and the
    // two triangular solves.  Called once per buildNodes() cycle.
    // h_nb is the nb32 array in face-major layout: h_nb[fi * nodeCount + i].
    void setupCSR(const uint32_t* h_nb, uint32_t nodeCount, int faceCount) {
        // Destroy any previous ILU state but keep the main CUDA allocations.
        iluReady = false;
        if (svDescrU)  { cusparseSpSV_destroyDescr(svDescrU);   svDescrU = nullptr; }
        if (svDescrL)  { cusparseSpSV_destroyDescr(svDescrL);   svDescrL = nullptr; }
        if (spMatU)    { cusparseDestroySpMat(spMatU);           spMatU   = nullptr; }
        if (spMatL)    { cusparseDestroySpMat(spMatL);           spMatL   = nullptr; }
        if (iluInfo)   { cusparseDestroyCsrilu02Info(iluInfo);   iluInfo  = nullptr; }
        if (descrA)    { cusparseDestroyMatDescr(descrA);        descrA   = nullptr; }
        if (csHandle)  { cusparseDestroy(csHandle);              csHandle = nullptr; }
        if (d_iluBuf)  { cudaFree(d_iluBuf);  d_iluBuf  = nullptr; }
        if (d_svBufL)  { cudaFree(d_svBufL);  d_svBufL  = nullptr; }
        if (d_svBufU)  { cudaFree(d_svBufU);  d_svBufU  = nullptr; }
        if (d_svTmp)   { cudaFree(d_svTmp);   d_svTmp   = nullptr; }
        if (d_iluVal)  { cudaFree(d_iluVal);  d_iluVal  = nullptr; }
        if (d_csrColInd) { cudaFree(d_csrColInd); d_csrColInd = nullptr; }
        if (d_csrRowPtr) { cudaFree(d_csrRowPtr); d_csrRowPtr = nullptr; }
        h_csrRowPtr.clear(); h_csrColInd.clear(); h_csrFaceIdx.clear();
        nnz = 0;

        // ── Build sorted CSR on CPU ────────────────────────────────────
        // For each row i, collect: the diagonal entry (col = i, faceIdx = -1)
        // and all valid off-diagonal entries (col = h_nb[fi*n+i], faceIdx = fi).
        // Then sort by column index.
        const uint32_t N = nodeCount;
        h_csrRowPtr.resize(static_cast<std::size_t>(N) + 1, 0);

        // First pass: count entries per row.
        for (uint32_t i = 0; i < N; ++i) {
            int cnt = 1; // diagonal always present
            for (int fi = 0; fi < faceCount; ++fi) {
                const uint32_t j = h_nb[static_cast<uint32_t>(fi) * N + i];
                if (j != 0xFFFFFFFFu)
                    ++cnt;
            }
            h_csrRowPtr[i + 1] = cnt;
        }
        // Prefix-sum into row pointers.
        for (uint32_t i = 0; i < N; ++i)
            h_csrRowPtr[i + 1] += h_csrRowPtr[i];

        const int totalNnz = h_csrRowPtr[N];
        h_csrColInd.resize(static_cast<std::size_t>(totalNnz));
        h_csrFaceIdx.resize(static_cast<std::size_t>(totalNnz));

        // Second pass: fill column indices and face indices.
        // Use a temporary write-pointer per row.
        std::vector<int> writePtr(h_csrRowPtr.begin(), h_csrRowPtr.begin() + N);
        for (uint32_t i = 0; i < N; ++i) {
            const int base = writePtr[i];
            int k = base;
            // Diagonal entry first (will sort to correct position).
            h_csrColInd[k]  = static_cast<int>(i);
            h_csrFaceIdx[k] = -1;
            ++k;
            // Off-diagonal entries.
            for (int fi = 0; fi < faceCount; ++fi) {
                const uint32_t j = h_nb[static_cast<uint32_t>(fi) * N + i];
                if (j != 0xFFFFFFFFu) {
                    h_csrColInd[k]  = static_cast<int>(j);
                    h_csrFaceIdx[k] = fi;
                    ++k;
                }
            }
            const int rowEnd = h_csrRowPtr[i + 1];
            const int rowLen = rowEnd - base;
            // Sort the row's entries by column index (insertion sort — rows are small).
            for (int a = 1; a < rowLen; ++a) {
                const int colA  = h_csrColInd[base + a];
                const int faceA = h_csrFaceIdx[base + a];
                int b = a - 1;
                while (b >= 0 && h_csrColInd[base + b] > colA) {
                    h_csrColInd[base + b + 1]  = h_csrColInd[base + b];
                    h_csrFaceIdx[base + b + 1] = h_csrFaceIdx[base + b];
                    --b;
                }
                h_csrColInd[base + b + 1]  = colA;
                h_csrFaceIdx[base + b + 1] = faceA;
            }
        }
        nnz = totalNnz;

        // ── Upload CSR index arrays to device ──────────────────────────
        LS_CUDA_CHECK(cudaMalloc(&d_csrRowPtr,
                                 (static_cast<std::size_t>(N) + 1) * sizeof(int)));
        LS_CUDA_CHECK(cudaMalloc(&d_csrColInd,
                                 static_cast<std::size_t>(nnz) * sizeof(int)));
        LS_CUDA_CHECK(cudaMalloc(&d_iluVal,
                                 static_cast<std::size_t>(nnz) * sizeof(double)));
        LS_CUDA_CHECK(cudaMalloc(&d_svTmp,
                                 static_cast<std::size_t>(N) * sizeof(double)));

        LS_CUDA_CHECK(cudaMemcpy(d_csrRowPtr, h_csrRowPtr.data(),
                                 (static_cast<std::size_t>(N) + 1) * sizeof(int),
                                 cudaMemcpyHostToDevice));
        LS_CUDA_CHECK(cudaMemcpy(d_csrColInd, h_csrColInd.data(),
                                 static_cast<std::size_t>(nnz) * sizeof(int),
                                 cudaMemcpyHostToDevice));

        // Fill placeholder values so that analysis can proceed:
        // diagonal = 1.0, off-diagonal = -1.0 / faceCount.
        {
            const double offVal = -1.0 / static_cast<double>(faceCount > 0 ? faceCount : 1);
            std::vector<double> hPlaceholder(static_cast<std::size_t>(nnz));
            for (uint32_t i = 0; i < N; ++i) {
                for (int k = h_csrRowPtr[i]; k < h_csrRowPtr[i + 1]; ++k) {
                    hPlaceholder[k] = (h_csrColInd[k] == static_cast<int>(i))
                                      ? 1.0 : offVal;
                }
            }
            LS_CUDA_CHECK(cudaMemcpy(d_iluVal, hPlaceholder.data(),
                                     static_cast<std::size_t>(nnz) * sizeof(double),
                                     cudaMemcpyHostToDevice));
        }

        // ── Create CUSPARSE handle and descriptors ────────────────────
        if (cusparseCreate(&csHandle) != CUSPARSE_STATUS_SUCCESS) {
            fprintf(stderr, "GpuBiCGSTABBuffers: cusparseCreate failed\n");
            return;
        }

        if (cusparseCreateMatDescr(&descrA) != CUSPARSE_STATUS_SUCCESS) return;
        cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);

        // ── Create ILU(0) info (legacy API, still present in CUDA 12) ────
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
        if (cusparseCreateCsrilu02Info(&iluInfo) != CUSPARSE_STATUS_SUCCESS) return;
#pragma GCC diagnostic pop

        // ── Create generic SpMat descriptors for L and U (CUDA 12 SpSV API) ─
        // Both L and U share d_iluVal (in-place ILU factorization).
        cusparseFillMode_t fillL = CUSPARSE_FILL_MODE_LOWER;
        cusparseDiagType_t diagL = CUSPARSE_DIAG_TYPE_UNIT;
        cusparseFillMode_t fillU = CUSPARSE_FILL_MODE_UPPER;
        cusparseDiagType_t diagU = CUSPARSE_DIAG_TYPE_NON_UNIT;

        if (cusparseCreateCsr(&spMatL, static_cast<int64_t>(N), static_cast<int64_t>(N),
                              static_cast<int64_t>(nnz),
                              d_csrRowPtr, d_csrColInd, d_iluVal,
                              CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                              CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F) != CUSPARSE_STATUS_SUCCESS) return;
        cusparseSpMatSetAttribute(spMatL, CUSPARSE_SPMAT_FILL_MODE,
                                  static_cast<void*>(&fillL), sizeof(fillL));
        cusparseSpMatSetAttribute(spMatL, CUSPARSE_SPMAT_DIAG_TYPE,
                                  static_cast<void*>(&diagL), sizeof(diagL));

        if (cusparseCreateCsr(&spMatU, static_cast<int64_t>(N), static_cast<int64_t>(N),
                              static_cast<int64_t>(nnz),
                              d_csrRowPtr, d_csrColInd, d_iluVal,
                              CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                              CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F) != CUSPARSE_STATUS_SUCCESS) return;
        cusparseSpMatSetAttribute(spMatU, CUSPARSE_SPMAT_FILL_MODE,
                                  static_cast<void*>(&fillU), sizeof(fillU));
        cusparseSpMatSetAttribute(spMatU, CUSPARSE_SPMAT_DIAG_TYPE,
                                  static_cast<void*>(&diagU), sizeof(diagU));

        if (cusparseSpSV_createDescr(&svDescrL) != CUSPARSE_STATUS_SUCCESS) return;
        if (cusparseSpSV_createDescr(&svDescrU) != CUSPARSE_STATUS_SUCCESS) return;

        const int iN = static_cast<int>(N);

        // ── ILU(0) analysis ───────────────────────────────────────────
        int iluBufSize = 0;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
        cusparseStatus_t st = cusparseDcsrilu02_bufferSize(
            csHandle, iN, nnz, descrA,
            d_iluVal, d_csrRowPtr, d_csrColInd,
            iluInfo, &iluBufSize);
#pragma GCC diagnostic pop
        if (st != CUSPARSE_STATUS_SUCCESS) {
            fprintf(stderr, "GpuBiCGSTABBuffers: cusparseDcsrilu02_bufferSize failed\n");
            return;
        }
        LS_CUDA_CHECK(cudaMalloc(&d_iluBuf, static_cast<std::size_t>(iluBufSize)));

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
        st = cusparseDcsrilu02_analysis(
            csHandle, iN, nnz, descrA,
            d_iluVal, d_csrRowPtr, d_csrColInd,
            iluInfo, CUSPARSE_SOLVE_POLICY_USE_LEVEL, d_iluBuf);
#pragma GCC diagnostic pop
        if (st != CUSPARSE_STATUS_SUCCESS) {
            fprintf(stderr, "GpuBiCGSTABBuffers: cusparseDcsrilu02_analysis failed\n");
            return;
        }

        // ── Triangular-solve analyses (L and U, done once per pattern) ─
        // Temporary dense-vector descriptors pointing to d_svTmp (analysis only).
        cusparseDnVecDescr_t vecTmp1 = nullptr, vecTmp2 = nullptr;
        cusparseCreateDnVec(&vecTmp1, static_cast<int64_t>(N), d_svTmp, CUDA_R_64F);
        cusparseCreateDnVec(&vecTmp2, static_cast<int64_t>(N), d_svTmp, CUDA_R_64F);

        const double one = 1.0;
        size_t svBufSizeL = 0, svBufSizeU = 0;

        st = cusparseSpSV_bufferSize(
            csHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
            spMatL, vecTmp1, vecTmp2, CUDA_R_64F,
            CUSPARSE_SPSV_ALG_DEFAULT, svDescrL, &svBufSizeL);
        if (st != CUSPARSE_STATUS_SUCCESS) {
            fprintf(stderr, "GpuBiCGSTABBuffers: cusparseSpSV_bufferSize (L) failed\n");
            cusparseDestroyDnVec(vecTmp1); cusparseDestroyDnVec(vecTmp2);
            return;
        }

        st = cusparseSpSV_bufferSize(
            csHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
            spMatU, vecTmp1, vecTmp2, CUDA_R_64F,
            CUSPARSE_SPSV_ALG_DEFAULT, svDescrU, &svBufSizeU);
        if (st != CUSPARSE_STATUS_SUCCESS) {
            fprintf(stderr, "GpuBiCGSTABBuffers: cusparseSpSV_bufferSize (U) failed\n");
            cusparseDestroyDnVec(vecTmp1); cusparseDestroyDnVec(vecTmp2);
            return;
        }

        LS_CUDA_CHECK(cudaMalloc(&d_svBufL, svBufSizeL ? svBufSizeL : 1));
        LS_CUDA_CHECK(cudaMalloc(&d_svBufU, svBufSizeU ? svBufSizeU : 1));

        st = cusparseSpSV_analysis(
            csHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
            spMatL, vecTmp1, vecTmp2, CUDA_R_64F,
            CUSPARSE_SPSV_ALG_DEFAULT, svDescrL, d_svBufL);
        if (st != CUSPARSE_STATUS_SUCCESS) {
            fprintf(stderr, "GpuBiCGSTABBuffers: cusparseSpSV_analysis (L) failed\n");
            cusparseDestroyDnVec(vecTmp1); cusparseDestroyDnVec(vecTmp2);
            return;
        }

        st = cusparseSpSV_analysis(
            csHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
            spMatU, vecTmp1, vecTmp2, CUDA_R_64F,
            CUSPARSE_SPSV_ALG_DEFAULT, svDescrU, d_svBufU);
        if (st != CUSPARSE_STATUS_SUCCESS) {
            fprintf(stderr, "GpuBiCGSTABBuffers: cusparseSpSV_analysis (U) failed\n");
            cusparseDestroyDnVec(vecTmp1); cusparseDestroyDnVec(vecTmp2);
            return;
        }

        cusparseDestroyDnVec(vecTmp1);
        cusparseDestroyDnVec(vecTmp2);
        iluReady = true;
    }

    // ── ILU(0) fill + factorize ────────────────────────────────────────
    // Assemble CSR values from diagH + coeffH, upload to d_iluVal, then
    // run in-place ILU(0) factorization.  Called from gpuUploadSolverArrays
    // every solve (coefficients change with pressure).
    // diagH[i]          = A[i,i]  (diagonal)
    // coeffH[fi*n+i]    = off-diagonal coupling for face fi at node i
    //                     (the SpMV computes A*v as: diag[i]*v[i] - coeff*v[nb],
    //                      so the matrix entry is A[i,j] = -coeff[fi*n+i])
    void fillAndFactorizeILU(const double* diagH,
                              const double* coeffH,
                              int           nFacesIn) {
        if (!iluReady) return;  // setupCSR not yet called or previously failed

        const int N = static_cast<int>(n);
        // Build CSR values on CPU.
        std::vector<double> hVal(static_cast<std::size_t>(nnz));
        for (int i = 0; i < N; ++i) {
            for (int k = h_csrRowPtr[i]; k < h_csrRowPtr[i + 1]; ++k) {
                const int fi = h_csrFaceIdx[k];
                if (fi < 0) {
                    // Diagonal entry.
                    hVal[k] = diagH[i];
                } else {
                    // Off-diagonal: A[i,j] = -coeff[fi*n+i].
                    hVal[k] = -coeffH[static_cast<std::size_t>(fi) * static_cast<std::size_t>(N) +
                                      static_cast<std::size_t>(i)];
                }
            }
        }
        LS_CUDA_CHECK(cudaMemcpy(d_iluVal, hVal.data(),
                                 static_cast<std::size_t>(nnz) * sizeof(double),
                                 cudaMemcpyHostToDevice));

        // In-place ILU(0) factorization: d_iluVal is overwritten with L\U factors.
        cusparseStatus_t st;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
        st = cusparseDcsrilu02(
            csHandle, N, nnz, descrA,
            d_iluVal, d_csrRowPtr, d_csrColInd,
            iluInfo, CUSPARSE_SOLVE_POLICY_USE_LEVEL, d_iluBuf);
#pragma GCC diagnostic pop
        if (st != CUSPARSE_STATUS_SUCCESS) {
            fprintf(stderr, "GpuBiCGSTABBuffers: cusparseDcsrilu02 failed\n");
            iluReady = false;
            return;
        }

        // Check for structural zero-pivot.
        int pivotPos = -1;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
        st = cusparseXcsrilu02_zeroPivot(csHandle, iluInfo, &pivotPos);
#pragma GCC diagnostic pop
        if (st == CUSPARSE_STATUS_ZERO_PIVOT) {
            fprintf(stderr, "GpuBiCGSTABBuffers: ILU(0) zero pivot at row %d\n",
                    pivotPos);
            iluReady = false;
        }
    }

    // ── Apply ILU(0) preconditioner: sol = (LU)^{-1} * rhs ───────────
    // Two-step triangular solve: L * tmp = rhs, U * sol = tmp.
    void applyILU(const double* rhs, double* sol) const {
        const double one = 1.0;
        const int64_t iN = static_cast<int64_t>(n);

        // Temporary dense-vector descriptors for this solve call.
        cusparseDnVecDescr_t vecRhs = nullptr, vecTmp = nullptr, vecSol = nullptr;
        cusparseCreateDnVec(&vecRhs, iN, const_cast<double*>(rhs), CUDA_R_64F);
        cusparseCreateDnVec(&vecTmp, iN, d_svTmp, CUDA_R_64F);
        cusparseCreateDnVec(&vecSol, iN, sol,      CUDA_R_64F);

        // Lower triangular solve: L * d_svTmp = rhs
        LS_CUSPARSE_CHECK(cusparseSpSV_solve(
            csHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
            spMatL, vecRhs, vecTmp, CUDA_R_64F,
            CUSPARSE_SPSV_ALG_DEFAULT, svDescrL));

        // Upper triangular solve: U * sol = d_svTmp
        LS_CUSPARSE_CHECK(cusparseSpSV_solve(
            csHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &one,
            spMatU, vecTmp, vecSol, CUDA_R_64F,
            CUSPARSE_SPSV_ALG_DEFAULT, svDescrU));

        cusparseDestroyDnVec(vecRhs);
        cusparseDestroyDnVec(vecTmp);
        cusparseDestroyDnVec(vecSol);
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
inline bool solveBiCGSTAB(
    GpuBiCGSTABBuffers&    gpu,
    std::vector<double>&   x,           // size n, host
    double                 diagEps,
    unsigned               maxIter,
    double                 tolerance,
    unsigned&              outIterations,
    double&                outResidual)
{
    assert(gpu.valid);
    const uint32_t n   = gpu.n;
    const int      blk = std::min(blocksFor(static_cast<int>(n)), 65535);

    // Upload initial guess to device
    LS_CUDA_CHECK(cudaMemcpy(gpu.d_x, x.data(), n * sizeof(double),
                             cudaMemcpyHostToDevice));
    // Zero all work vectors.  cudaMalloc does not zero GPU memory, and d_z in
    // particular is read (multiplied by omega=0) on the early-convergence path
    // before it has been written by jacobiKernel — giving 0*NaN = NaN if the
    // GPU heap happened to contain NaN from a prior allocation.
    LS_CUDA_CHECK(cudaMemset(gpu.d_p,   0, n * sizeof(double)));
    LS_CUDA_CHECK(cudaMemset(gpu.d_v,   0, n * sizeof(double)));
    LS_CUDA_CHECK(cudaMemset(gpu.d_y,   0, n * sizeof(double)));
    LS_CUDA_CHECK(cudaMemset(gpu.d_z,   0, n * sizeof(double)));
    LS_CUDA_CHECK(cudaMemset(gpu.d_s,   0, n * sizeof(double)));
    LS_CUDA_CHECK(cudaMemset(gpu.d_t,   0, n * sizeof(double)));
    LS_CUDA_CHECK(cudaMemset(gpu.d_Ax,  0, n * sizeof(double)));

    // r0 = b - A*x0
    gpu.spmv(gpu.d_x, gpu.d_Ax);
    residualKernel<<<blk, kThreadsPerBlock>>>(gpu.d_b, gpu.d_Ax, gpu.d_r, n);
    // r_hat = r0 (shadow residual; never changes)
    LS_CUDA_CHECK(cudaMemcpy(gpu.d_rhat, gpu.d_r, n * sizeof(double),
                             cudaMemcpyDeviceToDevice));

    // Normalise both convergence checks by ||b||_inf, matching the CPU path which
    // uses `residual < tolerance * b_norm`.  Without this the GPU checks the
    // raw (unnormalised) residual, which is 10–1000× tighter than the CPU
    // threshold.  BiCGSTAB then keeps iterating on an already-converged system
    // until breakdown (ρ → 0 → α = NaN) poisons the solution.
    // Both the s-check (early exit) and r-check (end-of-loop) use eff_tol.
    const double b_norm = gpu.maxAbs(gpu.d_b);
    if (!std::isfinite(b_norm) || !std::isfinite(tolerance) || tolerance <= 0.f)
        return false;
    const double eff_tol = (b_norm > 1e-100) ? tolerance * b_norm : tolerance;
    if (!std::isfinite(eff_tol) || eff_tol <= 0.f)
        return false;

    double rho   = 1.0;
    double alpha = 1.0;
    double omega = 1.0;
    bool converged = false;
    outIterations = 0;
    outResidual   = 0.0;

    // Restart counter: limits breakdown-restarts to avoid infinite loops.
    int restartCount = 0;
    static constexpr int kMaxRestarts = 3;

    outResidual = static_cast<double>(gpu.maxAbs(gpu.d_r));
    if (!std::isfinite(outResidual))
        return false;
    if (outResidual < eff_tol)
        return std::all_of(x.begin(), x.end(),
                           [](double value) { return std::isfinite(value); });

    for (unsigned iter = 0; iter < maxIter; ++iter) {

        // rho_new = <r_hat, r>
        const double rho_new =
            static_cast<double>(gpu.dotProduct(gpu.d_rhat, gpu.d_r));
        if (!std::isfinite(rho_new))
            return false;
        if (std::abs(rho_new) < 1e-100) {
            // ρ-breakdown: restart from current r.
            if (restartCount >= kMaxRestarts)
                return false;
            ++restartCount;
            // Reset shadow residual to current r and reinitialise scalars.
            LS_CUDA_CHECK(cudaMemcpy(gpu.d_rhat, gpu.d_r, n * sizeof(double),
                                     cudaMemcpyDeviceToDevice));
            LS_CUDA_CHECK(cudaMemset(gpu.d_p, 0, n * sizeof(double)));
            LS_CUDA_CHECK(cudaMemset(gpu.d_v, 0, n * sizeof(double)));
            rho = 1.0; alpha = 1.0; omega = 1.0;
            continue;
        }
        if (!std::isfinite(rho) || !std::isfinite(alpha) ||
            !std::isfinite(omega) || std::abs(omega) < 1e-100)
            return false;

        // p = r + beta*(p - omega*v),  beta = (rho_new/rho)*(alpha/omega)
        const double betaD = (rho_new / rho) * (alpha / omega);
        if (!std::isfinite(betaD))
            return false;
        const double beta = betaD;
        const double beta_omega = betaD * omega;
        if (!std::isfinite(beta) || !std::isfinite(beta_omega))
            return false;
        rho = rho_new;
        updatePKernel<<<blk, kThreadsPerBlock>>>(
            gpu.d_r, beta, beta_omega, gpu.d_v, gpu.d_p, n);

        // y = M^{-1} p
        if (gpu.iluReady)
            gpu.applyILU(gpu.d_p, gpu.d_y);
        else
            jacobiKernel<<<blk, kThreadsPerBlock>>>(
                gpu.d_p, gpu.d_diag, gpu.d_y, diagEps, n);

        // v = A y
        gpu.spmv(gpu.d_y, gpu.d_v);

        // alpha = rho_new / <r_hat, v>
        const double r_hat_v =
            static_cast<double>(gpu.dotProduct(gpu.d_rhat, gpu.d_v));
        if (!std::isfinite(r_hat_v))
            return false;
        if (std::abs(r_hat_v) < 1e-100) {
            // Lanczos breakdown: restart from current r.
            if (restartCount >= kMaxRestarts)
                return false;
            ++restartCount;
            LS_CUDA_CHECK(cudaMemcpy(gpu.d_rhat, gpu.d_r, n * sizeof(double),
                                     cudaMemcpyDeviceToDevice));
            LS_CUDA_CHECK(cudaMemset(gpu.d_p, 0, n * sizeof(double)));
            LS_CUDA_CHECK(cudaMemset(gpu.d_v, 0, n * sizeof(double)));
            rho = 1.0; alpha = 1.0; omega = 1.0;
            continue;
        }
        alpha = rho_new / r_hat_v;
        if (!std::isfinite(alpha))
            return false;

        // s = r - alpha*v
        subtractScaledKernel<<<blk, kThreadsPerBlock>>>(
            gpu.d_r, alpha, gpu.d_v, gpu.d_s, n);

        // Check convergence on s
        outResidual = static_cast<double>(gpu.maxAbs(gpu.d_s));
        if (!std::isfinite(outResidual))
            return false;
        if (outResidual < eff_tol) {
            // x += alpha*y  (z term is zero since we converged early)
            updateXKernel<<<blk, kThreadsPerBlock>>>(
                gpu.d_x, alpha, gpu.d_y,
                0.0, gpu.d_z, n);
            ++outIterations;
            converged = true;
            break;
        }

        // z = M^{-1} s
        if (gpu.iluReady)
            gpu.applyILU(gpu.d_s, gpu.d_z);
        else
            jacobiKernel<<<blk, kThreadsPerBlock>>>(
                gpu.d_s, gpu.d_diag, gpu.d_z, diagEps, n);

        // t = A z
        gpu.spmv(gpu.d_z, gpu.d_t);

        // omega = <t,s> / <t,t>
        const double t_s = static_cast<double>(gpu.dotProduct(gpu.d_t, gpu.d_s));
        const double t_t = static_cast<double>(gpu.dotProduct(gpu.d_t, gpu.d_t));
        if (!std::isfinite(t_s) || !std::isfinite(t_t) || t_t <= 1e-100)
            return false;
        omega = t_s / t_t;
        if (!std::isfinite(omega) || std::abs(omega) < 1e-100)
            return false;

        // x += alpha*y + omega*z
        updateXKernel<<<blk, kThreadsPerBlock>>>(
            gpu.d_x, alpha, gpu.d_y, omega, gpu.d_z, n);

        // r = s - omega*t
        updateRKernel<<<blk, kThreadsPerBlock>>>(
            gpu.d_s, omega, gpu.d_t, gpu.d_r, n);

        // Check convergence on r
        outResidual = static_cast<double>(gpu.maxAbs(gpu.d_r));
        if (!std::isfinite(outResidual))
            return false;
        ++outIterations;
        if (outResidual < eff_tol) {
            converged = true;
            break;
        }
    }

    if (!converged)
        return false;

    // Download solution to host
    LS_CUDA_CHECK(cudaMemcpy(x.data(), gpu.d_x, n * sizeof(double),
                             cudaMemcpyDeviceToHost));
    return std::all_of(x.begin(), x.end(),
                       [](double value) { return std::isfinite(value); });
}

} // namespace gpu
} // namespace viennals

#endif // VIENNALS_GPU_BICGSTAB
