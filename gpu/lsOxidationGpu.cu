// Compiled by nvcc into the ViennaLS_GPU shared library: the linear solves of
// the oxidation solvers on the GPU, exported through the C ABI in
// lsOxidationGpuAbi.hpp.
//
// Each solve is its CPU counterpart with AMGCL's CUDA backend instead of the
// builtin one: the same solver, preconditioner, parameters and stopping test
// on the same CSR system. AMGCL builds the AMG hierarchy on the host (OpenMP)
// and moves it to the device, where the Krylov iterations and V-cycles run.
// Its ILU0 smoother applies the factors with a few Jacobi sweeps on the
// device, where the CPU solves the triangular systems exactly: the iteration
// counts can differ slightly, the converged solution meets the same test.
//
// This library is the only binary that links against the CUDA runtime.
// ViennaLS opens it lazily with dlopen() (see lsOxidationGpuInterface.hpp)
// so that a GPU-enabled build still loads where no CUDA runtime is installed.

#ifndef AMGCL_NO_BOOST
#define AMGCL_NO_BOOST
#endif

#include <lsOxidationGpuAbi.hpp>

#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/cuda.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/preconditioner/schur_pressure_correction.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/fgmres.hpp>
#include <amgcl/solver/preonly.hpp>

#include <thrust/copy.h>
#include <thrust/device_vector.h>

#include <cmath>
#include <cstdio>
#include <exception>
#include <algorithm>
#include <memory>
#include <type_traits>
#include <new>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

typedef amgcl::backend::cuda<double> Backend;

thread_local std::string lastGpuError;

void setLastGpuError(const std::string &message) {
  lastGpuError = message;
  fprintf(stderr, "ViennaLS_GPU: %s\n", message.c_str());
}

// A CUDA device and context exist. cudaGetDeviceCount only enumerates
// devices; cudaFree(nullptr) forces the context (driver mismatch,
// permissions).
bool usableDevice() {
  int devices = 0;
  const cudaError_t countErr = cudaGetDeviceCount(&devices);
  if (countErr != cudaSuccess || devices == 0) {
    setLastGpuError(countErr != cudaSuccess
                        ? std::string("cudaGetDeviceCount failed: ") +
                              cudaGetErrorString(countErr)
                        : std::string("no CUDA device"));
    cudaGetLastError();
    return false;
  }
  const cudaError_t ctxErr = cudaFree(nullptr);
  if (ctxErr != cudaSuccess) {
    setLastGpuError(std::string("CUDA context initialization failed: ") +
                    cudaGetErrorString(ctxErr));
    cudaGetLastError();
    return false;
  }
  return true;
}

// One cuSPARSE handle for the diffusion solves of the process. It is never
// destroyed: at exit the CUDA runtime may be gone before static destructors.
cusparseHandle_t diffusionHandle() {
  static cusparseHandle_t handle = [] {
    cusparseHandle_t h = nullptr;
    if (cusparseCreate(&h) != CUSPARSE_STATUS_SUCCESS)
      h = nullptr;
    return h;
  }();
  return handle;
}

auto hostMatrix(int64_t n, const int64_t *ptr, const int64_t *col,
                const double *val) {
  return std::make_tuple(
      static_cast<ptrdiff_t>(n), amgcl::make_iterator_range(ptr, ptr + n + 1),
      amgcl::make_iterator_range(col, col + ptr[n]),
      amgcl::make_iterator_range(val, val + ptr[n]));
}

bool allFinite(const double *x, int64_t n) {
  for (int64_t i = 0; i < n; ++i)
    if (!std::isfinite(x[i]))
      return false;
  return true;
}

// ── diffusion (OxidationDiffusion::solveDiffusion) ──────────────────────────
typedef amgcl::make_solver<
    amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation,
               amgcl::relaxation::ilu0>,
    amgcl::solver::bicgstab<Backend>>
    DiffusionSolver;

// ── 3-D velocity-pressure (OxidationDeformation::macAmgSolve) ───────────────
typedef amgcl::make_solver<
    amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation,
               amgcl::relaxation::ilu0>,
    amgcl::solver::preonly<Backend>>
    USolver;
typedef amgcl::make_solver<
    amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::spai0>,
    amgcl::solver::preonly<Backend>>
    PSolver;
typedef amgcl::preconditioner::schur_pressure_correction<USolver, PSolver>
    MacPrecond;
typedef amgcl::solver::fgmres<Backend> MacKrylov;

// The system matrix of the 3-D solve on the device, kept between solves. On
// one node set its pattern rarely changes, and then only the values are
// copied (AMGCL's copy_matrix built a host CRS and allocated and converted
// ~10^7 entries on every solve). The product is AMGCL's cuda_matrix one.
struct MacDeviceMatrix {
  typedef double value_type;
  std::size_t nrows = 0, nnz = 0;
  cusparseHandle_t handle = nullptr;
  thrust::device_vector<int> ptr, col;
  thrust::device_vector<double> val;
  std::shared_ptr<std::remove_pointer<cusparseSpMatDescr_t>::type> desc;
  mutable thrust::device_vector<char> buf;
  std::vector<int64_t> hostPtr, hostCol; // the pattern on the device

  void set(int64_t n, const int64_t *p, const int64_t *c, const double *v,
           cusparseHandle_t h) {
    const auto rows = static_cast<std::size_t>(n);
    const auto entries = static_cast<std::size_t>(p[n]);
    if (desc && rows == nrows && entries == nnz && handle == h &&
        std::equal(p, p + n + 1, hostPtr.begin()) &&
        std::equal(c, c + entries, hostCol.begin())) {
      thrust::copy(v, v + entries, val.begin());
      return;
    }
    desc.reset(); // it points into the arrays replaced below
    hostPtr.assign(p, p + n + 1);
    hostCol.assign(c, c + entries);
    const std::vector<int> ptr32(p, p + n + 1), col32(c, c + entries);
    ptr.assign(ptr32.begin(), ptr32.end());
    col.assign(col32.begin(), col32.end());
    val.assign(v, v + entries);
    nrows = rows;
    nnz = entries;
    handle = h;
    desc.reset(amgcl::backend::detail::cuda_matrix_description(
                   nrows, nrows, nnz, ptr, col, val),
               amgcl::backend::detail::cuda_deleter());
  }

  void spmv(double alpha, const thrust::device_vector<double> &x, double beta,
            thrust::device_vector<double> &y) const {
    using namespace amgcl::backend::detail;
    std::shared_ptr<std::remove_pointer<cusparseDnVecDescr_t>::type> xdesc(
        cuda_vector_description(const_cast<thrust::device_vector<double> &>(x)),
        cuda_deleter());
    std::shared_ptr<std::remove_pointer<cusparseDnVecDescr_t>::type> ydesc(
        cuda_vector_description(y), cuda_deleter());
    std::size_t bufSize = 0;
    AMGCL_CALL_CUDA(cusparseSpMV_bufferSize(
        handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, desc.get(),
        xdesc.get(), &beta, ydesc.get(), cuda_datatype<double>(),
        CUSPARSE_SPMV_CSR_ALG1, &bufSize));
    if (buf.size() < bufSize)
      buf.resize(bufSize);
    AMGCL_CALL_CUDA(cusparseSpMV(
        handle, CUSPARSE_OPERATION_NON_TRANSPOSE, &alpha, desc.get(),
        xdesc.get(), &beta, ydesc.get(), cuda_datatype<double>(),
        CUSPARSE_SPMV_CSR_ALG1, thrust::raw_pointer_cast(buf.data())));
  }
};

} // namespace

// AMGCL's iterative solvers apply the system matrix through these two hooks
// (the same as for its cuda_matrix).
namespace amgcl {
namespace backend {
template <typename Alpha, typename Beta>
struct spmv_impl<Alpha, MacDeviceMatrix, thrust::device_vector<double>, Beta,
                 thrust::device_vector<double>> {
  static void apply(Alpha alpha, const MacDeviceMatrix &A,
                    const thrust::device_vector<double> &x, Beta beta,
                    thrust::device_vector<double> &y) {
    A.spmv(alpha, x, beta, y);
  }
};
template <>
struct residual_impl<MacDeviceMatrix, thrust::device_vector<double>,
                     thrust::device_vector<double>,
                     thrust::device_vector<double>> {
  static void apply(const thrust::device_vector<double> &rhs,
                    const MacDeviceMatrix &A,
                    const thrust::device_vector<double> &x,
                    thrust::device_vector<double> &r) {
    thrust::copy(rhs.begin(), rhs.end(), r.begin());
    A.spmv(-1, x, 1, r);
  }
};
} // namespace backend
} // namespace amgcl

namespace {

struct MacAmgGpu {
  Backend::params bprm;
  std::unique_ptr<MacPrecond> precond;
  std::unique_ptr<MacKrylov> krylov;
  MacDeviceMatrix matrix;
  ~MacAmgGpu() {
    krylov.reset();
    precond.reset();
    matrix.desc.reset();
    if (bprm.cusparse_handle)
      cusparseDestroy(bprm.cusparse_handle);
  }
};

} // namespace

extern "C" {

int viennalsGpuAbiVersion(void) { return VIENNALS_GPU_ABI_VERSION; }

const char *viennalsGpuGetLastErrorMessage(void) { return lastGpuError.c_str(); }

int viennalsGpuDiffusionSolve(int64_t n, const int64_t *ptr,
                              const int64_t *col, const double *val,
                              const double *rhs, double *x,
                              unsigned maxIterations, double absTolerance,
                              unsigned *outIterations) {
  lastGpuError.clear();
  if (!usableDevice())
    return 0;
  try {
    Backend::params bprm(diffusionHandle());
    if (!bprm.cusparse_handle)
      throw std::runtime_error("cusparseCreate failed");
    DiffusionSolver::params prm;
    prm.solver.maxiter = maxIterations;
    prm.solver.tol = 0.;
    prm.solver.abstol = absTolerance;
    DiffusionSolver solve(hostMatrix(n, ptr, col, val), prm, bprm);
    thrust::device_vector<double> f(rhs, rhs + n), u(x, x + n);
    std::size_t iterations = 0;
    double estimate = 0.;
    std::tie(iterations, estimate) = solve(f, u);
    thrust::copy(u.begin(), u.end(), x);
    *outIterations = static_cast<unsigned>(iterations);
  } catch (const std::exception &e) {
    setLastGpuError(std::string("diffusion solve failed: ") + e.what());
    cudaGetLastError();
    return -1;
  }
  if (!allFinite(x, n)) {
    setLastGpuError("diffusion solve gave non-finite values");
    return -1;
  }
  return 1;
}

void *viennalsGpuMacSetup(int64_t n, const int64_t *ptr, const int64_t *col,
                          const double *val, int64_t nvel, int dim,
                          const double *nullspace) {
  lastGpuError.clear();
  if (!usableDevice())
    return nullptr;
  auto *h = new (std::nothrow) MacAmgGpu();
  if (!h) {
    setLastGpuError("3-D mechanics setup: out of host memory");
    return nullptr;
  }
  try {
    if (cusparseCreate(&h->bprm.cusparse_handle) != CUSPARSE_STATUS_SUCCESS) {
      h->bprm.cusparse_handle = nullptr;
      throw std::runtime_error("cusparseCreate failed");
    }
    MacPrecond::params prm;
    prm.pmask.assign(static_cast<std::size_t>(n), 0);
    for (int64_t i = nvel; i < n; ++i)
      prm.pmask[static_cast<std::size_t>(i)] = 1;
    prm.type = 2;
    prm.approx_schur = true;
    prm.adjust_p = 2;
    auto &ns = prm.usolver.precond.coarsening.nullspace;
    ns.cols = dim;
    ns.B.assign(nullspace, nullspace + nvel * dim);
    h->precond = std::make_unique<MacPrecond>(hostMatrix(n, ptr, col, val),
                                              prm, h->bprm);
    MacKrylov::params sprm;
    sprm.M = 80;
    h->krylov = std::make_unique<MacKrylov>(static_cast<std::size_t>(n), sprm,
                                            h->bprm);
  } catch (const std::exception &e) {
    setLastGpuError(std::string("3-D mechanics setup failed: ") + e.what());
    delete h;
    cudaGetLastError();
    return nullptr;
  }
  return h;
}

int viennalsGpuMacSolve(void *handle, int64_t n, const int64_t *ptr,
                        const int64_t *col, const double *val,
                        const double *rhs, double *x, unsigned maxIterations,
                        double tolerance, unsigned *outIterations,
                        double *outResidual) {
  lastGpuError.clear();
  auto *h = static_cast<MacAmgGpu *>(handle);
  if (!h || !h->precond) {
    setLastGpuError("3-D mechanics solve: no setup");
    return 0;
  }
  try {
    h->matrix.set(n, ptr, col, val, h->bprm.cusparse_handle);
    thrust::device_vector<double> f(rhs, rhs + n), u(x, x + n);
    h->krylov->prm.maxiter = maxIterations;
    h->krylov->prm.tol = tolerance;
    std::size_t iterations = 0;
    double residual = 0.;
    std::tie(iterations, residual) =
        (*h->krylov)(h->matrix, *h->precond, f, u);
    thrust::copy(u.begin(), u.end(), x);
    *outIterations = static_cast<unsigned>(iterations);
    *outResidual = residual;
  } catch (const std::exception &e) {
    setLastGpuError(std::string("3-D mechanics solve failed: ") + e.what());
    cudaGetLastError();
    return 0;
  }
  if (!allFinite(x, n)) {
    setLastGpuError("3-D mechanics solve gave non-finite values");
    return 0;
  }
  return 1;
}

void viennalsGpuMacFree(void *handle) {
  delete static_cast<MacAmgGpu *>(handle);
}

} // extern "C"
