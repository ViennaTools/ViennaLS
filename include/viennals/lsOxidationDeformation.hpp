#pragma once

#include <lsOxidationDiffusion.hpp>
#include <lsOxidationSolverBase.hpp>

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include <omp.h>
#include <vcTimer.hpp>

namespace viennals {

/// Parameters for the Cartesian-grid oxide deformation model.
template <class T> struct OxidationDeformationParameters {
  T viscosity = 1.;
  T bulkModulus = 1.;
  T ambientPressure = 0.;
  T pressureTolerance = 1e-8;
  T minMechanicsBoundaryDistance = 0.05;
  T shearModulus = 0.;
  T stressRelaxationTime = 0.;
  T stressTimeStep = 1.;
  unsigned harmonicIterations = 500;
  unsigned mechanicsIterations = 5;
  unsigned pressureIterations = 10000;
  unsigned stokesIterations = 200;
  T mechanicsTolerance = 1e-8;
  T stokesTolerance = 1e-8;
  T tolerance = 1e-8;
  T relaxation = 0.7;         // SIMPLE velocity under-relaxation (0 < α ≤ 1)
  T pressureRelaxation = 0.5; // SIMPLE pressure under-relaxation (0 < β ≤ 1)
  std::size_t maxGridPoints = 5000000;
  int material = -1;
};

/// Propagates the volume expansion generated at the Si/SiO2 interface through
/// the oxide as a Cartesian-grid deformation velocity field.
///
/// This class implements the deformation part of the oxidation workflow in the
/// same fixed-grid spirit as the diffusion solver above. The Si/SiO2 boundary
/// is driven by the volume expansion velocity
///
///   (gamma - 1) / gamma * k C / N
///
/// along the reaction-interface normal. The initial velocity field in the oxide
/// is obtained from a component-wise harmonic extension and is then relaxed
/// with a Cartesian-grid quasi-static Stokes solve. The mechanical update
/// solves a pressure equation from the current velocity divergence and a
/// velocity momentum equation, including pressure-gradient and viscoelastic
/// deviatoric stress terms. The Si/SiO2 interface uses the oxidation expansion
/// velocity, the oxide/ambient interface uses a traction-free boundary, and
/// optional mask contacts use the mask velocity field.
template <class T, int D>
class OxidationDeformation final : public VelocityField<T>,
                                   public OxidationSolverBase<T, D> {
  using IndexType = viennahrle::Index<D>;
  using ConstSparseIterator =
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>;

private:
  // bring base members into scope
  using OxidationSolverBase<T, D>::noNode;
  using OxidationSolverBase<T, D>::nodeLookupFlat;
  using OxidationSolverBase<T, D>::initNodeLookup;
  using OxidationSolverBase<T, D>::lookupNode;
  using OxidationSolverBase<T, D>::minIndex;
  using OxidationSolverBase<T, D>::maxIndex;
  using OxidationSolverBase<T, D>::extents;
  using OxidationSolverBase<T, D>::strides;
  using OxidationSolverBase<T, D>::gridDelta;
  using OxidationSolverBase<T, D>::crosses;
  using OxidationSolverBase<T, D>::valueAt;
  using OxidationSolverBase<T, D>::inBounds;
  using OxidationSolverBase<T, D>::linearIndex;
  using OxidationSolverBase<T, D>::increment;
  using OxidationSolverBase<T, D>::findNearbyNode;
  using OxidationSolverBase<T, D>::initializeGridFromInterfaces;

  enum class Boundary { NONE, REACTION, AMBIENT, MASK };

  struct BoundaryIntersection {
    Boundary boundary = Boundary::NONE;
    T distance = 0.;
  };

  template <class ValueType> struct StencilPoint {
    ValueType value{};
    T distance = 1.;
  };

  struct Node {
    IndexType index;
    Vec3D<T> velocity{0., 0., 0.};
    T pressure = 0.;
    T strainTrace = 0.;
    std::array<T, 9> strainRateTensor{};
    std::array<T, 9> stressTensor{};
    T vonMisesStress = 0.;
  };

  SmartPointer<Domain<T, D>> reactionInterface = nullptr;
  SmartPointer<Domain<T, D>> ambientInterface = nullptr;
  SmartPointer<Domain<T, D>> maskInterface = nullptr;
  SmartPointer<OxidationDiffusion<T, D>> diffusionField = nullptr;
  SmartPointer<VelocityField<T>> maskVelocityField = nullptr;
  OxidationDeformationParameters<T> deformationParameters;
  OxidationParameters<T> oxidationParameters;
  int reactionSign = 1;
  int ambientSign = -1;
  int maskSign = 1;

  IndexType requestedMinIndex{};
  IndexType requestedMaxIndex{};
  unsigned iterations = 0;
  T residual = std::numeric_limits<T>::max();
  // Last achieved iteration counts and residuals for pressure and Stokes
  // solves.
  unsigned lastPressureIters_ = 0;
  T lastPressureResidual_ = 0.;
  unsigned lastStokesIters_ = 0;
  T lastStokesResidual_ = 0.;
  T avgExpansionSpeed_ = 0.;
  bool avgExpansionSpeedComputed = false;
  bool solved = false;
  bool nodesDirty_ = true;
  std::array<T, D> maxVelocity_{};
  bool useRequestedBounds = false;
  std::unordered_map<IndexType, std::array<T, 9>, detail::IndexTypeHasher<D>>
      deviatoricStressHistory;

  // Warm-start storage: solutions from previous time step used as initial guess
  std::vector<Vec3D<T>> previousVelocity_;
  std::vector<T> previousPressure_;
  bool hasPreviousSolution_ = false;

  static bool isFiniteVec(const Vec3D<T> &value) {
    for (unsigned i = 0; i < 3; ++i)
      if (!std::isfinite(value[i]))
        return false;
    return true;
  }

  static bool isFiniteTensor(const std::array<T, 9> &value) {
    for (const auto component : value)
      if (!std::isfinite(component))
        return false;
    return true;
  }

  // Face-major flat BC arrays: index = fi * n + nodeId, fi in [0, 2*D).
  std::vector<Boundary> faceBCTypes_;
  std::vector<T> faceBCDists_;
  std::vector<uint8_t>
      touchesAmbient_; // 1 if node touches the ambient (free) surface

  // GPU solver selection. Semantics match OxidationDiffusion.
  GpuMode gpuMode_ = GpuMode::Cpu;
  GpuPreconditioner gpuPreconditioner_ = GpuPreconditioner::Jacobi;
  mutable std::string
      lastLoggedBackend_; // suppresses repeated "using X" messages

#ifdef VIENNALS_GPU_BICGSTAB
  // Geometry-fixed (per buildNodes()) arrays for GPU pressure solve.
  // Layout matches pressCoeff/pressNeighId in solvePressure() so the
  // same spmvKernel can be reused without any CPU-side reformatting.
  std::vector<double> pressCoeffGpu_;    // face-major [2D * n]
  std::vector<uint32_t> pressNeighId32_; // face-major [2D * n]
  std::vector<double> actualDiagGpu_;    // [n], effective GPU matrix diagonal
  gpu::GpuBiCGSTABBuffers *gpuPressBufs_ = nullptr;

  // Geometry-fixed arrays for GPU Stokes velocity solve.  The diagonal is
  // component-major because mixed MASK contact is Dirichlet in the normal
  // component and Neumann/self-canceling in tangential components.
  std::vector<double> stokesCoeffGpu_;    // face-major [2D * n]
  std::vector<uint32_t> stokesNeighId32_; // face-major [2D * n]
  std::vector<double> stokesDiagGpu_;     // component-major [D * n]
  gpu::GpuBiCGSTABBuffers *gpuStokesBufs_ = nullptr;

  // Harmonic velocity solver GPU arrays.  The neighbor IDs are identical to
  // Stokes (stokesNeighId32_ is reused); only the coefficients and diagonal
  // differ (all interior coefficients = 1.0, diagonal = interior-face count).
  std::vector<double>
      harmonicCoeffGpu_; // face-major [2D * n], 1.0 for interior
  std::vector<double>
      harmonicDiagGpu_; // [n], = count of interior faces per node
  gpu::GpuBiCGSTABBuffers *gpuHarmonicBufs_ = nullptr;
#endif

public:
  std::vector<Node> nodes;

  OxidationDeformation() = default;

  OxidationDeformation(
      SmartPointer<Domain<T, D>> passedReactionInterface,
      SmartPointer<Domain<T, D>> passedAmbientInterface,
      SmartPointer<OxidationDiffusion<T, D>> passedDiffusionField,
      OxidationParameters<T> passedOxidationParameters,
      OxidationDeformationParameters<T> passedDeformationParameters = {})
      : reactionInterface(passedReactionInterface),
        ambientInterface(passedAmbientInterface),
        diffusionField(passedDiffusionField),
        deformationParameters(passedDeformationParameters),
        oxidationParameters(passedOxidationParameters) {}

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<OxidationDeformation>::New(std::forward<Args>(args)...);
  }

  ~OxidationDeformation() {
#ifdef VIENNALS_GPU_BICGSTAB
    gpu::freeGpuBuffers(gpuPressBufs_);
    gpu::freeGpuBuffers(gpuStokesBufs_);
    gpu::freeGpuBuffers(gpuHarmonicBufs_);
#endif
  }

  void setGpuMode(GpuMode mode) { gpuMode_ = mode; }
  void setGpuPreconditioner(GpuPreconditioner prec) {
    gpuPreconditioner_ = prec;
  }

  void setReactionInterface(SmartPointer<Domain<T, D>> passedInterface) {
    reactionInterface = passedInterface;
    nodesDirty_ = true;
    solved = false;
  }

  void setAmbientInterface(SmartPointer<Domain<T, D>> passedInterface) {
    ambientInterface = passedInterface;
    nodesDirty_ = true;
    solved = false;
  }

  void setMaskInterface(SmartPointer<Domain<T, D>> passedInterface,
                        int passedMaskSign = 1) {
    maskInterface = passedInterface;
    maskSign = (passedMaskSign < 0) ? -1 : 1;
    nodesDirty_ = true;
    solved = false;
  }

  void clearMaskInterface() {
    maskInterface = nullptr;
    nodesDirty_ = true;
    solved = false;
  }

  void
  setMaskVelocityField(SmartPointer<VelocityField<T>> passedVelocityField) {
    maskVelocityField = passedVelocityField;
    solved = false;
  }

  void clearMaskVelocityField() {
    maskVelocityField = nullptr;
    solved = false;
  }

  void setDiffusionField(
      SmartPointer<OxidationDiffusion<T, D>> passedDiffusionField) {
    diffusionField = passedDiffusionField;
    solved = false;
  }

  void setOxidationParameters(OxidationParameters<T> passedParameters) {
    oxidationParameters = passedParameters;
    solved = false;
  }

  void
  setDeformationParameters(OxidationDeformationParameters<T> passedParameters) {
    deformationParameters = passedParameters;
    solved = false;
  }

  void setOxideSigns(int passedReactionSign, int passedAmbientSign) {
    reactionSign = (passedReactionSign < 0) ? -1 : 1;
    ambientSign = (passedAmbientSign < 0) ? -1 : 1;
    solved = false;
  }

  void setSolveBounds(const IndexType &passedMinIndex,
                      const IndexType &passedMaxIndex) {
    requestedMinIndex = passedMinIndex;
    requestedMaxIndex = passedMaxIndex;
    useRequestedBounds = true;
    nodesDirty_ = true;
    solved = false;
  }

  void clearSolveBounds() {
    useRequestedBounds = false;
    nodesDirty_ = true;
    solved = false;
  }

  void markGeometryChanged() {
    nodesDirty_ = true;
    solved = false;
  }

  void apply() {
    if (reactionInterface == nullptr || ambientInterface == nullptr ||
        diffusionField == nullptr) {
      Logger::getInstance()
          .addError("OxidationDeformation: Missing interface or "
                    "diffusion field.")
          .print();
      return;
    }

    if (nodesDirty_) {
      if (!initialiseGrid())
        return; // base class already logged the error
      buildNodes();
      nodesDirty_ = false;
      if (nodes.empty())
        Logger::getInstance()
            .addWarning("OxidationDeformation: no oxide nodes found after "
                        "buildNodes(). Verify that the reaction and ambient "
                        "level sets enclose a non-empty oxide band.")
            .print();
      hasPreviousSolution_ =
          false; // Geometry changed, invalidate in-memory warm-start
      // Try restoring velocity, pressure, and stress history from level set
      // pointData (written by writeFieldsToLevelSet() before the previous
      // advection and remapped+filled by lsAdvect + lsInterior).
      seedFromLevelSet();
    } else if (hasPreviousSolution_ &&
               previousVelocity_.size() == nodes.size() &&
               previousPressure_.size() == nodes.size()) {
      // Warm-start: restore previous solution as initial guess for solver.
      // Geometry stability check: only warm-start if node count matches
      // (prevents using stale solutions after grid refinement/coarsening). This
      // typically reduces solver iterations by 30-50% since the previous step's
      // solution is close to the new one.
      for (std::size_t i = 0; i < nodes.size(); ++i) {
        nodes[i].velocity = previousVelocity_[i];
        nodes[i].pressure = previousPressure_[i];
      }
    }

    const std::size_t nn = nodes.size();
    Timer<> tHarmonic, tMechanics;
    tHarmonic.start();
    solveVelocity();
    tHarmonic.finish();
    tMechanics.start();
    solveMechanics();
    tMechanics.finish();
    Logger::getInstance()
        .addTiming("      deformation n=" + std::to_string(nn) + " harmonic",
                   tHarmonic)
        .addTiming("      deformation n=" + std::to_string(nn) +
                       " mechanics-total",
                   tMechanics)
        .print();
    avgExpansionSpeedComputed = false;

    maxVelocity_.fill(T(0));
    for (const auto &node : nodes) {
      for (unsigned d = 0; d < D; ++d)
        maxVelocity_[d] = std::max(maxVelocity_[d], std::abs(node.velocity[d]));
    }
    const auto unresolvedMax = estimateMaxUnresolvedAmbientVelocity();
    for (unsigned d = 0; d < D; ++d)
      maxVelocity_[d] = std::max(maxVelocity_[d], unresolvedMax[d]);

    // Save current solution for warm-start on next apply()
    previousVelocity_.resize(nodes.size());
    previousPressure_.resize(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      previousVelocity_[i] = nodes[i].velocity;
      previousPressure_[i] = nodes[i].pressure;
    }
    hasPreviousSolution_ = true;

    solved = true;
  }

  Vec3D<T> getVectorVelocity(const Vec3D<T> &coordinate, int material,
                             const Vec3D<T> & /*normalVector*/,
                             unsigned long /*pointId*/) final {
    if (!solved)
      apply();

    if (deformationParameters.material >= 0 &&
        material != deformationParameters.material)
      return {0., 0., 0.};

    const auto velocity = getVelocity(coordinate);
    T norm2 = 0.;
    for (unsigned d = 0; d < D; ++d)
      norm2 += velocity[d] * velocity[d];
    if (norm2 > std::numeric_limits<T>::epsilon())
      return velocity;

    return unresolvedAmbientVelocity(coordinate);
  }

  T getScalarVelocity(const Vec3D<T> &coordinate, int material,
                      const Vec3D<T> &normalVector,
                      unsigned long /*pointId*/) final {
    return 0.;
  }

  T getDissipationAlpha(int direction, int material,
                        const Vec3D<T> & /*centralDifferences*/) final {
    if (deformationParameters.material >= 0 &&
        material != deformationParameters.material)
      return 0.;
    return maxVelocity_[direction];
  }

private:
  template <class ValueType, class NodeAccessor>
  ValueType getField(const IndexType &index, ValueType fallback,
                     NodeAccessor accessor) const {
    const std::size_t nodeId = lookupNode(index);
    if (nodeId != noNode)
      return accessor(nodes[nodeId]);

    const auto nearby = findNearbyNode(index);
    if (nearby == noNode)
      return fallback;
    return accessor(nodes[nearby]);
  }

  template <class ValueType, class NodeAccessor>
  ValueType getField(const Vec3D<T> &coordinate, ValueType fallback,
                     NodeAccessor accessor) const {
    // D-linear interpolation for Vec3D<T> (velocity) and scalar T (pressure).
    // These are the types queried during level-set advection; interpolating at
    // the exact interface coordinate eliminates the O(gridDelta) nearest-node
    // error that caused the SiO2/mask interface to shift with grid delta.
    // Other types (e.g. stress tensor std::array<T,9>) use nearest-node because
    // they lack scalar multiply and are not used for advection.
    if constexpr (!std::is_same_v<ValueType, Vec3D<T>> &&
                  !std::is_same_v<ValueType, T>) {
      IndexType index;
      for (unsigned i = 0; i < D; ++i)
        index[i] = std::llround(coordinate[i] / gridDelta);
      return getField(index, fallback, accessor);
    } else {
      using IdxScalar = std::decay_t<decltype(std::declval<IndexType>()[0])>;
      IndexType lo;
      T frac[D];
      for (unsigned d = 0; d < D; ++d) {
        const T c = coordinate[d] / gridDelta;
        const T c_flo = std::floor(c);
        lo[d] = static_cast<IdxScalar>(c_flo);
        frac[d] = c - c_flo;
      }

      ValueType result{};
      T totalWeight = T(0);
      for (int corner = 0; corner < (1 << D); ++corner) {
        IndexType idx = lo;
        T w = T(1);
        for (unsigned d = 0; d < D; ++d) {
          if ((corner >> d) & 1) {
            idx[d]++;
            w *= frac[d];
          } else {
            w *= T(1) - frac[d];
          }
        }
        if (w < T(1e-14))
          continue;
        const std::size_t nodeId = lookupNode(idx);
        if (nodeId == noNode)
          continue;
        result = result + accessor(nodes[nodeId]) * w;
        totalWeight += w;
      }

      if (totalWeight < T(1e-14))
        return fallback;
      if (totalWeight < T(1) - T(1e-6))
        result = result * (T(1) / totalWeight);
      return result;
    }
  }

public:
  Vec3D<T> getVelocity(const Vec3D<T> &coordinate) const {
    return getField(coordinate, Vec3D<T>{0., 0., 0.},
                    [](const Node &n) { return n.velocity; });
  }

  Vec3D<T> getVelocity(const IndexType &index) const {
    return getField(index, Vec3D<T>{0., 0., 0.},
                    [](const Node &n) { return n.velocity; });
  }

  T getPressure(const Vec3D<T> &coordinate) const {
    return getField(coordinate, T(0), [](const Node &n) { return n.pressure; });
  }

  T getPressure(const IndexType &index) const {
    return getField(index, T(0), [](const Node &n) { return n.pressure; });
  }

  T getStrainTrace(const Vec3D<T> &coordinate) const {
    return getField(coordinate, T(0),
                    [](const Node &n) { return n.strainTrace; });
  }

  T getStrainTrace(const IndexType &index) const {
    return getField(index, T(0), [](const Node &n) { return n.strainTrace; });
  }

  std::array<T, 9> getStrainRateTensor(const Vec3D<T> &coordinate) const {
    return getField(coordinate, std::array<T, 9>{},
                    [](const Node &n) { return n.strainRateTensor; });
  }

  std::array<T, 9> getStrainRateTensor(const IndexType &index) const {
    return getField(index, std::array<T, 9>{},
                    [](const Node &n) { return n.strainRateTensor; });
  }

  std::array<T, 9> getStressTensor(const Vec3D<T> &coordinate) const {
    return getField(coordinate, std::array<T, 9>{},
                    [](const Node &n) { return n.stressTensor; });
  }

  std::array<T, 9> getStressTensor(const IndexType &index) const {
    return getField(index, std::array<T, 9>{},
                    [](const Node &n) { return n.stressTensor; });
  }

  T getVonMisesStress(const Vec3D<T> &coordinate) const {
    return getField(coordinate, T(0),
                    [](const Node &n) { return n.vonMisesStress; });
  }

  T getVonMisesStress(const IndexType &index) const {
    return getField(index, T(0),
                    [](const Node &n) { return n.vonMisesStress; });
  }

  unsigned getIterations() const { return iterations; }
  T getResidual() const { return residual; }
  T getLastPressureResidual() const { return lastPressureResidual_; }
  T getLastStokesResidual() const { return lastStokesResidual_; }
  bool lastSolveConverged() const {
    return std::isfinite(residual) &&
           residual <= deformationParameters.mechanicsTolerance &&
           std::isfinite(lastPressureResidual_) &&
           lastPressureResidual_ <= deformationParameters.pressureTolerance &&
           std::isfinite(lastStokesResidual_) &&
           lastStokesResidual_ <= deformationParameters.stokesTolerance &&
           hasFiniteSolution();
  }
  bool hasFiniteSolution() const {
    for (const auto &node : nodes) {
      if (!isFiniteVec(node.velocity) || !std::isfinite(node.pressure) ||
          !std::isfinite(node.strainTrace) ||
          !isFiniteTensor(node.strainRateTensor) ||
          !isFiniteTensor(node.stressTensor) ||
          !std::isfinite(node.vonMisesStress))
        return false;
    }
    return true;
  }
  std::size_t getNumberOfSolutionNodes() const { return nodes.size(); }
  T avgExpansionSpeed() {
    if (!avgExpansionSpeedComputed) {
      computeAvgExpansionSpeed();
      avgExpansionSpeedComputed = true;
    }
    return avgExpansionSpeed_;
  }
  template <class Callback> void forEachSolutionNode(Callback callback) const {
    for (const auto &node : nodes)
      callback(node.index, node.pressure);
  }

  /// Write velocity (Vec3D) and viscoelastic stress history (3 tensor-row
  /// vectors) into ambientInterface->getPointData(). Called alongside
  /// OxidationDiffusion::writePersistentFields() so a single lsInterior pass
  /// fills the oxide interior with all warm-start data.
  void writeFieldsToLevelSet() {
    if (nodes.empty() || ambientInterface == nullptr)
      return;

    using VD = typename PointData<T>::VectorDataType;
    VD velocity, stressR0, stressR1, stressR2;

    ConstSparseIterator it(ambientInterface->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const IndexType idx = it.getStartIndices();
      const std::size_t nId = lookupNode(idx);

      if (nId != noNode) {
        const auto &n = nodes[nId];
        velocity.push_back(
            isFiniteVec(n.velocity) ? n.velocity : Vec3D<T>{T(0), T(0), T(0)});
        const auto sIt = deviatoricStressHistory.find(idx);
        if (sIt != deviatoricStressHistory.end() &&
            isFiniteTensor(sIt->second)) {
          const auto &s = sIt->second;
          stressR0.push_back({s[0], s[1], s[2]});
          stressR1.push_back({s[3], s[4], s[5]});
          stressR2.push_back({s[6], s[7], s[8]});
        } else {
          stressR0.push_back({T(0), T(0), T(0)});
          stressR1.push_back({T(0), T(0), T(0)});
          stressR2.push_back({T(0), T(0), T(0)});
        }
      } else {
        velocity.push_back({T(0), T(0), T(0)});
        stressR0.push_back({T(0), T(0), T(0)});
        stressR1.push_back({T(0), T(0), T(0)});
        stressR2.push_back({T(0), T(0), T(0)});
      }
    }

    auto &pd = ambientInterface->getPointData();
    pd.insertReplaceVectorData(std::move(velocity), "OxVelocity");
    pd.insertReplaceVectorData(std::move(stressR0), "OxStressR0");
    pd.insertReplaceVectorData(std::move(stressR1), "OxStressR1");
    pd.insertReplaceVectorData(std::move(stressR2), "OxStressR2");
  }

private:
  /// Reads velocity, pressure (shared with diffusion's "OxPressure"), and
  /// stress history from ambientInterface pointData into the warm-start state.
  /// Called in apply() after buildNodes() when hasPreviousSolution_=false.
  void seedFromLevelSet() {
    if (ambientInterface == nullptr || nodes.empty())
      return;

    auto &pd = ambientInterface->getPointData();
    const int vIdx = pd.getVectorDataIndex("OxVelocity");
    const int r0Idx = pd.getVectorDataIndex("OxStressR0");
    const int r1Idx = pd.getVectorDataIndex("OxStressR1");
    const int r2Idx = pd.getVectorDataIndex("OxStressR2");
    const int pIdx = pd.getScalarDataIndex("OxPressure");

    const bool hasVelocity = (vIdx != -1);
    const bool hasStress = (r0Idx != -1 && r1Idx != -1 && r2Idx != -1);
    const bool hasPressure = (pIdx != -1);

    if (!hasVelocity && !hasStress && !hasPressure)
      return;

    const auto *vd = hasVelocity ? pd.getVectorData(vIdx) : nullptr;
    const auto *r0d = hasStress ? pd.getVectorData(r0Idx) : nullptr;
    const auto *r1d = hasStress ? pd.getVectorData(r1Idx) : nullptr;
    const auto *r2d = hasStress ? pd.getVectorData(r2Idx) : nullptr;
    const auto *ppd = hasPressure ? pd.getScalarData(pIdx) : nullptr;

    previousVelocity_.assign(nodes.size(), Vec3D<T>{});
    previousPressure_.assign(nodes.size(), T(0));

    ConstSparseIterator it(ambientInterface->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const auto ptId = it.getPointId();
      const IndexType idx = it.getStartIndices();
      const std::size_t ni = lookupNode(idx);

      if (ni != noNode) {
        if (vd && ptId < static_cast<decltype(ptId)>(vd->size()) &&
            isFiniteVec((*vd)[ptId]))
          previousVelocity_[ni] = (*vd)[ptId];
        if (ppd && ptId < static_cast<decltype(ptId)>(ppd->size()) &&
            std::isfinite((*ppd)[ptId]))
          previousPressure_[ni] = (*ppd)[ptId];
      }

      if (hasStress && ptId < static_cast<decltype(ptId)>(r0d->size()) &&
          ptId < static_cast<decltype(ptId)>(r1d->size()) &&
          ptId < static_cast<decltype(ptId)>(r2d->size())) {
        const auto &row0 = (*r0d)[ptId];
        const auto &row1 = (*r1d)[ptId];
        const auto &row2 = (*r2d)[ptId];
        std::array<T, 9> s{row0[0], row0[1], row0[2], row1[0], row1[1],
                           row1[2], row2[0], row2[1], row2[2]};
        if (isFiniteTensor(s))
          deviatoricStressHistory[idx] = s;
      }
    }

    if (hasVelocity || hasPressure) {
      hasPreviousSolution_ = true;
      // Apply the restored state immediately to the current nodes so
      // the Stokes solve warm-starts on this (nodesDirty_=true) apply() call.
      for (std::size_t i = 0; i < nodes.size(); ++i) {
        nodes[i].velocity = previousVelocity_[i];
        nodes[i].pressure = previousPressure_[i];
      }
    }
  }

#ifdef VIENNALS_GPU_BICGSTAB
  // Builds the face-major pressure matrix geometry for GPU reuse.
  // The logic mirrors the explicit precomputation in solvePressure() so the
  // same actualDiagGpu_ / pressCoeffGpu_ arrays can feed both the CPU ILU(0)
  // factorization and the GPU SpMV without recomputation.
  void buildPressureGpuGeometry() {
    const std::size_t n = nodes.size();
    const T eps = std::numeric_limits<T>::epsilon();
    pressCoeffGpu_.assign(2 * D * n, 0.0);
    pressNeighId32_.assign(2 * D * n, gpu::kNoNode);
    actualDiagGpu_.assign(n, 0.0);

    for (std::size_t id = 0; id < n; ++id) {
      if (touchesAmbient_[id]) {
        actualDiagGpu_[id] = 1.0;
        continue;
      }
      for (unsigned dir = 0; dir < D; ++dir) {
        const unsigned fiNeg = dir * 2u;
        const unsigned fiPos = dir * 2u + 1u;
        IndexType nbNeg = nodes[id].index;
        nbNeg[dir] -= 1;
        IndexType nbPos = nodes[id].index;
        nbPos[dir] += 1;
        const std::size_t jNeg =
            inBounds(nbNeg) ? nodeLookupFlat[linearIndex(nbNeg)] : noNode;
        const std::size_t jPos =
            inBounds(nbPos) ? nodeLookupFlat[linearIndex(nbPos)] : noNode;

        auto effDist = [&](unsigned fi, std::size_t j) -> T {
          if (j != noNode)
            return gridDelta;
          const Boundary bt = faceBCTypes_[fi * n + id];
          return (bt != Boundary::NONE) ? faceBCDists_[fi * n + id] : gridDelta;
        };

        const T dNeg = effDist(fiNeg, jNeg);
        const T dPos = effDist(fiPos, jPos);
        const T dSum = dNeg + dPos;
        if (dSum <= eps)
          continue;

        auto processFace = [&](unsigned fi, std::size_t j, T d) {
          const T c = T(2) / (d * dSum);
          if (j != noNode && !touchesAmbient_[j]) {
            pressCoeffGpu_[fi * n + id] = static_cast<double>(c);
            pressNeighId32_[fi * n + id] = static_cast<uint32_t>(j);
            actualDiagGpu_[id] += c;
          } else if (j != noNode ||
                     faceBCTypes_[fi * n + id] == Boundary::AMBIENT) {
            // j is an ambient-only neighbour (identity-row Dirichlet p=0), OR
            // this face crosses the free surface directly (AMBIENT Dirichlet
            // p=0 at the sub-grid crossing distance). REACTION faces are
            // solid-wall Neumann ∂p/∂n=0: no contribution.
            actualDiagGpu_[id] += c;
          }
        };
        processFace(fiNeg, jNeg, dNeg);
        processFace(fiPos, jPos, dPos);
      }
      if (actualDiagGpu_[id] <= eps)
        actualDiagGpu_[id] = 1.0;
    }
  }

  // Builds the face-major Stokes velocity matrix geometry for GPU reuse.
  // The effective diagonal excludes OOB, AMBIENT, and traction-coupled MASK
  // face self-coupling because those cancel exactly with the vBC correction in
  // the Stokes matvec.  Kinematic MASK faces remain Dirichlet-like and
  // contribute to the diagonal.
  void buildStokesGpuGeometry() {
    const std::size_t n = nodes.size();
    const T eps = std::numeric_limits<T>::epsilon();
    stokesCoeffGpu_.assign(2 * D * n, 0.0);
    stokesNeighId32_.assign(2 * D * n, gpu::kNoNode);
    stokesDiagGpu_.assign(D * n, 0.0);

    auto addDiag = [&](std::size_t id, T c) {
      for (unsigned comp = 0; comp < D; ++comp)
        stokesDiagGpu_[comp * n + id] += static_cast<double>(c);
    };

    for (std::size_t id = 0; id < n; ++id) {
      for (unsigned dir = 0; dir < D; ++dir) {
        const unsigned fiNeg = dir * 2u;
        const unsigned fiPos = dir * 2u + 1u;
        IndexType nbNeg = nodes[id].index;
        nbNeg[dir] -= 1;
        IndexType nbPos = nodes[id].index;
        nbPos[dir] += 1;
        const bool negInBounds = inBounds(nbNeg);
        const bool posInBounds = inBounds(nbPos);
        const std::size_t jNeg =
            negInBounds ? nodeLookupFlat[linearIndex(nbNeg)] : noNode;
        const std::size_t jPos =
            posInBounds ? nodeLookupFlat[linearIndex(nbPos)] : noNode;

        // Distance matching velocityStencilPoint: gridDelta for OOB/interior,
        // faceBCDists_ for boundary-crossing faces.
        auto stokesFaceDist = [&](unsigned fi, bool inb, std::size_t j) -> T {
          if (!inb || j != noNode)
            return gridDelta;
          const Boundary bt = faceBCTypes_[fi * n + id];
          return (bt != Boundary::NONE) ? faceBCDists_[fi * n + id] : gridDelta;
        };

        const T dNeg = stokesFaceDist(fiNeg, negInBounds, jNeg);
        const T dPos = stokesFaceDist(fiPos, posInBounds, jPos);
        const T dSum = dNeg + dPos;
        if (dSum <= eps)
          continue;

        auto processFace = [&](unsigned fi, bool inb, std::size_t j, T d) {
          const T c = T(2) / (d * dSum);
          if (inb && j != noNode) {
            // Interior neighbor: off-diagonal and diagonal contribution.
            stokesCoeffGpu_[fi * n + id] = static_cast<double>(c);
            stokesNeighId32_[fi * n + id] = static_cast<uint32_t>(j);
            addDiag(id, c);
          } else if (inb) {
            // Boundary-crossing face.
            // REACTION / MASK: Dirichlet → adds to diagonal.
            // AMBIENT / OOB: excluded (self-coupling cancels with vBC
            // correction).
            const Boundary bt = faceBCTypes_[fi * n + id];
            if (bt == Boundary::REACTION || bt == Boundary::MASK)
              addDiag(id, c);
          }
          // !inb (OOB): self-coupling, excluded.
        };
        processFace(fiNeg, negInBounds, jNeg, dNeg);
        processFace(fiPos, posInBounds, jPos, dPos);
      }
      for (unsigned comp = 0; comp < D; ++comp)
        if (stokesDiagGpu_[comp * n + id] <= static_cast<double>(eps))
          stokesDiagGpu_[comp * n + id] = 1.0;
    }
  }

  // Builds face-major harmonic velocity geometry for GPU reuse.
  // The neighbor IDs are identical to Stokes (stokesNeighId32_ is shared).
  //
  // The CPU harmonicMatvec is  Av[i][c] = 2*D*v[i] - stencil_sum(v)[i] + b[i]
  // where stencil_sum includes self-coupling (v[nodeId]) for OOB/AMBIENT/NONE
  // faces and constants for REACTION/MASK faces.  After the b/constant terms
  // cancel, the effective matrix diagonal is:
  //   2*D - n_OOB - n_AMBIENT - n_NONE  =  n_interior + n_REACTION + n_MASK
  //
  // REACTION and MASK faces have no off-diagonal entry (their velocity is a
  // constant in b) but they DO count toward the diagonal.  Excluding them
  // gives a matrix that is not diagonally dominant near the Si/SiO2 boundary,
  // which causes Jacobi-preconditioned BiCGSTAB to diverge.
  void buildHarmonicGpuGeometry() {
    const std::size_t n = nodes.size();
    harmonicCoeffGpu_.assign(2 * D * n, 0.0);
    harmonicDiagGpu_.assign(n, 0.0);
    for (std::size_t id = 0; id < n; ++id) {
      for (unsigned fi = 0; fi < 2 * D; ++fi) {
        if (stokesNeighId32_[fi * n + id] != gpu::kNoNode) {
          // Interior neighbor: off-diagonal coefficient 1.0 + diagonal 1.0.
          harmonicCoeffGpu_[fi * n + id] = 1.0;
          harmonicDiagGpu_[id] += 1.0;
        } else {
          // Non-interior face.  REACTION/MASK contribute a Dirichlet constant
          // to b but NOT self-coupling, so they add 1.0 to the diagonal just
          // like an interior neighbor would (matching the CPU matrix row).
          // OOB/AMBIENT/NONE contribute v[nodeId] (self-coupling), which
          // cancels with the 2*D diagonal term and must be excluded here.
          const Boundary bt = faceBCTypes_[fi * n + id];
          if (bt == Boundary::REACTION || bt == Boundary::MASK)
            harmonicDiagGpu_[id] += 1.0;
        }
      }
      if (harmonicDiagGpu_[id] < 1e-10)
        harmonicDiagGpu_[id] = 1.0;
    }
  }

  // Allocates GPU buffers for pressure and Stokes solves and uploads the
  // geometry-fixed CSR pattern.  Called from buildNodes() after the face BC
  // arrays and the two geometry helper arrays are populated.
  void setupDeformationGpuBuffers() {
    const std::size_t n = nodes.size();
    if (n == 0) {
      gpu::freeGpuBuffers(gpuPressBufs_);
      gpuPressBufs_ = nullptr;
      gpu::freeGpuBuffers(gpuStokesBufs_);
      gpuStokesBufs_ = nullptr;
      return;
    }
    const bool tryGpu = (gpuMode_ == GpuMode::Gpu);
    const bool useIlu0 = (gpuPreconditioner_ == GpuPreconditioner::ILU0);

    gpu::freeGpuBuffers(gpuPressBufs_);
    gpuPressBufs_ = nullptr;
    if (tryGpu) {
      gpuPressBufs_ =
          gpu::allocGpuBuffers(static_cast<uint32_t>(n), 2 * D, useIlu0);
      if (!gpuPressBufs_) {
        VIENNACORE_LOG_ERROR("OxidationDeformation: GPU mode was selected, but "
                             "pressure solver CUDA buffers could not be "
                             "allocated or the CUDA context could not be "
                             "initialized." +
                             gpuErrorDetail());
      }
      if (!gpu::gpuUploadNeighborIds(gpuPressBufs_, pressNeighId32_.data(),
                                     2u * D * n)) {
        gpu::freeGpuBuffers(gpuPressBufs_);
        gpuPressBufs_ = nullptr;
        VIENNACORE_LOG_ERROR("OxidationDeformation: GPU mode was selected, but "
                             "uploading pressure GPU neighbor IDs failed." +
                             gpuErrorDetail());
      }
      if (useIlu0 && !gpu::gpuSetupCSR(gpuPressBufs_, pressNeighId32_.data(),
                                       static_cast<uint32_t>(n), 2 * D)) {
        gpu::freeGpuBuffers(gpuPressBufs_);
        gpuPressBufs_ = nullptr;
        VIENNACORE_LOG_ERROR("OxidationDeformation: GPU mode was selected, but "
                             "CUSPARSE setup for the pressure GPU BiCGSTAB "
                             "solver failed." +
                             gpuErrorDetail());
      }
    }

    gpu::freeGpuBuffers(gpuStokesBufs_);
    gpuStokesBufs_ = nullptr;
    if (tryGpu) {
      gpuStokesBufs_ =
          gpu::allocGpuBuffers(static_cast<uint32_t>(n), 2 * D, useIlu0);
      if (!gpuStokesBufs_) {
        VIENNACORE_LOG_ERROR("OxidationDeformation: GPU mode was selected, but "
                             "Stokes solver CUDA buffers could not be "
                             "allocated or the CUDA context could not be "
                             "initialized." +
                             gpuErrorDetail());
      }
      if (!gpu::gpuUploadNeighborIds(gpuStokesBufs_, stokesNeighId32_.data(),
                                     2u * D * n)) {
        gpu::freeGpuBuffers(gpuStokesBufs_);
        gpuStokesBufs_ = nullptr;
        VIENNACORE_LOG_ERROR("OxidationDeformation: GPU mode was selected, but "
                             "uploading Stokes GPU neighbor IDs failed." +
                             gpuErrorDetail());
      }
      if (useIlu0 && !gpu::gpuSetupCSR(gpuStokesBufs_, stokesNeighId32_.data(),
                                       static_cast<uint32_t>(n), 2 * D)) {
        gpu::freeGpuBuffers(gpuStokesBufs_);
        gpuStokesBufs_ = nullptr;
        VIENNACORE_LOG_ERROR("OxidationDeformation: GPU mode was selected, but "
                             "CUSPARSE setup for the Stokes GPU BiCGSTAB "
                             "solver failed." +
                             gpuErrorDetail());
      }
    }

    // Harmonic velocity: same neighbor IDs as Stokes, different coefficients.
    gpu::freeGpuBuffers(gpuHarmonicBufs_);
    gpuHarmonicBufs_ = nullptr;
    if (tryGpu) {
      gpuHarmonicBufs_ =
          gpu::allocGpuBuffers(static_cast<uint32_t>(n), 2 * D, useIlu0);
      if (!gpuHarmonicBufs_) {
        VIENNACORE_LOG_ERROR("OxidationDeformation: GPU mode was selected, but "
                             "harmonic solver CUDA buffers could not be "
                             "allocated." +
                             gpuErrorDetail());
      }
      if (!gpu::gpuUploadNeighborIds(gpuHarmonicBufs_, stokesNeighId32_.data(),
                                     2u * D * n)) {
        gpu::freeGpuBuffers(gpuHarmonicBufs_);
        gpuHarmonicBufs_ = nullptr;
        VIENNACORE_LOG_ERROR("OxidationDeformation: GPU mode was selected, but "
                             "uploading harmonic GPU neighbor IDs failed." +
                             gpuErrorDetail());
      }
      if (useIlu0 &&
          !gpu::gpuSetupCSR(gpuHarmonicBufs_, stokesNeighId32_.data(),
                            static_cast<uint32_t>(n), 2 * D)) {
        gpu::freeGpuBuffers(gpuHarmonicBufs_);
        gpuHarmonicBufs_ = nullptr;
        VIENNACORE_LOG_ERROR("OxidationDeformation: GPU mode was selected, but "
                             "CUSPARSE setup for the harmonic GPU BiCGSTAB "
                             "solver failed." +
                             gpuErrorDetail());
      }
    }

    if (tryGpu)
      logDeformationBackend("GPU BiCGSTAB",
                            "pressure/Stokes/harmonic, preconditioner=" +
                                std::string(useIlu0 ? "ILU0" : "Jacobi"));
  }
#endif // VIENNALS_GPU_BICGSTAB

#ifdef VIENNALS_GPU_BICGSTAB
  static std::string gpuErrorDetail() {
    const char *detail = gpu::gpuGetLastErrorMessage();
    if (detail && detail[0] != '\0')
      return std::string(" Detail: ") + detail;
    return {};
  }
#endif

  void logDeformationBackend(const std::string &backend,
                             const std::string &detail) const {
    if (!Logger::hasInfo())
      return;
    const std::string msg = "OxidationDeformation: using " + backend +
                            " for mechanics pressure/Stokes solves (nodes=" +
                            std::to_string(nodes.size()) +
                            (detail.empty() ? std::string() : ", " + detail) +
                            ").";
    if (msg == lastLoggedBackend_)
      return;
    lastLoggedBackend_ = msg;
    Logger::getInstance().addInfo(msg).print();
  }

public:
  bool initialiseGrid() {
    return initializeGridFromInterfaces(
        reactionInterface, ambientInterface, maskInterface, useRequestedBounds,
        requestedMinIndex, requestedMaxIndex,
        deformationParameters.maxGridPoints, "OxidationDeformation");
  }

  void buildNodes() {
    nodes.clear();
    initNodeLookup();

    ConstSparseIterator reactionIt(reactionInterface->getDomain());
    ConstSparseIterator ambientIt(ambientInterface->getDomain());
    auto maskIt = makeMaskIterator();

    IndexType index = minIndex;
    while (true) {
      const T reactionPhi = valueAt(reactionIt, index);
      const T ambientPhi = valueAt(ambientIt, index);
      if (isInsideOxide(reactionPhi, ambientPhi) &&
          !isInsideMask(maskIt, index)) {
        const std::size_t id = nodes.size();
        nodeLookupFlat[linearIndex(index)] = id;
        nodes.push_back({index});
      }

      if (!increment(index))
        break;
    }

    // Precompute per-face boundary intersections into flat face-major arrays.
    const std::size_t n = nodes.size();
    faceBCTypes_.assign(2 * D * n, Boundary::NONE);
    faceBCDists_.assign(2 * D * n, T(1));
    touchesAmbient_.assign(n, uint8_t(0));
    for (std::size_t id = 0; id < n; ++id) {
      const auto &node = nodes[id];
      bool touchesAmbient = false;
      bool touchesSolidBoundary = false;
      for (unsigned dir = 0; dir < D; ++dir) {
        for (int off : {-1, 1}) {
          const unsigned fi = dir * 2u + (off == 1 ? 1u : 0u);
          IndexType nb = node.index;
          nb[dir] += off;
          if (!inBounds(nb) || lookupNode(nb) != noNode)
            continue; // NONE/1.0 already set
          const auto bi = boundaryIntersection(reactionIt, ambientIt, maskIt,
                                               node.index, nb);
          faceBCTypes_[fi * n + id] = bi.boundary;
          faceBCDists_[fi * n + id] = bi.distance;
          if (bi.boundary == Boundary::AMBIENT)
            touchesAmbient = true;
          else if (bi.boundary == Boundary::MASK ||
                   bi.boundary == Boundary::REACTION)
            touchesSolidBoundary = true;
        }
      }
      // Do not collapse mixed mask/reaction/ambient corner nodes to a single
      // ambient pressure identity row.  Those nodes need their per-face
      // boundary conditions; otherwise the mask-edge triple point injects an
      // artificial free-surface pressure singularity.
      touchesAmbient_[id] =
          (touchesAmbient && !touchesSolidBoundary) ? uint8_t(1) : uint8_t(0);
    }

#ifdef VIENNALS_GPU_BICGSTAB
    buildPressureGpuGeometry();
    buildStokesGpuGeometry();
    buildHarmonicGpuGeometry();
    setupDeformationGpuBuffers();
#else
    if (gpuMode_ == GpuMode::Gpu) {
      VIENNACORE_LOG_ERROR("OxidationDeformation: explicit GPU mode was "
                           "requested, but ViennaLS was built without "
                           "VIENNALS_GPU_BICGSTAB.");
    }
#endif
  }

  // Evaluates the harmonic stencil at one node.
  // sum = sum of neighbor/BC contributions (interior neighbors, reaction/mask
  // BCs, and self-coupling for OOB/NONE/AMBIENT faces). count is always 2*D
  // (every face is counted regardless of type).
  template <class SolverT>
  void computeHarmonicStencilAt(std::size_t nodeId,
                                const std::vector<Vec3D<SolverT>> &v,
                                Vec3D<T> &sum) const {
    const auto &node = nodes[nodeId];
    sum = {T(0), T(0), T(0)};

    const auto toT = [](const Vec3D<SolverT> &w) -> Vec3D<T> {
      return {static_cast<T>(w[0]), static_cast<T>(w[1]), static_cast<T>(w[2])};
    };

    for (unsigned direction = 0; direction < D; ++direction) {
      for (int offset : {-1, 1}) {
        IndexType neighbor = node.index;
        neighbor[direction] += offset;

        if (!inBounds(neighbor)) {
          detail::vecAddTo(sum, toT(v[nodeId])); // zero-flux: ghost = self
          continue;
        }

        const std::size_t neighborId = lookupNode(neighbor);
        if (neighborId != noNode) {
          detail::vecAddTo(sum, toT(v[neighborId]));
          continue;
        }

        const unsigned fi = direction * 2u + (offset == 1 ? 1u : 0u);
        const Boundary boundary = faceBCTypes_[fi * nodes.size() + nodeId];
        if (boundary == Boundary::REACTION) {
          detail::vecAddTo(sum, reactionBoundaryVelocity(node.index));
        } else if (boundary == Boundary::MASK) {
          detail::vecAddTo(sum,
                           maskVelocityBoundary(node.index, toT(v[nodeId])));
        } else {
          detail::vecAddTo(sum, toT(v[nodeId])); // AMBIENT/NONE: zero-flux
        }
      }
    }
  }

  // (Av)[i] = (2*D) * v[i] - sum_at_v[i] + b[i]
  template <class SolverT>
  void harmonicMatvec(const std::vector<Vec3D<SolverT>> &v,
                      const std::vector<Vec3D<T>> &b,
                      std::vector<Vec3D<SolverT>> &Av) const {
    const T diagVal = static_cast<T>(2 * D);
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      Vec3D<T> sum;
      computeHarmonicStencilAt(i, v, sum);
      for (unsigned c = 0; c < D; ++c)
        Av[i][c] = static_cast<SolverT>(diagVal * v[i][c] - sum[c] + b[i][c]);
    }
  }

  void solveVelocity() {
    iterations = 0;
    residual = 0.;
    if (nodes.empty())
      return;

    using SolverT = T;

    const std::size_t n = nodes.size();
    const T diagVal = static_cast<T>(2 * D); // constant for all nodes

    // b[i] = BC constants (reaction + mask velocities), computed at v = zeros.
    // OOB/NONE/AMBIENT faces contribute v[i] = 0 at zeros, so only Dirichlet
    // BCs survive — correctly isolating the RHS constant vector.
    std::vector<Vec3D<T>> b(n);
    {
      const std::vector<Vec3D<SolverT>> zeros(
          n, Vec3D<SolverT>{SolverT(0), SolverT(0), SolverT(0)});
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < n; ++i)
        computeHarmonicStencilAt(i, zeros, b[i]);
    }

    // Warm-start from previous substep's velocity field.
    std::vector<Vec3D<SolverT>> x(n);
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned c = 0; c < D; ++c) {
        const T value = nodes[i].velocity[c];
        x[i][c] = static_cast<SolverT>(std::isfinite(value) ? value : T(0));
      }

#ifdef VIENNALS_GPU_BICGSTAB
    if (gpu::gpuIsValid(gpuHarmonicBufs_)) {
      const std::size_t nf = 2u * D * n;
      if (harmonicDiagGpu_.size() != n || harmonicCoeffGpu_.size() != nf) {
        VIENNACORE_LOG_ERROR("OxidationDeformation: harmonic GPU geometry has "
                             "the wrong size for the current node set.");
      }

      Timer<> tUpload, tSolve;
      std::vector<Vec3D<SolverT>> xSolved(n);
      unsigned maxGpuIterations = 0;
      double maxGpuResidual = 0.0;

      for (unsigned c = 0; c < D; ++c) {
        std::vector<double> bGpu(n), xGpu(n);
        for (std::size_t i = 0; i < n; ++i) {
          bGpu[i] = static_cast<double>(b[i][c]);
          xGpu[i] = static_cast<double>(x[i][c]);
        }

        tUpload.start();
        const bool gpuUploadOk =
            (c == 0) ? gpu::gpuUploadSolverArrays(
                           gpuHarmonicBufs_, harmonicDiagGpu_.data(),
                           bGpu.data(), harmonicCoeffGpu_.data(),
                           static_cast<uint32_t>(n), harmonicCoeffGpu_.size())
                     : gpu::gpuUploadRhs(gpuHarmonicBufs_, bGpu.data(),
                                         static_cast<uint32_t>(n));
        tUpload.finish();
        if (!gpuUploadOk) {
          VIENNACORE_LOG_ERROR("OxidationDeformation: GPU mode was selected, "
                               "but uploading harmonic solver arrays failed." +
                               gpuErrorDetail());
        }

        unsigned gpuIterations = 0;
        double gpuResidual = 0.0;
        tSolve.start();
        const bool gpuConverged = gpu::gpuSolveBiCGSTAB(
            gpuHarmonicBufs_, xGpu.data(),
            static_cast<double>(std::numeric_limits<SolverT>::epsilon()),
            deformationParameters.harmonicIterations,
            static_cast<double>(deformationParameters.tolerance), gpuIterations,
            gpuResidual);
        tSolve.finish();

        if (!gpuConverged || !std::isfinite(gpuResidual)) {
          VIENNACORE_LOG_ERROR(
              "OxidationDeformation: harmonic GPU BiCGSTAB failed or produced "
              "a non-finite residual for component " +
              std::to_string(c) + " (iters=" + std::to_string(gpuIterations) +
              ", residual=" + std::to_string(gpuResidual) + ").");
        }

        maxGpuIterations = std::max(maxGpuIterations, gpuIterations);
        maxGpuResidual = std::max(maxGpuResidual, gpuResidual);
        for (std::size_t i = 0; i < n; ++i)
          xSolved[i][c] = static_cast<SolverT>(xGpu[i]);
      }

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          nodes[i].velocity[c] = static_cast<T>(xSolved[i][c]);
      iterations = maxGpuIterations;
      residual = maxGpuResidual;

      if (Logger::hasTiming()) {
        Logger::getInstance()
            .addTiming("harmonic n=" + std::to_string(n) +
                           " iters=" + std::to_string(iterations) + " res=" +
                           std::to_string(residual) + " [GPU] GPU BiCGSTAB",
                       tSolve)
            .print();
      }
      if (Logger::hasDebug()) {
        Logger::getInstance()
            .addTiming("harmonic n=" + std::to_string(n) + " [GPU] GPU upload",
                       tUpload)
            .print();
      }
      return;
    }
#endif

    // r = b - A*x
    const Vec3D<SolverT> zero3{SolverT(0), SolverT(0), SolverT(0)};
    std::vector<Vec3D<SolverT>> Ax(n);
    harmonicMatvec(x, b, Ax);
    std::vector<Vec3D<SolverT>> r(n), r_hat(n);
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned c = 0; c < D; ++c) {
        r[i][c] = static_cast<SolverT>(b[i][c] - Ax[i][c]);
        r_hat[i][c] = r[i][c];
      }

    // BiCGSTAB with diagonal preconditioner (diag = 2*D, constant).
    std::vector<Vec3D<SolverT>> pv(n, zero3), sv(n, zero3), y(n), z(n), s(n),
        t(n);
    T rho = T(1), alpha = T(1), omega = T(1);

    auto vecDot = [&](const std::vector<Vec3D<SolverT>> &a,
                      const std::vector<Vec3D<SolverT>> &bv) {
      T sum = T(0);
      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          sum += static_cast<T>(a[i][c]) * static_cast<T>(bv[i][c]);
      return sum;
    };

    auto vecMaxAbs = [&](const std::vector<Vec3D<SolverT>> &vin) {
      T m = T(0);
      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          m = std::max(m, std::abs(static_cast<T>(vin[i][c])));
      return m;
    };

    const T b_norm = [&] {
      T m = T(0);
      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          m = std::max(m, std::abs(b[i][c]));
      return (m < T(1e-100)) ? T(1) : m;
    }();

    for (; iterations < deformationParameters.harmonicIterations;
         ++iterations) {
      const T rho_new = vecDot(r_hat, r);
      if (!std::isfinite(rho_new) || std::abs(rho_new) < T(1e-100))
        break;
      if (!std::isfinite(rho) || !std::isfinite(alpha) ||
          !std::isfinite(omega) || std::abs(omega) < T(1e-100))
        break;

      const T beta = (rho_new / rho) * (alpha / omega);
      if (!std::isfinite(beta))
        break;
      rho = rho_new;

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          pv[i][c] = static_cast<SolverT>(r[i][c] +
                                          beta * (pv[i][c] - omega * sv[i][c]));

      // y = M^{-1} p = p / (2*D)
      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          y[i][c] = static_cast<SolverT>(static_cast<T>(pv[i][c]) / diagVal);

      harmonicMatvec(y, b, sv);

      const T r_hat_v = vecDot(r_hat, sv);
      if (!std::isfinite(r_hat_v) || std::abs(r_hat_v) < T(1e-100))
        break;

      alpha = rho_new / r_hat_v;
      if (!std::isfinite(alpha))
        break;

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          s[i][c] = static_cast<SolverT>(r[i][c] - alpha * sv[i][c]);

      residual = vecMaxAbs(s);
      if (!std::isfinite(residual))
        break;
      if (residual < deformationParameters.tolerance * b_norm) {
        for (std::size_t i = 0; i < n; ++i)
          for (unsigned c = 0; c < D; ++c)
            x[i][c] = static_cast<SolverT>(x[i][c] + alpha * y[i][c]);
        ++iterations;
        break;
      }

      // z = M^{-1} s
      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          z[i][c] = static_cast<SolverT>(static_cast<T>(s[i][c]) / diagVal);

      harmonicMatvec(z, b, t);

      const T t_s = vecDot(t, s);
      const T t_t = vecDot(t, t);
      if (!std::isfinite(t_s) || !std::isfinite(t_t))
        break;
      omega = (t_t > T(1e-100)) ? t_s / t_t : T(0);
      if (!std::isfinite(omega))
        break;

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c) {
          x[i][c] =
              static_cast<SolverT>(x[i][c] + alpha * y[i][c] + omega * z[i][c]);
          r[i][c] = static_cast<SolverT>(s[i][c] - omega * t[i][c]);
        }

      residual = vecMaxAbs(r);
      if (!std::isfinite(residual))
        break;
      if (residual < deformationParameters.tolerance * b_norm) {
        ++iterations;
        break;
      }
    }

    bool finiteSolution = true;
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned c = 0; c < D; ++c)
        if (!std::isfinite(static_cast<T>(x[i][c])))
          finiteSolution = false;

    if (finiteSolution) {
      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          nodes[i].velocity[c] = static_cast<T>(x[i][c]);
    } else {
      residual = std::numeric_limits<T>::infinity();
    }
    if (residual > deformationParameters.tolerance * b_norm)
      VIENNACORE_LOG_WARNING(
          "solveVelocity (harmonic): BiCGSTAB did not converge after " +
          std::to_string(iterations) + "/" +
          std::to_string(deformationParameters.harmonicIterations) +
          " iterations (residual=" + std::to_string(residual / b_norm) +
          ", tolerance=" + std::to_string(deformationParameters.tolerance) +
          ")");
  }

  // Returns component-wise diagonal entries of the Stokes operator A_v.
  // Geometry-fixed within a mechanics solve; computed once and reused by the
  // SIMPLE velocity-correction step: v_c^{k+1}=v_c* - grad_c(dp)/(eta*a_ic).
  //
  // With traction-coupled MASK contact, ghost=v_node+const for every component,
  // so the MASK face self-coupling cancels and the face coefficient is removed.
  std::vector<Vec3D<T>> computeVelocityDiagonals() const {
    const std::size_t n = nodes.size();
    std::vector<Vec3D<T>> diagV(n, Vec3D<T>{T(0), T(0), T(0)});
    if (n == 0)
      return diagV;
    const std::vector<Vec3D<T>> zeros(n, Vec3D<T>{T(0), T(0), T(0)});
    std::vector<Vec3D<T>> tmp(n);
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < n; ++i) {
      T diag{};
      computeVelocityStencilAt(i, zeros, diag, tmp[i]);
      for (unsigned comp = 0; comp < D; ++comp)
        diagV[i][comp] = diag;
    }

    return diagV;
  }

  void solveMechanics() {
    T mechanicsResidual = 0.;

    // SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) coupling.
    // The Gauss-Seidel p→v→p loop has spectral radius > 1 on thin geometries,
    // causing divergence that worsens with more iterations.  SIMPLE adds a
    // velocity-correction step after the pressure update that provably
    // eliminates the oscillation mode:
    //
    //   1. Momentum predictor:  A_v * v* = vBC - (∇p^k - ∇·σ'(v^k)) / η
    //   2. Pressure update:     A_p * p^{k+1} = pBC + K · div(v*)
    //   3. Velocity correction: v^{k+1} = v* - ∇δp / (η · a_i)
    //      where δp = p^{k+1} - p^k,  a_i = diag(A_v)[i]
    //
    // Step 3 ensures the corrected velocity is consistent with the new
    // pressure without re-solving the full momentum equation.  Unlike the
    // Aitken clamp (which can only damp, not stabilise, ρ > 1 iterations),
    // this correction is unconditionally convergent for steady Stokes.

    const std::vector<Vec3D<T>> diagV =
        computeVelocityDiagonals(); // geometry-fixed within this call

    for (unsigned iteration = 0;
         iteration < deformationParameters.mechanicsIterations; ++iteration) {
      const auto previousVelocity = collectVelocities(); // v^k
      const auto previousPressure = collectPressures();  // p^k

      computeDiagnostics();
      computeStressTensors();

      // Step 1: momentum predictor uses current p^k (in nodes[i].pressure).
      Timer<> tStokes, tPressure;
      tStokes.start();
      solveStokesVelocity(); // produces v* in nodes[i].velocity
      tStokes.finish();
      if (!std::isfinite(lastStokesResidual_)) {
        mechanicsResidual = std::numeric_limits<T>::infinity();
        break;
      }

      // Step 2: pressure solve uses divergence of v*.
      tPressure.start();
      solvePressure(); // produces p^{k+1} in nodes[i].pressure
      tPressure.finish();
      if (!std::isfinite(lastPressureResidual_)) {
        mechanicsResidual = std::numeric_limits<T>::infinity();
        break;
      }

      // Step 3: SIMPLE velocity correction: v^{k+1} = v* - ∇δp / (η · a_i).
      applySimpleVelocityCorrection(previousPressure, diagV);

      mechanicsResidual = std::max(maxVelocityChange(previousVelocity),
                                   maxPressureChange(previousPressure));
      if (!std::isfinite(mechanicsResidual)) {
        mechanicsResidual = std::numeric_limits<T>::infinity();
        break;
      }

      if (Logger::hasDebug())
        Logger::getInstance()
            .addTiming(
                "        mechanics[" + std::to_string(iteration) +
                    "] stokes   iters=" + std::to_string(lastStokesIters_) +
                    "/" +
                    std::to_string(deformationParameters.stokesIterations) +
                    " res=" + std::to_string(lastStokesResidual_),
                tStokes)
            .addTiming(
                "        mechanics[" + std::to_string(iteration) +
                    "] pressure iters=" + std::to_string(lastPressureIters_) +
                    "/" +
                    std::to_string(deformationParameters.pressureIterations) +
                    " res=" + std::to_string(lastPressureResidual_) +
                    " coupling=" + std::to_string(mechanicsResidual),
                tPressure)
            .print();

      if (mechanicsResidual < deformationParameters.mechanicsTolerance)
        break;
    }

    computeDiagnostics();
    computeStressTensors();
    residual = mechanicsResidual;
    if (residual > deformationParameters.mechanicsTolerance)
      VIENNACORE_LOG_WARNING(
          "solveMechanics: did not converge after " +
          std::to_string(deformationParameters.mechanicsIterations) +
          " iterations (residual=" + std::to_string(residual) + ", tolerance=" +
          std::to_string(deformationParameters.mechanicsTolerance) + ")");
  }

  // SIMPLE velocity correction: v^{k+1} = v* - ∇(p^{k+1} - p^k) / (η · a_i)
  //
  // δp gradient uses homogeneous Neumann at all boundary faces (δp ghost = 0).
  // The boundary pressure correction is re-enforced by the next pressure solve,
  // so this approximation only affects the current-iteration correction, not
  // the converged solution.
  void applySimpleVelocityCorrection(const std::vector<T> &pressureOld,
                                     const std::vector<Vec3D<T>> &diagV) {
    if (nodes.empty())
      return;
    if (deformationParameters.viscosity <= std::numeric_limits<T>::epsilon())
      return;

    const std::size_t n = nodes.size();
    const T invEta = T(1) / deformationParameters.viscosity;

#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < n; ++i) {
      for (unsigned dir = 0; dir < D; ++dir) {
        const T ai = diagV[i][dir];
        if (ai <= std::numeric_limits<T>::epsilon())
          continue;

        // Negative-offset face (fi = dir*2).
        T dpMinus, dMinus;
        {
          const unsigned fi = dir * 2u;
          IndexType nb = nodes[i].index;
          nb[dir] -= 1;
          if (inBounds(nb)) {
            const std::size_t j = nodeLookupFlat[linearIndex(nb)];
            if (j != noNode) {
              dpMinus = nodes[j].pressure - pressureOld[j];
              dMinus = gridDelta;
            } else {
              dpMinus = T(0);
              dMinus = faceBCDists_[fi * n + i];
            }
          } else {
            dpMinus = T(0);
            dMinus = gridDelta;
          }
        }

        // Positive-offset face (fi = dir*2+1).
        T dpPlus, dPlus;
        {
          const unsigned fi = dir * 2u + 1u;
          IndexType nb = nodes[i].index;
          nb[dir] += 1;
          if (inBounds(nb)) {
            const std::size_t j = nodeLookupFlat[linearIndex(nb)];
            if (j != noNode) {
              dpPlus = nodes[j].pressure - pressureOld[j];
              dPlus = gridDelta;
            } else {
              dpPlus = T(0);
              dPlus = faceBCDists_[fi * n + i];
            }
          } else {
            dpPlus = T(0);
            dPlus = gridDelta;
          }
        }

        const T dpCenter = nodes[i].pressure - pressureOld[i];
        const T gradDP =
            firstDerivative(dpMinus, dpCenter, dpPlus, dMinus, dPlus);
        const T correction =
            gradDP * invEta / ai * deformationParameters.relaxation;
        if (std::isfinite(correction) && std::isfinite(nodes[i].velocity[dir]))
          nodes[i].velocity[dir] -= correction;
      }
    }
  }

  // Fills diag = centerCoefficient and rhs = pressureSum for one node.
  // Dirichlet (ambient) nodes are encoded as identity rows: diag=1,
  // rhs=ambientBP.
  template <class SolverT>
  void computePressureStencilAt(std::size_t nodeId,
                                const std::vector<SolverT> &p,
                                const std::vector<T> &ambientBP, T &diag,
                                T &rhs) const {
    if (touchesAmbient_[nodeId]) {
      diag = T(1);
      rhs = ambientBP[nodeId];
      return;
    }
    diag = T(0);
    rhs = T(0);
    for (unsigned direction = 0; direction < D; ++direction) {
      const auto plus =
          pressureStencilPoint(p, ambientBP, nodeId, direction, 1);
      const auto minus =
          pressureStencilPoint(p, ambientBP, nodeId, direction, -1);
      const T dSum = plus.distance + minus.distance;
      const T plusCoeff = T(2) / (plus.distance * dSum);
      const T minusCoeff = T(2) / (minus.distance * dSum);
      rhs += plusCoeff * plus.value + minusCoeff * minus.value;
      diag += plusCoeff + minusCoeff;
    }
  }

  // (Av)[i] = precomputedDiag[i]*v[i] - rhs_at_v[i] + pBC[i]
  template <class SolverT>
  void
  pressureMatvec(const std::vector<SolverT> &v, const std::vector<T> &ambientBP,
                 const std::vector<T> &precomputedDiag,
                 const std::vector<T> &pBC, std::vector<SolverT> &Av) const {
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      T diag, rhs;
      computePressureStencilAt(i, v, ambientBP, diag, rhs);
      Av[i] = static_cast<SolverT>(precomputedDiag[i] * v[i] - rhs + pBC[i]);
    }
  }

  void solvePressure() {
    if (nodes.empty())
      return;

    using SolverT = T;

    const std::size_t n = nodes.size();
    const T eps = std::numeric_limits<T>::epsilon();

    std::vector<T> divergence(n), ambientBP(n);
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < n; ++i) {
      divergence[i] = divergenceAt(nodes[i].index);
      ambientBP[i] = freeSurfacePressureBoundary(nodes[i].index);
    }

    auto warnBadPressureAssembly = [](const std::string &stage,
                                      std::size_t nodeId, const IndexType &idx,
                                      T value) {
      VIENNACORE_LOG_WARNING(
          "solvePressure: non-finite/overflow " + stage +
          " at node=" + std::to_string(nodeId) + " index=(" +
          std::to_string(idx[0]) + "," + std::to_string(idx[1]) +
          (D == 3 ? "," + std::to_string(idx[2]) : std::string()) +
          ") value=" + std::to_string(value));
    };

    const T solverMax = static_cast<T>(std::numeric_limits<SolverT>::max());
    for (std::size_t i = 0; i < n; ++i) {
      if (!std::isfinite(divergence[i]) ||
          std::abs(divergence[i]) > solverMax) {
        warnBadPressureAssembly("divergence", i, nodes[i].index, divergence[i]);
        break;
      }
      if (!std::isfinite(ambientBP[i]) || std::abs(ambientBP[i]) > solverMax) {
        warnBadPressureAssembly("ambient pressure boundary", i, nodes[i].index,
                                ambientBP[i]);
        break;
      }
    }

    // Geometry-fixed diagonal and BC constants (kept in T for full precision).
    std::vector<T> diag(n), pBC(n);
    {
      const std::vector<SolverT> zeros(n, SolverT(0));
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < n; ++i)
        computePressureStencilAt(i, zeros, ambientBP, diag[i], pBC[i]);
    }

    for (std::size_t i = 0; i < n; ++i) {
      if (!std::isfinite(diag[i]) || std::abs(diag[i]) > solverMax) {
        warnBadPressureAssembly("pressure diagonal", i, nodes[i].index,
                                diag[i]);
        break;
      }
      if (!std::isfinite(pBC[i]) || std::abs(pBC[i]) > solverMax) {
        warnBadPressureAssembly("pressure boundary rhs", i, nodes[i].index,
                                pBC[i]);
        break;
      }
    }

    std::vector<T> b(n);
    T b_norm = T(0);
    for (std::size_t i = 0; i < n; ++i) {
      b[i] = pBC[i] + deformationParameters.bulkModulus * divergence[i];
      b_norm = std::max(b_norm, std::abs(b[i]));
    }
    for (std::size_t i = 0; i < n; ++i) {
      if (!std::isfinite(b[i]) || std::abs(b[i]) > solverMax) {
        warnBadPressureAssembly("pressure rhs", i, nodes[i].index, b[i]);
        break;
      }
    }
    if (b_norm < T(1e-100))
      b_norm = T(1);

    std::vector<SolverT> x(n);
    for (std::size_t i = 0; i < n; ++i) {
      T guess = touchesAmbient_[i] ? ambientBP[i] : nodes[i].pressure;
      if (!std::isfinite(guess))
        guess = deformationParameters.ambientPressure;
      x[i] = static_cast<SolverT>(guess);
    }

#ifdef VIENNALS_GPU_BICGSTAB
    if (gpu::gpuIsValid(gpuPressBufs_)) {
      const std::size_t nf = 2u * D * n;
      if (actualDiagGpu_.size() != n || pressCoeffGpu_.size() != nf) {
        VIENNACORE_LOG_ERROR("OxidationDeformation: pressure GPU geometry has "
                             "the wrong size for the current node set.");
      }

      Timer<> tUpload, tSolve;
      std::vector<double> bGpu(n), xGpu(n);
      for (std::size_t i = 0; i < n; ++i) {
        bGpu[i] = static_cast<double>(b[i]);
        xGpu[i] = static_cast<double>(x[i]);
      }

      tUpload.start();
      const bool gpuUploadOk = gpu::gpuUploadSolverArrays(
          gpuPressBufs_, actualDiagGpu_.data(), bGpu.data(),
          pressCoeffGpu_.data(), static_cast<uint32_t>(n),
          pressCoeffGpu_.size());
      tUpload.finish();
      if (!gpuUploadOk) {
        VIENNACORE_LOG_ERROR("OxidationDeformation: GPU mode was selected, but "
                             "uploading pressure solver arrays or factorizing "
                             "ILU failed." +
                             gpuErrorDetail());
      }

      unsigned gpuIterations = 0;
      double gpuResidual = 0.0;
      tSolve.start();
      const bool gpuConverged = gpu::gpuSolveBiCGSTAB(
          gpuPressBufs_, xGpu.data(), static_cast<double>(eps),
          deformationParameters.pressureIterations,
          static_cast<double>(deformationParameters.pressureTolerance),
          gpuIterations, gpuResidual);
      tSolve.finish();

      // gpuResidual is the GPU true residual ||b - A*x||_inf recomputed at
      // convergence (not the recursive BiCGSTAB residual), so no separate CPU
      // stencil evaluation is needed.
      if (!gpuConverged || !std::isfinite(gpuResidual)) {
        VIENNACORE_LOG_ERROR(
            "OxidationDeformation: pressure GPU BiCGSTAB failed or produced "
            "a non-finite residual (iters=" +
            std::to_string(gpuIterations) +
            ", residual=" + std::to_string(gpuResidual) + ").");
      }

      {
        const T beta = deformationParameters.pressureRelaxation;
        const T oneMinB = T(1) - beta;
        for (std::size_t i = 0; i < n; ++i)
          nodes[i].pressure =
              oneMinB * nodes[i].pressure + beta * static_cast<T>(xGpu[i]);
      }
      lastPressureIters_ = gpuIterations;
      lastPressureResidual_ = gpuResidual / b_norm;

      if (Logger::hasDebug()) {
        const std::string tag =
            "pressure n=" + std::to_string(n) +
            " iters=" + std::to_string(lastPressureIters_) +
            " res=" + std::to_string(lastPressureResidual_) + " [GPU]";
        Logger::getInstance()
            .addTiming(tag + " GPU upload", tUpload)
            .addTiming(tag + " GPU BiCGSTAB", tSolve)
            .print();
      }
      return;
    }
#endif

    // Precompute off-diagonal structure and the CORRECT matrix diagonal for
    // SSOR.
    //
    // Key insight: diag[i] from computePressureStencilAt includes self-coupling
    // contributions from REACTION/MASK/OOB faces (those return v[nodeId]
    // itself). The ACTUAL matrix diagonal A[i,i] = sum of off-diagonal
    // (interior-neighbor) coefficients only.  Using the wrong diagonal in the
    // SSOR sweeps makes the preconditioner invalid near boundaries.
    //
    // Also: NONE-type non-interior faces (OOB or no crossing) use gridDelta in
    // pressureStencilPoint, NOT faceBCDists_ (which defaults to T(1)).
    //
    // Face-major layout: fi = dir*2 + (offset==+1 ? 1 : 0)
    //   Even fi (offset=-1): lower-index neighbor → forward sweep
    //   Odd  fi (offset=+1): higher-index neighbor → backward sweep
    std::vector<T> pressCoeff(2 * D * n, T(0));
    std::vector<std::size_t> pressNeighId(2 * D * n, noNode);
    std::vector<T> actualDiag(n, T(0)); // A[i,i] = sum of interior coefficients

    for (std::size_t id = 0; id < n; ++id) {
      if (touchesAmbient_[id]) {
        actualDiag[id] = T(1);
        continue;
      } // identity row
      for (unsigned dir = 0; dir < D; ++dir) {
        const unsigned fiNeg = dir * 2u;
        const unsigned fiPos = dir * 2u + 1u;
        IndexType nbNeg = nodes[id].index;
        nbNeg[dir] -= 1;
        IndexType nbPos = nodes[id].index;
        nbPos[dir] += 1;
        const std::size_t jNeg =
            inBounds(nbNeg) ? nodeLookupFlat[linearIndex(nbNeg)] : noNode;
        const std::size_t jPos =
            inBounds(nbPos) ? nodeLookupFlat[linearIndex(nbPos)] : noNode;

        // Effective distance matching pressureStencilPoint:
        //   interior neighbour  → gridDelta
        //   AMBIENT/REACTION/MASK crossing → faceBCDists_ (actual sub-grid
        //   distance) NONE (OOB or no crossing) → gridDelta
        //   (pressureStencilPoint fallthrough)
        auto effectiveDist = [&](unsigned fi, std::size_t j) -> T {
          if (j != noNode)
            return gridDelta;
          const Boundary bt = faceBCTypes_[fi * n + id];
          if (bt != Boundary::NONE)
            return faceBCDists_[fi * n + id];
          return gridDelta;
        };

        const T dNeg = effectiveDist(fiNeg, jNeg);
        const T dPos = effectiveDist(fiPos, jPos);
        const T dSum = dNeg + dPos;
        if (dSum <= eps)
          continue;

        if (jNeg != noNode && !touchesAmbient_[jNeg]) {
          const T c = T(2) / (dNeg * dSum);
          pressCoeff[fiNeg * n + id] = c;
          pressNeighId[fiNeg * n + id] = jNeg;
          actualDiag[id] += c; // A[i,i] += interior off-diagonal coefficient
        } else if (jNeg != noNode ||
                   faceBCTypes_[fiNeg * n + id] == Boundary::AMBIENT) {
          // j is an ambient-only neighbour (identity-row Dirichlet p=0), OR
          // this face directly crosses the free surface (AMBIENT Dirichlet p=0
          // at the sub-grid crossing distance).  RHS contribution is c·0=0.
          // REACTION faces are solid-wall Neumann ∂p/∂n=0: no contribution.
          actualDiag[id] += T(2) / (dNeg * dSum);
        }
        if (jPos != noNode && !touchesAmbient_[jPos]) {
          const T c = T(2) / (dPos * dSum);
          pressCoeff[fiPos * n + id] = c;
          pressNeighId[fiPos * n + id] = jPos;
          actualDiag[id] += c;
        } else if (jPos != noNode ||
                   faceBCTypes_[fiPos * n + id] == Boundary::AMBIENT) {
          actualDiag[id] += T(2) / (dPos * dSum);
        }
      }
      // Guard against fully-isolated nodes (surrounded by boundaries on every
      // face)
      if (actualDiag[id] <= eps)
        actualDiag[id] = T(1);
    }

    // ILU(0) preconditioner for the (non-symmetric) pressure matrix.
    //
    // The sub-grid interface distances make A[i,j] ≠ A[j,i] in general, so
    // SSOR is not guaranteed to converge.  ILU(0) handles non-symmetric
    // matrices robustly.
    //
    // Factorisation A ≈ L * U (zero fill-in, natural node ordering):
    //   L  – unit lower triangular: L[i,j] = A[i,j] / U[j,j]  for j < i
    //   U  – upper triangular:      U[i,j] = A[i,j]            for j > i
    //   U[i,i] = A[i,i] - Σ_{k<i} L[i,k] * A[k,i]
    //
    // With A[i,j] = -pressCoeff[fi*n+i] and A[j,i] = -pressCoeff[(fi^1)*n+j]:
    //   U[i,i] = actualDiag[i] - Σ_{lower j} pressCoeff[fi_L*n+i]
    //                                         * pressCoeff[fi_U*n+j]
    //                                         / ilu_diag[j]
    //
    // Preconditioner application M_ILU^{-1} r = z:
    //   Forward (L y = r, unit lower triangular, no diagonal divide):
    //     y[i] = r[i] + Σ_{j<i} (pressCoeff[fi_L*n+i] / ilu_diag[j]) * y[j]
    //   Backward (U z = y):
    //     z[i] = (y[i] + Σ_{j>i} pressCoeff[fi_U*n+i] * z[j]) / ilu_diag[i]
    std::vector<T> ilu_diag(n);
    for (std::size_t id = 0; id < n; ++id) {
      if (touchesAmbient_[id]) {
        ilu_diag[id] = T(1);
        continue;
      }
      ilu_diag[id] = actualDiag[id];
      for (unsigned dir = 0; dir < D; ++dir) {
        const unsigned fi_L = dir * 2u; // lower face (offset=-1)
        const unsigned fi_U =
            fi_L + 1u; // upper face (offset=+1, j's face toward i)
        const std::size_t j = pressNeighId[fi_L * n + id];
        if (j == noNode || ilu_diag[j] <= eps)
          continue;
        // L[id,j] = A[id,j] / U[j,j] = (-pressCoeff_L) / ilu_diag[j]
        // A[j,id] = -pressCoeff[fi_U * n + j]   (j's upper-face coefficient
        // toward id) ilu_diag[id] -= L[id,j] * A[j,id]
        //              = (-pressCoeff_L / ilu_diag[j]) * (-pressCoeff_fi_U[j])
        //              = pressCoeff_L * pressCoeff_fi_U[j] / ilu_diag[j]
        //              (positive drop)
        ilu_diag[id] -=
            pressCoeff[fi_L * n + id] * pressCoeff[fi_U * n + j] / ilu_diag[j];
      }
      if (ilu_diag[id] <= eps)
        ilu_diag[id] = actualDiag[id]; // guard non-positive pivot
    }

    auto applyIlu = [&](const std::vector<SolverT> &in,
                        std::vector<SolverT> &out) {
      std::vector<T> y(n);
      // Forward solve: L * y = in  (L is unit lower triangular)
      for (std::size_t i = 0; i < n; ++i) {
        T val = static_cast<T>(in[i]);
        for (unsigned dir = 0; dir < D; ++dir) {
          const unsigned fi_L = dir * 2u;
          const std::size_t j = pressNeighId[fi_L * n + i];
          if (j != noNode)
            // L[i,j] = -pressCoeff[fi_L*n+i] / ilu_diag[j], subtract
            // A[i,j]*y[j]: y[i] -= L[i,j] * y[j] = -(-pressCoeff/ilu_diag[j]) *
            // y[j] = +(coeff/ilu) * y[j]
            val += (pressCoeff[fi_L * n + i] / ilu_diag[j]) * y[j];
        }
        y[i] = val; // no diagonal divide (unit lower triangular)
      }
      // Backward solve: U * z = y
      for (std::size_t i = n; i-- > 0;) {
        T val = y[i];
        for (unsigned dir = 0; dir < D; ++dir) {
          const unsigned fi_U = dir * 2u + 1u;
          const std::size_t j = pressNeighId[fi_U * n + i];
          if (j != noNode)
            // U[i,j] = -pressCoeff[fi_U*n+i], subtract U[i,j]*z[j]:
            // val -= U[i,j] * z[j] = -(-pressCoeff) * z[j] = +(pressCoeff) *
            // z[j]
            val += pressCoeff[fi_U * n + i] * static_cast<T>(out[j]);
        }
        out[i] = static_cast<SolverT>(val / ilu_diag[i]);
      }
    };

    std::vector<SolverT> Ax(n);
    pressureMatvec(x, ambientBP, diag, pBC, Ax);
    for (std::size_t i = 0; i < n; ++i) {
      if (!std::isfinite(static_cast<T>(Ax[i])) ||
          std::abs(static_cast<T>(Ax[i])) > solverMax) {
        warnBadPressureAssembly("initial pressure matvec", i, nodes[i].index,
                                static_cast<T>(Ax[i]));
        break;
      }
    }
    std::vector<SolverT> r(n), r_hat(n), p(n, SolverT(0)), v(n, SolverT(0)),
        y(n), z(n), s(n), t(n);
    for (std::size_t i = 0; i < n; ++i) {
      r[i] = static_cast<SolverT>(b[i] - Ax[i]);
      r_hat[i] = r[i];
    }

    T rho = T(1), alpha = T(1), omega = T(1);
    T pressureResidual = T(0);
    unsigned pressureIter = 0;
    bool pressureBreakdown = false;
    for (std::size_t i = 0; i < n; ++i) {
      const T ri = static_cast<T>(r[i]);
      if (!std::isfinite(ri)) {
        pressureBreakdown = true;
        break;
      }
      pressureResidual = std::max(pressureResidual, std::abs(ri));
    }

    for (; !pressureBreakdown &&
           pressureIter < deformationParameters.pressureIterations;
         ++pressureIter) {
      T rho_new = T(0);
      for (std::size_t i = 0; i < n; ++i)
        rho_new += static_cast<T>(r_hat[i]) * static_cast<T>(r[i]);

      if (!std::isfinite(rho_new)) {
        pressureBreakdown = true;
        break;
      }
      if (std::abs(rho_new) < T(1e-100))
        break;
      if (!std::isfinite(rho) || !std::isfinite(alpha) ||
          !std::isfinite(omega) || std::abs(omega) < T(1e-100)) {
        pressureBreakdown = true;
        break;
      }

      const T beta = (rho_new / rho) * (alpha / omega);
      if (!std::isfinite(beta)) {
        pressureBreakdown = true;
        break;
      }
      rho = rho_new;

      for (std::size_t i = 0; i < n; ++i)
        p[i] = static_cast<SolverT>(r[i] + beta * (p[i] - omega * v[i]));

      applyIlu(p, y);

      pressureMatvec(y, ambientBP, diag, pBC, v);

      T r_hat_v = T(0);
      for (std::size_t i = 0; i < n; ++i)
        r_hat_v += static_cast<T>(r_hat[i]) * static_cast<T>(v[i]);
      if (!std::isfinite(r_hat_v)) {
        pressureBreakdown = true;
        break;
      }
      if (std::abs(r_hat_v) < T(1e-100))
        break;

      alpha = rho_new / r_hat_v;
      if (!std::isfinite(alpha)) {
        pressureBreakdown = true;
        break;
      }

      for (std::size_t i = 0; i < n; ++i)
        s[i] = static_cast<SolverT>(r[i] - alpha * v[i]);

      pressureResidual = T(0);
      for (std::size_t i = 0; i < n; ++i)
        pressureResidual =
            std::max(pressureResidual, std::abs(static_cast<T>(s[i])));
      if (!std::isfinite(pressureResidual)) {
        pressureBreakdown = true;
        break;
      }
      if (pressureResidual < deformationParameters.pressureTolerance * b_norm) {
        for (std::size_t i = 0; i < n; ++i)
          x[i] = static_cast<SolverT>(x[i] + alpha * y[i]);
        break;
      }

      applyIlu(s, z);

      pressureMatvec(z, ambientBP, diag, pBC, t);

      T t_s = T(0), t_t = T(0);
      for (std::size_t i = 0; i < n; ++i) {
        t_s += static_cast<T>(t[i]) * static_cast<T>(s[i]);
        t_t += static_cast<T>(t[i]) * static_cast<T>(t[i]);
      }
      if (!std::isfinite(t_s) || !std::isfinite(t_t)) {
        pressureBreakdown = true;
        break;
      }
      omega = (t_t > T(1e-100)) ? t_s / t_t : T(0);
      if (!std::isfinite(omega)) {
        pressureBreakdown = true;
        break;
      }

      for (std::size_t i = 0; i < n; ++i) {
        x[i] = static_cast<SolverT>(x[i] + alpha * y[i] + omega * z[i]);
        r[i] = static_cast<SolverT>(s[i] - omega * t[i]);
      }

      pressureResidual = T(0);
      for (std::size_t i = 0; i < n; ++i)
        pressureResidual =
            std::max(pressureResidual, std::abs(static_cast<T>(r[i])));
      if (!std::isfinite(pressureResidual)) {
        pressureBreakdown = true;
        break;
      }
      if (pressureResidual < deformationParameters.pressureTolerance * b_norm)
        break;
    }

    if (pressureBreakdown)
      pressureResidual = std::numeric_limits<T>::infinity();

    bool finiteSolution = !pressureBreakdown;
    for (std::size_t i = 0; i < n; ++i)
      if (!std::isfinite(static_cast<T>(x[i])))
        finiteSolution = false;

    lastPressureIters_ = pressureIter;
    lastPressureResidual_ = pressureResidual / b_norm;
    if (finiteSolution) {
      const T beta = deformationParameters.pressureRelaxation;
      const T oneMinB = T(1) - beta;
      for (std::size_t i = 0; i < n; ++i)
        nodes[i].pressure =
            oneMinB * nodes[i].pressure + beta * static_cast<T>(x[i]);
    } else {
      lastPressureResidual_ = std::numeric_limits<T>::infinity();
    }
    if (lastPressureResidual_ > deformationParameters.pressureTolerance)
      VIENNACORE_LOG_WARNING(
          "solvePressure: BiCGSTAB did not converge after " +
          std::to_string(lastPressureIters_) + "/" +
          std::to_string(deformationParameters.pressureIterations) +
          " iterations (residual=" + std::to_string(lastPressureResidual_) +
          ", tolerance=" +
          std::to_string(deformationParameters.pressureTolerance) + ")");
  }

  // Fills scalar diag = centerCoefficient and Vec3D rhs = velocitySum for one
  // node.
  template <class SolverT>
  void computeVelocityStencilAt(std::size_t nodeId,
                                const std::vector<Vec3D<SolverT>> &v, T &diag,
                                Vec3D<T> &rhs) const {
    diag = T(0);
    rhs = {T(0), T(0), T(0)};
    for (unsigned direction = 0; direction < D; ++direction) {
      const auto plus = velocityStencilPoint(v, nodeId, direction, 1);
      const auto minus = velocityStencilPoint(v, nodeId, direction, -1);
      const T dSum = plus.distance + minus.distance;
      const T plusCoeff = T(2) / (plus.distance * dSum);
      const T minusCoeff = T(2) / (minus.distance * dSum);
      detail::vecAddTo(rhs, detail::vecScaled(plus.value, plusCoeff));
      detail::vecAddTo(rhs, detail::vecScaled(minus.value, minusCoeff));
      diag += plusCoeff + minusCoeff;
    }
  }

  void solveStokesVelocity() {
    if (deformationParameters.viscosity <= std::numeric_limits<T>::epsilon())
      return;
    if (nodes.empty())
      return;

    using SolverT = T;

    const std::size_t n = nodes.size();
    const T eps = std::numeric_limits<T>::epsilon();

    // Geometry-fixed diagonal, BC constants, and forcing (all in T).
    std::vector<T> diag(n);
    std::vector<Vec3D<T>> vBC(n), forcing(n);
    {
      const std::vector<Vec3D<SolverT>> zeros(
          n, Vec3D<SolverT>{SolverT(0), SolverT(0), SolverT(0)});
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < n; ++i) {
        computeVelocityStencilAt(i, zeros, diag[i], vBC[i]);
        forcing[i] = momentumForcing(nodes[i].index);
      }
    }
    const auto precondDiag = computeVelocityDiagonals();

    std::vector<Vec3D<T>> b(n);
    T b_norm = T(0);
    for (std::size_t i = 0; i < n; ++i) {
      for (unsigned c = 0; c < D; ++c) {
        b[i][c] = vBC[i][c] - forcing[i][c] / deformationParameters.viscosity;
        b_norm = std::max(b_norm, std::abs(b[i][c]));
      }
    }
    if (b_norm < T(1e-100))
      b_norm = T(1);

    // Initial guess from current node velocities (warm-start), converted to
    // SolverT.
    std::vector<Vec3D<SolverT>> x(n);
    {
      const auto vel = collectVelocities();
      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c) {
          const T value = vel[i][c];
          x[i][c] = static_cast<SolverT>(std::isfinite(value) ? value : T(0));
        }
    }

#ifdef VIENNALS_GPU_BICGSTAB
    if (gpu::gpuIsValid(gpuStokesBufs_)) {
      const std::size_t nf = 2u * D * n;
      if (stokesDiagGpu_.size() != D * n || stokesCoeffGpu_.size() != nf) {
        VIENNACORE_LOG_ERROR("OxidationDeformation: Stokes GPU geometry has "
                             "the wrong size for the current node set.");
      }

      Timer<> tUpload, tSolve;
      std::vector<Vec3D<SolverT>> xSolved(n);
      unsigned maxGpuIterations = 0;
      double maxGpuResidual = 0.0;

      for (unsigned c = 0; c < D; ++c) {
        std::vector<double> bGpu(n), xGpu(n);
        for (std::size_t i = 0; i < n; ++i) {
          bGpu[i] = static_cast<double>(b[i][c]);
          xGpu[i] = static_cast<double>(x[i][c]);
        }

        tUpload.start();
        const bool gpuUploadOk = gpu::gpuUploadSolverArrays(
            gpuStokesBufs_, stokesDiagGpu_.data() + c * n, bGpu.data(),
            stokesCoeffGpu_.data(), static_cast<uint32_t>(n),
            stokesCoeffGpu_.size());
        tUpload.finish();
        if (!gpuUploadOk) {
          VIENNACORE_LOG_ERROR("OxidationDeformation: GPU mode was selected, "
                               "but uploading Stokes solver arrays failed." +
                               gpuErrorDetail());
        }

        unsigned gpuIterations = 0;
        double gpuResidual = 0.0;
        tSolve.start();
        const bool gpuConverged = gpu::gpuSolveBiCGSTAB(
            gpuStokesBufs_, xGpu.data(), static_cast<double>(eps),
            deformationParameters.stokesIterations,
            static_cast<double>(deformationParameters.stokesTolerance),
            gpuIterations, gpuResidual);
        tSolve.finish();

        if (!gpuConverged || !std::isfinite(gpuResidual)) {
          VIENNACORE_LOG_ERROR(
              "OxidationDeformation: Stokes GPU BiCGSTAB failed or produced "
              "a non-finite residual for component " +
              std::to_string(c) + " (iters=" + std::to_string(gpuIterations) +
              ", residual=" + std::to_string(gpuResidual) + ").");
        }

        maxGpuIterations = std::max(maxGpuIterations, gpuIterations);
        maxGpuResidual = std::max(maxGpuResidual, gpuResidual);
        for (std::size_t i = 0; i < n; ++i)
          xSolved[i][c] = static_cast<SolverT>(xGpu[i]);
      }

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          nodes[i].velocity[c] = static_cast<T>(xSolved[i][c]);

      lastStokesIters_ = maxGpuIterations;
      lastStokesResidual_ = maxGpuResidual / b_norm;

      if (Logger::hasDebug()) {
        const std::string tag = "stokes n=" + std::to_string(n) +
                                " iters=" + std::to_string(lastStokesIters_) +
                                " res=" + std::to_string(lastStokesResidual_) +
                                " [GPU]";
        Logger::getInstance()
            .addTiming(tag + " GPU upload", tUpload)
            .addTiming(tag + " GPU BiCGSTAB", tSolve)
            .print();
      }
      return;
    }
#endif

    // Stokes SpMV: (Av)[i] = diag[i]*vin[i] - rhs_at_vin[i] + vBC[i], stored as
    // SolverT.
    auto stokesMatvec = [&](const std::vector<Vec3D<SolverT>> &vin,
                            std::vector<Vec3D<SolverT>> &Av) {
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < n; ++i) {
        T d;
        Vec3D<T> rhs;
        computeVelocityStencilAt(i, vin, d, rhs);
        for (unsigned c = 0; c < D; ++c)
          Av[i][c] =
              static_cast<SolverT>(diag[i] * vin[i][c] - rhs[c] + vBC[i][c]);
      }
    };

    // Dot product accumulated in T for numerical stability.
    auto vecDot = [&](const std::vector<Vec3D<SolverT>> &a,
                      const std::vector<Vec3D<SolverT>> &bv) {
      T sum = T(0);
      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c) {
          const T av = static_cast<T>(a[i][c]);
          const T bvVal = static_cast<T>(bv[i][c]);
          if (!std::isfinite(av) || !std::isfinite(bvVal))
            return std::numeric_limits<T>::quiet_NaN();
          sum += av * bvVal;
        }
      return sum;
    };

    auto vecMaxAbs = [&](const std::vector<Vec3D<SolverT>> &vin) {
      T m = T(0);
      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c) {
          const T value = static_cast<T>(vin[i][c]);
          if (!std::isfinite(value))
            return std::numeric_limits<T>::infinity();
          m = std::max(m, std::abs(value));
        }
      return m;
    };

    // r = b - A*x
    const Vec3D<SolverT> zero3{SolverT(0), SolverT(0), SolverT(0)};
    std::vector<Vec3D<SolverT>> Ax(n), r(n), r_hat(n);
    std::vector<Vec3D<SolverT>> pv(n, zero3), sv(n, zero3), y(n), z(n), s(n),
        t(n);
    stokesMatvec(x, Ax);
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned c = 0; c < D; ++c) {
        r[i][c] = static_cast<SolverT>(b[i][c] - Ax[i][c]);
        r_hat[i][c] = r[i][c];
      }

    T rho = T(1), alpha = T(1), omega = T(1);
    T velocityResidual = T(0);
    unsigned stokesIter = 0;
    bool stokesBreakdown = false;
    velocityResidual = vecMaxAbs(r);

    for (; stokesIter < deformationParameters.stokesIterations; ++stokesIter) {
      const T rho_new = vecDot(r_hat, r);
      if (!std::isfinite(rho_new)) {
        stokesBreakdown = true;
        break;
      }
      if (std::abs(rho_new) < T(1e-100))
        break;
      if (!std::isfinite(rho) || !std::isfinite(alpha) ||
          !std::isfinite(omega) || std::abs(omega) < T(1e-100)) {
        stokesBreakdown = true;
        break;
      }

      const T beta = (rho_new / rho) * (alpha / omega);
      if (!std::isfinite(beta)) {
        stokesBreakdown = true;
        break;
      }
      rho = rho_new;

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          pv[i][c] = static_cast<SolverT>(r[i][c] +
                                          beta * (pv[i][c] - omega * sv[i][c]));

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c) {
          const T pvc = pv[i][c];
          const T pcDiag = precondDiag[i][c];
          y[i][c] = static_cast<SolverT>((pcDiag > eps) ? pvc / pcDiag : pvc);
        }

      stokesMatvec(y, sv);

      const T r_hat_v = vecDot(r_hat, sv);
      if (!std::isfinite(r_hat_v)) {
        stokesBreakdown = true;
        break;
      }
      if (std::abs(r_hat_v) < T(1e-100))
        break;

      alpha = rho_new / r_hat_v;
      if (!std::isfinite(alpha)) {
        stokesBreakdown = true;
        break;
      }

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          s[i][c] = static_cast<SolverT>(r[i][c] - alpha * sv[i][c]);

      velocityResidual = vecMaxAbs(s);
      if (!std::isfinite(velocityResidual)) {
        stokesBreakdown = true;
        break;
      }
      if (velocityResidual < deformationParameters.stokesTolerance * b_norm) {
        for (std::size_t i = 0; i < n; ++i)
          for (unsigned c = 0; c < D; ++c)
            x[i][c] = static_cast<SolverT>(x[i][c] + alpha * y[i][c]);
        break;
      }

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c) {
          const T sc = s[i][c];
          const T pcDiag = precondDiag[i][c];
          z[i][c] = static_cast<SolverT>((pcDiag > eps) ? sc / pcDiag : sc);
        }

      stokesMatvec(z, t);

      const T t_s = vecDot(t, s);
      const T t_t = vecDot(t, t);
      if (!std::isfinite(t_s) || !std::isfinite(t_t)) {
        stokesBreakdown = true;
        break;
      }
      omega = (t_t > T(1e-100)) ? t_s / t_t : T(0);
      if (!std::isfinite(omega)) {
        stokesBreakdown = true;
        break;
      }

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c) {
          x[i][c] =
              static_cast<SolverT>(x[i][c] + alpha * y[i][c] + omega * z[i][c]);
          r[i][c] = static_cast<SolverT>(s[i][c] - omega * t[i][c]);
        }

      velocityResidual = vecMaxAbs(r);
      if (!std::isfinite(velocityResidual)) {
        stokesBreakdown = true;
        break;
      }
      if (velocityResidual < deformationParameters.stokesTolerance * b_norm)
        break;
    }

    if (stokesBreakdown)
      velocityResidual = std::numeric_limits<T>::infinity();

    bool finiteSolution = !stokesBreakdown;
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned c = 0; c < D; ++c)
        if (!std::isfinite(static_cast<T>(x[i][c])))
          finiteSolution = false;

    if (finiteSolution) {
      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          nodes[i].velocity[c] = static_cast<T>(x[i][c]);
    }

    lastStokesIters_ = stokesIter;
    lastStokesResidual_ = finiteSolution ? velocityResidual / b_norm
                                         : std::numeric_limits<T>::infinity();
    if (lastStokesResidual_ > deformationParameters.stokesTolerance)
      VIENNACORE_LOG_WARNING(
          "solveStokesVelocity: BiCGSTAB did not converge after " +
          std::to_string(lastStokesIters_) + "/" +
          std::to_string(deformationParameters.stokesIterations) +
          " iterations (residual=" + std::to_string(lastStokesResidual_) +
          ", tolerance=" +
          std::to_string(deformationParameters.stokesTolerance) + ")");
  }

  std::vector<Vec3D<T>> collectVelocities() const {
    std::vector<Vec3D<T>> velocities;
    velocities.reserve(nodes.size());
    for (const auto &node : nodes)
      velocities.push_back(node.velocity);
    return velocities;
  }

  std::vector<T> collectPressures() const {
    std::vector<T> pressures;
    pressures.reserve(nodes.size());
    for (const auto &node : nodes)
      pressures.push_back(node.pressure);
    return pressures;
  }

  template <class SolverT>
  StencilPoint<T>
  pressureStencilPoint(const std::vector<SolverT> &pressure,
                       const std::vector<T> &ambientBoundaryPressure,
                       std::size_t nodeId, unsigned direction,
                       int offset) const {
    const auto &node = nodes[nodeId];
    IndexType neighbor = node.index;
    neighbor[direction] += offset;

    if (!inBounds(neighbor))
      return {static_cast<T>(pressure[nodeId]), gridDelta};

    const std::size_t neighborId = nodeLookupFlat[linearIndex(neighbor)];
    if (neighborId != noNode) {
      if (touchesAmbient_[neighborId])
        return {ambientBoundaryPressure[neighborId], gridDelta};
      return {static_cast<T>(pressure[neighborId]), gridDelta};
    }

    const unsigned fi = direction * 2u + (offset == 1 ? 1u : 0u);
    const std::size_t nn = nodes.size();
    const Boundary faceType = faceBCTypes_[fi * nn + nodeId];
    const T faceDist = faceBCDists_[fi * nn + nodeId];
    if (faceType == Boundary::AMBIENT)
      return {ambientBoundaryPressure[nodeId], faceDist};
    // Reaction interface: Neumann ∂p/∂n=0 (solid-wall BC for pressure).
    // The ghost node takes the same value as the interior node, giving zero
    // contribution to the Laplacian stencil.  The pressure is anchored only by
    // the AMBIENT Dirichlet (p=0 at the free surface), which is always present
    // for any connected oxide region.
    if (faceType == Boundary::REACTION)
      return {static_cast<T>(pressure[nodeId]), faceDist};
    if (faceType == Boundary::MASK)
      return {maskPressureBoundary(node.index, direction, offset,
                                   static_cast<T>(pressure[nodeId])),
              faceDist};

    return {static_cast<T>(pressure[nodeId]), gridDelta};
  }

  template <class SolverT>
  StencilPoint<Vec3D<T>>
  velocityStencilPoint(const std::vector<Vec3D<SolverT>> &velocity,
                       std::size_t nodeId, unsigned direction,
                       int offset) const {
    const auto &node = nodes[nodeId];
    IndexType neighbor = node.index;
    neighbor[direction] += offset;

    const auto toT = [](const Vec3D<SolverT> &v) -> Vec3D<T> {
      return {static_cast<T>(v[0]), static_cast<T>(v[1]), static_cast<T>(v[2])};
    };

    if (!inBounds(neighbor))
      return {toT(velocity[nodeId]), gridDelta};

    const std::size_t neighborId = nodeLookupFlat[linearIndex(neighbor)];
    if (neighborId != noNode)
      return {toT(velocity[neighborId]), gridDelta};

    const unsigned fi = direction * 2u + (offset == 1 ? 1u : 0u);
    const std::size_t nn = nodes.size();
    const Boundary faceType = faceBCTypes_[fi * nn + nodeId];
    const T faceDist = faceBCDists_[fi * nn + nodeId];
    if (faceType == Boundary::REACTION)
      return {reactionBoundaryVelocity(node.index), faceDist};
    if (faceType == Boundary::AMBIENT)
      return {freeSurfaceVelocityBoundary(node.index, direction, offset,
                                          faceDist, toT(velocity[nodeId])),
              faceDist};
    if (faceType == Boundary::MASK)
      return {maskVelocityBoundary(node.index, toT(velocity[nodeId])),
              faceDist};

    return {toT(velocity[nodeId]), gridDelta};
  }

  StencilPoint<T> currentPressureStencilPoint(std::size_t nodeId,
                                              unsigned direction,
                                              int offset) const {
    const auto &node = nodes[nodeId];
    IndexType neighbor = node.index;
    neighbor[direction] += offset;

    if (!inBounds(neighbor))
      return {node.pressure, gridDelta};

    const std::size_t neighborId = nodeLookupFlat[linearIndex(neighbor)];
    if (neighborId != noNode)
      return {nodes[neighborId].pressure, gridDelta};

    const unsigned fi = direction * 2u + (offset == 1 ? 1u : 0u);
    const std::size_t nn = nodes.size();
    const Boundary faceType = faceBCTypes_[fi * nn + nodeId];
    const T faceDist = faceBCDists_[fi * nn + nodeId];
    if (faceType == Boundary::AMBIENT)
      return {freeSurfacePressureBoundary(node.index), faceDist};
    if (faceType == Boundary::REACTION)
      return {node.pressure, faceDist}; // Neumann ∂p/∂n=0
    if (faceType == Boundary::MASK)
      return {
          maskPressureBoundary(node.index, direction, offset, node.pressure),
          faceDist};

    return {node.pressure, gridDelta};
  }

  StencilPoint<Vec3D<T>> currentVelocityStencilPoint(std::size_t nodeId,
                                                     unsigned direction,
                                                     int offset) const {
    const auto &node = nodes[nodeId];
    IndexType neighbor = node.index;
    neighbor[direction] += offset;

    if (!inBounds(neighbor))
      return {node.velocity, gridDelta};

    const std::size_t neighborId = nodeLookupFlat[linearIndex(neighbor)];
    if (neighborId != noNode)
      return {nodes[neighborId].velocity, gridDelta};

    const unsigned fi = direction * 2u + (offset == 1 ? 1u : 0u);
    const std::size_t nn = nodes.size();
    const Boundary faceType = faceBCTypes_[fi * nn + nodeId];
    const T faceDist = faceBCDists_[fi * nn + nodeId];
    if (faceType == Boundary::REACTION)
      return {reactionBoundaryVelocity(node.index), faceDist};
    if (faceType == Boundary::AMBIENT)
      return {freeSurfaceVelocityBoundary(node.index, direction, offset,
                                          faceDist, node.velocity),
              faceDist};
    if (faceType == Boundary::MASK)
      return {maskVelocityBoundary(node.index, node.velocity), faceDist};

    return {node.velocity, gridDelta};
  }

  T maxVelocityChange(const std::vector<Vec3D<T>> &previous) const {
    T maxChange = 0.;
    T maxVelocity = 0.;
    const auto count = std::min(previous.size(), nodes.size());
    for (std::size_t i = 0; i < count; ++i) {
      for (unsigned j = 0; j < D; ++j) {
        maxChange = std::max(maxChange,
                             std::abs(nodes[i].velocity[j] - previous[i][j]));
        maxVelocity = std::max(maxVelocity, std::abs(nodes[i].velocity[j]));
      }
    }

    if (maxVelocity <= std::numeric_limits<T>::epsilon())
      return maxChange;
    return maxChange / maxVelocity;
  }

  T maxPressureChange(const std::vector<T> &previous) const {
    T maxChange = 0.;
    T maxPressure = 0.;
    const auto count = std::min(previous.size(), nodes.size());
    for (std::size_t i = 0; i < count; ++i) {
      maxChange =
          std::max(maxChange, std::abs(nodes[i].pressure - previous[i]));
      maxPressure = std::max(maxPressure, std::abs(nodes[i].pressure));
    }

    if (maxPressure <= std::numeric_limits<T>::epsilon())
      return maxChange;
    return maxChange / maxPressure;
  }

  std::array<T, 9>
  currentBoundaryDeviatoricStress(const IndexType &index) const {
    const auto strainRate = strainRateTensorAt(index);
    const auto deviatoricRate =
        deviatoricTensor(strainRate, divergenceAt(index));
    const auto previousStress = previousDeviatoricStress(index);
    const T relaxationTime = effectiveStressRelaxationTime();
    const T decay =
        (relaxationTime <= std::numeric_limits<T>::epsilon())
            ? T(0)
            : std::exp(-deformationParameters.stressTimeStep / relaxationTime);

    std::array<T, 9> deviatoricStress{};
    for (unsigned i = 0; i < 9; ++i) {
      const T viscousStress =
          T(2) * deformationParameters.viscosity * deviatoricRate[i];
      deviatoricStress[i] =
          decay * previousStress[i] + (T(1) - decay) * viscousStress;
    }

    return deviatoricStress;
  }

  T freeSurfacePressureBoundary(const IndexType &index) const {
    const auto normal = interfaceNormal(index, Boundary::AMBIENT);
    const auto deviatoricStress = currentBoundaryDeviatoricStress(index);

    return deformationParameters.ambientPressure +
           normalStress(deviatoricStress, normal);
  }

  T maskPressureBoundary(const IndexType & /*index*/, unsigned /*direction*/,
                         int /*offset*/, T fallbackPressure) const {
    return fallbackPressure;
  }

  Vec3D<T> freeSurfaceVelocityBoundary(const IndexType &index,
                                       unsigned direction, int offset,
                                       T distance,
                                       const Vec3D<T> &interiorVelocity) const {
    Vec3D<T> boundaryVelocity = interiorVelocity;
    const auto normal = interfaceNormal(index, Boundary::AMBIENT);
    const auto deviatoricStress = deviatoricStressAt(index);
    const T pressure = pressureAt(index);

    Vec3D<T> deviatoricTraction{0., 0., 0.};
    for (unsigned component = 0; component < D; ++component) {
      for (unsigned j = 0; j < D; ++j)
        deviatoricTraction[component] +=
            deviatoricStress[tensorIndex(component, j)] * normal[j];
    }

    for (unsigned component = 0; component < D; ++component) {
      const T normalTraction =
          pressure * normal[component] - deviatoricTraction[component];
      const T faceDerivative = normalTraction * normal[direction] /
                               std::max(deformationParameters.viscosity,
                                        std::numeric_limits<T>::epsilon());
      boundaryVelocity[component] +=
          static_cast<T>(offset) * distance * faceDerivative;
    }

    return boundaryVelocity;
  }

  Vec3D<T> maskVelocityBoundary(const IndexType &index,
                                const Vec3D<T> &interiorVelocity) const {
    if (maskVelocityField != nullptr) {
      Vec3D<T> coordinate{0., 0., 0.};
      for (unsigned i = 0; i < D; ++i)
        coordinate[i] = index[i] * gridDelta;
      return maskVelocityField->getVectorVelocity(
          coordinate, deformationParameters.material, {0., 0., 0.}, 0);
    }
    return {0., 0., 0.};
  }

  void computeAvgExpansionSpeed() {
    avgExpansionSpeed_ = 0.;
    if (nodes.empty())
      return;

    ConstSparseIterator reactionIt(reactionInterface->getDomain());
    ConstSparseIterator ambientIt(ambientInterface->getDomain());
    auto maskIt = makeMaskIterator();
    std::size_t count = 0;

    for (const auto &node : nodes) {
      bool touchesReactionBoundary = false;
      for (unsigned direction = 0; direction < D; ++direction) {
        for (int offset : {-1, 1}) {
          IndexType neighbor = node.index;
          neighbor[direction] += offset;
          if (!inBounds(neighbor))
            continue;

          if (lookupNode(neighbor) != noNode)
            continue;

          if (classifyBoundary(reactionIt, ambientIt, maskIt, node.index,
                               neighbor) == Boundary::REACTION) {
            touchesReactionBoundary = true;
            break;
          }
        }
        if (touchesReactionBoundary)
          break;
      }

      if (touchesReactionBoundary) {
        Vec3D<T> coordinate{0., 0., 0.};
        for (unsigned i = 0; i < D; ++i)
          coordinate[i] = node.index[i] * gridDelta;
        avgExpansionSpeed_ +=
            (oxidationParameters.expansionCoefficient - T(1)) *
            std::abs(diffusionField->getScalarVelocity(coordinate, 0,
                                                       {0., 0., 0.}, 0));
        ++count;
      }
    }

    if (count > 0)
      avgExpansionSpeed_ /= static_cast<T>(count);
  }

  Vec3D<T> reactionBoundaryVelocity(const IndexType &index) const {
    Vec3D<T> coordinate{0., 0., 0.};
    for (unsigned i = 0; i < D; ++i)
      coordinate[i] = index[i] * gridDelta;
    const T expansionVelocity = localExpansionSpeed(coordinate);
    return detail::vecScaled(reactionNormal(index),
                             reactionSign * expansionVelocity);
  }

  Vec3D<T> unresolvedAmbientVelocity(const Vec3D<T> &coordinate) const {
    if (diffusionField == nullptr || ambientInterface == nullptr)
      return {0., 0., 0.};

    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);

    ConstSparseIterator ambientIt(ambientInterface->getDomain());
    const auto normal = levelSetNormal(ambientIt, index);
    return detail::vecScaled(normal, localExpansionSpeed(coordinate));
  }

  Vec3D<T> estimateMaxUnresolvedAmbientVelocity() const {
    Vec3D<T> maxVelocity{0., 0., 0.};
    if (ambientInterface == nullptr || diffusionField == nullptr)
      return maxVelocity;

    ConstSparseIterator ambientIt(ambientInterface->getDomain());
    for (; !ambientIt.isFinished(); ++ambientIt) {
      if (!ambientIt.isDefined())
        continue;

      Vec3D<T> coordinate{0., 0., 0.};
      const auto &index = ambientIt.getStartIndices();
      for (unsigned d = 0; d < D; ++d)
        coordinate[d] = index[d] * gridDelta;

      const auto velocity = unresolvedAmbientVelocity(coordinate);
      for (unsigned d = 0; d < D; ++d)
        maxVelocity[d] = std::max(maxVelocity[d], std::abs(velocity[d]));
    }
    return maxVelocity;
  }

  T divergenceAt(const IndexType &index) const {
    T divergence = 0.;
    for (unsigned i = 0; i < D; ++i) {
      divergence += velocityDerivative(index, i, i);
    }
    return divergence;
  }

  Vec3D<T> pressureGradient(const IndexType &index) const {
    Vec3D<T> gradient{0., 0., 0.};
    for (unsigned i = 0; i < D; ++i)
      gradient[i] = pressureDerivative(index, i);
    return gradient;
  }

  Vec3D<T> momentumForcing(const IndexType &index) const {
    Vec3D<T> forcing = pressureGradient(index);
    const auto stressDivergence = deviatoricStressDivergence(index);
    for (unsigned i = 0; i < D; ++i)
      forcing[i] -= stressDivergence[i];
    return forcing;
  }

  Vec3D<T> deviatoricStressDivergence(const IndexType &index) const {
    Vec3D<T> divergence{0., 0., 0.};
    for (unsigned component = 0; component < D; ++component) {
      for (unsigned direction = 0; direction < D; ++direction) {
        IndexType pos = index;
        IndexType neg = index;
        pos[direction] += 1;
        neg[direction] -= 1;
        const auto posStress = deviatoricStressAt(pos);
        const auto negStress = deviatoricStressAt(neg);
        divergence[component] +=
            (posStress[tensorIndex(component, direction)] -
             negStress[tensorIndex(component, direction)]) /
            (T(2) * gridDelta);
      }
    }
    return divergence;
  }

  std::array<T, 9> deviatoricStressAt(const IndexType &index) const {
    if (!inBounds(index))
      return {};

    const std::size_t nodeId = lookupNode(index);
    if (nodeId == noNode)
      return {};

    std::array<T, 9> deviatoric = nodes[nodeId].stressTensor;
    for (unsigned i = 0; i < 3; ++i)
      deviatoric[tensorIndex(i, i)] += nodes[nodeId].pressure;
    return deviatoric;
  }

  T pressureAt(const IndexType &index) const {
    if (!inBounds(index))
      return deformationParameters.ambientPressure;

    const std::size_t nodeId = lookupNode(index);
    if (nodeId == noNode)
      return deformationParameters.ambientPressure;
    return nodes[nodeId].pressure;
  }

  T localExpansionSpeed(const Vec3D<T> &coordinate) const {
    return (oxidationParameters.expansionCoefficient - T(1)) *
           std::abs(diffusionField->getScalarVelocity(coordinate, 0,
                                                      {0., 0., 0.}, 0));
  }

  Vec3D<T> reactionNormal(const IndexType &index) const {
    ConstSparseIterator reactionIt(reactionInterface->getDomain());
    return levelSetNormal(reactionIt, index);
  }

  Vec3D<T> interfaceNormal(const IndexType &index, Boundary boundary) const {
    if (boundary == Boundary::AMBIENT) {
      ConstSparseIterator ambientIt(ambientInterface->getDomain());
      return levelSetNormal(ambientIt, index);
    }
    if (boundary == Boundary::MASK && maskInterface != nullptr) {
      ConstSparseIterator maskIt(maskInterface->getDomain());
      return levelSetNormal(maskIt, index);
    }

    ConstSparseIterator reactionIt(reactionInterface->getDomain());
    return levelSetNormal(reactionIt, index);
  }

  Vec3D<T> levelSetNormal(ConstSparseIterator &levelSetIt,
                          const IndexType &index) const {
    Vec3D<T> normal{0., 0., 0.};
    T norm = 0.;

    for (unsigned i = 0; i < D; ++i) {
      IndexType pos = index;
      IndexType neg = index;
      pos[i] += 1;
      neg[i] -= 1;
      if (!inBounds(pos))
        pos = index;
      if (!inBounds(neg))
        neg = index;
      normal[i] = detail::clampLevelSetPhi(valueAt(levelSetIt, pos)) -
                  detail::clampLevelSetPhi(valueAt(levelSetIt, neg));
      norm += normal[i] * normal[i];
    }

    if (norm <= std::numeric_limits<T>::epsilon()) {
      normal = Vec3D<T>{0., 0., 0.};
      normal[D - 1] = 1.;
      return normal;
    }

    norm = std::sqrt(norm);
    for (unsigned i = 0; i < D; ++i)
      normal[i] /= norm;
    return normal;
  }

  void computeDiagnostics() {
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < nodes.size(); ++i)
      nodes[i].strainTrace = divergenceAt(nodes[i].index);
  }

  void computeStressTensors() {
    const T relaxationTime = effectiveStressRelaxationTime();
    const T decay =
        (relaxationTime <= std::numeric_limits<T>::epsilon())
            ? T(0)
            : std::exp(-deformationParameters.stressTimeStep / relaxationTime);

    // Per-node computation is independent; collect history keys into a vector
    // to avoid concurrent map writes, then build the map sequentially below.
    std::vector<std::pair<IndexType, std::array<T, 9>>> historyEntries(
        nodes.size());
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      auto &node = nodes[i];
      node.strainRateTensor = strainRateTensorAt(node.index);
      const auto deviatoricRate =
          deviatoricTensor(node.strainRateTensor, node.strainTrace);
      const auto previousStress = previousDeviatoricStress(node.index);

      std::array<T, 9> deviatoricStress{};
      for (unsigned j = 0; j < 9; ++j) {
        const T viscousStress =
            T(2) * deformationParameters.viscosity * deviatoricRate[j];
        deviatoricStress[j] =
            decay * previousStress[j] + (T(1) - decay) * viscousStress;
      }

      node.stressTensor = deviatoricStress;
      for (unsigned j = 0; j < 3; ++j)
        node.stressTensor[tensorIndex(j, j)] -= node.pressure;

      node.vonMisesStress = vonMisesFromDeviatoric(deviatoricStress);
      historyEntries[i] = {node.index, deviatoricStress};
    }

    std::unordered_map<IndexType, std::array<T, 9>, detail::IndexTypeHasher<D>>
        nextHistory;
    nextHistory.reserve(nodes.size());
    for (const auto &entry : historyEntries)
      nextHistory[entry.first] = entry.second;
    deviatoricStressHistory.swap(nextHistory);
  }

  std::array<T, 9> strainRateTensorAt(const IndexType &index) const {
    std::array<T, 9> tensor{};
    for (unsigned i = 0; i < D; ++i) {
      for (unsigned j = 0; j < D; ++j) {
        tensor[tensorIndex(i, j)] = T(0.5) * (velocityDerivative(index, i, j) +
                                              velocityDerivative(index, j, i));
      }
    }
    return tensor;
  }

  T velocityDerivative(const IndexType &index, unsigned component,
                       unsigned direction) const {
    const std::size_t nodeId = lookupNode(index);
    if (nodeId == noNode)
      return 0.;

    const auto plus = currentVelocityStencilPoint(nodeId, direction, 1);
    const auto minus = currentVelocityStencilPoint(nodeId, direction, -1);
    return firstDerivative(
        minus.value[component], nodes[nodeId].velocity[component],
        plus.value[component], minus.distance, plus.distance);
  }

  T pressureDerivative(const IndexType &index, unsigned direction) const {
    const std::size_t nodeId = lookupNode(index);
    if (nodeId == noNode)
      return 0.;

    const auto plus = currentPressureStencilPoint(nodeId, direction, 1);
    const auto minus = currentPressureStencilPoint(nodeId, direction, -1);
    return firstDerivative(minus.value, nodes[nodeId].pressure, plus.value,
                           minus.distance, plus.distance);
  }

  std::array<T, 9> deviatoricTensor(const std::array<T, 9> &tensor,
                                    T trace) const {
    std::array<T, 9> result = tensor;
    const T mean = trace / T(3);
    for (unsigned i = 0; i < 3; ++i)
      result[tensorIndex(i, i)] -= mean;
    return result;
  }

  std::array<T, 9> previousDeviatoricStress(const IndexType &index) const {
    const auto found = deviatoricStressHistory.find(index);
    if (found == deviatoricStressHistory.end())
      return {};
    return found->second;
  }

  T effectiveStressRelaxationTime() const {
    if (deformationParameters.stressRelaxationTime > T(0))
      return deformationParameters.stressRelaxationTime;
    if (deformationParameters.shearModulus > std::numeric_limits<T>::epsilon())
      return deformationParameters.viscosity /
             deformationParameters.shearModulus;
    return T(0);
  }

  T vonMisesFromDeviatoric(const std::array<T, 9> &deviatoricStress) const {
    T doubleContraction = 0.;
    for (unsigned i = 0; i < 3; ++i) {
      for (unsigned j = 0; j < 3; ++j) {
        const T value = T(0.5) * (deviatoricStress[tensorIndex(i, j)] +
                                  deviatoricStress[tensorIndex(j, i)]);
        doubleContraction += value * value;
      }
    }
    return std::sqrt(T(1.5) * doubleContraction);
  }

  T normalStress(const std::array<T, 9> &tensor, const Vec3D<T> &normal) const {
    T result = 0.;
    for (unsigned i = 0; i < 3; ++i) {
      for (unsigned j = 0; j < 3; ++j)
        result += normal[i] * tensor[tensorIndex(i, j)] * normal[j];
    }
    return result;
  }

  Boundary classifyBoundary(ConstSparseIterator &reactionIt,
                            ConstSparseIterator &ambientIt,
                            ConstSparseIterator &maskIt,
                            const IndexType &inside,
                            const IndexType &outside) const {
    return boundaryIntersection(reactionIt, ambientIt, maskIt, inside, outside)
        .boundary;
  }

  BoundaryIntersection boundaryIntersection(ConstSparseIterator &reactionIt,
                                            ConstSparseIterator &ambientIt,
                                            ConstSparseIterator &maskIt,
                                            const IndexType &inside,
                                            const IndexType &outside) const {
    const T reactionInside = valueAt(reactionIt, inside);
    const T reactionOutside = valueAt(reactionIt, outside);
    const T ambientInside = valueAt(ambientIt, inside);
    const T ambientOutside = valueAt(ambientIt, outside);
    const T maskInside = valueAtMask(maskIt, inside);
    const T maskOutside = valueAtMask(maskIt, outside);

    const bool reactionCrosses = crosses(reactionInside, reactionOutside);
    const bool ambientCrosses = crosses(ambientInside, ambientOutside);
    const bool maskCrosses =
        maskInterface != nullptr && crosses(maskInside, maskOutside);

    if (!reactionCrosses && !ambientCrosses && !maskCrosses)
      return {Boundary::NONE, gridDelta};
    if (reactionCrosses && !ambientCrosses && !maskCrosses)
      return {Boundary::REACTION,
              crossingDistance(reactionInside, reactionOutside)};
    if (!reactionCrosses && ambientCrosses && !maskCrosses)
      return ambientCrossingInsideMask(
          maskInside, maskOutside,
          crossingDistance(ambientInside, ambientOutside));
    if (!reactionCrosses && !ambientCrosses && maskCrosses)
      return {Boundary::MASK, crossingDistance(maskInside, maskOutside)};

    const T reactionDistance =
        reactionCrosses ? crossingDistance(reactionInside, reactionOutside)
                        : std::numeric_limits<T>::max();
    const T ambientDistance =
        ambientCrosses ? crossingDistance(ambientInside, ambientOutside)
                       : std::numeric_limits<T>::max();
    const T maskDistance = maskCrosses
                               ? crossingDistance(maskInside, maskOutside)
                               : std::numeric_limits<T>::max();
    if (reactionDistance <= ambientDistance && reactionDistance <= maskDistance)
      return {Boundary::REACTION, reactionDistance};
    if (ambientDistance != std::numeric_limits<T>::max()) {
      const auto maskedAmbient =
          ambientCrossingInsideMask(maskInside, maskOutside, ambientDistance);
      if (maskedAmbient.boundary == Boundary::MASK)
        return maskedAmbient;
    }
    if (maskDistance <= ambientDistance)
      return {Boundary::MASK, maskDistance};
    return {Boundary::AMBIENT, ambientDistance};
  }

  bool touchesBoundary(ConstSparseIterator &reactionIt,
                       ConstSparseIterator &ambientIt,
                       ConstSparseIterator &maskIt, const IndexType &index,
                       Boundary requestedBoundary) const {
    for (unsigned direction = 0; direction < D; ++direction) {
      for (int offset : {-1, 1}) {
        IndexType neighbor = index;
        neighbor[direction] += offset;
        if (!inBounds(neighbor))
          continue;

        if (lookupNode(neighbor) != noNode)
          continue;

        if (classifyBoundary(reactionIt, ambientIt, maskIt, index, neighbor) ==
            requestedBoundary)
          return true;
      }
    }
    return false;
  }

  bool isInsideOxide(T reactionPhi, T ambientPhi) const {
    // GeometricAdvect can leave a tiny positive residual (~4*epsilon) when the
    // interface lands exactly on a grid point at non-zero coordinates, because
    // k*gridDelta is not exactly representable in floating point.  Allow a
    // tolerance of 1e-9 grid units so that grid points on the surface (phi≈0)
    // are correctly classified as inside the oxide.
    constexpr T eps = T(1e-9);
    return reactionSign * reactionPhi >= -eps &&
           ambientSign * ambientPhi >= -eps;
  }

  ConstSparseIterator makeMaskIterator() const {
    if (maskInterface == nullptr)
      return ConstSparseIterator(reactionInterface->getDomain());
    return ConstSparseIterator(maskInterface->getDomain());
  }

  bool isInsideMask(ConstSparseIterator &maskIt, const IndexType &index) const {
    if (maskInterface == nullptr)
      return false;
    return maskSign * valueAt(maskIt, index) >= 0.;
  }

  T valueAtMask(ConstSparseIterator &maskIt, const IndexType &index) const {
    if (maskInterface == nullptr)
      return std::numeric_limits<T>::max();
    return valueAt(maskIt, index);
  }

  BoundaryIntersection ambientCrossingInsideMask(T maskInside, T maskOutside,
                                                 T distance) const {
    if (isMaskAtCrossing(maskInside, maskOutside, distance))
      return {Boundary::MASK, distance};
    // Outer node is wholly inside the mask body: the oxide/gas surface has
    // drifted into the nitride. Apply mask Dirichlet BC (not traction-free)
    // so the deformation solver does not advance the surface further in.
    // Mirrors the equivalent check in lsOxidationDiffusion::classifyBoundary.
    if (maskInterface != nullptr &&
        static_cast<T>(maskSign) * maskOutside >= T(0))
      return {Boundary::MASK, distance};
    return {Boundary::AMBIENT, distance};
  }

  bool isMaskAtCrossing(T maskInside, T maskOutside, T distance) const {
    if (maskInterface == nullptr)
      return false;
    const T fraction = std::clamp(distance / gridDelta, T(0), T(1));
    const T insidePhi = detail::clampLevelSetPhi(maskInside);
    const T outsidePhi = detail::clampLevelSetPhi(maskOutside);
    const T maskPhi = insidePhi + fraction * (outsidePhi - insidePhi);
    return static_cast<T>(maskSign) * maskPhi >= T(0);
  }

  T crossingDistance(T insidePhi, T outsidePhi) const {
    return detail::levelSetCrossingDistance(
        insidePhi, outsidePhi,
        deformationParameters.minMechanicsBoundaryDistance, gridDelta);
  }

  static T firstDerivative(T minusValue, T centerValue, T plusValue,
                           T minusDistance, T plusDistance) {
    const T denominator =
        minusDistance * plusDistance * (minusDistance + plusDistance);
    if (denominator <= std::numeric_limits<T>::epsilon())
      return 0.;

    return (-plusDistance * plusDistance * minusValue +
            (plusDistance * plusDistance - minusDistance * minusDistance) *
                centerValue +
            minusDistance * minusDistance * plusValue) /
           denominator;
  }

  static constexpr unsigned tensorIndex(unsigned row, unsigned column) {
    return 3 * row + column;
  }
};

} // namespace viennals
