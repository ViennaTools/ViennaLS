#pragma once

#include <lsOxidationSolverBase.hpp>
#include <lsVelocityField.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <vcTimer.hpp>

#include <omp.h>

#ifdef VIENNALS_GPU_BICGSTAB
#include <lsOxidationBiCGSTABInterface.hpp>
#endif

namespace viennals {

/// Selects the BiCGSTAB back-end for the diffusion solve.
/// Cpu always uses CPU. Auto uses CPU below the GPU threshold and GPU above it.
/// Once GPU is selected or requested, GPU failures are reported instead of
/// falling back to CPU.
enum class GpuMode {
  Auto, ///< GPU when node count >= kGpuThreshold (default)
  Gpu,  ///< Always use GPU; fail if unavailable or unsuccessful
  Cpu   ///< Always use CPU
};

/// Selects the preconditioner used by the GPU BiCGSTAB solver.
/// Jacobi matches the CPU solver's preconditioner and is the default. ILU0 is
/// experimental because cuSPARSE SpSV value-analysis behavior differs between
/// CUDA versions and can introduce structured solver error if stale.
enum class GpuPreconditioner {
  Jacobi,
  ILU0
};

/// Parameters for the steady oxidant diffusion model used by
/// OxidationDiffusion.
template <class T> struct OxidationParameters {
  T diffusionCoefficient = 1.;
  T reactionRate = 1.;
  T transferCoefficient = 1.;
  T equilibriumConcentration = 1.;
  T oxidantMoleculeDensity = 1.;
  T expansionCoefficient = 1.;
  T velocitySign = 1.;

  // Stress coupling for reaction rate:
  // k_eff = k * exp(-(p - p_ref) * V_k / (k_B * T)).
  // Activation volumes are in m^3, pressure is in Pa, temperature is in K.
  T temperature = 1273.15;
  T reactionActivationVolume = 0.;
  T referencePressure = 0.;

  // Stress coupling for diffusion coefficient:
  // D_eff = D * exp(-(p - p_ref) * V_D / (k_B * T)).
  T diffusionActivationVolume = 0.;

  // Crystal orientation factor on reaction rate.
  // k(n) = k * [1 + (reactionRateRatio111 - 1) * (1 - (n . crystalAxis)^2)]
  // reactionRateRatio111 = 1 disables the correction (isotropic).
  T reactionRateRatio111 = 1.;
  Vec3D<T> crystalAxis = {0., 1., 0.};

  T maskTransferCoefficient = 0.;
  T maskConcentration = 0.;
  T minBoundaryDistance = 1e-6;
  unsigned maxIterations = 10000;
  T tolerance = 1e-8;
  T relaxation = 1.;
  std::size_t maxGridPoints = 5000000;
  int material = -1;
};

/// Solves the oxidant diffusion step of the Suvorov et al. (10.1007/s10825-006-0003-z)
/// oxidation model on the Cartesian grid carrying two level sets.
///
/// The oxide is assumed to lie between the Si/SiO2 reaction interface and the
/// SiO2/O2 ambient interface. Regular Cartesian grid nodes inside this band are
/// solved with finite differences. Whenever an axis-aligned grid edge leaves
/// the oxide, a cross-point is inserted at the zero of the crossed level-set
/// function and a Robin boundary condition is applied there. The diffusion
/// update uses the resulting sub-grid distances in a nonuniform three-point
/// stencil along each Cartesian axis.
///
/// The default sign convention matches the multilayer convention described in
/// the paper: oxide is above the reaction interface and below the ambient
/// interface, i.e. reactionPhi >= 0 and ambientPhi <= 0.
template <class T, int D>
class OxidationDiffusion final : public VelocityField<T>,
                                              public OxidationSolverBase<T, D> {
  using IndexType = viennahrle::Index<D>;
  using ConstSparseIterator =
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>;
private:
  static constexpr T boltzmannConstant = T(1.380649e-23);
  static constexpr T minStressFactor = T(1.e-6);
  static constexpr T maxStressFactor = T(1.e6);

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

  struct Node {
    IndexType index;
    T concentration = 0.;
    Vec3D<T> siNormal = {0., 0., 0.}; // unit outward normal of Si surface (into oxide)
  };

  struct StencilSide {
    T distance = 1.;
    T nodeCoefficient = 1.;
    T constant = 0.;
  };

  SmartPointer<Domain<T, D>> reactionInterface = nullptr;
  SmartPointer<Domain<T, D>> ambientInterface = nullptr;
  SmartPointer<Domain<T, D>> maskInterface = nullptr;
  OxidationParameters<T> parameters;
  int reactionSign = 1;
  int ambientSign = -1;
  int maskSign = 1;
  IndexType requestedMinIndex{};
  IndexType requestedMaxIndex{};
  unsigned iterations = 0;
  T residual = std::numeric_limits<T>::max();
  T normalizedResidual_ = std::numeric_limits<T>::max();
  bool lastSolveConverged_ = false;
  T maxScalarVelocity_ = 0.;
  bool solved = false;
  bool nodesDirty_ = true;  // true → rebuild grid/nodes on next apply()
  mutable std::string lastLoggedBackend_; // suppresses repeated "using X" messages
  bool useRequestedBounds = false;
  bool warmStartable_ = false; // true when nodes[i].concentration holds a prior solution
  std::unordered_map<std::size_t, T> pressureLookup;
  std::unordered_map<std::size_t, T> concentrationCache_;
  std::vector<Node> nodes;
  // Face-major flat BC arrays: index = fi * n + nodeId, where fi in [0, 2*D).
  // Face-major layout gives coalesced GPU reads when all warp threads access
  // the same face of consecutive nodes.
  std::vector<Boundary> faceBCTypes_;
  std::vector<T> faceBCDists_;
  // Neighbor node IDs for every face (geometry-only, rebuilt with nodes).
  // Entry == noNode for boundary / out-of-bounds faces.
  std::vector<std::size_t> neighborIds_;

  GpuMode gpuMode_ = GpuMode::Auto;
  GpuPreconditioner gpuPreconditioner_ = GpuPreconditioner::Jacobi;

#ifdef VIENNALS_GPU_BICGSTAB
  // Opaque pointer to device-side buffers (GpuBiCGSTABBuffers is defined only
  // in the CUDA translation unit; we hold a raw pointer here so g++ never sees
  // the CUDA internals).
  gpu::GpuBiCGSTABBuffers* gpuBufs_ = nullptr;
  // Minimum node count required to engage the GPU solver in Auto mode.
  static constexpr std::size_t kGpuThreshold = 20000;
#endif

public:
  /// Sub-grid accurate sample of the reaction boundary crossing closest to a
  /// given grid node. Carries the interpolated concentration at the boundary
  /// point, plus the axis and offset of the edge on which the crossing lies so
  /// that callers can reconstruct the boundary-point position and normal.
  struct ReactionBoundarySample {
    bool found = false;
    IndexType nodeIndex{};
    T distance = 1.;
    T concentration = 0.;
    unsigned crossingAxis = 0; // Cartesian axis of the crossing edge
    int crossingOffset = 0;    // +1 or -1: which neighbour the crossing faces
  };

  OxidationDiffusion() = default;

  OxidationDiffusion(
      SmartPointer<Domain<T, D>> passedReactionInterface,
      SmartPointer<Domain<T, D>> passedAmbientInterface,
      OxidationParameters<T> passedParameters = {})
      : reactionInterface(passedReactionInterface),
        ambientInterface(passedAmbientInterface), parameters(passedParameters) {
  }

  ~OxidationDiffusion() {
#ifdef VIENNALS_GPU_BICGSTAB
    gpu::freeGpuBuffers(gpuBufs_);
    gpuBufs_ = nullptr;
#endif
  }

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<OxidationDiffusion>::New(
        std::forward<Args>(args)...);
  }

  /// Call after any in-place modification of the level sets (e.g. after
  /// ls::Advect) so that the next apply() rebuilds the Cartesian node grid.
  void markGeometryChanged() { nodesDirty_ = true; solved = false; }

  void setReactionInterface(SmartPointer<Domain<T, D>> passedInterface) {
    reactionInterface = passedInterface;
    nodesDirty_ = true; solved = false;
  }

  void setAmbientInterface(SmartPointer<Domain<T, D>> passedInterface) {
    ambientInterface = passedInterface;
    nodesDirty_ = true; solved = false;
  }

  void setMaskInterface(SmartPointer<Domain<T, D>> passedInterface,
                        int passedMaskSign = 1) {
    maskInterface = passedInterface;
    maskSign = (passedMaskSign < 0) ? -1 : 1;
    nodesDirty_ = true; solved = false;
  }

  void clearMaskInterface() {
    maskInterface = nullptr;
    nodesDirty_ = true; solved = false;
  }

  void setParameters(OxidationParameters<T> passedParameters) {
    parameters = passedParameters;
    solved = false;
  }

  OxidationParameters<T> getParameters() const { return parameters; }

  T getEffectiveReactionRate(const Vec3D<T> &coordinate) const {
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return getEffectiveReactionRate(index);
  }

  void clearPressureField() {
    pressureLookup.clear();
    solved = false;
  }

  void setPressure(const IndexType &index, T pressure) {
    const auto key = detail::gridIndexHash<D>(index);
    if (std::isfinite(pressure))
      pressureLookup[key] = pressure;
    else
      pressureLookup.erase(key);
    solved = false;
  }

  void setPressure(const Vec3D<T> &coordinate, T pressure) {
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    setPressure(index, pressure);
  }

  /// Set signs defining the oxide band. A node is inside oxide if
  /// reactionSign * reactionPhi >= 0 and ambientSign * ambientPhi >= 0.
  void setOxideSigns(int passedReactionSign, int passedAmbientSign) {
    reactionSign = (passedReactionSign < 0) ? -1 : 1;
    ambientSign = (passedAmbientSign < 0) ? -1 : 1;
    solved = false;
  }

  /// Restrict the dense Cartesian diffusion solve to a finite index box.
  /// This is useful for level sets with infinite boundary conditions.
  void setSolveBounds(const IndexType &passedMinIndex,
                      const IndexType &passedMaxIndex) {
    requestedMinIndex = passedMinIndex;
    requestedMaxIndex = passedMaxIndex;
    useRequestedBounds = true;
    nodesDirty_ = true; solved = false;
  }

  void clearSolveBounds() {
    useRequestedBounds = false;
    nodesDirty_ = true; solved = false;
  }

  void apply() {
    if (reactionInterface == nullptr || ambientInterface == nullptr) {
      Logger::getInstance()
          .addError("OxidationDiffusion: Missing level-set "
                    "interface.")
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
            .addWarning("OxidationDiffusion: no oxide nodes found after "
                        "buildNodes(). Verify that the reaction and ambient "
                        "level sets enclose a non-empty oxide band.")
            .print();
    }
    solveDiffusion();

    concentrationCache_.clear();
    for (const auto &node : nodes)
      concentrationCache_[detail::gridIndexHash<D>(node.index)] = node.concentration;

    maxScalarVelocity_ = 0.;
    ConstSparseIterator reactionIt(reactionInterface->getDomain());
    for (const auto &node : nodes) {
      const auto sample = reactionBoundarySampleFromNode(reactionIt, node);
      if (!sample.found)
        continue;
      const T rate = getEffectiveReactionRate(node.index);
      const T vel =
          std::abs(parameters.velocitySign) * rate * sample.concentration /
          (parameters.oxidantMoleculeDensity * parameters.expansionCoefficient);
      maxScalarVelocity_ = std::max(maxScalarVelocity_, vel);
    }

    solved = true;
  }

  T getScalarVelocity(const Vec3D<T> &coordinate, int material,
                      const Vec3D<T> &normalVector,
                      unsigned long /*pointId*/) final {
    if (!solved)
      apply();

    if (parameters.material >= 0 && material != parameters.material)
      return 0.;

    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);

    const auto boundarySample = reactionBoundarySample(index);
    const IndexType rateIndex =
        boundarySample.found ? boundarySample.nodeIndex : index;
    const T concentration =
        boundarySample.found ? boundarySample.concentration
                             : getReactionBoundaryConcentration(index);
    return parameters.velocitySign * getEffectiveReactionRate(rateIndex) *
           concentration /
           (parameters.oxidantMoleculeDensity * parameters.expansionCoefficient);
  }

  T getDissipationAlpha(int /*direction*/, int material,
                        const Vec3D<T> & /*centralDifferences*/) final {
    if (parameters.material >= 0 && material != parameters.material)
      return 0.;

    T bulkVel = std::abs(parameters.velocitySign) * parameters.reactionRate *
                parameters.equilibriumConcentration /
                (parameters.oxidantMoleculeDensity * parameters.expansionCoefficient);
    
    return std::max(maxScalarVelocity_, bulkVel);
  }

  T getConcentration(const Vec3D<T> &coordinate) const {
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return getConcentration(index);
  }

  T getConcentration(const IndexType &index) const {
    const std::size_t nodeId = lookupNode(index);
    if (nodeId == noNode) {
      const auto nearby = findNearbyNode(index);
      if (nearby == noNode)
        return 0.;
      return nodes[nearby].concentration;
    }
    return nodes[nodeId].concentration;
  }

  T getReactionBoundaryConcentration(const Vec3D<T> &coordinate) const {
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return getReactionBoundaryConcentration(index);
  }

  T getReactionBoundaryConcentration(const IndexType &index) const {
    const auto sample = reactionBoundarySample(index);
    return sample.found ? sample.concentration : getConcentration(index);
  }

  unsigned getIterations() const { return iterations; }
  T getResidual() const { return residual; }
  T getNormalizedResidual() const { return normalizedResidual_; }
  bool lastSolveConverged() const { return lastSolveConverged_; }
  std::size_t getNumberOfSolutionNodes() const { return nodes.size(); }
  bool hasFiniteConcentrationField() const {
    for (const auto &node : nodes)
      if (!std::isfinite(node.concentration))
        return false;
    return true;
  }

  const std::unordered_map<std::size_t, T> &getConcentrationCache() const {
    return concentrationCache_;
  }

  void setConcentrationCache(std::unordered_map<std::size_t, T> cache) {
    concentrationCache_ = std::move(cache);
  }

  /// Set the GPU solver selection mode.  See GpuMode for the three options.
  /// On CPU-only builds (VIENNALS_GPU_BICGSTAB not defined) this is a no-op.
  void setGpuMode(GpuMode mode) { gpuMode_ = mode; }
  /// Set the GPU BiCGSTAB preconditioner. Jacobi matches the CPU solver.
  void setGpuPreconditioner(GpuPreconditioner preconditioner) {
    gpuPreconditioner_ = preconditioner;
  }

  /// Mark the current solution as valid without re-solving. Call this before
  /// any parallel advection (lsAdvect) that uses this field as a velocity
  /// source to prevent concurrent apply() calls inside getScalarVelocity().
  void markSolved() { solved = true; }

  /// Write per-node concentration into ambientInterface->getPointData() so
  /// that lsInterior + lsAdvect can carry it across timestep boundaries.
  /// Safe to call regardless of the solved flag; uses the most recent nodes.
  void writeConcentrationToLevelSet() {
    if (nodes.empty() || ambientInterface == nullptr)
      return;

    std::vector<T> concentrations;
    ConstSparseIterator it(ambientInterface->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const IndexType idx = it.getStartIndices();
      const std::size_t nodeId = lookupNode(idx);
      // Non-solve-grid narrow-band points (mask interior, gas-phase boundary)
      // get concentration 0.  Using equilibriumConcentration here contaminates
      // the output: mask-interior points are impermeable (C≈0), and newly
      // appeared bird's-beak points inherit C=1 spuriously after advection.
      // The in-memory concentrationCache_ is the warm-start source for the
      // next substep; the level-set value is only the fallback when that cache
      // is cold, and C=0 converges just as quickly as C=1 from a cold start.
      const T value = nodeId != noNode ? nodes[nodeId].concentration : T(0);
      concentrations.push_back(std::isfinite(value) ? value : T(0));
    }
    ambientInterface->getPointData().insertReplaceScalarData(
        std::move(concentrations), "OxConcentration");
  }

  /// Write per-node pressure into ambientInterface->getPointData() so that
  /// it survives advection and can warm-start the coupling loop next step.
  void writePressureToLevelSet() {
    if (pressureLookup.empty() || ambientInterface == nullptr)
      return;

    std::vector<T> pressures;
    ConstSparseIterator it(ambientInterface->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const IndexType idx = it.getStartIndices();
      const auto pIt = pressureLookup.find(detail::gridIndexHash<D>(idx));
      const T value = pIt != pressureLookup.end() ? pIt->second : T(0);
      pressures.push_back(std::isfinite(value) ? value : T(0));
    }
    ambientInterface->getPointData().insertReplaceScalarData(
        std::move(pressures), "OxPressure");
  }

  /// Convenience wrapper: persist both concentration and pressure in one call.
  void writePersistentFields() {
    writeConcentrationToLevelSet();
    writePressureToLevelSet();
  }

  /// Return the reaction boundary sample for the grid node nearest to
  /// `coordinate`. Used by the deformation solver to get both the sub-grid
  /// concentration and the crossing edge so it can compute a sub-grid normal.
  ReactionBoundarySample
  getReactionBoundarySample(const Vec3D<T> &coordinate) const {
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return reactionBoundarySample(index);
  }

  /// Absolute scalar velocity derived from an already-computed boundary sample,
  /// avoiding a second call to reactionBoundarySample inside getScalarVelocity.
  T getScalarVelocityFromSample(const ReactionBoundarySample &sample) const {
    if (!sample.found)
      return T(0);
    return std::abs(getEffectiveReactionRate(sample.nodeIndex) *
                    sample.concentration /
                    (parameters.oxidantMoleculeDensity *
                     parameters.expansionCoefficient));
  }

private:
  bool initialiseGrid() {
    return initializeGridFromInterfaces(reactionInterface, ambientInterface,
                                        maskInterface, useRequestedBounds,
                                        requestedMinIndex, requestedMaxIndex,
                                        parameters.maxGridPoints,
                                        "OxidationDiffusion");
  }

  void buildNodes() {
    nodes.clear();
    initNodeLookup();

    // On the first substep after an outer-step boundary the in-memory cache is
    // empty. Fall back to the concentration stored in the level set's pointData
    // (written by writeConcentrationToLevelSet() before the previous advection
    // and remapped by lsAdvect + lsInterior).
    warmStartable_ = !concentrationCache_.empty();
    if (!warmStartable_) {
      // Single HRLE pass restores both concentration and pressure from the
      // pointData written by writePersistentFields() before the last advection.
      const int cIdx =
          ambientInterface->getPointData().getScalarDataIndex("OxConcentration");
      const int pIdx =
          ambientInterface->getPointData().getScalarDataIndex("OxPressure");
      const auto *cd = (cIdx != -1)
                           ? ambientInterface->getPointData().getScalarData(cIdx)
                           : nullptr;
      const auto *pd = (pIdx != -1)
                           ? ambientInterface->getPointData().getScalarData(pIdx)
                           : nullptr;
      if (cd != nullptr || pd != nullptr) {
        ConstSparseIterator it(ambientInterface->getDomain());
        for (; !it.isFinished(); ++it) {
          if (!it.isDefined())
            continue;
          const auto ptId = it.getPointId();
          const std::size_t key =
              detail::gridIndexHash<D>(it.getStartIndices());
          if (cd != nullptr &&
              ptId < static_cast<decltype(ptId)>(cd->size()) &&
              std::isfinite((*cd)[ptId]))
            concentrationCache_[key] = (*cd)[ptId];
          if (pd != nullptr &&
              ptId < static_cast<decltype(ptId)>(pd->size()) &&
              std::isfinite((*pd)[ptId]))
            pressureLookup[key] = (*pd)[ptId];
        }
        warmStartable_ = !concentrationCache_.empty();
      }
    }

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
        T seedConc = parameters.equilibriumConcentration;
        auto cacheIt = concentrationCache_.find(detail::gridIndexHash<D>(index));
        if (cacheIt != concentrationCache_.end())
          seedConc = cacheIt->second;
        if (!std::isfinite(seedConc))
          seedConc = parameters.equilibriumConcentration;
        Node newNode{index, seedConc};
        if (parameters.reactionRateRatio111 != T(1))
          newNode.siNormal = computeSiNormal(index, reactionIt);
        nodes.push_back(newNode);
      }

      if (!increment(index))
        break;
    }

    // Precompute per-face boundary intersections into flat face-major arrays.
    // Level-set positions are fixed during the inner solve; one pass per apply().
    const std::size_t n = nodes.size();
    faceBCTypes_.assign(2 * D * n, Boundary::NONE);
    faceBCDists_.assign(2 * D * n, T(1));
    ConstSparseIterator faceReactionIt(reactionInterface->getDomain());
    ConstSparseIterator faceAmbientIt(ambientInterface->getDomain());
    auto faceMaskIt = makeMaskIterator();
    for (std::size_t id = 0; id < n; ++id) {
      const auto &node = nodes[id];
      for (unsigned dir = 0; dir < D; ++dir) {
        for (int off : {-1, 1}) {
          const unsigned fi = dir * 2u + (off == 1 ? 1u : 0u);
          IndexType nb = node.index;
          nb[dir] += off;
          if (!inBounds(nb) || lookupNode(nb) != noNode)
            continue; // NONE/1.0 already set by assign()
          const auto bc = classifyBoundary(faceReactionIt, faceAmbientIt,
                                           faceMaskIt, node.index, nb);
          faceBCTypes_[fi * n + id] = bc.first;
          faceBCDists_[fi * n + id] = bc.second;
        }
      }
    }

    // Precompute per-face neighbor node IDs (geometry-only).
    // Used by computeFaceCoeffs() and uploaded once to GPU.
    neighborIds_.assign(2 * D * n, noNode);
    for (std::size_t id = 0; id < n; ++id) {
      for (unsigned dir = 0; dir < D; ++dir) {
        for (int off : {-1, 1}) {
          const unsigned fi = dir * 2u + (off == 1 ? 1u : 0u);
          IndexType nb = nodes[id].index;
          nb[dir] += off;
          if (inBounds(nb))
            neighborIds_[fi * n + id] = lookupNode(nb);
        }
      }
    }

    bool loggedBackend = false;
#ifdef VIENNALS_GPU_BICGSTAB
    // Free any stale allocation from a previous buildNodes() call.
    gpu::freeGpuBuffers(gpuBufs_);
    gpuBufs_ = nullptr;

    const bool tryGpu = (gpuMode_ == GpuMode::Gpu) ||
                        (gpuMode_ == GpuMode::Auto && n >= kGpuThreshold);
    if (tryGpu) {
      const bool useIlu0 =
          gpuPreconditioner_ == GpuPreconditioner::ILU0;
      gpuBufs_ = gpu::allocGpuBuffers(static_cast<uint32_t>(n), 2 * D,
                                      useIlu0);

      if (gpuBufs_) {
        // Convert std::size_t neighbor IDs to uint32_t for the device
        const std::size_t nf = 2u * D * n;
        std::vector<uint32_t> nb32(nf);
        for (std::size_t k = 0; k < nf; ++k)
          nb32[k] = (neighborIds_[k] == noNode)
                        ? gpu::kNoNode
                        : static_cast<uint32_t>(neighborIds_[k]);
        if (!gpu::gpuUploadNeighborIds(gpuBufs_, nb32.data(), nf)) {
          gpu::freeGpuBuffers(gpuBufs_);
          gpuBufs_ = nullptr;
          throwStrictGpuError("OxidationDiffusion: GPU mode was selected, but "
                              "uploading GPU neighbor IDs failed." +
                              gpuErrorDetail());
        }
        if (useIlu0 &&
            !gpu::gpuSetupCSR(gpuBufs_, nb32.data(),
                              static_cast<uint32_t>(n), 2 * D)) {
          gpu::freeGpuBuffers(gpuBufs_);
          gpuBufs_ = nullptr;
          throwStrictGpuError("OxidationDiffusion: GPU mode was selected, but "
                              "CUSPARSE setup for the GPU BiCGSTAB solver "
                              "failed." + gpuErrorDetail());
        }
        logDiffusionBackend(
            "GPU BiCGSTAB",
            "preconditioner=" +
                std::string(useIlu0 ? "ILU0" : "Jacobi"));
        loggedBackend = true;
      } else {
        throwStrictGpuError("OxidationDiffusion: GPU mode was selected, but "
                            "CUDA buffers could not be "
                            "allocated or the CUDA context could not be "
                            "initialized." + gpuErrorDetail());
      }
    }
#endif
    if (!loggedBackend) {
      std::string reason;
      if (gpuMode_ == GpuMode::Cpu) {
        reason = "GPU disabled by configuration";
      } else {
#ifdef VIENNALS_GPU_BICGSTAB
        reason = "node count below GPU threshold " +
                 std::to_string(kGpuThreshold);
#else
        reason = "ViennaLS was built without GPU BiCGSTAB support";
#endif
      }
      logDiffusionBackend("CPU BiCGSTAB", reason);
    }
  }

  // Precompute per-face coupling coefficients for the GPU SpMV.
  //
  //   faceCoeffs[fi * n + id] = 2*D_eff / (dist_fi * distSum_axis)
  //
  // Only interior faces (neighbor != noNode) get a nonzero entry;
  // boundary and out-of-bounds faces contribute to diag/b but not to
  // the off-diagonal coupling stored here.  Must be called after diag/b
  // are computed (or at least after D_eff is known), because D_eff
  // depends on pressure when diffusionActivationVolume != 0.
  void computeFaceCoeffs(std::vector<T> &faceCoeffs) const {
    const std::size_t n = nodes.size();
    faceCoeffs.assign(2 * D * n, T(0));
    const T eps = std::numeric_limits<T>::epsilon();
    const std::vector<T> zeros(n, T(0));

    for (std::size_t id = 0; id < n; ++id) {
      const T D_eff = getEffectiveDiffusionCoefficient(nodes[id].index);
      for (unsigned dir = 0; dir < D; ++dir) {
        const unsigned fiNeg = dir * 2u;
        const unsigned fiPos = dir * 2u + 1u;

        // Use the exact side construction used by the CPU stencil so the GPU
        // SpMV cannot drift from computeStencilAt() at cut-cell boundaries.
        const auto neg = makeStencilSide(id, zeros, dir, -1, D_eff);
        const auto pos = makeStencilSide(id, zeros, dir,  1, D_eff);
        const T distNeg = neg.distance;
        const T distPos = pos.distance;
        const T distSum = distNeg + distPos;
        if (distSum <= eps) continue;

        if (neighborIds_[fiNeg * n + id] != noNode && distNeg > eps)
          faceCoeffs[fiNeg * n + id] = T(2) * D_eff / (distNeg * distSum);
        if (neighborIds_[fiPos * n + id] != noNode && distPos > eps)
          faceCoeffs[fiPos * n + id] = T(2) * D_eff / (distPos * distSum);
      }
    }
  }

  // Evaluates the stencil at one node: fills diag = A[i,i] and
  // rhs = sum_j(A_off[i,j] * x[j]) + bc_constants[i].
  // Called with x = zeros to precompute the geometry-fixed diagonal and b.
  template <class SolverT>
  void computeStencilAt(std::size_t nodeId, const std::vector<SolverT> &x,
                        T &diag, T &rhs) const {
    diag = T(0);
    rhs  = T(0);
    const T D_eff = getEffectiveDiffusionCoefficient(nodes[nodeId].index);
    for (unsigned direction = 0; direction < D; ++direction) {
      const auto neg = makeStencilSide(nodeId, x, direction, -1, D_eff);
      const auto pos = makeStencilSide(nodeId, x, direction,  1, D_eff);
      addAxisContribution(rhs, diag, neg, pos, D_eff);
    }
  }

  // Computes Av = A * v using precomputed diagonal and b (RHS constants).
  // (Av)[i] = precomputedDiag[i]*v[i] - rhs_at_v[i] + b[i]
  // Stencil arithmetic stays in T; only storage uses SolverT.
  template <class SolverT>
  void matvec(const std::vector<SolverT> &v,
              const std::vector<T> &precomputedDiag,
              const std::vector<T> &b,
              std::vector<SolverT> &Av) const {
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      T diag, rhs;
      computeStencilAt(i, v, diag, rhs);
      Av[i] = static_cast<SolverT>(precomputedDiag[i] * v[i] - rhs + b[i]);
    }
  }

  void solveDiffusion() {
    iterations = 0;
    residual = 0.;
    normalizedResidual_ = 0.;
    lastSolveConverged_ = false;
    if (nodes.empty()) {
      lastSolveConverged_ = true;
      return;
    }

    // Work vectors use float to halve memory bandwidth in the SpMV hot path.
    // Stencil arithmetic and dot-product accumulation remain in T (double).
    using SolverT = float;

    const std::size_t n = nodes.size();
    const T eps = std::numeric_limits<T>::epsilon();

    // Geometry-fixed diagonal and BC source vector (kept in T for full precision).
    Timer<> tDiag;
    tDiag.start();
    std::vector<T> diag(n), b(n);
    {
      const std::vector<SolverT> zeros(n, SolverT(0));
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < n; ++i)
        computeStencilAt(i, zeros, diag[i], b[i]);
    }
    tDiag.finish();

    T b_norm = T(0);
    bool finiteSystem = true;
    for (std::size_t i = 0; i < n; ++i) {
      finiteSystem = finiteSystem && std::isfinite(diag[i]) &&
                     std::isfinite(b[i]);
      b_norm = std::max(b_norm, std::abs(b[i]));
    }
    if (b_norm < T(1e-100))
      b_norm = T(1);
    if (!finiteSystem) {
      residual = std::numeric_limits<T>::infinity();
      normalizedResidual_ = residual;
      Logger::getInstance()
          .addWarning("solveDiffusion: assembled non-finite matrix/RHS; "
                      "rejecting this coupled trial.")
          .print();
      return;
    }

    // Initial guess: warm-start or diagonal-preconditioned b.  Sanitize the
    // warm start so a previous failed solve cannot seed NaNs into the next one.
    auto fallbackInitialGuess = [&](std::size_t i) -> T {
      if (std::isfinite(diag[i]) && std::isfinite(b[i]) && diag[i] > eps) {
        const T guess = b[i] / diag[i];
        if (std::isfinite(guess))
          return guess;
      }
      return T(0);
    };

    std::vector<SolverT> x(n);
    for (std::size_t i = 0; i < n; ++i) {
      T guess = warmStartable_ ? nodes[i].concentration : fallbackInitialGuess(i);
      if (!std::isfinite(guess))
        guess = fallbackInitialGuess(i);
      x[i] = static_cast<SolverT>(guess);
      if (!std::isfinite(x[i]))
        x[i] = SolverT(0);
    }
    const auto initialGuess = x;

#ifdef VIENNALS_GPU_BICGSTAB
    // ── GPU BiCGSTAB path ──────────────────────────────────────────────
    // Engaged when buildNodes() allocated the GPU buffers (n >= kGpuThreshold).
    // Uploads diag, b, and face coefficients fresh each call (they are
    // pressure-dependent via D_eff), then runs the full BiCGSTAB loop on
    // the device.  Only dot-product scalars and the max-abs residual are
    // downloaded; work vectors stay on the GPU between iterations.
    if (gpu::gpuIsValid(gpuBufs_)) {
      Timer<> tPrep, tUpload, tSolve;
      tPrep.start();
      std::vector<double> diagGpu(n), bGpu(n);
      for (std::size_t i = 0; i < n; ++i) {
        diagGpu[i] = static_cast<double>(diag[i]);
        bGpu[i]    = static_cast<double>(b[i]);
      }

      // Face coupling coefficients (pressure-dependent; recomputed each call)
      std::vector<T> faceCoeffs;
      computeFaceCoeffs(faceCoeffs);
      std::vector<double> coeffGpu(faceCoeffs.size());
      for (std::size_t k = 0; k < faceCoeffs.size(); ++k)
        coeffGpu[k] = static_cast<double>(faceCoeffs[k]);

      std::vector<double> xGpu(n);
      for (std::size_t i = 0; i < n; ++i)
        xGpu[i] = static_cast<double>(x[i]);
      tPrep.finish();

      tUpload.start();
      const bool gpuUploadOk =
          gpu::gpuUploadSolverArrays(gpuBufs_, diagGpu.data(), bGpu.data(),
                                     coeffGpu.data(),
                                     static_cast<uint32_t>(n),
                                     faceCoeffs.size());
      tUpload.finish();
      if (!gpuUploadOk) {
        if (Logger::hasDebug()) {
          const std::string tag = "diffusion n=" + std::to_string(n) +
                                  " [GPU upload failed]";
          Logger::getInstance()
              .addTiming(tag + " diag/b precompute", tDiag)
              .addTiming(tag + " GPU prep+faceCoeffs", tPrep)
              .addTiming(tag + " GPU upload", tUpload)
              .print();
        }
        throwStrictGpuError("OxidationDiffusion: GPU mode was selected, but "
                            "uploading GPU solver arrays or factorizing ILU "
                            "failed." + gpuErrorDetail());
      }

      // Solve on GPU; xGpu holds the initial guess and receives the solution.
      tSolve.start();
      const bool gpuConverged =
          gpu::gpuSolveBiCGSTAB(gpuBufs_, xGpu.data(),
                                static_cast<double>(eps),
                                parameters.maxIterations,
                                static_cast<double>(parameters.tolerance),
                                iterations, residual);
      tSolve.finish();

      const bool finiteGpuSolution =
          gpuConverged &&
          std::all_of(xGpu.begin(), xGpu.end(),
                      [](double value) { return std::isfinite(value); });
      const unsigned gpuIterations = iterations;
      const T gpuResidual = residual;

      if (finiteGpuSolution) {
        // Write solution back to nodes only after convergence and finite checks.
        for (std::size_t i = 0; i < n; ++i)
          nodes[i].concentration = static_cast<T>(xGpu[i]);
        normalizedResidual_ = gpuResidual / b_norm;
        lastSolveConverged_ = std::isfinite(normalizedResidual_) &&
                              normalizedResidual_ <= parameters.tolerance;

        if (Logger::hasTiming()) {
          const std::string tag = "diffusion n=" + std::to_string(n) +
                                  " iters=" + std::to_string(iterations) +
                                  " res=" + std::to_string(residual) +
                                  " [GPU]";
          Logger::getInstance()
              .addTiming(tag + " GPU BiCGSTAB", tSolve)
              .print();
        }
        if (Logger::hasDebug()) {
          const std::string tag = "diffusion n=" + std::to_string(n) + " [GPU]";
          Logger::getInstance()
              .addTiming(tag + " diag/b precompute", tDiag)
              .addTiming(tag + " GPU prep+faceCoeffs", tPrep)
              .addTiming(tag + " GPU upload", tUpload)
              .print();
        }
        return;
      }

      if (Logger::hasDebug()) {
        const std::string tag = "diffusion n=" + std::to_string(n) +
                                " iters=" + std::to_string(gpuIterations) +
                                " res=" + std::to_string(gpuResidual) +
                                " [GPU failed]";
        Logger::getInstance()
            .addTiming(tag + " diag/b precompute", tDiag)
            .addTiming(tag + " GPU prep+faceCoeffs", tPrep)
            .addTiming(tag + " GPU upload", tUpload)
            .addTiming(tag + " GPU BiCGSTAB", tSolve)
            .print();
      }
      throwStrictGpuError(
          "OxidationDiffusion: GPU mode was selected, but GPU BiCGSTAB "
          "failed, did not converge, or produced non-finite concentrations "
          "(iters=" + std::to_string(gpuIterations) +
          ", residual=" + std::to_string(gpuResidual) + ").");
    }
#else
    if (gpuMode_ == GpuMode::Gpu) {
      throwStrictGpuError("OxidationDiffusion: explicit GPU mode was "
                          "requested, but ViennaLS was built without "
                          "VIENNALS_GPU_BICGSTAB.");
    }
#endif

    // r = b - A*x
    std::vector<SolverT> Ax(n);
    matvec(x, diag, b, Ax);
    std::vector<SolverT> r(n), r_hat(n);
    for (std::size_t i = 0; i < n; ++i) {
      r[i]     = static_cast<SolverT>(b[i] - Ax[i]);
      r_hat[i] = r[i];
    }

    // BiCGSTAB with diagonal (Jacobi) preconditioner.
    // Scalars (rho, alpha, omega, beta) and dot products stay in T for stability.
    Timer<> tCpuSolve;
    tCpuSolve.start();
    std::vector<SolverT> p(n, SolverT(0)), v(n, SolverT(0)), y(n), z(n), s(n), t(n);
    T rho = T(1), alpha = T(1), omega = T(1);
    bool bicgstabBreakdown = false;
    for (std::size_t i = 0; i < n; ++i) {
      const T ri = static_cast<T>(r[i]);
      if (!std::isfinite(ri)) {
        bicgstabBreakdown = true;
        break;
      }
      residual = std::max(residual, std::abs(ri));
    }

    for (; !bicgstabBreakdown &&
           iterations < parameters.maxIterations; ++iterations) {
      T rho_new = T(0);
      for (std::size_t i = 0; i < n; ++i)
        rho_new += static_cast<T>(r_hat[i]) * static_cast<T>(r[i]);

      if (!std::isfinite(rho_new)) {
        bicgstabBreakdown = true;
        break;
      }
      if (std::abs(rho_new) < T(1e-100))
        break;
      if (!std::isfinite(rho) || !std::isfinite(alpha) ||
          !std::isfinite(omega) || std::abs(omega) < T(1e-100)) {
        bicgstabBreakdown = true;
        break;
      }

      const T beta = (rho_new / rho) * (alpha / omega);
      if (!std::isfinite(beta)) {
        bicgstabBreakdown = true;
        break;
      }
      rho = rho_new;

      for (std::size_t i = 0; i < n; ++i)
        p[i] = static_cast<SolverT>(r[i] + beta * (p[i] - omega * v[i]));

      for (std::size_t i = 0; i < n; ++i) {
        const T pi = p[i];
        y[i] = static_cast<SolverT>((diag[i] > eps) ? pi / diag[i] : pi);
      }

      matvec(y, diag, b, v);

      T r_hat_v = T(0);
      for (std::size_t i = 0; i < n; ++i)
        r_hat_v += static_cast<T>(r_hat[i]) * static_cast<T>(v[i]);
      if (!std::isfinite(r_hat_v)) {
        bicgstabBreakdown = true;
        break;
      }
      if (std::abs(r_hat_v) < T(1e-100))
        break;

      alpha = rho_new / r_hat_v;
      if (!std::isfinite(alpha)) {
        bicgstabBreakdown = true;
        break;
      }

      for (std::size_t i = 0; i < n; ++i)
        s[i] = static_cast<SolverT>(r[i] - alpha * v[i]);

      residual = T(0);
      for (std::size_t i = 0; i < n; ++i)
        residual = std::max(residual, std::abs(static_cast<T>(s[i])));
      if (!std::isfinite(residual)) {
        bicgstabBreakdown = true;
        break;
      }
      if (residual < parameters.tolerance * b_norm) {
        for (std::size_t i = 0; i < n; ++i)
          x[i] = static_cast<SolverT>(x[i] + alpha * y[i]);
        ++iterations;
        break;
      }

      for (std::size_t i = 0; i < n; ++i) {
        const T si = s[i];
        z[i] = static_cast<SolverT>((diag[i] > eps) ? si / diag[i] : si);
      }

      matvec(z, diag, b, t);

      T t_s = T(0), t_t = T(0);
      for (std::size_t i = 0; i < n; ++i) {
        t_s += static_cast<T>(t[i]) * static_cast<T>(s[i]);
        t_t += static_cast<T>(t[i]) * static_cast<T>(t[i]);
      }
      if (!std::isfinite(t_s) || !std::isfinite(t_t)) {
        bicgstabBreakdown = true;
        break;
      }
      omega = (t_t > T(1e-100)) ? t_s / t_t : T(0);
      if (!std::isfinite(omega)) {
        bicgstabBreakdown = true;
        break;
      }

      for (std::size_t i = 0; i < n; ++i) {
        x[i] = static_cast<SolverT>(x[i] + alpha * y[i] + omega * z[i]);
        r[i]  = static_cast<SolverT>(s[i] - omega * t[i]);
      }

      residual = T(0);
      for (std::size_t i = 0; i < n; ++i)
        residual = std::max(residual, std::abs(static_cast<T>(r[i])));
      if (!std::isfinite(residual)) {
        bicgstabBreakdown = true;
        break;
      }
      if (residual < parameters.tolerance * b_norm) {
        ++iterations;
        break;
      }
    }

    if (bicgstabBreakdown)
      residual = std::numeric_limits<T>::infinity();

    const bool finiteCpuSolution =
        !bicgstabBreakdown &&
        std::all_of(x.begin(), x.end(),
                    [](SolverT value) { return std::isfinite(value); });
    if (finiteCpuSolution) {
      for (std::size_t i = 0; i < n; ++i)
        nodes[i].concentration = x[i];
    } else {
      for (std::size_t i = 0; i < n; ++i)
        nodes[i].concentration = initialGuess[i];
      residual = std::numeric_limits<T>::infinity();
      Logger::getInstance()
          .addWarning("solveDiffusion: CPU BiCGSTAB produced non-finite "
                      "concentrations; keeping the sanitized initial guess.")
          .print();
    }
    normalizedResidual_ = residual / b_norm;
    lastSolveConverged_ = finiteCpuSolution &&
                          std::isfinite(normalizedResidual_) &&
                          normalizedResidual_ <= parameters.tolerance;

    tCpuSolve.finish();
    if (Logger::hasTiming()) {
#ifdef VIENNALS_GPU_BICGSTAB
      const std::string path =
          gpuMode_ == GpuMode::Cpu
              ? " [CPU]"
              : " [CPU n<" + std::to_string(kGpuThreshold) + "]";
#else
      const std::string path = " [CPU]";
#endif
      const std::string tag = "diffusion n=" + std::to_string(n) +
                              " iters=" + std::to_string(iterations) +
                              " res=" + std::to_string(residual) + path;
      Logger::getInstance()
          .addTiming(tag + " CPU BiCGSTAB", tCpuSolve)
          .print();
    }
    if (Logger::hasDebug()) {
#ifdef VIENNALS_GPU_BICGSTAB
      const std::string path =
          gpuMode_ == GpuMode::Cpu
              ? " [CPU]"
              : " [CPU n<" + std::to_string(kGpuThreshold) + "]";
#else
      const std::string path = " [CPU]";
#endif
      Logger::getInstance()
          .addTiming("diffusion n=" + std::to_string(n) + path +
                     " diag/b precompute", tDiag)
          .print();
    }
    if (residual > parameters.tolerance * b_norm)
      Logger::getInstance()
          .addWarning("solveDiffusion: BiCGSTAB did not converge after " +
                      std::to_string(iterations) + "/" +
                      std::to_string(parameters.maxIterations) +
                      " iterations (residual=" + std::to_string(residual / b_norm) +
                      ", tolerance=" + std::to_string(parameters.tolerance) + ")")
          .print();
  }

  [[noreturn]] static void throwStrictGpuError(const std::string &message) {
    Logger::getInstance().addError(message).print();
    throw std::runtime_error(message);
  }

#ifdef VIENNALS_GPU_BICGSTAB
  static std::string gpuErrorDetail() {
    const char *detail = gpu::gpuGetLastErrorMessage();
    if (detail && detail[0] != '\0')
      return std::string(" Detail: ") + detail;
    return {};
  }
#endif

  void logDiffusionBackend(const std::string &backend,
                           const std::string &detail) const {
    if (!Logger::hasInfo())
      return;
    const std::string msg = "OxidationDiffusion: using " + backend +
                            " for diffusion solve (nodes=" +
                            std::to_string(nodes.size()) +
                            (detail.empty() ? std::string() : ", " + detail) +
                            ").";
    if (msg == lastLoggedBackend_)
      return;
    lastLoggedBackend_ = msg;
    Logger::getInstance().addInfo(msg).print();
  }

  // Uses precomputed flat faceBC arrays — no HRLE access, safe for parallel execution.
  template <class SolverT>
  StencilSide makeStencilSide(std::size_t nodeId, const std::vector<SolverT> &previous,
                              unsigned direction, int offset, T diffusion) const {
    const IndexType &nodeIndex = nodes[nodeId].index;
    IndexType neighbor = nodeIndex;
    neighbor[direction] += offset;

    if (!inBounds(neighbor))
      return zeroFluxSide();

    const std::size_t neighborId = lookupNode(neighbor);
    if (neighborId != noNode)
      return {gridDelta, 0., previous[neighborId]};

    const unsigned fi = direction * 2u + (offset == 1 ? 1u : 0u);
    const std::size_t n = nodes.size();
    const Boundary faceType = faceBCTypes_[fi * n + nodeId];
    const T faceDist       = faceBCDists_[fi * n + nodeId];
    if (faceType == Boundary::REACTION)
      return reactionBoundarySide(nodeIndex, faceDist, diffusion);
    if (faceType == Boundary::AMBIENT)
      return ambientBoundarySide(faceDist, diffusion);
    if (faceType == Boundary::MASK)
      return maskBoundarySide(faceDist, diffusion);
    return zeroFluxSide();
  }

  void addAxisContribution(T &rightHandSide, T &diagonal,
                           const StencilSide &negativeSide,
                           const StencilSide &positiveSide, T diffusion) const {
    const T distanceSum = negativeSide.distance + positiveSide.distance;
    if (distanceSum <= std::numeric_limits<T>::epsilon())
      return;

    addSideContribution(rightHandSide, diagonal, negativeSide, distanceSum, diffusion);
    addSideContribution(rightHandSide, diagonal, positiveSide, distanceSum, diffusion);
  }

  void addSideContribution(T &rightHandSide, T &diagonal,
                           const StencilSide &side, T distanceSum,
                           T diffusion) const {
    if (side.distance <= std::numeric_limits<T>::epsilon())
      return;

    const T coefficient = T(2) * diffusion / (side.distance * distanceSum);
    rightHandSide += coefficient * side.constant;
    diagonal += coefficient * (T(1) - side.nodeCoefficient);
  }

  StencilSide zeroFluxSide() const { return {gridDelta, 1., 0.}; }

  StencilSide reactionBoundarySide(const IndexType &nodeIndex, T distance,
                                   T diffusion) const {
    const T conductance = diffusion / distance;
    const T reactionRate = getEffectiveReactionRate(nodeIndex);
    const T denominator = conductance + reactionRate;
    if (denominator <= std::numeric_limits<T>::epsilon())
      return zeroFluxSide();

    return {distance, conductance / denominator, 0.};
  }

  ReactionBoundarySample reactionBoundarySample(const IndexType &index) const {
    ConstSparseIterator reactionIt(reactionInterface->getDomain());

    const std::size_t directId = lookupNode(index);
    if (directId != noNode) {
      const auto sample =
          reactionBoundarySampleFromNode(reactionIt, nodes[directId]);
      if (sample.found)
        return sample;
    }

    // The reaction velocity is a local interface quantity.  A global nearest
    // search can smear open-window concentrations deep under a mask if the
    // local oxide band is temporarily missing or clipped near contact.  Keep
    // the fallback strictly local, consistent with the velocity-extension
    // radius used elsewhere in the oxidation fields.
    std::size_t bestNode = std::numeric_limits<std::size_t>::max();
    T bestDistance2 = std::numeric_limits<T>::max();

    for (int radius = 1; radius <= 4; ++radius) {
      IndexType offset{};
      offset.fill(-radius);
      while (true) {
        IndexType candidate = index;
        T distance2 = 0.;
        for (unsigned d = 0; d < D; ++d) {
          candidate[d] += offset[d];
          distance2 += static_cast<T>(offset[d] * offset[d]);
        }

        if (distance2 > T(0) && inBounds(candidate)) {
          const std::size_t foundId = nodeLookupFlat[linearIndex(candidate)];
          if (foundId != noNode && distance2 < bestDistance2) {
            const auto sample =
                reactionBoundarySampleFromNode(reactionIt, nodes[foundId]);
            if (sample.found) {
              bestDistance2 = distance2;
              bestNode = foundId;
            }
          }
        }

        unsigned dim = 0;
        for (; dim < D; ++dim) {
          if (offset[dim] < radius) {
            ++offset[dim];
            break;
          }
          offset[dim] = -radius;
        }
        if (dim == D)
          break;
      }

      if (bestNode != std::numeric_limits<std::size_t>::max())
        break;
    }

    if (bestNode == std::numeric_limits<std::size_t>::max())
      return {};
    return reactionBoundarySampleFromNode(reactionIt, nodes[bestNode]);
  }

  ReactionBoundarySample
  reactionBoundarySampleFromNode(ConstSparseIterator &reactionIt,
                                 const Node &node) const {
    ReactionBoundarySample best;
    best.nodeIndex = node.index;
    T bestDistance = std::numeric_limits<T>::max();
    const T insidePhi = valueAt(reactionIt, node.index);

    for (unsigned direction = 0; direction < D; ++direction) {
      for (const int offset : {-1, 1}) {
        IndexType neighbor = node.index;
        neighbor[direction] += offset;
        if (!inBounds(neighbor))
          continue;

        const T outsidePhi = valueAt(reactionIt, neighbor);
        if (!crosses(insidePhi, outsidePhi))
          continue;

        const T distance = crossingDistance(insidePhi, outsidePhi);
        if (distance >= bestDistance)
          continue;

        const T diffusion = getEffectiveDiffusionCoefficient(node.index);
        const auto side = reactionBoundarySide(node.index, distance, diffusion);
        best.found = true;
        best.distance = distance;
        best.concentration =
            side.nodeCoefficient * node.concentration + side.constant;
        best.crossingAxis = direction;
        best.crossingOffset = offset;
        bestDistance = distance;
      }
    }

    return best;
  }

  StencilSide ambientBoundarySide(T distance, T diffusion) const {
    const T conductance = diffusion / distance;
    const T denominator = conductance + parameters.transferCoefficient;
    if (denominator <= std::numeric_limits<T>::epsilon())
      return zeroFluxSide();

    return {distance, conductance / denominator,
            parameters.transferCoefficient *
                parameters.equilibriumConcentration / denominator};
  }

  StencilSide maskBoundarySide(T distance, T diffusion) const {
    if (parameters.maskTransferCoefficient <= std::numeric_limits<T>::epsilon())
      return zeroFluxSide();

    const T conductance = diffusion / distance;
    const T denominator = conductance + parameters.maskTransferCoefficient;
    if (denominator <= std::numeric_limits<T>::epsilon())
      return zeroFluxSide();

    return {distance, conductance / denominator,
            parameters.maskTransferCoefficient * parameters.maskConcentration /
                denominator};
  }

  T getEffectiveReactionRate(const IndexType &index) const {
    T rate = parameters.reactionRate;

    if (parameters.reactionActivationVolume != T(0)) {
      T pressure = parameters.referencePressure;
      const auto foundPressure = pressureLookup.find(detail::gridIndexHash<D>(index));
      if (foundPressure != pressureLookup.end())
        pressure = foundPressure->second;
      if (!std::isfinite(pressure))
        pressure = parameters.referencePressure;
      const T exponent = stressExponent(pressure,
                                        parameters.reactionActivationVolume);
      rate *= stressFactor(exponent);
    }

    if (parameters.reactionRateRatio111 != T(1)) {
      const std::size_t nodeId = lookupNode(index);
      if (nodeId != noNode) {
        const auto &normal = nodes[nodeId].siNormal;
        T dot = T(0);
        for (unsigned d = 0; d < D; ++d)
          dot += normal[d] * parameters.crystalAxis[d];
        rate *= T(1) + (parameters.reactionRateRatio111 - T(1)) * (T(1) - dot * dot);
      }
    }

    return rate;
  }

  T getEffectiveDiffusionCoefficient(const IndexType &index) const {
    if (parameters.diffusionActivationVolume == T(0))
      return parameters.diffusionCoefficient;

    T pressure = parameters.referencePressure;
    const auto found = pressureLookup.find(detail::gridIndexHash<D>(index));
    if (found != pressureLookup.end())
      pressure = found->second;
    if (!std::isfinite(pressure))
      pressure = parameters.referencePressure;

    const T exponent = stressExponent(pressure,
                                      parameters.diffusionActivationVolume);
    return parameters.diffusionCoefficient * stressFactor(exponent);
  }

  T stressExponent(T pressure, T activationVolume) const {
    const T thermalEnergy =
        boltzmannConstant * std::max(parameters.temperature, T(1.));
    return -(pressure - parameters.referencePressure) * activationVolume /
           thermalEnergy;
  }

  static T stressFactor(T exponent) {
    if (!std::isfinite(exponent))
      return T(1);
    if (exponent <= std::log(minStressFactor))
      return minStressFactor;
    if (exponent >= std::log(maxStressFactor))
      return maxStressFactor;
    return std::exp(exponent);
  }

  Vec3D<T> computeSiNormal(const IndexType &index,
                           ConstSparseIterator &reactionIt) const {
    // Reflect neighbor indices that fall outside the HRLE grid.  This handles
    // REFLECTIVE boundary conditions correctly: phi(b-k) = phi(b+k), so the
    // centered difference at b gives zero lateral gradient as expected.
    auto reflectToGrid = [&](IndexType idx) {
      auto &g = reactionInterface->getGrid();
      for (unsigned d2 = 0; d2 < D; ++d2) {
        const auto lo = g.getMinGridPoint(d2);
        const auto hi = g.getMaxGridPoint(d2);
        if (idx[d2] < lo) idx[d2] = 2 * lo - idx[d2];
        if (idx[d2] > hi) idx[d2] = 2 * hi - idx[d2];
      }
      return idx;
    };
    Vec3D<T> gradient{0., 0., 0.};
    for (unsigned d = 0; d < D; ++d) {
      IndexType plus = index, minus = index;
      plus[d] += 1;
      minus[d] -= 1;
      gradient[d] =
          (detail::clampLevelSetPhi(valueAt(reactionIt, reflectToGrid(plus))) -
           detail::clampLevelSetPhi(valueAt(reactionIt, reflectToGrid(minus)))) /
          (T(2) * gridDelta);
    }
    T len = T(0);
    for (unsigned d = 0; d < D; ++d)
      len += gradient[d] * gradient[d];
    len = std::sqrt(len);
    if (len > std::numeric_limits<T>::epsilon()) {
      for (unsigned d = 0; d < D; ++d)
        gradient[d] /= len;
    } else {
      gradient = Vec3D<T>{0., 0., 0.};
      gradient[D - 1] = T(1);
    }
    return gradient;
  }

  std::pair<Boundary, T>
  classifyBoundary(ConstSparseIterator &reactionIt, ConstSparseIterator &ambientIt,
                   ConstSparseIterator &maskIt,
                   const IndexType &inside, const IndexType &outside) const {
    const T reactionInside = valueAt(reactionIt, inside);
    const T reactionOutside = valueAt(reactionIt, outside);
    const T ambientInside = valueAt(ambientIt, inside);
    const T ambientOutside = valueAt(ambientIt, outside);
    const T maskInside = valueAtMask(maskIt, inside);
    const T maskOutside = valueAtMask(maskIt, outside);

    const T reactionDistance =
        crosses(reactionInside, reactionOutside)
            ? crossingDistance(reactionInside, reactionOutside)
            : std::numeric_limits<T>::max();
    const T ambientDistance =
        crosses(ambientInside, ambientOutside)
            ? crossingDistance(ambientInside, ambientOutside)
            : std::numeric_limits<T>::max();
    const T maskDistance =
        maskInterface != nullptr && crosses(maskInside, maskOutside)
            ? crossingDistance(maskInside, maskOutside)
            : std::numeric_limits<T>::max();

    if (reactionDistance == std::numeric_limits<T>::max() &&
        ambientDistance == std::numeric_limits<T>::max() &&
        maskDistance == std::numeric_limits<T>::max())
      return {Boundary::NONE, gridDelta};

    if (reactionDistance <= ambientDistance && reactionDistance <= maskDistance)
      return {Boundary::REACTION, reactionDistance};

    // Ambient crossing heads into mask-occupied space: the oxide/gas surface
    // has drifted under the nitride. Classify as MASK regardless of whether
    // the crossings are coincident (isMaskAtCrossing) or the outer node is
    // wholly inside the mask body (valueAtMask check). This is equivalent to
    // sealing the diffusion domain with UNION(ambientInterface, maskInterface)
    // without mutating any level set.
    if (ambientDistance != std::numeric_limits<T>::max() &&
        (isMaskAtCrossing(maskInside, maskOutside, ambientDistance) ||
         (maskInterface != nullptr &&
          static_cast<T>(maskSign) * valueAtMask(maskIt, outside) >= T(0))))
      return {Boundary::MASK, ambientDistance};

    if (maskDistance <= ambientDistance)
      return {Boundary::MASK, maskDistance};
    return {Boundary::AMBIENT, ambientDistance};
  }

  bool isInsideOxide(T reactionPhi, T ambientPhi) const {
    // GeometricAdvect can leave a tiny positive residual (~4*epsilon) when the
    // interface lands exactly on a grid point at non-zero coordinates, because
    // k*gridDelta is not exactly representable in floating point.  Allow a
    // tolerance of 1e-9 grid units so that grid points on the surface (phi≈0)
    // are correctly classified as inside the oxide.
    constexpr T eps = T(1e-9);
    return reactionSign * reactionPhi >= -eps && ambientSign * ambientPhi >= -eps;
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
    return detail::levelSetCrossingDistance(insidePhi, outsidePhi,
                                            parameters.minBoundaryDistance,
                                            gridDelta);
  }
};

} // namespace viennals
