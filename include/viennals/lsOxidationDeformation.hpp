#pragma once

#include <lsOxidationSolverBase.hpp>
#include <lsOxidationDiffusion.hpp>

#include <algorithm>
#include <unordered_map>

#include <omp.h>

namespace viennals {

/// Parameters for the Cartesian-grid oxide deformation model.
template <class T> struct OxidationDeformationParameters {
  T viscosity = 1.;
  T bulkModulus = 1.;
  T ambientPressure = 0.;
  T pressureRelaxation = 1.;
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
  T relaxation = 1.;
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
/// is obtained from a component-wise harmonic extension and is then relaxed with
/// a Cartesian-grid quasi-static Stokes solve. The mechanical update solves a
/// pressure equation from the current velocity divergence and a velocity
/// momentum equation, including pressure-gradient and viscoelastic deviatoric
/// stress terms. The Si/SiO2 interface uses the oxidation expansion velocity,
/// the oxide/ambient interface uses a traction-free boundary, and optional mask
/// contacts use the mask velocity field.
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
    // Per-face boundary precomputation: avoids HRLE goToIndices in hot solver loops.
    // Populated once per apply() in buildNodes(); level-set positions are fixed during the solve.
    struct FaceBC { Boundary type = Boundary::NONE; T distance = 1.; };
    std::array<FaceBC, 2 * D> faceBC{};
    bool touchesAmbient = false;
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
  T avgExpansionSpeed_ = 0.;
  bool avgExpansionSpeedComputed = false;
  bool solved = false;
  bool nodesDirty_ = true;
  std::array<T, D> maxVelocity_{};
  bool useRequestedBounds = false;
  std::unordered_map<std::size_t, std::array<T, 9>> deviatoricStressHistory;

  // Warm-start storage: solutions from previous time step used as initial guess
  std::vector<Vec3D<T>> previousVelocity_;
  std::vector<T> previousPressure_;
  bool hasPreviousSolution_ = false;

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
    return SmartPointer<OxidationDeformation>::New(
        std::forward<Args>(args)...);
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

  void setMaskVelocityField(SmartPointer<VelocityField<T>> passedVelocityField) {
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

  void setDeformationParameters(
      OxidationDeformationParameters<T> passedParameters) {
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

  void markGeometryChanged() { nodesDirty_ = true; solved = false; }

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
      hasPreviousSolution_ = false; // Geometry changed, invalidate in-memory warm-start
      // Try restoring velocity, pressure, and stress history from level set
      // pointData (written by writeFieldsToLevelSet() before the previous
      // advection and remapped+filled by lsAdvect + lsInterior).
      seedFromLevelSet();
    } else if (hasPreviousSolution_ && previousVelocity_.size() == nodes.size() &&
               previousPressure_.size() == nodes.size()) {
      // Warm-start: restore previous solution as initial guess for solver.
      // Geometry stability check: only warm-start if node count matches (prevents using
      // stale solutions after grid refinement/coarsening). This typically reduces solver
      // iterations by 30-50% since the previous step's solution is close to the new one.
      for (std::size_t i = 0; i < nodes.size(); ++i) {
        nodes[i].velocity = previousVelocity_[i];
        nodes[i].pressure = previousPressure_[i];
      }
    }

    solveVelocity();
    solveMechanics();
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
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return getField(index, fallback, accessor);
  }

public:
  Vec3D<T> getVelocity(const Vec3D<T> &coordinate) const {
    return getField(coordinate, Vec3D<T>{0., 0., 0.}, [](const Node &n) { return n.velocity; });
  }

  Vec3D<T> getVelocity(const IndexType &index) const {
    return getField(index, Vec3D<T>{0., 0., 0.}, [](const Node &n) { return n.velocity; });
  }

  T getPressure(const Vec3D<T> &coordinate) const {
    return getField(coordinate, T(0), [](const Node &n) { return n.pressure; });
  }

  T getPressure(const IndexType &index) const {
    return getField(index, T(0), [](const Node &n) { return n.pressure; });
  }

  T getStrainTrace(const Vec3D<T> &coordinate) const {
    return getField(coordinate, T(0), [](const Node &n) { return n.strainTrace; });
  }

  T getStrainTrace(const IndexType &index) const {
    return getField(index, T(0), [](const Node &n) { return n.strainTrace; });
  }

  std::array<T, 9> getStrainRateTensor(const Vec3D<T> &coordinate) const {
    return getField(coordinate, std::array<T, 9>{}, [](const Node &n) { return n.strainRateTensor; });
  }

  std::array<T, 9> getStrainRateTensor(const IndexType &index) const {
    return getField(index, std::array<T, 9>{}, [](const Node &n) { return n.strainRateTensor; });
  }

  std::array<T, 9> getStressTensor(const Vec3D<T> &coordinate) const {
    return getField(coordinate, std::array<T, 9>{}, [](const Node &n) { return n.stressTensor; });
  }

  std::array<T, 9> getStressTensor(const IndexType &index) const {
    return getField(index, std::array<T, 9>{}, [](const Node &n) { return n.stressTensor; });
  }

  T getVonMisesStress(const Vec3D<T> &coordinate) const {
    return getField(coordinate, T(0), [](const Node &n) { return n.vonMisesStress; });
  }

  T getVonMisesStress(const IndexType &index) const {
    return getField(index, T(0), [](const Node &n) { return n.vonMisesStress; });
  }

  unsigned getIterations() const { return iterations; }
  T getResidual() const { return residual; }
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
      const auto key = detail::gridIndexHash<D>(idx);
      const std::size_t nId = lookupNode(idx);

      if (nId != noNode) {
        const auto &n = nodes[nId];
        velocity.push_back(n.velocity);
        const auto sIt = deviatoricStressHistory.find(key);
        if (sIt != deviatoricStressHistory.end()) {
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
    pd.insertReplaceVectorData(std::move(velocity),  "OxVelocity");
    pd.insertReplaceVectorData(std::move(stressR0),  "OxStressR0");
    pd.insertReplaceVectorData(std::move(stressR1),  "OxStressR1");
    pd.insertReplaceVectorData(std::move(stressR2),  "OxStressR2");
  }

private:
  /// Reads velocity, pressure (shared with diffusion's "OxPressure"), and
  /// stress history from ambientInterface pointData into the warm-start state.
  /// Called in apply() after buildNodes() when hasPreviousSolution_=false.
  void seedFromLevelSet() {
    if (ambientInterface == nullptr || nodes.empty())
      return;

    auto &pd = ambientInterface->getPointData();
    const int vIdx  = pd.getVectorDataIndex("OxVelocity");
    const int r0Idx = pd.getVectorDataIndex("OxStressR0");
    const int r1Idx = pd.getVectorDataIndex("OxStressR1");
    const int r2Idx = pd.getVectorDataIndex("OxStressR2");
    const int pIdx  = pd.getScalarDataIndex("OxPressure");

    const bool hasVelocity = (vIdx != -1);
    const bool hasStress   = (r0Idx != -1 && r1Idx != -1 && r2Idx != -1);
    const bool hasPressure = (pIdx  != -1);

    if (!hasVelocity && !hasStress && !hasPressure)
      return;

    const auto *vd  = hasVelocity ? pd.getVectorData(vIdx)  : nullptr;
    const auto *r0d = hasStress   ? pd.getVectorData(r0Idx) : nullptr;
    const auto *r1d = hasStress   ? pd.getVectorData(r1Idx) : nullptr;
    const auto *r2d = hasStress   ? pd.getVectorData(r2Idx) : nullptr;
    const auto *ppd = hasPressure ? pd.getScalarData(pIdx)  : nullptr;

    previousVelocity_.assign(nodes.size(), Vec3D<T>{});
    previousPressure_.assign(nodes.size(), T(0));

    ConstSparseIterator it(ambientInterface->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const auto ptId = it.getPointId();
      const IndexType idx = it.getStartIndices();
      const auto key = detail::gridIndexHash<D>(idx);
      const std::size_t ni = lookupNode(idx);

      if (ni != noNode) {
        if (vd  && ptId < static_cast<decltype(ptId)>(vd->size()))
          previousVelocity_[ni] = (*vd)[ptId];
        if (ppd && ptId < static_cast<decltype(ptId)>(ppd->size()))
          previousPressure_[ni] = (*ppd)[ptId];
      }

      if (hasStress &&
          ptId < static_cast<decltype(ptId)>(r0d->size())) {
        const auto &row0 = (*r0d)[ptId];
        const auto &row1 = (*r1d)[ptId];
        const auto &row2 = (*r2d)[ptId];
        std::array<T, 9> s{row0[0], row0[1], row0[2],
                           row1[0], row1[1], row1[2],
                           row2[0], row2[1], row2[2]};
        deviatoricStressHistory[key] = s;
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

public:
  bool initialiseGrid() {
    return initializeGridFromInterfaces(reactionInterface, ambientInterface,
                                        maskInterface, useRequestedBounds,
                                        requestedMinIndex, requestedMaxIndex,
                                        deformationParameters.maxGridPoints,
                                        "OxidationDeformation");
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

    // Precompute per-face boundary intersections. Level-set positions don't
    // change during the inner solver, so one computation per apply() suffices.
    for (auto &node : nodes) {
      node.touchesAmbient = false;
      for (unsigned dir = 0; dir < D; ++dir) {
        for (int off : {-1, 1}) {
          const unsigned fi = dir * 2u + (off == 1 ? 1u : 0u);
          IndexType nb = node.index;
          nb[dir] += off;
          if (!inBounds(nb) || lookupNode(nb) != noNode) {
            node.faceBC[fi] = {};
            continue;
          }
          const auto bi =
              boundaryIntersection(reactionIt, ambientIt, maskIt, node.index, nb);
          node.faceBC[fi] = {bi.boundary, bi.distance};
          if (bi.boundary == Boundary::AMBIENT)
            node.touchesAmbient = true;
        }
      }
    }

  }

  // Evaluates the harmonic stencil at one node.
  // sum = sum of neighbor/BC contributions (interior neighbors, reaction/mask BCs,
  // and self-coupling for OOB/NONE/AMBIENT faces).
  // count is always 2*D (every face is counted regardless of type).
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
        const auto boundary = node.faceBC[fi].type;
        if (boundary == Boundary::REACTION) {
          detail::vecAddTo(sum, reactionBoundaryVelocity(node.index));
        } else if (boundary == Boundary::MASK) {
          detail::vecAddTo(sum, maskVelocityBoundary(node.index, toT(v[nodeId])));
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

    using SolverT = float;

    const std::size_t n = nodes.size();
    const T diagVal = static_cast<T>(2 * D); // constant for all nodes

    // b[i] = BC constants (reaction + mask velocities), computed at v = zeros.
    // OOB/NONE/AMBIENT faces contribute v[i] = 0 at zeros, so only Dirichlet
    // BCs survive — correctly isolating the RHS constant vector.
    std::vector<Vec3D<T>> b(n);
    {
      const std::vector<Vec3D<SolverT>> zeros(n, Vec3D<SolverT>{SolverT(0), SolverT(0), SolverT(0)});
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < n; ++i)
        computeHarmonicStencilAt(i, zeros, b[i]);
    }

    // Warm-start from previous substep's velocity field.
    std::vector<Vec3D<SolverT>> x(n);
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned c = 0; c < D; ++c)
        x[i][c] = static_cast<SolverT>(nodes[i].velocity[c]);

    // r = b - A*x
    const Vec3D<SolverT> zero3{SolverT(0), SolverT(0), SolverT(0)};
    std::vector<Vec3D<SolverT>> Ax(n);
    harmonicMatvec(x, b, Ax);
    std::vector<Vec3D<SolverT>> r(n), r_hat(n);
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned c = 0; c < D; ++c) {
        r[i][c]     = static_cast<SolverT>(b[i][c] - Ax[i][c]);
        r_hat[i][c] = r[i][c];
      }

    // BiCGSTAB with diagonal preconditioner (diag = 2*D, constant).
    std::vector<Vec3D<SolverT>> pv(n, zero3), sv(n, zero3), y(n), z(n), s(n), t(n);
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

    for (; iterations < deformationParameters.harmonicIterations; ++iterations) {
      const T rho_new = vecDot(r_hat, r);
      if (std::abs(rho_new) < T(1e-100))
        break;

      const T beta = (rho_new / rho) * (alpha / omega);
      rho = rho_new;

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          pv[i][c] = static_cast<SolverT>(r[i][c] + beta * (pv[i][c] - omega * sv[i][c]));

      // y = M^{-1} p = p / (2*D)
      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          y[i][c] = static_cast<SolverT>(static_cast<T>(pv[i][c]) / diagVal);

      harmonicMatvec(y, b, sv);

      const T r_hat_v = vecDot(r_hat, sv);
      if (std::abs(r_hat_v) < T(1e-100))
        break;

      alpha = rho_new / r_hat_v;

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          s[i][c] = static_cast<SolverT>(r[i][c] - alpha * sv[i][c]);

      residual = vecMaxAbs(s);
      if (residual < deformationParameters.tolerance) {
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
      omega = (t_t > T(1e-100)) ? t_s / t_t : T(0);

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c) {
          x[i][c] = static_cast<SolverT>(x[i][c] + alpha * y[i][c] + omega * z[i][c]);
          r[i][c]  = static_cast<SolverT>(s[i][c] - omega * t[i][c]);
        }

      residual = vecMaxAbs(r);
      if (residual < deformationParameters.tolerance) {
        ++iterations;
        break;
      }
    }

    for (std::size_t i = 0; i < n; ++i)
      for (unsigned c = 0; c < D; ++c)
        nodes[i].velocity[c] = static_cast<T>(x[i][c]);
  }

  void solveMechanics() {
    T mechanicsResidual = 0.;
    for (unsigned iteration = 0;
         iteration < deformationParameters.mechanicsIterations; ++iteration) {
      const auto previousVelocity = collectVelocities();
      const auto previousPressure = collectPressures();

      computeDiagnostics();
      computeStressTensors();
      solvePressure();
      solveStokesVelocity();

      mechanicsResidual =
          std::max(maxVelocityChange(previousVelocity),
                   maxPressureChange(previousPressure));
      if (mechanicsResidual < deformationParameters.mechanicsTolerance)
        break;
    }

    computeDiagnostics();
    computeStressTensors();
    residual = mechanicsResidual;
  }

  // Fills diag = centerCoefficient and rhs = pressureSum for one node.
  // Dirichlet (ambient) nodes are encoded as identity rows: diag=1, rhs=ambientBP.
  template <class SolverT>
  void computePressureStencilAt(std::size_t nodeId,
                                 const std::vector<SolverT> &p,
                                 const std::vector<T> &ambientBP,
                                 T &diag, T &rhs) const {
    if (nodes[nodeId].touchesAmbient) {
      diag = T(1);
      rhs  = ambientBP[nodeId];
      return;
    }
    diag = T(0);
    rhs  = T(0);
    for (unsigned direction = 0; direction < D; ++direction) {
      const auto plus  = pressureStencilPoint(p, ambientBP, nodeId, direction,  1);
      const auto minus = pressureStencilPoint(p, ambientBP, nodeId, direction, -1);
      const T dSum = plus.distance + minus.distance;
      const T plusCoeff  = T(2) / (plus.distance  * dSum);
      const T minusCoeff = T(2) / (minus.distance * dSum);
      rhs  += plusCoeff * plus.value + minusCoeff * minus.value;
      diag += plusCoeff + minusCoeff;
    }
  }

  // (Av)[i] = precomputedDiag[i]*v[i] - rhs_at_v[i] + pBC[i]
  template <class SolverT>
  void pressureMatvec(const std::vector<SolverT> &v,
                      const std::vector<T> &ambientBP,
                      const std::vector<T> &precomputedDiag,
                      const std::vector<T> &pBC,
                      std::vector<SolverT> &Av) const {
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

    using SolverT = float;

    const std::size_t n = nodes.size();
    const T eps = std::numeric_limits<T>::epsilon();

    std::vector<T> divergence(n), ambientBP(n);
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < n; ++i) {
      divergence[i] = divergenceAt(nodes[i].index);
      ambientBP[i]  = freeSurfacePressureBoundary(nodes[i].index);
    }

    // Geometry-fixed diagonal and BC constants (kept in T for full precision).
    std::vector<T> diag(n), pBC(n);
    {
      const std::vector<SolverT> zeros(n, SolverT(0));
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < n; ++i)
        computePressureStencilAt(i, zeros, ambientBP, diag[i], pBC[i]);
    }

    std::vector<T> b(n);
    for (std::size_t i = 0; i < n; ++i)
      b[i] = pBC[i] + deformationParameters.bulkModulus * divergence[i];

    std::vector<SolverT> x(n);
    for (std::size_t i = 0; i < n; ++i)
      x[i] = static_cast<SolverT>(nodes[i].touchesAmbient ? ambientBP[i] : nodes[i].pressure);

    std::vector<SolverT> Ax(n);
    pressureMatvec(x, ambientBP, diag, pBC, Ax);
    std::vector<SolverT> r(n), r_hat(n), p(n, SolverT(0)), v(n, SolverT(0)), y(n), z(n), s(n), t(n);
    for (std::size_t i = 0; i < n; ++i) {
      r[i]     = static_cast<SolverT>(b[i] - Ax[i]);
      r_hat[i] = r[i];
    }

    T rho = T(1), alpha = T(1), omega = T(1);
    T pressureResidual = T(0);

    for (unsigned iteration = 0;
         iteration < deformationParameters.pressureIterations; ++iteration) {
      T rho_new = T(0);
      for (std::size_t i = 0; i < n; ++i)
        rho_new += static_cast<T>(r_hat[i]) * static_cast<T>(r[i]);

      if (std::abs(rho_new) < T(1e-100))
        break;

      const T beta = (rho_new / rho) * (alpha / omega);
      rho = rho_new;

      for (std::size_t i = 0; i < n; ++i)
        p[i] = static_cast<SolverT>(r[i] + beta * (p[i] - omega * v[i]));

      for (std::size_t i = 0; i < n; ++i) {
        const T pi = p[i];
        y[i] = static_cast<SolverT>((diag[i] > eps) ? pi / diag[i] : pi);
      }

      pressureMatvec(y, ambientBP, diag, pBC, v);

      T r_hat_v = T(0);
      for (std::size_t i = 0; i < n; ++i)
        r_hat_v += static_cast<T>(r_hat[i]) * static_cast<T>(v[i]);
      if (std::abs(r_hat_v) < T(1e-100))
        break;

      alpha = rho_new / r_hat_v;

      for (std::size_t i = 0; i < n; ++i)
        s[i] = static_cast<SolverT>(r[i] - alpha * v[i]);

      pressureResidual = T(0);
      for (std::size_t i = 0; i < n; ++i)
        pressureResidual = std::max(pressureResidual, std::abs(static_cast<T>(s[i])));
      if (pressureResidual < deformationParameters.pressureTolerance) {
        for (std::size_t i = 0; i < n; ++i)
          x[i] = static_cast<SolverT>(x[i] + alpha * y[i]);
        break;
      }

      for (std::size_t i = 0; i < n; ++i) {
        const T si = s[i];
        z[i] = static_cast<SolverT>((diag[i] > eps) ? si / diag[i] : si);
      }

      pressureMatvec(z, ambientBP, diag, pBC, t);

      T t_s = T(0), t_t = T(0);
      for (std::size_t i = 0; i < n; ++i) {
        t_s += static_cast<T>(t[i]) * static_cast<T>(s[i]);
        t_t += static_cast<T>(t[i]) * static_cast<T>(t[i]);
      }
      omega = (t_t > T(1e-100)) ? t_s / t_t : T(0);

      for (std::size_t i = 0; i < n; ++i) {
        x[i] = static_cast<SolverT>(x[i] + alpha * y[i] + omega * z[i]);
        r[i]  = static_cast<SolverT>(s[i] - omega * t[i]);
      }

      pressureResidual = T(0);
      for (std::size_t i = 0; i < n; ++i)
        pressureResidual = std::max(pressureResidual, std::abs(static_cast<T>(r[i])));
      if (pressureResidual < deformationParameters.pressureTolerance)
        break;
    }

    for (std::size_t i = 0; i < n; ++i)
      nodes[i].pressure = x[i];
  }

  // Fills scalar diag = centerCoefficient and Vec3D rhs = velocitySum for one node.
  template <class SolverT>
  void computeVelocityStencilAt(std::size_t nodeId,
                                 const std::vector<Vec3D<SolverT>> &v,
                                 T &diag, Vec3D<T> &rhs) const {
    diag = T(0);
    rhs  = {T(0), T(0), T(0)};
    for (unsigned direction = 0; direction < D; ++direction) {
      const auto plus  = velocityStencilPoint(v, nodeId, direction,  1);
      const auto minus = velocityStencilPoint(v, nodeId, direction, -1);
      const T dSum = plus.distance + minus.distance;
      const T plusCoeff  = T(2) / (plus.distance  * dSum);
      const T minusCoeff = T(2) / (minus.distance * dSum);
      detail::vecAddTo(rhs, detail::vecScaled(plus.value,  plusCoeff));
      detail::vecAddTo(rhs, detail::vecScaled(minus.value, minusCoeff));
      diag += plusCoeff + minusCoeff;
    }
  }

  void solveStokesVelocity() {
    if (deformationParameters.viscosity <= std::numeric_limits<T>::epsilon())
      return;
    if (nodes.empty())
      return;

    using SolverT = float;

    const std::size_t n = nodes.size();
    const T eps = std::numeric_limits<T>::epsilon();

    // Geometry-fixed diagonal, BC constants, and forcing (all in T).
    std::vector<T> diag(n);
    std::vector<Vec3D<T>> vBC(n), forcing(n);
    {
      const std::vector<Vec3D<SolverT>> zeros(n, Vec3D<SolverT>{SolverT(0), SolverT(0), SolverT(0)});
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < n; ++i) {
        computeVelocityStencilAt(i, zeros, diag[i], vBC[i]);
        forcing[i] = momentumForcing(nodes[i].index);
      }
    }

    std::vector<Vec3D<T>> b(n);
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned c = 0; c < D; ++c)
        b[i][c] = vBC[i][c] - forcing[i][c] / deformationParameters.viscosity;

    // Initial guess from current node velocities (warm-start), converted to SolverT.
    std::vector<Vec3D<SolverT>> x(n);
    {
      const auto vel = collectVelocities();
      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          x[i][c] = static_cast<SolverT>(vel[i][c]);
    }

    // Stokes SpMV: (Av)[i] = diag[i]*vin[i] - rhs_at_vin[i] + vBC[i], stored as SolverT.
    auto stokesMatvec = [&](const std::vector<Vec3D<SolverT>> &vin,
                             std::vector<Vec3D<SolverT>> &Av) {
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < n; ++i) {
        T d; Vec3D<T> rhs;
        computeVelocityStencilAt(i, vin, d, rhs);
        for (unsigned c = 0; c < D; ++c)
          Av[i][c] = static_cast<SolverT>(diag[i] * vin[i][c] - rhs[c] + vBC[i][c]);
      }
    };

    // Dot product accumulated in T for numerical stability.
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

    // r = b - A*x
    const Vec3D<SolverT> zero3{SolverT(0), SolverT(0), SolverT(0)};
    std::vector<Vec3D<SolverT>> Ax(n), r(n), r_hat(n);
    std::vector<Vec3D<SolverT>> pv(n, zero3), sv(n, zero3), y(n), z(n), s(n), t(n);
    stokesMatvec(x, Ax);
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned c = 0; c < D; ++c) {
        r[i][c]     = static_cast<SolverT>(b[i][c] - Ax[i][c]);
        r_hat[i][c] = r[i][c];
      }

    T rho = T(1), alpha = T(1), omega = T(1);
    T velocityResidual = T(0);

    for (unsigned iteration = 0;
         iteration < deformationParameters.stokesIterations; ++iteration) {
      const T rho_new = vecDot(r_hat, r);
      if (std::abs(rho_new) < T(1e-100))
        break;

      const T beta = (rho_new / rho) * (alpha / omega);
      rho = rho_new;

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          pv[i][c] = static_cast<SolverT>(r[i][c] + beta * (pv[i][c] - omega * sv[i][c]));

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c) {
          const T pvc = pv[i][c];
          y[i][c] = static_cast<SolverT>((diag[i] > eps) ? pvc / diag[i] : pvc);
        }

      stokesMatvec(y, sv);

      const T r_hat_v = vecDot(r_hat, sv);
      if (std::abs(r_hat_v) < T(1e-100))
        break;

      alpha = rho_new / r_hat_v;

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          s[i][c] = static_cast<SolverT>(r[i][c] - alpha * sv[i][c]);

      velocityResidual = vecMaxAbs(s);
      if (velocityResidual < deformationParameters.stokesTolerance) {
        for (std::size_t i = 0; i < n; ++i)
          for (unsigned c = 0; c < D; ++c)
            x[i][c] = static_cast<SolverT>(x[i][c] + alpha * y[i][c]);
        break;
      }

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c) {
          const T sc = s[i][c];
          z[i][c] = static_cast<SolverT>((diag[i] > eps) ? sc / diag[i] : sc);
        }

      stokesMatvec(z, t);

      const T t_s = vecDot(t, s);
      const T t_t = vecDot(t, t);
      omega = (t_t > T(1e-100)) ? t_s / t_t : T(0);

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c) {
          x[i][c] = static_cast<SolverT>(x[i][c] + alpha * y[i][c] + omega * z[i][c]);
          r[i][c]  = static_cast<SolverT>(s[i][c] - omega * t[i][c]);
        }

      velocityResidual = vecMaxAbs(r);
      if (velocityResidual < deformationParameters.stokesTolerance)
        break;
    }

    for (std::size_t i = 0; i < n; ++i)
      for (unsigned c = 0; c < D; ++c)
        nodes[i].velocity[c] = static_cast<T>(x[i][c]);
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
  StencilPoint<T> pressureStencilPoint(const std::vector<SolverT> &pressure,
                                       const std::vector<T> &ambientBoundaryPressure,
                                       std::size_t nodeId, unsigned direction,
                                       int offset) const {
    const auto &node = nodes[nodeId];
    IndexType neighbor = node.index;
    neighbor[direction] += offset;

    if (!inBounds(neighbor))
      return {static_cast<T>(pressure[nodeId]), gridDelta};

    const std::size_t neighborId = nodeLookupFlat[linearIndex(neighbor)];
    if (neighborId != noNode)
      return {static_cast<T>(pressure[neighborId]), gridDelta};

    const unsigned fi = direction * 2u + (offset == 1 ? 1u : 0u);
    const auto &intersection = node.faceBC[fi];
    if (intersection.type == Boundary::AMBIENT)
      return {ambientBoundaryPressure[nodeId], intersection.distance};
    if (intersection.type == Boundary::REACTION ||
        intersection.type == Boundary::MASK)
      return {solidInterfacePressureBoundary(node.index, static_cast<T>(pressure[nodeId])),
              intersection.distance};

    return {static_cast<T>(pressure[nodeId]), gridDelta};
  }

  template <class SolverT>
  StencilPoint<Vec3D<T>> velocityStencilPoint(const std::vector<Vec3D<SolverT>> &velocity,
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
    const auto &intersection = node.faceBC[fi];
    if (intersection.type == Boundary::REACTION)
      return {reactionBoundaryVelocity(node.index), intersection.distance};
    if (intersection.type == Boundary::AMBIENT)
      return {freeSurfaceVelocityBoundary(node.index, direction, offset,
                                          intersection.distance, toT(velocity[nodeId])),
              intersection.distance};
    if (intersection.type == Boundary::MASK)
      return {maskVelocityBoundary(node.index, toT(velocity[nodeId])),
              intersection.distance};

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
    const auto &intersection = node.faceBC[fi];
    if (intersection.type == Boundary::AMBIENT)
      return {freeSurfacePressureBoundary(node.index), intersection.distance};
    if (intersection.type == Boundary::REACTION ||
        intersection.type == Boundary::MASK)
      return {solidInterfacePressureBoundary(node.index, node.pressure),
              intersection.distance};

    return {node.pressure, gridDelta};
  }

  StencilPoint<Vec3D<T>>
  currentVelocityStencilPoint(std::size_t nodeId, unsigned direction,
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
    const auto &intersection = node.faceBC[fi];
    if (intersection.type == Boundary::REACTION)
      return {reactionBoundaryVelocity(node.index), intersection.distance};
    if (intersection.type == Boundary::AMBIENT)
      return {freeSurfaceVelocityBoundary(node.index, direction, offset,
                                          intersection.distance, node.velocity),
              intersection.distance};
    if (intersection.type == Boundary::MASK)
      return {maskVelocityBoundary(node.index, node.velocity),
              intersection.distance};

    return {node.velocity, gridDelta};
  }

  T maxVelocityChange(const std::vector<Vec3D<T>> &previous) const {
    T maxChange = 0.;
    T maxVelocity = 0.;
    const auto count = std::min(previous.size(), nodes.size());
    for (std::size_t i = 0; i < count; ++i) {
      for (unsigned j = 0; j < D; ++j) {
        maxChange =
            std::max(maxChange, std::abs(nodes[i].velocity[j] - previous[i][j]));
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
      maxChange = std::max(maxChange,
                           std::abs(nodes[i].pressure - previous[i]));
      maxPressure = std::max(maxPressure, std::abs(nodes[i].pressure));
    }

    if (maxPressure <= std::numeric_limits<T>::epsilon())
      return maxChange;
    return maxChange / maxPressure;
  }

  T freeSurfacePressureBoundary(const IndexType &index) const {
    const auto normal = interfaceNormal(index, Boundary::AMBIENT);
    const auto strainRate = strainRateTensorAt(index);
    const auto deviatoricRate = deviatoricTensor(strainRate, divergenceAt(index));
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

    return deformationParameters.ambientPressure +
           normalStress(deviatoricStress, normal);
  }

  T solidInterfacePressureBoundary(const IndexType &,
                                   T fallbackPressure) const {
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
      const T faceDerivative =
          normalTraction * normal[direction] /
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
      return maskVelocityField->getVectorVelocity(coordinate,
                                                  deformationParameters.material,
                                                  {0., 0., 0.}, 0);
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
    std::vector<std::pair<std::size_t, std::array<T, 9>>> historyEntries(nodes.size());
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
      historyEntries[i] = {detail::gridIndexHash<D>(node.index), deviatoricStress};
    }

    std::unordered_map<std::size_t, std::array<T, 9>> nextHistory;
    nextHistory.reserve(nodes.size());
    for (const auto &entry : historyEntries)
      nextHistory[entry.first] = entry.second;
    deviatoricStressHistory.swap(nextHistory);
  }

  std::array<T, 9> strainRateTensorAt(const IndexType &index) const {
    std::array<T, 9> tensor{};
    for (unsigned i = 0; i < D; ++i) {
      for (unsigned j = 0; j < D; ++j) {
        tensor[tensorIndex(i, j)] =
            T(0.5) * (velocityDerivative(index, i, j) +
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
    return firstDerivative(minus.value[component],
                           nodes[nodeId].velocity[component],
                           plus.value[component], minus.distance,
                           plus.distance);
  }

  T pressureDerivative(const IndexType &index, unsigned direction) const {
    const std::size_t nodeId = lookupNode(index);
    if (nodeId == noNode)
      return 0.;

    const auto plus = currentPressureStencilPoint(nodeId, direction, 1);
    const auto minus = currentPressureStencilPoint(nodeId, direction, -1);
    return firstDerivative(minus.value, nodes[nodeId].pressure,
                           plus.value, minus.distance, plus.distance);
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
    const auto found =
        deviatoricStressHistory.find(detail::gridIndexHash<D>(index));
    if (found == deviatoricStressHistory.end())
      return {};
    return found->second;
  }

  T effectiveStressRelaxationTime() const {
    if (deformationParameters.stressRelaxationTime > T(0))
      return deformationParameters.stressRelaxationTime;
    if (deformationParameters.shearModulus > std::numeric_limits<T>::epsilon())
      return deformationParameters.viscosity / deformationParameters.shearModulus;
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

  T normalStress(const std::array<T, 9> &tensor,
                 const Vec3D<T> &normal) const {
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
      return ambientCrossingInsideMask(maskInside, maskOutside,
                                       crossingDistance(ambientInside,
                                                        ambientOutside));
    if (!reactionCrosses && !ambientCrosses && maskCrosses)
      return {Boundary::MASK, crossingDistance(maskInside, maskOutside)};

    const T reactionDistance =
        reactionCrosses ? crossingDistance(reactionInside, reactionOutside)
                        : std::numeric_limits<T>::max();
    const T ambientDistance =
        ambientCrosses ? crossingDistance(ambientInside, ambientOutside)
                       : std::numeric_limits<T>::max();
    const T maskDistance =
        maskCrosses ? crossingDistance(maskInside, maskOutside)
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
    return reactionSign * reactionPhi >= 0. && ambientSign * ambientPhi >= 0.;
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
    const T denominator = minusDistance * plusDistance *
                          (minusDistance + plusDistance);
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
