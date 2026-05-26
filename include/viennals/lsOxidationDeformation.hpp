#pragma once

#include <lsOxidationSolverBase.hpp>
#include <lsOxidationDiffusion.hpp>

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
  using OxidationSolverBase<T, D>::nodeLookup;
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
  std::array<T, D> maxVelocity_{};
  bool useRequestedBounds = false;
  std::unordered_map<std::size_t, std::array<T, 9>> deviatoricStressHistory;

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
    solved = false;
  }

  void setAmbientInterface(SmartPointer<Domain<T, D>> passedInterface) {
    ambientInterface = passedInterface;
    solved = false;
  }

  void setMaskInterface(SmartPointer<Domain<T, D>> passedInterface,
                        int passedMaskSign = 1) {
    maskInterface = passedInterface;
    maskSign = (passedMaskSign < 0) ? -1 : 1;
    solved = false;
  }

  void clearMaskInterface() {
    maskInterface = nullptr;
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
    solved = false;
  }

  void clearSolveBounds() {
    useRequestedBounds = false;
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

    initialiseGrid();
    buildNodes();
    solveVelocity();
    solveMechanics();
    avgExpansionSpeedComputed = false;

    maxVelocity_.fill(T(0));
    for (const auto &node : nodes) {
      for (unsigned d = 0; d < D; ++d)
        maxVelocity_[d] = std::max(maxVelocity_[d], std::abs(node.velocity[d]));
    }

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

    return getVelocity(coordinate);
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
    const auto found = nodeLookup.find(linearIndex(index));
    if (found != nodeLookup.end())
      return accessor(nodes[found->second]);

    const auto nearby = findNearbyNode(index);
    if (nearby == std::numeric_limits<std::size_t>::max())
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

private:
  void initialiseGrid() {
    initializeGridFromInterfaces(reactionInterface, ambientInterface,
                                 maskInterface, useRequestedBounds,
                                 requestedMinIndex, requestedMaxIndex,
                                 deformationParameters.maxGridPoints,
                                 "OxidationDeformation");
  }

  void buildNodes() {
    nodes.clear();
    nodeLookup.clear();

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
        nodeLookup.emplace(linearIndex(index), id);
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
          if (!inBounds(nb) || nodeLookup.count(linearIndex(nb))) {
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

  void solveVelocity() {
    iterations = 0;
    residual = 0.;
    if (nodes.empty())
      return;

    std::vector<Vec3D<T>> previous(nodes.size(), {0., 0., 0.});
    std::vector<Vec3D<T>> next = previous;

    for (; iterations < deformationParameters.harmonicIterations; ++iterations) {
      residual = 0.;
      for (std::size_t nodeId = 0; nodeId < nodes.size(); ++nodeId) {
        const auto &node = nodes[nodeId];
        Vec3D<T> sum{0., 0., 0.};
        unsigned count = 0;

        for (unsigned direction = 0; direction < D; ++direction) {
          for (int offset : {-1, 1}) {
            IndexType neighbor = node.index;
            neighbor[direction] += offset;
            if (!inBounds(neighbor)) {
              detail::vecAddTo(sum, previous[nodeId]);
              ++count;
              continue;
            }

            const auto foundNeighbor = nodeLookup.find(linearIndex(neighbor));
            if (foundNeighbor != nodeLookup.end()) {
              detail::vecAddTo(sum, previous[foundNeighbor->second]);
              ++count;
              continue;
            }

            const unsigned fi = direction * 2u + (offset == 1 ? 1u : 0u);
            const auto boundary = node.faceBC[fi].type;
            if (boundary == Boundary::REACTION) {
              detail::vecAddTo(sum, reactionBoundaryVelocity(node.index));
            } else if (boundary == Boundary::MASK) {
              detail::vecAddTo(sum,
                               maskVelocityBoundary(node.index, previous[nodeId]));
            } else {
              detail::vecAddTo(sum, previous[nodeId]);
            }
            ++count;
          }
        }

        const Vec3D<T> updated =
            (count == 0) ? previous[nodeId] : detail::vecScaled(sum, T(1) / count);
        next[nodeId] =
            detail::vecAdd(detail::vecScaled(updated, deformationParameters.relaxation),
                           detail::vecScaled(previous[nodeId],
                                             T(1) - deformationParameters.relaxation));
        Vec3D<T> delta = detail::vecSubtract(next[nodeId], previous[nodeId]);
        for (unsigned i = 0; i < D; ++i)
          residual = std::max(residual, std::abs(delta[i]));
      }

      previous.swap(next);
      if (residual < deformationParameters.tolerance)
        break;
    }

    for (std::size_t i = 0; i < nodes.size(); ++i)
      nodes[i].velocity = previous[i];
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

  void solvePressure() {
    if (nodes.empty())
      return;

    std::vector<T> divergence(nodes.size(), 0.);
    for (std::size_t i = 0; i < nodes.size(); ++i)
      divergence[i] = divergenceAt(nodes[i].index);

    std::vector<T> previous = collectPressures();
    std::vector<T> next = previous;
    std::vector<T> ambientBoundaryPressure(nodes.size(), 0.);
    for (std::size_t i = 0; i < nodes.size(); ++i)
      ambientBoundaryPressure[i] = freeSurfacePressureBoundary(nodes[i].index);

    T pressureResidual = 0.;
    for (unsigned iteration = 0;
         iteration < deformationParameters.pressureIterations; ++iteration) {
      pressureResidual = 0.;

      for (std::size_t nodeId = 0; nodeId < nodes.size(); ++nodeId) {
        const auto &node = nodes[nodeId];
        if (node.touchesAmbient) {
          next[nodeId] = ambientBoundaryPressure[nodeId];
          continue;
        }

        T pressureSum = 0.;
        T centerCoefficient = 0.;
        for (unsigned direction = 0; direction < D; ++direction) {
          const auto plus = pressureStencilPoint(
              previous, ambientBoundaryPressure, nodeId, direction, 1);
          const auto minus = pressureStencilPoint(
              previous, ambientBoundaryPressure, nodeId, direction, -1);
          const T plusCoefficient =
              T(2) / (plus.distance * (plus.distance + minus.distance));
          const T minusCoefficient =
              T(2) / (minus.distance * (plus.distance + minus.distance));
          pressureSum += plusCoefficient * plus.value +
                         minusCoefficient * minus.value;
          centerCoefficient += plusCoefficient + minusCoefficient;
        }

        const T source = -deformationParameters.bulkModulus * divergence[nodeId];
        const T updated =
            (centerCoefficient <= std::numeric_limits<T>::epsilon())
                ? previous[nodeId]
                : (pressureSum - source) / centerCoefficient;
        next[nodeId] =
            deformationParameters.pressureRelaxation * updated +
            (T(1) - deformationParameters.pressureRelaxation) *
                previous[nodeId];
        pressureResidual =
            std::max(pressureResidual, std::abs(next[nodeId] - previous[nodeId]));
      }

      previous.swap(next);
      if (pressureResidual < deformationParameters.pressureTolerance)
        break;
    }

    for (std::size_t i = 0; i < nodes.size(); ++i)
      nodes[i].pressure = previous[i];
  }

  void solveStokesVelocity() {
    if (deformationParameters.viscosity <= std::numeric_limits<T>::epsilon())
      return;

    std::vector<Vec3D<T>> previous = collectVelocities();
    std::vector<Vec3D<T>> next = previous;

    T velocityResidual = 0.;
    for (unsigned iteration = 0;
         iteration < deformationParameters.stokesIterations; ++iteration) {
      velocityResidual = 0.;

      for (std::size_t nodeId = 0; nodeId < nodes.size(); ++nodeId) {
        const auto &node = nodes[nodeId];
        Vec3D<T> sum{0., 0., 0.};
        T centerCoefficient = 0.;

        for (unsigned direction = 0; direction < D; ++direction) {
          const auto plus = velocityStencilPoint(previous, nodeId, direction, 1);
          const auto minus = velocityStencilPoint(previous, nodeId, direction, -1);
          const T plusCoefficient =
              T(2) / (plus.distance * (plus.distance + minus.distance));
          const T minusCoefficient =
              T(2) / (minus.distance * (plus.distance + minus.distance));
          detail::vecAddTo(sum, detail::vecScaled(plus.value, plusCoefficient));
          detail::vecAddTo(sum, detail::vecScaled(minus.value, minusCoefficient));
          centerCoefficient += plusCoefficient + minusCoefficient;
        }

        if (centerCoefficient <= std::numeric_limits<T>::epsilon())
          continue;

        const auto forcing = momentumForcing(node.index);
        Vec3D<T> updated = sum;
        for (unsigned component = 0; component < D; ++component) {
          updated[component] =
              (sum[component] -
               forcing[component] / deformationParameters.viscosity) /
              centerCoefficient;
        }

        next[nodeId] =
            detail::vecAdd(detail::vecScaled(updated, deformationParameters.relaxation),
                           detail::vecScaled(previous[nodeId],
                                             T(1) - deformationParameters.relaxation));

        const Vec3D<T> delta = detail::vecSubtract(next[nodeId], previous[nodeId]);
        for (unsigned component = 0; component < D; ++component)
          velocityResidual =
              std::max(velocityResidual, std::abs(delta[component]));
      }

      previous.swap(next);
      if (velocityResidual < deformationParameters.stokesTolerance)
        break;
    }

    for (std::size_t i = 0; i < nodes.size(); ++i)
      nodes[i].velocity = previous[i];
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

  StencilPoint<T> pressureStencilPoint(const std::vector<T> &pressure,
                                       const std::vector<T> &ambientBoundaryPressure,
                                       std::size_t nodeId, unsigned direction,
                                       int offset) const {
    const auto &node = nodes[nodeId];
    IndexType neighbor = node.index;
    neighbor[direction] += offset;

    if (!inBounds(neighbor))
      return {pressure[nodeId], gridDelta};

    const auto foundNeighbor = nodeLookup.find(linearIndex(neighbor));
    if (foundNeighbor != nodeLookup.end())
      return {pressure[foundNeighbor->second], gridDelta};

    const unsigned fi = direction * 2u + (offset == 1 ? 1u : 0u);
    const auto &intersection = node.faceBC[fi];
    if (intersection.type == Boundary::AMBIENT)
      return {ambientBoundaryPressure[nodeId], intersection.distance};
    if (intersection.type == Boundary::REACTION ||
        intersection.type == Boundary::MASK)
      return {solidInterfacePressureBoundary(node.index, pressure[nodeId]),
              intersection.distance};

    return {pressure[nodeId], gridDelta};
  }

  StencilPoint<Vec3D<T>> velocityStencilPoint(const std::vector<Vec3D<T>> &velocity,
                                              std::size_t nodeId, unsigned direction,
                                              int offset) const {
    const auto &node = nodes[nodeId];
    IndexType neighbor = node.index;
    neighbor[direction] += offset;

    if (!inBounds(neighbor))
      return {velocity[nodeId], gridDelta};

    const auto foundNeighbor = nodeLookup.find(linearIndex(neighbor));
    if (foundNeighbor != nodeLookup.end())
      return {velocity[foundNeighbor->second], gridDelta};

    const unsigned fi = direction * 2u + (offset == 1 ? 1u : 0u);
    const auto &intersection = node.faceBC[fi];
    if (intersection.type == Boundary::REACTION)
      return {reactionBoundaryVelocity(node.index), intersection.distance};
    if (intersection.type == Boundary::AMBIENT)
      return {freeSurfaceVelocityBoundary(node.index, direction, offset,
                                          intersection.distance,
                                          velocity[nodeId]),
              intersection.distance};
    if (intersection.type == Boundary::MASK)
      return {maskVelocityBoundary(node.index, velocity[nodeId]),
              intersection.distance};

    return {velocity[nodeId], gridDelta};
  }

  StencilPoint<T> currentPressureStencilPoint(std::size_t nodeId,
                                              unsigned direction,
                                              int offset) const {
    const auto &node = nodes[nodeId];
    IndexType neighbor = node.index;
    neighbor[direction] += offset;

    if (!inBounds(neighbor))
      return {node.pressure, gridDelta};

    const auto foundNeighbor = nodeLookup.find(linearIndex(neighbor));
    if (foundNeighbor != nodeLookup.end())
      return {nodes[foundNeighbor->second].pressure, gridDelta};

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

    const auto foundNeighbor = nodeLookup.find(linearIndex(neighbor));
    if (foundNeighbor != nodeLookup.end())
      return {nodes[foundNeighbor->second].velocity, gridDelta};

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

          if (nodeLookup.find(linearIndex(neighbor)) != nodeLookup.end())
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

    const auto found = nodeLookup.find(linearIndex(index));
    if (found == nodeLookup.end())
      return {};

    std::array<T, 9> deviatoric = nodes[found->second].stressTensor;
    for (unsigned i = 0; i < 3; ++i)
      deviatoric[tensorIndex(i, i)] += nodes[found->second].pressure;
    return deviatoric;
  }

  T pressureAt(const IndexType &index) const {
    if (!inBounds(index))
      return deformationParameters.ambientPressure;

    const auto found = nodeLookup.find(linearIndex(index));
    if (found == nodeLookup.end())
      return deformationParameters.ambientPressure;
    return nodes[found->second].pressure;
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
      normal[D - 1] = 1.;
      return normal;
    }

    norm = std::sqrt(norm);
    for (unsigned i = 0; i < D; ++i)
      normal[i] /= norm;
    return normal;
  }

  void computeDiagnostics() {
    for (auto &node : nodes) {
      node.strainTrace = divergenceAt(node.index);
    }
  }

  void computeStressTensors() {
    std::unordered_map<std::size_t, std::array<T, 9>> nextHistory;
    const T relaxationTime = effectiveStressRelaxationTime();
    const T decay =
        (relaxationTime <= std::numeric_limits<T>::epsilon())
            ? T(0)
            : std::exp(-deformationParameters.stressTimeStep / relaxationTime);

    for (auto &node : nodes) {
      node.strainRateTensor = strainRateTensorAt(node.index);
      const auto deviatoricRate =
          deviatoricTensor(node.strainRateTensor, node.strainTrace);
      const auto previousStress =
          previousDeviatoricStress(node.index);

      std::array<T, 9> deviatoricStress{};
      for (unsigned i = 0; i < 9; ++i) {
        const T viscousStress =
            T(2) * deformationParameters.viscosity * deviatoricRate[i];
        deviatoricStress[i] =
            decay * previousStress[i] + (T(1) - decay) * viscousStress;
      }

      node.stressTensor = deviatoricStress;
      for (unsigned i = 0; i < 3; ++i)
        node.stressTensor[tensorIndex(i, i)] -= node.pressure;

      node.vonMisesStress = vonMisesFromDeviatoric(deviatoricStress);
      nextHistory[detail::gridIndexHash<D>(node.index)] = deviatoricStress;
    }

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
    const auto found = nodeLookup.find(linearIndex(index));
    if (found == nodeLookup.end())
      return 0.;

    const auto plus = currentVelocityStencilPoint(found->second, direction, 1);
    const auto minus = currentVelocityStencilPoint(found->second, direction, -1);
    return firstDerivative(minus.value[component],
                           nodes[found->second].velocity[component],
                           plus.value[component], minus.distance,
                           plus.distance);
  }

  T pressureDerivative(const IndexType &index, unsigned direction) const {
    const auto found = nodeLookup.find(linearIndex(index));
    if (found == nodeLookup.end())
      return 0.;

    const auto plus = currentPressureStencilPoint(found->second, direction, 1);
    const auto minus = currentPressureStencilPoint(found->second, direction, -1);
    return firstDerivative(minus.value, nodes[found->second].pressure,
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
      return {Boundary::AMBIENT, crossingDistance(ambientInside, ambientOutside)};
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

        if (nodeLookup.find(linearIndex(neighbor)) != nodeLookup.end())
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

  T crossingDistance(T insidePhi, T outsidePhi) const {
    return detail::levelSetCrossingDistance(
        insidePhi, outsidePhi,
        deformationParameters.minMechanicsBoundaryDistance, gridDelta);
  }

  // Vec3D<T> getNearbyVelocity(const IndexType &index) const {
  //   const auto nearby = findNearbyNode(index);
  //   if (nearby == std::numeric_limits<std::size_t>::max())
  //     return {0., 0., 0.};
  //   return nodes[nearby].velocity;
  // }


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
