#pragma once

#include <lsOxidationSolverBase.hpp>
#include <lsVelocityField.hpp>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include <omp.h>

namespace viennals {

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

  struct Node {
    IndexType index;
    T concentration = 0.;
    Vec3D<T> siNormal = {0., 0., 0.}; // unit outward normal of Si surface (into oxide)
    // Per-face boundary cache: eliminates goToIndices() from the parallel solve loop.
    struct FaceBC { Boundary type = Boundary::NONE; T distance = 1.; };
    std::array<FaceBC, 2 * D> faceBC{};
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
  T maxScalarVelocity_ = 0.;
  bool solved = false;
  bool nodesDirty_ = true;  // true → rebuild grid/nodes on next apply()
  bool useRequestedBounds = false;
  std::unordered_map<std::size_t, T> pressureLookup;
  std::vector<Node> nodes;

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
    pressureLookup[detail::gridIndexHash<D>(index)] = pressure;
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
    auto it = nodeLookup.find(linearIndex(index));
    if (it == nodeLookup.end()) {
      const auto nearby = findNearbyNode(index);
      if (nearby == std::numeric_limits<std::size_t>::max())
        return 0.;
      return nodes[nearby].concentration;
    }
    return nodes[it->second].concentration;
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
  std::size_t getNumberOfSolutionNodes() const { return nodes.size(); }

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
        Node newNode{index, parameters.equilibriumConcentration};
        if (parameters.reactionRateRatio111 != T(1))
          newNode.siNormal = computeSiNormal(index, reactionIt);
        nodes.push_back(newNode);
      }

      if (!increment(index))
        break;
    }

    // Precompute per-face boundary intersections. Level-set positions are fixed
    // during the inner solve, so one computation per apply() suffices and the
    // hot loop becomes iterator-free (enabling parallelism).
    ConstSparseIterator faceReactionIt(reactionInterface->getDomain());
    ConstSparseIterator faceAmbientIt(ambientInterface->getDomain());
    auto faceMaskIt = makeMaskIterator();
    for (auto &node : nodes) {
      for (unsigned dir = 0; dir < D; ++dir) {
        for (int off : {-1, 1}) {
          const unsigned fi = dir * 2u + (off == 1 ? 1u : 0u);
          IndexType nb = node.index;
          nb[dir] += off;
          if (!inBounds(nb) || nodeLookup.count(linearIndex(nb))) {
            node.faceBC[fi] = {};
            continue;
          }
          const auto bc = classifyBoundary(faceReactionIt, faceAmbientIt,
                                           faceMaskIt, node.index, nb);
          node.faceBC[fi] = {bc.first, bc.second};
        }
      }
    }
  }

  void solveDiffusion() {
    iterations = 0;
    residual = 0.;
    if (nodes.empty())
      return;

    std::vector<T> previous(nodes.size(), 0.);
    std::vector<T> next = previous;

    // Face BCs are precomputed in buildNodes(); no per-iteration iterators needed.
    for (; iterations < parameters.maxIterations; ++iterations) {
      residual = 0.;
#pragma omp parallel for schedule(static) reduction(max : residual)
      for (std::size_t nodeId = 0; nodeId < nodes.size(); ++nodeId) {
        const auto &node = nodes[nodeId];
        T rightHandSide = 0.;
        T diagonal = 0.;

        const T D_eff = getEffectiveDiffusionCoefficient(node.index);
        for (unsigned direction = 0; direction < D; ++direction) {
          const auto negativeSide =
              makeStencilSide(node, previous, direction, -1, D_eff);
          const auto positiveSide =
              makeStencilSide(node, previous, direction, 1, D_eff);
          addAxisContribution(rightHandSide, diagonal, negativeSide,
                              positiveSide, D_eff);
        }

        const T updated =
            (diagonal <= std::numeric_limits<T>::epsilon())
                ? previous[nodeId]
                : rightHandSide / diagonal;
        next[nodeId] = parameters.relaxation * updated +
                       (T(1) - parameters.relaxation) * previous[nodeId];
        residual = std::max(residual, std::abs(next[nodeId] - previous[nodeId]));
      }

      previous.swap(next);
      if (residual < parameters.tolerance)
        break;
    }

    for (std::size_t i = 0; i < nodes.size(); ++i)
      nodes[i].concentration = previous[i];
  }

  // Uses precomputed node.faceBC — no HRLE iterator access, safe for parallel execution.
  StencilSide makeStencilSide(const Node &node, const std::vector<T> &previous,
                              unsigned direction, int offset, T diffusion) const {
    IndexType neighbor = node.index;
    neighbor[direction] += offset;

    if (!inBounds(neighbor))
      return zeroFluxSide();

    const auto foundNeighbor = nodeLookup.find(linearIndex(neighbor));
    if (foundNeighbor != nodeLookup.end())
      return {gridDelta, 0., previous[foundNeighbor->second]};

    const unsigned fi = direction * 2u + (offset == 1 ? 1u : 0u);
    const auto &face = node.faceBC[fi];
    if (face.type == Boundary::REACTION)
      return reactionBoundarySide(node.index, face.distance, diffusion);
    if (face.type == Boundary::AMBIENT)
      return ambientBoundarySide(face.distance, diffusion);
    if (face.type == Boundary::MASK)
      return maskBoundarySide(face.distance, diffusion);
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

    const auto direct = nodeLookup.find(linearIndex(index));
    if (direct != nodeLookup.end()) {
      const auto sample =
          reactionBoundarySampleFromNode(reactionIt, nodes[direct->second]);
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
          const auto found = nodeLookup.find(linearIndex(candidate));
          if (found != nodeLookup.end() && distance2 < bestDistance2) {
            const auto sample =
                reactionBoundarySampleFromNode(reactionIt, nodes[found->second]);
            if (sample.found) {
              bestDistance2 = distance2;
              bestNode = found->second;
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
      const T exponent = stressExponent(pressure,
                                        parameters.reactionActivationVolume);
      rate *= stressFactor(exponent);
    }

    if (parameters.reactionRateRatio111 != T(1)) {
      const auto nodeIt = nodeLookup.find(linearIndex(index));
      if (nodeIt != nodeLookup.end()) {
        const auto &normal = nodes[nodeIt->second].siNormal;
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
    if (exponent <= std::log(minStressFactor))
      return minStressFactor;
    if (exponent >= std::log(maxStressFactor))
      return maxStressFactor;
    return std::exp(exponent);
  }

  Vec3D<T> computeSiNormal(const IndexType &index,
                           ConstSparseIterator &reactionIt) const {
    Vec3D<T> gradient{0., 0., 0.};
    for (unsigned d = 0; d < D; ++d) {
      IndexType plus = index, minus = index;
      plus[d] += 1;
      minus[d] -= 1;
      gradient[d] =
          (detail::clampLevelSetPhi(valueAt(reactionIt, plus)) -
           detail::clampLevelSetPhi(valueAt(reactionIt, minus))) /
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
