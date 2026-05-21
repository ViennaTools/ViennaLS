#pragma once

#include <lsDomain.hpp>
#include <lsVelocityField.hpp>

#include <hrleSparseIterator.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace viennals {

using namespace viennacore;

/// Parameters for the steady oxidant diffusion model used by
/// OxidationDiffusionVelocityField.
template <class T> struct OxidationParameters {
  T diffusionCoefficient = 1.;
  T reactionRate = 1.;
  T transferCoefficient = 1.;
  T equilibriumConcentration = 1.;
  T oxidantMoleculeDensity = 1.;
  T expansionCoefficient = 1.;
  T velocitySign = 1.;
  T stressCouplingCoefficient = 0.;
  T referencePressure = 0.;
  T minStressRateFactor = 0.01;
  T maxStressRateFactor = 100.;
  T maskTransferCoefficient = 0.;
  T maskConcentration = 0.;
  unsigned maxIterations = 10000;
  T tolerance = 1e-8;
  T relaxation = 1.;
  std::size_t maxGridPoints = 5000000;
  int material = -1;
};

/// Parameters for the Cartesian-grid oxide deformation model.
template <class T> struct OxidationDeformationParameters {
  T viscosity = 1.;
  T bulkModulus = 1.;
  T ambientPressure = 0.;
  T pressureRelaxation = 1.;
  T pressureTolerance = 1e-8;
  T pressureGradientScale = 1e-3;
  T freeSurfaceTractionScale = 1.;
  T substrateNormalStiffness = 0.;
  T minMechanicsBoundaryDistance = 0.05;
  T maskNormalStiffness = 0.;
  T maskVelocityScale = 0.;
  T maskPressure = 0.;
  T shearModulus = 0.;
  T stressRelaxationTime = 0.;
  T stressTimeStep = 1.;
  T freeSurfaceVelocityScale = 0.;
  T vectorVelocityScale = 1.;
  unsigned maxIterations = 10000;
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

template <class T> struct OxidationCouplingParameters {
  unsigned maxIterations = 10;
  T tolerance = 1e-6;
  T relaxation = 1.;
};

template <class T> struct OxidationMaskParameters {
  T youngModulus = 2.5e11;
  T poissonRatio = 0.27;
  T thickness = 0.1;
  T referenceThickness = 0.1;
  T velocityScale = 1.;
  T pressureVelocityScale = 0.;
  T referencePressure = 0.;
  T maxVelocity = std::numeric_limits<T>::max();
  T relaxation = 1.;
  T tolerance = 1e-8;
  T minBoundaryDistance = 0.05;
  unsigned maxIterations = 10000;
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
class OxidationDiffusionVelocityField final : public VelocityField<T> {
  using IndexType = viennahrle::Index<D>;
  using ConstSparseIterator =
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>;

  enum class Boundary { NONE, REACTION, AMBIENT, MASK };

  struct Node {
    IndexType index;
    T concentration = 0.;
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

  std::vector<Node> nodes;
  std::unordered_map<std::size_t, std::size_t> nodeLookup;
  IndexType minIndex{};
  IndexType maxIndex{};
  std::array<std::size_t, D> extents{};
  std::array<std::size_t, D> strides{};
  IndexType requestedMinIndex{};
  IndexType requestedMaxIndex{};
  T gridDelta = 1.;
  unsigned iterations = 0;
  T residual = std::numeric_limits<T>::max();
  bool solved = false;
  bool useRequestedBounds = false;
  std::unordered_map<std::size_t, T> pressureLookup;

public:
  OxidationDiffusionVelocityField() = default;

  OxidationDiffusionVelocityField(
      SmartPointer<Domain<T, D>> passedReactionInterface,
      SmartPointer<Domain<T, D>> passedAmbientInterface,
      OxidationParameters<T> passedParameters = {})
      : reactionInterface(passedReactionInterface),
        ambientInterface(passedAmbientInterface), parameters(passedParameters) {
  }

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<OxidationDiffusionVelocityField>::New(
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
    pressureLookup[pressureKey(index)] = pressure;
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
    solved = false;
  }

  void clearSolveBounds() {
    useRequestedBounds = false;
    solved = false;
  }

  void apply() {
    if (reactionInterface == nullptr || ambientInterface == nullptr) {
      Logger::getInstance()
          .addError("OxidationDiffusionVelocityField: Missing level-set "
                    "interface.")
          .print();
      return;
    }

    initialiseGrid();
    buildNodes();
    solveDiffusion();
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

    const T concentration = getConcentration(index);
    return parameters.velocitySign * getEffectiveReactionRate(index) * concentration /
           (parameters.oxidantMoleculeDensity * parameters.expansionCoefficient);
  }

  T getDissipationAlpha(int /*direction*/, int material,
                        const Vec3D<T> & /*centralDifferences*/) final {
    if (parameters.material >= 0 && material != parameters.material)
      return 0.;
    return std::abs(parameters.velocitySign) * parameters.reactionRate *
           parameters.equilibriumConcentration /
           (parameters.oxidantMoleculeDensity * parameters.expansionCoefficient);
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

  unsigned getIterations() const { return iterations; }
  T getResidual() const { return residual; }
  std::size_t getNumberOfSolutionNodes() const { return nodes.size(); }

private:
  void initialiseGrid() {
    auto &reactionGrid = reactionInterface->getGrid();
    auto &ambientGrid = ambientInterface->getGrid();
    gridDelta = reactionGrid.getGridDelta();

    if (std::abs(gridDelta - ambientGrid.getGridDelta()) >
            std::numeric_limits<T>::epsilon() ||
        (maskInterface != nullptr &&
         std::abs(gridDelta - maskInterface->getGrid().getGridDelta()) >
             std::numeric_limits<T>::epsilon())) {
      Logger::getInstance()
          .addError("OxidationDiffusionVelocityField: Interface grid deltas "
                    "must match.")
          .print();
      return;
    }

    minIndex = reactionGrid.getMinGridPoint();
    maxIndex = reactionGrid.getMaxGridPoint();
    std::size_t numGridPoints = 1;
    for (unsigned i = 0; i < D; ++i) {
      minIndex[i] = std::max(minIndex[i], ambientGrid.getMinGridPoint(i));
      maxIndex[i] = std::min(maxIndex[i], ambientGrid.getMaxGridPoint(i));
      if (maskInterface != nullptr) {
        minIndex[i] =
            std::max(minIndex[i], maskInterface->getGrid().getMinGridPoint(i));
        maxIndex[i] =
            std::min(maxIndex[i], maskInterface->getGrid().getMaxGridPoint(i));
      }
      if (useRequestedBounds) {
        minIndex[i] = std::max(minIndex[i], requestedMinIndex[i]);
        maxIndex[i] = std::min(maxIndex[i], requestedMaxIndex[i]);
      }
      extents[i] = static_cast<std::size_t>(maxIndex[i] - minIndex[i] + 1);
      numGridPoints *= extents[i];
      strides[i] = (i == 0) ? 1 : strides[i - 1] * extents[i - 1];
    }

    if (numGridPoints > parameters.maxGridPoints) {
      Logger::getInstance()
          .addError("OxidationDiffusionVelocityField: Cartesian solve region "
                    "exceeds maxGridPoints.")
          .print();
    }
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
        nodes.push_back({index, parameters.equilibriumConcentration});
      }

      if (!increment(index))
        break;
    }
  }

  void solveDiffusion() {
    iterations = 0;
    residual = 0.;
    if (nodes.empty())
      return;

    std::vector<T> previous(nodes.size(), parameters.equilibriumConcentration);
    std::vector<T> next = previous;

    ConstSparseIterator reactionIt(reactionInterface->getDomain());
    ConstSparseIterator ambientIt(ambientInterface->getDomain());
    auto maskIt = makeMaskIterator();

    for (; iterations < parameters.maxIterations; ++iterations) {
      residual = 0.;
      for (std::size_t nodeId = 0; nodeId < nodes.size(); ++nodeId) {
        const auto &node = nodes[nodeId];
        T rightHandSide = 0.;
        T diagonal = 0.;

        for (unsigned direction = 0; direction < D; ++direction) {
          const auto negativeSide = makeStencilSide(
              reactionIt, ambientIt, maskIt, previous, node.index, direction, -1);
          const auto positiveSide = makeStencilSide(
              reactionIt, ambientIt, maskIt, previous, node.index, direction, 1);
          addAxisContribution(rightHandSide, diagonal, negativeSide,
                              positiveSide);
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

  StencilSide makeStencilSide(ConstSparseIterator &reactionIt,
                              ConstSparseIterator &ambientIt,
                              ConstSparseIterator &maskIt,
                              const std::vector<T> &previous,
                              const IndexType &nodeIndex, unsigned direction,
                              int offset) const {
    IndexType neighbor = nodeIndex;
    neighbor[direction] += offset;

    if (!inBounds(neighbor))
      return zeroFluxSide();

    const auto foundNeighbor = nodeLookup.find(linearIndex(neighbor));
    if (foundNeighbor != nodeLookup.end())
      return {gridDelta, 0., previous[foundNeighbor->second]};

    const auto boundary =
        classifyBoundary(reactionIt, ambientIt, maskIt, nodeIndex, neighbor);
    if (boundary.first == Boundary::REACTION)
      return reactionBoundarySide(nodeIndex, boundary.second);
    if (boundary.first == Boundary::AMBIENT)
      return ambientBoundarySide(boundary.second);
    if (boundary.first == Boundary::MASK)
      return maskBoundarySide(boundary.second);
    return zeroFluxSide();
  }

  void addAxisContribution(T &rightHandSide, T &diagonal,
                           const StencilSide &negativeSide,
                           const StencilSide &positiveSide) const {
    const T distanceSum = negativeSide.distance + positiveSide.distance;
    if (distanceSum <= std::numeric_limits<T>::epsilon())
      return;

    addSideContribution(rightHandSide, diagonal, negativeSide, distanceSum);
    addSideContribution(rightHandSide, diagonal, positiveSide, distanceSum);
  }

  void addSideContribution(T &rightHandSide, T &diagonal,
                           const StencilSide &side, T distanceSum) const {
    if (side.distance <= std::numeric_limits<T>::epsilon())
      return;

    const T coefficient =
        T(2) * parameters.diffusionCoefficient / (side.distance * distanceSum);
    rightHandSide += coefficient * side.constant;
    diagonal += coefficient * (T(1) - side.nodeCoefficient);
  }

  StencilSide zeroFluxSide() const { return {gridDelta, 1., 0.}; }

  StencilSide reactionBoundarySide(const IndexType &nodeIndex, T distance) const {
    const T conductance = parameters.diffusionCoefficient / distance;
    const T reactionRate = getEffectiveReactionRate(nodeIndex);
    const T denominator = conductance + reactionRate;
    if (denominator <= std::numeric_limits<T>::epsilon())
      return zeroFluxSide();

    return {distance, conductance / denominator, 0.};
  }

  StencilSide ambientBoundarySide(T distance) const {
    const T conductance = parameters.diffusionCoefficient / distance;
    const T denominator = conductance + parameters.transferCoefficient;
    if (denominator <= std::numeric_limits<T>::epsilon())
      return zeroFluxSide();

    return {distance, conductance / denominator,
            parameters.transferCoefficient *
                parameters.equilibriumConcentration / denominator};
  }

  StencilSide maskBoundarySide(T distance) const {
    if (parameters.maskTransferCoefficient <= std::numeric_limits<T>::epsilon())
      return zeroFluxSide();

    const T conductance = parameters.diffusionCoefficient / distance;
    const T denominator = conductance + parameters.maskTransferCoefficient;
    if (denominator <= std::numeric_limits<T>::epsilon())
      return zeroFluxSide();

    return {distance, conductance / denominator,
            parameters.maskTransferCoefficient * parameters.maskConcentration /
                denominator};
  }

  T getEffectiveReactionRate(const IndexType &index) const {
    if (parameters.stressCouplingCoefficient == T(0))
      return parameters.reactionRate;

    T pressure = parameters.referencePressure;
    const auto foundPressure = pressureLookup.find(pressureKey(index));
    if (foundPressure != pressureLookup.end())
      pressure = foundPressure->second;

    const T exponent = -parameters.stressCouplingCoefficient *
                       (pressure - parameters.referencePressure);
    const T unclampedFactor = std::exp(exponent);
    const T factor = std::min(parameters.maxStressRateFactor,
                              std::max(parameters.minStressRateFactor,
                                       unclampedFactor));
    return parameters.reactionRate * factor;
  }

  std::size_t pressureKey(const IndexType &index) const {
    std::size_t seed = 0;
    for (unsigned i = 0; i < D; ++i) {
      const auto value = std::hash<long long>{}(static_cast<long long>(index[i]));
      seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
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

  bool crosses(T a, T b) const {
    return (a <= 0. && b >= 0.) || (a >= 0. && b <= 0.);
  }

  T crossingDistance(T insidePhi, T outsidePhi) const {
    const T denom = std::abs(insidePhi) + std::abs(outsidePhi);
    if (denom <= std::numeric_limits<T>::epsilon())
      return gridDelta;
    return std::max(T(1e-6) * gridDelta,
                    gridDelta * std::abs(insidePhi) / denom);
  }

  T valueAt(ConstSparseIterator &it, const IndexType &index) const {
    it.goToIndices(index);
    return it.getValue();
  }

  bool inBounds(const IndexType &index) const {
    for (unsigned i = 0; i < D; ++i) {
      if (index[i] < minIndex[i] || index[i] > maxIndex[i])
        return false;
    }
    return true;
  }

  std::size_t linearIndex(const IndexType &index) const {
    if (!inBounds(index))
      return std::numeric_limits<std::size_t>::max();

    std::size_t result = 0;
    for (unsigned i = 0; i < D; ++i) {
      result += static_cast<std::size_t>(index[i] - minIndex[i]) * strides[i];
    }
    return result;
  }

  bool increment(IndexType &index) const {
    for (unsigned i = 0; i < D; ++i) {
      if (index[i] < maxIndex[i]) {
        ++index[i];
        return true;
      }
      index[i] = minIndex[i];
    }
    return false;
  }

  std::size_t findNearbyNode(const IndexType &index) const {
    T bestDistance2 = std::numeric_limits<T>::max();
    std::size_t bestNode = std::numeric_limits<std::size_t>::max();

    IndexType offset{};
    offset.fill(-1);
    while (true) {
      IndexType candidate = index;
      T distance2 = 0.;
      for (unsigned i = 0; i < D; ++i) {
        candidate[i] += offset[i];
        distance2 += static_cast<T>(offset[i] * offset[i]);
      }

      if (inBounds(candidate)) {
        const auto found = nodeLookup.find(linearIndex(candidate));
        if (found != nodeLookup.end() && distance2 < bestDistance2) {
          bestDistance2 = distance2;
          bestNode = found->second;
        }
      }

      unsigned dim = 0;
      for (; dim < D; ++dim) {
        if (offset[dim] < 1) {
          ++offset[dim];
          break;
        }
        offset[dim] = -1;
      }
      if (dim == D)
        break;
    }

    return bestNode;
  }
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
/// stress terms, with approximate traction/free-surface and elastic-substrate
/// boundary conditions.
template <class T, int D>
class OxidationDeformationVelocityField final : public VelocityField<T> {
  using IndexType = viennahrle::Index<D>;
  using ConstSparseIterator =
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>;

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
  SmartPointer<OxidationDiffusionVelocityField<T, D>> diffusionField = nullptr;
  OxidationDeformationParameters<T> deformationParameters;
  OxidationParameters<T> oxidationParameters;
  int reactionSign = 1;
  int ambientSign = -1;
  int maskSign = 1;

  std::vector<Node> nodes;
  std::unordered_map<std::size_t, std::size_t> nodeLookup;
  IndexType minIndex{};
  IndexType maxIndex{};
  std::array<std::size_t, D> extents{};
  std::array<std::size_t, D> strides{};
  IndexType requestedMinIndex{};
  IndexType requestedMaxIndex{};
  T gridDelta = 1.;
  unsigned iterations = 0;
  T residual = std::numeric_limits<T>::max();
  T averageBoundaryExpansionVelocity = 0.;
  bool solved = false;
  bool useRequestedBounds = false;
  std::unordered_map<std::size_t, std::array<T, 9>> deviatoricStressHistory;

public:
  OxidationDeformationVelocityField() = default;

  OxidationDeformationVelocityField(
      SmartPointer<Domain<T, D>> passedReactionInterface,
      SmartPointer<Domain<T, D>> passedAmbientInterface,
      SmartPointer<OxidationDiffusionVelocityField<T, D>> passedDiffusionField,
      OxidationParameters<T> passedOxidationParameters,
      OxidationDeformationParameters<T> passedDeformationParameters = {})
      : reactionInterface(passedReactionInterface),
        ambientInterface(passedAmbientInterface),
        diffusionField(passedDiffusionField),
        deformationParameters(passedDeformationParameters),
        oxidationParameters(passedOxidationParameters) {}

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<OxidationDeformationVelocityField>::New(
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

  void setDiffusionField(
      SmartPointer<OxidationDiffusionVelocityField<T, D>> passedDiffusionField) {
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
          .addError("OxidationDeformationVelocityField: Missing interface or "
                    "diffusion field.")
          .print();
      return;
    }

    initialiseGrid();
    buildNodes();
    computeAverageBoundaryExpansionVelocity();
    solveVelocity();
    solveMechanics();
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

    return scaled(getVelocity(coordinate), deformationParameters.vectorVelocityScale);
  }

  T getScalarVelocity(const Vec3D<T> &coordinate, int material,
                      const Vec3D<T> &normalVector,
                      unsigned long /*pointId*/) final {
    if (!solved)
      apply();

    if (deformationParameters.material >= 0 &&
        material != deformationParameters.material)
      return 0.;

    return deformationParameters.freeSurfaceVelocityScale *
           localFreeSurfaceExpansionSpeed(coordinate, normalVector);
  }

  T getDissipationAlpha(int direction, int material,
                        const Vec3D<T> & /*centralDifferences*/) final {
    if (deformationParameters.material >= 0 &&
        material != deformationParameters.material)
      return 0.;

    T alpha = 0.;
    for (const auto &node : nodes)
      alpha = std::max(alpha, std::abs(node.velocity[direction]));
    return alpha;
  }

  Vec3D<T> getVelocity(const Vec3D<T> &coordinate) const {
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return getVelocity(index);
  }

  Vec3D<T> getVelocity(const IndexType &index) const {
    auto it = nodeLookup.find(linearIndex(index));
    if (it == nodeLookup.end())
      return getNearbyVelocity(index);
    return nodes[it->second].velocity;
  }

  T getPressure(const Vec3D<T> &coordinate) const {
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return getPressure(index);
  }

  T getPressure(const IndexType &index) const {
    auto it = nodeLookup.find(linearIndex(index));
    if (it == nodeLookup.end()) {
      const auto nearby = findNearbyNode(index);
      if (nearby == std::numeric_limits<std::size_t>::max())
        return 0.;
      return nodes[nearby].pressure;
    }
    return nodes[it->second].pressure;
  }

  T getStrainTrace(const Vec3D<T> &coordinate) const {
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return getStrainTrace(index);
  }

  T getStrainTrace(const IndexType &index) const {
    auto it = nodeLookup.find(linearIndex(index));
    if (it == nodeLookup.end()) {
      const auto nearby = findNearbyNode(index);
      if (nearby == std::numeric_limits<std::size_t>::max())
        return 0.;
      return nodes[nearby].strainTrace;
    }
    return nodes[it->second].strainTrace;
  }

  std::array<T, 9> getStrainRateTensor(const Vec3D<T> &coordinate) const {
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return getStrainRateTensor(index);
  }

  std::array<T, 9> getStrainRateTensor(const IndexType &index) const {
    auto it = nodeLookup.find(linearIndex(index));
    if (it == nodeLookup.end()) {
      const auto nearby = findNearbyNode(index);
      if (nearby == std::numeric_limits<std::size_t>::max())
        return {};
      return nodes[nearby].strainRateTensor;
    }
    return nodes[it->second].strainRateTensor;
  }

  std::array<T, 9> getStressTensor(const Vec3D<T> &coordinate) const {
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return getStressTensor(index);
  }

  std::array<T, 9> getStressTensor(const IndexType &index) const {
    auto it = nodeLookup.find(linearIndex(index));
    if (it == nodeLookup.end()) {
      const auto nearby = findNearbyNode(index);
      if (nearby == std::numeric_limits<std::size_t>::max())
        return {};
      return nodes[nearby].stressTensor;
    }
    return nodes[it->second].stressTensor;
  }

  T getVonMisesStress(const Vec3D<T> &coordinate) const {
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return getVonMisesStress(index);
  }

  T getVonMisesStress(const IndexType &index) const {
    auto it = nodeLookup.find(linearIndex(index));
    if (it == nodeLookup.end()) {
      const auto nearby = findNearbyNode(index);
      if (nearby == std::numeric_limits<std::size_t>::max())
        return 0.;
      return nodes[nearby].vonMisesStress;
    }
    return nodes[it->second].vonMisesStress;
  }

  unsigned getIterations() const { return iterations; }
  T getResidual() const { return residual; }
  std::size_t getNumberOfSolutionNodes() const { return nodes.size(); }
  T getAverageBoundaryExpansionVelocity() const {
    return averageBoundaryExpansionVelocity;
  }
  template <class Callback> void forEachSolutionNode(Callback callback) const {
    for (const auto &node : nodes)
      callback(node.index, node.pressure);
  }

private:
  void initialiseGrid() {
    auto &reactionGrid = reactionInterface->getGrid();
    auto &ambientGrid = ambientInterface->getGrid();
    gridDelta = reactionGrid.getGridDelta();

    if (std::abs(gridDelta - ambientGrid.getGridDelta()) >
            std::numeric_limits<T>::epsilon() ||
        (maskInterface != nullptr &&
         std::abs(gridDelta - maskInterface->getGrid().getGridDelta()) >
             std::numeric_limits<T>::epsilon())) {
      Logger::getInstance()
          .addError("OxidationDeformationVelocityField: Interface grid deltas "
                    "must match.")
          .print();
      return;
    }

    minIndex = reactionGrid.getMinGridPoint();
    maxIndex = reactionGrid.getMaxGridPoint();
    std::size_t numGridPoints = 1;
    for (unsigned i = 0; i < D; ++i) {
      minIndex[i] = std::max(minIndex[i], ambientGrid.getMinGridPoint(i));
      maxIndex[i] = std::min(maxIndex[i], ambientGrid.getMaxGridPoint(i));
      if (maskInterface != nullptr) {
        minIndex[i] =
            std::max(minIndex[i], maskInterface->getGrid().getMinGridPoint(i));
        maxIndex[i] =
            std::min(maxIndex[i], maskInterface->getGrid().getMaxGridPoint(i));
      }
      if (useRequestedBounds) {
        minIndex[i] = std::max(minIndex[i], requestedMinIndex[i]);
        maxIndex[i] = std::min(maxIndex[i], requestedMaxIndex[i]);
      }
      extents[i] = static_cast<std::size_t>(maxIndex[i] - minIndex[i] + 1);
      numGridPoints *= extents[i];
      strides[i] = (i == 0) ? 1 : strides[i - 1] * extents[i - 1];
    }

    if (numGridPoints > deformationParameters.maxGridPoints) {
      Logger::getInstance()
          .addError("OxidationDeformationVelocityField: Cartesian solve "
                    "region exceeds maxGridPoints.")
          .print();
    }
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
  }

  void solveVelocity() {
    iterations = 0;
    residual = 0.;
    if (nodes.empty())
      return;

    std::vector<Vec3D<T>> previous(nodes.size(), {0., 0., 0.});
    std::vector<Vec3D<T>> next = previous;

    ConstSparseIterator reactionIt(reactionInterface->getDomain());
    ConstSparseIterator ambientIt(ambientInterface->getDomain());
    auto maskIt = makeMaskIterator();

    for (; iterations < deformationParameters.maxIterations; ++iterations) {
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
              addTo(sum, previous[nodeId]);
              ++count;
              continue;
            }

            const auto foundNeighbor = nodeLookup.find(linearIndex(neighbor));
            if (foundNeighbor != nodeLookup.end()) {
              addTo(sum, previous[foundNeighbor->second]);
              ++count;
              continue;
            }

            const auto boundary =
                classifyBoundary(reactionIt, ambientIt, maskIt, node.index,
                                 neighbor);
            if (boundary == Boundary::REACTION) {
              addTo(sum, reactionBoundaryVelocity(node.index));
            } else if (boundary == Boundary::MASK) {
              addTo(sum, maskVelocityBoundary(previous[nodeId]));
            } else {
              addTo(sum, previous[nodeId]);
            }
            ++count;
          }
        }

        const Vec3D<T> updated =
            (count == 0) ? previous[nodeId] : scaled(sum, T(1) / count);
        next[nodeId] =
            add(scaled(updated, deformationParameters.relaxation),
                scaled(previous[nodeId],
                       T(1) - deformationParameters.relaxation));
        Vec3D<T> delta = subtract(next[nodeId], previous[nodeId]);
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

    ConstSparseIterator reactionIt(reactionInterface->getDomain());
    ConstSparseIterator ambientIt(ambientInterface->getDomain());
    auto maskIt = makeMaskIterator();

    T pressureResidual = 0.;
    for (unsigned iteration = 0;
         iteration < deformationParameters.pressureIterations; ++iteration) {
      pressureResidual = 0.;

      for (std::size_t nodeId = 0; nodeId < nodes.size(); ++nodeId) {
        const auto &node = nodes[nodeId];
        if (touchesBoundary(reactionIt, ambientIt, maskIt, node.index,
                            Boundary::AMBIENT)) {
          next[nodeId] = ambientBoundaryPressure[nodeId];
          continue;
        }

        T pressureSum = 0.;
        T centerCoefficient = 0.;
        for (unsigned direction = 0; direction < D; ++direction) {
          const auto plus = pressureStencilPoint(
              reactionIt, ambientIt, maskIt, previous, ambientBoundaryPressure,
              nodeId, direction, 1);
          const auto minus = pressureStencilPoint(
              reactionIt, ambientIt, maskIt, previous, ambientBoundaryPressure,
              nodeId, direction, -1);
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

    ConstSparseIterator reactionIt(reactionInterface->getDomain());
    ConstSparseIterator ambientIt(ambientInterface->getDomain());
    auto maskIt = makeMaskIterator();

    T velocityResidual = 0.;
    for (unsigned iteration = 0;
         iteration < deformationParameters.stokesIterations; ++iteration) {
      velocityResidual = 0.;

      for (std::size_t nodeId = 0; nodeId < nodes.size(); ++nodeId) {
        const auto &node = nodes[nodeId];
        Vec3D<T> sum{0., 0., 0.};
        T centerCoefficient = 0.;

        for (unsigned direction = 0; direction < D; ++direction) {
          const auto plus = velocityStencilPoint(reactionIt, ambientIt, maskIt,
                                                 previous, nodeId, direction, 1);
          const auto minus = velocityStencilPoint(reactionIt, ambientIt, maskIt,
                                                  previous, nodeId, direction, -1);
          const T plusCoefficient =
              T(2) / (plus.distance * (plus.distance + minus.distance));
          const T minusCoefficient =
              T(2) / (minus.distance * (plus.distance + minus.distance));
          addTo(sum, scaled(plus.value, plusCoefficient));
          addTo(sum, scaled(minus.value, minusCoefficient));
          centerCoefficient += plusCoefficient + minusCoefficient;
        }

        if (centerCoefficient <= std::numeric_limits<T>::epsilon())
          continue;

        const auto forcing = momentumForcing(node.index);
        Vec3D<T> updated = sum;
        for (unsigned component = 0; component < D; ++component) {
          updated[component] =
              (sum[component] -
               deformationParameters.pressureGradientScale * forcing[component] /
                   deformationParameters.viscosity) /
              centerCoefficient;
        }

        next[nodeId] =
            add(scaled(updated, deformationParameters.relaxation),
                scaled(previous[nodeId],
                       T(1) - deformationParameters.relaxation));

        const Vec3D<T> delta = subtract(next[nodeId], previous[nodeId]);
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

  StencilPoint<T> pressureStencilPoint(ConstSparseIterator &reactionIt,
                                       ConstSparseIterator &ambientIt,
                                       ConstSparseIterator &maskIt,
                                       const std::vector<T> &pressure,
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

    const auto intersection =
        boundaryIntersection(reactionIt, ambientIt, maskIt, node.index, neighbor);
    if (intersection.boundary == Boundary::AMBIENT)
      return {ambientBoundaryPressure[nodeId], intersection.distance};
    if (intersection.boundary == Boundary::REACTION)
      return {reactionPressureBoundary(node.index, pressure[nodeId]),
              intersection.distance};
    if (intersection.boundary == Boundary::MASK)
      return {maskPressureBoundary(node.index, pressure[nodeId]),
              intersection.distance};

    return {pressure[nodeId], gridDelta};
  }

  StencilPoint<Vec3D<T>> velocityStencilPoint(
      ConstSparseIterator &reactionIt, ConstSparseIterator &ambientIt,
      ConstSparseIterator &maskIt,
      const std::vector<Vec3D<T>> &velocity, std::size_t nodeId,
      unsigned direction, int offset) const {
    const auto &node = nodes[nodeId];
    IndexType neighbor = node.index;
    neighbor[direction] += offset;

    if (!inBounds(neighbor))
      return {velocity[nodeId], gridDelta};

    const auto foundNeighbor = nodeLookup.find(linearIndex(neighbor));
    if (foundNeighbor != nodeLookup.end())
      return {velocity[foundNeighbor->second], gridDelta};

    const auto intersection =
        boundaryIntersection(reactionIt, ambientIt, maskIt, node.index, neighbor);
    if (intersection.boundary == Boundary::REACTION)
      return {reactionBoundaryVelocity(node.index), intersection.distance};
    if (intersection.boundary == Boundary::AMBIENT)
      return {freeSurfaceVelocityBoundary(node.index, direction, offset,
                                          intersection.distance,
                                          velocity[nodeId]),
              intersection.distance};
    if (intersection.boundary == Boundary::MASK)
      return {maskVelocityBoundary(velocity[nodeId]), intersection.distance};

    return {velocity[nodeId], gridDelta};
  }

  StencilPoint<T> currentPressureStencilPoint(ConstSparseIterator &reactionIt,
                                              ConstSparseIterator &ambientIt,
                                              ConstSparseIterator &maskIt,
                                              std::size_t nodeId,
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

    const auto intersection =
        boundaryIntersection(reactionIt, ambientIt, maskIt, node.index, neighbor);
    if (intersection.boundary == Boundary::AMBIENT)
      return {freeSurfacePressureBoundary(node.index), intersection.distance};
    if (intersection.boundary == Boundary::REACTION)
      return {reactionPressureBoundary(node.index, node.pressure),
              intersection.distance};
    if (intersection.boundary == Boundary::MASK)
      return {maskPressureBoundary(node.index, node.pressure),
              intersection.distance};

    return {node.pressure, gridDelta};
  }

  StencilPoint<Vec3D<T>>
  currentVelocityStencilPoint(ConstSparseIterator &reactionIt,
                              ConstSparseIterator &ambientIt,
                              ConstSparseIterator &maskIt,
                              std::size_t nodeId, unsigned direction,
                              int offset) const {
    const auto &node = nodes[nodeId];
    IndexType neighbor = node.index;
    neighbor[direction] += offset;

    if (!inBounds(neighbor))
      return {node.velocity, gridDelta};

    const auto foundNeighbor = nodeLookup.find(linearIndex(neighbor));
    if (foundNeighbor != nodeLookup.end())
      return {nodes[foundNeighbor->second].velocity, gridDelta};

    const auto intersection =
        boundaryIntersection(reactionIt, ambientIt, maskIt, node.index, neighbor);
    if (intersection.boundary == Boundary::REACTION)
      return {reactionBoundaryVelocity(node.index), intersection.distance};
    if (intersection.boundary == Boundary::AMBIENT)
      return {freeSurfaceVelocityBoundary(node.index, direction, offset,
                                          intersection.distance, node.velocity),
              intersection.distance};
    if (intersection.boundary == Boundary::MASK)
      return {maskVelocityBoundary(node.velocity), intersection.distance};

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
    if (deformationParameters.freeSurfaceTractionScale == T(0))
      return deformationParameters.ambientPressure;

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
           deformationParameters.freeSurfaceTractionScale *
               normalStress(deviatoricStress, normal);
  }

  T reactionPressureBoundary(const IndexType &index,
                             T fallbackPressure) const {
    if (deformationParameters.substrateNormalStiffness <= T(0))
      return fallbackPressure;

    const auto normal = interfaceNormal(index, Boundary::REACTION);
    const auto boundaryVelocity = reactionBoundaryVelocity(index);
    T normalVelocity = 0.;
    for (unsigned i = 0; i < D; ++i)
      normalVelocity += boundaryVelocity[i] * normal[i];

    return fallbackPressure +
           deformationParameters.substrateNormalStiffness *
               deformationParameters.stressTimeStep * normalVelocity;
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
          deformationParameters.freeSurfaceTractionScale * normalTraction *
          normal[direction] /
          std::max(deformationParameters.viscosity,
                   std::numeric_limits<T>::epsilon());
      boundaryVelocity[component] +=
          static_cast<T>(offset) * distance * faceDerivative;
    }

    return boundaryVelocity;
  }

  T maskPressureBoundary(const IndexType &index, T fallbackPressure) const {
    if (deformationParameters.maskNormalStiffness <= T(0))
      return deformationParameters.maskPressure;

    const auto normal = interfaceNormal(index, Boundary::MASK);
    const auto velocity = getVelocity(index);
    T normalVelocity = 0.;
    for (unsigned i = 0; i < D; ++i)
      normalVelocity += velocity[i] * normal[i];

    return fallbackPressure +
           deformationParameters.maskNormalStiffness *
               deformationParameters.stressTimeStep * normalVelocity;
  }

  Vec3D<T> maskVelocityBoundary(const Vec3D<T> &interiorVelocity) const {
    return scaled(interiorVelocity, deformationParameters.maskVelocityScale);
  }

  void computeAverageBoundaryExpansionVelocity() {
    averageBoundaryExpansionVelocity = 0.;
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
        averageBoundaryExpansionVelocity +=
            (oxidationParameters.expansionCoefficient - T(1)) *
            std::abs(diffusionField->getScalarVelocity(coordinate, 0,
                                                       {0., 0., 0.}, 0));
        ++count;
      }
    }

    if (count > 0)
      averageBoundaryExpansionVelocity /= static_cast<T>(count);
  }

  Vec3D<T> reactionBoundaryVelocity(const IndexType &index) const {
    Vec3D<T> coordinate{0., 0., 0.};
    for (unsigned i = 0; i < D; ++i)
      coordinate[i] = index[i] * gridDelta;

    const T expansionVelocity = localExpansionSpeed(coordinate);
    return scaled(reactionNormal(index), reactionSign * expansionVelocity);
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

  T localFreeSurfaceExpansionSpeed(const Vec3D<T> &coordinate,
                                   const Vec3D<T> &normalVector) const {
    ConstSparseIterator reactionIt(reactionInterface->getDomain());
    const auto reactionCoordinate =
        projectToReactionInterface(reactionIt, coordinate, normalVector);

    Vec3D<T> sampleCoordinate = reactionCoordinate;
    Vec3D<T> oxideDirection = subtract(coordinate, reactionCoordinate);
    T directionNorm = 0.;
    for (unsigned i = 0; i < D; ++i)
      directionNorm += oxideDirection[i] * oxideDirection[i];
    directionNorm = std::sqrt(directionNorm);
    if (directionNorm > std::numeric_limits<T>::epsilon()) {
      for (unsigned i = 0; i < D; ++i)
        sampleCoordinate[i] += oxideDirection[i] / directionNorm * gridDelta;
    }

    return (oxidationParameters.expansionCoefficient - T(1)) *
           std::abs(diffusionField->getScalarVelocity(sampleCoordinate, 0,
                                                      {0., 0., 0.}, 0));
  }

  Vec3D<T> projectToReactionInterface(ConstSparseIterator &reactionIt,
                                      const Vec3D<T> &coordinate,
                                      const Vec3D<T> &normalVector) const {
    Vec3D<T> normal = normalVector;
    T normalNorm = 0.;
    for (unsigned i = 0; i < D; ++i)
      normalNorm += normal[i] * normal[i];
    normalNorm = std::sqrt(normalNorm);
    if (normalNorm <= std::numeric_limits<T>::epsilon())
      return coordinate;
    for (unsigned i = 0; i < D; ++i)
      normal[i] /= normalNorm;

    const T originPhi = valueAtCoordinate(reactionIt, coordinate);
    T bestDistance = std::numeric_limits<T>::max();
    Vec3D<T> bestCoordinate = coordinate;

    for (T sign : {T(-1), T(1)}) {
      Vec3D<T> previousCoordinate = coordinate;
      T previousPhi = originPhi;

      for (unsigned step = 1; step <= 200; ++step) {
        Vec3D<T> currentCoordinate = coordinate;
        for (unsigned i = 0; i < D; ++i)
          currentCoordinate[i] += sign * normal[i] * gridDelta *
                                  static_cast<T>(step);

        const T currentPhi = valueAtCoordinate(reactionIt, currentCoordinate);
        if (crosses(previousPhi, currentPhi)) {
          const T denominator = std::abs(previousPhi) + std::abs(currentPhi);
          const T fraction =
              (denominator <= std::numeric_limits<T>::epsilon())
                  ? T(0)
                  : std::abs(previousPhi) / denominator;
          Vec3D<T> crossingCoordinate = previousCoordinate;
          for (unsigned i = 0; i < D; ++i) {
            crossingCoordinate[i] +=
                (currentCoordinate[i] - previousCoordinate[i]) * fraction;
          }

          const T distance = gridDelta * static_cast<T>(step - 1) +
                             gridDelta * fraction;
          if (distance < bestDistance) {
            bestDistance = distance;
            bestCoordinate = crossingCoordinate;
          }
          break;
        }

        previousCoordinate = currentCoordinate;
        previousPhi = currentPhi;
      }
    }

    return bestCoordinate;
  }

  T valueAtCoordinate(ConstSparseIterator &reactionIt,
                      const Vec3D<T> &coordinate) const {
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return valueAt(reactionIt, index);
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
      normal[i] = valueAt(levelSetIt, pos) - valueAt(levelSetIt, neg);
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
      nextHistory[nodeKey(node.index)] = deviatoricStress;
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

    auto reactionIt = ConstSparseIterator(reactionInterface->getDomain());
    auto ambientIt = ConstSparseIterator(ambientInterface->getDomain());
    auto maskIt = makeMaskIterator();
    const auto plus = currentVelocityStencilPoint(reactionIt, ambientIt, maskIt,
                                                  found->second, direction, 1);
    const auto minus = currentVelocityStencilPoint(reactionIt, ambientIt, maskIt,
                                                   found->second, direction, -1);
    return firstDerivative(minus.value[component],
                           nodes[found->second].velocity[component],
                           plus.value[component], minus.distance,
                           plus.distance);
  }

  T pressureDerivative(const IndexType &index, unsigned direction) const {
    const auto found = nodeLookup.find(linearIndex(index));
    if (found == nodeLookup.end())
      return 0.;

    auto reactionIt = ConstSparseIterator(reactionInterface->getDomain());
    auto ambientIt = ConstSparseIterator(ambientInterface->getDomain());
    auto maskIt = makeMaskIterator();
    const auto plus = currentPressureStencilPoint(reactionIt, ambientIt, maskIt,
                                                  found->second, direction, 1);
    const auto minus = currentPressureStencilPoint(reactionIt, ambientIt, maskIt,
                                                   found->second, direction, -1);
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
    const auto found = deviatoricStressHistory.find(nodeKey(index));
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

  bool crosses(T a, T b) const {
    return (a <= 0. && b >= 0.) || (a >= 0. && b <= 0.);
  }

  T crossingDistance(T insidePhi, T outsidePhi) const {
    const T denom = std::abs(insidePhi) + std::abs(outsidePhi);
    if (denom <= std::numeric_limits<T>::epsilon())
      return gridDelta;
    return std::max(deformationParameters.minMechanicsBoundaryDistance * gridDelta,
                    gridDelta * std::abs(insidePhi) / denom);
  }

  T valueAt(ConstSparseIterator &it, const IndexType &index) const {
    it.goToIndices(index);
    return it.getValue();
  }

  bool inBounds(const IndexType &index) const {
    for (unsigned i = 0; i < D; ++i) {
      if (index[i] < minIndex[i] || index[i] > maxIndex[i])
        return false;
    }
    return true;
  }

  std::size_t linearIndex(const IndexType &index) const {
    if (!inBounds(index))
      return std::numeric_limits<std::size_t>::max();

    std::size_t result = 0;
    for (unsigned i = 0; i < D; ++i)
      result += static_cast<std::size_t>(index[i] - minIndex[i]) * strides[i];
    return result;
  }

  bool increment(IndexType &index) const {
    for (unsigned i = 0; i < D; ++i) {
      if (index[i] < maxIndex[i]) {
        ++index[i];
        return true;
      }
      index[i] = minIndex[i];
    }
    return false;
  }

  Vec3D<T> getNearbyVelocity(const IndexType &index) const {
    const auto nearby = findNearbyNode(index);
    if (nearby == std::numeric_limits<std::size_t>::max())
      return {0., 0., 0.};
    return nodes[nearby].velocity;
  }

  std::size_t findNearbyNode(const IndexType &index) const {
    T bestDistance2 = std::numeric_limits<T>::max();
    std::size_t bestNode = std::numeric_limits<std::size_t>::max();

    IndexType offset{};
    offset.fill(-1);
    while (true) {
      IndexType candidate = index;
      T distance2 = 0.;
      for (unsigned i = 0; i < D; ++i) {
        candidate[i] += offset[i];
        distance2 += static_cast<T>(offset[i] * offset[i]);
      }

      if (inBounds(candidate)) {
        const auto found = nodeLookup.find(linearIndex(candidate));
        if (found != nodeLookup.end() && distance2 < bestDistance2) {
          bestDistance2 = distance2;
          bestNode = found->second;
        }
      }

      unsigned dim = 0;
      for (; dim < D; ++dim) {
        if (offset[dim] < 1) {
          ++offset[dim];
          break;
        }
        offset[dim] = -1;
      }
      if (dim == D)
        break;
    }

    return bestNode;
  }

  std::size_t nodeKey(const IndexType &index) const {
    std::size_t seed = 0;
    for (unsigned i = 0; i < D; ++i) {
      const auto value = std::hash<long long>{}(static_cast<long long>(index[i]));
      seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }

  static void addTo(Vec3D<T> &target, const Vec3D<T> &source) {
    for (unsigned i = 0; i < 3; ++i)
      target[i] += source[i];
  }

  static Vec3D<T> scaled(const Vec3D<T> &source, T factor) {
    Vec3D<T> result{0., 0., 0.};
    for (unsigned i = 0; i < 3; ++i)
      result[i] = source[i] * factor;
    return result;
  }

  static Vec3D<T> add(const Vec3D<T> &a, const Vec3D<T> &b) {
    Vec3D<T> result{0., 0., 0.};
    for (unsigned i = 0; i < 3; ++i)
      result[i] = a[i] + b[i];
    return result;
  }

  static Vec3D<T> subtract(const Vec3D<T> &a, const Vec3D<T> &b) {
    Vec3D<T> result{0., 0., 0.};
    for (unsigned i = 0; i < 3; ++i)
      result[i] = a[i] - b[i];
    return result;
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

/// Vector velocity field for a compliant oxidation mask driven by the solved
/// oxide mechanics field. A Cartesian elasticity solve is built inside the mask
/// level set, oxide-contact motion is applied as a boundary displacement rate,
/// and the resulting vector field is used to move the entire mask body.
template <class T, int D>
class OxidationMaskBendingVelocityField final : public VelocityField<T> {
  using IndexType = viennahrle::Index<D>;
  using ConstSparseIterator =
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>;

  struct Node {
    IndexType index;
    Vec3D<T> velocity{0., 0., 0.};
    bool contact = false;
  };

  SmartPointer<OxidationDeformationVelocityField<T, D>> deformationField =
      nullptr;
  SmartPointer<Domain<T, D>> maskInterface = nullptr;
  OxidationMaskParameters<T> parameters;
  int maskSign = 1;

  std::vector<Node> nodes;
  std::unordered_map<std::size_t, std::size_t> nodeLookup;
  IndexType minIndex{};
  IndexType maxIndex{};
  IndexType requestedMinIndex{};
  IndexType requestedMaxIndex{};
  std::array<std::size_t, D> extents{};
  std::array<std::size_t, D> strides{};
  T gridDelta = 1.;
  bool solved = false;
  bool useRequestedBounds = false;
  T residual = std::numeric_limits<T>::max();
  unsigned iterations = 0;
  std::size_t contactNodes = 0;

public:
  OxidationMaskBendingVelocityField() = default;

  OxidationMaskBendingVelocityField(
      SmartPointer<OxidationDeformationVelocityField<T, D>> passedDeformation,
      OxidationMaskParameters<T> passedParameters = {})
      : deformationField(passedDeformation), parameters(passedParameters) {}

  OxidationMaskBendingVelocityField(
      SmartPointer<OxidationDeformationVelocityField<T, D>> passedDeformation,
      SmartPointer<Domain<T, D>> passedMaskInterface,
      OxidationMaskParameters<T> passedParameters = {},
      int passedMaskSign = 1)
      : deformationField(passedDeformation), maskInterface(passedMaskInterface),
        parameters(passedParameters),
        maskSign((passedMaskSign < 0) ? -1 : 1) {}

  static SmartPointer<OxidationMaskBendingVelocityField>
  New(SmartPointer<OxidationDeformationVelocityField<T, D>> passedDeformation,
      OxidationMaskParameters<T> passedParameters = {}) {
    return SmartPointer<OxidationMaskBendingVelocityField>::New(
        passedDeformation, passedParameters);
  }

  static SmartPointer<OxidationMaskBendingVelocityField>
  New(SmartPointer<OxidationDeformationVelocityField<T, D>> passedDeformation,
      SmartPointer<Domain<T, D>> passedMaskInterface,
      OxidationMaskParameters<T> passedParameters = {},
      int passedMaskSign = 1) {
    return SmartPointer<OxidationMaskBendingVelocityField>::New(
        passedDeformation, passedMaskInterface, passedParameters,
        passedMaskSign);
  }

  void setMaskInterface(SmartPointer<Domain<T, D>> passedMaskInterface,
                        int passedMaskSign = 1) {
    maskInterface = passedMaskInterface;
    maskSign = (passedMaskSign < 0) ? -1 : 1;
    solved = false;
  }

  void setParameters(OxidationMaskParameters<T> passedParameters) {
    parameters = passedParameters;
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

  OxidationMaskParameters<T> getParameters() const { return parameters; }
  unsigned getIterations() const { return iterations; }
  T getResidual() const { return residual; }
  std::size_t getNumberOfSolutionNodes() const { return nodes.size(); }
  std::size_t getNumberOfContactNodes() const { return contactNodes; }

  void apply() {
    if (deformationField == nullptr || maskInterface == nullptr) {
      solved = true;
      return;
    }

    initialiseGrid();
    buildNodes();
    solveElasticVelocity();
    solved = true;
  }

  Vec3D<T> getVectorVelocity(const Vec3D<T> &coordinate, int material,
                             const Vec3D<T> &normalVector,
                             unsigned long /*pointId*/) final {
    if (deformationField == nullptr)
      return {0., 0., 0.};

    if (parameters.material >= 0 && material != parameters.material)
      return {0., 0., 0.};

    if (!solved)
      apply();

    if (maskInterface == nullptr || nodes.empty()) {
      Vec3D<T> velocity =
          scaled(deformationField->getVelocity(coordinate),
                 elasticCompliance() * parameters.velocityScale);
      const T pressure =
          std::abs(deformationField->getPressure(coordinate) -
                   parameters.referencePressure);
      if (parameters.pressureVelocityScale != T(0) && pressure > T(0))
        velocity = add(velocity,
                       scaled(normalVector, elasticCompliance() *
                                                parameters.pressureVelocityScale *
                                                pressure));
      return clampMagnitude(velocity);
    }

    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return clampMagnitude(getVelocity(index));
  }

  T getDissipationAlpha(int /*direction*/, int /*material*/,
                        const Vec3D<T> & /*centralDifferences*/) final {
    return parameters.maxVelocity == std::numeric_limits<T>::max()
               ? T(0)
               : parameters.maxVelocity;
  }

private:
  void initialiseGrid() {
    auto &grid = maskInterface->getGrid();
    gridDelta = grid.getGridDelta();

    minIndex = grid.getMinGridPoint();
    maxIndex = grid.getMaxGridPoint();
    std::size_t numGridPoints = 1;
    for (unsigned i = 0; i < D; ++i) {
      if (useRequestedBounds) {
        minIndex[i] = std::max(minIndex[i], requestedMinIndex[i]);
        maxIndex[i] = std::min(maxIndex[i], requestedMaxIndex[i]);
      }
      extents[i] = static_cast<std::size_t>(maxIndex[i] - minIndex[i] + 1);
      numGridPoints *= extents[i];
      strides[i] = (i == 0) ? 1 : strides[i - 1] * extents[i - 1];
    }

    if (numGridPoints > parameters.maxGridPoints) {
      Logger::getInstance()
          .addError("OxidationMaskBendingVelocityField: Cartesian solve region "
                    "exceeds maxGridPoints.")
          .print();
    }
  }

  void buildNodes() {
    nodes.clear();
    nodeLookup.clear();

    ConstSparseIterator maskIt(maskInterface->getDomain());
    IndexType index = minIndex;
    while (true) {
      if (isInsideMask(maskIt, index)) {
        const std::size_t id = nodes.size();
        nodeLookup.emplace(linearIndex(index), id);
        nodes.push_back({index});
      }

      if (!increment(index))
        break;
    }

    contactNodes = 0;
    viennahrle::IndexType bottomIndex = std::numeric_limits<viennahrle::IndexType>::max();
    for (const auto &node : nodes)
      bottomIndex = std::min(bottomIndex, node.index[D - 1]);

    for (auto &node : nodes) {
      node.contact = touchesContactBoundary(maskIt, node.index) ||
                     node.index[D - 1] <= bottomIndex + 1;
      if (node.contact)
        ++contactNodes;
    }
  }

  void solveElasticVelocity() {
    iterations = 0;
    residual = 0.;
    if (nodes.empty())
      return;

    std::vector<Vec3D<T>> previous(nodes.size(), {0., 0., 0.});
    std::vector<Vec3D<T>> next = previous;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      if (nodes[i].contact)
        previous[i] = contactVelocity(nodes[i].index);
    }
    next = previous;

    ConstSparseIterator maskIt(maskInterface->getDomain());
    const T lambda = lameLambda();
    const T mu = lameMu();
    const T gradDivWeight =
        (lambda + mu) / std::max(lambda + T(2) * mu,
                                 std::numeric_limits<T>::epsilon());

    for (; iterations < parameters.maxIterations; ++iterations) {
      residual = 0.;
      for (std::size_t nodeId = 0; nodeId < nodes.size(); ++nodeId) {
        if (nodes[nodeId].contact) {
          next[nodeId] = contactVelocity(nodes[nodeId].index);
        } else {
          Vec3D<T> laplaceAverage{0., 0., 0.};
          unsigned count = 0;
          for (unsigned direction = 0; direction < D; ++direction) {
            for (int offset : {-1, 1}) {
              addTo(laplaceAverage, neighborVelocity(maskIt, previous, nodeId,
                                                     direction, offset));
              ++count;
            }
          }

          Vec3D<T> updated =
              (count == 0) ? previous[nodeId]
                           : scaled(laplaceAverage, T(1) / count);
          const Vec3D<T> gradDivCorrection =
              scaled(divergenceGradient(maskIt, previous, nodes[nodeId].index),
                     gradDivWeight * gridDelta * gridDelta /
                         (T(2) * static_cast<T>(D)));
          updated = add(updated, gradDivCorrection);

          next[nodeId] =
              add(scaled(updated, parameters.relaxation),
                  scaled(previous[nodeId], T(1) - parameters.relaxation));
        }

        const Vec3D<T> delta = subtract(next[nodeId], previous[nodeId]);
        for (unsigned component = 0; component < D; ++component)
          residual = std::max(residual, std::abs(delta[component]));
      }

      previous.swap(next);
      if (residual < parameters.tolerance)
        break;
    }

    for (std::size_t i = 0; i < nodes.size(); ++i)
      nodes[i].velocity = previous[i];
  }

  Vec3D<T> neighborVelocity(ConstSparseIterator &maskIt,
                            const std::vector<Vec3D<T>> &velocity,
                            std::size_t nodeId, unsigned direction,
                            int offset) const {
    IndexType neighbor = nodes[nodeId].index;
    neighbor[direction] += offset;
    if (!inBounds(neighbor))
      return velocity[nodeId];

    const auto found = nodeLookup.find(linearIndex(neighbor));
    if (found != nodeLookup.end())
      return velocity[found->second];

    const T insidePhi = valueAt(maskIt, nodes[nodeId].index);
    const T outsidePhi = valueAt(maskIt, neighbor);
    const T distance = crossingDistance(insidePhi, outsidePhi);
    if (isContactBoundary(nodes[nodeId].index, direction, offset))
      return contactVelocity(nodes[nodeId].index);

    (void)distance;
    return velocity[nodeId];
  }

  Vec3D<T> divergenceGradient(ConstSparseIterator &maskIt,
                              const std::vector<Vec3D<T>> &velocity,
                              const IndexType &index) const {
    Vec3D<T> gradient{0., 0., 0.};
    for (unsigned component = 0; component < D; ++component) {
      IndexType plus = index;
      IndexType minus = index;
      plus[component] += 1;
      minus[component] -= 1;
      const T plusDiv = divergence(maskIt, velocity, plus);
      const T minusDiv = divergence(maskIt, velocity, minus);
      gradient[component] = (plusDiv - minusDiv) / (T(2) * gridDelta);
    }
    return gradient;
  }

  T divergence(ConstSparseIterator &maskIt,
               const std::vector<Vec3D<T>> &velocity,
               const IndexType &index) const {
    if (!inBounds(index))
      return 0.;
    const auto found = nodeLookup.find(linearIndex(index));
    if (found == nodeLookup.end())
      return 0.;

    T result = 0.;
    for (unsigned direction = 0; direction < D; ++direction) {
      const auto plus =
          neighborVelocity(maskIt, velocity, found->second, direction, 1);
      const auto minus =
          neighborVelocity(maskIt, velocity, found->second, direction, -1);
      result += (plus[direction] - minus[direction]) / (T(2) * gridDelta);
    }
    return result;
  }

  bool touchesContactBoundary(ConstSparseIterator &maskIt,
                              const IndexType &index) const {
    for (unsigned direction = 0; direction < D; ++direction) {
      for (int offset : {-1, 1}) {
        IndexType neighbor = index;
        neighbor[direction] += offset;
        if (!inBounds(neighbor))
          continue;
        if (nodeLookup.find(linearIndex(neighbor)) != nodeLookup.end())
          continue;
        if (!crosses(valueAt(maskIt, index), valueAt(maskIt, neighbor)))
          continue;
        if (isContactBoundary(index, direction, offset))
          return true;
      }
    }
    return false;
  }

  bool isContactBoundary(const IndexType &index, unsigned direction,
                         int offset) const {
    const auto normal = maskNormal(index);
    return direction == D - 1 && static_cast<T>(offset) * normal[direction] < T(0);
  }

  Vec3D<T> contactVelocity(const IndexType &index) const {
    Vec3D<T> coordinate{0., 0., 0.};
    for (unsigned i = 0; i < D; ++i)
      coordinate[i] = index[i] * gridDelta;
    Vec3D<T> oxideSampleCoordinate = coordinate;
    oxideSampleCoordinate[D - 1] -= gridDelta;

    Vec3D<T> velocity =
        scaled(deformationField->getVelocity(oxideSampleCoordinate),
               elasticCompliance() * parameters.velocityScale);
    velocity = add(velocity,
                   scaled(maskNormal(index),
                          elasticCompliance() * parameters.velocityScale *
                              deformationField
                                  ->getAverageBoundaryExpansionVelocity()));
    const T pressure =
        std::abs(deformationField->getPressure(oxideSampleCoordinate) -
                 parameters.referencePressure);
    if (parameters.pressureVelocityScale != T(0) && pressure > T(0)) {
      velocity = add(velocity,
                     scaled(maskNormal(index), elasticCompliance() *
                                                   parameters.pressureVelocityScale *
                                                   pressure));
    }
    return clampMagnitude(velocity);
  }

  Vec3D<T> maskNormal(const IndexType &index) const {
    ConstSparseIterator maskIt(maskInterface->getDomain());
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
      normal[i] = valueAt(maskIt, pos) - valueAt(maskIt, neg);
      norm += normal[i] * normal[i];
    }
    if (norm <= std::numeric_limits<T>::epsilon()) {
      normal[D - 1] = -1.;
      return normal;
    }

    norm = std::sqrt(norm);
    for (unsigned i = 0; i < D; ++i)
      normal[i] /= norm;
    return normal;
  }

  T elasticCompliance() const {
    const T nu = std::clamp(parameters.poissonRatio, T(-0.95), T(0.49));
    const T thickness = std::max(parameters.thickness,
                                 std::numeric_limits<T>::epsilon());
    const T referenceThickness =
        std::max(parameters.referenceThickness,
                 std::numeric_limits<T>::epsilon());
    const T bendingRigidity =
        parameters.youngModulus * thickness * thickness * thickness /
        (T(12) * (T(1) - nu * nu));
    const T referenceRigidity =
        parameters.youngModulus * referenceThickness * referenceThickness *
        referenceThickness / (T(12) * (T(1) - nu * nu));

    if (bendingRigidity <= std::numeric_limits<T>::epsilon())
      return T(1);
    return std::clamp(referenceRigidity / bendingRigidity, T(0), T(1));
  }

  T lameMu() const {
    const T nu = std::clamp(parameters.poissonRatio, T(-0.95), T(0.49));
    return parameters.youngModulus / (T(2) * (T(1) + nu));
  }

  T lameLambda() const {
    const T nu = std::clamp(parameters.poissonRatio, T(-0.95), T(0.49));
    return parameters.youngModulus * nu /
           ((T(1) + nu) * (T(1) - T(2) * nu));
  }

  Vec3D<T> getVelocity(const IndexType &index) const {
    const auto found = nodeLookup.find(linearIndex(index));
    if (found != nodeLookup.end())
      return nodes[found->second].velocity;

    T bestDistance2 = std::numeric_limits<T>::max();
    std::size_t bestNode = std::numeric_limits<std::size_t>::max();
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      T distance2 = 0.;
      for (unsigned d = 0; d < D; ++d) {
        const T delta = static_cast<T>(nodes[i].index[d] - index[d]);
        distance2 += delta * delta;
      }
      if (distance2 < bestDistance2) {
        bestDistance2 = distance2;
        bestNode = i;
      }
    }

    if (bestNode == std::numeric_limits<std::size_t>::max())
      return {0., 0., 0.};
    return nodes[bestNode].velocity;
  }

  bool isInsideMask(ConstSparseIterator &maskIt, const IndexType &index) const {
    return maskSign * valueAt(maskIt, index) >= 0.;
  }

  bool crosses(T a, T b) const {
    return (a <= 0. && b >= 0.) || (a >= 0. && b <= 0.);
  }

  T crossingDistance(T insidePhi, T outsidePhi) const {
    const T denom = std::abs(insidePhi) + std::abs(outsidePhi);
    if (denom <= std::numeric_limits<T>::epsilon())
      return gridDelta;
    return std::max(parameters.minBoundaryDistance * gridDelta,
                    gridDelta * std::abs(insidePhi) / denom);
  }

  T valueAt(ConstSparseIterator &it, const IndexType &index) const {
    it.goToIndices(index);
    return it.getValue();
  }

  bool inBounds(const IndexType &index) const {
    for (unsigned i = 0; i < D; ++i) {
      if (index[i] < minIndex[i] || index[i] > maxIndex[i])
        return false;
    }
    return true;
  }

  std::size_t linearIndex(const IndexType &index) const {
    if (!inBounds(index))
      return std::numeric_limits<std::size_t>::max();

    std::size_t result = 0;
    for (unsigned i = 0; i < D; ++i)
      result += static_cast<std::size_t>(index[i] - minIndex[i]) * strides[i];
    return result;
  }

  bool increment(IndexType &index) const {
    for (unsigned i = 0; i < D; ++i) {
      if (index[i] < maxIndex[i]) {
        ++index[i];
        return true;
      }
      index[i] = minIndex[i];
    }
    return false;
  }

  Vec3D<T> clampMagnitude(const Vec3D<T> &velocity) const {
    if (parameters.maxVelocity == std::numeric_limits<T>::max())
      return velocity;

    T magnitudeSquared = 0.;
    for (unsigned i = 0; i < 3; ++i)
      magnitudeSquared += velocity[i] * velocity[i];
    const T magnitude = std::sqrt(magnitudeSquared);
    if (magnitude <= parameters.maxVelocity || magnitude == T(0))
      return velocity;

    return scaled(velocity, parameters.maxVelocity / magnitude);
  }

  static Vec3D<T> scaled(const Vec3D<T> &source, T factor) {
    Vec3D<T> result{0., 0., 0.};
    for (unsigned i = 0; i < 3; ++i)
      result[i] = source[i] * factor;
    return result;
  }

  static Vec3D<T> add(const Vec3D<T> &a, const Vec3D<T> &b) {
    Vec3D<T> result{0., 0., 0.};
    for (unsigned i = 0; i < 3; ++i)
      result[i] = a[i] + b[i];
    return result;
  }

  static Vec3D<T> subtract(const Vec3D<T> &a, const Vec3D<T> &b) {
    Vec3D<T> result{0., 0., 0.};
    for (unsigned i = 0; i < 3; ++i)
      result[i] = a[i] - b[i];
    return result;
  }

  static void addTo(Vec3D<T> &target, const Vec3D<T> &source) {
    for (unsigned i = 0; i < 3; ++i)
      target[i] += source[i];
  }
};

/// Velocity field for the oxide outer interface when part of that interface is
/// in contact with a mask. Open oxide/ambient regions use the oxide deformation
/// field. Mask-contact regions are kinematically constrained to the solved mask
/// vector velocity and do not receive additional free-surface scalar growth.
template <class T, int D>
class OxidationConstrainedAmbientVelocityField final : public VelocityField<T> {
  using IndexType = viennahrle::Index<D>;
  using ConstSparseIterator =
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>;

  SmartPointer<OxidationDeformationVelocityField<T, D>> deformationField =
      nullptr;
  SmartPointer<OxidationMaskBendingVelocityField<T, D>> maskVelocityField =
      nullptr;
  SmartPointer<Domain<T, D>> maskInterface = nullptr;
  int maskSign = 1;
  T contactTolerance = 1.;

public:
  OxidationConstrainedAmbientVelocityField() = default;

  OxidationConstrainedAmbientVelocityField(
      SmartPointer<OxidationDeformationVelocityField<T, D>> passedDeformation,
      SmartPointer<OxidationMaskBendingVelocityField<T, D>> passedMaskVelocity,
      SmartPointer<Domain<T, D>> passedMaskInterface, int passedMaskSign = 1,
      T passedContactTolerance = 1.)
      : deformationField(passedDeformation),
        maskVelocityField(passedMaskVelocity), maskInterface(passedMaskInterface),
        maskSign((passedMaskSign < 0) ? -1 : 1),
        contactTolerance(passedContactTolerance) {}

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<OxidationConstrainedAmbientVelocityField>::New(
        std::forward<Args>(args)...);
  }

  void setContactTolerance(T passedContactTolerance) {
    contactTolerance = passedContactTolerance;
  }

  Vec3D<T> getVectorVelocity(const Vec3D<T> &coordinate, int material,
                             const Vec3D<T> &normalVector,
                             unsigned long pointId) final {
    if (isMaskContact(coordinate) && maskVelocityField != nullptr)
      return maskVelocityField->getVectorVelocity(coordinate, material,
                                                  normalVector, pointId);
    if (deformationField == nullptr)
      return {0., 0., 0.};
    return deformationField->getVectorVelocity(coordinate, material,
                                               normalVector, pointId);
  }

  T getScalarVelocity(const Vec3D<T> &coordinate, int material,
                      const Vec3D<T> &normalVector,
                      unsigned long pointId) final {
    if (isMaskContact(coordinate))
      return 0.;
    if (deformationField == nullptr)
      return 0.;
    return deformationField->getScalarVelocity(coordinate, material,
                                               normalVector, pointId);
  }

  T getDissipationAlpha(int direction, int material,
                        const Vec3D<T> &centralDifferences) final {
    T alpha = 0.;
    if (deformationField != nullptr)
      alpha = std::max(alpha, deformationField->getDissipationAlpha(
                                  direction, material, centralDifferences));
    if (maskVelocityField != nullptr)
      alpha = std::max(alpha, maskVelocityField->getDissipationAlpha(
                                  direction, material, centralDifferences));
    return alpha;
  }

private:
  bool isMaskContact(const Vec3D<T> &coordinate) const {
    if (maskInterface == nullptr)
      return false;

    ConstSparseIterator maskIt(maskInterface->getDomain());
    IndexType index;
    const T gridDelta = maskInterface->getGrid().getGridDelta();
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);

    maskIt.goToIndices(index);
    return maskSign * maskIt.getValue() >= -contactTolerance * gridDelta;
  }
};

/// Iterates diffusion, oxide deformation, and pressure-dependent reaction-rate
/// feedback on the shared Cartesian solve grid.
template <class T, int D> class OxidationCoupledModel {
  using IndexType = viennahrle::Index<D>;

  SmartPointer<OxidationDiffusionVelocityField<T, D>> diffusionField = nullptr;
  SmartPointer<OxidationDeformationVelocityField<T, D>> deformationField =
      nullptr;
  OxidationCouplingParameters<T> parameters;
  IndexType minIndex{};
  IndexType maxIndex{};
  bool useRequestedBounds = false;
  unsigned iterations = 0;
  T residual = std::numeric_limits<T>::max();

public:
  OxidationCoupledModel() = default;

  OxidationCoupledModel(
      SmartPointer<OxidationDiffusionVelocityField<T, D>> passedDiffusionField,
      SmartPointer<OxidationDeformationVelocityField<T, D>>
          passedDeformationField,
      OxidationCouplingParameters<T> passedParameters = {})
      : diffusionField(passedDiffusionField),
        deformationField(passedDeformationField), parameters(passedParameters) {}

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<OxidationCoupledModel>::New(std::forward<Args>(args)...);
  }

  void setDiffusionField(
      SmartPointer<OxidationDiffusionVelocityField<T, D>> passedDiffusionField) {
    diffusionField = passedDiffusionField;
  }

  void setDeformationField(
      SmartPointer<OxidationDeformationVelocityField<T, D>>
          passedDeformationField) {
    deformationField = passedDeformationField;
  }

  void setParameters(OxidationCouplingParameters<T> passedParameters) {
    parameters = passedParameters;
  }

  void setSolveBounds(const IndexType &passedMinIndex,
                      const IndexType &passedMaxIndex) {
    minIndex = passedMinIndex;
    maxIndex = passedMaxIndex;
    useRequestedBounds = true;
  }

  void clearSolveBounds() { useRequestedBounds = false; }

  void apply() {
    iterations = 0;
    residual = 0.;
    if (diffusionField == nullptr || deformationField == nullptr) {
      Logger::getInstance()
          .addError("OxidationCoupledModel: Missing diffusion or deformation "
                    "field.")
          .print();
      return;
    }

    if (!useRequestedBounds) {
      Logger::getInstance()
          .addError("OxidationCoupledModel: Solve bounds must be set.")
          .print();
      return;
    }

    diffusionField->setSolveBounds(minIndex, maxIndex);
    deformationField->setSolveBounds(minIndex, maxIndex);

    std::unordered_map<std::size_t, T> previousPressures;
    for (; iterations < parameters.maxIterations; ++iterations) {
      diffusionField->apply();
      deformationField->apply();

      residual = updateDiffusionPressureField(previousPressures);
      if (residual < parameters.tolerance)
        break;
    }

    diffusionField->apply();
    deformationField->apply();
  }

  unsigned getIterations() const { return iterations; }
  T getResidual() const { return residual; }

private:
  T updateDiffusionPressureField(
      std::unordered_map<std::size_t, T> &previousPressures) {
    T maxChange = 0.;
    T maxPressure = 0.;
    deformationField->forEachSolutionNode(
        [&](const IndexType &index, T newPressure) {
      const auto key = pressureKey(index);
      const auto previousIt = previousPressures.find(key);
      const T oldPressure =
          (previousIt == previousPressures.end()) ? T(0) : previousIt->second;
      const T relaxedPressure =
          parameters.relaxation * newPressure +
          (T(1) - parameters.relaxation) * oldPressure;

      previousPressures[key] = relaxedPressure;
      diffusionField->setPressure(index, relaxedPressure);
      maxChange = std::max(maxChange, std::abs(relaxedPressure - oldPressure));
      maxPressure = std::max(maxPressure, std::abs(relaxedPressure));
        });
    if (maxPressure <= std::numeric_limits<T>::epsilon())
      return maxChange;
    return maxChange / maxPressure;
  }

  std::size_t pressureKey(const IndexType &index) const {
    std::size_t seed = 0;
    for (unsigned i = 0; i < D; ++i) {
      const auto value = std::hash<long long>{}(static_cast<long long>(index[i]));
      seed ^= value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

} // namespace viennals
