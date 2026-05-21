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

namespace detail {

template <int D>
inline std::size_t gridIndexHash(const viennahrle::Index<D> &index) {
  std::size_t seed = 0;
  for (unsigned i = 0; i < static_cast<unsigned>(D); ++i) {
    seed ^= std::hash<long long>{}(static_cast<long long>(index[i])) +
             0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}

template <class T>
inline Vec3D<T> vecScaled(const Vec3D<T> &source, T factor) {
  Vec3D<T> result{0., 0., 0.};
  for (unsigned i = 0; i < 3; ++i)
    result[i] = source[i] * factor;
  return result;
}

template <class T>
inline Vec3D<T> vecAdd(const Vec3D<T> &a, const Vec3D<T> &b) {
  Vec3D<T> result{0., 0., 0.};
  for (unsigned i = 0; i < 3; ++i)
    result[i] = a[i] + b[i];
  return result;
}

template <class T>
inline Vec3D<T> vecSubtract(const Vec3D<T> &a, const Vec3D<T> &b) {
  Vec3D<T> result{0., 0., 0.};
  for (unsigned i = 0; i < 3; ++i)
    result[i] = a[i] - b[i];
  return result;
}

template <class T>
inline void vecAddTo(Vec3D<T> &target, const Vec3D<T> &source) {
  for (unsigned i = 0; i < 3; ++i)
    target[i] += source[i];
}

} // namespace detail

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
    return detail::gridIndexHash<D>(index);
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

      if (distance2 > 0 && inBounds(candidate)) {
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

} // namespace viennals
