#pragma once

#include <lsOxidationDeformation.hpp>

namespace viennals {

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
  T minBoundaryDistance = 1e-6;
  unsigned maxIterations = 10000;
  std::size_t maxGridPoints = 5000000;
  int material = -1;
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
    Vec3D<T> outwardNormal{0., 0., -1.}; // unit normal from mask toward oxide
    bool contact = false;
  };

  SmartPointer<OxidationDeformationVelocityField<T, D>> deformationField =
      nullptr;
  SmartPointer<Domain<T, D>> maskInterface = nullptr;
  SmartPointer<Domain<T, D>> ambientInterface = nullptr;
  OxidationMaskParameters<T> parameters;
  int maskSign = 1;
  int ambientSign = -1;

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

  /// Provide the SiO₂/ambient interface so that contact faces can be detected
  /// on any mask face that borders the oxide, not only the bottom face.
  /// ambientSign follows the same convention as OxidationDeformationVelocityField:
  /// ambientSign * φ_ambient >= 0 selects nodes inside the oxide band.
  void setAmbientInterface(SmartPointer<Domain<T, D>> passedAmbientInterface,
                           int passedAmbientSign = -1) {
    ambientInterface = passedAmbientInterface;
    ambientSign = (passedAmbientSign < 0) ? -1 : 1;
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
    for (auto &node : nodes) {
      node.contact = touchesContactBoundary(maskIt, node.index);
      if (node.contact) {
        ++contactNodes;
        // maskNormal is the gradient of phi_mask, pointing from the mask
        // interior (phi < 0) toward the oxide exterior (phi > 0). This is
        // the direction from the mask material toward the oxide, which is
        // where we sample the deformation velocity in contactVelocity.
        const auto normal = maskNormal(node.index, maskIt);
        for (unsigned i = 0; i < D; ++i)
          node.outwardNormal[i] = normal[i];
      }
    }
  }

  void solveElasticVelocity() {
    iterations = 0;
    residual = 0.;
    if (nodes.empty())
      return;

    ConstSparseIterator maskIt(maskInterface->getDomain());
    std::vector<Vec3D<T>> previous(nodes.size(), {0., 0., 0.});
    std::vector<Vec3D<T>> next = previous;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      if (nodes[i].contact)
        previous[i] = contactVelocity(nodes[i]);
    }
    next = previous;

    const T lambda = lameLambda();
    const T mu = lameMu();
    const T gradDivWeight =
        (lambda + mu) / std::max(lambda + T(2) * mu,
                                 std::numeric_limits<T>::epsilon());

    for (; iterations < parameters.maxIterations; ++iterations) {
      residual = 0.;
      for (std::size_t nodeId = 0; nodeId < nodes.size(); ++nodeId) {
        if (nodes[nodeId].contact) {
          next[nodeId] = contactVelocity(nodes[nodeId]);
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

    if (isContactBoundary(nodes[nodeId].index, direction, offset, maskIt))
      return contactVelocity(nodes[nodeId]);

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
        if (!inBounds(neighbor)) {
          // Solve region is clipped to the mask/oxide interface: an
          // out-of-bounds neighbor is by definition outside the mask.
          if (isContactBoundary(index, direction, offset, maskIt))
            return true;
          continue;
        }
        if (nodeLookup.find(linearIndex(neighbor)) != nodeLookup.end())
          continue;
        if (!crosses(valueAt(maskIt, index), valueAt(maskIt, neighbor)))
          continue;
        if (isContactBoundary(index, direction, offset, maskIt))
          return true;
      }
    }
    return false;
  }

  bool isContactBoundary(const IndexType &index, unsigned direction,
                         int offset, ConstSparseIterator &maskIt) const {
    const auto normal = maskNormal(index, maskIt);
    if (static_cast<T>(offset) * normal[direction] >= T(0))
      return false; // face does not exit the mask bending domain
    if (ambientInterface == nullptr)
      return direction == D - 1; // no ambient LS: fall back to bottom-face only
    // Check that the current node (oxide-side, inside the bending domain) is
    // below the ambient surface, not in the gas above the mask.
    return isInsideOxide(index);
  }

  bool isInsideOxide(const IndexType &index) const {
    ConstSparseIterator ambientIt(ambientInterface->getDomain());
    ambientIt.goToIndices(index);
    return ambientSign * ambientIt.getValue() >= T(0);
  }

  Vec3D<T> contactVelocity(const Node &node) const {
    Vec3D<T> oxideSampleCoordinate{0., 0., 0.};
    for (unsigned i = 0; i < D; ++i)
      oxideSampleCoordinate[i] =
          node.index[i] * gridDelta + node.outwardNormal[i] * gridDelta;

    return clampMagnitude(
        scaled(deformationField->getVelocity(oxideSampleCoordinate),
               elasticCompliance() * parameters.velocityScale));
  }

  Vec3D<T> maskNormal(const IndexType &index, ConstSparseIterator &maskIt) const {
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

  T clampedPoissonRatio() const {
    return std::clamp(parameters.poissonRatio, T(-0.95), T(0.49));
  }

  T elasticCompliance() const {
    const T nu = clampedPoissonRatio();
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
    const T nu = clampedPoissonRatio();
    return parameters.youngModulus / (T(2) * (T(1) + nu));
  }

  T lameLambda() const {
    const T nu = clampedPoissonRatio();
    return parameters.youngModulus * nu /
           ((T(1) + nu) * (T(1) - T(2) * nu));
  }

  Vec3D<T> getVelocity(const IndexType &index) const {
    const auto found = nodeLookup.find(linearIndex(index));
    if (found != nodeLookup.end())
      return nodes[found->second].velocity;

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

  static Vec3D<T> scaled(const Vec3D<T> &v, T f) { return detail::vecScaled(v, f); }
  static Vec3D<T> add(const Vec3D<T> &a, const Vec3D<T> &b) { return detail::vecAdd(a, b); }
  static Vec3D<T> subtract(const Vec3D<T> &a, const Vec3D<T> &b) { return detail::vecSubtract(a, b); }
  static void addTo(Vec3D<T> &t, const Vec3D<T> &s) { detail::vecAddTo(t, s); }
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

public:
  OxidationConstrainedAmbientVelocityField() = default;

  OxidationConstrainedAmbientVelocityField(
      SmartPointer<OxidationDeformationVelocityField<T, D>> passedDeformation,
      SmartPointer<OxidationMaskBendingVelocityField<T, D>> passedMaskVelocity,
      SmartPointer<Domain<T, D>> passedMaskInterface, int passedMaskSign = 1)
      : deformationField(passedDeformation),
        maskVelocityField(passedMaskVelocity), maskInterface(passedMaskInterface),
        maskSign((passedMaskSign < 0) ? -1 : 1) {}

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<OxidationConstrainedAmbientVelocityField>::New(
        std::forward<Args>(args)...);
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
    return maskSign * maskIt.getValue() >= T(0);
  }
};

} // namespace viennals
