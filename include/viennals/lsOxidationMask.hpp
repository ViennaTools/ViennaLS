#pragma once

#include <lsOxidationSolverBase.hpp>
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
class OxidationMaskBendingVelocityField final
    : public VelocityField<T>,
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
  using OxidationSolverBase<T, D>::initializeGridFromMask;

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

  bool solved = false;
  bool useRequestedBounds = false;
  T residual = std::numeric_limits<T>::max();
  unsigned iterations = 0;
  std::array<T, D> maxVelocity_{};
  std::size_t contactNodes = 0;

public:
  std::vector<Node> nodes;
  IndexType requestedMinIndex{};
  IndexType requestedMaxIndex{};

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

    maxVelocity_.fill(T(0));
    for (const auto &node : nodes) {
      for (unsigned d = 0; d < D; ++d)
        maxVelocity_[d] = std::max(maxVelocity_[d], std::abs(node.velocity[d]));
    }

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
          detail::vecScaled(deformationField->getVelocity(coordinate),
                            elasticCompliance() * parameters.velocityScale);
      const T pressure =
          std::abs(deformationField->getPressure(coordinate) -
                   parameters.referencePressure);
      if (parameters.pressureVelocityScale != T(0) && pressure > T(0))
        velocity = detail::vecAdd(velocity,
                                  detail::vecScaled(normalVector, elasticCompliance() *
                                                                    parameters.pressureVelocityScale *
                                                                    pressure));
      return clampMagnitude(velocity);
    }

    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return clampMagnitude(getVelocity(index));
  }

  T getDissipationAlpha(int direction, int material,
                        const Vec3D<T> &centralDifferences) final {
    T alpha = maxVelocity_[direction];
    if (nodes.empty() && deformationField != nullptr) {
      alpha = std::max(alpha, deformationField->getDissipationAlpha(
                                  direction, material, centralDifferences) *
                                  elasticCompliance() *
                                  parameters.velocityScale);
    }
    if (parameters.maxVelocity != std::numeric_limits<T>::max()) {
      alpha = std::min(alpha, parameters.maxVelocity);
    }
    return alpha;
  }

private:
  void initialiseGrid() {
    initializeGridFromMask(maskInterface, useRequestedBounds, requestedMinIndex,
                           requestedMaxIndex, parameters.maxGridPoints,
                           "OxidationMaskBendingVelocityField");
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
              detail::vecAddTo(laplaceAverage, neighborVelocity(maskIt, previous, nodeId,
                                                                direction, offset));
              ++count;
            }
          }

          Vec3D<T> updated =
              (count == 0) ? previous[nodeId]
                           : detail::vecScaled(laplaceAverage, T(1) / count);
          const Vec3D<T> gradDivCorrection =
              detail::vecScaled(divergenceGradient(maskIt, previous, nodes[nodeId].index),
                                gradDivWeight * gridDelta * gridDelta /
                                    (T(2) * static_cast<T>(D)));
          updated = detail::vecAdd(updated, gradDivCorrection);

          next[nodeId] =
              detail::vecAdd(detail::vecScaled(updated, parameters.relaxation),
                             detail::vecScaled(previous[nodeId], T(1) - parameters.relaxation));
        }

        const Vec3D<T> delta = detail::vecSubtract(next[nodeId], previous[nodeId]);
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
        detail::vecScaled(deformationField->getVelocity(oxideSampleCoordinate),
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

  bool isInsideMask(ConstSparseIterator &maskIt, const IndexType &index) const {
    return maskSign * valueAt(maskIt, index) >= 0.;
  }

  T crossingDistance(T insidePhi, T outsidePhi) const {
    return detail::levelSetCrossingDistance(
        insidePhi, outsidePhi, parameters.minBoundaryDistance, gridDelta);
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

    return detail::vecScaled(velocity, parameters.maxVelocity / magnitude);
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
