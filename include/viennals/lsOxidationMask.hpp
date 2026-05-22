#pragma once

#include <lsOxidationSolverBase.hpp>
#include <lsOxidationDeformation.hpp>

namespace viennals {

enum class OxidationMaskAnchorMode : int {
  NONE = 0,
  MIN_BOUNDARY = 1,
  MAX_BOUNDARY = 2,
  BOTH_BOUNDARIES = 3
};

template <class T> struct OxidationMaskParameters {
  // Effective viscosity of the mask material (Pa·hr). Controls how fast the
  // mask creeps under the oxide traction: larger value → stiffer → less
  // deflection per unit oxide pressure. This is an effective stack viscosity;
  // fitted process values may be higher than intrinsic high-temperature Si3N4.
  T maskViscosity = 5e8;
  T poissonRatio = 0.27;
  T maxVelocity = std::numeric_limits<T>::max();
  bool unilateralContact = true;
  // Maximum oxide/mask separation, in level-set grid units, still considered
  // closed contact. 0 requires the oxide-side ghost node to be inside oxide.
  T contactGapTolerance = 0.;
  OxidationMaskAnchorMode anchorMode =
      OxidationMaskAnchorMode::MIN_BOUNDARY;
  unsigned anchorDirection = 0;
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
    bool contact = false;
    bool anchor = false;
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
  T lastApplyVelocityChange = std::numeric_limits<T>::max();

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
  T getLastApplyVelocityChange() const { return lastApplyVelocityChange; }

  void apply() {
    if (deformationField == nullptr || maskInterface == nullptr) {
      solved = true;
      return;
    }

    const auto previousVelocities = collectVelocitiesByGridPoint();
    initialiseGrid();
    buildNodes();
    solveElasticVelocity();
    lastApplyVelocityChange = maxVelocityChange(previousVelocities);

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

    if (maskInterface == nullptr || nodes.empty())
      return {0., 0., 0.};

    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return clampMagnitude(getVelocity(index));
  }

  T getDissipationAlpha(int direction, int /*material*/,
                        const Vec3D<T> & /*centralDifferences*/) final {
    T alpha = maxVelocity_[direction];
    if (parameters.maxVelocity != std::numeric_limits<T>::max())
      alpha = std::min(alpha, parameters.maxVelocity);
    return alpha;
  }

private:
  void initialiseGrid() {
    initializeGridFromMask(maskInterface, useRequestedBounds, requestedMinIndex,
                           requestedMaxIndex, parameters.maxGridPoints,
                           "OxidationMaskBendingVelocityField");
  }

  std::unordered_map<std::size_t, Vec3D<T>> collectVelocitiesByGridPoint() const {
    std::unordered_map<std::size_t, Vec3D<T>> result;
    result.reserve(nodes.size());
    for (const auto &node : nodes)
      if (inBounds(node.index))
        result.emplace(linearIndex(node.index), node.velocity);
    return result;
  }

  T maxVelocityChange(
      const std::unordered_map<std::size_t, Vec3D<T>> &previous) const {
    if (previous.empty())
      return std::numeric_limits<T>::max();

    T maxChange = 0.;
    T maxMagnitude = 0.;
    for (const auto &node : nodes) {
      if (!node.contact)
        continue;
      const auto found = previous.find(linearIndex(node.index));
      if (found == previous.end())
        continue;

      T changeSquared = 0.;
      T magnitudeSquared = 0.;
      for (unsigned i = 0; i < D; ++i) {
        const T delta = node.velocity[i] - found->second[i];
        changeSquared += delta * delta;
        magnitudeSquared += node.velocity[i] * node.velocity[i];
      }
      maxChange = std::max(maxChange, std::sqrt(changeSquared));
      maxMagnitude = std::max(maxMagnitude, std::sqrt(magnitudeSquared));
    }

    if (maxMagnitude <= std::numeric_limits<T>::epsilon())
      return maxChange;
    return maxChange / maxMagnitude;
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
      node.anchor = isRigidModeAnchor(node.index);
      if (node.contact)
        ++contactNodes;
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

    const T lambda = lameLambda();
    const T mu = lameMu();
    const T gradDivWeight =
        (lambda + mu) / std::max(lambda + T(2) * mu,
                                 std::numeric_limits<T>::epsilon());

    for (; iterations < parameters.maxIterations; ++iterations) {
      residual = 0.;
      for (std::size_t nodeId = 0; nodeId < nodes.size(); ++nodeId) {
        if (nodes[nodeId].anchor) {
          next[nodeId] = {0., 0., 0.};
          continue;
        }

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
      return tractionGhostVelocity(velocity, nodeId, direction, offset);

    return velocity[nodeId];
  }

  bool isRigidModeAnchor(const IndexType &index) const {
    if constexpr (D <= 1) {
      return false;
    } else {
      if (parameters.anchorMode == OxidationMaskAnchorMode::NONE)
        return false;

      const unsigned direction =
          std::min(parameters.anchorDirection, static_cast<unsigned>(D - 1));
      const bool minAnchor =
          parameters.anchorMode == OxidationMaskAnchorMode::MIN_BOUNDARY ||
          parameters.anchorMode == OxidationMaskAnchorMode::BOTH_BOUNDARIES;
      const bool maxAnchor =
          parameters.anchorMode == OxidationMaskAnchorMode::MAX_BOUNDARY ||
          parameters.anchorMode == OxidationMaskAnchorMode::BOTH_BOUNDARIES;
      return (minAnchor && index[direction] == minIndex[direction]) ||
             (maxAnchor && index[direction] == maxIndex[direction]);
    }
  }

  // Ghost velocity for a Neumann traction BC at the oxide/mask contact face.
  // The face outward normal (mask → oxide) is offset*ê_direction.
  // Oxide traction: t = σ_oxide · n̂  (full Cauchy stress, row-major 3×3)
  // Ghost: v_ghost = v_node + h/(λ+2μ)·t_n·n̂ + h/μ·t_tangential
  Vec3D<T> tractionGhostVelocity(const std::vector<Vec3D<T>> &velocity,
                                  std::size_t nodeId, unsigned direction,
                                  int offset) const {
    // Outward face normal: from mask into oxide
    Vec3D<T> faceNormal{0., 0., 0.};
    faceNormal[direction] = static_cast<T>(offset);

    // Sample oxide stress one grid cell into the oxide from the contact face
    Vec3D<T> oxidePt{0., 0., 0.};
    for (unsigned i = 0; i < D; ++i)
      oxidePt[i] = nodes[nodeId].index[i] * gridDelta;
    oxidePt[direction] += offset * gridDelta;

    // Full Cauchy stress σ = s_dev − p·I, stored row-major: σ[3i+j] = σ_ij
    const auto sigma = deformationField->getStressTensor(oxidePt);

    // Traction on the mask face: t_i = Σ_j σ_ij · n̂_j
    Vec3D<T> t{0., 0., 0.};
    for (unsigned i = 0; i < D; ++i)
      for (unsigned j = 0; j < D; ++j)
        t[i] += sigma[3 * i + j] * faceNormal[j];

    // Decompose into normal and tangential components
    T tn = 0.;
    for (unsigned i = 0; i < D; ++i)
      tn += t[i] * faceNormal[i];

    if (parameters.unilateralContact && tn >= T(0))
      return velocity[nodeId];

    const T mu = lameMu();
    const T lam = lameLambda();
    const T normalScale = gridDelta / std::max(lam + T(2) * mu,
                                               std::numeric_limits<T>::epsilon());
    const T tangentialScale = gridDelta / std::max(mu,
                                                   std::numeric_limits<T>::epsilon());

    // v_ghost = v_node + h/(λ+2μ)·t_n·n̂ + h/μ·t_tangential
    const Vec3D<T> &vNode = velocity[nodeId];
    Vec3D<T> vGhost{0., 0., 0.};
    for (unsigned i = 0; i < D; ++i) {
      const T tTangential = t[i] - tn * faceNormal[i];
      vGhost[i] = vNode[i] + normalScale * tn * faceNormal[i]
                            + tangentialScale * tTangential;
    }
    return vGhost;
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
    const T grad = maskGradientComponent(index, direction, maskIt);

    // The face exits the selected bending domain only if it points from an
    // inside node toward the outside of that domain.  Use the signed level-set
    // gradient instead of a normalized normal so HRLE far-field sentinels cannot
    // produce a spurious zero normal.
    if (static_cast<T>(maskSign) * static_cast<T>(offset) * grad >= T(0))
      return false;

    if (ambientInterface == nullptr)
      return direction == D - 1; // no ambient LS: fall back to bottom-face only

    IndexType ghostIndex = index;
    ghostIndex[direction] += offset;
    return isInsideOxide(ghostIndex);
  }

  bool isInsideOxide(const IndexType &index) const {
    ConstSparseIterator ambientIt(ambientInterface->getDomain());
    ambientIt.goToIndices(index);
    return ambientSign * ambientIt.getValue() >=
           -std::max(parameters.contactGapTolerance, T(0));
  }

  T maskGradientComponent(const IndexType &index, unsigned direction,
                          ConstSparseIterator &maskIt) const {
    IndexType pos = index;
    IndexType neg = index;
    pos[direction] += 1;
    neg[direction] -= 1;
    if (!inBounds(pos))
      pos = index;
    if (!inBounds(neg))
      neg = index;
    return detail::clampLevelSetPhi(valueAt(maskIt, pos)) -
           detail::clampLevelSetPhi(valueAt(maskIt, neg));
  }

  T clampedPoissonRatio() const {
    return std::clamp(parameters.poissonRatio, T(-0.95), T(0.49));
  }

  // Lamé viscosity parameters: same structure as elastic Lamé constants but
  // with maskViscosity (Pa·hr) in place of Young's modulus (Pa), making the
  // governing equation a viscous Stokes flow rather than elastic equilibrium.
  T lameMu() const {
    const T nu = clampedPoissonRatio();
    return parameters.maskViscosity / (T(2) * (T(1) + nu));
  }

  T lameLambda() const {
    const T nu = clampedPoissonRatio();
    return parameters.maskViscosity * nu /
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
