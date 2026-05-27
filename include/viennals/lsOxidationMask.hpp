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
  // Arrhenius creep viscosity law:
  // eta(T) = eta_ref * exp(E/R * (1/T - 1/T_ref)).
  // Viscosity is in Pa·hr, activation energy in J/mol, temperature in K.
  T temperature = 1273.15;
  T referenceTemperature = 1273.15;
  T referenceViscosity = 5e8;
  T creepActivationEnergy = 0.;
  T poissonRatio = 0.27;
  bool unilateralContact = true;
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

/// Vector velocity field for a compliant oxidation mask driven by solved oxide
/// traction. A Cartesian viscous elasticity solve is built inside the mask
/// level set, oxide-contact faces use a traction ghost velocity, and the
/// resulting vector field moves the entire mask body. Repeated applications use
/// Aitken relaxation on the contact-interface velocity update.
template <class T, int D>
class OxidationMaskBending final
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
    // contactFaceActive[direction*2 + (offset==+1?1:0)]:
    //   true  → traction BC applies; ghost = velocity[node] + contactFaceTraction[face]
    //   false → face is free or tensile-contact (use velocity[node] as ghost)
    // Populated once in buildNodes(); oxide stress doesn't change during inner solve.
    std::array<bool, 2 * D> contactFaceActive{};
    std::array<Vec3D<T>, 2 * D> contactFaceTraction{};
  };

  SmartPointer<OxidationDeformation<T, D>> deformationField =
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
  T aitkenOmega = 1.;
  std::vector<T> previousAitkenResidual;
  IndexType requestedMinIndex{};
  IndexType requestedMaxIndex{};
  std::vector<Node> nodes;

public:
  OxidationMaskBending() = default;

  OxidationMaskBending(
      SmartPointer<OxidationDeformation<T, D>> passedDeformation,
      OxidationMaskParameters<T> passedParameters = {})
      : deformationField(passedDeformation), parameters(passedParameters) {}

  OxidationMaskBending(
      SmartPointer<OxidationDeformation<T, D>> passedDeformation,
      SmartPointer<Domain<T, D>> passedMaskInterface,
      OxidationMaskParameters<T> passedParameters = {},
      int passedMaskSign = 1)
      : deformationField(passedDeformation), maskInterface(passedMaskInterface),
        parameters(passedParameters),
        maskSign((passedMaskSign < 0) ? -1 : 1) {}

  static SmartPointer<OxidationMaskBending>
  New(SmartPointer<OxidationDeformation<T, D>> passedDeformation,
      OxidationMaskParameters<T> passedParameters = {}) {
    return SmartPointer<OxidationMaskBending>::New(
        passedDeformation, passedParameters);
  }

  static SmartPointer<OxidationMaskBending>
  New(SmartPointer<OxidationDeformation<T, D>> passedDeformation,
      SmartPointer<Domain<T, D>> passedMaskInterface,
      OxidationMaskParameters<T> passedParameters = {},
      int passedMaskSign = 1) {
    return SmartPointer<OxidationMaskBending>::New(
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
  /// ambientSign follows the same convention as OxidationDeformation:
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
    correctAnchorFlags();
    solveElasticVelocity();

    const auto fixedPointResidual =
        contactResidualVector(previousVelocities);
    const T omega = aitkenRelaxation(fixedPointResidual);
    relaxVelocities(previousVelocities, omega);
    lastApplyVelocityChange = maxVelocityChange(previousVelocities);

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
    return getVelocity(index);
  }

  T getDissipationAlpha(int direction, int /*material*/,
                        const Vec3D<T> & /*centralDifferences*/) final {
    return maxVelocity_[direction];
  }

private:
  void initialiseGrid() {
    initializeGridFromMask(maskInterface, useRequestedBounds, requestedMinIndex,
                           requestedMaxIndex, parameters.maxGridPoints,
                           "OxidationMaskBending");
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

    T changeSquaredSum = 0.;
    T magnitudeSquaredSum = std::numeric_limits<T>::epsilon();
    for (const auto &node : nodes) {
      if (!node.contact)
        continue;
      const auto found = previous.find(linearIndex(node.index));
      if (found == previous.end())
        continue;

      for (unsigned i = 0; i < D; ++i) {
        const T delta = node.velocity[i] - found->second[i];
        changeSquaredSum += delta * delta;
        magnitudeSquaredSum += node.velocity[i] * node.velocity[i];
        magnitudeSquaredSum += found->second[i] * found->second[i];
      }
    }

    return std::sqrt(changeSquaredSum / magnitudeSquaredSum);
  }

  std::vector<T> contactResidualVector(
      const std::unordered_map<std::size_t, Vec3D<T>> &previous) const {
    std::vector<T> residualVector;
    residualVector.reserve(contactNodes * D);
    if (previous.empty())
      return residualVector;

    for (const auto &node : nodes) {
      if (!node.contact)
        continue;
      const auto found = previous.find(linearIndex(node.index));
      if (found == previous.end())
        continue;
      for (unsigned i = 0; i < D; ++i)
        residualVector.push_back(node.velocity[i] - found->second[i]);
    }
    return residualVector;
  }

  T aitkenRelaxation(const std::vector<T> &residualVector) {
    if (residualVector.empty()) {
      previousAitkenResidual.clear();
      aitkenOmega = 1.;
      return 1.;
    }

    T omega = aitkenOmega;
    if (previousAitkenResidual.size() == residualVector.size()) {
      T numerator = 0.;
      T denominator = 0.;
      for (std::size_t i = 0; i < residualVector.size(); ++i) {
        const T delta = residualVector[i] - previousAitkenResidual[i];
        numerator += previousAitkenResidual[i] * delta;
        denominator += delta * delta;
      }

      if (denominator > std::numeric_limits<T>::epsilon()) {
        omega = -aitkenOmega * numerator / denominator;
        omega = std::clamp(omega, T(0.05), T(1.5));
      }
    }

    previousAitkenResidual = residualVector;
    aitkenOmega = omega;
    return omega;
  }

  void relaxVelocities(const std::unordered_map<std::size_t, Vec3D<T>> &previous,
                       T omega) {
    if (previous.empty() || omega == T(1))
      return;

    for (auto &node : nodes) {
      const auto found = previous.find(linearIndex(node.index));
      if (found == previous.end())
        continue;
      for (unsigned i = 0; i < D; ++i)
        node.velocity[i] =
            found->second[i] + omega * (node.velocity[i] - found->second[i]);
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
      node.anchor = isRigidModeAnchor(node.index);
      if (node.contact)
        ++contactNodes;

      // Precompute traction BCs so the inner GS loop has no HRLE lookups or
      // stress-tensor queries; oxide stress is fixed during the inner solve.
      const T mu = lameMu();
      const T lam = lameLambda();
      const T normalScale =
          gridDelta / std::max(lam + T(2) * mu, std::numeric_limits<T>::epsilon());
      const T tangentialScale =
          gridDelta / std::max(mu, std::numeric_limits<T>::epsilon());

      for (unsigned dir = 0; dir < D; ++dir) {
        for (int offset : {-1, 1}) {
          const unsigned faceIdx = dir * 2u + (offset == 1 ? 1u : 0u);
          IndexType neighbor = node.index;
          neighbor[dir] += offset;
          const bool neighborIsNode =
              inBounds(neighbor) && nodeLookup.count(linearIndex(neighbor));
          if (neighborIsNode ||
              !isContactBoundary(node.index, dir, offset, maskIt)) {
            node.contactFaceActive[faceIdx] = false;
            node.contactFaceTraction[faceIdx] = {};
            continue;
          }

          // Sample oxide stress at the ghost node (one cell into the oxide).
          Vec3D<T> faceNormal{0., 0., 0.};
          faceNormal[dir] = static_cast<T>(offset);
          Vec3D<T> oxidePt{0., 0., 0.};
          for (unsigned i = 0; i < D; ++i)
            oxidePt[i] = node.index[i] * gridDelta;
          oxidePt[dir] += offset * gridDelta;

          const auto sigma = deformationField->getStressTensor(oxidePt);
          Vec3D<T> t{0., 0., 0.};
          for (unsigned i = 0; i < D; ++i)
            for (unsigned j = 0; j < D; ++j)
              t[i] += sigma[3 * i + j] * faceNormal[j];

          T tn = 0.;
          for (unsigned i = 0; i < D; ++i)
            tn += t[i] * faceNormal[i];

          if (parameters.unilateralContact && tn >= T(0)) {
            node.contactFaceActive[faceIdx] = false;
            node.contactFaceTraction[faceIdx] = {};
            continue;
          }

          // Fixed traction addition: ghost = vNode + contactFaceTraction
          Vec3D<T> traction{};
          for (unsigned i = 0; i < D; ++i) {
            const T tTangential = t[i] - tn * faceNormal[i];
            traction[i] = normalScale * tn * faceNormal[i]
                        + tangentialScale * tTangential;
          }
          node.contactFaceActive[faceIdx] = true;
          node.contactFaceTraction[faceIdx] = traction;
        }
      }
    }
  }

  // Re-tag anchor nodes based on actual node bounds (not padded solve bounds).
  // buildNodes() calls isRigidModeAnchor() which compares against the padded
  // solve-bounds minIndex/maxIndex; those are 4+ cells beyond actual mask nodes
  // so no node ever matches → no anchor → rigid-body floating.
  // This method fixes that by scanning the built node list for the true extremes.
  void correctAnchorFlags() {
    if (parameters.anchorMode == OxidationMaskAnchorMode::NONE || nodes.empty())
      return;

    const unsigned dir =
        std::min(parameters.anchorDirection, static_cast<unsigned>(D - 1));
    const bool applyMin =
        parameters.anchorMode == OxidationMaskAnchorMode::MIN_BOUNDARY ||
        parameters.anchorMode == OxidationMaskAnchorMode::BOTH_BOUNDARIES;
    const bool applyMax =
        parameters.anchorMode == OxidationMaskAnchorMode::MAX_BOUNDARY ||
        parameters.anchorMode == OxidationMaskAnchorMode::BOTH_BOUNDARIES;

    auto nodeMin = nodes[0].index[dir];
    auto nodeMax = nodes[0].index[dir];
    for (const auto &node : nodes) {
      if (node.index[dir] < nodeMin) nodeMin = node.index[dir];
      if (node.index[dir] > nodeMax) nodeMax = node.index[dir];
    }

    for (auto &node : nodes)
      node.anchor =
          (applyMin && node.index[dir] == nodeMin) ||
          (applyMax && node.index[dir] == nodeMax);
  }

  void solveElasticVelocity() {
    iterations = 0;
    residual = 0.;
    if (nodes.empty())
      return;

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
            detail::vecAddTo(laplaceAverage, neighborVelocity(previous, nodeId,
                                                              direction, offset));
            ++count;
          }
        }

        Vec3D<T> updated =
            (count == 0) ? previous[nodeId]
                         : detail::vecScaled(laplaceAverage, T(1) / count);
        const Vec3D<T> gradDivCorrection =
            detail::vecScaled(divergenceGradient(previous, nodes[nodeId].index),
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

  // Returns neighbor velocity without ghost extrapolation (used inside ghost
  // computations to avoid recursion).
  Vec3D<T> simpleNeighborVelocity(const std::vector<Vec3D<T>> &velocity,
                                  std::size_t nodeId, unsigned direction,
                                  int offset) const {
    IndexType neighbor = nodes[nodeId].index;
    neighbor[direction] += offset;
    if (!inBounds(neighbor))
      return velocity[nodeId];
    const auto found = nodeLookup.find(linearIndex(neighbor));
    return (found != nodeLookup.end()) ? velocity[found->second]
                                       : velocity[nodeId];
  }

  // Ghost velocity for a traction-free face (σ·n = 0).
  // The formula drops gridDelta because the derivative scale and the ghost
  // extension scale cancel: ghost[n] = node[n] - coeff * Σ_t (v+t[t]-v-t[t])/2,
  //                          ghost[t] = node[t] - (v+t[n]-v-t[n])/2.
  // This is independent of which side (offset) the ghost sits on.
  Vec3D<T> tractionFreeGhost(const std::vector<Vec3D<T>> &velocity,
                              std::size_t nodeId, unsigned normalDir) const {
    const T coeff = lameLambda() /
                    std::max(lameLambda() + T(2) * lameMu(),
                             std::numeric_limits<T>::epsilon());
    Vec3D<T> ghost = velocity[nodeId];
    for (unsigned tanDir = 0; tanDir < D; ++tanDir) {
      if (tanDir == normalDir)
        continue;
      const Vec3D<T> vPlus  = simpleNeighborVelocity(velocity, nodeId, tanDir, +1);
      const Vec3D<T> vMinus = simpleNeighborVelocity(velocity, nodeId, tanDir, -1);
      const T dvTanTan  = (vPlus[tanDir]    - vMinus[tanDir])   * T(0.5);
      const T dvNormTan = (vPlus[normalDir] - vMinus[normalDir]) * T(0.5);
      ghost[normalDir] -= coeff * dvTanTan;
      ghost[tanDir]    -= dvNormTan;
    }
    return ghost;
  }

  Vec3D<T> neighborVelocity(const std::vector<Vec3D<T>> &velocity,
                            std::size_t nodeId, unsigned direction,
                            int offset) const {
    IndexType neighbor = nodes[nodeId].index;
    neighbor[direction] += offset;

    if (inBounds(neighbor)) {
      const auto found = nodeLookup.find(linearIndex(neighbor));
      if (found != nodeLookup.end())
        return velocity[found->second];
    }

    const unsigned faceIdx = direction * 2u + (offset == 1 ? 1u : 0u);
    if (nodes[nodeId].contactFaceActive[faceIdx])
      return detail::vecAdd(velocity[nodeId],
                            nodes[nodeId].contactFaceTraction[faceIdx]);

    return tractionFreeGhost(velocity, nodeId, direction);
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

  Vec3D<T> divergenceGradient(const std::vector<Vec3D<T>> &velocity,
                              const IndexType &index) const {
    Vec3D<T> gradient{0., 0., 0.};
    for (unsigned component = 0; component < D; ++component) {
      IndexType plus = index;
      IndexType minus = index;
      plus[component] += 1;
      minus[component] -= 1;
      const T plusDiv = divergence(velocity, plus);
      const T minusDiv = divergence(velocity, minus);
      gradient[component] = (plusDiv - minusDiv) / (T(2) * gridDelta);
    }
    return gradient;
  }

  T divergence(const std::vector<Vec3D<T>> &velocity,
               const IndexType &index) const {
    if (!inBounds(index))
      return 0.;
    const auto found = nodeLookup.find(linearIndex(index));
    if (found == nodeLookup.end())
      return 0.;

    T result = 0.;
    for (unsigned direction = 0; direction < D; ++direction) {
      const auto plus =
          neighborVelocity(velocity, found->second, direction, 1);
      const auto minus =
          neighborVelocity(velocity, found->second, direction, -1);
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
    return ambientSign * ambientIt.getValue() >= T(0);
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

  static constexpr T gasConstant = T(8.31446261815324);

  T effectiveMaskViscosity() const {
    const T temperature = std::max(parameters.temperature, T(1.));
    const T referenceTemperature =
        std::max(parameters.referenceTemperature, T(1.));
    return parameters.referenceViscosity *
           std::exp(parameters.creepActivationEnergy / gasConstant *
                    (T(1) / temperature - T(1) / referenceTemperature));
  }

  // Lamé viscosity parameters: same structure as elastic Lamé constants but
  // with eta(T) in place of Young's modulus, making the governing equation a
  // viscous Stokes flow rather than elastic equilibrium.
  T lameMu() const {
    const T nu = clampedPoissonRatio();
    return effectiveMaskViscosity() / (T(2) * (T(1) + nu));
  }

  T lameLambda() const {
    const T nu = clampedPoissonRatio();
    return effectiveMaskViscosity() * nu /
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

};

/// Velocity field for the oxide outer interface when part of that interface is
/// in contact with a mask. Open oxide/ambient regions use the oxide deformation
/// field. Mask-contact regions are kinematically constrained to the solved mask
/// vector velocity and do not receive additional free-surface scalar growth.
template <class T, int D>
class OxidationConstrainedAmbient final : public VelocityField<T> {
  using IndexType = viennahrle::Index<D>;
  using ConstSparseIterator =
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>;

  SmartPointer<OxidationDeformation<T, D>> deformationField =
      nullptr;
  SmartPointer<OxidationMaskBending<T, D>> maskVelocityField =
      nullptr;
  SmartPointer<Domain<T, D>> maskInterface = nullptr;
  int maskSign = 1;

public:
  OxidationConstrainedAmbient() = default;

  OxidationConstrainedAmbient(
      SmartPointer<OxidationDeformation<T, D>> passedDeformation,
      SmartPointer<OxidationMaskBending<T, D>> passedMaskVelocity,
      SmartPointer<Domain<T, D>> passedMaskInterface, int passedMaskSign = 1)
      : deformationField(passedDeformation),
        maskVelocityField(passedMaskVelocity), maskInterface(passedMaskInterface),
        maskSign((passedMaskSign < 0) ? -1 : 1) {}

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<OxidationConstrainedAmbient>::New(
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
