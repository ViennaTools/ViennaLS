#pragma once

#include <lsOxidationSolverBase.hpp>
#include <lsOxidationDeformation.hpp>

#include <omp.h>
#include <stdexcept>
#include <unordered_map>

namespace viennals {

template <class T> struct OxidationMaskParameters {
  // Contact mode:
  //   0 = kinematic, prescribed mask/oxide contact velocity (legacy).
  //   1 = traction, oxide stress drives the mask through a viscous-elastic
  //       Neumann boundary condition.
  int contactMode = 1;
  // Arrhenius creep viscosity law:
  // eta(T) = eta_ref * exp(E/R * (1/T - 1/T_ref)).
  // Viscosity is in Pa·hr, activation energy in J/mol, temperature in K.
  T temperature = 1273.15;
  T referenceTemperature = 1273.15;
  T referenceViscosity = 5e8;
  T creepActivationEnergy = 0.;
  T poissonRatio = 0.27;
  bool unilateralContact = true;
  T relaxation = 1.;
  T tolerance = 1e-8;
  T minBoundaryDistance = 1e-6;
  unsigned maxIterations = 10000;
  std::size_t maxGridPoints = 5000000;
  int material = -1;
  // Optional far-field clamp used to remove rigid-body mask drift and provide
  // the support that makes mask thickness mechanically meaningful.  A side of
  // -1 clamps the lower index side, +1 clamps the upper side, and 0 disables
  // the clamp.  Direction defaults to x in a LOCOS cross-section.
  int anchorBoundaryDirection = 0;
  int anchorBoundarySide = -1;
  unsigned anchorBoundaryLayers = 1;
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
  using OxidationSolverBase<T, D>::initializeGridFromMask;

  struct Node {
    IndexType index;
    Vec3D<T> velocity{0., 0., 0.};
    bool contact = false;
    bool fixed = false;
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
  std::size_t fixedNodes = 0;
  T lastApplyVelocityChange = std::numeric_limits<T>::max();
  T aitkenOmega = 1.;
  std::vector<T> previousAitkenResidual;
  IndexType requestedMinIndex{};
  IndexType requestedMaxIndex{};
  std::vector<Node> nodes;
  // Face-major flat contact BC arrays: index = faceIdx * n + nodeId.
  std::vector<uint8_t> contactFaceActive_;   // 1 if contact BC applies
  std::vector<Vec3D<T>> contactFaceVelocity_; // Legacy Dirichlet velocity
  std::vector<Vec3D<T>> contactFaceTraction_; // Oxide traction on mask face
  std::vector<T> contactFaceDistance_;        // Node-to-interface distance
  std::unordered_map<std::size_t, T> ambientPhiCache_;

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
  std::size_t getNumberOfFixedNodes() const { return fixedNodes; }
  T getLastApplyVelocityChange() const { return lastApplyVelocityChange; }

  void apply() {
    if (maskInterface == nullptr) {
      solved = true;
      return;
    }
    if (deformationField == nullptr) {
      Logger::getInstance()
          .addWarning("OxidationMaskBending: deformation field is null; "
                      "mask bending will produce zero velocities.")
          .print();
      solved = true;
      return;
    }

    const auto previousVelocities = collectVelocitiesByGridPoint();
    if (!initialiseGrid())
      return; // base class already logged the error
    buildNodes();
    seedFromLevelSet(); // warm-start velocity from previous substep's pointData
    seedFromPrevious(previousVelocities);
    if (nodes.empty()) {
      Logger::getInstance()
          .addWarning("OxidationMaskBending: no mask nodes found after "
                      "buildNodes(). Verify that the mask level set defines "
                      "a non-empty interior region within the solve bounds.")
          .print();
      solved = true;
      return;
    }
    if (contactNodes == 0)
      Logger::getInstance()
          .addDebug("OxidationMaskBending: mask has no oxide-contact nodes; "
                    "no traction will be applied this step.")
          .print();
    solveElasticVelocity();
    validateNodeVelocities("solveElasticVelocity");

    const auto fixedPointResidual =
        contactResidualVector(previousVelocities);
    const T omega = aitkenRelaxation(fixedPointResidual);
    relaxVelocities(previousVelocities, omega);
    validateNodeVelocities("Aitken relaxation");
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

  /// Write mask bending velocity into maskInterface->getPointData() so that
  /// lsInterior + lsAdvect carry it across timestep boundaries.
  void writeFieldsToLevelSet() {
    if (nodes.empty() || maskInterface == nullptr)
      return;

    using VD = typename PointData<T>::VectorDataType;
    VD velocity;

    ConstSparseIterator it(maskInterface->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const IndexType idx = it.getStartIndices();
      const std::size_t nId = lookupNode(idx);
      velocity.push_back(nId != noNode
                             ? nodes[nId].velocity
                             : Vec3D<T>{T(0), T(0), T(0)});
    }
    maskInterface->getPointData().insertReplaceVectorData(
        std::move(velocity), "MaskVelocity");
  }

private:
  bool initialiseGrid() {
    return initializeGridFromMask(maskInterface, useRequestedBounds,
                                  requestedMinIndex, requestedMaxIndex,
                                  parameters.maxGridPoints,
                                  "OxidationMaskBending");
  }

  /// Seed nodes[i].velocity from maskInterface->getPointData()["MaskVelocity"]
  /// so solveElasticVelocity() warm-starts from the previous substep's field.
  void seedFromLevelSet() {
    if (maskInterface == nullptr || nodes.empty())
      return;
    const int vIdx = maskInterface->getPointData().getVectorDataIndex("MaskVelocity");
    if (vIdx == -1)
      return;
    const auto *vd = maskInterface->getPointData().getVectorData(vIdx);
    if (vd == nullptr)
      return;

    ConstSparseIterator it(maskInterface->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const auto ptId = it.getPointId();
      if (ptId >= static_cast<decltype(ptId)>(vd->size()))
        continue;
      const std::size_t nId = lookupNode(it.getStartIndices());
      if (nId == noNode)
        continue;
      if (!isFinite((*vd)[ptId]))
        throwNonFinite("stored mask velocity point data");
      else
        nodes[nId].velocity = (*vd)[ptId];
    }
  }

  void seedFromPrevious(
      const std::unordered_map<std::size_t, Vec3D<T>> &previous) {
    if (previous.empty())
      return;
    for (auto &node : nodes) {
      const auto found = previous.find(linearIndex(node.index));
      if (found != previous.end() && isFinite(found->second))
        node.velocity = found->second;
      if (node.fixed)
        node.velocity = {T(0), T(0), T(0)};
    }
  }

  std::unordered_map<std::size_t, Vec3D<T>> collectVelocitiesByGridPoint() const {
    std::unordered_map<std::size_t, Vec3D<T>> result;
    result.reserve(nodes.size());
    for (const auto &node : nodes)
      if (inBounds(node.index) && isFinite(node.velocity))
        result.emplace(linearIndex(node.index), node.velocity);
    return result;
  }

  static bool isFinite(const Vec3D<T> &v) {
    for (unsigned i = 0; i < 3; ++i)
      if (!std::isfinite(v[i]))
        return false;
    return true;
  }

  static bool isFiniteTensor(const std::array<T, 9> &tensor) {
    for (T value : tensor)
      if (!std::isfinite(value))
        return false;
    return true;
  }

  [[noreturn]] void throwNonFinite(const std::string &stage) const {
    const std::string message =
        "OxidationMaskBending: " + stage + " produced non-finite values.";
    Logger::getInstance().addError(message).print();
    throw std::runtime_error(message);
  }

  void validateNodeVelocities(const std::string &stage) const {
    for (const auto &node : nodes)
      for (unsigned i = 0; i < D; ++i)
        if (!std::isfinite(node.velocity[i]))
          throwNonFinite(stage);
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

    const T change = std::sqrt(changeSquaredSum / magnitudeSquaredSum);
    if (!std::isfinite(change))
      throwNonFinite("mask velocity coupling residual");
    return change;
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
      for (unsigned i = 0; i < D; ++i) {
        const T residualComponent = node.velocity[i] - found->second[i];
        if (!std::isfinite(residualComponent))
          throwNonFinite("mask contact residual");
        residualVector.push_back(residualComponent);
      }
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
        if (!std::isfinite(residualVector[i]) ||
            !std::isfinite(previousAitkenResidual[i]))
          throwNonFinite("mask Aitken residual");
        const T delta = residualVector[i] - previousAitkenResidual[i];
        numerator += previousAitkenResidual[i] * delta;
        denominator += delta * delta;
      }

      if (!std::isfinite(numerator) || !std::isfinite(denominator))
        throwNonFinite("mask Aitken coefficient");

      if (denominator > std::numeric_limits<T>::epsilon()) {
        omega = -aitkenOmega * numerator / denominator;
        if (std::isfinite(omega))
          omega = std::clamp(omega, T(0.05), T(1.5));
        else
          throwNonFinite("mask Aitken coefficient");
      }
    }

    previousAitkenResidual = residualVector;
    aitkenOmega = omega;
    return omega;
  }

  void relaxVelocities(const std::unordered_map<std::size_t, Vec3D<T>> &previous,
                       T omega) {
    if (!std::isfinite(omega))
      throwNonFinite("mask Aitken coefficient");
    if (previous.empty() || omega == T(1))
      return;

    for (auto &node : nodes) {
      const auto found = previous.find(linearIndex(node.index));
      if (found == previous.end())
        continue;
      for (unsigned i = 0; i < D; ++i) {
        const T relaxed =
            found->second[i] + omega * (node.velocity[i] - found->second[i]);
        if (!std::isfinite(relaxed))
          throwNonFinite("mask velocity relaxation");
        node.velocity[i] = relaxed;
      }
    }
  }

  void buildAmbientPhiCache() {
    ambientPhiCache_.clear();
    if (ambientInterface == nullptr)
      return;
    for (ConstSparseIterator it(ambientInterface->getDomain());
         !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      ambientPhiCache_[detail::gridIndexHash<D>(it.getStartIndices())] =
          static_cast<T>(ambientSign) * it.getValue();
    }
  }

  void buildNodes() {
    nodes.clear();
    initNodeLookup();
    buildAmbientPhiCache();

    ConstSparseIterator maskIt(maskInterface->getDomain());
    IndexType index = minIndex;
    while (true) {
      if (isInsideMask(maskIt, index)) {
        const std::size_t id = nodes.size();
        nodeLookupFlat[linearIndex(index)] = id;
        nodes.push_back({index});
      }

      if (!increment(index))
        break;
    }

    const std::size_t n = nodes.size();
    contactFaceActive_.assign(2 * D * n, uint8_t(0));
    contactFaceVelocity_.assign(2 * D * n, Vec3D<T>{T(0), T(0), T(0)});
    contactFaceTraction_.assign(2 * D * n, Vec3D<T>{T(0), T(0), T(0)});
    contactFaceDistance_.assign(2 * D * n, gridDelta);

    contactNodes = 0;
    for (std::size_t id = 0; id < n; ++id) {
      auto &node = nodes[id];
      node.contact = touchesContactBoundary(maskIt, node.index);
      if (node.contact)
        ++contactNodes;

      for (unsigned dir = 0; dir < D; ++dir) {
        for (int offset : {-1, 1}) {
          const unsigned faceIdx = dir * 2u + (offset == 1 ? 1u : 0u);
          IndexType neighbor = node.index;
          neighbor[dir] += offset;
          const bool neighborIsNode =
              inBounds(neighbor) && lookupNode(neighbor) != noNode;
          if (neighborIsNode ||
              !isContactBoundary(node.index, dir, offset, maskIt))
            continue; // inactive already set by assign()

          const T faceDistance = maskFaceDistance(maskIt, node.index, neighbor);
          Vec3D<T> faceNormal{0., 0., 0.};
          faceNormal[dir] = static_cast<T>(offset);
          Vec3D<T> oxidePt{0., 0., 0.};
          for (unsigned i = 0; i < D; ++i)
            oxidePt[i] = node.index[i] * gridDelta;
          oxidePt[dir] += offset * gridDelta;

          const auto sigma = deformationField->getStressTensor(oxidePt);
          if (!isFiniteTensor(sigma))
            throwNonFinite("oxide stress at mask contact");
          Vec3D<T> t{0., 0., 0.};
          for (unsigned i = 0; i < D; ++i)
            for (unsigned j = 0; j < D; ++j)
              t[i] += sigma[3 * i + j] * faceNormal[j];

          T tn = 0.;
          for (unsigned i = 0; i < D; ++i)
            tn += t[i] * faceNormal[i];

          if (!isFinite(t) || !std::isfinite(tn))
            throwNonFinite("oxide traction at mask contact");

          if (parameters.unilateralContact && tn >= T(0))
            continue; // tensile contact — leave inactive

          contactFaceActive_[faceIdx * n + id]    = uint8_t(1);
          contactFaceTraction_[faceIdx * n + id]  = t;
          contactFaceDistance_[faceIdx * n + id]  = faceDistance;
          if (parameters.contactMode == 0) {
            contactFaceVelocity_[faceIdx * n + id] =
                deformationField->getVectorVelocity(oxidePt, parameters.material,
                                                    faceNormal, 0);
          }
        }
      }
    }
    markFixedNodes();
  }

  // Returns neighbor velocity without ghost extrapolation.
  // All arithmetic is in T; reads from the SolverT work vector are widened.
  template <class SolverT>
  Vec3D<T> simpleNeighborVelocity(const std::vector<Vec3D<SolverT>> &velocity,
                                   std::size_t nodeId, unsigned direction,
                                   int offset) const {
    IndexType neighbor = nodes[nodeId].index;
    neighbor[direction] += offset;
    const auto toT = [](const Vec3D<SolverT> &v) -> Vec3D<T> {
      return {static_cast<T>(v[0]), static_cast<T>(v[1]), static_cast<T>(v[2])};
    };
    if (!inBounds(neighbor))
      return toT(velocity[nodeId]);
    const std::size_t foundId = nodeLookupFlat[linearIndex(neighbor)];
    return (foundId != noNode) ? toT(velocity[foundId]) : toT(velocity[nodeId]);
  }

  // Ghost velocity enforcing sigma_mask * n = traction at a mask boundary face.
  template <class SolverT>
  Vec3D<T> stressBoundaryGhost(const std::vector<Vec3D<SolverT>> &velocity,
                               std::size_t nodeId, unsigned normalDir,
                               int offset, T distance,
                               const Vec3D<T> &traction) const {
    const T lambda = lameLambda();
    const T mu = lameMu();
    const T denom = std::max(lambda + T(2) * mu,
                             std::numeric_limits<T>::epsilon());
    const T signedOffset = static_cast<T>(offset);
    Vec3D<T> ghost{static_cast<T>(velocity[nodeId][0]),
                   static_cast<T>(velocity[nodeId][1]),
                   static_cast<T>(velocity[nodeId][2])};
    T tangentialDivergence = T(0);
    std::array<T, D> normalTangentialDerivative{};
    for (unsigned tanDir = 0; tanDir < D; ++tanDir) {
      if (tanDir == normalDir)
        continue;
      const Vec3D<T> vPlus  = simpleNeighborVelocity(velocity, nodeId, tanDir, +1);
      const Vec3D<T> vMinus = simpleNeighborVelocity(velocity, nodeId, tanDir, -1);
      tangentialDivergence +=
          (vPlus[tanDir] - vMinus[tanDir]) / (T(2) * gridDelta);
      normalTangentialDerivative[tanDir] =
          (vPlus[normalDir] - vMinus[normalDir]) / (T(2) * gridDelta);
    }

    const T normalTraction = traction[normalDir] * signedOffset;
    const T normalDerivative =
        (normalTraction - lambda * tangentialDivergence) / denom;
    ghost[normalDir] += signedOffset * distance * normalDerivative;

    for (unsigned tanDir = 0; tanDir < D; ++tanDir) {
      if (tanDir == normalDir)
        continue;
      const T normalDerivativeTangential =
          signedOffset * traction[tanDir] /
              std::max(mu, std::numeric_limits<T>::epsilon()) -
          normalTangentialDerivative[tanDir];
      ghost[tanDir] += signedOffset * distance * normalDerivativeTangential;
    }
    return ghost;
  }

  template <class SolverT>
  Vec3D<T> tractionFreeGhost(const std::vector<Vec3D<SolverT>> &velocity,
                              std::size_t nodeId, unsigned normalDir,
                              int offset) const {
    return stressBoundaryGhost(velocity, nodeId, normalDir, offset, gridDelta,
                               Vec3D<T>{T(0), T(0), T(0)});
  }

  template <class SolverT>
  Vec3D<T> neighborVelocity(const std::vector<Vec3D<SolverT>> &velocity,
                             std::size_t nodeId, unsigned direction,
                             int offset) const {
    IndexType neighbor = nodes[nodeId].index;
    neighbor[direction] += offset;

    if (inBounds(neighbor)) {
      const std::size_t foundId = nodeLookupFlat[linearIndex(neighbor)];
      if (foundId != noNode)
        return {static_cast<T>(velocity[foundId][0]),
                static_cast<T>(velocity[foundId][1]),
                static_cast<T>(velocity[foundId][2])};
    }

    const unsigned faceIdx = direction * 2u + (offset == 1 ? 1u : 0u);
    const std::size_t nn = nodes.size();
    if (contactFaceActive_[faceIdx * nn + nodeId]) {
      if (parameters.contactMode == 0)
        return detail::vecSubtract(
            detail::vecScaled(contactFaceVelocity_[faceIdx * nn + nodeId], T(2)),
            Vec3D<T>{static_cast<T>(velocity[nodeId][0]),
                     static_cast<T>(velocity[nodeId][1]),
                     static_cast<T>(velocity[nodeId][2])});

      return stressBoundaryGhost(
          velocity, nodeId, direction, offset,
          contactFaceDistance_[faceIdx * nn + nodeId],
          contactFaceTraction_[faceIdx * nn + nodeId]);
    }

    return tractionFreeGhost(velocity, nodeId, direction, offset);
  }

  template <class SolverT>
  Vec3D<T> divergenceGradient(const std::vector<Vec3D<SolverT>> &velocity,
                               const IndexType &index) const {
    Vec3D<T> gradient{T(0), T(0), T(0)};
    for (unsigned component = 0; component < D; ++component) {
      IndexType plus = index;
      IndexType minus = index;
      plus[component]  += 1;
      minus[component] -= 1;
      const T plusDiv  = divergence(velocity, plus);
      const T minusDiv = divergence(velocity, minus);
      gradient[component] = (plusDiv - minusDiv) / (T(2) * gridDelta);
    }
    return gradient;
  }

  template <class SolverT>
  T divergence(const std::vector<Vec3D<SolverT>> &velocity,
               const IndexType &index) const {
    if (!inBounds(index))
      return T(0);
    const std::size_t nodeId = nodeLookupFlat[linearIndex(index)];
    if (nodeId == noNode)
      return T(0);

    T result = T(0);
    for (unsigned direction = 0; direction < D; ++direction) {
      const Vec3D<T> plus  = neighborVelocity(velocity, nodeId, direction,  1);
      const Vec3D<T> minus = neighborVelocity(velocity, nodeId, direction, -1);
      result += (plus[direction] - minus[direction]) / (T(2) * gridDelta);
    }
    return result;
  }

  // Evaluates the elastic stencil F(v)[i] = laplaceAverage + gradDivCorrection.
  // (A·v)[i] = v[i] - F(v)[i] + b[i],  where b[i] = F(0)[i] (contact BC constants).
  template <class SolverT>
  Vec3D<T> computeElasticStencilAt(std::size_t nodeId,
                                    const std::vector<Vec3D<SolverT>> &v,
                                    T gradDivWeight) const {
    Vec3D<T> laplaceAverage{T(0), T(0), T(0)};
    for (unsigned direction = 0; direction < D; ++direction)
      for (int offset : {-1, 1})
        detail::vecAddTo(laplaceAverage, neighborVelocity(v, nodeId, direction, offset));

    const T count = static_cast<T>(2 * D);
    const Vec3D<T> lapAvg = detail::vecScaled(laplaceAverage, T(1) / count);
    const Vec3D<T> gradDivCorr = detail::vecScaled(
        divergenceGradient(v, nodes[nodeId].index),
        gradDivWeight * gridDelta * gridDelta / (T(2) * static_cast<T>(D)));

    return detail::vecAdd(lapAvg, gradDivCorr);
  }

  // (Av)[i] = v[i] - F(v)[i] + b[i], stored as SolverT.
  template <class SolverT>
  void elasticMatvec(const std::vector<Vec3D<SolverT>> &v,
                     const std::vector<Vec3D<T>> &b,
                     T gradDivWeight,
                     std::vector<Vec3D<SolverT>> &Av) const {
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      if (nodes[i].fixed) {
        for (unsigned c = 0; c < D; ++c)
          Av[i][c] = v[i][c];
        continue;
      }
      const Vec3D<T> Fv = computeElasticStencilAt(i, v, gradDivWeight);
      for (unsigned c = 0; c < D; ++c)
        Av[i][c] = static_cast<SolverT>(static_cast<T>(v[i][c]) - Fv[c] + b[i][c]);
    }
  }

  void solveElasticVelocity() {
    if (parameters.contactMode == 1) {
      solveElasticVelocityRelaxation();
      return;
    }
    solveElasticVelocityBiCGSTAB();
  }

  void solveElasticVelocityRelaxation() {
    iterations = 0;
    residual = 0.;
    if (nodes.empty())
      return;

    const T lambda = lameLambda();
    const T mu = lameMu();
    const T gradDivWeight =
        (lambda + mu) / std::max(lambda + T(2) * mu,
                                 std::numeric_limits<T>::epsilon());
    const T relaxation =
        std::clamp(parameters.relaxation, T(0.01), T(1));

    std::vector<Vec3D<T>> current(nodes.size());
    std::vector<Vec3D<T>> next(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i)
      current[i] = nodes[i].fixed ? Vec3D<T>{T(0), T(0), T(0)}
                                  : nodes[i].velocity;

    for (; iterations < parameters.maxIterations; ++iterations) {
      T maxDelta = T(0);
      T maxMagnitude = std::numeric_limits<T>::epsilon();
      int finiteFlag = 1;

#pragma omp parallel for schedule(static) reduction(max:maxDelta,maxMagnitude) reduction(min:finiteFlag)
      for (std::size_t i = 0; i < nodes.size(); ++i) {
        Vec3D<T> candidate =
            nodes[i].fixed ? Vec3D<T>{T(0), T(0), T(0)}
                           : computeElasticStencilAt(i, current, gradDivWeight);

        for (unsigned c = 0; c < D; ++c) {
          if (!std::isfinite(candidate[c]) || !std::isfinite(current[i][c])) {
            finiteFlag = 0;
            next[i][c] = current[i][c];
            continue;
          }
          next[i][c] =
              current[i][c] + relaxation * (candidate[c] - current[i][c]);
          maxDelta =
              std::max(maxDelta, std::abs(next[i][c] - current[i][c]));
          maxMagnitude = std::max(maxMagnitude, std::abs(next[i][c]));
          maxMagnitude = std::max(maxMagnitude, std::abs(current[i][c]));
        }
      }

      if (!finiteFlag) {
        residual = std::numeric_limits<T>::infinity();
        throwNonFinite("traction mask solve");
      }

      residual = maxDelta / maxMagnitude;
      current.swap(next);
      if (residual < parameters.tolerance) {
        ++iterations;
        break;
      }
    }

    for (std::size_t i = 0; i < nodes.size(); ++i)
      nodes[i].velocity = current[i];

    if (residual > parameters.tolerance)
      Logger::getInstance()
          .addWarning("solveElasticVelocity: traction relaxation did not "
                      "converge after " + std::to_string(iterations) + "/" +
                      std::to_string(parameters.maxIterations) +
                      " iterations (residual=" + std::to_string(residual) +
                      ", tolerance=" + std::to_string(parameters.tolerance) + ")")
          .print();
  }

  void solveElasticVelocityBiCGSTAB() {
    iterations = 0;
    residual = 0.;
    if (nodes.empty())
      return;

    using SolverT = float;

    const T lambda = lameLambda();
    const T mu = lameMu();
    const T gradDivWeight =
        (lambda + mu) / std::max(lambda + T(2) * mu,
                                 std::numeric_limits<T>::epsilon());

    const std::size_t n = nodes.size();
    const Vec3D<SolverT> zero3{SolverT(0), SolverT(0), SolverT(0)};

    // b[i] = F(0)[i]: contact BC constants (reaction/contact velocities).
    // OOB and traction-free faces contribute v[i]=0 at zeros.
    std::vector<Vec3D<T>> b(n);
    {
      const std::vector<Vec3D<SolverT>> zeros(n, zero3);
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < n; ++i) {
        if (nodes[i].fixed)
          b[i] = Vec3D<T>{T(0), T(0), T(0)};
        else
          b[i] = computeElasticStencilAt(i, zeros, gradDivWeight);
      }
    }

    // Cold start from zeros (matching original Jacobi behavior).
    std::vector<Vec3D<SolverT>> x(n, zero3);

    // r = b - A*x. With x=0: A*0 = 0 - F(0) + b = 0, so r = b.
    std::vector<Vec3D<SolverT>> r(n), r_hat(n);
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned c = 0; c < D; ++c) {
        r[i][c]     = static_cast<SolverT>(b[i][c]);
        r_hat[i][c] = r[i][c];
      }

    // BiCGSTAB with identity preconditioner.
    // For most interior nodes the diagonal of A is ≈ 1; the identity preconditioner
    // is exact there and a reasonable approximation at boundary nodes.
    std::vector<Vec3D<SolverT>> pv(n, zero3), sv(n, zero3), y(n), z(n), s(n), t(n);
    T rho = T(1), alpha = T(1), omega = T(1);

    auto vecDot = [&](const std::vector<Vec3D<SolverT>> &a,
                      const std::vector<Vec3D<SolverT>> &bv) {
      T sum = T(0);
      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c) {
          const T av = static_cast<T>(a[i][c]);
          const T bvVal = static_cast<T>(bv[i][c]);
          if (!std::isfinite(av) || !std::isfinite(bvVal))
            return std::numeric_limits<T>::quiet_NaN();
          sum += av * bvVal;
        }
      return sum;
    };

    auto vecMaxAbs = [&](const std::vector<Vec3D<SolverT>> &vin) {
      T m = T(0);
      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c) {
          const T value = static_cast<T>(vin[i][c]);
          if (!std::isfinite(value))
            return std::numeric_limits<T>::infinity();
          m = std::max(m, std::abs(value));
        }
      return m;
    };

    const T b_norm = [&]{ T m = T(0); for (std::size_t i = 0; i < n; ++i) for (unsigned c = 0; c < D; ++c) m = std::max(m, std::abs(b[i][c])); return (m < T(1e-100)) ? T(1) : m; }();

    for (; iterations < parameters.maxIterations; ++iterations) {
      const T rho_new = vecDot(r_hat, r);
      if (!std::isfinite(rho_new) || std::abs(rho_new) < T(1e-100))
        break;
      if (!std::isfinite(rho) || !std::isfinite(alpha) ||
          !std::isfinite(omega) || std::abs(omega) < T(1e-100))
        break;

      const T beta = (rho_new / rho) * (alpha / omega);
      if (!std::isfinite(beta))
        break;
      rho = rho_new;

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          pv[i][c] = static_cast<SolverT>(r[i][c] + beta * (pv[i][c] - omega * sv[i][c]));

      // Identity preconditioner: y = p
      y = pv;
      elasticMatvec(y, b, gradDivWeight, sv);

      const T r_hat_v = vecDot(r_hat, sv);
      if (!std::isfinite(r_hat_v) || std::abs(r_hat_v) < T(1e-100))
        break;

      alpha = rho_new / r_hat_v;
      if (!std::isfinite(alpha))
        break;

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          s[i][c] = static_cast<SolverT>(r[i][c] - alpha * sv[i][c]);

      residual = vecMaxAbs(s);
      if (!std::isfinite(residual))
        break;
      if (residual < parameters.tolerance * b_norm) {
        for (std::size_t i = 0; i < n; ++i)
          for (unsigned c = 0; c < D; ++c)
            x[i][c] = static_cast<SolverT>(x[i][c] + alpha * y[i][c]);
        ++iterations;
        break;
      }

      // Identity preconditioner: z = s
      z = s;
      elasticMatvec(z, b, gradDivWeight, t);

      const T t_s = vecDot(t, s);
      const T t_t = vecDot(t, t);
      if (!std::isfinite(t_s) || !std::isfinite(t_t) || t_t <= T(1e-100))
        break;
      omega = t_s / t_t;
      if (!std::isfinite(omega))
        break;

      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c) {
          x[i][c] = static_cast<SolverT>(x[i][c] + alpha * y[i][c] + omega * z[i][c]);
          r[i][c]  = static_cast<SolverT>(s[i][c] - omega * t[i][c]);
        }

      residual = vecMaxAbs(r);
      if (!std::isfinite(residual))
        break;
      if (residual < parameters.tolerance * b_norm) {
        ++iterations;
        break;
      }
    }

    bool finiteSolution = true;
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned c = 0; c < D; ++c)
        if (!std::isfinite(static_cast<T>(x[i][c])))
          finiteSolution = false;

    if (finiteSolution) {
      for (std::size_t i = 0; i < n; ++i)
        for (unsigned c = 0; c < D; ++c)
          nodes[i].velocity[c] = static_cast<T>(x[i][c]);
    } else {
      residual = std::numeric_limits<T>::infinity();
      throwNonFinite("legacy kinematic mask solve");
    }
    if (residual > parameters.tolerance * b_norm)
      Logger::getInstance()
          .addWarning("solveElasticVelocity: BiCGSTAB did not converge after " +
                      std::to_string(iterations) + "/" +
                      std::to_string(parameters.maxIterations) +
                      " iterations (residual=" + std::to_string(residual / b_norm) +
                      ", tolerance=" + std::to_string(parameters.tolerance) + ")")
          .print();
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
        if (lookupNode(neighbor) != noNode)
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
    if (ambientInterface == nullptr)
      return false;
    const auto it = ambientPhiCache_.find(detail::gridIndexHash<D>(index));
    return it != ambientPhiCache_.end() && it->second >= T(0);
  }

  T maskFaceDistance(ConstSparseIterator &maskIt, const IndexType &inside,
                     const IndexType &outside) const {
    if (!inBounds(outside))
      return gridDelta;
    return crossingDistance(valueAt(maskIt, inside), valueAt(maskIt, outside));
  }

  void markFixedNodes() {
    fixedNodes = 0;
    if (nodes.empty() || parameters.anchorBoundarySide == 0)
      return;
    if (parameters.anchorBoundaryDirection < 0 ||
        parameters.anchorBoundaryDirection >= D)
      return;

    const unsigned dir =
        static_cast<unsigned>(parameters.anchorBoundaryDirection);
    IndexType nodeMin = nodes.front().index;
    IndexType nodeMax = nodes.front().index;
    for (const auto &node : nodes) {
      for (unsigned d = 0; d < D; ++d) {
        nodeMin[d] = std::min(nodeMin[d], node.index[d]);
        nodeMax[d] = std::max(nodeMax[d], node.index[d]);
      }
    }

    const auto layers =
        static_cast<viennahrle::IndexType>(
            std::max(1u, parameters.anchorBoundaryLayers));
    for (auto &node : nodes) {
      const bool onLower =
          parameters.anchorBoundarySide < 0 &&
          node.index[dir] <= nodeMin[dir] + layers - 1;
      const bool onUpper =
          parameters.anchorBoundarySide > 0 &&
          node.index[dir] >= nodeMax[dir] - layers + 1;
      if (onLower || onUpper) {
        node.fixed = true;
        node.velocity = {T(0), T(0), T(0)};
        ++fixedNodes;
      }
    }
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
    const std::size_t nodeId = lookupNode(index);
    if (nodeId != noNode)
      return nodes[nodeId].velocity;

    const auto nearby = findNearbyNode(index);
    if (nearby == noNode)
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
  std::unordered_map<std::size_t, T> maskPhiCache_;
  T maskGridDelta_ = 1.;

public:
  OxidationConstrainedAmbient() = default;

  OxidationConstrainedAmbient(
      SmartPointer<OxidationDeformation<T, D>> passedDeformation,
      SmartPointer<OxidationMaskBending<T, D>> passedMaskVelocity,
      SmartPointer<Domain<T, D>> passedMaskInterface, int passedMaskSign = 1)
      : deformationField(passedDeformation),
        maskVelocityField(passedMaskVelocity), maskInterface(passedMaskInterface),
        maskSign((passedMaskSign < 0) ? -1 : 1) {
    buildMaskPhiCache();
  }

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<OxidationConstrainedAmbient>::New(
        std::forward<Args>(args)...);
  }

  Vec3D<T> getVectorVelocity(const Vec3D<T> &coordinate, int material,
                             const Vec3D<T> &normalVector,
                             unsigned long pointId) final {
    const T signedPhi = getSignedMaskPhi(coordinate);

    // At/inside mask surface: kinematic constraint — oxide tracks mask exactly.
    if (signedPhi >= T(0) && maskVelocityField != nullptr)
      return maskVelocityField->getVectorVelocity(coordinate, material,
                                                  normalVector, pointId);
    if (deformationField == nullptr)
      return {0., 0., 0.};

    auto v_def = deformationField->getVectorVelocity(coordinate, material,
                                                     normalVector, pointId);

    // Gap zone (0 < b < Δx, i.e. signedPhi ∈ (−Δx, 0)): enforce no-separation.
    // Project both velocities onto the oxide outward normal; if the mask is
    // moving away faster than the deformation field, boost the deformation
    // velocity in the normal direction to match.  This prevents the gap from
    // widening without touching the level-set geometry, so there are no
    // corner-snap artefacts at the bird's beak edge.
    if (maskVelocityField != nullptr) {
      // Use a 3-cell zone so gaps that have grown to 1–2 grid cells are still
      // caught.  The boost is safe at this radius: at the bird's beak
      // def_n >> mask_n, so target_n < def_n and the branch is never taken.
      if (signedPhi > -T(3) * maskGridDelta_) {
        auto v_mask = maskVelocityField->getVectorVelocity(
            coordinate, material, normalVector, pointId);
        T def_n = T(0), mask_n = T(0);
        for (int k = 0; k < D; ++k) {
          def_n  += v_def[k]  * normalVector[k];
          mask_n += v_mask[k] * normalVector[k];
        }
        // Target: oxide normal velocity exceeds mask by a closing fraction of
        // the local velocity scale.  This drives the oxide into the mask so
        // the RELATIVE_COMPLEMENT clip always snaps the surfaces co-planar
        // rather than leaving any persistent void.
        const T ref_v = std::max(std::abs(mask_n), std::abs(def_n));
        const T target_n = mask_n + T(1.) * ref_v;
        if (target_n > def_n) {
          const T boost = target_n - def_n;
          Vec3D<T> v_out = v_def;
          for (int k = 0; k < D; ++k)
            v_out[k] += boost * normalVector[k];
          return v_out;
        }
      }
    }
    return v_def;
  }

  T getScalarVelocity(const Vec3D<T> &coordinate, int material,
                      const Vec3D<T> &normalVector,
                      unsigned long pointId) final {
    if (getSignedMaskPhi(coordinate) >= T(0))
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
  // Returns maskSign * maskPhi at the grid node nearest to coordinate.
  // Positive → inside mask (contact), negative → outside mask (gap or free).
  // Uses a pre-built cache so no HRLE iterator is constructed in the hot path.
  // Nodes outside the narrow band return lowest() (unambiguously outside mask).
  T getSignedMaskPhi(const Vec3D<T> &coordinate) const {
    if (maskInterface == nullptr)
      return std::numeric_limits<T>::lowest();
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / maskGridDelta_);
    const auto it = maskPhiCache_.find(detail::gridIndexHash<D>(index));
    return (it != maskPhiCache_.end()) ? it->second
                                       : std::numeric_limits<T>::lowest();
  }

  void buildMaskPhiCache() {
    maskPhiCache_.clear();
    if (maskInterface == nullptr)
      return;
    maskGridDelta_ = maskInterface->getGrid().getGridDelta();
    for (ConstSparseIterator it(maskInterface->getDomain());
         !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const auto key = detail::gridIndexHash<D>(it.getStartIndices());
      maskPhiCache_[key] = static_cast<T>(maskSign) * it.getValue();
    }
  }

};

} // namespace viennals
