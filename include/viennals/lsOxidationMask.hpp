#pragma once

#include <lsOxidationSolverBase.hpp>
#include <lsOxidationDeformation.hpp>

#include <hrleSparseStarIterator.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <omp.h>
#include <stdexcept>
#include <unordered_map>

namespace viennals {

template <class T> struct OxidationMaskParameters {
  // Contact mode:
  //   0 = kinematic: mask velocity at oxide-contact faces equals the solved
  //       oxide velocity (legacy Dirichlet BC). No traction computation.
  //   1 = oneway: mask solved with oxide stress traction as Neumann BC
  //       (multigrid GMRES); oxide uses kinematic Dirichlet at mask faces.
  //       Stable; captures mask bending driven by oxide force. (Config
  //       aliases: "oneway", "traction", "1", "2".)
  //   2 = elastic: multigrid-GMRES elastic solve with youngModulus and
  //       poissonRatio. Oxide contact prescribes elastic displacement
  //       (v_oxide × dt). The outer coupling loop in lsOxidation feeds the
  //       solved mask displacement back into the next oxide solve. (Config
  //       aliases: "elastic", "twoway" and variants, "3", "4".)
  int contactMode = 1;
  // Arrhenius creep viscosity law:
  // eta(T) = eta_ref * exp(E/R * (1/T - 1/T_ref)).
  // Viscosity is in Pa·hr, activation energy in J/mol, temperature in K.
  T temperature = 1273.15;
  T referenceTemperature = 1273.15;
  T referenceViscosity = 5e8;
  T creepActivationEnergy = 0.;
  // Young's modulus for elastic contact mode (2), in Pa.
  // Si₃N₄ at 1000 °C: ~250–270 GPa.
  T youngModulus = T(270e9);
  // Current solve timestep in hours; set by lsOxidation before each apply().
  // Elastic modes use it to scale oxide velocity to contact displacement
  // (v * dt) and to convert the solved displacement back to advection velocity.
  T stressTimeStep = T(1);
  T poissonRatio = 0.27;
  // false = bonded mask/oxide interface (traction continuity in compression
  // and tension); true = optional unilateral contact/release model.
  bool unilateralContact = true;
  // Outer Aitken relaxation factor applied to the mask/oxide coupling residual.
  // Independent of the multigrid smoother below.
  T relaxation = 1.;
  // Under-relaxation for the unilateral contact load active set.  A hard
  // compressive/tensile switch makes the mask/oxide fixed point oscillate when
  // contact faces release; relaxing the load keeps the complementarity limit
  // but makes the iteration continuous.
  T contactLoadRelaxation = 0.25;
  // Relative pressure floor for releasing a relaxed contact face.  Machine
  // epsilon is far too small for Pa-scale contact stresses: a tensile face with
  // a decaying old compressive load would otherwise remain "active" for many
  // nonlinear iterations.  The scale is the previous/current normal traction.
  T contactReleaseFraction = 5e-3;
  // SOR omega for the multigrid V-cycle smoother (both forward and backward
  // sweeps).  1.0 is standard Gauss-Seidel; values in (1, 1.4] add over-
  // relaxation.  Do not share this with the Aitken relaxation above.
  T multigridSmootherOmega = 1.0;
  T tolerance = 1e-8;
  // Minimum sub-grid boundary distance as a fraction of gridDelta. This is a
  // mechanical derivative length scale, so keep it comparable to the oxide
  // deformation solver; near-zero distances make elastic contact stresses blow
  // up as E * u / d.
  T minBoundaryDistance = 0.05;
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

  struct SparseMatrix {
    std::size_t nodeCount = 0;
    std::vector<std::size_t> rowPtr;
    std::vector<std::size_t> colIndex;
    std::vector<T> values;
    std::vector<T> invDiagonal;
  };

  struct MultigridLevel {
    std::vector<IndexType> indices;
    // For level l > 0, children maps each coarse node to level l-1 nodes.
    std::vector<std::vector<std::size_t>> children;
    std::vector<std::size_t> fineToCoarse;
    SparseMatrix matrix;
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
  T lastApplyAbsoluteVelocityChange = std::numeric_limits<T>::max();
  T aitkenOmega = 1.;
  std::size_t candidateContactFaces_ = 0;
  std::size_t activeContactFaces_ = 0;
  std::size_t tensileContactFaces_ = 0;
  T minContactNormalTraction_ = std::numeric_limits<T>::max();
  T maxContactNormalTraction_ = std::numeric_limits<T>::lowest();
  T contactReleaseThreshold_ = 0.;
  std::unordered_map<std::size_t, Vec3D<T>> previousContactTraction_;
  std::unordered_map<std::size_t, T> previousContactReleaseScale_;
  std::vector<T> previousAitkenResidual;
  // Elastic mode (contactMode==2): snapshot of u_new (elastic equilibrium
  // displacement in µm, stored as µm/hr with implicit dt_ref=1 hr) taken by
  // finalizeElasticAdvectionVelocity(). Used by writeFieldsToLevelSet() to
  // write the "MaskVelocity" warm-start for the next substep's solver.
  std::vector<Vec3D<T>> elasticU_;
  IndexType requestedMinIndex{};
  IndexType requestedMaxIndex{};
  std::vector<Node> nodes;
  // Face-major flat contact BC arrays: index = faceIdx * n + nodeId.
  std::vector<uint8_t> contactFaceActive_;      // 1 if contact BC applies
  std::vector<Vec3D<T>> contactFaceVelocity_;   // Dirichlet velocity at contact
  std::vector<Vec3D<T>> contactFaceTraction_;   // Oxide traction on mask face
  std::vector<T> contactFaceDistance_;          // Node-to-interface distance
  std::unordered_map<std::size_t, T> ambientPhiCache_;

  // Cached multigrid hierarchy.  The stiffness matrix depends on the node
  // geometry and the contact face classification (active/inactive) but NOT on
  // the traction magnitudes, which only affect the load vector b.  When the
  // contact pattern is unchanged from the previous apply() call, the hierarchy
  // can be reused and the O(n) matrix build and Galerkin coarsening skipped.
  std::vector<MultigridLevel> cachedMultigridLevels_;
  std::vector<uint8_t>        cachedContactFaceActive_; // snapshot at last build
  std::size_t                 cachedNodeCount_ = 0;

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

  ~OxidationMaskBending() = default;

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
  T getLastApplyAbsoluteVelocityChange() const {
    return lastApplyAbsoluteVelocityChange;
  }

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
    elasticU_.clear(); // reset before coupling so getVectorVelocity divides by dt
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
    smoothVelocityField();

    const auto fixedPointResidual =
        contactResidualVector(previousVelocities);
    const T omega = aitkenRelaxation(fixedPointResidual);
    relaxVelocities(previousVelocities, omega);
    validateNodeVelocities("Aitken relaxation");
    lastApplyVelocityChange = maxVelocityChange(previousVelocities);
    lastApplyAbsoluteVelocityChange =
        absoluteVelocityChange(previousVelocities);

    maxVelocity_.fill(T(0));
    for (const auto &node : nodes) {
      for (unsigned d = 0; d < D; ++d)
        maxVelocity_[d] = std::max(maxVelocity_[d], std::abs(node.velocity[d]));
    }
    if (Logger::hasDebug()) {
      Logger::getInstance()
          .addDebug("OxidationMaskBending: contact faces active=" +
                    std::to_string(activeContactFaces_) + "/" +
                    std::to_string(candidateContactFaces_) +
                    ", tensileSkipped=" +
                    std::to_string(tensileContactFaces_) +
                    ", normalTraction=[" +
                    std::to_string(candidateContactFaces_ == 0 ? T(0)
                                                              : minContactNormalTraction_) +
                    ", " +
                    std::to_string(candidateContactFaces_ == 0 ? T(0)
                                                              : maxContactNormalTraction_) +
                    "], aitkenOmega=" + std::to_string(omega) +
                    ", contactLoadRelaxation=" +
                    std::to_string(std::clamp(parameters.contactLoadRelaxation,
                                              T(0.02), T(1))) +
                    ", contactReleaseThreshold=" +
                    std::to_string(contactReleaseThreshold_) +
                    ", residual=" + std::to_string(lastApplyVelocityChange) +
                    ", absVelocityChange=" +
                    std::to_string(lastApplyAbsoluteVelocityChange))
          .print();
    }

    solved = true;
  }

  // Called from lsOxidation after writeFieldsToLevelSet() and before advect().
  //
  // Elastic mode (2): the contact faces use kinematic (Dirichlet) BC so the
  // mask stays bonded to the oxide surface.  The solver output u_new is the
  // elastic bending displacement (in µm, stored as µm/hr with implicit
  // dt_ref=1 hr).  Convert to advection velocity u_new/dt and update
  // maxVelocity_ for CFL subcycling.  No-op for viscous modes.
  void finalizeElasticAdvectionVelocity() {
    if (!isElasticContactMode() || nodes.empty())
      return;
    const T dt = parameters.stressTimeStep;
    if (dt <= T(0)) {
      for (auto &node : nodes)
        node.velocity = {T(0), T(0), T(0)};
      return;
    }

    const std::size_t n = nodes.size();
    elasticU_.resize(n);
    for (std::size_t i = 0; i < n; ++i)
      elasticU_[i] = nodes[i].velocity;

    T maxUNew = T(0);
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned d = 0; d < D; ++d)
        maxUNew = std::max(maxUNew, std::abs(elasticU_[i][d]));

    const T invDt = T(1) / dt;
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned d = 0; d < D; ++d)
        nodes[i].velocity[d] = elasticU_[i][d] * invDt;

    Logger::getInstance()
        .addInfo("ElasticFinalize: dt=" + std::to_string(dt) +
                 " max|u_new|=" + std::to_string(maxUNew) +
                 " max|u_new/dt|=" + std::to_string(maxUNew * invDt))
        .print();

    maxVelocity_.fill(T(0));
    for (const auto &node : nodes)
      for (unsigned d = 0; d < D; ++d)
        maxVelocity_[d] = std::max(maxVelocity_[d], std::abs(node.velocity[d]));
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

    // Elastic mode: the solver stores u_new (displacement in µm) as if it were
    // a velocity (µm/hr, with implicit dt_ref=1hr).  The oxide needs the actual
    // advection velocity u_new/dt.  Before finalization elasticU_ is empty and
    // nodes[i].velocity = u_new; we divide by dt here.  After finalization
    // nodes[i].velocity = u_new/dt already, so we return it directly.
    if (isElasticContactMode() && elasticU_.empty()) {
      const T dt = parameters.stressTimeStep;
      if (dt <= T(0))
        return {T(0), T(0), T(0)};
      return detail::vecScaled(getVelocity(index), T(1) / dt);
    }

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

    // Elastic mode: after finalization, elasticU_ holds u_new (the displacement
    // in µm, stored as µm/hr). Use it as the "MaskVelocity" warm-start so the
    // next substep's solver begins near its expected solution.
    // Before finalization (elasticU_ empty), fall back to nodes[nId].velocity
    // which still equals u_new from the current solve.
    const bool useElasticU =
        isElasticContactMode() && !elasticU_.empty();

    VD velocity;
    ConstSparseIterator it(maskInterface->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const std::size_t nId = lookupNode(it.getStartIndices());
      Vec3D<T> v{T(0), T(0), T(0)};
      if (nId != noNode)
        v = (useElasticU && nId < elasticU_.size()) ? elasticU_[nId]
                                                     : nodes[nId].velocity;
      velocity.push_back(v);
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

  T absoluteVelocityChange(
      const std::unordered_map<std::size_t, Vec3D<T>> &previous) const {
    if (previous.empty())
      return std::numeric_limits<T>::max();

    T changeSquaredSum = 0.;
    std::size_t components = 0;
    for (const auto &node : nodes) {
      if (!node.contact)
        continue;
      const auto found = previous.find(linearIndex(node.index));
      if (found == previous.end())
        continue;

      for (unsigned i = 0; i < D; ++i) {
        const T delta = node.velocity[i] - found->second[i];
        if (!std::isfinite(delta))
          throwNonFinite("mask absolute velocity coupling residual");
        changeSquaredSum += delta * delta;
        ++components;
      }
    }

    if (components == 0)
      return std::numeric_limits<T>::max();
    const T change =
        std::sqrt(changeSquaredSum / static_cast<T>(components));
    if (!std::isfinite(change))
      throwNonFinite("mask absolute velocity coupling residual");
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
    const T baseOmega = std::clamp(parameters.relaxation, T(0.01), T(1));
    if (residualVector.empty()) {
      previousAitkenResidual.clear();
      aitkenOmega = baseOmega;
      return 1.;
    }

    T omega = std::clamp(aitkenOmega, T(0.01), baseOmega);
    if (previousAitkenResidual.size() != residualVector.size())
      omega = baseOmega;
    if (previousAitkenResidual.size() == residualVector.size()) {
      T numerator = 0.;
      T denominator = 0.;
      T residualNorm2 = 0.;
      T previousNorm2 = 0.;
      for (std::size_t i = 0; i < residualVector.size(); ++i) {
        if (!std::isfinite(residualVector[i]) ||
            !std::isfinite(previousAitkenResidual[i]))
          throwNonFinite("mask Aitken residual");
        const T delta = residualVector[i] - previousAitkenResidual[i];
        numerator += previousAitkenResidual[i] * delta;
        denominator += delta * delta;
        residualNorm2 += residualVector[i] * residualVector[i];
        previousNorm2 += previousAitkenResidual[i] * previousAitkenResidual[i];
      }

      if (!std::isfinite(numerator) || !std::isfinite(denominator))
        throwNonFinite("mask Aitken coefficient");

      if (denominator > std::numeric_limits<T>::epsilon()) {
        omega = -aitkenOmega * numerator / denominator;
        if (std::isfinite(omega)) {
          // Traction contact mode: cap at 1.0 (no extrapolation).  The
          // compressive/tensile contact state can alternate across coupling
          // iterations, and omega > 1 extrapolates through the fixed point
          // for those nodes, producing sign-alternating velocities that
          // appear as zigzag kinks after level-set advection.
          const T omegaMax =
              (parameters.contactMode > 0) ? baseOmega : T(1.5);
          omega = std::clamp(omega, T(0.05), omegaMax);
        } else {
          throwNonFinite("mask Aitken coefficient");
        }
      }

      if (std::isfinite(residualNorm2) && std::isfinite(previousNorm2) &&
          residualNorm2 > previousNorm2 * T(1.05)) {
        omega = std::min(omega, std::max(T(0.05), aitkenOmega * T(0.5)));
      }
    }

    previousAitkenResidual = residualVector;
    aitkenOmega = omega;
    return omega;
  }

  // One pass of Laplacian (neighbour-average) smoothing on the mask velocity
  // field.  Damps grid-scale oscillations that arise near the contact-to-free
  // boundary transition without materially changing the bulk velocity.
  // Fixed (anchor) nodes are skipped — they must stay at zero.
  void smoothVelocityField() {
    const std::size_t n = nodes.size();
    if (n == 0)
      return;

    std::vector<Vec3D<T>> smoothed(n);
    for (std::size_t id = 0; id < n; ++id) {
      if (nodes[id].fixed) {
        smoothed[id] = {T(0), T(0), T(0)};
        continue;
      }
      Vec3D<T> sum = nodes[id].velocity;
      int count = 1;
      for (unsigned dir = 0; dir < D; ++dir) {
        for (int off : {-1, 1}) {
          IndexType nb = nodes[id].index;
          nb[dir] += off;
          if (!inBounds(nb))
            continue;
          const std::size_t nbId = nodeLookupFlat[linearIndex(nb)];
          if (nbId == noNode)
            continue;
          for (unsigned c = 0; c < D; ++c)
            sum[c] += nodes[nbId].velocity[c];
          ++count;
        }
      }
      const T w = T(1) / static_cast<T>(count);
      for (unsigned c = 0; c < D; ++c)
        smoothed[id][c] = sum[c] * w;
    }

    for (std::size_t id = 0; id < n; ++id)
      nodes[id].velocity = smoothed[id];
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
    auto oldContactTraction = std::move(previousContactTraction_);
    auto oldContactReleaseScale = std::move(previousContactReleaseScale_);
    previousContactTraction_.clear();
    previousContactReleaseScale_.clear();
    contactReleaseThreshold_ = T(0);
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
    candidateContactFaces_ = 0;
    activeContactFaces_ = 0;
    tensileContactFaces_ = 0;
    minContactNormalTraction_ = std::numeric_limits<T>::max();
    maxContactNormalTraction_ = std::numeric_limits<T>::lowest();
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
          ++candidateContactFaces_;
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
          minContactNormalTraction_ =
              std::min(minContactNormalTraction_, tn);
          maxContactNormalTraction_ =
              std::max(maxContactNormalTraction_, tn);

          Vec3D<T> contactLoad = t;
          if (parameters.unilateralContact && tn >= T(0)) {
            ++tensileContactFaces_;
            contactLoad = {T(0), T(0), T(0)};
          }

          if (parameters.contactMode > 0) {
            const auto key = contactFaceKey(node.index, dir, offset);
            const auto prevIt = oldContactTraction.find(key);
            if (prevIt != oldContactTraction.end()) {
              const T alpha =
                  std::clamp(parameters.contactLoadRelaxation, T(0.02), T(1));
              for (unsigned c = 0; c < D; ++c)
                contactLoad[c] =
                    prevIt->second[c] +
                    alpha * (contactLoad[c] - prevIt->second[c]);
            }
          }

          T loadNormal = T(0);
          for (unsigned i = 0; i < D; ++i)
            loadNormal += contactLoad[i] * faceNormal[i];

          T releaseThreshold = std::numeric_limits<T>::epsilon();
          T releaseScale = std::max(std::abs(tn), T(1));
          if (parameters.contactMode > 0 && parameters.unilateralContact) {
            const auto scaleIt =
                oldContactReleaseScale.find(contactFaceKey(node.index, dir,
                                                           offset));
            if (scaleIt != oldContactReleaseScale.end())
              releaseScale = std::max(releaseScale, scaleIt->second);
            releaseThreshold =
                std::clamp(parameters.contactReleaseFraction, T(0), T(0.25)) *
                releaseScale;
            contactReleaseThreshold_ =
                std::max(contactReleaseThreshold_, releaseThreshold);
          }

          if (parameters.unilateralContact &&
              loadNormal >= -releaseThreshold)
            continue;

          if (!isFinite(contactLoad))
            throwNonFinite("relaxed oxide contact load");

          if (parameters.contactMode > 0) {
            const auto key = contactFaceKey(node.index, dir, offset);
            previousContactTraction_[key] = contactLoad;
            if (parameters.unilateralContact)
              previousContactReleaseScale_[key] =
                  std::max(releaseScale, std::abs(loadNormal));
          }

          contactFaceActive_[faceIdx * n + id]    = uint8_t(1);
          ++activeContactFaces_;
          contactFaceTraction_[faceIdx * n + id]  = contactLoad;
          contactFaceDistance_[faceIdx * n + id]  = faceDistance;
          if (usesKinematicContactBoundary()) {
            auto oxVel = deformationField->getVectorVelocity(oxidePt,
                                                             parameters.material,
                                                             faceNormal, 0);
            if (isElasticContactMode()) {
              // Scale by dt: elastic solver uses dt_ref=1hr so the Dirichlet BC
              // must be the physical displacement (v_oxide × dt), not the velocity.
              oxVel = detail::vecScaled(oxVel, parameters.stressTimeStep);
            }
            contactFaceVelocity_[faceIdx * n + id] = oxVel;
          }
        }
      }
    }
    markFixedNodes();
  }

  std::size_t contactFaceKey(const IndexType &index, unsigned direction,
                             int offset) const {
    std::size_t seed = detail::gridIndexHash<D>(index);
    seed ^= std::hash<unsigned>{}(direction) +
            std::size_t(0x9e3779b97f4a7c15ULL) + (seed << 6) + (seed >> 2);
    seed ^= std::hash<int>{}(offset) +
            std::size_t(0x9e3779b97f4a7c15ULL) + (seed << 6) + (seed >> 2);
    return seed;
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
      if (usesKinematicContactBoundary())
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

  static constexpr std::size_t mgNoNode =
      std::numeric_limits<std::size_t>::max();

  static Vec3D<T> zeroVec() { return Vec3D<T>{T(0), T(0), T(0)}; }

  static IndexType coarsenIndex(const IndexType &index) {
    IndexType coarse{};
    for (unsigned d = 0; d < D; ++d) {
      const auto value = index[d];
      coarse[d] = (value >= 0) ? value / 2 : -((-value + 1) / 2);
    }
    return coarse;
  }

  std::vector<MultigridLevel> buildMultigridHierarchy() const {
    std::vector<MultigridLevel> levels;
    levels.emplace_back();
    auto &fine = levels.back();
    fine.indices.reserve(nodes.size());
    for (const auto &node : nodes)
      fine.indices.push_back(node.index);

    constexpr unsigned maxLevels = 12;
    constexpr std::size_t minCoarsenNodes = 48;
    while (levels.size() < maxLevels &&
           levels.back().indices.size() > minCoarsenNodes) {
      const auto &previous = levels.back();
      MultigridLevel coarse;
      coarse.indices.reserve(previous.indices.size() / 2 + 1);
      coarse.children.reserve(previous.indices.size() / 2 + 1);
      coarse.fineToCoarse.assign(previous.indices.size(), mgNoNode);

      std::unordered_map<std::size_t, std::size_t> coarseLookup;
      coarseLookup.reserve(previous.indices.size());
      for (std::size_t fineId = 0; fineId < previous.indices.size(); ++fineId) {
        const IndexType coarseIndex = coarsenIndex(previous.indices[fineId]);
        const std::size_t key = detail::gridIndexHash<D>(coarseIndex);
        auto found = coarseLookup.find(key);
        if (found == coarseLookup.end()) {
          const std::size_t coarseId = coarse.indices.size();
          coarseLookup.emplace(key, coarseId);
          coarse.indices.push_back(coarseIndex);
          coarse.children.emplace_back();
          found = coarseLookup.find(key);
        }

        const std::size_t coarseId = found->second;
        coarse.children[coarseId].push_back(fineId);
        coarse.fineToCoarse[fineId] = coarseId;
      }

      if (coarse.indices.empty() ||
          coarse.indices.size() >= previous.indices.size())
        break;

      // Stop if any grid dimension collapses to a single cell on the coarse
      // level.  For thin masks this happens quickly in the thickness direction,
      // and a 1-cell-thick level loses all stress-gradient information across
      // that dimension, making the smoother degenerate.
      {
        IndexType lo = coarse.indices.front(), hi = coarse.indices.front();
        for (const auto &idx : coarse.indices)
          for (unsigned d = 0; d < D; ++d) {
            lo[d] = std::min(lo[d], idx[d]);
            hi[d] = std::max(hi[d], idx[d]);
          }
        bool degenerate = false;
        for (unsigned d = 0; d < D; ++d)
          if (hi[d] == lo[d])
            degenerate = true;
        if (degenerate)
          break;
      }

      levels.push_back(std::move(coarse));
    }
    return levels;
  }

  static T vectorDot(const std::vector<Vec3D<T>> &a,
                     const std::vector<Vec3D<T>> &b) {
    T sum = T(0);
    for (std::size_t i = 0; i < a.size(); ++i)
      for (unsigned c = 0; c < D; ++c)
        sum += a[i][c] * b[i][c];
    return sum;
  }

  static T vectorNorm(const std::vector<Vec3D<T>> &v) {
    return std::sqrt(std::max(vectorDot(v, v), T(0)));
  }

  static void vectorScale(std::vector<Vec3D<T>> &v, T scale) {
    for (auto &entry : v)
      for (unsigned c = 0; c < D; ++c)
        entry[c] *= scale;
  }

  static void vectorAxpy(std::vector<Vec3D<T>> &y, T alpha,
                         const std::vector<Vec3D<T>> &x) {
    for (std::size_t i = 0; i < y.size(); ++i)
      for (unsigned c = 0; c < D; ++c)
        y[i][c] += alpha * x[i][c];
  }

  static void vectorSubtractInPlace(std::vector<Vec3D<T>> &y, T alpha,
                                    const std::vector<Vec3D<T>> &x) {
    for (std::size_t i = 0; i < y.size(); ++i)
      for (unsigned c = 0; c < D; ++c)
        y[i][c] -= alpha * x[i][c];
  }

  static T flatValue(const std::vector<Vec3D<T>> &v, std::size_t row) {
    return v[row / D][row % D];
  }

  static void flatSet(std::vector<Vec3D<T>> &v, std::size_t row, T value) {
    v[row / D][row % D] = value;
  }

  void collectLocalCandidatesRecursive(unsigned dim, const IndexType &center,
                                       IndexType &candidate,
                                       std::vector<std::size_t> &out) const {
    if (dim == D) {
      if (!inBounds(candidate))
        return;
      const std::size_t nodeId = lookupNode(candidate);
      if (nodeId != noNode)
        out.push_back(nodeId);
      return;
    }

    for (int offset = -2; offset <= 2; ++offset) {
      candidate[dim] = center[dim] + offset;
      collectLocalCandidatesRecursive(dim + 1, center, candidate, out);
    }
  }

  void collectLocalCandidates(std::size_t rowNode,
                              std::vector<std::size_t> &out) const {
    out.clear();
    IndexType candidate = nodes[rowNode].index;
    collectLocalCandidatesRecursive(0, nodes[rowNode].index, candidate, out);
    out.push_back(rowNode);
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
  }

  SparseMatrix compressSparseRows(
      std::vector<std::unordered_map<std::size_t, T>> &rows,
      std::size_t nodeCount) const {
    SparseMatrix matrix;
    matrix.nodeCount = nodeCount;
    matrix.rowPtr.assign(rows.size() + 1, 0);
    matrix.invDiagonal.assign(rows.size(), T(1));

    for (std::size_t row = 0; row < rows.size(); ++row) {
      std::vector<std::pair<std::size_t, T>> entries(rows[row].begin(),
                                                     rows[row].end());
      std::sort(entries.begin(), entries.end(),
                [](const auto &a, const auto &b) { return a.first < b.first; });

      matrix.rowPtr[row] = matrix.values.size();
      T diagonal = T(0);
      for (const auto &[col, value] : entries) {
        if (std::abs(value) <= std::numeric_limits<T>::epsilon() * T(100))
          continue;
        if (!std::isfinite(value))
          throwNonFinite("traction sparse matrix assembly");
        matrix.colIndex.push_back(col);
        matrix.values.push_back(value);
        if (col == row)
          diagonal += value;
      }

      if (std::abs(diagonal) > std::numeric_limits<T>::epsilon())
        matrix.invDiagonal[row] = T(1) / diagonal;
      else
        matrix.invDiagonal[row] = T(1);
    }
    matrix.rowPtr[rows.size()] = matrix.values.size();
    return matrix;
  }

  SparseMatrix buildFineElasticMatrix(const std::vector<Vec3D<T>> &b,
                                      T gradDivWeight) const {
    const std::size_t n = nodes.size();
    std::vector<std::unordered_map<std::size_t, T>> rows(n * D);
    std::vector<Vec3D<T>> basis(n, zeroVec());
    std::vector<std::size_t> candidates;

    for (std::size_t rowNode = 0; rowNode < n; ++rowNode) {
      if (nodes[rowNode].fixed) {
        for (unsigned rowComponent = 0; rowComponent < D; ++rowComponent) {
          const std::size_t row = rowNode * D + rowComponent;
          rows[row][row] = T(1);
        }
        continue;
      }

      // Build A column-by-column via probing.  The operator is (A·v)[i][c] =
      // v[i][c] - F(v)[i][c], where F is affine: F(v) = F_linear(v) + F(0).
      // b[rowNode] = F(0)[rowNode] (the traction load vector).  The matrix
      // entry is the linear part only:
      //   A[row, col] = δ(row,col) - (F(e_col)[rowComp] - F(0)[rowComp])
      //               = δ(row,col) - (Fbasis[rowComp] - b[rowNode][rowComp])
      // b is constant with respect to col so it cancels between probes; it
      // is included here to isolate the linear part of F from its affine offset.
      collectLocalCandidates(rowNode, candidates);
      for (std::size_t colNode : candidates) {
        for (unsigned colComponent = 0; colComponent < D; ++colComponent) {
          basis[colNode][colComponent] = T(1);
          const Vec3D<T> Fbasis =
              computeElasticStencilAt(rowNode, basis, gradDivWeight);
          basis[colNode][colComponent] = T(0);

          for (unsigned rowComponent = 0; rowComponent < D; ++rowComponent) {
            const std::size_t row = rowNode * D + rowComponent;
            const std::size_t col = colNode * D + colComponent;
            const T identity =
                (rowNode == colNode && rowComponent == colComponent) ? T(1)
                                                                     : T(0);
            const T value = identity - (Fbasis[rowComponent] -
                                        b[rowNode][rowComponent]);
            rows[row][col] += value;
          }
        }
      }
    }

    return compressSparseRows(rows, n);
  }

  SparseMatrix buildGalerkinMatrix(const SparseMatrix &fine,
                                   const MultigridLevel &coarseLevel) const {
    const std::size_t coarseNodes = coarseLevel.indices.size();
    std::vector<std::unordered_map<std::size_t, T>> rows(coarseNodes * D);

    for (std::size_t fineRow = 0; fineRow < fine.nodeCount * D; ++fineRow) {
      const std::size_t fineRowNode = fineRow / D;
      if (fineRowNode >= coarseLevel.fineToCoarse.size())
        continue;
      const std::size_t coarseRowNode = coarseLevel.fineToCoarse[fineRowNode];
      if (coarseRowNode == mgNoNode)
        continue;

      const T restrictionWeight =
          T(1) / static_cast<T>(
                     std::max<std::size_t>(coarseLevel.children[coarseRowNode].size(),
                                           1));
      const std::size_t coarseRow = coarseRowNode * D + fineRow % D;

      for (std::size_t nz = fine.rowPtr[fineRow]; nz < fine.rowPtr[fineRow + 1];
           ++nz) {
        const std::size_t fineCol = fine.colIndex[nz];
        const std::size_t fineColNode = fineCol / D;
        if (fineColNode >= coarseLevel.fineToCoarse.size())
          continue;
        const std::size_t coarseColNode = coarseLevel.fineToCoarse[fineColNode];
        if (coarseColNode == mgNoNode)
          continue;
        const std::size_t coarseCol = coarseColNode * D + fineCol % D;
        rows[coarseRow][coarseCol] += restrictionWeight * fine.values[nz];
      }
    }

    return compressSparseRows(rows, coarseNodes);
  }

  void sparseMatvec(const SparseMatrix &matrix,
                    const std::vector<Vec3D<T>> &x,
                    std::vector<Vec3D<T>> &Ax) const {
    Ax.assign(matrix.nodeCount, zeroVec());
#pragma omp parallel for schedule(static)
    for (std::size_t row = 0; row < matrix.nodeCount * D; ++row) {
      T sum = T(0);
      for (std::size_t nz = matrix.rowPtr[row]; nz < matrix.rowPtr[row + 1];
           ++nz)
        sum += matrix.values[nz] * flatValue(x, matrix.colIndex[nz]);
      flatSet(Ax, row, sum);
    }
  }

  void multigridProlong(const std::vector<MultigridLevel> &levels,
                        std::size_t coarseLevelId,
                        const std::vector<Vec3D<T>> &coarse,
                        std::vector<Vec3D<T>> &fine) const {
    const auto &fineLevel = levels[coarseLevelId - 1];
    const auto &coarseLevel = levels[coarseLevelId];
    fine.assign(fineLevel.indices.size(), zeroVec());
    for (std::size_t coarseId = 0; coarseId < coarseLevel.children.size();
         ++coarseId)
      for (std::size_t fineId : coarseLevel.children[coarseId])
        for (unsigned c = 0; c < D; ++c)
          fine[fineId][c] = coarse[coarseId][c];
  }

  void multigridSmooth(const SparseMatrix &matrix,
                       const std::vector<Vec3D<T>> &rhs,
                       std::vector<Vec3D<T>> &x,
                       unsigned sweeps,
                       T omega) const {
    const std::size_t rows = matrix.nodeCount * D;
    auto relaxRow = [&](std::size_t row) {
      T offDiagonal = T(0);
      for (std::size_t nz = matrix.rowPtr[row]; nz < matrix.rowPtr[row + 1];
           ++nz) {
        const std::size_t col = matrix.colIndex[nz];
        if (col == row)
          continue;
        offDiagonal += matrix.values[nz] * flatValue(x, col);
      }
      const T updated =
          (flatValue(rhs, row) - offDiagonal) * matrix.invDiagonal[row];
      flatSet(x, row, flatValue(x, row) + omega * (updated - flatValue(x, row)));
    };

    for (unsigned sweep = 0; sweep < sweeps; ++sweep) {
      for (std::size_t row = 0; row < rows; ++row)
        relaxRow(row);
      for (std::size_t row = rows; row-- > 0;)
        relaxRow(row);
    }
  }

  void multigridResidual(const SparseMatrix &matrix,
                         const std::vector<Vec3D<T>> &rhs,
                         const std::vector<Vec3D<T>> &x,
                         std::vector<Vec3D<T>> &residualOut) const {
    sparseMatvec(matrix, x, residualOut);
    for (std::size_t i = 0; i < rhs.size(); ++i)
      for (unsigned c = 0; c < D; ++c)
        residualOut[i][c] = rhs[i][c] - residualOut[i][c];
  }

  void multigridRestrict(const std::vector<Vec3D<T>> &fineResidual,
                         const MultigridLevel &coarseLevel,
                         std::vector<Vec3D<T>> &coarseRhs) const {
    coarseRhs.assign(coarseLevel.indices.size(), zeroVec());
    for (std::size_t coarseId = 0; coarseId < coarseLevel.children.size();
         ++coarseId) {
      const auto &children = coarseLevel.children[coarseId];
      if (children.empty())
        continue;
      for (std::size_t fineId : children)
        for (unsigned c = 0; c < D; ++c)
          coarseRhs[coarseId][c] += fineResidual[fineId][c];
      const T scale = T(1) / static_cast<T>(children.size());
      for (unsigned c = 0; c < D; ++c)
        coarseRhs[coarseId][c] *= scale;
    }
  }

  void multigridProlongAdd(const std::vector<Vec3D<T>> &coarseCorrection,
                            const MultigridLevel &coarseLevel,
                            std::vector<Vec3D<T>> &fineCorrection) const {
    for (std::size_t coarseId = 0; coarseId < coarseLevel.children.size();
         ++coarseId) {
      for (std::size_t fineId : coarseLevel.children[coarseId])
        for (unsigned c = 0; c < D; ++c)
          fineCorrection[fineId][c] += coarseCorrection[coarseId][c];
    }
  }

  void multigridVCycle(const std::vector<MultigridLevel> &levels,
                       std::size_t levelId,
                       const std::vector<Vec3D<T>> &rhs,
                       std::vector<Vec3D<T>> &x,
                       T smootherOmega) const {
    const auto &level = levels[levelId];
    if (levelId + 1 >= levels.size() || level.indices.size() <= 32) {
      multigridSmooth(level.matrix, rhs, x, 40, smootherOmega);
      return;
    }

    multigridSmooth(level.matrix, rhs, x, 2, smootherOmega);

    std::vector<Vec3D<T>> fineResidual;
    multigridResidual(level.matrix, rhs, x, fineResidual);

    std::vector<Vec3D<T>> coarseRhs;
    const auto &coarseLevel = levels[levelId + 1];
    multigridRestrict(fineResidual, coarseLevel, coarseRhs);

    std::vector<Vec3D<T>> coarseCorrection(coarseRhs.size(), zeroVec());
    multigridVCycle(levels, levelId + 1, coarseRhs, coarseCorrection,
                    smootherOmega);
    multigridProlongAdd(coarseCorrection, coarseLevel, x);

    multigridSmooth(level.matrix, rhs, x, 2, smootherOmega);
  }

  std::vector<Vec3D<T>>
  multigridPrecondition(const std::vector<MultigridLevel> &levels,
                        const std::vector<Vec3D<T>> &rhs,
                        T smootherOmega) const {
    std::vector<Vec3D<T>> correction(rhs.size(), zeroVec());
    if (levels.empty())
      return correction;
    multigridVCycle(levels, 0, rhs, correction, smootherOmega);
    return correction;
  }

  void exactElasticResidual(const SparseMatrix &matrix,
                            const std::vector<Vec3D<T>> &x,
                            const std::vector<Vec3D<T>> &b,
                            std::vector<Vec3D<T>> &r) const {
    std::vector<Vec3D<T>> Ax(x.size(), zeroVec());
    sparseMatvec(matrix, x, Ax);
    r.assign(x.size(), zeroVec());
    for (std::size_t i = 0; i < x.size(); ++i)
      for (unsigned c = 0; c < D; ++c)
        r[i][c] = b[i][c] - Ax[i][c];
  }

  static std::vector<T>
  solveUpperTriangular(const std::vector<std::vector<T>> &h,
                       const std::vector<T> &g,
                       unsigned usedColumns) {
    std::vector<T> y(usedColumns, T(0));
    for (int row = static_cast<int>(usedColumns) - 1; row >= 0; --row) {
      T sum = g[static_cast<std::size_t>(row)];
      for (unsigned col = static_cast<unsigned>(row) + 1;
           col < usedColumns; ++col)
        sum -= h[static_cast<std::size_t>(row)][col] * y[col];
      const T diag = h[static_cast<std::size_t>(row)]
                      [static_cast<std::size_t>(row)];
      if (std::abs(diag) > std::numeric_limits<T>::epsilon())
        y[static_cast<std::size_t>(row)] = sum / diag;
    }
    return y;
  }

  void solveElasticVelocity() {
    if (parameters.contactMode > 0) {
      solveElasticVelocityMultigridGMRES();
      return;
    }
    solveElasticVelocityBiCGSTAB();
  }

  void solveElasticVelocityMultigridGMRES() {
    iterations = 0;
    residual = 0.;
    if (nodes.empty())
      return;

    const T lambda = lameLambda();
    const T mu = lameMu();
    const T gradDivWeight =
        (lambda + mu) / std::max(lambda + T(2) * mu,
                                 std::numeric_limits<T>::epsilon());
    const T smootherOmega =
        std::clamp(parameters.multigridSmootherOmega, T(0.2), T(1.4));

    const std::size_t n = nodes.size();
    const std::vector<Vec3D<T>> zeros(n, zeroVec());
    std::vector<Vec3D<T>> b(n, zeroVec());
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < n; ++i) {
      if (nodes[i].fixed)
        b[i] = zeroVec();
      else
        b[i] = computeElasticStencilAt(i, zeros, gradDivWeight);
    }

    std::vector<Vec3D<T>> x(n, zeroVec());
    for (std::size_t i = 0; i < n; ++i)
      x[i] = nodes[i].fixed ? zeroVec() : nodes[i].velocity;

    // Rebuild the multigrid hierarchy only when the node count or the contact
    // face classification changes.  The stiffness matrix depends only on node
    // geometry and which faces are active (compressive), not on the traction
    // magnitudes — those only enter the load vector b computed above.  Within
    // a single time step the coupling iterations run without advection, so the
    // mask geometry and contact pattern are typically stable after the first
    // call and the matrix can be reused for all subsequent iterations.
    const bool hierarchyDirty =
        cachedNodeCount_ != n ||
        cachedContactFaceActive_ != contactFaceActive_;

    if (hierarchyDirty) {
      cachedMultigridLevels_ = buildMultigridHierarchy();
      if (cachedMultigridLevels_.empty())
        return;
      cachedMultigridLevels_[0].matrix =
          buildFineElasticMatrix(b, gradDivWeight);
      for (std::size_t level = 1; level < cachedMultigridLevels_.size();
           ++level)
        cachedMultigridLevels_[level].matrix =
            buildGalerkinMatrix(cachedMultigridLevels_[level - 1].matrix,
                                cachedMultigridLevels_[level]);
      cachedNodeCount_          = n;
      cachedContactFaceActive_  = contactFaceActive_;
    }

    if (cachedMultigridLevels_.empty())
      return;
    const auto &multigridLevels = cachedMultigridLevels_;

    std::vector<Vec3D<T>> r;
    exactElasticResidual(multigridLevels[0].matrix, x, b, r);

    const T rhsNorm = vectorNorm(b);
    const T absTolerance =
        std::max(parameters.tolerance * rhsNorm,
                 std::numeric_limits<T>::epsilon() * T(100) *
                     std::sqrt(static_cast<T>(
                         std::max<std::size_t>(std::size_t(1), n * D))));
    const T residualNormDenom = std::max(rhsNorm, absTolerance);
    auto updateGmresResidual = [&](T absResidual) {
      residual = (absResidual <= absTolerance)
                     ? T(0)
                     : absResidual / residualNormDenom;
    };

    T absResidual = vectorNorm(r);
    updateGmresResidual(absResidual);
    if (absResidual <= absTolerance || residual < parameters.tolerance) {
      for (std::size_t i = 0; i < n; ++i)
        nodes[i].velocity = x[i];
      return;
    }

    constexpr unsigned restart = 32;
    while (iterations < parameters.maxIterations &&
           residual > parameters.tolerance) {
      const T beta = vectorNorm(r);
      if (!std::isfinite(beta))
        throwNonFinite("traction multigrid GMRES residual");
      if (beta <= absTolerance) {
        residual = T(0);
        break;
      }

      const unsigned innerLimit =
          std::min<unsigned>(restart, parameters.maxIterations - iterations);
      std::vector<std::vector<Vec3D<T>>> v(innerLimit + 1,
                                           std::vector<Vec3D<T>>(n, zeroVec()));
      std::vector<std::vector<Vec3D<T>>> z(innerLimit,
                                           std::vector<Vec3D<T>>(n, zeroVec()));
      v[0] = r;
      vectorScale(v[0], T(1) / beta);

      std::vector<std::vector<T>> h(innerLimit + 1,
                                    std::vector<T>(innerLimit, T(0)));
      std::vector<T> cs(innerLimit, T(0));
      std::vector<T> sn(innerLimit, T(0));
      std::vector<T> g(innerLimit + 1, T(0));
      g[0] = beta;

      unsigned usedColumns = 0;
      for (unsigned j = 0; j < innerLimit; ++j) {
        z[j] = multigridPrecondition(multigridLevels, v[j], smootherOmega);

        std::vector<Vec3D<T>> w(n, zeroVec());
        sparseMatvec(multigridLevels[0].matrix, z[j], w);

        for (unsigned i = 0; i <= j; ++i) {
          h[i][j] = vectorDot(w, v[i]);
          vectorSubtractInPlace(w, h[i][j], v[i]);
        }

        h[j + 1][j] = vectorNorm(w);
        if (h[j + 1][j] > std::numeric_limits<T>::epsilon()) {
          v[j + 1] = w;
          vectorScale(v[j + 1], T(1) / h[j + 1][j]);
        }

        for (unsigned i = 0; i < j; ++i) {
          const T h0 = h[i][j];
          const T h1 = h[i + 1][j];
          h[i][j] = cs[i] * h0 + sn[i] * h1;
          h[i + 1][j] = -sn[i] * h0 + cs[i] * h1;
        }

        const T h0 = h[j][j];
        const T h1 = h[j + 1][j];
        const T denom = std::hypot(h0, h1);
        if (denom <= std::numeric_limits<T>::epsilon()) {
          cs[j] = T(1);
          sn[j] = T(0);
        } else {
          cs[j] = h0 / denom;
          sn[j] = h1 / denom;
        }
        h[j][j] = cs[j] * h0 + sn[j] * h1;
        h[j + 1][j] = T(0);

        const T g0 = g[j];
        g[j] = cs[j] * g0;
        g[j + 1] = -sn[j] * g0;

        ++iterations;
        usedColumns = j + 1;
        const T projectedAbsResidual = std::abs(g[j + 1]);
        updateGmresResidual(projectedAbsResidual);
        if (!std::isfinite(residual))
          throwNonFinite("traction multigrid GMRES residual");
        if (projectedAbsResidual <= absTolerance ||
            residual < parameters.tolerance)
          break;
      }

      if (usedColumns == 0)
        break;

      const auto y = solveUpperTriangular(h, g, usedColumns);
      for (unsigned col = 0; col < usedColumns; ++col)
        vectorAxpy(x, y[col], z[col]);

      exactElasticResidual(multigridLevels[0].matrix, x, b, r);
      absResidual = vectorNorm(r);
      updateGmresResidual(absResidual);
      if (!std::isfinite(residual))
        throwNonFinite("traction multigrid GMRES residual");
    }

    for (std::size_t i = 0; i < n; ++i)
      nodes[i].velocity = x[i];

    if (residual > parameters.tolerance)
      Logger::getInstance()
          .addWarning("solveElasticVelocity: traction multigrid GMRES did not "
                      "converge after " + std::to_string(iterations) + "/" +
                      std::to_string(parameters.maxIterations) +
                      " iterations (residual=" + std::to_string(residual) +
                      ", tolerance=" + std::to_string(parameters.tolerance) + ")")
          .print();
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

    // Warn if the anchor is placed along the oxide-growth direction (D-1 in a
    // standard LOCOS cross-section is the vertical/y axis).  Clamping the top
    // or bottom of the mask instead of a far lateral edge removes the degrees
    // of freedom the bending solve needs and will suppress the bird's-beak
    // deflection entirely.
    if (dir == static_cast<unsigned>(D - 1))
      Logger::getInstance()
          .addWarning("OxidationMaskBending: anchorBoundaryDirection=" +
                      std::to_string(parameters.anchorBoundaryDirection) +
                      " points along the oxide-growth axis (direction D-1=" +
                      std::to_string(D - 1) + ").  The anchor is normally "
                      "placed at a lateral edge (direction 0) so bending "
                      "degrees of freedom are preserved.  Set "
                      "anchorBoundaryDirection=0 for a LOCOS cross-section.")
          .print();
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

  bool isElasticContactMode() const {
    return parameters.contactMode == 2;
  }

  bool usesKinematicContactBoundary() const {
    return parameters.contactMode == 0 || isElasticContactMode();
  }

  T effectiveMaskViscosity() const {
    // Elastic mode: η_eff = E (Pa·hr, with implicit dt_ref = 1 hr).
    // lameMu/lameLambda equal the standard elastic Lamé constants G and λ.
    // Contact faces use kinematic (Dirichlet) BC: v_contact = v_oxide × dt,
    // so the solver gives u_new ≈ v_oxide × dt (displacement in µm stored as
    // µm/hr).  finalizeElasticAdvectionVelocity() converts u_new → u_new/dt
    // (the actual advection velocity µm/hr).  Displacement feedback into the
    // oxide comes through the outer coupling loop in lsOxidation.
    if (isElasticContactMode())
      return parameters.youngModulus; // Pa·hr with implicit dt_ref = 1 hr
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
  SmartPointer<Domain<T, D>> ambientInterface = nullptr;
  int maskSign = 1;
  std::unordered_map<std::size_t, T> maskPhiCache_;
  T maskGridDelta_ = 1.;
  std::array<T, D> maxVelocity_{};

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
    buildEffectiveVelocityCache();
  }

  OxidationConstrainedAmbient(
      SmartPointer<OxidationDeformation<T, D>> passedDeformation,
      SmartPointer<OxidationMaskBending<T, D>> passedMaskVelocity,
      SmartPointer<Domain<T, D>> passedMaskInterface,
      SmartPointer<Domain<T, D>> passedAmbientInterface,
      int passedMaskSign = 1)
      : deformationField(passedDeformation),
        maskVelocityField(passedMaskVelocity), maskInterface(passedMaskInterface),
        ambientInterface(passedAmbientInterface),
        maskSign((passedMaskSign < 0) ? -1 : 1) {
    buildMaskPhiCache();
    buildEffectiveVelocityCache();
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

    // Near-contact gap zone: smoothly approach the mask velocity instead of
    // letting stress spikes just outside the mask edge advect the ambient level
    // set with the full oxide deformation velocity.
    if (maskVelocityField != nullptr) {
      if (signedPhi > -T(3) * maskGridDelta_) {
        auto v_mask = maskVelocityField->getVectorVelocity(
            coordinate, material, normalVector, pointId);
        const T blend =
            std::max(T(0), std::min(T(1),
                                    (signedPhi + T(3) * maskGridDelta_) /
                                        (T(3) * maskGridDelta_)));
        for (int k = 0; k < D; ++k)
          v_def[k] = (T(1) - blend) * v_def[k] + blend * v_mask[k];

        T def_n = T(0), mask_n = T(0);
        for (int k = 0; k < D; ++k) {
          def_n  += v_def[k]  * normalVector[k];
          mask_n += v_mask[k] * normalVector[k];
        }
        if (mask_n > def_n) {
          const T boost = mask_n - def_n;
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
    if (direction < 0 || direction >= static_cast<int>(D))
      return T(0);
    (void)material;
    (void)centralDifferences;
    return maxVelocity_[direction];
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

  void buildEffectiveVelocityCache() {
    maxVelocity_.fill(T(0));
    if (ambientInterface == nullptr) {
      for (unsigned d = 0; d < D; ++d) {
        if (deformationField != nullptr)
          maxVelocity_[d] = std::max(
              maxVelocity_[d],
              deformationField->getDissipationAlpha(d, -1, {}));
        if (maskVelocityField != nullptr)
          maxVelocity_[d] = std::max(
              maxVelocity_[d],
              maskVelocityField->getDissipationAlpha(d, -1, {}));
      }
      return;
    }

    const T gridDelta = ambientInterface->getGrid().getGridDelta();
    bool foundInterfacePoint = false;
    viennahrle::ConstSparseStarIterator<typename Domain<T, D>::DomainType, 1>
        neighborIterator(ambientInterface->getDomain());
    for (ConstSparseIterator it(ambientInterface->getDomain());
         !it.isFinished(); ++it) {
      if (!it.isDefined() || std::abs(it.getValue()) > T(1))
        continue;
      foundInterfacePoint = true;

      const auto index = it.getStartIndices();
      neighborIterator.goToIndicesSequential(index);

      Vec3D<T> coordinate{0., 0., 0.};
      Vec3D<T> normal{0., 0., 0.};
      T normalNorm2 = T(0);
      for (unsigned d = 0; d < D; ++d) {
        coordinate[d] = static_cast<T>(index[d]) * gridDelta;
        normal[d] = neighborIterator.getNeighbor(d).getValue() -
                    neighborIterator.getNeighbor(d + D).getValue();
        normalNorm2 += normal[d] * normal[d];
      }
      if (normalNorm2 > std::numeric_limits<T>::epsilon()) {
        const T invNorm = T(1) / std::sqrt(normalNorm2);
        for (unsigned d = 0; d < D; ++d)
          normal[d] *= invNorm;
      }

      const auto vectorVelocity =
          getVectorVelocity(coordinate, -1, normal, 0);
      const T scalarVelocity =
          getScalarVelocity(coordinate, -1, normal, 0);
      for (unsigned d = 0; d < D; ++d) {
        maxVelocity_[d] =
            std::max(maxVelocity_[d],
                     std::abs(vectorVelocity[d] + scalarVelocity * normal[d]));
      }
    }

    if (!foundInterfacePoint) {
      for (unsigned d = 0; d < D; ++d) {
        if (deformationField != nullptr)
          maxVelocity_[d] = std::max(
              maxVelocity_[d],
              deformationField->getDissipationAlpha(d, -1, {}));
        if (maskVelocityField != nullptr)
          maxVelocity_[d] = std::max(
              maxVelocity_[d],
              maskVelocityField->getDissipationAlpha(d, -1, {}));
      }
    }
  }

};

} // namespace viennals
