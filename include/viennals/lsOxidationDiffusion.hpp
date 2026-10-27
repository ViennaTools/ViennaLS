#pragma once

#include <lsOxidationSolverBase.hpp>
#include <lsVelocityField.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <vcTimer.hpp>

#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/solver/bicgstab.hpp>

#include <omp.h>

#ifdef VIENNALS_GPU_BICGSTAB
#include <lsOxidationBiCGSTABInterface.hpp>
#endif

namespace viennals {

/// Selects the BiCGSTAB back-end for the diffusion solve.
/// GPU failures are reported and not silently fallen back to CPU.
enum class GpuMode {
  Cpu, ///< Always use CPU (default)
  Gpu, ///< Always use GPU; fail if unavailable or unsuccessful
  /// Use the GPU when it is usable, otherwise warn and fall back to the CPU.
  /// This covers the GPU being *unavailable*: no CUDA runtime installed, no
  /// device, or device setup failing. A GPU that is present and initialises
  /// but then fails mid-solve is still reported as an error, because that
  /// indicates a malfunction rather than an absent GPU.
  Auto
};

/// Selects the preconditioner used by the GPU BiCGSTAB solver.
/// Jacobi is the default. ILU0 is experimental because cuSPARSE SpSV
/// value-analysis behavior differs between CUDA versions and can introduce
/// structured solver error if stale.
enum class GpuPreconditioner { Jacobi, ILU0 };

/// Parameters for the steady oxidant diffusion model used by
/// OxidationDiffusion.
struct OxidationParameters {
  double diffusionCoefficient = 1.;
  double reactionRate = 1.;
  double transferCoefficient = 1.;
  double equilibriumConcentration = 1.;
  double oxidantMoleculeDensity = 1.;
  double expansionCoefficient = 1.;
  double velocitySign = 1.;

  // Stress coupling for reaction rate (Kao et al., IEEE TED 35, 25 (1988);
  // the "Vr" of SUPREM-IV and TCAD): the new oxide expands against the
  // stress NORMAL to the Si/SiO2 interface,
  // k_eff = k * exp((sigma_nn + p_ref) * V_k / (k_B * T)),
  // with sigma_nn = n.sigma.n the oxide's Cauchy stress (tension positive) at
  // the interface, so normal compression suppresses the rate.  This is not
  // the pressure: at a convex surface the interface is in normal compression
  // while the oxide pressure there can be tensile.
  // Activation volumes are in m^3, stresses in Pa, temperature in K.
  double temperature = 1273.15;
  double reactionActivationVolume = 0.;
  double referencePressure = 0.;

  // Stress coupling for diffusion coefficient:
  // D_eff = D * exp(-(p - p_ref) * V_D / (k_B * T)).
  double diffusionActivationVolume = 0.;

  // Crystal orientation factor on reaction rate.
  // k(n) = k * [1 + (reactionRateRatio111 - 1) * (1 - (n . crystalAxis)^2)]
  // reactionRateRatio111 = 1 disables the correction (isotropic).
  // crystalAxis is the wafer normal: the default is the 2-D convention (y
  // up); a 3-D domain with z up needs (0, 0, 1) (ViennaPS sets the last
  // axis).
  double reactionRateRatio111 = 1.;
  Vec3D<double> crystalAxis = {0., 1., 0.};

  double maskTransferCoefficient = 0.;
  double maskConcentration = 0.;
  /// Meshed mask diffusion: gamma = D_SiO2 / D_mask.  When > 0 the mask BODY
  /// becomes part of the diffusion domain — oxidant enters its outer surface
  /// through the ambient Robin BC, diffuses through the film with D/gamma,
  /// crosses into the oxide, and reacts at the Si interface.  The film's
  /// resistance then comes from its actual geometry (a thinned or stretched
  /// mask leaks more), not from a lumped coefficient.  0 keeps the mask
  /// impermeable (or lumped via maskTransferCoefficient, if that is set).
  double maskPermeabilityRatio = 0.;
  double minBoundaryDistance = 1e-6;
  unsigned maxIterations = 10000;
  double tolerance = 1e-8;
  double relaxation = 1.;
  std::size_t maxGridPoints = 5000000;
  int material = -1;
};

/// Solves the oxidant diffusion step of the Suvorov et al.
/// (10.1007/s10825-006-0003-z) oxidation model on the Cartesian grid carrying
/// two level sets.
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
  using OxidationSolverBase<T, D>::stencilOffset;
  using OxidationSolverBase<T, D>::linearIndex;
  using OxidationSolverBase<T, D>::increment;
  using OxidationSolverBase<T, D>::findNearbyNode;
  using OxidationSolverBase<T, D>::initializeGridFromInterfaces;

  enum class Boundary { NONE, REACTION, AMBIENT, MASK };

  struct Node {
    IndexType index;
    T concentration = 0.;
    Vec3D<T> siNormal = {0., 0.,
                         0.}; // unit outward normal of Si surface (into oxide)
    // Node lies inside the mask body (meshed mask diffusion): its diffusion
    // coefficient is D / maskPermeabilityRatio instead of D.
    bool inMask = false;
  };

  struct StencilSide {
    T distance = 1.;
    T nodeCoefficient = 1.;
    T constant = 0.;
  };

  SmartPointer<Domain<T, D>> reactionInterface = nullptr;
  SmartPointer<Domain<T, D>> ambientInterface = nullptr;
  SmartPointer<Domain<T, D>> maskInterface = nullptr;
  OxidationParameters parameters;
  int reactionSign = 1;
  int ambientSign = -1;
  int maskSign = 1;
  IndexType requestedMinIndex{};
  IndexType requestedMaxIndex{};
  unsigned iterations = 0;
  T residual = std::numeric_limits<T>::max();
  T normalizedResidual_ = std::numeric_limits<T>::max();
  bool lastSolveConverged_ = false;
  T maxScalarVelocity_ = 0.;
  bool solved = false;
  bool nodesDirty_ = true; // true → rebuild grid/nodes on next apply()
  mutable std::string
      lastLoggedBackend_; // suppresses repeated "using X" messages
  bool useRequestedBounds = false;
  bool warmStartable_ =
      false; // true when nodes[i].concentration holds a prior solution
  std::unordered_map<std::size_t, T> pressureLookup;
  // Oxide Cauchy stress per node (row-major 3x3), fed back by OxidationModel
  // for the normal-stress reaction law.
  std::unordered_map<std::size_t, std::array<T, 9>> stressLookup;
  // Interface normal stress per node for the reaction law (NaN = none);
  // see updateInterfaceNormalStress.
  std::vector<T> interfaceNormalStress_;
  std::vector<T> lawPressure_;
  std::unordered_map<std::size_t, T> concentrationCache_;
  std::vector<Node> nodes;
  // Face-major flat BC arrays: index = fi * n + nodeId, where fi in [0, 2*D).
  // Face-major layout gives coalesced GPU reads when all warp threads access
  // the same face of consecutive nodes.
  std::vector<Boundary> faceBCTypes_;
  std::vector<T> faceBCDists_;
  // |n_axis| of the Si interface normal on REACTION faces (see
  // reactionBoundarySide); 1 elsewhere.
  std::vector<T> faceReactionNormal_;
  // The Si front's interface points are the crossings of the REACTION faces
  // (node + faceBCDists_ along the face's axis).  Everything the front does
  // is evaluated there, never per node: the normal (the level-set gradient
  // interpolated to the crossing), the reaction rate k (stress and
  // orientation laws), and the front speed k C_i / (N alpha) with the
  // crossing's own interface concentration C_i.  The oxidant sink, the
  // advection speed of every level-set point and the expansion velocity of
  // every mechanics ghost interpolate these along the front.  Per-node values
  // copied to nearby points (the earlier form) made the speed depend on where
  // the front sat in its cell, and the stress feedback grew that stair-step
  // pattern into a corrugated front (2026-10-26, B9).
  std::vector<Vec3D<T>> faceCrossingNormal_;
  std::vector<T> faceReactionRate_;
  std::vector<T> faceFrontSpeed_;
  // Point on the front (index units) a reaction node's law stress is sampled
  // at: the node moved back along its normal by its front distance.
  std::vector<Vec3D<T>> interfaceAnchor_;
  // Neighbor node IDs for every face (geometry-only, rebuilt with nodes).
  // Entry == noNode for boundary / out-of-bounds faces.
  std::vector<std::size_t> neighborIds_;

  GpuMode gpuMode_ = GpuMode::Cpu;
  GpuPreconditioner gpuPreconditioner_ = GpuPreconditioner::Jacobi;

#ifdef VIENNALS_GPU_BICGSTAB
  // Opaque pointer to device-side buffers (GpuBiCGSTABBuffers is defined only
  // in the CUDA translation unit; we hold a raw pointer here so g++ never sees
  // the CUDA internals).
  gpu::GpuBiCGSTABBuffers *gpuBufs_ = nullptr;
#endif

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

  OxidationDiffusion(SmartPointer<Domain<T, D>> passedReactionInterface,
                     SmartPointer<Domain<T, D>> passedAmbientInterface,
                     OxidationParameters passedParameters = {})
      : reactionInterface(passedReactionInterface),
        ambientInterface(passedAmbientInterface), parameters(passedParameters) {
  }

  ~OxidationDiffusion() {
#ifdef VIENNALS_GPU_BICGSTAB
    gpu::freeGpuBuffers(gpuBufs_);
    gpuBufs_ = nullptr;
#endif
  }

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<OxidationDiffusion>::New(std::forward<Args>(args)...);
  }

  /// Call after any in-place modification of the level sets (e.g. after
  /// ls::Advect) so that the next apply() rebuilds the Cartesian node grid.
  void markGeometryChanged() {
    nodesDirty_ = true;
    solved = false;
  }

  void setReactionInterface(SmartPointer<Domain<T, D>> passedInterface) {
    reactionInterface = passedInterface;
    nodesDirty_ = true;
    solved = false;
  }

  void setAmbientInterface(SmartPointer<Domain<T, D>> passedInterface) {
    ambientInterface = passedInterface;
    nodesDirty_ = true;
    solved = false;
  }

  void setMaskInterface(SmartPointer<Domain<T, D>> passedInterface,
                        int passedMaskSign = 1) {
    maskInterface = passedInterface;
    maskSign = (passedMaskSign < 0) ? -1 : 1;
    nodesDirty_ = true;
    solved = false;
  }

  void clearMaskInterface() {
    maskInterface = nullptr;
    nodesDirty_ = true;
    solved = false;
  }

  void setParameters(OxidationParameters passedParameters) {
    parameters = passedParameters;
    solved = false;
  }

  OxidationParameters getParameters() const { return parameters; }

  T getEffectiveReactionRate(const Vec3D<T> &coordinate) const {
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    return getEffectiveReactionRate(index);
  }

  void clearPressureField() {
    pressureLookup.clear();
    stressLookup.clear();
    solved = false;
  }

  /// Oxide Cauchy stress (Pa, tension positive, row-major 3x3) at a node,
  /// for the normal-stress reaction law.
  void setStress(const IndexType &index, const std::array<T, 9> &stress) {
    const auto key = detail::gridIndexHash<D>(index);
    if (std::all_of(stress.begin(), stress.end(),
                    [](T value) { return std::isfinite(value); }))
      stressLookup[key] = stress;
    else
      stressLookup.erase(key);
    solved = false;
  }

  /// The pressure and stress the last solve used at a node (reference
  /// pressure / zero stress where nothing was fed back).
  void getFeedback(const IndexType &index, T &pressure,
                   std::array<T, 9> &stress) const {
    const auto key = detail::gridIndexHash<D>(index);
    const auto p = pressureLookup.find(key);
    pressure = p != pressureLookup.end() ? p->second
                                         : T(parameters.referencePressure);
    const auto s = stressLookup.find(key);
    if (s != stressLookup.end())
      stress = s->second;
    else
      stress.fill(T(0));
  }

  using StressLookup = std::unordered_map<std::size_t, std::array<T, 9>>;

  using PressureLookup = std::unordered_map<std::size_t, T>;

  /// Pressure and stress feedback carried from one solve to the next, like
  /// the concentration cache: lsOxidation builds a new OxidationDiffusion
  /// for every coupled solve, and without the carry-over every solve began
  /// its stress coupling from zero stress (~7 coupling iterations per solve
  /// on a Kao cylinder, residual 1.0 at the first).
  const PressureLookup &getPressureCache() const { return pressureLookup; }
  const StressLookup &getStressCache() const { return stressLookup; }
  void setFeedbackCache(PressureLookup pressures, StressLookup stresses) {
    pressureLookup = std::move(pressures);
    stressLookup = std::move(stresses);
    solved = false;
  }

  /// The interface normal stress the reaction law used in the last solve
  /// (per node; NaN where the law falls back to the unstressed rate).
  const std::vector<T> &getInterfaceNormalStress() const {
    return interfaceNormalStress_;
  }

  /// The interface normal stress the reaction law would use with `lookup`:
  /// the median sample at every node with a REACTION face, smoothed along the
  /// interface (see updateInterfaceNormalStress for why).  Non-mutating, so
  /// the coupling loop can measure convergence without changing the state.
  std::vector<T> interfaceNormalStressFor(const StressLookup &lookup) const {
    const std::size_t n = nodes.size();
    constexpr T nan = std::numeric_limits<T>::quiet_NaN();
    std::vector<T> value(n, nan);
    if (faceBCTypes_.size() != 2 * D * n)
      return value;
    for (std::size_t i = 0; i < n; ++i) {
      if (nodes[i].inMask || !reactionNode(i))
        continue;
      const auto &normal = nodes[i].siNormal;
      Vec3D<T> into = normal; // into the oxide
      for (unsigned d = 0; d < D; ++d)
        into[d] *= static_cast<T>(reactionSign);
      const auto normalStress = [&](std::size_t j) {
        const auto found =
            lookup.find(detail::gridIndexHash<D>(nodes[j].index));
        if (found == lookup.end())
          return std::numeric_limits<T>::quiet_NaN();
        T v = T(0);
        for (unsigned r = 0; r < D; ++r)
          for (unsigned c = 0; c < D; ++c)
            v += normal[r] * found->second[3 * r + c] * normal[c];
        return v;
      };
      T sample;
      if (frontStressAt(i, into, frontDistance(nodes[i].index), T(0),
                        normalStress, sample) ||
          medianNormalStress(nodes[i].index, lookup, sample))
        value[i] = sample;
    }
    std::vector<T> next(n);
    for (unsigned pass = 0; pass < 6; ++pass) {
      for (std::size_t i = 0; i < n; ++i) {
        next[i] = value[i];
        if (!std::isfinite(value[i]))
          continue;
        T sum = value[i];
        unsigned count = 1;
        IndexType offset;
        offset.fill(-1);
        while (true) {
          bool centre = true;
          IndexType neighbor = nodes[i].index;
          for (unsigned d = 0; d < D; ++d) {
            neighbor[d] += offset[d];
            centre = centre && offset[d] == 0;
          }
          if (!centre && inBounds(neighbor)) {
            const std::size_t j = lookupNode(neighbor);
            if (j != noNode && std::isfinite(value[j])) {
              sum += value[j];
              ++count;
            }
          }
          unsigned d = 0;
          while (d < D && offset[d] == 1)
            offset[d++] = -1;
          if (d == D)
            break;
          ++offset[d];
        }
        next[i] = sum / static_cast<T>(count);
      }
      value.swap(next);
    }
    return value;
  }

  void setStress(const Vec3D<T> &coordinate, const std::array<T, 9> &stress) {
    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);
    setStress(index, stress);
  }

  void setPressure(const IndexType &index, T pressure) {
    const auto key = detail::gridIndexHash<D>(index);
    if (std::isfinite(pressure))
      pressureLookup[key] = pressure;
    else
      pressureLookup.erase(key);
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
    nodesDirty_ = true;
    solved = false;
  }

  void clearSolveBounds() {
    useRequestedBounds = false;
    nodesDirty_ = true;
    solved = false;
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
    if (parameters.reactionActivationVolume != T(0))
      updateInterfaceNormalStress();
    if (parameters.diffusionActivationVolume != T(0))
      updateLawPressure();
    updateFaceRates();
    solveDiffusion();

    concentrationCache_.clear();
    for (const auto &node : nodes)
      concentrationCache_[detail::gridIndexHash<D>(node.index)] =
          node.concentration;

    updateFaceSpeeds();
    maxScalarVelocity_ = 0.;
    for (const T speed : faceFrontSpeed_)
      maxScalarVelocity_ = std::max(maxScalarVelocity_, speed);

    solved = true;

    // TEST INSTRUMENTATION (env-gated): dump one vertical column of node
    // concentrations — including meshed-mask nodes — so the steady profile
    // through the film is observable.  Locates where a resistance error sits:
    // a linear drop across the film means D_mask itself is off; a jump at one
    // face means that face's coefficient or BC is wrong.
    if (std::getenv("OX_DUMP_CONC") != nullptr && !nodes.empty()) {
      const auto column = nodes[nodes.size() / 2].index[0];
      if (FILE *dump = std::fopen("/tmp/ox_conc_dump.txt", "w")) {
        for (const auto &node : nodes)
          if (node.index[0] == column)
            std::fprintf(dump, "%.6f %d %.8f\n",
                         static_cast<double>(node.index[D - 1]) * gridDelta,
                         node.inMask ? 1 : 0,
                         static_cast<double>(node.concentration));
        std::fclose(dump);
      }
    }
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

    // Meshed mask: mask nodes carry oxidant concentration but DO NOT react.
    // A reaction speed queried at or nearest to a mask point is zero — never
    // k*C from the film's near-C_eq concentration.  Without this, velocity
    // fallbacks that sample the reaction speed at mask coordinates invent a
    // full oxidation rate inside the mask body.
    {
      std::size_t nodeId = lookupNode(index);
      if (nodeId == noNode)
        nodeId = findNearbyNode(index);
      if (nodeId != noNode && nodes[nodeId].inMask)
        return 0.;
    }

    // The front speed at the point of the front closest to this one,
    // x - phi grad(phi) / |grad(phi)|.
    {
      ConstSparseIterator reactionIt(reactionInterface->getDomain());
      const T phi = valueAt(reactionIt, index);
      if (std::isfinite(phi) && std::abs(phi) <= T(2)) {
        const auto g = computeSiNormal(index, reactionIt);
        Vec3D<T> p{T(0), T(0), T(0)};
        for (unsigned d = 0; d < D; ++d)
          p[d] = (static_cast<T>(index[d]) - phi * g[d]) * gridDelta;
        const auto front = frontVelocityAt(p);
        if (front.found)
          return parameters.velocitySign * front.speed;
      }
    }

    const auto boundarySample = reactionBoundarySample(index);
    const IndexType rateIndex =
        boundarySample.found ? boundarySample.nodeIndex : index;
    const T concentration = boundarySample.found
                                ? boundarySample.concentration
                                : getReactionBoundaryConcentration(index);
    return parameters.velocitySign * getEffectiveReactionRate(rateIndex) *
           concentration /
           (parameters.oxidantMoleculeDensity *
            parameters.expansionCoefficient);
  }

  /// Si front speed (k C_i / (N alpha), >= 0) and normal (level-set gradient
  /// direction) at a point on the front: the area-weighted average over the
  /// face crossings within two cells (Gaussian weights, 0.75 cells wide).
  /// Not found before the first solve or away from the front.
  struct FrontVelocity {
    bool found = false;
    T speed = 0.;
    Vec3D<T> normal{0., 0., 0.};
  };

  FrontVelocity frontVelocityAt(const Vec3D<T> &coordinate) const {
    FrontVelocity result;
    const std::size_t n = nodes.size();
    if (faceFrontSpeed_.size() != 2 * D * n)
      return result;
    Vec3D<T> p{T(0), T(0), T(0)};
    for (unsigned d = 0; d < D; ++d)
      p[d] = coordinate[d] / gridDelta;
    constexpr T width = T(0.75);
    T weightSum = T(0), speed = T(0);
    Vec3D<T> normal{T(0), T(0), T(0)};
    forEachNodeNear(p, 2, [&](std::size_t j) {
      for (unsigned fi = 0; fi < 2 * D; ++fi) {
        if (faceBCTypes_[fi * n + j] != Boundary::REACTION)
          continue;
        const auto c = crossingPoint(j, fi);
        T r2 = T(0);
        for (unsigned d = 0; d < D; ++d)
          r2 += (c[d] - p[d]) * (c[d] - p[d]);
        // weighted by the crossing's share of the front area, |n_axis|
        // (the factor its oxidant sink carries): an area average, so the
        // extra crossings at stair steps do not count twice, and grazing
        // arms -- the least accurate interface concentration -- hardly count
        const T w = std::exp(-r2 / (T(2) * width * width)) *
                    faceReactionNormal_[fi * n + j];
        weightSum += w;
        speed += w * faceFrontSpeed_[fi * n + j];
        for (unsigned d = 0; d < D; ++d)
          normal[d] += w * faceCrossingNormal_[fi * n + j][d];
      }
    });
    if (weightSum <= T(1e-12))
      return result;
    T len = T(0);
    for (unsigned d = 0; d < D; ++d)
      len += normal[d] * normal[d];
    len = std::sqrt(len);
    if (len <= std::numeric_limits<T>::epsilon())
      return result;
    result.found = true;
    result.speed = speed / weightSum;
    for (unsigned d = 0; d < D; ++d)
      result.normal[d] = normal[d] / len;
    return result;
  }

  T getDissipationAlpha(int /*direction*/, int material,
                        const Vec3D<T> & /*centralDifferences*/) final {
    if (parameters.material >= 0 && material != parameters.material)
      return 0.;

    T bulkVel =
        std::abs(parameters.velocitySign) * parameters.reactionRate *
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
    const std::size_t nodeId = lookupNode(index);
    if (nodeId == noNode) {
      const auto nearby = findNearbyNode(index);
      if (nearby == noNode)
        return 0.;
      return nodes[nearby].concentration;
    }
    return nodes[nodeId].concentration;
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
  T getNormalizedResidual() const { return normalizedResidual_; }
  bool lastSolveConverged() const { return lastSolveConverged_; }
  std::size_t getNumberOfSolutionNodes() const { return nodes.size(); }
  bool hasFiniteConcentrationField() const {
    for (const auto &node : nodes)
      if (!std::isfinite(node.concentration))
        return false;
    return true;
  }

  const std::unordered_map<std::size_t, T> &getConcentrationCache() const {
    return concentrationCache_;
  }

  void setConcentrationCache(std::unordered_map<std::size_t, T> cache) {
    concentrationCache_ = std::move(cache);
  }

  /// Set the GPU solver selection mode.  See GpuMode for the two options.
  /// On CPU-only builds (VIENNALS_GPU_BICGSTAB not defined) this is a no-op.
  void setGpuMode(GpuMode mode) { gpuMode_ = mode; }
  /// Set the GPU BiCGSTAB preconditioner (Jacobi by default).
  void setGpuPreconditioner(GpuPreconditioner preconditioner) {
    gpuPreconditioner_ = preconditioner;
  }

  /// Mark the current solution as valid without re-solving. Call this before
  /// any parallel advection (lsAdvect) that uses this field as a velocity
  /// source to prevent concurrent apply() calls inside getScalarVelocity().
  void markSolved() { solved = true; }

  /// Write per-node concentration into ambientInterface->getPointData() so
  /// that lsInterior + lsAdvect can carry it across timestep boundaries.
  /// Safe to call regardless of the solved flag; uses the most recent nodes.
  void writeConcentrationToLevelSet() {
    if (nodes.empty() || ambientInterface == nullptr)
      return;

    std::vector<T> concentrations;
    ConstSparseIterator it(ambientInterface->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const IndexType idx = it.getStartIndices();
      const std::size_t nodeId = lookupNode(idx);
      // Non-solve-grid narrow-band points (mask interior, gas-phase boundary)
      // get concentration 0.  Using equilibriumConcentration here contaminates
      // the output: mask-interior points are impermeable (C≈0), and newly
      // appeared bird's-beak points inherit C=1 spuriously after advection.
      // The in-memory concentrationCache_ is the warm-start source for the
      // next substep; the level-set value is only the fallback when that cache
      // is cold, and C=0 converges just as quickly as C=1 from a cold start.
      const T value = nodeId != noNode ? nodes[nodeId].concentration : T(0);
      concentrations.push_back(std::isfinite(value) ? value : T(0));
    }
    ambientInterface->getPointData().insertReplaceScalarData(
        std::move(concentrations), "OxConcentration");
  }

  /// Write per-node pressure into ambientInterface->getPointData() so that
  /// it survives advection and can warm-start the coupling loop next step.
  void writePressureToLevelSet() {
    if (pressureLookup.empty() || ambientInterface == nullptr)
      return;

    std::vector<T> pressures;
    ConstSparseIterator it(ambientInterface->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const IndexType idx = it.getStartIndices();
      const auto pIt = pressureLookup.find(detail::gridIndexHash<D>(idx));
      const T value = pIt != pressureLookup.end() ? pIt->second : T(0);
      pressures.push_back(std::isfinite(value) ? value : T(0));
    }
    ambientInterface->getPointData().insertReplaceScalarData(
        std::move(pressures), "OxPressure");
  }

  /// Convenience wrapper: persist both concentration and pressure in one call.
  void writePersistentFields() {
    writeConcentrationToLevelSet();
    writePressureToLevelSet();
  }

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
    return std::abs(
        getEffectiveReactionRate(sample.nodeIndex) * sample.concentration /
        (parameters.oxidantMoleculeDensity * parameters.expansionCoefficient));
  }

private:
  bool initialiseGrid() {
    return initializeGridFromInterfaces(
        reactionInterface, ambientInterface, maskInterface, useRequestedBounds,
        requestedMinIndex, requestedMaxIndex, parameters.maxGridPoints,
        "OxidationDiffusion");
  }

  void buildNodes() {
    nodes.clear();
    initNodeLookup();

    // On the first substep after an outer-step boundary the in-memory cache is
    // empty. Fall back to the concentration stored in the level set's pointData
    // (written by writeConcentrationToLevelSet() before the previous advection
    // and remapped by lsAdvect + lsInterior).
    warmStartable_ = !concentrationCache_.empty();
    if (!warmStartable_) {
      // Single HRLE pass restores both concentration and pressure from the
      // pointData written by writePersistentFields() before the last advection.
      const int cIdx = ambientInterface->getPointData().getScalarDataIndex(
          "OxConcentration");
      const int pIdx =
          ambientInterface->getPointData().getScalarDataIndex("OxPressure");
      const auto *cd =
          (cIdx != -1) ? ambientInterface->getPointData().getScalarData(cIdx)
                       : nullptr;
      const auto *pd =
          (pIdx != -1) ? ambientInterface->getPointData().getScalarData(pIdx)
                       : nullptr;
      // The normal-stress reaction law starts from the previous step's oxide
      // stress: the deviatoric end-of-step stress the deformation solver
      // carries with the flow (OxStressR0..2) minus the pressure.
      const auto &pointData = ambientInterface->getPointData();
      const std::array<int, 3> sIdx = {
          pointData.getVectorDataIndex("OxStressR0"),
          pointData.getVectorDataIndex("OxStressR1"),
          pointData.getVectorDataIndex("OxStressR2")};
      const bool stressWarmStart = parameters.reactionActivationVolume != 0. &&
                                   pd != nullptr && sIdx[0] != -1 &&
                                   sIdx[1] != -1 && sIdx[2] != -1;
      if (cd != nullptr || pd != nullptr) {
        ConstSparseIterator it(ambientInterface->getDomain());
        for (; !it.isFinished(); ++it) {
          if (!it.isDefined())
            continue;
          const auto ptId = it.getPointId();
          const std::size_t key =
              detail::gridIndexHash<D>(it.getStartIndices());
          if (cd != nullptr && ptId < static_cast<decltype(ptId)>(cd->size()) &&
              std::isfinite((*cd)[ptId]))
            concentrationCache_[key] = (*cd)[ptId];
          if (pd != nullptr && ptId < static_cast<decltype(ptId)>(pd->size()) &&
              std::isfinite((*pd)[ptId]))
            pressureLookup[key] = (*pd)[ptId];
          if (stressWarmStart &&
              ptId < static_cast<decltype(ptId)>(pd->size())) {
            std::array<T, 9> stress{};
            bool valid = std::isfinite((*pd)[ptId]);
            for (unsigned r = 0; r < 3 && valid; ++r) {
              const auto *rows = pointData.getVectorData(sIdx[r]);
              valid = rows != nullptr &&
                      ptId < static_cast<decltype(ptId)>(rows->size());
              for (unsigned c = 0; c < 3 && valid; ++c)
                stress[3 * r + c] = T((*rows)[ptId][c]);
            }
            if (valid) {
              for (unsigned d = 0; d < 3; ++d)
                stress[4 * d] -= T((*pd)[ptId]);
              stressLookup[key] = stress;
            }
          }
        }
        warmStartable_ = !concentrationCache_.empty();
      }
    }

    ConstSparseIterator reactionIt(reactionInterface->getDomain());
    ConstSparseIterator ambientIt(ambientInterface->getDomain());
    auto maskIt = makeMaskIterator();

    // Meshed mask diffusion: with gamma > 0 the mask body joins the diffusion
    // domain as nodes of their own (D/gamma), so the oxidant concentration is
    // solved THROUGH the film instead of being lumped into a boundary
    // coefficient.
    const bool meshMask =
        maskInterface != nullptr && parameters.maskPermeabilityRatio > 0.;

    IndexType index = minIndex;
    while (true) {
      const T reactionPhi = valueAt(reactionIt, index);
      const T ambientPhi = valueAt(ambientIt, index);
      const bool insideMask = isInsideMask(maskIt, index);
      constexpr T eps = T(1e-9);
      const bool oxideNode =
          isInsideOxide(reactionPhi, ambientPhi) && !insideMask;
      const bool maskNode =
          meshMask && insideMask && reactionSign * reactionPhi >= -eps;
      if (oxideNode || maskNode) {
        const std::size_t id = nodes.size();
        nodeLookupFlat[linearIndex(index)] = id;
        T seedConc = parameters.equilibriumConcentration;
        auto cacheIt =
            concentrationCache_.find(detail::gridIndexHash<D>(index));
        if (cacheIt != concentrationCache_.end())
          seedConc = cacheIt->second;
        if (!std::isfinite(seedConc))
          seedConc = parameters.equilibriumConcentration;
        Node newNode{index, seedConc};
        newNode.inMask = maskNode;
        // The orientation factor and the normal-stress reaction law both
        // need the interface normal.
        if (!maskNode && (parameters.reactionRateRatio111 != T(1) ||
                          parameters.reactionActivationVolume != T(0)))
          newNode.siNormal = computeSiNormal(index, reactionIt);
        nodes.push_back(newNode);
      }

      if (!increment(index))
        break;
    }

    // Precompute per-face boundary intersections into flat face-major arrays.
    // Level-set positions are fixed during the inner solve; one pass per
    // apply().
    const std::size_t n = nodes.size();
    faceBCTypes_.assign(2 * D * n, Boundary::NONE);
    faceBCDists_.assign(2 * D * n, T(1));
    faceReactionNormal_.assign(2 * D * n, T(1));
    faceCrossingNormal_.assign(2 * D * n, Vec3D<T>{T(0), T(0), T(0)});
    ConstSparseIterator faceReactionIt(reactionInterface->getDomain());
    ConstSparseIterator faceAmbientIt(ambientInterface->getDomain());
    auto faceMaskIt = makeMaskIterator();
    for (std::size_t id = 0; id < n; ++id) {
      const auto &node = nodes[id];
      for (unsigned dir = 0; dir < D; ++dir) {
        for (int off : {-1, 1}) {
          const unsigned fi = dir * 2u + (off == 1 ? 1u : 0u);
          IndexType nb = node.index;
          nb[dir] += off;
          if (!inBounds(nb) || lookupNode(nb) != noNode)
            continue; // NONE/1.0 already set by assign()
          auto bc = classifyBoundary(faceReactionIt, faceAmbientIt,
                                     faceMaskIt, node.index, nb);
          if (node.inMask) {
            // Faces of a meshed mask node.  Crossing the mask level set here
            // is the film's OUTER surface — the gas inlet — so it takes the
            // ambient Robin BC.  Crossing the Si level set is mask resting on
            // silicon: the mask does not react with Si, so zero flux.
            if (bc.first == Boundary::MASK)
              bc.first = Boundary::AMBIENT;
            else if (bc.first == Boundary::REACTION)
              bc.first = Boundary::NONE;
          }
          faceBCTypes_[fi * n + id] = bc.first;
          faceBCDists_[fi * n + id] = bc.second;
          if (bc.first == Boundary::REACTION) {
            // the normal at the crossing: the node's and the neighbour's
            // level-set gradients, interpolated to the crossing point
            const auto a = computeSiNormal(node.index, faceReactionIt);
            const auto b = computeSiNormal(nb, faceReactionIt);
            const T t = std::clamp(bc.second / gridDelta, T(0), T(1));
            Vec3D<T> m{T(0), T(0), T(0)};
            T len = T(0);
            for (unsigned e = 0; e < D; ++e) {
              m[e] = (T(1) - t) * a[e] + t * b[e];
              len += m[e] * m[e];
            }
            len = std::sqrt(len);
            if (len > std::numeric_limits<T>::epsilon())
              for (unsigned e = 0; e < D; ++e)
                m[e] /= len;
            else
              m = a;
            faceCrossingNormal_[fi * n + id] = m;
            faceReactionNormal_[fi * n + id] = std::abs(m[dir]);
          }
        }
      }
    }

    interfaceAnchor_.assign(n, Vec3D<T>{T(0), T(0), T(0)});
    if (parameters.reactionActivationVolume != T(0))
      for (std::size_t id = 0; id < n; ++id) {
        if (nodes[id].inMask || !reactionNode(id))
          continue;
        const T phi = frontDistance(nodes[id].index);
        for (unsigned d = 0; d < D; ++d)
          interfaceAnchor_[id][d] =
              static_cast<T>(nodes[id].index[d]) -
              static_cast<T>(reactionSign) * nodes[id].siNormal[d] * phi;
      }

    // Precompute per-face neighbor node IDs (geometry-only).
    // Used by computeFaceCoeffs() and uploaded once to GPU.
    neighborIds_.assign(2 * D * n, noNode);
    for (std::size_t id = 0; id < n; ++id) {
      for (unsigned dir = 0; dir < D; ++dir) {
        for (int off : {-1, 1}) {
          const unsigned fi = dir * 2u + (off == 1 ? 1u : 0u);
          IndexType nb = nodes[id].index;
          nb[dir] += off;
          if (inBounds(nb))
            neighborIds_[fi * n + id] = lookupNode(nb);
        }
      }
    }

    bool loggedBackend = false;
#ifdef VIENNALS_GPU_BICGSTAB
    // Free any stale allocation from a previous buildNodes() call.
    gpu::freeGpuBuffers(gpuBufs_);
    gpuBufs_ = nullptr;

    const bool tryGpu = (gpuMode_ == GpuMode::Gpu || gpuMode_ == GpuMode::Auto);
    if (tryGpu) {
      const bool useIlu0 = gpuPreconditioner_ == GpuPreconditioner::ILU0;
      gpuBufs_ = gpu::allocGpuBuffers(static_cast<uint32_t>(n), 2 * D, useIlu0);

      // Each setup step is chained with else-if: once one fails the handle is
      // released and set to null, and every later step must be skipped rather
      // than called with a null handle. (In GpuMode::Gpu reportGpuUnavailable
      // throws, but in GpuMode::Auto it returns and execution continues here.)
      if (!gpuBufs_) {
        reportGpuUnavailable("OxidationDiffusion: GPU mode was selected, but "
                             "CUDA buffers could not be allocated or the CUDA "
                             "context could not be initialized.");
      } else {
        // Convert std::size_t neighbor IDs to uint32_t for the device
        const std::size_t nf = 2u * D * n;
        std::vector<uint32_t> nb32(nf);
        for (std::size_t k = 0; k < nf; ++k)
          nb32[k] = (neighborIds_[k] == noNode)
                        ? gpu::kNoNode
                        : static_cast<uint32_t>(neighborIds_[k]);
        if (!gpu::gpuUploadNeighborIds(gpuBufs_, nb32.data(), nf)) {
          gpu::freeGpuBuffers(gpuBufs_);
          gpuBufs_ = nullptr;
          reportGpuUnavailable("OxidationDiffusion: GPU mode was selected, but "
                               "uploading GPU neighbor IDs failed.");
        } else if (useIlu0 &&
                   !gpu::gpuSetupCSR(gpuBufs_, nb32.data(),
                                     static_cast<uint32_t>(n), 2 * D)) {
          gpu::freeGpuBuffers(gpuBufs_);
          gpuBufs_ = nullptr;
          reportGpuUnavailable("OxidationDiffusion: GPU mode was selected, but "
                               "CUSPARSE setup for the GPU BiCGSTAB solver "
                               "failed.");
        } else {
          // Only claim the GPU backend once every setup step has succeeded.
          logDiffusionBackend("GPU BiCGSTAB",
                              "preconditioner=" +
                                  std::string(useIlu0 ? "ILU0" : "Jacobi"));
          loggedBackend = true;
        }
      }
    }
#endif
    if (!loggedBackend) {
#ifdef VIENNALS_GPU_BICGSTAB
      logDiffusionBackend("CPU AMG-BiCGSTAB",
                          gpuMode_ == GpuMode::Cpu
                              ? "GPU mode not selected"
                              : "GPU requested but unavailable");
#else
      logDiffusionBackend("CPU AMG-BiCGSTAB",
                          "ViennaLS was built without GPU BiCGSTAB support");
#endif
    }
  }

  // Per-face coupling coefficients: the off-diagonal entries of the
  // diffusion matrix (CPU solve and GPU SpMV).
  //
  //   faceCoeffs[fi * n + id] = 2*D_eff / (dist_fi * distSum_axis)
  //
  // Only interior faces (neighbor != noNode) get a nonzero entry;
  // boundary and out-of-bounds faces contribute to diag/b but not to
  // the off-diagonal coupling stored here.  Must be called after diag/b
  // are computed (or at least after D_eff is known), because D_eff
  // depends on pressure when diffusionActivationVolume != 0.
  void computeFaceCoeffs(std::vector<T> &faceCoeffs) const {
    const std::size_t n = nodes.size();
    faceCoeffs.assign(2 * D * n, T(0));
    const T eps = std::numeric_limits<T>::epsilon();
    const std::vector<T> zeros(n, T(0));

    for (std::size_t id = 0; id < n; ++id) {
      const T D_eff = nodeDiffusionCoefficient(id);
      for (unsigned dir = 0; dir < D; ++dir) {
        const unsigned fiNeg = dir * 2u;
        const unsigned fiPos = dir * 2u + 1u;

        // Use the exact side construction used by the CPU stencil so the GPU
        // SpMV cannot drift from computeStencilAt() at cut-cell boundaries.
        const auto neg = makeStencilSide(id, zeros, dir, -1, D_eff);
        const auto pos = makeStencilSide(id, zeros, dir, 1, D_eff);
        const T distNeg = neg.distance;
        const T distPos = pos.distance;
        const T distSum = distNeg + distPos;
        if (distSum <= eps)
          continue;

        // Per-face diffusion, matching computeStencilAt: the harmonic mean
        // of the two nodes' D.  A mirror side couples to the opposite
        // neighbour: its coefficient goes there.
        const unsigned fiNegEff =
            stencilOffset(nodes[id].index, dir, -1) == 1 ? fiPos : fiNeg;
        const unsigned fiPosEff =
            stencilOffset(nodes[id].index, dir, 1) == -1 ? fiNeg : fiPos;
        if (neighborIds_[fiNegEff * n + id] != noNode && distNeg > eps)
          faceCoeffs[fiNegEff * n + id] +=
              T(2) * faceDiffusionCoefficient(id, dir, -1, D_eff) /
              (distNeg * distSum);
        if (neighborIds_[fiPosEff * n + id] != noNode && distPos > eps)
          faceCoeffs[fiPosEff * n + id] +=
              T(2) * faceDiffusionCoefficient(id, dir, 1, D_eff) /
              (distPos * distSum);
      }
    }
  }

  // Evaluates the stencil at one node: fills diag = A[i,i] and
  // rhs = sum_j(A_off[i,j] * x[j]) + bc_constants[i].
  // Called with x = zeros to precompute the geometry-fixed diagonal and b.
  template <class SolverT>
  void computeStencilAt(std::size_t nodeId, const std::vector<SolverT> &x,
                        T &diag, T &rhs) const {
    diag = T(0);
    rhs = T(0);
    const T D_eff = nodeDiffusionCoefficient(nodeId);
    for (unsigned direction = 0; direction < D; ++direction) {
      const auto neg = makeStencilSide(nodeId, x, direction, -1, D_eff);
      const auto pos = makeStencilSide(nodeId, x, direction, 1, D_eff);
      // Each side carries its own face diffusion, the harmonic mean of the
      // two nodes' D.
      const T distSum = neg.distance + pos.distance;
      if (distSum <= std::numeric_limits<T>::epsilon())
        continue;
      addSideContribution(rhs, diag, neg, distSum,
                          faceDiffusionCoefficient(nodeId, direction, -1,
                                                   D_eff));
      addSideContribution(rhs, diag, pos, distSum,
                          faceDiffusionCoefficient(nodeId, direction, 1,
                                                   D_eff));
    }
  }

  void solveDiffusion() {
    iterations = 0;
    residual = 0.;
    normalizedResidual_ = 0.;
    lastSolveConverged_ = false;
    if (nodes.empty()) {
      lastSolveConverged_ = true;
      return;
    }

    const std::size_t n = nodes.size();
    const T eps = std::numeric_limits<T>::epsilon();

    // Geometry-fixed diagonal and BC source vector (kept in T for full
    // precision).
    Timer<> tDiag;
    tDiag.start();
    std::vector<T> diag(n), b(n);
    {
      const std::vector<T> zeros(n, T(0));
#pragma omp parallel for schedule(static)
      for (std::size_t i = 0; i < n; ++i)
        computeStencilAt(i, zeros, diag[i], b[i]);
    }
    tDiag.finish();

    bool finiteSystem = true;
    for (std::size_t i = 0; i < n; ++i)
      finiteSystem =
          finiteSystem && std::isfinite(diag[i]) && std::isfinite(b[i]);
    if (!finiteSystem) {
      residual = std::numeric_limits<T>::infinity();
      normalizedResidual_ = residual;
      Logger::getInstance()
          .addWarning("solveDiffusion: assembled non-finite matrix/RHS; "
                      "rejecting this coupled trial.")
          .print();
      return;
    }

    // Scale of the residual test: a typical row (the median diagonal) times
    // the concentration scale (max |b_i / diag_i|).  Not max|b|: a node on
    // the ambient surface with the ambient on both sides along one axis has
    // a row of 1 / distance (2e10 for a node exactly on it), and the solve
    // then accepted its initial guess, C = 0 in the oxide (a dry Kao pillar
    // whose surface passed through grid points did not oxidise).
    T b_norm = T(0), bMax = T(0);
    {
      std::vector<T> sorted(diag);
      const auto median = sorted.begin() + n / 2;
      std::nth_element(sorted.begin(), median, sorted.end());
      T concentrationScale = T(0);
      for (std::size_t i = 0; i < n; ++i) {
        bMax = std::max(bMax, std::abs(b[i]));
        if (diag[i] > eps)
          concentrationScale =
              std::max(concentrationScale, std::abs(b[i] / diag[i]));
      }
      b_norm = *median * concentrationScale;
    }
    if (b_norm < T(1e-100))
      b_norm = T(1);

    // Initial guess: warm-start or diagonal-preconditioned b.  Sanitize the
    // warm start so a previous failed solve cannot seed NaNs into the next one.
    auto fallbackInitialGuess = [&](std::size_t i) -> T {
      if (std::isfinite(diag[i]) && std::isfinite(b[i]) && diag[i] > eps) {
        const T guess = b[i] / diag[i];
        if (std::isfinite(guess))
          return guess;
      }
      return T(0);
    };

    std::vector<T> x(n);
    for (std::size_t i = 0; i < n; ++i) {
      x[i] = warmStartable_ ? nodes[i].concentration : fallbackInitialGuess(i);
      if (!std::isfinite(x[i]))
        x[i] = fallbackInitialGuess(i);
    }
    const auto initialGuess = x;

#ifdef VIENNALS_GPU_BICGSTAB
    // ── GPU BiCGSTAB path ──────────────────────────────────────────────
    // Engaged when GpuMode::Gpu is set and buildNodes() allocated GPU buffers.
    // Uploads diag, b, and face coefficients fresh each call (they are
    // pressure-dependent via D_eff), then runs the full BiCGSTAB loop on
    // the device.  Only dot-product scalars and the max-abs residual are
    // downloaded; work vectors stay on the GPU between iterations.
    if (gpu::gpuIsValid(gpuBufs_)) {
      Timer<> tPrep, tUpload, tSolve;
      tPrep.start();
      std::vector<double> diagGpu(n), bGpu(n);
      for (std::size_t i = 0; i < n; ++i) {
        diagGpu[i] = static_cast<double>(diag[i]);
        bGpu[i] = static_cast<double>(b[i]);
      }

      // Face coupling coefficients (pressure-dependent; recomputed each call)
      std::vector<T> faceCoeffs;
      computeFaceCoeffs(faceCoeffs);
      std::vector<double> coeffGpu(faceCoeffs.size());
      for (std::size_t k = 0; k < faceCoeffs.size(); ++k)
        coeffGpu[k] = static_cast<double>(faceCoeffs[k]);

      std::vector<double> xGpu(n);
      for (std::size_t i = 0; i < n; ++i)
        xGpu[i] = static_cast<double>(x[i]);
      tPrep.finish();

      tUpload.start();
      const bool gpuUploadOk = gpu::gpuUploadSolverArrays(
          gpuBufs_, diagGpu.data(), bGpu.data(), coeffGpu.data(),
          static_cast<uint32_t>(n), faceCoeffs.size());
      tUpload.finish();
      if (!gpuUploadOk) {
        if (Logger::hasDebug()) {
          const std::string tag =
              "diffusion n=" + std::to_string(n) + " [GPU upload failed]";
          Logger::getInstance()
              .addTiming(tag + " diag/b precompute", tDiag)
              .addTiming(tag + " GPU prep+faceCoeffs", tPrep)
              .addTiming(tag + " GPU upload", tUpload)
              .print();
        }
        VIENNACORE_LOG_ERROR("OxidationDiffusion: GPU mode was selected, but "
                             "uploading GPU solver arrays or factorizing ILU "
                             "failed." +
                             gpuErrorDetail());
      }

      // Solve on GPU; xGpu holds the initial guess and receives the solution.
      tSolve.start();
      // the GPU solver tests max|r| < tolerance x max|b| itself
      const T gpuTolerance = bMax > T(1e-100)
                                 ? parameters.tolerance * b_norm / bMax
                                 : T(parameters.tolerance);
      const bool gpuConverged = gpu::gpuSolveBiCGSTAB(
          gpuBufs_, xGpu.data(), static_cast<double>(eps),
          parameters.maxIterations, static_cast<double>(gpuTolerance),
          iterations, residual);
      tSolve.finish();

      const bool finiteGpuSolution =
          gpuConverged &&
          std::all_of(xGpu.begin(), xGpu.end(),
                      [](double value) { return std::isfinite(value); });
      const unsigned gpuIterations = iterations;
      const T gpuResidual = residual;

      if (finiteGpuSolution) {
        // Write solution back to nodes only after convergence and finite
        // checks.
        for (std::size_t i = 0; i < n; ++i)
          nodes[i].concentration = static_cast<T>(xGpu[i]);
        normalizedResidual_ = gpuResidual / b_norm;
        lastSolveConverged_ = std::isfinite(normalizedResidual_) &&
                              normalizedResidual_ <= parameters.tolerance;

        if (Logger::hasTiming()) {
          const std::string tag = "diffusion n=" + std::to_string(n) +
                                  " iters=" + std::to_string(iterations) +
                                  " res=" + std::to_string(residual) + " [GPU]";
          Logger::getInstance()
              .addTiming(tag + " GPU BiCGSTAB", tSolve)
              .print();
        }
        if (Logger::hasDebug()) {
          const std::string tag = "diffusion n=" + std::to_string(n) + " [GPU]";
          Logger::getInstance()
              .addTiming(tag + " diag/b precompute", tDiag)
              .addTiming(tag + " GPU prep+faceCoeffs", tPrep)
              .addTiming(tag + " GPU upload", tUpload)
              .print();
        }
        return;
      }

      if (Logger::hasDebug()) {
        const std::string tag = "diffusion n=" + std::to_string(n) +
                                " iters=" + std::to_string(gpuIterations) +
                                " res=" + std::to_string(gpuResidual) +
                                " [GPU failed]";
        Logger::getInstance()
            .addTiming(tag + " diag/b precompute", tDiag)
            .addTiming(tag + " GPU prep+faceCoeffs", tPrep)
            .addTiming(tag + " GPU upload", tUpload)
            .addTiming(tag + " GPU BiCGSTAB", tSolve)
            .print();
      }
      // Under GpuMode::Auto this warns and execution continues into the CPU
      // solve below, which recomputes the solution from the initial guess.
      reportGpuUnavailable(
          "OxidationDiffusion: GPU mode was selected, but GPU BiCGSTAB "
          "failed, did not converge, or produced non-finite concentrations "
          "(iters=" +
          std::to_string(gpuIterations) +
          ", residual=" + std::to_string(gpuResidual) + ").");
    }
#else
    if (gpuMode_ == GpuMode::Gpu) {
      VIENNACORE_LOG_ERROR("OxidationDiffusion: explicit GPU mode was "
                           "requested, but ViennaLS was built without "
                           "VIENNALS_GPU_BICGSTAB.");
    } else if (gpuMode_ == GpuMode::Auto) {
      VIENNACORE_LOG_WARNING("OxidationDiffusion: GPU mode Auto was requested, "
                             "but ViennaLS was built without "
                             "VIENNALS_GPU_BICGSTAB. Using the CPU solver.");
    }
#endif

    // CPU: BiCGSTAB preconditioned by smoothed-aggregation AMG with an ILU0
    // smoother (AMGCL) on the matrix the GPU path uploads (diag plus
    // computeFaceCoeffs).  A Jacobi-preconditioned BiCGSTAB (the former CPU
    // solver, float work vectors) stagnated for thousands of iterations and
    // broke down (omega -> 0) on LOCOS bird's beaks, where the pressure-
    // dependent D spans five decades and the cut-cell rows make the matrix
    // strongly non-symmetric; the breakdown repeated on every CFL retry (the
    // quasi-steady system has no dt in it) until the run aborted.  On those
    // systems this solver takes 3-4 iterations.
    Timer<> tCpuSolve;
    tCpuSolve.start();
    std::vector<T> faceCoeffs;
    computeFaceCoeffs(faceCoeffs);
    std::vector<std::ptrdiff_t> rowPtr(n + 1, 0), cols;
    std::vector<T> values;
    cols.reserve((2 * D + 1) * n);
    values.reserve((2 * D + 1) * n);
    std::vector<std::pair<std::ptrdiff_t, T>> row;
    for (std::size_t i = 0; i < n; ++i) {
      row.assign(1, {static_cast<std::ptrdiff_t>(i), diag[i]});
      for (unsigned f = 0; f < 2 * D; ++f) {
        const std::size_t j = neighborIds_[f * n + i];
        if (j != noNode && faceCoeffs[f * n + i] != T(0))
          row.emplace_back(static_cast<std::ptrdiff_t>(j),
                           -faceCoeffs[f * n + i]);
      }
      std::sort(row.begin(), row.end());
      for (const auto &entry : row) {
        cols.push_back(entry.first);
        values.push_back(entry.second);
      }
      rowPtr[i + 1] = static_cast<std::ptrdiff_t>(cols.size());
    }

    typedef amgcl::backend::builtin<T> Backend;
    typedef amgcl::make_solver<
        amgcl::amg<Backend, amgcl::coarsening::smoothed_aggregation,
                   amgcl::relaxation::ilu0>,
        amgcl::solver::bicgstab<Backend>>
        Solver;
    typename Solver::params prm;
    prm.solver.maxiter = parameters.maxIterations;
    // ||r||_2 <= tolerance b_norm bounds the max-norm residual tested below
    prm.solver.tol = T(0);
    prm.solver.abstol = parameters.tolerance * b_norm;
    const auto rows = static_cast<std::ptrdiff_t>(n);
    auto A = std::tie(rows, rowPtr, cols, values);
    bool solverFailed = false;
    try {
      Solver solve(A, prm);
      std::size_t count = 0;
      T estimate = T(0);
      std::tie(count, estimate) = solve(A, b, x);
      iterations = static_cast<unsigned>(count);
    } catch (const std::exception &e) {
      solverFailed = true;
      VIENNACORE_LOG_WARNING(
          std::string("solveDiffusion: AMGCL failed: ") + e.what());
    }

    // the max-norm residual of the result
    residual = T(0);
    bool finiteCpuSolution = !solverFailed;
    for (std::size_t i = 0; i < n && finiteCpuSolution; ++i) {
      T ax = T(0);
      for (auto k = rowPtr[i]; k < rowPtr[i + 1]; ++k)
        ax += values[k] * x[cols[k]];
      finiteCpuSolution = std::isfinite(x[i]) && std::isfinite(ax);
      residual = std::max(residual, std::abs(b[i] - ax));
    }
    if (finiteCpuSolution) {
      for (std::size_t i = 0; i < n; ++i)
        nodes[i].concentration = x[i];
    } else {
      for (std::size_t i = 0; i < n; ++i)
        nodes[i].concentration = initialGuess[i];
      residual = std::numeric_limits<T>::infinity();
      Logger::getInstance()
          .addWarning("solveDiffusion: the CPU solve failed or produced "
                      "non-finite concentrations; keeping the sanitized "
                      "initial guess.")
          .print();
    }
    normalizedResidual_ = residual / b_norm;
    lastSolveConverged_ = finiteCpuSolution &&
                          std::isfinite(normalizedResidual_) &&
                          normalizedResidual_ <= parameters.tolerance;

    tCpuSolve.finish();
    if (Logger::hasTiming()) {
      const std::string path = " [CPU]";
      const std::string tag = "diffusion n=" + std::to_string(n) +
                              " iters=" + std::to_string(iterations) +
                              " res=" + std::to_string(residual) + path;
      Logger::getInstance()
          .addTiming(tag + " CPU AMG-BiCGSTAB", tCpuSolve)
          .print();
    }
    if (Logger::hasDebug()) {
      const std::string path = " [CPU]";
      Logger::getInstance()
          .addTiming("diffusion n=" + std::to_string(n) + path +
                         " diag/b precompute",
                     tDiag)
          .print();
    }
    if (!lastSolveConverged_) {
      char numbers[96];
      std::snprintf(numbers, sizeof(numbers), "residual=%.3e, tolerance=%.1e",
                    double(normalizedResidual_), double(parameters.tolerance));
      VIENNACORE_LOG_WARNING(
          "solveDiffusion: AMG-BiCGSTAB did not converge after " +
          std::to_string(iterations) + "/" +
          std::to_string(parameters.maxIterations) + " iterations (" +
          numbers + ")");
    }
  }

#ifdef VIENNALS_GPU_BICGSTAB
  static std::string gpuErrorDetail() {
    const char *detail = gpu::gpuGetLastErrorMessage();
    if (detail && detail[0] != '\0')
      return std::string(" Detail: ") + detail;
    return {};
  }

  /// Reports that the GPU solver cannot be used.
  /// GpuMode::Gpu means "GPU or nothing", so this raises an error (which
  /// aborts). GpuMode::Auto asked for the GPU only if it happens to work, so
  /// it warns and lets the caller continue on to the CPU solver.
  void reportGpuUnavailable(const std::string &message) const {
    if (gpuMode_ == GpuMode::Auto) {
      VIENNACORE_LOG_WARNING(message + gpuErrorDetail() +
                             " Falling back to the CPU solver.");
    } else {
      VIENNACORE_LOG_ERROR(message + gpuErrorDetail());
    }
  }
#endif

  void logDiffusionBackend(const std::string &backend,
                           const std::string &detail) const {
    if (!Logger::hasInfo())
      return;
    const std::string msg =
        "OxidationDiffusion: using " + backend +
        " for diffusion solve (nodes=" + std::to_string(nodes.size()) +
        (detail.empty() ? std::string() : ", " + detail) + ").";
    if (msg == lastLoggedBackend_)
      return;
    lastLoggedBackend_ = msg;
    Logger::getInstance().addInfo(msg).print();
  }

  // Uses precomputed flat faceBC arrays — no HRLE access, safe for parallel
  // execution.
  template <class SolverT>
  StencilSide
  makeStencilSide(std::size_t nodeId, const std::vector<SolverT> &previous,
                  unsigned direction, int offset, T diffusion) const {
    const IndexType &nodeIndex = nodes[nodeId].index;
    offset = stencilOffset(nodeIndex, direction, offset); // mirror sides
    IndexType neighbor = nodeIndex;
    neighbor[direction] += offset;

    if (!inBounds(neighbor))
      return zeroFluxSide();

    const std::size_t neighborId = lookupNode(neighbor);
    if (neighborId != noNode)
      return {gridDelta, 0., previous[neighborId]};

    const unsigned fi = direction * 2u + (offset == 1 ? 1u : 0u);
    const std::size_t n = nodes.size();
    const Boundary faceType = faceBCTypes_[fi * n + nodeId];
    const T faceDist = faceBCDists_[fi * n + nodeId];
    if (faceType == Boundary::REACTION)
      return reactionBoundarySide(faceRate(nodeId, fi), faceDist, diffusion,
                                  faceReactionNormal_[fi * n + nodeId]);
    if (faceType == Boundary::AMBIENT)
      return ambientBoundarySide(faceDist, diffusion);
    if (faceType == Boundary::MASK)
      return maskBoundarySide(faceDist, diffusion);
    return zeroFluxSide();
  }

  void addAxisContribution(T &rightHandSide, T &diagonal,
                           const StencilSide &negativeSide,
                           const StencilSide &positiveSide, T diffusion) const {
    const T distanceSum = negativeSide.distance + positiveSide.distance;
    if (distanceSum <= std::numeric_limits<T>::epsilon())
      return;

    addSideContribution(rightHandSide, diagonal, negativeSide, distanceSum,
                        diffusion);
    addSideContribution(rightHandSide, diagonal, positiveSide, distanceSum,
                        diffusion);
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

  // Robin condition D dC/dn = k C on one axis arm that crosses the Si
  // interface (distance = the crossing along that axis).  Along the arm only
  // the axis component of the normal flux is available, D dC/dx_a =
  // |n_a| k C, so the rate is scaled by |n_a|.  Imposing the full k on every
  // crossing arm sums to a normal sink of k (|n_x| + |n_y| (+ |n_z|)): exact
  // on an axis-aligned front, sqrt(2) too strong at 45 degrees, 4/pi on
  // average around a cylinder — measured as 3-7% too little Si consumed by
  // the Kao cylinders with no stress coupling (2026-10-21).  With the
  // scaling the reconstructed interface concentration is also the same on
  // every arm, C_node (D/d) / (D/d + k) for the normal distance d.
  StencilSide reactionBoundarySide(T rate, T distance, T diffusion,
                                   T normalComponent) const {
    const T conductance = diffusion / distance;
    const T reactionRate = rate * normalComponent;
    const T denominator = conductance + reactionRate;
    if (denominator <= std::numeric_limits<T>::epsilon())
      return zeroFluxSide();

    return {distance, conductance / denominator, 0.};
  }

  ReactionBoundarySample reactionBoundarySample(const IndexType &index) const {
    ConstSparseIterator reactionIt(reactionInterface->getDomain());

    const std::size_t directId = lookupNode(index);
    if (directId != noNode) {
      const auto sample =
          reactionBoundarySampleFromNode(reactionIt, nodes[directId]);
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
          const std::size_t foundId = nodeLookupFlat[linearIndex(candidate)];
          if (foundId != noNode && distance2 < bestDistance2) {
            const auto sample =
                reactionBoundarySampleFromNode(reactionIt, nodes[foundId]);
            if (sample.found) {
              bestDistance2 = distance2;
              bestNode = foundId;
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
    // A meshed mask node never reacts with silicon: its reaction faces are
    // zero-flux, so it must not supply a reaction-boundary sample either.
    if (node.inMask)
      return best;
    best.nodeIndex = node.index;
    T bestDistance = std::numeric_limits<T>::max();
    const T insidePhi = valueAt(reactionIt, node.index);
    const std::size_t nodeId = lookupNode(node.index);
    const std::size_t n = nodes.size();
    std::optional<Vec3D<T>> siNormal;

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
        const unsigned fi = direction * 2u + (offset == 1 ? 1u : 0u);
        T normalComponent, rate;
        if (nodeId != noNode && faceReactionNormal_.size() == 2 * D * n &&
            faceBCTypes_[fi * n + nodeId] == Boundary::REACTION) {
          normalComponent = faceReactionNormal_[fi * n + nodeId];
          rate = faceRate(nodeId, fi);
        } else {
          if (!siNormal)
            siNormal = computeSiNormal(node.index, reactionIt);
          normalComponent = std::abs((*siNormal)[direction]);
          rate = getEffectiveReactionRate(node.index);
        }
        const auto side =
            reactionBoundarySide(rate, distance, diffusion, normalComponent);
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

  // Crossing of REACTION face fi of node i, in index units.
  Vec3D<T> crossingPoint(std::size_t i, unsigned fi) const {
    const std::size_t n = nodes.size();
    Vec3D<T> c{T(0), T(0), T(0)};
    for (unsigned d = 0; d < D; ++d)
      c[d] = static_cast<T>(nodes[i].index[d]);
    c[fi / 2] += (fi % 2 == 1 ? T(1) : T(-1)) * faceBCDists_[fi * n + i] /
                 gridDelta;
    return c;
  }

  // Calls f(j) for every node j in the box of half-width `radius` cells
  // around p (index units).
  template <class F> void forEachNodeNear(const Vec3D<T> &p, int radius,
                                          F f) const {
    IndexType centre, offset;
    for (unsigned d = 0; d < D; ++d)
      centre[d] = static_cast<viennahrle::IndexType>(std::llround(p[d]));
    offset.fill(-radius);
    while (true) {
      IndexType q = centre;
      for (unsigned d = 0; d < D; ++d)
        q[d] += offset[d];
      if (inBounds(q))
        if (const std::size_t j = lookupNode(q); j != noNode)
          f(j);
      unsigned d = 0;
      while (d < D && offset[d] == radius)
        offset[d++] = -radius;
      if (d == D)
        break;
      ++offset[d];
    }
  }

  // The law's normal stress at a front point p (index units): the reaction
  // nodes' smoothed samples (interfaceNormalStress_, each taken at its anchor
  // on the front), interpolated with Gaussian weights one cell wide.  NaN
  // where no reaction node within two cells carries one.
  T interfaceNormalStressAt(const Vec3D<T> &p) const {
    if (interfaceNormalStress_.size() != nodes.size() ||
        interfaceAnchor_.size() != nodes.size())
      return std::numeric_limits<T>::quiet_NaN();
    T weightSum = T(0), sum = T(0);
    forEachNodeNear(p, 2, [&](std::size_t j) {
      const T v = interfaceNormalStress_[j];
      if (!std::isfinite(v))
        return;
      T r2 = T(0);
      for (unsigned d = 0; d < D; ++d)
        r2 += (interfaceAnchor_[j][d] - p[d]) * (interfaceAnchor_[j][d] - p[d]);
      const T w = std::exp(-r2 / T(2));
      weightSum += w;
      sum += w * v;
    });
    return weightSum > T(1e-12) ? sum / weightSum
                                : std::numeric_limits<T>::quiet_NaN();
  }

  // Reaction rate at the crossing of REACTION face fi of node i: the
  // normal-stress law with the stress at the crossing, the orientation factor
  // with the crossing's normal.
  T interfaceReactionRate(std::size_t i, unsigned fi) const {
    const std::size_t n = nodes.size();
    T rate = parameters.reactionRate;
    if (parameters.reactionActivationVolume != T(0)) {
      T normalStress = interfaceNormalStressAt(crossingPoint(i, fi));
      if (!std::isfinite(normalStress)) {
        normalStress = -parameters.referencePressure;
        medianNormalStress(nodes[i].index, normalStress);
      }
      const T exponent =
          stressExponent(-normalStress, parameters.reactionActivationVolume);
      // floor 0.1 and ceiling 1: see getEffectiveReactionRate
      rate *= std::clamp(stressFactor(exponent), T(0.1), T(1.));
    }
    if (parameters.reactionRateRatio111 != T(1)) {
      const auto &normal = faceCrossingNormal_[fi * n + i];
      T dot = T(0);
      for (unsigned d = 0; d < D; ++d)
        dot += normal[d] * parameters.crystalAxis[d];
      rate *= T(1) +
              (parameters.reactionRateRatio111 - T(1)) * (T(1) - dot * dot);
    }
    return rate;
  }

  void updateFaceRates() {
    const std::size_t n = nodes.size();
    faceReactionRate_.assign(2 * D * n, T(0));
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned fi = 0; fi < 2 * D; ++fi)
        if (faceBCTypes_[fi * n + i] == Boundary::REACTION)
          faceReactionRate_[fi * n + i] = interfaceReactionRate(i, fi);
  }

  T faceRate(std::size_t i, unsigned fi) const {
    const std::size_t n = nodes.size();
    return faceReactionRate_.size() == 2 * D * n
               ? faceReactionRate_[fi * n + i]
               : getEffectiveReactionRate(nodes[i].index);
  }

  // Front speed k C_i / (N alpha) at every REACTION face crossing, C_i the
  // crossing's interface concentration from the solved node value (the Robin
  // arm's own reconstruction).
  void updateFaceSpeeds() {
    const std::size_t n = nodes.size();
    faceFrontSpeed_.assign(2 * D * n, T(0));
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned fi = 0; fi < 2 * D; ++fi) {
        if (faceBCTypes_[fi * n + i] != Boundary::REACTION)
          continue;
        const T rate = faceRate(i, fi);
        const auto side = reactionBoundarySide(
            rate, faceBCDists_[fi * n + i], nodeDiffusionCoefficient(i),
            faceReactionNormal_[fi * n + i]);
        const T concentration =
            side.nodeCoefficient * nodes[i].concentration + side.constant;
        faceFrontSpeed_[fi * n + i] =
            std::abs(parameters.velocitySign) * rate * concentration /
            (parameters.oxidantMoleculeDensity *
             parameters.expansionCoefficient);
      }
  }

  T getEffectiveReactionRate(const IndexType &index) const {
    T rate = parameters.reactionRate;

    if (parameters.reactionActivationVolume != T(0)) {
      T normalStress = -parameters.referencePressure;
      const std::size_t nodeId = lookupNode(index);
      if (nodeId != noNode && nodeId < interfaceNormalStress_.size() &&
          std::isfinite(interfaceNormalStress_[nodeId]))
        normalStress = interfaceNormalStress_[nodeId];
      else
        medianNormalStress(index, normalStress);
      // -sigma_nn is the normal compression, measured like a pressure.
      const T exponent =
          stressExponent(-normalStress, parameters.reactionActivationVolume);
      // Reaction floor 0.1 (vs the generic 1e-6 clamp): the source term
      // scales with this factor, and single-cell stress spikes (e.g.
      // -15 GPa at a contact corner vs -2 GPa a cell away) otherwise make
      // the volume source orders-of-magnitude rough on the grid scale,
      // stalling the SIMPLE mechanics solve (2026-08-10: residual plateau
      // 0.03-0.10 vs 5e-3, dt collapse).
      // Ceiling 1 (2026-10-22): interface stress only RETARDS the reaction.
      // In Kao's cylinders the interface normal stress is compressive for
      // convex and concave Si alike, and there is no data for the reaction
      // under interface tension.  Locally, tension appears where the Si
      // front dents: a dent in a convex front under the oxide's tensile hoop
      // stress is pulled open (membrane effect), so a tension-ENHANCED rate
      // deepens it — measured: Kao pillar (r = 0.45 um, Vr = 15 A^3) dented
      // 12.7 nm in 0.1 hr with ceiling 10, round within 0.6 nm with ceiling 1.
      // (The earlier ceiling of 10 was itself a guard against tensile
      // single-cell spikes at degenerate cut cells, 2026-08-12/13.)
      rate *= std::clamp(stressFactor(exponent), T(0.1), T(1.));
    }

    if (parameters.reactionRateRatio111 != T(1)) {
      const std::size_t nodeId = lookupNode(index);
      if (nodeId != noNode) {
        const auto &normal = nodes[nodeId].siNormal;
        T dot = T(0);
        for (unsigned d = 0; d < D; ++d)
          dot += normal[d] * parameters.crystalAxis[d];
        rate *= T(1) +
                (parameters.reactionRateRatio111 - T(1)) * (T(1) - dot * dot);
      }
    }

    return rate;
  }

  // Median of n.sigma.n over the node and its oxide face neighbours, all
  // projected on the node's own interface normal (2026-08-12 option A): cut
  // cells carry the least trustworthy stress in the field (single-cell spikes
  // of -15 GPa next to -2 GPa), and with up to 2D+1 samples a lone spike
  // cannot move the median, while a genuine gradient passes unchanged.
  bool medianNormalStress(const IndexType &index, T &result) const {
    return medianNormalStress(index, stressLookup, result);
  }

  // Oxide stress at the Si front for the stress laws, extrapolated along the
  // front normal n (into the oxide) from the points 2 and 3 cells into the
  // oxide, each D-linear over the nodes around it.  The first node layer
  // (the nodes with a reaction face) is left out: its MAC stress is O(1)
  // wrong next to the front -- on a Kao cylinder the pressure 0.2-0.4 GPa too
  // high and sigma_rr 0.1-0.2 GPa too compressive, at any grid spacing --
  // while the second layer is good to ~0.05 GPa.  Fixed distances along the
  // normal keep the estimate independent of the front's orientation on the
  // grid (a least-squares fit over the nodes behind the front was not).
  // phi: the node's distance from the front in cells; returns the value at
  // the distance sAt.  sample(j) gives the quantity at node j, NaN where
  // missing.  False where a sample point lacks oxide nodes (a film thinner
  // than ~4 cells); the caller then falls back to the node layer.
  template <class Sample>
  bool frontStressAt(std::size_t i, const Vec3D<T> &n, T phi, T sAt,
                     Sample sample, T &out) const {
    const std::size_t count = nodes.size();
    T q[2];
    for (unsigned k = 0; k < 2; ++k) {
      const T s = T(2 + k);
      std::array<T, D> base{}, frac{};
      for (unsigned d = 0; d < D; ++d) {
        const T y = static_cast<T>(nodes[i].index[d]) + n[d] * (s - phi);
        base[d] = std::floor(y);
        frac[d] = y - base[d];
      }
      T sum = T(0), weightSum = T(0);
      for (unsigned corner = 0; corner < (1u << D); ++corner) {
        IndexType idx;
        T w = T(1);
        for (unsigned d = 0; d < D; ++d) {
          const bool up = (corner >> d) & 1u;
          idx[d] = static_cast<viennahrle::IndexType>(base[d]) + (up ? 1 : 0);
          w *= up ? frac[d] : T(1) - frac[d];
        }
        if (w <= T(0) || !inBounds(idx))
          continue;
        const std::size_t j = lookupNode(idx);
        if (j == noNode || nodes[j].inMask)
          continue;
        bool firstLayer = false;
        if (faceBCTypes_.size() == 2 * D * count)
          for (unsigned f = 0; f < 2 * D && !firstLayer; ++f)
            firstLayer = faceBCTypes_[f * count + j] == Boundary::REACTION;
        if (firstLayer)
          continue;
        const T v = sample(j);
        if (!std::isfinite(v))
          continue;
        sum += w * v;
        weightSum += w;
      }
      if (weightSum < T(0.5))
        return false;
      q[k] = sum / weightSum;
    }
    out = q[0] + (q[1] - q[0]) * (sAt - T(2));
    return std::isfinite(out);
  }

  // Distance (cells) from the Si front into the oxide at a node next to it.
  T frontDistance(const IndexType &index) const {
    ConstSparseIterator it(reactionInterface->getDomain());
    it.goToIndices(index);
    const T v = it.isDefined() ? static_cast<T>(reactionSign) * it.getValue()
                               : T(0);
    return std::isfinite(v) ? std::clamp(v, T(0), T(1.5)) : T(0);
  }

  bool reactionNode(std::size_t i) const {
    const std::size_t n = nodes.size();
    if (faceBCTypes_.size() != 2 * D * n)
      return false;
    for (unsigned f = 0; f < 2 * D; ++f)
      if (faceBCTypes_[f * n + i] == Boundary::REACTION)
        return true;
    return false;
  }

  bool medianNormalStress(const IndexType &index, const StressLookup &lookup,
                          T &result) const {
    Vec3D<T> normal{T(0), T(0), T(0)};
    normal[D - 1] = T(1);
    if (const std::size_t centerNode = lookupNode(index); centerNode != noNode)
      normal = nodes[centerNode].siNormal;
    T samples[2 * D + 1];
    unsigned count = 0;
    const auto sampleAt = [&](const IndexType &sampleIndex) {
      // Only genuine oxide nodes carry an oxide stress (not the gas-side
      // band, not the mask body).
      const std::size_t sampleNode = lookupNode(sampleIndex);
      if (sampleNode == noNode || nodes[sampleNode].inMask)
        return;
      const auto found = lookup.find(detail::gridIndexHash<D>(sampleIndex));
      if (found == lookup.end())
        return;
      T value = T(0);
      for (unsigned r = 0; r < D; ++r)
        for (unsigned c = 0; c < D; ++c)
          value += normal[r] * found->second[3 * r + c] * normal[c];
      if (std::isfinite(value))
        samples[count++] = value;
    };
    sampleAt(index);
    for (unsigned d = 0; d < D; ++d) {
      IndexType neighbor = index;
      neighbor[d] -= 1;
      sampleAt(neighbor);
      neighbor[d] += 2;
      sampleAt(neighbor);
    }
    if (count == 0)
      return false;
    std::nth_element(samples, samples + count / 2, samples + count);
    result = samples[count / 2];
    return true;
  }

  // Interface normal stress for the reaction law: the median sample at every
  // node with a REACTION face, then smoothed along the interface.
  //
  // Why smoothed (2026-10-22): the reaction law k = k0 exp(sigma_nn V/kT)
  // feeds back through the mechanics, and that loop is stabilising but STIFF.
  // A rate ripple of wavenumber q makes a normal stress ~ 2 eta_impl q v
  // (eta_impl = (1 - e^(-dt/tau)) eta ~ G dt), so the loop gain is
  // g(q) ~ (V/kT) 2 G dt q v ~ (V/kT) G q h at the CFL step: ~80 for the
  // grid ripple at V = 15 A^3, G = 30 GPa, independent of h.  The partitioned
  // coupling cannot converge gains much above 10 (measured on a flat film:
  // undamped iterates alternate in sign and grow 4-5x per iteration; with the
  // loop exiting after one iteration the ripple doubled every step until the
  // mechanics failed).  Six passes of neighbour averaging along the interface
  // (a ~2h Gaussian) cap the gain at ~8 for every q.  It only removes stress
  // structure shorter than ~2 cells, which the converged solution suppresses
  // anyway (the stiffness is the physics damping such ripples), and it
  // vanishes with grid refinement.
  void updateInterfaceNormalStress() {
    interfaceNormalStress_ = interfaceNormalStressFor(stressLookup);
  }

  // Pressure the diffusivity law uses at the first node layer (the nodes
  // with a reaction face), extrapolated from behind it (frontStressAt); NaN
  // elsewhere, where the node's own pressure is used.
  void updateLawPressure() {
    const std::size_t n = nodes.size();
    lawPressure_.assign(n, std::numeric_limits<T>::quiet_NaN());
    ConstSparseIterator reactionIt(reactionInterface->getDomain());
    for (std::size_t i = 0; i < n; ++i) {
      if (nodes[i].inMask || !reactionNode(i))
        continue;
      // the node normal exists where a reaction law or the orientation
      // factor needs it
      Vec3D<T> into = (parameters.reactionRateRatio111 != T(1) ||
                       parameters.reactionActivationVolume != T(0))
                          ? nodes[i].siNormal
                          : computeSiNormal(nodes[i].index, reactionIt);
      for (unsigned d = 0; d < D; ++d)
        into[d] *= static_cast<T>(reactionSign);
      const T phi = frontDistance(nodes[i].index);
      const auto pressureAt = [&](std::size_t j) {
        const auto found =
            pressureLookup.find(detail::gridIndexHash<D>(nodes[j].index));
        return found == pressureLookup.end()
                   ? std::numeric_limits<T>::quiet_NaN()
                   : found->second;
      };
      T value;
      if (frontStressAt(i, into, phi, phi, pressureAt, value))
        lawPressure_[i] = value;
    }
  }

  /// Diffusion coefficient of a node's own material: D in the oxide,
  /// D / gamma inside a meshed mask.
  T nodeDiffusionCoefficient(std::size_t nodeId) const {
    T value = getEffectiveDiffusionCoefficient(nodes[nodeId].index);
    if (nodes[nodeId].inMask)
      value /= static_cast<T>(parameters.maskPermeabilityRatio);
    return value;
  }

  /// Diffusion coefficient on the face between a node and its neighbour: the
  /// HARMONIC mean of the two nodes' D, the value that makes the flux
  /// continuous across a jump (the mask/oxide interface) and conservative
  /// where D varies with pressure.  With each node's own D on its faces (the
  /// former rule inside one material) the flux leaving a node was not the
  /// flux entering its neighbour; a Kao cylinder then captured 73 % of the
  /// tension enhancement of the oxidant diffusivity (the 1-D model's).
  T faceDiffusionCoefficient(std::size_t nodeId, unsigned direction, int offset,
                             T ownDiffusion) const {
    offset = stencilOffset(nodes[nodeId].index, direction, offset);
    IndexType neighbor = nodes[nodeId].index;
    neighbor[direction] += offset;
    if (!inBounds(neighbor))
      return ownDiffusion;
    const std::size_t neighborId = lookupNode(neighbor);
    if (neighborId == noNode)
      return ownDiffusion;
    const T neighborDiffusion = nodeDiffusionCoefficient(neighborId);
    const T sum = ownDiffusion + neighborDiffusion;
    if (sum <= std::numeric_limits<T>::epsilon())
      return T(0);
    return T(2) * ownDiffusion * neighborDiffusion / sum;
  }

  T getEffectiveDiffusionCoefficient(const IndexType &index) const {
    if (parameters.diffusionActivationVolume == T(0))
      return parameters.diffusionCoefficient;

    T pressure = parameters.referencePressure;
    const auto found = pressureLookup.find(detail::gridIndexHash<D>(index));
    if (found != pressureLookup.end())
      pressure = found->second;
    if (const std::size_t id = lookupNode(index);
        id != noNode && id < lawPressure_.size() &&
        std::isfinite(lawPressure_[id]))
      pressure = lawPressure_[id];
    if (!std::isfinite(pressure))
      pressure = parameters.referencePressure;

    const T exponent =
        stressExponent(pressure, parameters.diffusionActivationVolume);
    return parameters.diffusionCoefficient * stressFactor(exponent);
  }

  T stressExponent(T pressure, T activationVolume) const {
    const T thermalEnergy =
        boltzmannConstant * std::max(parameters.temperature, T(1.));
    // Compression (p > 0) suppresses, tension enhances.
    return -(pressure - parameters.referencePressure) * activationVolume /
           thermalEnergy;
  }

  static T stressFactor(T exponent) {
    if (!std::isfinite(exponent))
      return T(1);
    if (exponent <= std::log(minStressFactor))
      return minStressFactor;
    if (exponent >= std::log(maxStressFactor))
      return maxStressFactor;
    return std::exp(exponent);
  }

  Vec3D<T> computeSiNormal(const IndexType &index,
                           ConstSparseIterator &reactionIt) const {
    // Reflect neighbor indices that fall outside the HRLE grid.  This handles
    // REFLECTIVE boundary conditions correctly: phi(b-k) = phi(b+k), so the
    // centered difference at b gives zero lateral gradient as expected.
    auto reflectToGrid = [&](IndexType idx) {
      auto &g = reactionInterface->getGrid();
      for (unsigned d2 = 0; d2 < D; ++d2) {
        const auto lo = g.getMinGridPoint(d2);
        const auto hi = g.getMaxGridPoint(d2);
        if (idx[d2] < lo)
          idx[d2] = 2 * lo - idx[d2];
        if (idx[d2] > hi)
          idx[d2] = 2 * hi - idx[d2];
      }
      return idx;
    };
    Vec3D<T> gradient = this->definedLevelSetGradient(
        reactionIt, index, [&](const IndexType &q, bool &ok) {
          ok = true;
          return reflectToGrid(q);
        });
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

  std::pair<Boundary, T> classifyBoundary(ConstSparseIterator &reactionIt,
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
    // GeometricAdvect can leave a tiny positive residual (~4*epsilon) when the
    // interface lands exactly on a grid point at non-zero coordinates, because
    // k*gridDelta is not exactly representable in floating point.  Allow a
    // tolerance of 1e-9 grid units so that grid points on the surface (phi≈0)
    // are correctly classified as inside the oxide.
    constexpr T eps = T(1e-9);
    return reactionSign * reactionPhi >= -eps &&
           ambientSign * ambientPhi >= -eps;
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
    return detail::levelSetCrossingDistance(
        insidePhi, outsidePhi, parameters.minBoundaryDistance, gridDelta);
  }
};

} // namespace viennals
