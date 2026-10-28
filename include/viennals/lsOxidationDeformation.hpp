#pragma once

#include <lsOxidationDiffusion.hpp>
#include <lsOxidationSolverBase.hpp>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>

#include <omp.h>
#include <vcTimer.hpp>

#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
// AMGCL's solver parameters as plain structs (no Boost property trees).
#ifndef AMGCL_NO_BOOST
#define AMGCL_NO_BOOST
#endif
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/preconditioner/schur_pressure_correction.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/fgmres.hpp>
#include <amgcl/solver/preonly.hpp>

namespace viennals {

/// Parameters for the Cartesian-grid oxide deformation model.
struct OxidationDeformationParameters {
  double viscosity = 1.;
  double ambientPressure = 0.;
  // Floor on the distance from a grid point to an interface crossing, in
  // grid spacings.  Interface ghost values extrapolate with 1/theta (theta
  // this distance), so the floor bounds their amplification where a crossing
  // lands next to a grid point.
  double minMechanicsBoundaryDistance = 0.20;
  // Maxwell oxide: shear modulus G (Pa), or the relaxation time tau (hr,
  // overrides eta/G).  Per substep dt the deviatoric stress follows
  //   S_end = decay S_start + (1 - decay) 2 eta D',  decay = exp(-dt/tau)
  // (exact for constant D'); S_start is the stored stress, carried with the
  // oxide.  The momentum solve treats (1 - decay) eta implicitly.  G = tau =
  // 0: purely viscous oxide.
  double shearModulus = 0.;
  double stressRelaxationTime = 0.;
  double stressTimeStep = 1.;
  // Stress-dependent oxide viscosity (Sutardja & Oldham, IEEE TED 36, 2415
  // (1989)): eta(tau) = eta0 x / sinh(x), x = tau Vc / (2 k T), with tau the
  // shear stress sqrt(J2) of the oxide's deviatoric stress and T the
  // oxidation temperature.  Below the critical shear stress 2kT/Vc (0.117 GPa
  // for Vc = 300 A^3 at 1000 C) the viscosity stays near eta0; above it, it
  // falls exponentially.  This bounds the stress at sharp concave corners,
  // where a linear viscous incompressible oxide has a pressure singularity
  // that grows as the grid is refined.  The Maxwell relaxation time scales
  // with it (tau = eta / G).  The momentum equation carries the full
  // div(2 eta D') with the local viscosity.  Vc in m^3; 0 = linear viscosity.
  double shearActivationVolume = 0.;
  // Viscosity iteration of the stress-dependent oxide (one exact velocity-
  // pressure solve per iteration, Anderson-accelerated): at most
  // mechanicsIterations, until the relative velocity change between two
  // iterations is below mechanicsTolerance.
  unsigned mechanicsIterations = 200;
  double mechanicsTolerance = 1e-4;
  // Monolithic oxide + nitride mechanics.  nitrideYoungModulus
  // > 0 together with a mask interface makes the grid points inside the mask
  // join the mechanics solve as a nitride body with stress memory
  // (linear-elastic, or Maxwell with nitrideViscosity).  The two materials
  // share nodes and interface faces, so velocity and traction are continuous
  // across the oxide/nitride interface by construction and no coupling
  // iteration is needed.  At a reflective domain boundary (the symmetry plane)
  // the nitride, like the oxide, has no normal velocity and no shear.  A mask
  // interface requires it.
  double nitrideYoungModulus = 0.;
  double nitridePoissonRatio = 0.27;
  // Maxwell nitride (Senez et al., IEEE TED 43, 720 (1996)): the deviatoric
  // stress relaxes over tau = eta_N / mu, the mean stress stays elastic.
  // nitrideViscosity in Pa hr; 0 = linear-elastic nitride.  With
  // nitrideShearActivationVolume (m^3) > 0 the viscosity falls with the
  // nitride's shear stress like the oxide's (Eyring, eta_N x / sinh x).
  double nitrideViscosity = 0.;
  double nitrideShearActivationVolume = 0.;
  std::size_t maxGridPoints = 5000000;
  int material = -1;
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
/// along the reaction-interface normal. The oxide (incompressible Maxwell
/// fluid with a stress-dependent viscosity) and an optional nitride mask
/// (incremental linear elasticity, monolithic mechanics) share one
/// quasi-static Stokes problem on a staggered (MAC) grid over the level-set
/// grid points, solved coupled in velocity and pressure: sparse LU in 2-D,
/// AMGCL in 3-D. The Si/SiO2 interface uses the oxidation expansion velocity,
/// the oxide/ambient and nitride/ambient interfaces are traction-free.
template <class T, int D>
class OxidationDeformation final : public VelocityField<T>,
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
  using OxidationSolverBase<T, D>::initializeGridFromInterfaces;

  enum class Boundary { NONE, REACTION, AMBIENT, ANCHOR };

  struct BoundaryIntersection {
    Boundary boundary = Boundary::NONE;
    T distance = 0.;
  };

  struct Node {
    IndexType index;
    Vec3D<T> velocity{0., 0., 0.};
    T pressure = 0.;
    T strainTrace = 0.;
    std::array<T, 9> strainRateTensor{};
    std::array<T, 9> stressTensor{};
    T vonMisesStress = 0.;
    bool nitride = false; // monolithic mode: node belongs to the nitride
  };

  SmartPointer<Domain<T, D>> reactionInterface = nullptr;
  SmartPointer<Domain<T, D>> ambientInterface = nullptr;
  SmartPointer<Domain<T, D>> maskInterface = nullptr;
  SmartPointer<OxidationDiffusion<T, D>> diffusionField = nullptr;
  OxidationDeformationParameters deformationParameters;
  OxidationParameters oxidationParameters;
  int reactionSign = 1;
  int ambientSign = -1;
  int maskSign = 1;

  IndexType requestedMinIndex{};
  IndexType requestedMaxIndex{};
  unsigned iterations = 0;
  T residual = std::numeric_limits<T>::max();
  // Volume error of the last solve (relative to the expansion inflow) and
  // backward error of its last linear solve.
  T lastContinuityResidual_ = 0.;
  T lastStokesResidual_ = 0.;
  T avgExpansionSpeed_ = 0.;
  bool avgExpansionSpeedComputed = false;
  bool solved = false;
  bool nodesDirty_ = true;
  std::array<T, D> maxVelocity_{};
  bool useRequestedBounds = false;
  std::unordered_map<IndexType, std::array<T, 9>, detail::IndexTypeHasher<D>>
      deviatoricStressHistory;
  // Stress-dependent oxide viscosity per node (shearActivationVolume > 0).
  std::vector<T> oxideEta_;
  // Stress-dependent nitride viscosity per node
  // (nitrideShearActivationVolume > 0).
  std::vector<T> nitrideEta_;
  // Monolithic nitride state: the step-start stress (from the mask level set;
  // per node in nitrideStartNode_), the end-of-step stress, and that stress
  // carried with the material.
  std::unordered_map<IndexType, std::array<T, 9>, detail::IndexTypeHasher<D>>
      nitrideStressStart_;
  std::vector<std::array<T, 9>> nitrideStartNode_;
  std::vector<std::array<T, 9>> nitrideStressEnd_;
  std::vector<std::array<T, 9>> nitrideStressAdvected_;
  bool nitrideStressAdvectedValid_ = false;
  // Maxwell oxide: the end-of-step deviatoric stress per node (shear
  // components the mean of the surrounding corners), and that stress carried
  // with the flow.
  std::vector<std::array<T, 9>> oxideStressEnd_;
  std::vector<std::array<T, 9>> oxideStressAdvected_;
  bool oxideStressAdvectedValid_ = false;
  // Per-node values that stay fixed during one apply(), cached because the
  // stencils need them for every row: the expansion velocity at Si faces and
  // the step-start oxide stress.
  // Valid only inside apply(); outside it they are recomputed per call.
  std::vector<Vec3D<T>> reactionVelocity_;
  std::vector<std::array<T, 9>> startStress_;
  bool boundaryCacheValid_ = false;

  static bool isFiniteVec(const Vec3D<T> &value) {
    for (unsigned i = 0; i < 3; ++i)
      if (!std::isfinite(value[i]))
        return false;
    return true;
  }

  static bool isFiniteTensor(const std::array<T, 9> &value) {
    for (const auto component : value)
      if (!std::isfinite(component))
        return false;
    return true;
  }

  // Face-major flat BC arrays: index = fi * n + nodeId, fi in [0, 2*D).
  std::vector<Boundary> faceBCTypes_;
  std::vector<T> faceBCDists_;

public:
  std::vector<Node> nodes;

  OxidationDeformation() = default;
#ifdef VIENNALS_GPU_BICGSTAB
  ~OxidationDeformation() { gpu::gpuMacFree(macGpu_); }
#endif

  OxidationDeformation(
      SmartPointer<Domain<T, D>> passedReactionInterface,
      SmartPointer<Domain<T, D>> passedAmbientInterface,
      SmartPointer<OxidationDiffusion<T, D>> passedDiffusionField,
      OxidationParameters passedOxidationParameters,
      OxidationDeformationParameters passedDeformationParameters = {})
      : reactionInterface(passedReactionInterface),
        ambientInterface(passedAmbientInterface),
        diffusionField(passedDiffusionField),
        deformationParameters(passedDeformationParameters),
        oxidationParameters(passedOxidationParameters) {}

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<OxidationDeformation>::New(std::forward<Args>(args)...);
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

  void setDiffusionField(
      SmartPointer<OxidationDiffusion<T, D>> passedDiffusionField) {
    diffusionField = passedDiffusionField;
    solved = false;
  }

  void setOxidationParameters(OxidationParameters passedParameters) {
    oxidationParameters = passedParameters;
    solved = false;
  }

  void
  setDeformationParameters(OxidationDeformationParameters passedParameters) {
    deformationParameters = passedParameters;
    solved = false;
  }

  /// While the stress coupling iterates (OxidationModel), a solve takes one
  /// pass of the viscosity iteration -- one linear solve -- and then updates
  /// the viscosity for the next call: the viscosities converge together with
  /// the stress feedback instead of inside every coupling iteration.  The
  /// solve after the coupling has converged runs the full iteration.
  void setCouplingIterate(bool iterate) { couplingIterate_ = iterate; }

  /// GPU for the 3-D velocity-pressure solve (systems of at least
  /// macGpuMinUnknowns unknowns): Cpu (default), Gpu, or Auto (the CPU when
  /// the GPU is unavailable). The 2-D sparse LU stays on the CPU.
  void setGpuMode(GpuMode mode) { gpuMode_ = mode; }

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
    nodesDirty_ = true;
    solved = false;
  }

  void clearSolveBounds() {
    useRequestedBounds = false;
    nodesDirty_ = true;
    solved = false;
  }

  void markGeometryChanged() {
    nodesDirty_ = true;
    solved = false;
  }

  // ===================== monolithic oxide + nitride =====================
  bool monolithicNitride() const {
    return deformationParameters.nitrideYoungModulus > 0. &&
           maskInterface != nullptr;
  }
  T stepTime() const {
    return std::max(T(deformationParameters.stressTimeStep),
                    std::numeric_limits<T>::min());
  }
  T nitridePoisson() const {
    return std::clamp(T(deformationParameters.nitridePoissonRatio), T(-0.95),
                      T(0.49));
  }
  T nitrideMu() const {
    return T(deformationParameters.nitrideYoungModulus) /
           (T(2) * (T(1) + nitridePoisson()));
  }
  T nitrideLambda() const {
    const T nu = nitridePoisson();
    return T(deformationParameters.nitrideYoungModulus) * nu /
           ((T(1) + nu) * (T(1) - T(2) * nu));
  }

  // Maxwell nitride: over one substep the deviator decays by
  // exp(-dt mu / eta) and responds with the implicit modulus
  // (1 - decay) eta (about mu dt for dt << tau); without a viscosity the
  // nitride is elastic (decay 1, modulus mu dt).
  bool viscousNitride() const {
    return deformationParameters.nitrideViscosity > 0.;
  }
  bool stressDependentNitride() const {
    return viscousNitride() &&
           deformationParameters.nitrideShearActivationVolume > 0.;
  }
  T nitrideDecayFor(T eta) const {
    if (!viscousNitride())
      return T(1);
    return std::exp(-stepTime() * nitrideMu() / eta);
  }
  T nitrideImplicitFor(T eta) const {
    if (!viscousNitride())
      return nitrideMu() * stepTime();
    return std::max(-std::expm1(-stepTime() * nitrideMu() / eta) * eta,
                    std::numeric_limits<T>::min() * T(1e10));
  }
  T nitrideNodeViscosity(std::size_t i) const {
    return (stressDependentNitride() && i < nitrideEta_.size())
               ? nitrideEta_[i]
               : T(deformationParameters.nitrideViscosity);
  }
  // Normal-stress moduli of nitride node i: sigma_dd = twoEta e_dd +
  // volumetric div v + memory (elastic: 2 mu dt and lambda dt).
  void nitrideModuli(std::size_t i, T &twoEta, T &volumetric) const {
    const T dt = stepTime();
    if (!viscousNitride()) {
      twoEta = T(2) * nitrideMu() * dt;
      volumetric = nitrideLambda() * dt;
      return;
    }
    const T etaImpl = nitrideImplicitFor(nitrideNodeViscosity(i));
    twoEta = T(2) * etaImpl;
    volumetric = (nitrideLambda() + T(2) * nitrideMu() / T(3)) * dt -
                 T(2) * etaImpl / T(3);
  }

  std::array<T, 9> nitrideStressStartAt(const IndexType &index) const {
    const auto it = nitrideStressStart_.find(index);
    return it == nitrideStressStart_.end() ? std::array<T, 9>{} : it->second;
  }
  const std::array<T, 9> &nitrideStartOf(std::size_t i) const {
    static const std::array<T, 9> zero{};
    return i < nitrideStartNode_.size() ? nitrideStartNode_[i] : zero;
  }

  BoundaryIntersection nitrideBoundaryIntersection(
      ConstSparseIterator &reactionIt, ConstSparseIterator &maskIt,
      const IndexType &inside, const IndexType &outside) const {
    const T maskInside = valueAtMask(maskIt, inside);
    const T maskOutside = valueAtMask(maskIt, outside);
    if (crosses(maskInside, maskOutside))
      return {Boundary::AMBIENT, crossingDistance(maskInside, maskOutside)};
    const T reactionInside = valueAt(reactionIt, inside);
    const T reactionOutside = valueAt(reactionIt, outside);
    if (crosses(reactionInside, reactionOutside)) // nitride on bare Si: weld
      return {Boundary::ANCHOR,
              crossingDistance(reactionInside, reactionOutside)};
    return {Boundary::AMBIENT, gridDelta};
  }

  // Nitride stress: the step-start stress (the momentum's memory term; the
  // nitride has no pressure) plus this step's increment.  Elastic, or
  // Maxwell (Senez 1996): the mean stress elastic, the deviator decayed,
  //   sigma_end = (m_start + K dt tr D) I + decay s_start + 2 eta_impl D'.
  void computeNitrideStress(std::size_t i) {
    auto &node = nodes[i];
    const auto &start = nitrideStartOf(i);
    const T dt = stepTime();
    const T mu = nitrideMu();
    const T lambda = nitrideLambda();
    T trace = T(0);
    for (unsigned d = 0; d < D; ++d)
      trace += node.strainRateTensor[tensorIndex(d, d)];
    std::array<T, 9> end = start;
    if (!viscousNitride()) {
      for (unsigned r = 0; r < 3; ++r) {
        for (unsigned c = 0; c < 3; ++c)
          end[tensorIndex(r, c)] +=
              dt * T(2) * mu * node.strainRateTensor[tensorIndex(r, c)];
        end[tensorIndex(r, r)] += dt * lambda * trace;
      }
    } else {
      const T decay = nitrideDecayFor(nitrideNodeViscosity(i));
      const T mean = (start[0] + start[4] + start[8]) / T(3);
      T twoEta, volumetric;
      nitrideModuli(i, twoEta, volumetric);
      for (unsigned r = 0; r < 3; ++r) {
        for (unsigned c = 0; c < 3; ++c)
          end[tensorIndex(r, c)] =
              decay * start[tensorIndex(r, c)] +
              twoEta * node.strainRateTensor[tensorIndex(r, c)];
        end[tensorIndex(r, r)] += (T(1) - decay) * mean + volumetric * trace;
      }
    }
    if (i < nitrideStressEnd_.size())
      nitrideStressEnd_[i] = end;
    nitrideStressAdvectedValid_ = false;
    node.stressTensor = end; // physical (output)
    const T meanStress = (end[0] + end[4] + end[8]) / T(3);
    auto deviator = end;
    for (unsigned d = 0; d < 3; ++d)
      deviator[tensorIndex(d, d)] -= meanStress;
    node.vonMisesStress = vonMisesFromDeviatoric(deviator);
  }

  // End-of-step stress carried with the material (semi-Lagrangian: the value
  // at x - v dt, D-linear over the nodes of the same material).  Stored at
  // grid points unadvected, the memory would slide relative to the moving
  // material.  Points whose upstream position has no node of the material
  // keep their own value.  Oxide whose upstream position lies in the Si was
  // formed during the step: no history (new oxide is born stress-free).
  // Renormalising over the oxide nodes there instead handed the first oxide
  // layer its own old stress every step, so it never took in new material
  // and its history grew without bound -- tension over the troughs of a
  // corrugated front that grew as 1/h.
  void advectStress(const std::vector<std::array<T, 9>> &end,
                    std::vector<std::array<T, 9>> &out, bool nitride) const {
    const std::size_t n = nodes.size();
    out = end;
    if (end.size() != n)
      return;
    const T dt = stepTime();
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < n; ++i) {
      if (nodes[i].nitride != nitride)
        continue;
      std::array<T, D> base{};
      std::array<T, D> frac{};
      std::array<T, D> departure{};
      for (unsigned d = 0; d < D; ++d) {
        const T x = static_cast<T>(nodes[i].index[d]) -
                    nodes[i].velocity[d] * dt / gridDelta;
        departure[d] = x;
        base[d] = std::floor(x);
        frac[d] = x - base[d];
      }
      if (!nitride && macPhiRAt(departure) < T(0)) {
        out[i] = std::array<T, 9>{};
        continue;
      }
      std::array<T, 9> sum{};
      T weightSum = T(0);
      for (unsigned corner = 0; corner < (1u << D); ++corner) {
        IndexType idx;
        T weight = T(1);
        for (unsigned d = 0; d < D; ++d) {
          const bool up = (corner >> d) & 1u;
          idx[d] = static_cast<viennahrle::IndexType>(base[d]) + (up ? 1 : 0);
          weight *= up ? frac[d] : (T(1) - frac[d]);
        }
        if (weight <= T(0) || !inBounds(idx))
          continue;
        const std::size_t j = nodeLookupFlat[linearIndex(idx)];
        if (j == noNode || nodes[j].nitride != nitride)
          continue;
        for (unsigned k = 0; k < 9; ++k)
          sum[k] += weight * end[j][k];
        weightSum += weight;
      }
      if (weightSum > T(1e-6))
        for (unsigned k = 0; k < 9; ++k)
          out[i][k] = sum[k] / weightSum;
    }
  }
  void advectNitrideStress() {
    advectStress(nitrideStressEnd_, nitrideStressAdvected_, true);
    nitrideStressAdvectedValid_ = true;
  }
  void advectOxideStress() {
    advectStress(oxideStressEnd_, oxideStressAdvected_, false);
    oxideStressAdvectedValid_ = true;
  }

  // Mean of field over the nodes of one material next to idx (face
  // neighbours, else the full 3^D box).  False if there are none.
  bool neighbourStress(const IndexType &idx,
                       const std::vector<std::array<T, 9>> &field, bool nitride,
                       std::array<T, 9> &out) const {
    if (field.size() != nodes.size())
      return false;
    std::array<T, 9> sum{};
    unsigned count = 0;
    auto add = [&](const IndexType &nb) {
      const std::size_t j = lookupNode(nb);
      if (j == noNode || nodes[j].nitride != nitride)
        return;
      for (unsigned k = 0; k < 9; ++k)
        sum[k] += field[j][k];
      ++count;
    };
    for (unsigned d = 0; d < D; ++d)
      for (int off : {-1, 1}) {
        IndexType nb = idx;
        nb[d] += off;
        add(nb);
      }
    if (count == 0) {
      unsigned corners = 1;
      for (unsigned d = 0; d < D; ++d)
        corners *= 3;
      for (unsigned code = 0; code < corners; ++code) {
        IndexType nb = idx;
        unsigned rest = code;
        bool centre = true;
        for (unsigned d = 0; d < D; ++d) {
          const int off = static_cast<int>(rest % 3) - 1;
          rest /= 3;
          nb[d] += off;
          centre = centre && off == 0;
        }
        if (!centre)
          add(nb);
      }
    }
    if (count == 0)
      return false;
    for (unsigned k = 0; k < 9; ++k)
      out[k] = sum[k] / static_cast<T>(count);
    return true;
  }

  // Persist nitride state on the mask level set: the advected stress as the
  // next step's history ("NitStressR*"), plus velocity and end-of-step stress
  // for visualisation ("MaskVelocity", "MaskStressR*", the names the elastic
  // mask solver writes).
  void writeNitrideFieldsToMask() {
    if (maskInterface == nullptr || nodes.empty())
      return;
    if (!nitrideStressAdvectedValid_)
      advectNitrideStress();
    using VD = typename PointData<T>::VectorDataType;
    VD hist0, hist1, hist2, vel, out0, out1, out2;
    ConstSparseIterator it(maskInterface->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const IndexType idx = it.getStartIndices();
      const std::size_t j = lookupNode(idx);
      if (j != noNode && nodes[j].nitride &&
          j < nitrideStressAdvected_.size()) {
        const auto &h = nitrideStressAdvected_[j];
        const auto &e = nitrideStressEnd_[j];
        hist0.push_back({h[0], h[1], h[2]});
        hist1.push_back({h[3], h[4], h[5]});
        hist2.push_back({h[6], h[7], h[8]});
        vel.push_back(isFiniteVec(nodes[j].velocity)
                          ? nodes[j].velocity
                          : Vec3D<T>{T(0), T(0), T(0)});
        out0.push_back({e[0], e[1], e[2]});
        out1.push_back({e[3], e[4], e[5]});
        out2.push_back({e[6], e[7], e[8]});
      } else {
        // Not a nitride node of this solve: a gas point the rising nitride
        // will cover, or a point outside the solve.  Carry the stress of the
        // adjacent nitride nodes so a newly covered surface layer keeps its
        // (bending) stress; else keep the stored value.
        std::array<T, 9> old{};
        if (!neighbourStress(idx, nitrideStressAdvected_, true, old))
          old = nitrideStressStartAt(idx);
        hist0.push_back({old[0], old[1], old[2]});
        hist1.push_back({old[3], old[4], old[5]});
        hist2.push_back({old[6], old[7], old[8]});
        vel.push_back({T(0), T(0), T(0)});
        out0.push_back({T(0), T(0), T(0)});
        out1.push_back({T(0), T(0), T(0)});
        out2.push_back({T(0), T(0), T(0)});
      }
    }
    auto &pd = maskInterface->getPointData();
    pd.insertReplaceVectorData(std::move(hist0), "NitStressR0");
    pd.insertReplaceVectorData(std::move(hist1), "NitStressR1");
    pd.insertReplaceVectorData(std::move(hist2), "NitStressR2");
    pd.insertReplaceVectorData(std::move(vel), "MaskVelocity");
    pd.insertReplaceVectorData(std::move(out0), "MaskStressR0");
    pd.insertReplaceVectorData(std::move(out1), "MaskStressR1");
    pd.insertReplaceVectorData(std::move(out2), "MaskStressR2");
  }

  void seedNitrideFromMask() {
    nitrideStressStart_.clear();
    if (!monolithicNitride())
      return;
    auto &pd = maskInterface->getPointData();
    const int r0 = pd.getVectorDataIndex("NitStressR0");
    const int r1 = pd.getVectorDataIndex("NitStressR1");
    const int r2 = pd.getVectorDataIndex("NitStressR2");
    if (r0 == -1 || r1 == -1 || r2 == -1)
      return; // first step: stress-free nitride
    const auto *d0 = pd.getVectorData(r0);
    const auto *d1 = pd.getVectorData(r1);
    const auto *d2 = pd.getVectorData(r2);
    ConstSparseIterator it(maskInterface->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const auto id = it.getPointId();
      if (id >= static_cast<decltype(id)>(d0->size()) ||
          id >= static_cast<decltype(id)>(d1->size()) ||
          id >= static_cast<decltype(id)>(d2->size()))
        continue;
      const auto &a = (*d0)[id];
      const auto &b = (*d1)[id];
      const auto &c = (*d2)[id];
      std::array<T, 9> st{a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]};
      if (isFiniteTensor(st))
        nitrideStressStart_[it.getStartIndices()] = st;
    }
    // Nitride nodes without a stored value: mean of stored face neighbours.
    std::vector<std::pair<IndexType, std::array<T, 9>>> fills;
    for (const auto &node : nodes) {
      if (!node.nitride || nitrideStressStart_.count(node.index))
        continue;
      std::array<T, 9> sum{};
      unsigned count = 0;
      for (unsigned d = 0; d < D; ++d)
        for (int off : {-1, 1}) {
          IndexType nb = node.index;
          nb[d] += off;
          const auto found = nitrideStressStart_.find(nb);
          if (found == nitrideStressStart_.end())
            continue;
          for (unsigned k = 0; k < 9; ++k)
            sum[k] += found->second[k];
          ++count;
        }
      if (count == 0)
        continue;
      for (unsigned k = 0; k < 9; ++k)
        sum[k] /= static_cast<T>(count);
      fills.push_back({node.index, sum});
    }
    for (const auto &fill : fills)
      nitrideStressStart_[fill.first] = fill.second;
  }
  // =======================================================================

  void apply() {
    if (reactionInterface == nullptr || ambientInterface == nullptr ||
        diffusionField == nullptr) {
      Logger::getInstance()
          .addError("OxidationDeformation: Missing interface or "
                    "diffusion field.")
          .print();
      return;
    }
    if (maskInterface != nullptr && !monolithicNitride())
      VIENNACORE_LOG_ERROR(
          "OxidationDeformation: a mask interface needs the monolithic "
          "nitride mechanics (nitrideYoungModulus > 0).");

    if (nodesDirty_) {
      if (!initialiseGrid())
        return; // base class already logged the error
      buildNodes();
      nodesDirty_ = false;
      if (nodes.empty())
        Logger::getInstance()
            .addWarning("OxidationDeformation: no oxide nodes found after "
                        "buildNodes(). Verify that the reaction and ambient "
                        "level sets enclose a non-empty oxide band.")
            .print();
      // Restore the stress history from the level-set point data (written
      // by writeFieldsToLevelSet() before the previous advection and
      // remapped + filled by lsAdvect + lsInterior).
      seedFromLevelSet();
      seedNitrideFromMask();
    }

    const std::size_t nn = nodes.size();
    cacheBoundaryValues();
    Timer<> tMechanics;
    tMechanics.start();
    solveMechanics();
    tMechanics.finish();
    boundaryCacheValid_ = false;
    Logger::getInstance()
        .addTiming("      deformation n=" + std::to_string(nn) +
                       " mechanics-total",
                   tMechanics)
        .print();
    avgExpansionSpeedComputed = false;

    maxVelocity_.fill(T(0));
    for (const auto &node : nodes) {
      for (unsigned d = 0; d < D; ++d)
        maxVelocity_[d] = std::max(maxVelocity_[d], std::abs(node.velocity[d]));
    }
    // The gas-side (unresolved) velocities enter the CFL step and the
    // advection's dissipation, i.e. only the solve the advection uses: a
    // stress-coupling iteration skips this walk over the ambient interface
    // (with stress kinetics it took ~30 % of a LOCOS run).
    if (!couplingIterate_) {
      const auto unresolvedMax = estimateMaxUnresolvedAmbientVelocity();
      for (unsigned d = 0; d < D; ++d)
        maxVelocity_[d] = std::max(maxVelocity_[d], unresolvedMax[d]);
    }

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

    auto velocity = getVelocity(coordinate);
    T norm2 = 0.;
    for (unsigned d = 0; d < D; ++d)
      norm2 += velocity[d] * velocity[d];
    bool resolved = norm2 > std::numeric_limits<T>::epsilon();

    // Gas-side advection queries land on grid points with no solver node at
    // any interpolation corner, so the D-linear lookup above comes back
    // empty.  Extend the solved velocity from the nearest oxide node instead
    // of dropping to the expansion-speed heuristic: that heuristic evaluates
    // the reaction speed from the LOCAL concentration, which at the oxide TOP
    // is ~C_eq — overestimating the true interface speed by C_eq/C_si.  Fed
    // into the upwind stencil, that inflated gas-side value advects the free
    // surface faster than the flow field and fabricates oxide volume — the
    // error grows with thickness as C_si drops (measured: expansion 2.35 at
    // 10 nm, 2.68 at 100 nm, target 2.27).
    if (!resolved) {
      IndexType index;
      for (unsigned d = 0; d < D; ++d)
        index[d] = std::llround(coordinate[d] / gridDelta);
      resolved =
          nearbyMeanVelocity(index, [](std::size_t) { return true; }, velocity);
    }

    const auto result =
        resolved ? velocity : unresolvedAmbientVelocity(coordinate);

    return result;
  }

  T getScalarVelocity(const Vec3D<T> &coordinate, int material,
                      const Vec3D<T> &normalVector,
                      unsigned long /*pointId*/) final {
    return 0.;
  }

  T getDissipationAlpha(int direction, int material,
                        const Vec3D<T> & /*centralDifferences*/) final {
    if (deformationParameters.material >= 0 &&
        material != deformationParameters.material)
      return 0.;
    return maxVelocity_[direction];
  }

private:
  template <class ValueType, class NodeAccessor>
  ValueType getField(const IndexType &index, ValueType fallback,
                     NodeAccessor accessor) const {
    const std::size_t nodeId = lookupNode(index);
    if (nodeId != noNode)
      return accessor(nodes[nodeId]);

    const auto nearby = findNearbyNode(index);
    if (nearby == noNode)
      return fallback;
    return accessor(nodes[nearby]);
  }

  template <class ValueType, class NodeAccessor>
  ValueType getField(const Vec3D<T> &coordinate, ValueType fallback,
                     NodeAccessor accessor) const {
    // D-linear interpolation for Vec3D<T> (velocity) and scalar T (pressure).
    // These are the types queried during level-set advection; interpolating at
    // the exact interface coordinate eliminates the O(gridDelta) nearest-node
    // error that caused the SiO2/mask interface to shift with grid delta.
    // Other types (e.g. stress tensor std::array<T,9>) use nearest-node because
    // they lack scalar multiply and are not used for advection.
    if constexpr (!std::is_same_v<ValueType, Vec3D<T>> &&
                  !std::is_same_v<ValueType, T>) {
      IndexType index;
      for (unsigned i = 0; i < D; ++i)
        index[i] = std::llround(coordinate[i] / gridDelta);
      return getField(index, fallback, accessor);
    } else {
      using IdxScalar = std::decay_t<decltype(std::declval<IndexType>()[0])>;
      IndexType lo;
      T frac[D];
      for (unsigned d = 0; d < D; ++d) {
        const T c = coordinate[d] / gridDelta;
        const T c_flo = std::floor(c);
        lo[d] = static_cast<IdxScalar>(c_flo);
        frac[d] = c - c_flo;
      }

      ValueType result{};
      T totalWeight = T(0);
      for (int corner = 0; corner < (1 << D); ++corner) {
        IndexType idx = lo;
        T w = T(1);
        for (unsigned d = 0; d < D; ++d) {
          if ((corner >> d) & 1) {
            idx[d]++;
            w *= frac[d];
          } else {
            w *= T(1) - frac[d];
          }
        }
        if (w < T(1e-14))
          continue;
        const std::size_t nodeId = lookupNode(idx);
        if (nodeId == noNode)
          continue;
        result = result + accessor(nodes[nodeId]) * w;
        totalWeight += w;
      }

      if (totalWeight < T(1e-14))
        return fallback;
      if (totalWeight < T(1) - T(1e-6))
        result = result * (T(1) / totalWeight);
      return result;
    }
  }

public:
  Vec3D<T> getVelocity(const Vec3D<T> &coordinate) const {
    return getField(coordinate, Vec3D<T>{0., 0., 0.},
                    [](const Node &n) { return n.velocity; });
  }

  /// Purely mechanical velocity: the interpolated Stokes solution, extended
  /// from the nearest solved node when the interpolation stencil has no node,
  /// and ZERO when nothing solved is nearby.  NEVER falls back to the
  /// concentration-derived expansion heuristic — consumers that move a
  /// carried mask must see only mechanics, since oxidant concentration has
  /// nothing to do with mask motion.
  Vec3D<T> getResolvedVectorVelocity(const Vec3D<T> &coordinate) const {
    auto velocity = getVelocity(coordinate);
    T norm2 = 0.;
    for (unsigned d = 0; d < D; ++d)
      norm2 += velocity[d] * velocity[d];
    if (norm2 > std::numeric_limits<T>::epsilon())
      return velocity;
    IndexType index;
    for (unsigned d = 0; d < D; ++d)
      index[d] = std::llround(coordinate[d] / gridDelta);
    if (nearbyMeanVelocity(index, [](std::size_t) { return true; }, velocity))
      return velocity;
    return {0., 0., 0.};
  }

  // Mean velocity of the nearest nodes passing accept(j) (findNearbyNodes);
  // false if there is none within 4 cells.
  template <class Accept>
  bool nearbyMeanVelocity(const IndexType &index, Accept accept,
                          Vec3D<T> &velocity) const {
    const auto nearby = this->findNearbyNodes(index, accept);
    if (nearby.empty())
      return false;
    velocity = Vec3D<T>{T(0), T(0), T(0)};
    for (const auto j : nearby)
      for (unsigned d = 0; d < D; ++d)
        velocity[d] += nodes[j].velocity[d] / static_cast<T>(nearby.size());
    return true;
  }

  Vec3D<T> getVelocity(const IndexType &index) const {
    return getField(index, Vec3D<T>{0., 0., 0.},
                    [](const Node &n) { return n.velocity; });
  }

  /// Monolithic contact: the velocity of one material only (oxide or
  /// nitride), for advecting that material's surface.  D-linear
  /// interpolation over that material's nodes, else its nearest node (within
  /// 4 cells), else zero for the nitride and the unresolved-ambient estimate
  /// for the oxide (as getVectorVelocity).  Interpolating across both
  /// materials advected the nitride's edge corner with the free oxide next
  /// to it (measured: corner thinned 50 -> 35 nm in 0.1 h).
  Vec3D<T> getMaterialVelocity(const Vec3D<T> &coordinate, bool nitride) {
    if (!solved)
      apply();
    using IdxScalar = std::decay_t<decltype(std::declval<IndexType>()[0])>;
    IndexType lo;
    T frac[D];
    for (unsigned d = 0; d < D; ++d) {
      const T c = coordinate[d] / gridDelta;
      const T cFloor = std::floor(c);
      lo[d] = static_cast<IdxScalar>(cFloor);
      frac[d] = c - cFloor;
    }
    Vec3D<T> result{T(0), T(0), T(0)};
    T totalWeight = T(0);
    for (int corner = 0; corner < (1 << D); ++corner) {
      IndexType idx = lo;
      T w = T(1);
      for (unsigned d = 0; d < D; ++d) {
        if ((corner >> d) & 1) {
          idx[d]++;
          w *= frac[d];
        } else {
          w *= T(1) - frac[d];
        }
      }
      if (w < T(1e-14))
        continue;
      const std::size_t j = lookupNode(idx);
      if (j == noNode || nodes[j].nitride != nitride)
        continue;
      for (unsigned d = 0; d < D; ++d)
        result[d] += w * nodes[j].velocity[d];
      totalWeight += w;
    }
    if (totalWeight >= T(1e-14)) {
      for (unsigned d = 0; d < D; ++d)
        result[d] /= totalWeight;
      return result;
    }
    IndexType index;
    for (unsigned d = 0; d < D; ++d)
      index[d] = std::llround(coordinate[d] / gridDelta);
    if (nearbyMeanVelocity(
            index, [&](std::size_t j) { return nodes[j].nitride == nitride; },
            result))
      return result;
    return nitride ? Vec3D<T>{T(0), T(0), T(0)}
                   : unresolvedAmbientVelocity(coordinate);
  }

  T getPressure(const Vec3D<T> &coordinate) const {
    return getField(coordinate, T(0), [](const Node &n) { return n.pressure; });
  }

  T getPressure(const IndexType &index) const {
    return getField(index, T(0), [](const Node &n) { return n.pressure; });
  }

  T getStrainTrace(const Vec3D<T> &coordinate) const {
    return getField(coordinate, T(0),
                    [](const Node &n) { return n.strainTrace; });
  }

  T getStrainTrace(const IndexType &index) const {
    return getField(index, T(0), [](const Node &n) { return n.strainTrace; });
  }

  std::array<T, 9> getStrainRateTensor(const Vec3D<T> &coordinate) const {
    return getField(coordinate, std::array<T, 9>{},
                    [](const Node &n) { return n.strainRateTensor; });
  }

  std::array<T, 9> getStrainRateTensor(const IndexType &index) const {
    return getField(index, std::array<T, 9>{},
                    [](const Node &n) { return n.strainRateTensor; });
  }

  std::array<T, 9> getStressTensor(const Vec3D<T> &coordinate) const {
    return getField(coordinate, std::array<T, 9>{},
                    [](const Node &n) { return n.stressTensor; });
  }

  std::array<T, 9> getStressTensor(const IndexType &index) const {
    return getField(index, std::array<T, 9>{},
                    [](const Node &n) { return n.stressTensor; });
  }

  T getVonMisesStress(const Vec3D<T> &coordinate) const {
    return getField(coordinate, T(0),
                    [](const Node &n) { return n.vonMisesStress; });
  }

  T getVonMisesStress(const IndexType &index) const {
    return getField(index, T(0),
                    [](const Node &n) { return n.vonMisesStress; });
  }

  unsigned getIterations() const { return iterations; }
  T getResidual() const { return residual; }
  // Volume the final flow field fails to conserve, relative to the expansion
  // inflow (see solveMechanics).
  T getLastContinuityResidual() const { return lastContinuityResidual_; }
  // Residual of the last velocity-pressure solve: componentwise backward
  // error (LU, 2-D), relative residual (AMGCL, 3-D).
  T getLastStokesResidual() const { return lastStokesResidual_; }
  // Largest velocity-pressure solve residual lastSolveConverged() accepts.
  static constexpr T stokesTolerance = T(1e-6);
  bool lastSolveConverged() const {
    return std::isfinite(residual) &&
           residual <= deformationParameters.mechanicsTolerance &&
           std::isfinite(lastStokesResidual_) &&
           lastStokesResidual_ <= stokesTolerance && hasFiniteSolution();
  }
  bool hasFiniteSolution() const {
    for (const auto &node : nodes) {
      if (!isFiniteVec(node.velocity) || !std::isfinite(node.pressure) ||
          !std::isfinite(node.strainTrace) ||
          !isFiniteTensor(node.strainRateTensor) ||
          !isFiniteTensor(node.stressTensor) ||
          !std::isfinite(node.vonMisesStress))
        return false;
    }
    return true;
  }
  std::size_t getNumberOfSolutionNodes() const { return nodes.size(); }
  T avgExpansionSpeed() {
    if (!avgExpansionSpeedComputed) {
      computeAvgExpansionSpeed();
      avgExpansionSpeedComputed = true;
    }
    return avgExpansionSpeed_;
  }
  /// callback(index, pressure, Cauchy stress) for every solution node.
  template <class Callback> void forEachSolutionNode(Callback callback) const {
    for (const auto &node : nodes)
      callback(node.index, node.pressure, node.stressTensor);
  }

  /// Write velocity (Vec3D) and viscoelastic stress history (3 tensor-row
  /// vectors) into ambientInterface->getPointData(). Called alongside
  /// OxidationDiffusion::writePersistentFields() so a single lsInterior pass
  /// fills the oxide interior with all warm-start data.
  void writeFieldsToLevelSet() {
    if (nodes.empty() || ambientInterface == nullptr)
      return;
    macAdvectCorners();

    using VD = typename PointData<T>::VectorDataType;
    VD velocity, stressR0, stressR1, stressR2;
    // The oxide's end-of-step stress, carried with the flow, is the next
    // substep's S_start.  Points that are not oxide nodes (gas band, newly
    // covered points, nitride) take the adjacent oxide nodes' value, else
    // keep their stored one.
    if (!oxideStressAdvectedValid_)
      advectOxideStress();

    ConstSparseIterator it(ambientInterface->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const IndexType idx = it.getStartIndices();
      const std::size_t nId = lookupNode(idx);

      velocity.push_back((nId != noNode && isFiniteVec(nodes[nId].velocity))
                             ? nodes[nId].velocity
                             : Vec3D<T>{T(0), T(0), T(0)});
      std::array<T, 9> stress{};
      if (nId != noNode && !nodes[nId].nitride &&
          nId < oxideStressAdvected_.size())
        stress = oxideStressAdvected_[nId];
      else if (nId == noNode && macPhiR(idx) < T(0))
        stress = {}; // Si: becomes oxide only as new, stress-free material
      else if (!neighbourStress(idx, oxideStressAdvected_, false, stress))
        stress = previousDeviatoricStress(idx);
      if (!isFiniteTensor(stress))
        stress = {};
      stressR0.push_back({stress[0], stress[1], stress[2]});
      stressR1.push_back({stress[3], stress[4], stress[5]});
      stressR2.push_back({stress[6], stress[7], stress[8]});
    }

    auto &pd = ambientInterface->getPointData();
    pd.insertReplaceVectorData(std::move(velocity), "OxVelocity");
    pd.insertReplaceVectorData(std::move(stressR0), "OxStressR0");
    pd.insertReplaceVectorData(std::move(stressR1), "OxStressR1");
    pd.insertReplaceVectorData(std::move(stressR2), "OxStressR2");
    if (stressDependentViscosity()) { // viscosity / eta0, for visualisation
      typename PointData<T>::ScalarDataType factor;
      for (ConstSparseIterator vit(ambientInterface->getDomain());
           !vit.isFinished(); ++vit) {
        if (!vit.isDefined())
          continue;
        const std::size_t nId = lookupNode(vit.getStartIndices());
        factor.push_back((nId != noNode && !nodes[nId].nitride)
                             ? nodeViscosity(nId) /
                                   T(deformationParameters.viscosity)
                             : T(1));
      }
      pd.insertReplaceScalarData(std::move(factor), "OxViscosity");
    }
    { // mechanics pressure (Pa, compression positive), for visualisation
      typename PointData<T>::ScalarDataType pressure;
      for (ConstSparseIterator pit(ambientInterface->getDomain());
           !pit.isFinished(); ++pit) {
        if (!pit.isDefined())
          continue;
        const std::size_t nId = lookupNode(pit.getStartIndices());
        pressure.push_back((nId != noNode && !nodes[nId].nitride &&
                            std::isfinite(nodes[nId].pressure))
                               ? nodes[nId].pressure
                               : T(0));
      }
      pd.insertReplaceScalarData(std::move(pressure), "OxMechanicsPressure");
    }

    if (monolithicNitride())
      writeNitrideFieldsToMask();
    macWriteCornerHistory();
  }

private:
  /// Reads the oxide stress history from ambientInterface pointData.  Called
  /// in apply() after buildNodes().
  void seedFromLevelSet() {
    if (ambientInterface == nullptr || nodes.empty())
      return;
    macReadCornerHistory();

    auto &pd = ambientInterface->getPointData();
    const int r0Idx = pd.getVectorDataIndex("OxStressR0");
    const int r1Idx = pd.getVectorDataIndex("OxStressR1");
    const int r2Idx = pd.getVectorDataIndex("OxStressR2");
    if (r0Idx == -1 || r1Idx == -1 || r2Idx == -1)
      return;
    const auto *r0d = pd.getVectorData(r0Idx);
    const auto *r1d = pd.getVectorData(r1Idx);
    const auto *r2d = pd.getVectorData(r2Idx);

    ConstSparseIterator it(ambientInterface->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const auto ptId = it.getPointId();
      const IndexType idx = it.getStartIndices();
      if (ptId < static_cast<decltype(ptId)>(r0d->size()) &&
          ptId < static_cast<decltype(ptId)>(r1d->size()) &&
          ptId < static_cast<decltype(ptId)>(r2d->size())) {
        const auto &row0 = (*r0d)[ptId];
        const auto &row1 = (*r1d)[ptId];
        const auto &row2 = (*r2d)[ptId];
        std::array<T, 9> s{row0[0], row0[1], row0[2], row1[0], row1[1],
                           row1[2], row2[0], row2[1], row2[2]};
        if (isFiniteTensor(s))
          deviatoricStressHistory[idx] = s;
      }
    }
  }

public:
  bool initialiseGrid() {
    return initializeGridFromInterfaces(
        reactionInterface, ambientInterface, maskInterface, useRequestedBounds,
        requestedMinIndex, requestedMaxIndex,
        deformationParameters.maxGridPoints, "OxidationDeformation");
  }

  void buildNodes() {
    ++nodesVersion_;
    nodes.clear();
    initNodeLookup();

    ConstSparseIterator reactionIt(reactionInterface->getDomain());
    ConstSparseIterator ambientIt(ambientInterface->getDomain());
    auto maskIt = makeMaskIterator();

    const bool mono = monolithicNitride();
    IndexType index = minIndex;
    while (true) {
      const T reactionPhi = valueAt(reactionIt, index);
      const T ambientPhi = valueAt(ambientIt, index);
      const bool inMask = isInsideMask(maskIt, index);
      if (isInsideOxide(reactionPhi, ambientPhi) && !inMask) {
        const std::size_t id = nodes.size();
        nodeLookupFlat[linearIndex(index)] = id;
        nodes.push_back({index});
      } else if (mono && inMask &&
                 static_cast<T>(reactionSign) * reactionPhi >= T(-1e-9)) {
        // Monolithic: the nitride body joins the mechanics solve.
        const std::size_t id = nodes.size();
        nodeLookupFlat[linearIndex(index)] = id;
        nodes.push_back({index});
        nodes.back().nitride = true;
      }

      if (!increment(index))
        break;
    }

    // Precompute per-face boundary intersections into flat face-major arrays.
    const std::size_t n = nodes.size();
    faceBCTypes_.assign(2 * D * n, Boundary::NONE);
    faceBCDists_.assign(2 * D * n, T(1));
    for (std::size_t id = 0; id < n; ++id) {
      const auto &node = nodes[id];
      for (unsigned dir = 0; dir < D; ++dir) {
        for (int off : {-1, 1}) {
          const unsigned fi = dir * 2u + (off == 1 ? 1u : 0u);
          IndexType nb = node.index;
          nb[dir] += off;
          if (!inBounds(nb))
            continue; // reflective domain boundary: the mirror plane
          if (lookupNode(nb) != noNode)
            continue; // NONE/1.0 already set
          BoundaryIntersection bi;
          if (node.nitride) {
            bi =
                nitrideBoundaryIntersection(reactionIt, maskIt, node.index, nb);
          } else {
            bi = boundaryIntersection(reactionIt, ambientIt, maskIt, node.index,
                                      nb);
          }
          faceBCTypes_[fi * n + id] = bi.boundary;
          faceBCDists_[fi * n + id] = bi.distance;
        }
      }
    }

    nitrideStressEnd_.assign(n, std::array<T, 9>{});
    nitrideStressAdvectedValid_ = false;
    oxideStressEnd_.assign(n, std::array<T, 9>{});
    oxideStressAdvectedValid_ = false;
  }

  // ============================ MAC mechanics ============================
  // Staggered (marker-and-cell) discretisation of the oxide + nitride
  // mechanics on the level-set grid: pressure at the grid points (oxide
  // nodes), velocity component d on the midpoint of the grid edge between a
  // and a + e_d ("face (d, a)"), shear stress sigma_ce at the corner (2-D) /
  // edge (3-D) a + (e_c + e_e)/2 ("corner (ce, a)").  Pressure differences
  // drive the face velocities directly, so no node-to-node pressure mode can
  // hide from the momentum equation (a collocated grid has one, which a
  // thick, stiff nitride drives to GPa amplitudes).  The full variable-
  // viscosity stress div(2 eta D') is discretised.  Interfaces come from the
  // level sets:
  //  - Si front: a face is an unknown if its midpoint lies on the oxide side
  //    (phi_R mean of its two grid points); a face in Si is linear through
  //    its partner face and the front velocity at the crossing, theta from
  //    phi_R, floor minMechanicsBoundaryDistance.
  //  - Free surface: every face next to a solid grid point is an unknown.
  //    The normal stress at a gas point is the solid point's stress,
  //    extrapolated to its traction-free value (P sigma P)_dd at the
  //    crossing (faceBCDists_), P = I - n n^T with the surface normal n from
  //    the defined level-set values: a free surface carries only tangential
  //    stress (extrapolating to zero strips it from the top cell row of every
  //    sloped surface).  A corner with two or more gas points lies on the
  //    surface and carries the traction-free shear (P sigma P)_ce of its
  //    solid points, plus a 1e-3 residual modulus that keeps a one-cell
  //    feature (a bump of the mask level set) attached; a gas face in its
  //    stencil takes its partner's value.  Full shear at corners with one gas
  //    point (steps of the surface), whose stencils need no gas faces.
  //  - Oxide / nitride: one velocity field; the corner modulus is the
  //    harmonic mean over the corner's solid points.
  //  - Domain boundary (reflective level-set boundary through the outermost
  //    grid points): normal velocity odd, tangential velocity even.
  // Shear-stress memory lives on the corners (macCornerStart_), carried with
  // the flow once per accepted step (macAdvectCorners) and stored on the
  // level sets between solves (macWriteCornerHistory).  Each viscosity
  // iteration solves the coupled velocity-pressure system exactly: sparse LU
  // in 2-D (ordering and symbolic analysis reused per node set), AMGCL in 3-D
  // (a direct solve does not scale there; in 2-D it is 2x faster than AMGCL
  // on the stiff-nitride LOCOS stack).

  struct MacLin { // affine expression sum t_k x_k + c
    std::vector<std::pair<std::int64_t, T>> t;
    T c = T(0);
    MacLin &add(const MacLin &o, T s = T(1)) {
      for (const auto &kv : o.t)
        t.emplace_back(kv.first, s * kv.second);
      c += s * o.c;
      return *this;
    }
    MacLin scaled(T s) const {
      MacLin r;
      r.t.reserve(t.size());
      for (const auto &kv : t)
        r.t.emplace_back(kv.first, s * kv.second);
      r.c = s * c;
      return r;
    }
  };

  static constexpr unsigned macPairs = static_cast<unsigned>(D * (D - 1) / 2);
  static unsigned macPairIndex(unsigned c, unsigned e) {
    if (c > e)
      std::swap(c, e);
    return (D == 2) ? 0u : (c == 0 ? (e == 1 ? 0u : 1u) : 2u);
  }
  static void macPairDirs(unsigned pair, unsigned &c, unsigned &e) {
    if (D == 2 || pair == 0) {
      c = 0;
      e = 1;
    } else if (pair == 1) {
      c = 0;
      e = 2;
    } else {
      c = 1;
      e = 2;
    }
  }

  enum : std::uint8_t { macGas = 0, macSi = 1, macOxide = 2, macNitride = 3 };

  // MAC state (rebuilt with the node set).
  std::size_t macNodesVersion_ = std::numeric_limits<std::size_t>::max();
  std::size_t nodesVersion_ = 0;
  std::array<std::vector<std::int64_t>, D> macPlusFace_; // face (d, node)
  std::array<std::vector<std::int64_t>, D>
      macMinusFace_;                 // face (d, node - e_d)
  std::vector<std::int64_t> macPid_; // per node
  std::int64_t macNVel_ = 0, macN_ = 0;
  std::unordered_map<std::size_t, T>
      macPhiR_; // phi_R (> 0 oxide side) near nodes
  // Inside-positive ambient / mask level sets near the nodes (NaN where the
  // level set is undefined): free-surface normals.
  std::unordered_map<std::size_t, T> macPsiA_, macPsiM_;
  // Corners: id per (pair, lower grid point); per corner data.
  std::vector<std::unordered_map<std::size_t, std::int64_t>> macCornerId_;
  std::vector<std::pair<unsigned, IndexType>>
      macCorners_; // (pair, lower point)
  // macCornerEta_ the oxide's corner viscosity, macCornerEtaN_ the nitride's
  std::vector<T> macCornerEta_, macCornerEtaN_, macCornerStartV_,
      macCornerEndV_, macCornerMod_, macCornerMem_;
  std::vector<std::uint8_t> macCornerGas_, macCornerOxide_;
  // Corner shear history across steps (pair, lower point) -> start stress.
  std::unordered_map<IndexType, std::array<T, 3>, detail::IndexTypeHasher<D>>
      macCornerStart_;
  std::unordered_map<IndexType, std::array<T, 3>, detail::IndexTypeHasher<D>>
      macCornerEndMap_;
  bool macCornerAdvectedValid_ = false;
  std::vector<T> macX_; // last solution
  // Front velocity of the Si ghost faces during one solveMechanics, per
  // (partner face unknown, axis and side of the ghost): it depends only on
  // the geometry and the diffusion solution, which the solve does not
  // change, but each ghost enters several stencils in every assembly, rate
  // and write-back pass.  The first evaluation is kept (state 0 empty,
  // 1 being written, 2 ready); the parallel passes read it without locks.
  mutable std::vector<T> macGhostV_;
  std::unique_ptr<std::atomic<std::uint8_t>[]> macGhostState_;
  std::size_t macGhostSlots_ = 0;
  bool macGhostOn_ = false;

  // Reflect a grid point into the solve box (reflective boundary through the
  // outermost grid points).
  IndexType macReflect(IndexType a) const {
    for (unsigned d = 0; d < D; ++d)
      while (a[d] < minIndex[d] || a[d] > maxIndex[d])
        a[d] = (a[d] < minIndex[d] ? 2 * minIndex[d] : 2 * maxIndex[d]) - a[d];
    return a;
  }

  T macPhiR(const IndexType &a) const {
    const IndexType r = macReflect(a);
    const auto it = macPhiR_.find(linearIndex(r));
    return it == macPhiR_.end() ? T(1) : it->second;
  }

  // phi_R (positive in the oxide) at a point x in index units, D-linear
  // over the grid points around it.
  T macPhiRAt(const std::array<T, D> &x) const {
    std::array<T, D> base{}, frac{};
    for (unsigned d = 0; d < D; ++d) {
      base[d] = std::floor(x[d]);
      frac[d] = x[d] - base[d];
    }
    T value = T(0);
    for (unsigned corner = 0; corner < (1u << D); ++corner) {
      IndexType idx;
      T w = T(1);
      for (unsigned d = 0; d < D; ++d) {
        const bool up = (corner >> d) & 1u;
        idx[d] = static_cast<viennahrle::IndexType>(base[d]) + (up ? 1 : 0);
        w *= up ? frac[d] : (T(1) - frac[d]);
      }
      if (w > T(0))
        value += w * std::clamp(macPhiR(idx), T(-1), T(1));
    }
    return value;
  }

  std::uint8_t macKind(const IndexType &a) const {
    const IndexType r = macReflect(a);
    const std::size_t j = nodeLookupFlat[linearIndex(r)];
    if (j != noNode)
      return nodes[j].nitride ? macNitride : macOxide;
    return macPhiR(r) < T(0) ? macSi : macGas;
  }

  std::size_t macNode(const IndexType &a) const {
    if (!inBounds(a))
      return noNode;
    return nodeLookupFlat[linearIndex(a)];
  }

  // Unknown id of the in-box face (d, a) or -1.
  std::int64_t macFaceId(unsigned d, const IndexType &a) const {
    const std::size_t ja = macNode(a);
    if (ja != noNode)
      return macPlusFace_[d][ja];
    IndexType b = a;
    b[d] += 1;
    const std::size_t jb = macNode(b);
    if (jb != noNode)
      return macMinusFace_[d][jb];
    return -1;
  }

  // Map face (d, a) into the box: returns the in-box face and the sign of
  // the velocity (odd in the normal direction of each reflection).
  T macFaceReflect(unsigned d, IndexType &a) const {
    T sign = T(1);
    for (unsigned e = 0; e < D; ++e) {
      if (e == d) {
        while (a[e] < minIndex[e] || a[e] > maxIndex[e] - 1) {
          if (a[e] < minIndex[e])
            a[e] = 2 * minIndex[e] - a[e] - 1;
          else
            a[e] = 2 * maxIndex[e] - a[e] - 1;
          sign = -sign;
        }
      } else {
        while (a[e] < minIndex[e] || a[e] > maxIndex[e]) {
          if (a[e] < minIndex[e])
            a[e] = 2 * minIndex[e] - a[e];
          else
            a[e] = 2 * maxIndex[e] - a[e];
        }
      }
    }
    return sign;
  }

  T macPhiFace(unsigned d, IndexType a) const {
    IndexType b = a;
    b[d] += 1;
    return T(0.5) * (macPhiR(a) + macPhiR(b));
  }

  // Front velocity component d for a Si ghost: the expansion velocity at the
  // front crossing between the in-box face (d, a) and the ghost face (d,
  // ghost), a fraction theta of the way (the oxidation front's own interface
  // point, see OxidationDiffusion::frontVelocityAt); zero for a nitride face
  // (weld).  Falls back to the mean reaction velocity of the face's oxide
  // nodes where no crossing is near.
  T macFrontVelocity(unsigned d, const IndexType &a, const IndexType &ghost,
                     T theta) const {
    IndexType b = a;
    b[d] += 1;
    T sum = T(0);
    int count = 0;
    bool nitride = false;
    for (const auto &q : {a, b}) {
      const std::size_t j = macNode(q);
      if (j == noNode)
        continue;
      if (nodes[j].nitride) {
        nitride = true;
        continue;
      }
      const auto v = reactionVelocityAt(j);
      if (v[0] != T(0) || v[1] != T(0) || (D == 3 && v[2] != T(0))) {
        sum += v[d];
        ++count;
      }
    }
    if (boundaryCacheValid_ && !(nitride && count == 0)) {
      Vec3D<T> p{T(0), T(0), T(0)};
      for (unsigned e = 0; e < D; ++e)
        p[e] = (static_cast<T>(a[e]) + (e == d ? T(0.5) : T(0)) +
                theta * static_cast<T>(ghost[e] - a[e])) *
               gridDelta;
      const auto front = diffusionField->frontVelocityAt(p);
      if (front.found)
        return static_cast<T>(reactionSign) *
               (oxidationParameters.expansionCoefficient - T(1)) * front.speed *
               front.normal[d];
    }
    if (count > 0)
      return sum / static_cast<T>(count);
    if (nitride)
      return T(0);
    for (const auto &q : {a, b}) {
      const std::size_t j = macNode(q);
      if (j != noNode)
        return reactionBoundaryVelocity(nodes[j].index)[d];
    }
    return T(0);
  }

  // Velocity component d at face (d, a) as an affine expression of the
  // unknowns.  partner: the adjacent face (same d, one step along another
  // or the same direction) that a ghost is built from.
  MacLin macFace(unsigned d, IndexType a, const IndexType *partner) const {
    const IndexType unreflected = a;
    const T sign = macFaceReflect(d, a);
    MacLin r;
    const std::int64_t id = macFaceId(d, a);
    if (id >= 0) {
      r.t.emplace_back(id, sign);
      return r;
    }
    // not an unknown: Si ghost, gas ghost, or nothing
    const T phiF = macPhiFace(d, a);
    if (partner == nullptr)
      return r; // zero
    IndexType pa = *partner;
    const T psign = macFaceReflect(d, pa);
    bool inBox = true;
    for (unsigned e = 0; e < D; ++e)
      inBox = inBox && a[e] == unreflected[e] && pa[e] == (*partner)[e];
    const std::int64_t pid = macFaceId(d, pa);
    MacLin pe;
    if (pid >= 0)
      pe.t.emplace_back(pid, psign);
    if (phiF >= T(0)) // gas face: zero gradient from the partner
      return pe;
    const T phiP = macPhiFace(d, pa);
    T theta = T(1);
    if (phiP > T(0) && phiP - phiF > T(0))
      theta = phiP / (phiP - phiF);
    const T thMin = std::clamp(
        T(deformationParameters.minMechanicsBoundaryDistance), T(0), T(1));
    theta = std::clamp(theta, thMin, T(1));
    // ghost = (1 - 1/theta) partner + v_front / theta (true orientation)
    const T vf =
        macGhostFrontVelocity(d, pa, a, theta, inBox ? pid : -1) * psign;
    MacLin g = pe.scaled(T(1) - T(1) / theta);
    g.c += vf / theta;
    return g;
  }

  // macFrontVelocity(d, pa, a, theta) of a ghost face a next to the
  // unknown pid = (d, pa), kept during solveMechanics (see macGhostV_).
  // pid < 0 (a reflected face or no unknown partner): not kept.
  T macGhostFrontVelocity(unsigned d, const IndexType &pa, const IndexType &a,
                          T theta, std::int64_t pid) const {
    std::size_t slot = macGhostSlots_;
    if (macGhostOn_ && pid >= 0) {
      int axis = -1, step = 0;
      for (unsigned e = 0; e < D; ++e) {
        const auto delta = a[e] - pa[e];
        if (delta == 0)
          continue;
        if (axis != -1 || (delta != 1 && delta != -1)) {
          axis = -2;
          break;
        }
        axis = static_cast<int>(e);
        step = static_cast<int>(delta);
      }
      if (axis >= 0)
        slot = static_cast<std::size_t>(pid) * 2 * D +
               2 * static_cast<std::size_t>(axis) + (step > 0 ? 1 : 0);
    }
    if (slot >= macGhostSlots_)
      return macFrontVelocity(d, pa, a, theta);
    auto &state = macGhostState_[slot];
    if (state.load(std::memory_order_acquire) == 2)
      return macGhostV_[slot];
    const T v = macFrontVelocity(d, pa, a, theta);
    std::uint8_t empty = 0;
    if (state.compare_exchange_strong(empty, 1, std::memory_order_acq_rel)) {
      macGhostV_[slot] = v;
      state.store(2, std::memory_order_release);
    }
    return v;
  }

  // Empty ghost front velocities for this solve's node set and boundary
  // values; on until the solve ends (macGhostOn_).
  void macGhostReset() {
    const std::size_t slots = static_cast<std::size_t>(macNVel_) * 2 * D;
    if (slots != macGhostSlots_) {
      macGhostV_.assign(slots, T(0));
      macGhostState_ = std::make_unique<std::atomic<std::uint8_t>[]>(slots);
      macGhostSlots_ = slots;
    }
    for (std::size_t k = 0; k < slots; ++k)
      macGhostState_[k].store(0, std::memory_order_relaxed);
    macGhostOn_ = true;
  }

  // dv_d/dx_d at grid point q (faces (d, q - e_d), (d, q)).
  MacLin macNormalRate(unsigned d, const IndexType &q) const {
    IndexType lo = q;
    lo[d] -= 1;
    const IndexType hi = q;
    MacLin r = macFace(d, hi, &lo);
    r.add(macFace(d, lo, &hi), T(-1));
    return r.scaled(T(1) / gridDelta);
  }

  // dv_c/dx_e + dv_e/dx_c at corner (ce, a).
  MacLin macShearRate(unsigned c, unsigned e, const IndexType &a) const {
    IndexType ae = a, ac = a;
    ae[e] += 1;
    ac[c] += 1;
    MacLin r = macFace(c, ae, &a);
    r.add(macFace(c, a, &ae), T(-1));
    r.add(macFace(e, ac, &a));
    r.add(macFace(e, a, &ac), T(-1));
    return r.scaled(T(1) / gridDelta);
  }

  std::int64_t macCornerLookup(unsigned pair, const IndexType &a) const {
    const auto &m = macCornerId_[pair];
    const auto it = m.find(linearIndex(a));
    return it == m.end() ? -1 : it->second;
  }

  // Corner (pair, a) mapped into the box: its id (or -1) and the sign of its
  // shear stress, odd under each reflection in c or e.
  std::int64_t macCornerReflected(unsigned pair, IndexType a, T &sign) const {
    unsigned c, e;
    macPairDirs(pair, c, e);
    sign = T(1);
    for (unsigned k : {c, e}) {
      while (a[k] < minIndex[k] || a[k] > maxIndex[k] - 1) {
        if (a[k] < minIndex[k])
          a[k] = 2 * minIndex[k] - a[k] - 1;
        else
          a[k] = 2 * maxIndex[k] - a[k] - 1;
        sign = -sign;
      }
    }
    return macCornerLookup(pair, a);
  }

  // Inside-positive level set of the oxide (ambient) or the nitride (mask)
  // at a grid point; NaN where it is undefined or not cached.
  T macPsi(bool nitride, const IndexType &a) const {
    const auto &m = nitride ? macPsiM_ : macPsiA_;
    const auto it = m.find(linearIndex(macReflect(a)));
    return it == m.end() ? std::numeric_limits<T>::quiet_NaN() : it->second;
  }

  // d psi / d x_e (per cell) at x: the mean of the one-sided differences
  // whose two ends are defined; NaN if there is none.
  T macPsiSlope(bool nitride, const IndexType &x, unsigned e) const {
    const T c = macPsi(nitride, x);
    T sum = T(0);
    int cnt = 0;
    if (std::isfinite(c))
      for (int sg : {-1, 1}) {
        IndexType y = x;
        y[e] += sg;
        const T v = macPsi(nitride, y);
        if (std::isfinite(v)) {
          sum += (v - c) * static_cast<T>(sg);
          ++cnt;
        }
      }
    return cnt ? sum / static_cast<T>(cnt)
               : std::numeric_limits<T>::quiet_NaN();
  }

  static std::array<T, 3> macUnitOutward(const std::array<T, 3> &grad) {
    T norm = T(0);
    for (unsigned e = 0; e < D; ++e)
      norm += grad[e] * grad[e];
    std::array<T, 3> n{T(0), T(0), T(0)};
    if (norm <= T(1e-12))
      return n; // degenerate
    norm = std::sqrt(norm);
    for (unsigned e = 0; e < D; ++e)
      n[e] = -grad[e] / norm; // psi grows into the material
    return n;
  }

  // Outward free-surface normal at the crossing between solid point s and
  // gas point q (their midpoint), from the material's level set; the axis
  // normal if the gradient is degenerate.
  std::array<T, 3> macSurfaceNormal(bool nitride, const IndexType &s,
                                    const IndexType &q) const {
    unsigned dq = 0;
    for (unsigned e = 0; e < D; ++e)
      if (q[e] != s[e])
        dq = e;
    const T step = static_cast<T>(q[dq] - s[dq]);
    std::array<T, 3> g{T(0), T(0), T(0)};
    g[dq] = (macPsi(nitride, q) - macPsi(nitride, s)) * step;
    if (!std::isfinite(g[dq])) {
      // the level set is undefined at the gas point: the surface passes
      // through the node, crossed along d
      g[dq] = step;
      return g;
    }
    for (unsigned e = 0; e < D; ++e) {
      if (e == dq)
        continue;
      T sum = T(0);
      int cnt = 0;
      for (const auto &x : {s, q}) {
        const T slope = macPsiSlope(nitride, x, e);
        if (std::isfinite(slope)) {
          sum += slope;
          ++cnt;
        }
      }
      g[e] = cnt ? sum / static_cast<T>(cnt) : T(0);
    }
    auto n = macUnitOutward(g);
    if (n[0] == T(0) && n[1] == T(0) && n[2] == T(0))
      n[dq] = step;
    return n;
  }

  // Outward free-surface normal at corner (ce, a): gradient over its points.
  std::array<T, 3> macCornerNormal(bool nitride, unsigned c, unsigned e,
                                   const IndexType &a) const {
    // edge differences of the corner square whose two ends are defined
    bool defined = true;
    auto edges = [&](unsigned along, unsigned across) {
      T sum = T(0);
      int cnt = 0;
      for (int b = 0; b <= 1; ++b) {
        IndexType p = a;
        p[across] += b;
        IndexType q = p;
        q[along] += 1;
        const T diff = macPsi(nitride, q) - macPsi(nitride, p);
        if (std::isfinite(diff)) {
          sum += diff;
          ++cnt;
        }
      }
      defined = defined && cnt > 0;
      return cnt ? sum / static_cast<T>(cnt) : T(0);
    };
    std::array<T, 3> g{T(0), T(0), T(0)};
    g[c] = edges(c, e);
    g[e] = edges(e, c);
    if (!defined) // too far from the band: no normal (zero shear)
      return {T(0), T(0), T(0)};
    if constexpr (D == 3) {
      const unsigned k = 3 - c - e;
      T sum = T(0);
      int cnt = 0;
      for (int bc = 0; bc <= 1; ++bc)
        for (int be = 0; be <= 1; ++be) {
          IndexType p = a;
          p[c] += bc;
          p[e] += be;
          const T slope = macPsiSlope(nitride, p, k);
          if (std::isfinite(slope)) {
            sum += slope;
            ++cnt;
          }
        }
      g[k] = cnt ? sum / static_cast<T>(cnt) : T(0);
    }
    return macUnitOutward(g);
  }

  // (P S P)_ab with P = I - n n^T: the part of the stress a free surface
  // with normal n can carry.
  static MacLin macProjected(const std::array<MacLin, 9> &S,
                             const std::array<T, 3> &n, unsigned a,
                             unsigned b) {
    auto P = [&](unsigned i, unsigned j) {
      return (i == j ? T(1) : T(0)) - n[i] * n[j];
    };
    MacLin r;
    for (unsigned c = 0; c < D; ++c)
      for (unsigned e = 0; e < D; ++e) {
        const T w = P(a, c) * P(e, b);
        if (w != T(0))
          r.add(S[tensorIndex(c, e)], w);
      }
    return r;
  }

  // Cauchy stress at solid grid point q, scaled by 1/eta_ref: normal
  // components from the material law, shear the mean of the surrounding
  // corners that are not on the surface.
  std::array<MacLin, 9> macPointStress(const IndexType &q, T etaRef) const {
    std::array<MacLin, 9> S{};
    for (unsigned c = 0; c < D; ++c)
      S[tensorIndex(c, c)] = macNormalStress(c, q, q, etaRef);
    for (unsigned pair = 0; pair < macPairs; ++pair) {
      unsigned c, e;
      macPairDirs(pair, c, e);
      MacLin sum;
      int cnt = 0;
      for (int dc = -1; dc <= 0; ++dc)
        for (int de = -1; de <= 0; ++de) {
          IndexType lo = q;
          lo[c] += dc;
          lo[e] += de;
          T sign;
          const std::int64_t k = macCornerReflected(pair, lo, sign);
          if (k < 0 || macCornerGas_[static_cast<std::size_t>(k)])
            continue;
          sum.add(macShearStress(c, e, lo, etaRef));
          ++cnt;
        }
      if (cnt > 0)
        S[tensorIndex(c, e)] = S[tensorIndex(e, c)] =
            sum.scaled(T(1) / static_cast<T>(cnt));
    }
    return S;
  }

  // Shear stress sigma_ce at corner (ce, a), scaled by 1/eta_ref.  A corner
  // on the surface (>= 2 gas points) carries the traction-free shear of its
  // solid points plus the residual modulus.
  MacLin macShearStress(unsigned c, unsigned e, const IndexType &a,
                        T etaRef) const {
    T sign;
    const std::int64_t k = macCornerReflected(macPairIndex(c, e), a, sign);
    if (k < 0)
      return MacLin{};
    const std::size_t kk = static_cast<std::size_t>(k);
    MacLin r = macShearRate(c, e, a).scaled(macCornerMod_[kk]);
    if (!macCornerGas_[kk]) {
      r.c += sign * macCornerMem_[kk];
      return r;
    }
    const auto n = macCornerNormal(!macCornerOxide_[kk], c, e, a);
    if (n[0] == T(0) && n[1] == T(0) && n[2] == T(0))
      return r;
    MacLin proj;
    int cnt = 0;
    for (int bc = 0; bc <= 1; ++bc)
      for (int be = 0; be <= 1; ++be) {
        IndexType q = a;
        q[c] += bc;
        q[e] += be;
        const auto kind = macKind(q);
        if (kind != macOxide && kind != macNitride)
          continue;
        proj.add(macProjected(macPointStress(q, etaRef), n, c, e));
        ++cnt;
      }
    if (cnt > 0)
      r.add(proj, T(1) / static_cast<T>(cnt));
    return r;
  }

  // Normal stress sigma_dd at grid point q seen from face f whose other end
  // is s, scaled by 1/eta_ref.
  MacLin macNormalStress(unsigned d, const IndexType &q, const IndexType &s,
                         T etaRef) const {
    const std::uint8_t kq = macKind(q);
    const std::size_t js = macNode(s);
    auto oxideLike = [&](std::size_t jp, const IndexType &at) {
      // -p + 2 eta_impl dv_d/dx_d + memory, properties of node jp
      MacLin r = macNormalRate(d, at).scaled(T(2) * nodeImplicitViscosity(jp) /
                                             etaRef);
      r.c += macNodeMemory(jp, d) / etaRef;
      if (macPid_[jp] >= 0)
        r.t.emplace_back(macPid_[jp], T(-1));
      return r;
    };
    auto nitrideLike = [&](std::size_t jp, const IndexType &at) {
      // 2 eta_impl dv_d/dx_d + (K dt - 2 eta_impl / 3) div v + memory
      T twoEta, volumetric;
      nitrideModuli(jp, twoEta, volumetric);
      MacLin r = macNormalRate(d, at).scaled(twoEta / etaRef);
      for (unsigned c = 0; c < D; ++c)
        r.add(macNormalRate(c, at), volumetric / etaRef);
      r.c += macNodeMemory(jp, d) / etaRef;
      return r;
    };
    if (kq == macOxide || kq == macNitride) {
      const IndexType qr = macReflect(q);
      const std::size_t jq = nodeLookupFlat[linearIndex(qr)];
      return kq == macOxide ? oxideLike(jq, q) : nitrideLike(jq, q);
    }
    if (js == noNode)
      return MacLin{};
    if (kq == macGas) {
      // extrapolated from s to its traction-free value (P sigma_s P)_dd at
      // the crossing: a free surface carries only tangential stress
      unsigned fi = 0;
      for (unsigned e = 0; e < D; ++e)
        if (q[e] != s[e])
          fi = e * 2u + (q[e] > s[e] ? 1u : 0u);
      const std::size_t n = nodes.size();
      T theta = faceBCDists_[fi * n + js] / gridDelta;
      if (faceBCTypes_[fi * n + js] != Boundary::AMBIENT)
        theta = T(1);
      theta = std::clamp(theta, T(1e-3), T(1));
      const auto S = macPointStress(s, etaRef);
      const auto normal = macSurfaceNormal(nodes[js].nitride, s, q);
      MacLin g = S[tensorIndex(d, d)].scaled(T(1) - T(1) / theta);
      g.add(macProjected(S, normal, d, d), T(1) / theta);
      return g;
    }
    // Si point: the node's material law with the ghost velocities
    return nodes[js].nitride ? nitrideLike(js, q) : oxideLike(js, q);
  }

  // Normal memory stress (component d) of node j, scaled: oxide decay S_start,
  // nitride sigma_start (Maxwell: decayed deviator, elastic mean stress).
  T macNodeMemory(std::size_t j, unsigned d) const {
    if (nodes[j].nitride) {
      const auto &start = nitrideStartOf(j);
      if (!viscousNitride())
        return start[tensorIndex(d, d)];
      const T decay = nitrideDecayFor(nitrideNodeViscosity(j));
      const T mean = (start[0] + start[4] + start[8]) / T(3);
      return decay * start[tensorIndex(d, d)] + (T(1) - decay) * mean;
    }
    return nodeDecay(j) * startStress_[j][tensorIndex(d, d)];
  }

  // Number faces, pressures and corners for the current node set.
  void macBuild() {
    const std::size_t n = nodes.size();
    for (unsigned d = 0; d < D; ++d) {
      macPlusFace_[d].assign(n, -1);
      macMinusFace_[d].assign(n, -1);
    }
    // phi_R, psi_A, psi_M on the nodes and their neighbours (two steps)
    macPhiR_.clear();
    macPsiA_.clear();
    macPsiM_.clear();
    {
      ConstSparseIterator reactionIt(reactionInterface->getDomain());
      ConstSparseIterator ambientIt(ambientInterface->getDomain());
      ConstSparseIterator maskIt = makeMaskIterator();
      auto put = [&](const IndexType &q) {
        if (!inBounds(q))
          return;
        const std::size_t li = linearIndex(q);
        if (macPhiR_.count(li))
          return;
        macPhiR_[li] = static_cast<T>(reactionSign) * valueAt(reactionIt, q);
        auto defined = [&](ConstSparseIterator &it) {
          it.goToIndices(q);
          return it.isDefined() ? it.getValue()
                                : std::numeric_limits<T>::quiet_NaN();
        };
        macPsiA_[li] = static_cast<T>(ambientSign) * defined(ambientIt);
        if (maskInterface != nullptr)
          macPsiM_[li] = static_cast<T>(maskSign) * defined(maskIt);
      };
      for (std::size_t i = 0; i < n; ++i) {
        unsigned span = 1;
        for (unsigned d = 0; d < D; ++d)
          span *= 5;
        for (unsigned code = 0; code < span; ++code) {
          IndexType q = nodes[i].index;
          unsigned rest = code;
          for (unsigned d = 0; d < D; ++d) {
            q[d] += static_cast<int>(rest % 5) - 2;
            rest /= 5;
          }
          put(q);
        }
      }
    }
    // faces: unknown if a solid end point and the midpoint on the oxide side
    std::int64_t count = 0;
    for (unsigned d = 0; d < D; ++d) {
      for (std::size_t i = 0; i < n; ++i) {
        for (int side = 0; side < 2; ++side) {
          IndexType a = nodes[i].index;
          if (side == 0)
            a[d] -= 1; // face (d, node - e_d)
          IndexType b = a;
          b[d] += 1;
          if (a[d] < minIndex[d] || b[d] > maxIndex[d])
            continue; // mirror face
          auto &slot = side == 0 ? macMinusFace_[d][i] : macPlusFace_[d][i];
          if (slot >= 0)
            continue;
          // shared with the other end point if that is a node
          const std::size_t other = macNode(side == 0 ? a : b);
          std::int64_t existing = -1;
          if (other != noNode)
            existing =
                side == 0 ? macPlusFace_[d][other] : macMinusFace_[d][other];
          if (existing >= 0) {
            slot = existing;
            continue;
          }
          if (macPhiFace(d, a) < T(0))
            continue; // midpoint in Si
          slot = count++;
          if (other != noNode) {
            if (side == 0)
              macPlusFace_[d][other] = slot;
            else
              macMinusFace_[d][other] = slot;
          }
        }
      }
    }
    macNVel_ = count;
    macPid_.assign(n, -1);
    bool gasContact = false;
    bool hasNitride = false;
    for (std::size_t i = 0; i < n; ++i) {
      if (nodes[i].nitride) {
        hasNitride = true;
        continue;
      }
      macPid_[i] = count++;
      for (unsigned fi = 0; fi < 2 * D; ++fi)
        gasContact =
            gasContact || faceBCTypes_[fi * n + i] == Boundary::AMBIENT;
    }
    macN_ = count;
    macPinned_ = !gasContact && !hasNitride; // closed oxide: pressure gauge
    // corners: every corner with a solid point
    macCornerId_.assign(macPairs, {});
    macCorners_.clear();
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned pair = 0; pair < macPairs; ++pair) {
        unsigned c, e;
        macPairDirs(pair, c, e);
        for (int dc = -1; dc <= 0; ++dc)
          for (int de = -1; de <= 0; ++de) {
            IndexType lo = nodes[i].index;
            lo[c] += dc;
            lo[e] += de;
            if (lo[c] < minIndex[c] || lo[c] > maxIndex[c] - 1 ||
                lo[e] < minIndex[e] || lo[e] > maxIndex[e] - 1)
              continue;
            auto &m = macCornerId_[pair];
            const std::size_t li = linearIndex(lo);
            if (m.count(li))
              continue;
            m[li] = static_cast<std::int64_t>(macCorners_.size());
            macCorners_.emplace_back(pair, lo);
          }
      }
    const std::size_t nk = macCorners_.size();
    macCornerGas_.assign(nk, 0);
    macCornerOxide_.assign(nk, 0);
    macCornerStartV_.assign(nk, T(0));
    macCornerEndV_.assign(nk, T(0));
    macCornerMod_.assign(nk, T(0));
    macCornerMem_.assign(nk, T(0));
    macCornerEta_.assign(nk, T(deformationParameters.viscosity));
    macCornerEtaN_.assign(nk, T(deformationParameters.nitrideViscosity));
    for (std::size_t k = 0; k < nk; ++k) {
      unsigned c, e;
      macPairDirs(macCorners_[k].first, c, e);
      const IndexType &lo = macCorners_[k].second;
      int oxide = 0, nitride = 0, gas = 0;
      for (int bc = 0; bc <= 1; ++bc)
        for (int be = 0; be <= 1; ++be) {
          IndexType q = lo;
          q[c] += bc;
          q[e] += be;
          const auto kind = macKind(q);
          gas += kind == macGas;
          oxide += kind == macOxide;
          nitride += kind == macNitride;
        }
      macCornerGas_[k] = gas >= 2;
      macCornerOxide_[k] = oxide >= nitride && oxide > 0;
    }
    macX_.assign(static_cast<std::size_t>(macN_), T(0));
    macNodesVersion_ = nodesVersion_;
  }
  bool macPinned_ = false;
  bool couplingIterate_ = false;
  // Node set version of the last converged viscosities (warm start).
  std::size_t macWarmVersion_ = std::numeric_limits<std::size_t>::max();

  // Corner start stresses: the carried history, else the mean of the
  // surrounding nodes' stored shear.
  void macLoadCornerStart() {
    for (std::size_t k = 0; k < macCorners_.size(); ++k) {
      const unsigned pair = macCorners_[k].first;
      unsigned c, e;
      macPairDirs(pair, c, e);
      const IndexType &lo = macCorners_[k].second;
      const auto it = macCornerStart_.find(lo);
      if (it != macCornerStart_.end() && macHasHistory(it->second[pair])) {
        macCornerStartV_[k] = it->second[pair];
        continue;
      }
      T sum = T(0);
      int cnt = 0;
      for (int bc = 0; bc <= 1; ++bc)
        for (int be = 0; be <= 1; ++be) {
          IndexType q = lo;
          q[c] += bc;
          q[e] += be;
          const std::size_t j = macNode(q);
          if (j == noNode)
            continue;
          sum += nodes[j].nitride ? nitrideStartOf(j)[tensorIndex(c, e)]
                                  : startStress_[j][tensorIndex(c, e)];
          ++cnt;
        }
      macCornerStartV_[k] = cnt ? sum / static_cast<T>(cnt) : T(0);
    }
  }

  T macCornerDecay(std::size_t k) const {
    const T tau0 = effectiveStressRelaxationTime();
    if (tau0 <= std::numeric_limits<T>::epsilon())
      return T(0);
    const T tau = tau0 * macCornerEta_[k] / T(deformationParameters.viscosity);
    return std::exp(-T(deformationParameters.stressTimeStep) / tau);
  }

  // Corner modulus and memory (scaled by 1/eta_ref) for the current viscosity.
  void macCornerCoefficients(T etaRef) {
#pragma omp parallel for schedule(static)
    for (std::size_t k = 0; k < macCorners_.size(); ++k) {
      unsigned c, e;
      macPairDirs(macCorners_[k].first, c, e);
      const IndexType &lo = macCorners_[k].second;
      const T decay = macCornerDecay(k);
      const T etaImpl = std::max((T(1) - decay) * macCornerEta_[k],
                                 std::numeric_limits<T>::min() * T(1e10));
      const T etaImplN = nitrideImplicitFor(macCornerEtaN_[k]);
      T inv = T(0);
      int cnt = 0;
      for (int bc = 0; bc <= 1; ++bc)
        for (int be = 0; be <= 1; ++be) {
          IndexType q = lo;
          q[c] += bc;
          q[e] += be;
          const auto kind = macKind(q);
          if (kind == macOxide) {
            inv += T(1) / etaImpl;
            ++cnt;
          } else if (kind == macNitride) {
            inv += T(1) / etaImplN;
            ++cnt;
          }
        }
      macCornerMod_[k] = cnt ? static_cast<T>(cnt) / inv / etaRef : T(0);
      if (macCornerGas_[k]) { // on the surface: residual modulus, no memory
        macCornerMod_[k] *= T(1e-3);
        macCornerMem_[k] = T(0);
        continue;
      }
      macCornerMem_[k] =
          (macCornerOxide_[k] ? decay : nitrideDecayFor(macCornerEtaN_[k])) *
          macCornerStartV_[k] / etaRef;
    }
  }

  // Assemble and solve the linear MAC system for the current viscosity.
  bool macSolveLinear(T etaRef, T &backward) {
    macCornerCoefficients(etaRef);
    const std::size_t n = nodes.size();
    const std::int64_t N = macN_;
    std::vector<MacLin> rows(static_cast<std::size_t>(N));
    // momentum rows: one per face unknown (visit from the node owning it)
#pragma omp parallel for schedule(dynamic, 64)
    for (std::size_t i = 0; i < n; ++i) {
      for (unsigned d = 0; d < D; ++d)
        for (int side = 0; side < 2; ++side) {
          const std::int64_t id =
              side == 0 ? macMinusFace_[d][i] : macPlusFace_[d][i];
          if (id < 0)
            continue;
          IndexType a = nodes[i].index;
          if (side == 0)
            a[d] -= 1;
          // the face's owner: the node at a if it is one, else this node
          const std::size_t ja = macNode(a);
          if (side == 0 && ja != noNode)
            continue; // visited from a
          IndexType b = a;
          b[d] += 1;
          MacLin r = macNormalStress(d, b, a, etaRef);
          r.add(macNormalStress(d, a, b, etaRef), T(-1));
          for (unsigned e = 0; e < D; ++e) {
            if (e == d)
              continue;
            IndexType lo = a;
            lo[e] -= 1;
            r.add(macShearStress(d, e, a, etaRef));
            r.add(macShearStress(d, e, lo, etaRef), T(-1));
          }
          rows[static_cast<std::size_t>(id)] = r.scaled(T(1) / gridDelta);
        }
    }
    // continuity rows
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < n; ++i) {
      const std::int64_t pid = macPid_[i];
      if (pid < 0)
        continue;
      MacLin r;
      for (unsigned d = 0; d < D; ++d)
        r.add(macNormalRate(d, nodes[i].index));
      rows[static_cast<std::size_t>(pid)] = std::move(r);
    }
    if (macPinned_) { // closed oxide: gauge at the first oxide node
      for (std::size_t i = 0; i < n; ++i)
        if (macPid_[i] >= 0) {
          MacLin g;
          g.t.emplace_back(macPid_[i], T(1));
          g.c = -T(deformationParameters.ambientPressure) / etaRef;
          rows[static_cast<std::size_t>(macPid_[i])] = g;
          break;
        }
    }
    if constexpr (D == 3)
      return macAmgSolve(rows, backward); // consumes the rows
    else
      return macDirectSolve(rows, backward);
  }

  // AMGCL: FGMRES, block-triangular velocity-pressure preconditioner (one
  // smoothed-aggregation V-cycle with ILU0 for the velocities, SPAI0 on
  // B diag(A)^-1 G for the pressure), row/column-equilibrated system,
  // warm-started from the previous iterate. The preconditioner works in
  // single precision (half the memory traffic of its V-cycles); FGMRES and
  // its residual stay double, so the solve meets the same tolerance.
  typedef amgcl::backend::builtin<double> MacAmgBackend;
  typedef amgcl::backend::builtin<float> MacAmgPrecondBackend;
  typedef amgcl::make_solver<
      amgcl::amg<MacAmgPrecondBackend, amgcl::coarsening::smoothed_aggregation,
                 amgcl::relaxation::ilu0>,
      amgcl::solver::preonly<MacAmgPrecondBackend>>
      MacAmgUSolver;
  typedef amgcl::make_solver<
      amgcl::relaxation::as_preconditioner<MacAmgPrecondBackend,
                                           amgcl::relaxation::spai0>,
      amgcl::solver::preonly<MacAmgPrecondBackend>>
      MacAmgPSolver;
  typedef amgcl::preconditioner::schur_pressure_correction<MacAmgUSolver,
                                                           MacAmgPSolver>
      MacAmgPrecond;
  typedef amgcl::solver::fgmres<MacAmgBackend> MacAmgKrylov;

  // The 3-D solve runs on the GPU from this many unknowns (with the GPU
  // mode). Per solve on a 3-D LOCOS corner, CPU (10 threads) / GPU: 15.6k
  // unknowns 0.11 / 0.11 s, 59k 0.47 / 0.17 s, 360k 3.9 / 0.40 s; the
  // setup stays on the host.
  static constexpr std::size_t macGpuMinUnknowns = 20000;
  GpuMode gpuMode_ = GpuMode::Cpu;
  bool macGpuOff_ = false; // Auto mode after the GPU was unavailable
  void *macGpu_ = nullptr; // GPU setup (ViennaLS_GPU), owned

  // The preconditioner of the last setup and the equilibration it was built
  // for. The next solve on the same node set (a stress-coupling iteration,
  // the converged solve, a CFL re-solve) changes only the viscosities and
  // the right-hand side: it reuses them, which saves the setup (~15 % of a
  // 3-D LOCOS solve) at the same iterations. A reused solve may take twice
  // the iterations of the solve after the setup (+20); beyond that, or when
  // it fails, a fresh setup solves it.
  std::unique_ptr<MacAmgPrecond> macAmgP_;
  std::unique_ptr<MacAmgKrylov> macAmgS_;
  std::size_t macAmgVersion_ = std::numeric_limits<std::size_t>::max();
  std::vector<double> macAmgRs_, macAmgCs_;
  std::size_t macAmgSetupIters_ = 0;
  // The CSR system and its equilibrated copy, kept between solves: on one
  // node set their sizes do not change, so a solve neither reallocates nor
  // zero-fills ~10^7 entries.
  std::vector<ptrdiff_t> macCsrPtr_, macCsrCol_;
  std::vector<double> macCsrVal_, macCsrRhs_, macSval_, macSrhs_, macY_;

  bool macAmgSolve(std::vector<MacLin> &rows, T &backward) {
    // CSR, per row the columns sorted and duplicates summed, in two parallel
    // passes: sort each row in place and count its entries, then sum them
    // again into place (and free the row) -- the matrix of a serial build.
    const std::size_t N = rows.size();
    auto &ptr = macCsrPtr_, &col = macCsrCol_;
    auto &val = macCsrVal_, &rhs = macCsrRhs_;
    ptr.resize(N + 1);
    rhs.resize(N);
    ptr[0] = 0;
    const auto byColumn = [](const auto &a, const auto &b) {
      return a.first < b.first;
    };
#pragma omp parallel for schedule(dynamic, 256)
    for (std::size_t i = 0; i < N; ++i) {
      auto &t = rows[i].t;
      std::sort(t.begin(), t.end(), byColumn);
      ptrdiff_t count = 0;
      for (std::size_t q = 0; q < t.size();) {
        const std::int64_t c = t[q].first;
        double v = 0.;
        for (; q < t.size() && t[q].first == c; ++q)
          v += static_cast<double>(t[q].second);
        if (v != 0.)
          ++count;
      }
      ptr[i + 1] = count;
      rhs[i] = -static_cast<double>(rows[i].c);
    }
    for (std::size_t i = 0; i < N; ++i)
      ptr[i + 1] += ptr[i];
    col.resize(static_cast<std::size_t>(ptr[N]));
    val.resize(static_cast<std::size_t>(ptr[N]));
#pragma omp parallel for schedule(dynamic, 256)
    for (std::size_t i = 0; i < N; ++i) {
      auto &t = rows[i].t;
      auto k = static_cast<std::size_t>(ptr[i]);
      for (std::size_t q = 0; q < t.size();) {
        const std::int64_t c = t[q].first;
        double v = 0.;
        for (; q < t.size() && t[q].first == c; ++q)
          v += static_cast<double>(t[q].second);
        if (v != 0.) {
          col[k] = static_cast<ptrdiff_t>(c);
          val[k++] = v;
        }
      }
      std::vector<std::pair<std::int64_t, T>>().swap(t);
    }
    if ((macAmgP_ || macGpu_) && macAmgVersion_ == macNodesVersion_ &&
        macAmgRs_.size() == N &&
        macAmgSolveWith(ptr, col, val, rhs, false, backward))
      return true;
    return macAmgSolveWith(ptr, col, val, rhs, true, backward);
  }

  // FGMRES target (relative residual of the equilibrated system): 1e-9 for
  // the velocities the interfaces move with; 1e-7 for the solves inside the
  // stress-coupling iteration, whose pressure and stress only feed the
  // coupling test (couplingTolerance, 2e-2..5e-2). 1e-7 is 10x below the
  // acceptance (stokesTolerance) and saves 40 % of their iterations; on a
  // 3-D LOCOS corner it changed the result ~400x less than tightening the
  // coupling test from 2e-2 to 2e-3 does.
  static constexpr double macTolerance = 1e-9;
  static constexpr double macCouplingTolerance = 1e-7;

  // One solve of the assembled system, with a fresh setup or with the kept
  // preconditioner and equilibration.
  bool macAmgSolveWith(const std::vector<ptrdiff_t> &ptr,
                       const std::vector<ptrdiff_t> &col,
                       const std::vector<double> &val,
                       const std::vector<double> &rhs, bool fresh,
                       T &backward) {
    const std::size_t N = rhs.size();
    if (fresh) {
      macAmgReset(); // the old hierarchy goes before the new one is built
      macAmgRs_.assign(N, 0.);
      macAmgCs_.assign(N, 0.);
      auto &rs = macAmgRs_, &cs = macAmgCs_;
      for (std::size_t i = 0; i < N; ++i)
        for (ptrdiff_t q = ptr[i]; q < ptr[i + 1]; ++q)
          rs[i] = std::max(rs[i], std::abs(val[q]));
      for (auto &v : rs)
        v = v > 0. ? 1. / v : 1.;
      for (std::size_t i = 0; i < N; ++i)
        for (ptrdiff_t q = ptr[i]; q < ptr[i + 1]; ++q)
          cs[col[q]] = std::max(cs[col[q]], std::abs(rs[i] * val[q]));
      for (auto &v : cs)
        v = v > 0. ? 1. / v : 1.;
    }
    const auto &rs = macAmgRs_, &cs = macAmgCs_;
    auto &sval = macSval_, &srhs = macSrhs_, &y = macY_;
    sval.resize(val.size());
    srhs.resize(N);
    y.resize(N);
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < N; ++i) {
      for (ptrdiff_t q = ptr[i]; q < ptr[i + 1]; ++q)
        sval[q] = rs[i] * val[q] * cs[col[q]];
      srhs[i] = rs[i] * rhs[i];
      y[i] = i < macX_.size() && std::isfinite(macX_[i])
                 ? static_cast<double>(macX_[i]) / cs[i]
                 : 0.;
    }
    const ptrdiff_t n = static_cast<ptrdiff_t>(N);
    auto A = std::tie(n, ptr, col, sval);
    std::size_t iters = 0;
    double err = 0.;
    const unsigned maxIterations =
        fresh ? 3000u : static_cast<unsigned>(2 * macAmgSetupIters_ + 20);
    const double target =
        couplingIterate_ ? macCouplingTolerance : macTolerance;
    if (fresh) {
      if (!macAmgSetup(ptr, col, sval, cs))
        return false;
      macAmgVersion_ = macNodesVersion_;
    }
#ifdef VIENNALS_GPU_BICGSTAB
    if (macGpu_) {
      static_assert(sizeof(ptrdiff_t) == sizeof(std::int64_t));
      unsigned gpuIters = 0;
      if (!gpu::gpuMacSolve(
              macGpu_, n, reinterpret_cast<const std::int64_t *>(ptr.data()),
              reinterpret_cast<const std::int64_t *>(col.data()), sval.data(),
              srhs.data(), y.data(), maxIterations, target, gpuIters, err)) {
        // a GPU that set up and then fails is a malfunction, not an absence
        VIENNACORE_LOG_ERROR(
            std::string(
                "OxidationDeformation: the GPU mechanics solve failed: ") +
            gpu::gpuGetLastErrorMessage());
        return false;
      }
      iters = gpuIters;
    } else
#endif
    {
      try {
        macAmgS_->prm.maxiter = maxIterations;
        macAmgS_->prm.tol = target;
        std::tie(iters, err) = (*macAmgS_)(A, *macAmgP_, srhs, y);
      } catch (const std::exception &e) {
        macAmgReset();
        VIENNACORE_LOG_WARNING(
            std::string("OxidationDeformation: AMGCL failed: ") + e.what());
        return false;
      }
    }
    if (fresh)
      macAmgSetupIters_ = iters;
    lastAmgIters_ = static_cast<unsigned>(iters);
    std::vector<double> x(N);
    for (std::size_t i = 0; i < N; ++i)
      x[i] = cs[i] * y[i];
    for (std::size_t i = 0; i < N; ++i)
      if (!std::isfinite(x[i]))
        return false;
    // the relative residual: a componentwise backward error is ~1 on rows
    // whose data is zero, however good the iterative solution
    backward = static_cast<T>(err);
    if (!(err < 1e-6)) {
      if (fresh)
        VIENNACORE_LOG_WARNING(
            "OxidationDeformation: AMGCL did not converge (iters=" +
            std::to_string(iters) + ", residual=" + std::to_string(err) + ")");
      return false;
    }
    macX_.assign(N, T(0));
    for (std::size_t i = 0; i < N; ++i)
      macX_[i] = static_cast<T>(x[i]);
    return true;
  }
  unsigned lastAmgIters_ = 0;

  void macAmgReset() {
    macAmgP_.reset();
    macAmgS_.reset();
#ifdef VIENNALS_GPU_BICGSTAB
    gpu::gpuMacFree(macGpu_);
#endif
    macGpu_ = nullptr;
  }

  // A fresh preconditioner for the equilibrated system: on the GPU with the
  // GPU mode and at least macGpuMinUnknowns unknowns, else on the CPU.
  bool macAmgSetup(const std::vector<ptrdiff_t> &ptr,
                   const std::vector<ptrdiff_t> &col,
                   const std::vector<double> &sval,
                   const std::vector<double> &cs) {
    const std::size_t N = cs.size();
    // Near nullspace of the velocity block: one translation per component
    // (in the equilibrated unknowns). AMGCL's default, the constant vector,
    // moves all components of an aggregate together, and the aggregates mix
    // components: 8x the iterations on 3-D LOCOS.
    std::vector<double> B(static_cast<std::size_t>(macNVel_) * D, 0.);
    for (std::size_t i = 0; i < nodes.size(); ++i)
      for (unsigned d = 0; d < D; ++d)
        for (const std::int64_t id : {macMinusFace_[d][i], macPlusFace_[d][i]})
          if (id >= 0)
            B[static_cast<std::size_t>(id) * D + d] = 1. / cs[id];
#ifdef VIENNALS_GPU_BICGSTAB
    if (gpuMode_ != GpuMode::Cpu && !macGpuOff_ && N >= macGpuMinUnknowns) {
      macGpu_ = gpu::gpuMacSetup(
          static_cast<std::int64_t>(N),
          reinterpret_cast<const std::int64_t *>(ptr.data()),
          reinterpret_cast<const std::int64_t *>(col.data()), sval.data(),
          static_cast<std::int64_t>(macNVel_), static_cast<int>(D), B.data());
      if (macGpu_)
        return true;
      const std::string message =
          std::string(
              "OxidationDeformation: GPU mode was selected, but the GPU "
              "mechanics setup failed: ") +
          gpu::gpuGetLastErrorMessage();
      if (gpuMode_ == GpuMode::Gpu) {
        VIENNACORE_LOG_ERROR(message);
        return false;
      }
      VIENNACORE_LOG_WARNING(message + " Falling back to the CPU solver.");
      macGpuOff_ = true;
    }
#endif
    try {
      typename MacAmgPrecond::params prm;
      prm.pmask.assign(N, 0);
      for (std::size_t i = static_cast<std::size_t>(macNVel_); i < N; ++i)
        prm.pmask[i] = 1;
      prm.type = 2;
      prm.approx_schur = true;
      prm.adjust_p = 2;
      auto &nullspace = prm.usolver.precond.coarsening.nullspace;
      nullspace.cols = static_cast<int>(D);
      nullspace.B = std::move(B);
      const ptrdiff_t n = static_cast<ptrdiff_t>(N);
      macAmgP_ =
          std::make_unique<MacAmgPrecond>(std::tie(n, ptr, col, sval), prm);
      typename MacAmgKrylov::params sprm; // tol: set per solve
      sprm.M = 80;
      macAmgS_ = std::make_unique<MacAmgKrylov>(N, sprm);
    } catch (const std::exception &e) {
      macAmgReset();
      VIENNACORE_LOG_WARNING(
          std::string("OxidationDeformation: AMGCL failed: ") + e.what());
      return false;
    }
    return true;
  }

  std::unique_ptr<
      Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>>
      macLU_;
  std::size_t macLUVersion_ = std::numeric_limits<std::size_t>::max();
  Eigen::Index macLUNnz_ = -1, macLUN_ = -1;
  // Equilibration of the system the LU holds (valid with macLUFactored_).
  Eigen::VectorXd macLURs_, macLUCs_;
  bool macLUFactored_ = false;

  // GMRES (right-preconditioned by the last factorisation, warm-started from
  // the last solution, <= 30 iterations, no restart) for a system whose
  // pattern the LU was factorised for. False when it does not converge, or
  // when it fails the componentwise test of lastSolveConverged() (the
  // normwise residual can be tiny while rows with small entries are off by
  // 1e-3); the fresh factorisation then solves it.
  bool macReuseSolve(const Eigen::SparseMatrix<double> &A,
                     const Eigen::VectorXd &b, T &backward) {
    const Eigen::Index N = A.rows();
    const auto precondition = [&](const Eigen::VectorXd &v) {
      return Eigen::VectorXd(
          macLUCs_.asDiagonal() *
          Eigen::VectorXd(macLU_->solve(macLURs_.asDiagonal() * v)));
    };
    Eigen::VectorXd x = Eigen::VectorXd::Zero(N);
    if (macX_.size() == static_cast<std::size_t>(N))
      for (Eigen::Index k = 0; k < N; ++k)
        x[k] = static_cast<double>(macX_[static_cast<std::size_t>(k)]);
    const double bnorm = b.norm();
    if (!(bnorm > 0.))
      return false;
    constexpr int m = 30;
    constexpr double tol = 1e-10;
    Eigen::VectorXd r = b - A * x;
    const double beta = r.norm();
    if (beta > tol * bnorm) {
      std::vector<Eigen::VectorXd> V, Z;
      Eigen::MatrixXd H = Eigen::MatrixXd::Zero(m + 1, m);
      Eigen::VectorXd g = Eigen::VectorXd::Zero(m + 1);
      std::array<double, m> c{}, s{};
      V.push_back(r / beta);
      g[0] = beta;
      int k = 0;
      bool converged = false;
      for (int j = 0; j < m; ++j) {
        Z.push_back(precondition(V[j]));
        Eigen::VectorXd w = A * Z[j];
        for (int i = 0; i <= j; ++i) {
          H(i, j) = w.dot(V[i]);
          w -= H(i, j) * V[i];
        }
        const double hNext = w.norm();
        H(j + 1, j) = hNext;
        for (int i = 0; i < j; ++i) {
          const double t = c[i] * H(i, j) + s[i] * H(i + 1, j);
          H(i + 1, j) = -s[i] * H(i, j) + c[i] * H(i + 1, j);
          H(i, j) = t;
        }
        const double den = std::hypot(H(j, j), H(j + 1, j));
        if (!(den > 0.))
          return false;
        c[j] = H(j, j) / den;
        s[j] = H(j + 1, j) / den;
        H(j, j) = den;
        H(j + 1, j) = 0.;
        g[j + 1] = -s[j] * g[j];
        g[j] = c[j] * g[j];
        k = j + 1;
        if (std::abs(g[j + 1]) <= tol * bnorm || hNext == 0.) {
          converged = true;
          break;
        }
        V.push_back(w / hNext);
      }
      if (!converged)
        return false;
      const Eigen::VectorXd y =
          H.topLeftCorner(k, k).triangularView<Eigen::Upper>().solve(g.head(k));
      for (int i = 0; i < k; ++i)
        x += y[i] * Z[i];
    }
    const Eigen::VectorXd res = b - A * x;
    if (!x.allFinite() || !(res.norm() <= 10. * tol * bnorm))
      return false;
    const Eigen::VectorXd terms = A.cwiseAbs() * x.cwiseAbs() + b.cwiseAbs();
    backward = T(0);
    for (Eigen::Index k = 0; k < N; ++k)
      if (terms[k] > 0.)
        backward =
            std::max(backward, static_cast<T>(std::abs(res[k]) / terms[k]));
    if (!(backward <= stokesTolerance))
      return false;
    macX_.assign(static_cast<std::size_t>(N), T(0));
    for (Eigen::Index k = 0; k < N; ++k)
      macX_[static_cast<std::size_t>(k)] = static_cast<T>(x[k]);
    return true;
  }

  // Sparse LU of the row/column-equilibrated system, two refinement passes.
  bool macDirectSolve(const std::vector<MacLin> &rows, T &backward) {
    const Eigen::Index N = static_cast<Eigen::Index>(rows.size());
    std::vector<Eigen::Triplet<double>> trip;
    std::size_t nnz = 0;
    for (const auto &r : rows)
      nnz += r.t.size();
    trip.reserve(nnz);
    Eigen::VectorXd rhs(N);
    for (Eigen::Index k = 0; k < N; ++k) {
      for (const auto &kv : rows[static_cast<std::size_t>(k)].t)
        trip.emplace_back(static_cast<int>(k), static_cast<int>(kv.first),
                          static_cast<double>(kv.second));
      rhs[k] = -static_cast<double>(rows[static_cast<std::size_t>(k)].c);
    }
    Eigen::SparseMatrix<double> A(N, N);
    A.setFromTriplets(trip.begin(), trip.end());
    A.makeCompressed();
    Eigen::VectorXd rs = Eigen::VectorXd::Zero(N),
                    cs = Eigen::VectorXd::Zero(N);
    for (Eigen::Index k = 0; k < A.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
        rs[it.row()] = std::max(rs[it.row()], std::abs(it.value()));
    for (Eigen::Index k = 0; k < N; ++k)
      rs[k] = rs[k] > 0. ? 1. / rs[k] : 1.;
    for (Eigen::Index k = 0; k < A.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it)
        cs[it.col()] =
            std::max(cs[it.col()], std::abs(rs[it.row()] * it.value()));
    for (Eigen::Index k = 0; k < N; ++k)
      cs[k] = cs[k] > 0. ? 1. / cs[k] : 1.;
    Eigen::SparseMatrix<double> As = rs.asDiagonal() * A * cs.asDiagonal();
    As.makeCompressed();
    // the pattern is fixed per node set: order and analyse once
    if (!macLU_ || macLUVersion_ != macNodesVersion_ ||
        macLUNnz_ != As.nonZeros() || macLUN_ != N) {
      macLU_ = std::make_unique<Eigen::SparseLU<Eigen::SparseMatrix<double>,
                                                Eigen::COLAMDOrdering<int>>>();
      macLU_->analyzePattern(As);
      macLUVersion_ = macNodesVersion_;
      macLUNnz_ = As.nonZeros();
      macLUN_ = N;
      macLUFactored_ = false;
    }
    // A coupling iteration on an unchanged node set changes the viscosity and
    // the front velocities only: GMRES on the last factorisation.
    if (couplingIterate_ && macLUFactored_ && macReuseSolve(A, rhs, backward))
      return true;
    auto &lu = *macLU_;
    lu.factorize(As);
    if (lu.info() != Eigen::Success) {
      VIENNACORE_LOG_WARNING(
          "OxidationDeformation: MAC factorisation failed: " +
          lu.lastErrorMessage());
      macLUVersion_ = std::numeric_limits<std::size_t>::max();
      macLUFactored_ = false;
      return false;
    }
    macLURs_ = rs;
    macLUCs_ = cs;
    macLUFactored_ = true;
    Eigen::VectorXd x =
        cs.asDiagonal() * Eigen::VectorXd(lu.solve(rs.asDiagonal() * rhs));
    for (int pass = 0; pass < 2; ++pass) {
      const Eigen::VectorXd r = rhs - A * x;
      x += cs.asDiagonal() * Eigen::VectorXd(lu.solve(rs.asDiagonal() * r));
    }
    const Eigen::VectorXd res = rhs - A * x;
    const Eigen::VectorXd terms = A.cwiseAbs() * x.cwiseAbs() + rhs.cwiseAbs();
    backward = T(0);
    for (Eigen::Index k = 0; k < N; ++k)
      if (terms[k] > 0.)
        backward =
            std::max(backward, static_cast<T>(std::abs(res[k]) / terms[k]));
    if (!x.allFinite())
      return false;
    macX_.assign(static_cast<std::size_t>(N), T(0));
    for (Eigen::Index k = 0; k < N; ++k)
      macX_[static_cast<std::size_t>(k)] = static_cast<T>(x[k]);
    return true;
  }

  T macEval(const MacLin &e) const {
    T v = e.c;
    for (const auto &kv : e.t)
      v += kv.second * macX_[static_cast<std::size_t>(kv.first)];
    return v;
  }

  // Node strain-rate tensor (normal rates at the node, shear rates the mean
  // of the surrounding corners) and corner shear rates.
  void macRates(std::vector<std::array<T, 9>> &nodeRate,
                std::vector<T> &cornerRate) const {
    const std::size_t n = nodes.size();
    nodeRate.assign(n, std::array<T, 9>{});
    cornerRate.assign(macCorners_.size(), T(0));
#pragma omp parallel for schedule(static)
    for (std::size_t k = 0; k < macCorners_.size(); ++k) {
      unsigned c, e;
      macPairDirs(macCorners_[k].first, c, e);
      cornerRate[k] = macEval(macShearRate(c, e, macCorners_[k].second));
    }
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < n; ++i) {
      auto &R = nodeRate[i];
      for (unsigned d = 0; d < D; ++d)
        R[tensorIndex(d, d)] = macEval(macNormalRate(d, nodes[i].index));
      for (unsigned pair = 0; pair < macPairs; ++pair) {
        unsigned c, e;
        macPairDirs(pair, c, e);
        T sum = T(0);
        int cnt = 0;
        for (int dc = -1; dc <= 0; ++dc)
          for (int de = -1; de <= 0; ++de) {
            IndexType lo = nodes[i].index;
            lo[c] += dc;
            lo[e] += de;
            T sign;
            const std::int64_t k = macCornerReflected(pair, lo, sign);
            if (k < 0)
              continue;
            sum += sign * cornerRate[static_cast<std::size_t>(k)];
            ++cnt;
          }
        const T g = cnt ? T(0.5) * sum / static_cast<T>(cnt) : T(0);
        R[tensorIndex(c, e)] = g;
        R[tensorIndex(e, c)] = g;
      }
    }
  }

  // Anderson acceleration (history 5) of the viscosity fixed point u = G(u),
  // u = log eta at the oxide nodes and corners, G the consistent viscosity
  // for the strain rate of the solution with u.  Restarted when the
  // residual more than doubles; the first step is the relaxed mean.  On the
  // stiff-nitride LOCOS stack it needs 2.7x fewer iterations than the
  // relaxed (1/2) fixed point, with the same result.
  std::vector<std::vector<T>> macAAu_, macAAg_;
  void macAndersonStep(std::vector<T> &u, const std::vector<T> &g) {
    constexpr int m = 5;
    const std::size_t L = u.size();
    macAAu_.push_back(u);
    macAAg_.push_back(g);
    if (static_cast<int>(macAAu_.size()) > m + 1) {
      macAAu_.erase(macAAu_.begin());
      macAAg_.erase(macAAg_.begin());
    }
    const std::size_t h = macAAu_.size();
    // residuals f_j = g_j - u_j; restart if the newest grew much
    auto fnorm = [&](std::size_t j) {
      T sq = T(0);
      for (std::size_t k = 0; k < L; ++k) {
        const T f = macAAg_[j][k] - macAAu_[j][k];
        sq += f * f;
      }
      return std::sqrt(sq);
    };
    if (h >= 2 && fnorm(h - 1) > T(2) * fnorm(h - 2)) {
      macAAu_.erase(macAAu_.begin(), macAAu_.end() - 1);
      macAAg_.erase(macAAg_.begin(), macAAg_.end() - 1);
    }
    const std::size_t hh = macAAu_.size();
    if (hh < 2) { // plain relaxed step
      for (std::size_t k = 0; k < L; ++k)
        u[k] = T(0.5) * (u[k] + g[k]);
      return;
    }
    // least squares: min || f_last - sum_j gamma_j (f_last - f_j) ||
    const std::size_t q = hh - 1;
    std::vector<std::vector<T>> dF(q, std::vector<T>(L));
    std::vector<T> fl(L);
    for (std::size_t k = 0; k < L; ++k)
      fl[k] = macAAg_[q][k] - macAAu_[q][k];
    for (std::size_t j = 0; j < q; ++j)
      for (std::size_t k = 0; k < L; ++k)
        dF[j][k] = fl[k] - (macAAg_[j][k] - macAAu_[j][k]);
    std::vector<T> M(q * q, T(0)), rhs(q, T(0)), gamma(q, T(0));
    for (std::size_t a = 0; a < q; ++a) {
      for (std::size_t b = 0; b < q; ++b) {
        T sum = T(0);
        for (std::size_t k = 0; k < L; ++k)
          sum += dF[a][k] * dF[b][k];
        M[a * q + b] = sum;
      }
      T sum = T(0);
      for (std::size_t k = 0; k < L; ++k)
        sum += dF[a][k] * fl[k];
      rhs[a] = sum;
      M[a * q + a] *= T(1) + T(1e-10); // regularise
    }
    // Gaussian elimination (q <= m small)
    for (std::size_t c = 0; c < q; ++c) {
      std::size_t piv = c;
      for (std::size_t r = c + 1; r < q; ++r)
        if (std::abs(M[r * q + c]) > std::abs(M[piv * q + c]))
          piv = r;
      if (std::abs(M[piv * q + c]) < T(1e-300))
        continue;
      if (piv != c) {
        for (std::size_t k = 0; k < q; ++k)
          std::swap(M[c * q + k], M[piv * q + k]);
        std::swap(rhs[c], rhs[piv]);
      }
      for (std::size_t r = c + 1; r < q; ++r) {
        const T f = M[r * q + c] / M[c * q + c];
        for (std::size_t k = c; k < q; ++k)
          M[r * q + k] -= f * M[c * q + k];
        rhs[r] -= f * rhs[c];
      }
    }
    for (std::size_t c = q; c-- > 0;) {
      T sum = rhs[c];
      for (std::size_t k = c + 1; k < q; ++k)
        sum -= M[c * q + k] * gamma[k];
      gamma[c] = std::abs(M[c * q + c]) > T(1e-300) ? sum / M[c * q + c] : T(0);
    }
    // u_new = g_last - sum_j gamma_j (g_last - g_j)
    for (std::size_t k = 0; k < L; ++k) {
      T val = macAAg_[q][k];
      for (std::size_t j = 0; j < q; ++j)
        val -= gamma[j] * (macAAg_[q][k] - macAAg_[j][k]);
      u[k] = val;
    }
  }

  void macUpdateViscosity(const std::vector<std::array<T, 9>> &nodeRate,
                          const std::vector<T> &cornerRate) {
    const bool oxide = stressDependentViscosity();
    const bool nitride = stressDependentNitride();
    if (!oxide && !nitride)
      return;
    const std::size_t n = nodes.size();
    std::vector<T> locNode(n, T(0)), locCorner(macCorners_.size(), T(0)),
        locCornerN(macCorners_.size(), T(0));
    // dynamic: the bisection is costly only where the stress is high
#pragma omp parallel for schedule(dynamic, 64)
    for (std::size_t i = 0; i < n; ++i) {
      if (nodes[i].nitride ? !nitride : !oxide)
        continue;
      std::array<T, 9> start =
          nodes[i].nitride ? deviatorOf(nitrideStartOf(i)) : startStress_[i];
      for (unsigned pair = 0; pair < macPairs; ++pair) { // corner mean
        unsigned c, e;
        macPairDirs(pair, c, e);
        T sum = T(0);
        int cnt = 0;
        for (int dc = -1; dc <= 0; ++dc)
          for (int de = -1; de <= 0; ++de) {
            IndexType lo = nodes[i].index;
            lo[c] += dc;
            lo[e] += de;
            T sign;
            const std::int64_t k = macCornerReflected(pair, lo, sign);
            if (k < 0)
              continue;
            sum += sign * macCornerStartV_[static_cast<std::size_t>(k)];
            ++cnt;
          }
        start[tensorIndex(c, e)] = start[tensorIndex(e, c)] =
            cnt ? sum / static_cast<T>(cnt) : T(0);
      }
      const auto &R = nodeRate[i];
      T tr = T(0);
      for (unsigned d = 0; d < D; ++d)
        tr += R[tensorIndex(d, d)];
      std::array<T, 9> dev = R;
      for (unsigned d = 0; d < 3; ++d)
        dev[tensorIndex(d, d)] -= tr / T(3);
      locNode[i] = nodes[i].nitride ? localNitrideViscosity(start, dev)
                                    : localViscosity(start, dev);
    }
    // corner viscosity of one material from its points' normal start stress
    // and rate and the corner's shear
    auto cornerViscosity = [&](std::size_t k, bool ofNitride) {
      unsigned c, e;
      macPairDirs(macCorners_[k].first, c, e);
      const IndexType &lo = macCorners_[k].second;
      std::array<T, 9> start{}, rate{};
      int cnt = 0;
      for (int bc = 0; bc <= 1; ++bc)
        for (int be = 0; be <= 1; ++be) {
          IndexType q = lo;
          q[c] += bc;
          q[e] += be;
          const std::size_t j = macNode(q);
          if (j == noNode || nodes[j].nitride != ofNitride)
            continue;
          const auto s =
              ofNitride ? deviatorOf(nitrideStartOf(j)) : startStress_[j];
          for (unsigned d = 0; d < 3; ++d) {
            start[tensorIndex(d, d)] += s[tensorIndex(d, d)];
            rate[tensorIndex(d, d)] += nodeRate[j][tensorIndex(d, d)];
          }
          ++cnt;
        }
      if (cnt == 0)
        return T(0);
      for (unsigned d = 0; d < 3; ++d) {
        start[tensorIndex(d, d)] /= static_cast<T>(cnt);
        rate[tensorIndex(d, d)] /= static_cast<T>(cnt);
      }
      start[tensorIndex(c, e)] = start[tensorIndex(e, c)] = macCornerStartV_[k];
      rate[tensorIndex(c, e)] = rate[tensorIndex(e, c)] =
          T(0.5) * cornerRate[k];
      T tr = T(0);
      for (unsigned d = 0; d < D; ++d)
        tr += rate[tensorIndex(d, d)];
      for (unsigned d = 0; d < 3; ++d)
        rate[tensorIndex(d, d)] -= tr / T(3);
      return ofNitride ? localNitrideViscosity(start, rate)
                       : localViscosity(start, rate);
    };
#pragma omp parallel for schedule(dynamic, 64)
    for (std::size_t k = 0; k < macCorners_.size(); ++k) {
      if (macCornerGas_[k])
        continue;
      if (oxide && macCornerOxide_[k])
        locCorner[k] = cornerViscosity(k, false);
      if (nitride)
        locCornerN[k] = cornerViscosity(k, true);
    }
    // Anderson-accelerated fixed-point update in log eta
    std::vector<T> u, g;
    for (std::size_t i = 0; i < n; ++i)
      if (locNode[i] > T(0)) {
        u.push_back(std::log(nodes[i].nitride ? nitrideEta_[i] : oxideEta_[i]));
        g.push_back(std::log(locNode[i]));
      }
    for (std::size_t k = 0; k < macCorners_.size(); ++k) {
      if (locCorner[k] > T(0)) {
        u.push_back(std::log(macCornerEta_[k]));
        g.push_back(std::log(locCorner[k]));
      }
      if (locCornerN[k] > T(0)) {
        u.push_back(std::log(macCornerEtaN_[k]));
        g.push_back(std::log(locCornerN[k]));
      }
    }
    macAndersonStep(u, g);
    const T hi = std::log(T(deformationParameters.viscosity));
    const T lo = hi + std::log(T(1e-6));
    const T hiN = std::log(std::max(T(deformationParameters.nitrideViscosity),
                                    std::numeric_limits<T>::min()));
    const T loN = hiN + std::log(T(1e-6));
    std::size_t q = 0;
    for (std::size_t i = 0; i < n; ++i)
      if (locNode[i] > T(0)) {
        if (nodes[i].nitride)
          nitrideEta_[i] = std::exp(std::clamp(u[q++], loN, hiN));
        else
          oxideEta_[i] = std::exp(std::clamp(u[q++], lo, hi));
      }
    for (std::size_t k = 0; k < macCorners_.size(); ++k) {
      if (locCorner[k] > T(0))
        macCornerEta_[k] = std::exp(std::clamp(u[q++], lo, hi));
      if (locCornerN[k] > T(0))
        macCornerEtaN_[k] = std::exp(std::clamp(u[q++], loN, hiN));
    }
  }

  void macInitViscosity() {
    if (stressDependentViscosity())
      initOxideViscosity(); // nodes, from the start stress
    if (stressDependentNitride())
      initNitrideViscosity();
    const T eta0 = T(deformationParameters.viscosity);
    const T etaN0 = T(deformationParameters.nitrideViscosity);
    const T oxideVolume = T(deformationParameters.shearActivationVolume);
    const T nitrideVolume =
        T(deformationParameters.nitrideShearActivationVolume);
    for (std::size_t k = 0; k < macCorners_.size(); ++k) {
      unsigned c, e;
      macPairDirs(macCorners_[k].first, c, e);
      std::array<T, 9> s{};
      s[tensorIndex(c, e)] = s[tensorIndex(e, c)] = macCornerStartV_[k];
      macCornerEta_[k] = stressDependentViscosity() && macCornerOxide_[k]
                             ? eta0 * eyringFactor(s, oxideVolume)
                             : eta0;
      macCornerEtaN_[k] = stressDependentNitride()
                              ? etaN0 * eyringFactor(s, nitrideVolume)
                              : etaN0;
    }
    macAAu_.clear();
    macAAg_.clear();
  }

  // Write the converged MAC solution into the node arrays the rest of the
  // class (output, stress history, level-set velocity) uses.
  void macWriteBack(const std::vector<std::array<T, 9>> &nodeRate,
                    const std::vector<T> &cornerRate, T etaRef) {
    const std::size_t n = nodes.size();
    if (oxideStressEnd_.size() != n)
      oxideStressEnd_.assign(n, std::array<T, 9>{});
    oxideStressAdvectedValid_ = false;
    // corner end stresses (the stress the momentum used); dynamic: the
    // costly surface corners are clustered
#pragma omp parallel for schedule(dynamic, 64)
    for (std::size_t k = 0; k < macCorners_.size(); ++k) {
      unsigned c, e;
      macPairDirs(macCorners_[k].first, c, e);
      macCornerEndV_[k] =
          macCornerGas_[k]
              ? macEval(macShearStress(c, e, macCorners_[k].second, etaRef)) *
                    etaRef
              : (macCornerMod_[k] * cornerRate[k] + macCornerMem_[k]) * etaRef;
    }
#pragma omp parallel for schedule(dynamic, 64)
    for (std::size_t i = 0; i < n; ++i) {
      auto &node = nodes[i];
      for (unsigned d = 0; d < D; ++d) {
        IndexType lo = node.index;
        lo[d] -= 1;
        const IndexType hi = node.index;
        node.velocity[d] = T(0.5) * (macEval(macFace(d, lo, &hi)) +
                                     macEval(macFace(d, hi, &lo)));
      }
      node.pressure = macPid_[i] >= 0
                          ? macX_[static_cast<std::size_t>(macPid_[i])] * etaRef
                          : T(0);
      node.strainRateTensor = nodeRate[i];
      T tr = T(0);
      for (unsigned d = 0; d < D; ++d)
        tr += nodeRate[i][tensorIndex(d, d)];
      node.strainTrace = tr;
      // node shear stress: mean of the surrounding corners' end stress
      std::array<T, 9> shear{};
      for (unsigned pair = 0; pair < macPairs; ++pair) {
        unsigned c, e;
        macPairDirs(pair, c, e);
        T sum = T(0);
        int cnt = 0;
        for (int dc = -1; dc <= 0; ++dc)
          for (int de = -1; de <= 0; ++de) {
            IndexType lo = node.index;
            lo[c] += dc;
            lo[e] += de;
            T sign;
            const std::int64_t k = macCornerReflected(pair, lo, sign);
            if (k < 0)
              continue;
            sum += sign * macCornerEndV_[static_cast<std::size_t>(k)];
            ++cnt;
          }
        shear[tensorIndex(c, e)] = shear[tensorIndex(e, c)] =
            cnt ? sum / static_cast<T>(cnt) : T(0);
      }
      if (node.nitride) {
        computeNitrideStress(i); // normal end stress from the strain rate
        auto &end = nitrideStressEnd_[i];
        for (unsigned c = 0; c < 3; ++c)
          for (unsigned e = 0; e < 3; ++e)
            if (c != e)
              end[tensorIndex(c, e)] = shear[tensorIndex(c, e)];
        node.stressTensor = end;
        continue;
      }
      const T decay = nodeDecay(i);
      const T eta = nodeViscosity(i);
      const auto &start = startStress_[i];
      std::array<T, 9> end = shear;
      for (unsigned d = 0; d < 3; ++d)
        end[tensorIndex(d, d)] =
            decay * start[tensorIndex(d, d)] +
            (T(1) - decay) * T(2) * eta *
                (nodeRate[i][tensorIndex(d, d)] - tr / T(3));
      oxideStressEnd_[i] = end;
      node.stressTensor = end;
      for (unsigned d = 0; d < 3; ++d)
        node.stressTensor[tensorIndex(d, d)] -= node.pressure;
      std::array<T, 9> dev = end;
      const T mean = (end[0] + end[4] + end[8]) / T(3);
      for (unsigned d = 0; d < 3; ++d)
        dev[tensorIndex(d, d)] -= mean;
      node.vonMisesStress = vonMisesFromDeviatoric(dev);
    }
    // carried shear history for the next step (advected on write)
    macCornerEndMap_.clear();
    for (std::size_t k = 0; k < macCorners_.size(); ++k) {
      auto &slot = macCornerEndMap_[macCorners_[k].second];
      slot[macCorners_[k].first] = macCornerEndV_[k];
    }
    macCornerAdvectedValid_ = false;
  }

  // The mechanics solve: the viscosity iteration of the stress-dependent
  // oxide, one exact velocity-pressure solve per iteration.
  void solveMechanics() {
    // node set version: rebuild the numbering when the nodes changed
    if (macNodesVersion_ != nodesVersion_ || macPid_.size() != nodes.size())
      macBuild();
    macGhostReset();
    struct GhostOff { // the kept front velocities end with this solve
      bool &on;
      ~GhostOff() { on = false; }
    } ghostOff{macGhostOn_};
    nitrideStartNode_.clear();
    if (monolithicNitride()) {
      nitrideStartNode_.assign(nodes.size(), std::array<T, 9>{});
      for (std::size_t i = 0; i < nodes.size(); ++i)
        if (nodes[i].nitride)
          nitrideStartNode_[i] = nitrideStressStartAt(nodes[i].index);
    }
    if (nitrideStressEnd_.size() != nodes.size())
      nitrideStressEnd_.assign(nodes.size(), std::array<T, 9>{});
    macLoadCornerStart();
    // A re-solve on the same node set (the next coupling iteration, a CFL
    // re-solve) starts from the last converged viscosities: the forcing
    // changes little, so the viscosity iteration needs a few steps instead
    // of a restart from the step-start stress.
    const bool nonlinear =
        stressDependentViscosity() || stressDependentNitride();
    if (nonlinear && macWarmVersion_ == macNodesVersion_ &&
        oxideEta_.size() == nodes.size() &&
        nitrideEta_.size() == nodes.size() &&
        macCornerEta_.size() == macCorners_.size()) {
      macAAu_.clear();
      macAAg_.clear();
    } else {
      macInitViscosity();
    }
    if (oxideEta_.size() != nodes.size())
      oxideEta_.assign(nodes.size(), T(deformationParameters.viscosity));
    if (nitrideEta_.size() != nodes.size())
      nitrideEta_.assign(nodes.size(),
                         T(deformationParameters.nitrideViscosity));
    const T etaRef = implicitViscosity();
    std::vector<std::array<T, 9>> nodeRate;
    std::vector<T> cornerRate;
    std::vector<T> previous;
    T change = std::numeric_limits<T>::infinity();
    T backward = T(0);
    unsigned it = 0;
    const unsigned maxIt =
        std::max(1u, deformationParameters.mechanicsIterations);
    const bool onePass = couplingIterate_ && nonlinear;
    bool ok = true;
    for (it = 1; it <= maxIt; ++it) {
      if (it > 1)
        macUpdateViscosity(nodeRate, cornerRate);
      previous = macX_;
      if (!macSolveLinear(etaRef, backward)) {
        ok = false;
        break;
      }
      macRates(nodeRate, cornerRate);
      T dmax = T(0), vmax = T(0);
      for (std::int64_t k = 0; k < macNVel_; ++k) {
        const std::size_t kk = static_cast<std::size_t>(k);
        vmax = std::max(vmax, std::abs(macX_[kk]));
        if (kk < previous.size())
          dmax = std::max(dmax, std::abs(macX_[kk] - previous[kk]));
      }
      change = vmax > T(0) ? dmax / vmax : dmax;
      if (onePass)
        break;
      if (!nonlinear ||
          (it > 1 && change < T(deformationParameters.mechanicsTolerance)))
        break;
    }
    iterations = it;
    if (!ok) {
      residual = std::numeric_limits<T>::infinity();
      lastStokesResidual_ = std::numeric_limits<T>::infinity();
      macWarmVersion_ = std::numeric_limits<std::size_t>::max();
      return;
    }
    macWriteBack(nodeRate, cornerRate, etaRef);
    residual = nonlinear ? change : T(0);
    if (onePass) { // the viscosity for the next coupling iteration
      macUpdateViscosity(nodeRate, cornerRate);
      residual = T(0);
    }
    lastStokesResidual_ = backward;
    { // volume the MAC flow fails to conserve, relative to the inflow
      T sumDiv = T(0), cellVolume = T(1);
      for (unsigned d = 0; d < D; ++d)
        cellVolume *= gridDelta;
      for (std::size_t i = 0; i < nodes.size(); ++i)
        if (!nodes[i].nitride)
          sumDiv += std::abs(nodes[i].strainTrace);
      const T inflow = reactionInflow();
      lastContinuityResidual_ = inflow > std::numeric_limits<T>::epsilon()
                                    ? sumDiv * cellVolume / inflow
                                    : sumDiv * cellVolume;
    }
    logContinuityDiagnostic(it);
    macWarmVersion_ = residual <= deformationParameters.mechanicsTolerance
                          ? macNodesVersion_
                          : std::numeric_limits<std::size_t>::max();
    if (residual > deformationParameters.mechanicsTolerance)
      VIENNACORE_LOG_WARNING(
          "solveMechanics: viscosity iteration did not converge after " +
          std::to_string(it - 1) +
          " iterations (change=" + std::to_string(residual) + ", tolerance=" +
          std::to_string(deformationParameters.mechanicsTolerance) + ")");
  }

  // Carry the corner shear history with the flow (semi-Lagrangian over the
  // corner lattice of the same pair), once per accepted step.
  void macAdvectCorners() {
    if (macCornerAdvectedValid_)
      return;
    macCornerAdvectedValid_ = true;
    if (macCornerEndMap_.empty())
      return;
    const T dt = stepTime();
    std::unordered_map<IndexType, std::array<T, 3>, detail::IndexTypeHasher<D>>
        out;
    for (std::size_t k = 0; k < macCorners_.size(); ++k) {
      const unsigned pair = macCorners_[k].first;
      unsigned c, e;
      macPairDirs(pair, c, e);
      const IndexType &lo = macCorners_[k].second;
      // corner velocity: mean of the nodes around it
      Vec3D<T> v{T(0), T(0), T(0)};
      int cnt = 0;
      for (int bc = 0; bc <= 1; ++bc)
        for (int be = 0; be <= 1; ++be) {
          IndexType q = lo;
          q[c] += bc;
          q[e] += be;
          const std::size_t j = macNode(q);
          if (j == noNode)
            continue;
          for (unsigned d = 0; d < D; ++d)
            v[d] += nodes[j].velocity[d];
          ++cnt;
        }
      T value = macCornerEndV_[k];
      if (cnt > 0) {
        std::array<T, D> base{}, frac{}, departure{};
        for (unsigned d = 0; d < D; ++d) {
          const T x = static_cast<T>(lo[d]) - v[d] / cnt * dt / gridDelta;
          base[d] = std::floor(x);
          frac[d] = x - base[d];
          // the corner sits half a cell along both directions of its pair
          departure[d] = x + (d == c || d == e ? T(0.5) : T(0));
        }
        // upstream of the Si front: oxide born this step, no history
        if (macCornerOxide_[k] && macPhiRAt(departure) < T(0)) {
          auto &slot = out.try_emplace(lo, macNoHistory()).first->second;
          slot[pair] = T(0);
          continue;
        }
        T sum = T(0), wsum = T(0);
        for (unsigned corner = 0; corner < (1u << D); ++corner) {
          IndexType idx;
          T w = T(1);
          for (unsigned d = 0; d < D; ++d) {
            const bool up = (corner >> d) & 1u;
            idx[d] = static_cast<viennahrle::IndexType>(base[d]) + (up ? 1 : 0);
            w *= up ? frac[d] : (T(1) - frac[d]);
          }
          if (w <= T(0))
            continue;
          T sign;
          const std::int64_t kk = macCornerReflected(pair, idx, sign);
          if (kk < 0 || macCornerGas_[static_cast<std::size_t>(kk)])
            continue;
          sum += w * sign * macCornerEndV_[static_cast<std::size_t>(kk)];
          wsum += w;
        }
        if (wsum > T(1e-6))
          value = sum / wsum;
      }
      auto &slot = out.try_emplace(lo, macNoHistory()).first->second;
      slot[pair] = value;
    }
    macCornerStart_ = std::move(out);
  }

  // "No history" marker for a corner pair: far outside any stress value and
  // finite in single precision -- the arrays go to the output files as
  // Float32, where NaN or an overflow to infinity makes the norm range VTK
  // stores with them unreadable.
  static constexpr T macNone = T(1e30);
  static std::array<T, 3> macNoHistory() { return {macNone, macNone, macNone}; }
  static bool macHasHistory(T value) { return std::abs(value) < T(1e29); }

  // The corner shear history (advected to the end of the step) is stored on
  // the level sets so that the next solver -- a new object every step, and in
  // a later process call -- starts from it.  Every defined point keeps the
  // corners it is a vertex of, one array per vertex position k (offset
  // (k & 1, k >> 1) from the corner's lower point in the pair's plane), one
  // component per pair: a corner survives if any of its points is defined,
  // also next to Si where neither level set is.  On the ambient level set
  // ("OxCornerShear0..3") and, with a nitride, on the mask level set
  // ("NitCornerShear0..3").
  IndexType macCornerLower(IndexType point, unsigned pair, unsigned k) const {
    unsigned c, e;
    macPairDirs(pair, c, e);
    point[c] -= static_cast<int>(k & 1u);
    point[e] -= static_cast<int>(k >> 1);
    return point;
  }

  void macWriteCornerHistory() {
    macAdvectCorners();
    auto store = [&](SmartPointer<Domain<T, D>> levelSet,
                     const std::string &name) {
      using VD = typename PointData<T>::VectorDataType;
      std::array<VD, 4> data;
      ConstSparseIterator it(levelSet->getDomain());
      for (; !it.isFinished(); ++it) {
        if (!it.isDefined())
          continue;
        const IndexType point = it.getStartIndices();
        for (unsigned k = 0; k < 4; ++k) {
          Vec3D<T> value{macNone, macNone, macNone};
          for (unsigned pair = 0; pair < macPairs; ++pair) {
            const auto found =
                macCornerStart_.find(macCornerLower(point, pair, k));
            if (found != macCornerStart_.end())
              value[pair] = found->second[pair];
          }
          data[k].push_back(value);
        }
      }
      for (unsigned k = 0; k < 4; ++k)
        levelSet->getPointData().insertReplaceVectorData(
            std::move(data[k]), name + std::to_string(k));
    };
    store(ambientInterface, "OxCornerShear");
    if (monolithicNitride())
      store(maskInterface, "NitCornerShear");
  }

  void macReadCornerHistory() {
    macCornerStart_.clear();
    auto load = [&](SmartPointer<Domain<T, D>> levelSet,
                    const std::string &name) {
      auto &pd = levelSet->getPointData();
      for (unsigned k = 0; k < 4; ++k) {
        const int index = pd.getVectorDataIndex(name + std::to_string(k));
        if (index == -1)
          return;
        const auto *data = pd.getVectorData(index);
        ConstSparseIterator it(levelSet->getDomain());
        for (; !it.isFinished(); ++it) {
          if (!it.isDefined())
            continue;
          const auto id = it.getPointId();
          if (id >= static_cast<decltype(id)>(data->size()))
            continue;
          const auto &v = (*data)[id];
          for (unsigned pair = 0; pair < macPairs; ++pair) {
            if (!macHasHistory(v[pair]))
              continue;
            auto &slot =
                macCornerStart_
                    .try_emplace(macCornerLower(it.getStartIndices(), pair, k),
                                 macNoHistory())
                    .first->second;
            if (!macHasHistory(slot[pair]))
              slot[pair] = v[pair];
          }
        }
      }
    };
    load(ambientInterface, "OxCornerShear");
    if (maskInterface != nullptr)
      load(maskInterface, "NitCornerShear");
  }

  // Oxide volume created per unit time at the Si/SiO2 interface (per unit
  // length in 2D): sum over REACTION faces of |v_exp . e_dir| * h^(D-1).
  T reactionInflow() const {
    const std::size_t n = nodes.size();
    T faceArea = T(1);
    for (unsigned d = 0; d + 1 < D; ++d)
      faceArea *= gridDelta;
    T inflow = T(0);
#pragma omp parallel for reduction(+ : inflow) schedule(static)
    for (std::size_t i = 0; i < n; ++i)
      for (unsigned fi = 0; fi < 2 * D; ++fi)
        if (faceBCTypes_[fi * n + i] == Boundary::REACTION) {
          const auto v = reactionVelocityAt(i);
          inflow += std::abs(v[fi / 2]) * faceArea;
        }
    return inflow;
  }

  // One INFO line per mechanics solve: how well the final velocity satisfies
  // continuity.  netVolumeRate = Σ div(v)·gridDelta^D over the oxide nodes
  // is the volume the discrete flow creates (+) or destroys (−) inside the
  // oxide per hour, in µm^D/hr.  For an incompressible oxide it must vanish.
  void logContinuityDiagnostic(unsigned iterationsUsed) const {
    const std::size_t n = nodes.size();
    if (n == 0)
      return;
    T maxDiv = T(0), sumSq = T(0), net = T(0);
    T pMin = std::numeric_limits<T>::max();
    T pMax = std::numeric_limits<T>::lowest();
    std::size_t count = 0;
    for (std::size_t i = 0; i < n; ++i) {
      if (nodes[i].nitride)
        continue; // compressible elastic body without a pressure
      const T p = nodes[i].pressure;
      if (std::isfinite(p)) {
        pMin = std::min(pMin, p);
        pMax = std::max(pMax, p);
      }
      const T div = nodes[i].strainTrace;
      if (!std::isfinite(div))
        continue;
      maxDiv = std::max(maxDiv, std::abs(div));
      sumSq += div * div;
      net += div;
      ++count;
    }
    T cellVolume = T(1);
    for (unsigned d = 0; d < D; ++d)
      cellVolume *= gridDelta;
    if (pMin > pMax)
      pMin = pMax = T(0);
    char buffer[384];
    std::snprintf(
        buffer, sizeof(buffer),
        "OxidationDeformation: mechanics iterations=%u "
        "residual=%.3e continuity=%.3e |div v|max=%.3e rms=%.3e /hr "
        "netVolumeRate=%.4e "
        "um^%d/hr p=[%.3f, %.3f] GPa",
        iterationsUsed, static_cast<double>(residual),
        static_cast<double>(lastContinuityResidual_),
        static_cast<double>(maxDiv),
        static_cast<double>(count ? std::sqrt(sumSq / T(count)) : T(0)),
        static_cast<double>(net * cellVolume), D,
        static_cast<double>(pMin * T(1e-9)),
        static_cast<double>(pMax * T(1e-9)));
    VIENNACORE_LOG_INFO(std::string(buffer));
  }

  void computeAvgExpansionSpeed() {
    avgExpansionSpeed_ = 0.;
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

          if (lookupNode(neighbor) != noNode)
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
        avgExpansionSpeed_ +=
            (oxidationParameters.expansionCoefficient - T(1)) *
            std::abs(diffusionField->getScalarVelocity(coordinate, 0,
                                                       {0., 0., 0.}, 0));
        ++count;
      }
    }

    if (count > 0)
      avgExpansionSpeed_ /= static_cast<T>(count);
  }

  void cacheBoundaryValues() {
    boundaryCacheValid_ = false;
    const std::size_t n = nodes.size();
    reactionVelocity_.assign(n, Vec3D<T>{T(0), T(0), T(0)});
    startStress_.assign(n, std::array<T, 9>{});
    if (diffusionField != nullptr) // a pending diffusion solve runs here,
      diffusionField->getScalarVelocity({T(0), T(0), T(0)}, 0,  // not in the
                                        {T(0), T(0), T(0)}, 0); // loop below
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < n; ++i) {
      bool reaction = false;
      for (unsigned fi = 0; fi < 2 * D; ++fi)
        reaction = reaction || faceBCTypes_[fi * n + i] == Boundary::REACTION;
      if (reaction)
        reactionVelocity_[i] = reactionBoundaryVelocity(nodes[i].index);
      startStress_[i] = previousDeviatoricStress(nodes[i].index);
    }
    boundaryCacheValid_ = true;
  }

  Vec3D<T> reactionVelocityAt(std::size_t i) const {
    return boundaryCacheValid_ ? reactionVelocity_[i]
                               : reactionBoundaryVelocity(nodes[i].index);
  }

  std::array<T, 9> startStressAt(const IndexType &index) const {
    if (boundaryCacheValid_) {
      const std::size_t i = lookupNode(index);
      if (i != noNode)
        return startStress_[i];
    }
    return previousDeviatoricStress(index);
  }

  Vec3D<T> reactionBoundaryVelocity(const IndexType &index) const {
    Vec3D<T> coordinate{0., 0., 0.};
    for (unsigned i = 0; i < D; ++i)
      coordinate[i] = index[i] * gridDelta;
    const T expansionVelocity = localExpansionSpeed(coordinate);
    return detail::vecScaled(reactionNormal(index),
                             reactionSign * expansionVelocity);
  }

  Vec3D<T> unresolvedAmbientVelocity(const Vec3D<T> &coordinate) const {
    if (diffusionField == nullptr || ambientInterface == nullptr)
      return {0., 0., 0.};

    IndexType index;
    for (unsigned i = 0; i < D; ++i)
      index[i] = std::llround(coordinate[i] / gridDelta);

    ConstSparseIterator ambientIt(ambientInterface->getDomain());
    const auto normal = levelSetNormal(ambientIt, index);
    return detail::vecScaled(normal, localExpansionSpeed(coordinate));
  }

  Vec3D<T> estimateMaxUnresolvedAmbientVelocity() const {
    Vec3D<T> maxVelocity{0., 0., 0.};
    if (ambientInterface == nullptr || diffusionField == nullptr)
      return maxVelocity;

    // Every defined point of the ambient interface, in parallel as in
    // cacheBoundaryValues (each query uses its own iterators; the max does
    // not depend on the order).  Far from the Si front a query searches
    // ~10^3 nodes: on one thread this walk took 1/4 of a 3-D LOCOS run.
    std::vector<IndexType> points;
    for (ConstSparseIterator ambientIt(ambientInterface->getDomain());
         !ambientIt.isFinished(); ++ambientIt)
      if (ambientIt.isDefined())
        points.push_back(ambientIt.getStartIndices());
    if (points.empty())
      return maxVelocity;
    diffusionField->getScalarVelocity({T(0), T(0), T(0)}, 0,  // a pending
                                      {T(0), T(0), T(0)}, 0); // solve runs here
#pragma omp parallel
    {
      Vec3D<T> local{0., 0., 0.};
#pragma omp for schedule(dynamic, 64) nowait
      for (std::size_t k = 0; k < points.size(); ++k) {
        Vec3D<T> coordinate{0., 0., 0.};
        for (unsigned d = 0; d < D; ++d)
          coordinate[d] = points[k][d] * gridDelta;

        const auto velocity = unresolvedAmbientVelocity(coordinate);
        for (unsigned d = 0; d < D; ++d)
          local[d] = std::max(local[d], std::abs(velocity[d]));
      }
#pragma omp critical
      for (unsigned d = 0; d < D; ++d)
        maxVelocity[d] = std::max(maxVelocity[d], local[d]);
    }
    return maxVelocity;
  }

  T localExpansionSpeed(const Vec3D<T> &coordinate) const {
    return (oxidationParameters.expansionCoefficient - T(1)) *
           std::abs(diffusionField->getScalarVelocity(coordinate, 0,
                                                      {0., 0., 0.}, 0));
  }

  Vec3D<T> reactionNormal(const IndexType &index) const {
    ConstSparseIterator reactionIt(reactionInterface->getDomain());
    return levelSetNormal(reactionIt, index);
  }

  Vec3D<T> levelSetNormal(ConstSparseIterator &levelSetIt,
                          const IndexType &index) const {
    Vec3D<T> normal = this->definedLevelSetGradient(
        levelSetIt, index, [&](const IndexType &q, bool &ok) {
          ok = inBounds(q);
          return q;
        });
    T norm = 0.;
    for (unsigned i = 0; i < D; ++i)
      norm += normal[i] * normal[i];

    if (norm <= std::numeric_limits<T>::epsilon()) {
      normal = Vec3D<T>{0., 0., 0.};
      normal[D - 1] = 1.;
      return normal;
    }

    norm = std::sqrt(norm);
    for (unsigned i = 0; i < D; ++i)
      normal[i] /= norm;
    return normal;
  }

  std::array<T, 9> previousDeviatoricStress(const IndexType &index) const {
    const auto found = deviatoricStressHistory.find(index);
    if (found == deviatoricStressHistory.end())
      return {};
    return found->second;
  }

  // exp(-dt/tau) of the Maxwell oxide; 0 for a purely viscous oxide.
  T maxwellDecay() const {
    const T tau = effectiveStressRelaxationTime();
    if (tau <= std::numeric_limits<T>::epsilon())
      return T(0);
    return std::exp(-T(deformationParameters.stressTimeStep) / tau);
  }

  // Instantaneous (implicit) modulus of the oxide over one substep:
  // (1 - decay) eta, about G dt for dt << tau, eta when purely viscous.
  T implicitViscosity() const {
    return std::max(T(deformationParameters.viscosity) *
                        (T(1) - maxwellDecay()),
                    std::numeric_limits<T>::min() * T(1e10));
  }

  bool stressDependentViscosity() const {
    return deformationParameters.shearActivationVolume > 0.;
  }

  // Oxide viscosity, Maxwell decay and implicit modulus at node i; the
  // constant-viscosity values unless the viscosity is stress-dependent.
  T nodeViscosity(std::size_t i) const {
    return (stressDependentViscosity() && i < oxideEta_.size())
               ? oxideEta_[i]
               : T(deformationParameters.viscosity);
  }

  T nodeDecay(std::size_t i) const {
    if (!stressDependentViscosity())
      return maxwellDecay();
    const T tau0 = effectiveStressRelaxationTime();
    if (tau0 <= std::numeric_limits<T>::epsilon())
      return T(0);
    const T tau = tau0 * nodeViscosity(i) / T(deformationParameters.viscosity);
    return std::exp(-T(deformationParameters.stressTimeStep) / tau);
  }

  T nodeImplicitViscosity(std::size_t i) const {
    if (!stressDependentViscosity())
      return implicitViscosity();
    return std::max(nodeViscosity(i) * (T(1) - nodeDecay(i)),
                    std::numeric_limits<T>::min() * T(1e10));
  }

  // x / sinh(x) for the shear stress of a deviatoric stress tensor, floored
  // at 1e-6 (the viscosity never drops below 1e-6 eta0).
  T eyringFactor(const std::array<T, 9> &deviatoric, T volume) const {
    T j2 = T(0);
    for (unsigned k = 0; k < 9; ++k)
      j2 += deviatoric[k] * deviatoric[k];
    const T shear = std::sqrt(T(0.5) * j2);
    const T kT =
        T(1.380649e-23) * std::max(T(oxidationParameters.temperature), T(1));
    const T x = shear * volume / (T(2) * kT);
    constexpr T floor = T(1e-6);
    if (!(x > T(1e-6)))
      return T(1);
    if (x > T(30))
      return floor;
    return std::max(x / std::sinh(x), floor);
  }

  // Start of a mechanics solve: viscosity from the step-start stress.
  void initOxideViscosity() {
    const T eta0 = T(deformationParameters.viscosity);
    oxideEta_.assign(nodes.size(), eta0);
    for (std::size_t i = 0; i < nodes.size(); ++i)
      if (!nodes[i].nitride)
        oxideEta_[i] =
            eta0 * eyringFactor(startStressAt(nodes[i].index),
                                T(deformationParameters.shearActivationVolume));
  }

  // Start of a mechanics solve: nitride viscosity from the step-start stress.
  void initNitrideViscosity() {
    const T eta0 = T(deformationParameters.nitrideViscosity);
    nitrideEta_.assign(nodes.size(), eta0);
    for (std::size_t i = 0; i < nodes.size(); ++i)
      if (nodes[i].nitride)
        nitrideEta_[i] =
            eta0 *
            eyringFactor(deviatorOf(nitrideStartOf(i)),
                         T(deformationParameters.nitrideShearActivationVolume));
  }

  static std::array<T, 9> deviatorOf(const std::array<T, 9> &stress) {
    auto s = stress;
    const T mean = (stress[0] + stress[4] + stress[8]) / T(3);
    for (unsigned d = 0; d < 3; ++d)
      s[d * 3 + d] -= mean;
    return s;
  }

  // Viscosity consistent with this node's strain rate: the root of
  //   eta = eta0 f(|S_end(eta)|),  S_end = d(eta) S_start + (1 - d(eta)) 2 eta
  //   D'
  // (f = x / sinh x, d = exp(-dt G / eta)).  |S_end| grows with eta and f
  // falls with |S_end|, so the root is unique; bisection in log(eta) between
  // the floor and eta0.  Iterating the stress -> viscosity map instead
  // oscillates for this strongly shear-thinning law.
  T localViscosity(const std::array<T, 9> &start,
                   const std::array<T, 9> &deviatoricRate) const {
    return localViscosity(start, deviatoricRate,
                          T(deformationParameters.viscosity),
                          effectiveStressRelaxationTime(),
                          T(deformationParameters.shearActivationVolume));
  }
  // The nitride's (start: the deviator of its stress).
  T localNitrideViscosity(const std::array<T, 9> &start,
                          const std::array<T, 9> &deviatoricRate) const {
    const T eta0 = T(deformationParameters.nitrideViscosity);
    return localViscosity(
        start, deviatoricRate, eta0, eta0 / nitrideMu(),
        T(deformationParameters.nitrideShearActivationVolume));
  }
  T localViscosity(const std::array<T, 9> &start,
                   const std::array<T, 9> &deviatoricRate, T eta0, T tau0,
                   T volume) const {
    const T dt = T(deformationParameters.stressTimeStep);
    const auto residual = [&](T logEta) {
      const T eta = std::exp(logEta);
      const T decay = (tau0 <= std::numeric_limits<T>::epsilon())
                          ? T(0)
                          : std::exp(-dt / (tau0 * eta / eta0));
      std::array<T, 9> end{};
      for (unsigned k = 0; k < 9; ++k)
        end[k] =
            decay * start[k] + (T(1) - decay) * T(2) * eta * deviatoricRate[k];
      return logEta - std::log(eta0 * eyringFactor(end, volume));
    };
    const T hi = std::log(eta0);
    const T rHi = residual(hi);
    if (rHi <= T(0))
      return eta0; // below the critical stress at eta0
    const T lo = hi + std::log(T(1e-6));
    const T rLo = residual(lo);
    if (rLo >= T(0))
      return std::exp(lo);
    // within 0.5e-10 of a root, as the midpoint of the bisection this
    // replaced (~40 evaluations; ~8 here)
    return std::exp(bracketedRoot(residual, lo, rLo, hi, rHi, T(0.5e-10)));
  }

  // A root of f in [a, b], f(a) and f(b) of opposite signs, within tolerance
  // (Brent: inverse quadratic and secant steps, a bisection step whenever
  // they would not shrink the bracket fast enough; the bracket always
  // holds a sign change).
  template <class Function>
  static T bracketedRoot(const Function &f, T a, T fa, T b, T fb, T tolerance) {
    T c = b, fc = fb, d = b - a, e = d;
    for (int iteration = 0; iteration < 100; ++iteration) {
      if ((fb > T(0) && fc > T(0)) || (fb < T(0) && fc < T(0))) {
        c = a;
        fc = fa;
        d = e = b - a;
      }
      if (std::abs(fc) < std::abs(fb)) {
        a = b;
        b = c;
        c = a;
        fa = fb;
        fb = fc;
        fc = fa;
      }
      const T tol = T(2) * std::numeric_limits<T>::epsilon() * std::abs(b) +
                    T(0.5) * tolerance;
      const T m = T(0.5) * (c - b);
      if (std::abs(m) <= tol || fb == T(0))
        return b;
      if (std::abs(e) >= tol && std::abs(fa) > std::abs(fb)) {
        const T s = fb / fa;
        T p, q;
        if (a == c) {
          p = T(2) * m * s;
          q = T(1) - s;
        } else {
          const T qa = fa / fc, r = fb / fc;
          p = s * (T(2) * m * qa * (qa - r) - (b - a) * (r - T(1)));
          q = (qa - T(1)) * (r - T(1)) * (s - T(1));
        }
        if (p > T(0))
          q = -q;
        else
          p = -p;
        if (T(2) * p <
            std::min(T(3) * m * q - std::abs(tol * q), std::abs(e * q))) {
          e = d;
          d = p / q;
        } else {
          d = m;
          e = m;
        }
      } else {
        d = m;
        e = m;
      }
      a = b;
      fa = fb;
      b += std::abs(d) > tol ? d : (m > T(0) ? tol : -tol);
      fb = f(b);
    }
    return b;
  }

  T effectiveStressRelaxationTime() const {
    if (deformationParameters.stressRelaxationTime > T(0))
      return deformationParameters.stressRelaxationTime;
    if (deformationParameters.shearModulus > std::numeric_limits<T>::epsilon())
      return deformationParameters.viscosity /
             deformationParameters.shearModulus;
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

    const T reactionDistance =
        reactionCrosses ? crossingDistance(reactionInside, reactionOutside)
                        : std::numeric_limits<T>::max();
    const T ambientDistance =
        ambientCrosses ? crossingDistance(ambientInside, ambientOutside)
                       : std::numeric_limits<T>::max();
    const T maskDistance = maskCrosses
                               ? crossingDistance(maskInside, maskOutside)
                               : std::numeric_limits<T>::max();
    if (reactionCrosses && reactionDistance <= ambientDistance &&
        reactionDistance <= maskDistance)
      return {Boundary::REACTION, reactionDistance};
    // Every grid point inside the mask is a nitride node, so a face that
    // leaves the oxide across the mask interface borders gas: free surface.
    // An oxide/gas crossing that lies inside the mask keeps its own distance.
    if (ambientCrosses &&
        ambientEndsInMask(maskInside, maskOutside, ambientDistance))
      return {Boundary::AMBIENT, ambientDistance};
    return {Boundary::AMBIENT, std::min(ambientDistance, maskDistance)};
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

  // The oxide/gas crossing lies inside the mask, or the outer point does.
  bool ambientEndsInMask(T maskInside, T maskOutside, T distance) const {
    return isMaskAtCrossing(maskInside, maskOutside, distance) ||
           (maskInterface != nullptr &&
            static_cast<T>(maskSign) * maskOutside >= T(0));
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
        insidePhi, outsidePhi,
        deformationParameters.minMechanicsBoundaryDistance, gridDelta);
  }

  static constexpr unsigned tensorIndex(unsigned row, unsigned column) {
    return 3 * row + column;
  }
};

/// Monolithic contact: one material's velocity from the shared mechanics
/// solve, as a VelocityField for advecting that material's surface.
template <class T, int D>
class OxidationMaterialVelocity final : public VelocityField<T> {
  SmartPointer<OxidationDeformation<T, D>> field_;
  bool nitride_ = false;

public:
  OxidationMaterialVelocity(SmartPointer<OxidationDeformation<T, D>> field,
                            bool nitride)
      : field_(field), nitride_(nitride) {}

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<OxidationMaterialVelocity>::New(
        std::forward<Args>(args)...);
  }

  Vec3D<T> getVectorVelocity(const Vec3D<T> &coordinate, int /*material*/,
                             const Vec3D<T> & /*normalVector*/,
                             unsigned long /*pointId*/) final {
    return field_->getMaterialVelocity(coordinate, nitride_);
  }

  T getDissipationAlpha(int direction, int material,
                        const Vec3D<T> &centralDifferences) final {
    return field_->getDissipationAlpha(direction, material, centralDifferences);
  }
};

} // namespace viennals
