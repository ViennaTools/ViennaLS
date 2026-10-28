#pragma once

#include <lsAdvect.hpp>
#include <lsAdvectIntegrationSchemes.hpp>
#include <lsBooleanOperation.hpp>
#include <lsExpand.hpp>
#include <lsInterior.hpp>
#include <lsOxidationModel.hpp>
#include <lsReduce.hpp>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>

#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vcTimer.hpp>

namespace viennals {

using namespace viennacore;

template <class T> struct LOCOSConservationDiagnostics {
  T siliconRecession = 0.;
  T ambientLift = 0.;
  T expectedAmbientLift = 0.;
  T ambientLiftRatio = 0.;
  T relativeError = 0.;
  unsigned samples = 0;
};

template <class T>
std::optional<T> findLOCOSInterfaceY(SmartPointer<Domain<T, 2>> levelSet,
                                     viennahrle::IndexType i,
                                     viennahrle::IndexType jMin,
                                     viennahrle::IndexType jMax) {
  using ConstIterator =
      viennahrle::ConstSparseIterator<typename Domain<T, 2>::DomainType>;
  ConstIterator it(levelSet->getDomain());

  const T gridDelta = levelSet->getGrid().getGridDelta();
  viennahrle::Index<2> previousIndex{i, jMin};
  it.goToIndices(previousIndex);
  T previous = it.getValue();

  for (auto j = jMin + 1; j <= jMax; ++j) {
    viennahrle::Index<2> currentIndex{i, j};
    it.goToIndices(currentIndex);
    const T current = it.getValue();
    if ((previous <= 0. && current >= 0.) ||
        (previous >= 0. && current <= 0.)) {
      const T denominator = std::abs(previous) + std::abs(current);
      const T fraction = denominator > std::numeric_limits<T>::epsilon()
                             ? std::abs(previous) / denominator
                             : 0.;
      return (static_cast<T>(j - 1) + fraction) * gridDelta;
    }
    previous = current;
  }

  return std::nullopt;
}

/// Measure the open-window volume balance after one 2D LOCOS step.
///
/// The diagnostic samples vertical columns, finds the Si/SiO2 and
/// oxide/ambient zero crossings, and compares ambient lift against the
/// Deal-Grove volume expansion expectation:
///   ambient lift = silicon consumed * (expansionCoefficient - 1)
template <class T>
LOCOSConservationDiagnostics<T> computeLOCOSOpenWindowConservation(
    SmartPointer<Domain<T, 2>> siInitial, SmartPointer<Domain<T, 2>> siAfter,
    SmartPointer<Domain<T, 2>> ambientInitial,
    SmartPointer<Domain<T, 2>> ambientAfter, T xMin, T xMax,
    viennahrle::IndexType jMin, viennahrle::IndexType jMax,
    T expansionCoefficient) {
  LOCOSConservationDiagnostics<T> result;
  const T gridDelta = siInitial->getGrid().getGridDelta();
  const auto iMin =
      static_cast<viennahrle::IndexType>(std::ceil(xMin / gridDelta));
  const auto iMax =
      static_cast<viennahrle::IndexType>(std::floor(xMax / gridDelta));

  for (auto i = iMin; i <= iMax; ++i) {
    const auto si0 = findLOCOSInterfaceY<T>(siInitial, i, jMin, jMax);
    const auto si1 = findLOCOSInterfaceY<T>(siAfter, i, jMin, jMax);
    const auto amb0 = findLOCOSInterfaceY<T>(ambientInitial, i, jMin, jMax);
    const auto amb1 = findLOCOSInterfaceY<T>(ambientAfter, i, jMin, jMax);
    if (!si0 || !si1 || !amb0 || !amb1)
      continue;

    result.siliconRecession += std::max(*si0 - *si1, T(0.)) * gridDelta;
    result.ambientLift += std::max(*amb1 - *amb0, T(0.)) * gridDelta;
    ++result.samples;
  }

  result.expectedAmbientLift =
      result.siliconRecession * (expansionCoefficient - T(1.));
  if (result.siliconRecession > 0.)
    result.ambientLiftRatio = result.ambientLift / result.siliconRecession;
  if (result.expectedAmbientLift > 0.)
    result.relativeError =
        std::abs(result.ambientLift - result.expectedAmbientLift) /
        result.expectedAmbientLift;
  return result;
}

/// Unified oxidation time-step orchestrator.
///
/// Couples diffusion and viscoelastic deformation into a single CFL-limited
/// step. When a mask interface is provided the solver operates in LOCOS mode;
/// without one it performs standard planar or step oxidation.
///
/// Standard mode (no mask):
///   1. Diffusion+deformation coupled solve (OxidationModel)
///   2. Ambient interface advection (deformation velocity)
///   3. Si/SiO2 reaction interface advection (diffusion velocity)
///   4. Interior fill for warm-start persistence
///
/// LOCOS mode (mask provided):
///   1. Diffusion+deformation coupled solve (OxidationModel); the nitride
///      joins the deformation solve as a linear-elastic body (monolithic
///      oxide + nitride mechanics, OxidationMaskParameters)
///   2. Pre-advection boolean clip (ambientInterface \ maskInterface)
///   3. Ambient interface advection (oxide velocity)
///   4. Si/SiO2 reaction interface advection (diffusion velocity)
///   5. Mask advection (nitride velocity), re-sampled from the mask's shape
///      at the start of the call and the nitride displacement
///      (OxidationMaskTracker: keeps its corners sharp)
///   6. Interior fills + post-advection boolean clip
/// With setMaskBarrierOnly(true) the mask is a diffusion barrier only: it is
/// withheld from the mechanics and carried along by the oxide.
///
/// Usage (standard):
/// @code
///   auto ox = ls::Oxidation<double, 2>::New(si, ambient);
///   ox->setOxidationParameters(oxParams);
///   ox->setDeformationParameters(defParams);
///   ox->apply(0.1);  // advection time in hours
/// @endcode
///
/// Usage (LOCOS):
/// @code
///   auto ox = ls::Oxidation<double, 2>::New(si, ambient, mask);
///   ox->setOxidationParameters(oxParams);
///   ox->setDeformationParameters(defParams);
///   ox->setMaskParameters(maskParams);
///   ox->setSolveBounds(minIdx, maxIdx);
///   ox->applyCFLLimited(0.35, 0.499);
///   ox->getDiffusionField()->getConcentration(pt);  // diagnostics
/// @endcode
/// Lagrangian nitride surface.  The nitride is an elastic solid that moves
/// almost rigidly, by up to hundreds of nm, and level-set advection rounds
/// its convex corners by a numerical diffusion that grows with the distance
/// travelled, for every scheme order: over a 0.1 um lift on a 5 nm grid the
/// side wall 2.5 nm from a corner moved in by 13 nm (Engquist-Osher 1st
/// order) or 6.5 nm (WENO5 + RK3), and a 50 nm strip lost its contact with
/// the oxide.  After each advection the mask is therefore re-sampled from
/// the shape it had when the tracker was created and the accumulated
/// displacement u of the nitride, phi(x) = phi_0(x - u(x)) (D-linear in
/// phi_0).  u follows the material: u(x) = u_old(x - v dt) + v dt with the
/// velocity v the mask was advected with (accumulating v at a fixed point
/// instead turns a rotation of angle a by an error ~ a^2 r / 2).  The
/// advection still moves the narrow band and carries the point data; the
/// re-sampled band is normalised like an advected one.  A tracker lives
/// for one Oxidation object, i.e. one process call: phi_0 is sampled once
/// per call (over a 0.1 um lift the corner moved by 0.9 nm in one call,
/// 2.5 nm in six; re-sampled every step, 17 nm).
template <class T, int D> class OxidationMaskTracker {
  using IndexType = viennahrle::Index<D>;
  using ConstSparseIterator =
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>;
  using DisplacementMap =
      std::unordered_map<IndexType, Vec3D<T>, detail::IndexTypeHasher<D>>;

  // phi_0 is expanded to |phi| <= 4.5: a departure point lies as far from
  // the reference surface as its band point from the current one.
  static constexpr int referenceWidth = 9;
  static constexpr T referenceLimit = T(5);

  const Domain<T, D> *mask_;
  SmartPointer<Domain<T, D>> reference_;
  DisplacementMap displacement_;
  T gridDelta_;

public:
  explicit OxidationMaskTracker(SmartPointer<Domain<T, D>> mask)
      : mask_(mask.get()), reference_(SmartPointer<Domain<T, D>>::New(mask)),
        gridDelta_(mask->getGrid().getGridDelta()) {
    Expand<T, D>(reference_, referenceWidth).apply();
  }

  /// False for another mask, or after the mask was resampled to another grid.
  bool tracks(const SmartPointer<Domain<T, D>> &mask) const {
    return mask.get() == mask_ && mask->getGrid().getGridDelta() == gridDelta_;
  }

  /// Re-sample the just advected mask: dt is the advection time, velocity
  /// (Vec3D<T>(const Vec3D<T> &coordinate)) the one it was advected with.
  template <class VelocityFunction>
  void correct(SmartPointer<Domain<T, D>> mask, VelocityFunction &&velocity,
               T dt) {
    const auto &grid = mask->getGrid();
    ConstSparseIterator referenceIt(reference_->getDomain());
    auto referenceValue = [&](const IndexType &index) {
      referenceIt.goToIndices(grid.globalIndices2LocalIndices(index));
      return std::clamp(referenceIt.getValue(), -referenceLimit,
                        referenceLimit);
    };

    struct Update {
      unsigned segment;
      std::size_t localId;
      T value;
    };
    std::vector<Update> updates;
    DisplacementMap next;
    auto &domain = mask->getDomain();
    for (ConstSparseIterator it(domain); !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const IndexType index = it.getStartIndices();
      Vec3D<T> coordinate{T(0), T(0), T(0)};
      for (unsigned d = 0; d < D; ++d)
        coordinate[d] = static_cast<T>(index[d]) * gridDelta_;
      const Vec3D<T> v = velocity(coordinate);
      T departure[D];
      for (unsigned d = 0; d < D; ++d)
        departure[d] = static_cast<T>(index[d]) - v[d] * dt / gridDelta_;
      Vec3D<T> u = departureDisplacement(departure);
      for (unsigned d = 0; d < D; ++d)
        u[d] += v[d] * dt;

      // phi_0 at the departure point x - u, D-linear.
      IndexType lower;
      T fraction[D];
      for (unsigned d = 0; d < D; ++d) {
        const T c = static_cast<T>(index[d]) - u[d] / gridDelta_;
        const T cFloor = std::floor(c);
        lower[d] = static_cast<viennahrle::IndexType>(cFloor);
        fraction[d] = c - cFloor;
      }
      T value = T(0);
      for (int corner = 0; corner < (1 << D); ++corner) {
        IndexType vertex = lower;
        T weight = T(1);
        for (unsigned d = 0; d < D; ++d) {
          const bool up = (corner >> d) & 1;
          vertex[d] += up ? 1 : 0;
          weight *= up ? fraction[d] : T(1) - fraction[d];
        }
        if (weight > T(0))
          value += weight * referenceValue(vertex);
      }
      const unsigned segment = it.getSegmentId();
      updates.push_back(
          {segment, it.getPointId() - domain.getPointIdOffset(segment), value});
      next.emplace(index, u);
    }
    // The re-sampled values are no ViennaLS narrow band: a stretched or bent
    // nitride gives neighbours like -0.504 / +0.527 across its surface,
    // which the level-set algorithms reject (they walked off the HRLE and
    // crashed).  Normalise as Advect does: such a pair (both beyond 0.5) is
    // clamped to +-0.5, the outer layers are rebuilt from the active points.
    // Only such pairs: clamping e.g. -0.88 next to +0.12 moves the surface
    // by 0.07 cells every step, which lifted the nitride off the oxide until
    // the mechanics lost the contact.
    std::unordered_map<IndexType, T, detail::IndexTypeHasher<D>> values;
    values.reserve(updates.size());
    std::size_t k = 0;
    for (ConstSparseIterator it(domain); !it.isFinished(); ++it)
      if (it.isDefined())
        values.emplace(it.getStartIndices(), updates[k++].value);
    k = 0;
    for (ConstSparseIterator it(domain); !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      auto &update = updates[k++];
      const IndexType index = it.getStartIndices();
      for (unsigned d = 0; d < D && std::abs(update.value) > T(0.5); ++d)
        for (const int side : {-1, 1}) {
          IndexType neighbor = index;
          neighbor[d] += side;
          const auto found = values.find(neighbor);
          if (found != values.end() && std::abs(found->second) > T(0.5) &&
              std::signbit(found->second) != std::signbit(update.value))
            update.value = std::clamp(update.value, T(-0.5), T(0.5));
        }
    }
    for (const auto &update : updates)
      domain.getDomainSegment(update.segment).definedValues[update.localId] =
          update.value;
    const int width = std::max(mask->getLevelSetWidth(), 2);
    Reduce<T, D>(mask, 1).apply();
    Expand<T, D>(mask, width).apply();
    displacement_ = std::move(next);
  }

private:
  // u_old at a departure point (grid units): D-linear over the stored band
  // points, skipping corners outside the old band; else the nearest stored
  // point within 2 cells (the band moves less than a cell per step and u
  // varies on the scale of the body); zero before the first step.
  Vec3D<T> departureDisplacement(const T (&position)[D]) const {
    IndexType lower;
    T fraction[D];
    for (unsigned d = 0; d < D; ++d) {
      const T cFloor = std::floor(position[d]);
      lower[d] = static_cast<viennahrle::IndexType>(cFloor);
      fraction[d] = position[d] - cFloor;
    }
    Vec3D<T> sum{T(0), T(0), T(0)};
    T weightSum = T(0);
    for (int corner = 0; corner < (1 << D); ++corner) {
      IndexType vertex = lower;
      T weight = T(1);
      for (unsigned d = 0; d < D; ++d) {
        const bool up = (corner >> d) & 1;
        vertex[d] += up ? 1 : 0;
        weight *= up ? fraction[d] : T(1) - fraction[d];
      }
      if (weight <= T(0))
        continue;
      const auto found = displacement_.find(vertex);
      if (found == displacement_.end())
        continue;
      for (unsigned d = 0; d < D; ++d)
        sum[d] += weight * found->second[d];
      weightSum += weight;
    }
    if (weightSum > T(1e-6)) {
      for (unsigned d = 0; d < D; ++d)
        sum[d] /= weightSum;
      return sum;
    }
    IndexType nearest;
    for (unsigned d = 0; d < D; ++d)
      nearest[d] =
          static_cast<viennahrle::IndexType>(std::llround(position[d]));
    return nearestDisplacement(nearest);
  }

  Vec3D<T> nearestDisplacement(const IndexType &index) const {
    const Vec3D<T> *best = nullptr;
    int bestDistance = std::numeric_limits<int>::max();
    IndexType offset;
    for (unsigned d = 0; d < D; ++d)
      offset[d] = -2;
    while (true) {
      int distance = 0;
      IndexType neighbor = index;
      for (unsigned d = 0; d < D; ++d) {
        neighbor[d] += offset[d];
        distance += offset[d] * offset[d];
      }
      if (distance < bestDistance) {
        const auto it = displacement_.find(neighbor);
        if (it != displacement_.end()) {
          best = &it->second;
          bestDistance = distance;
        }
      }
      unsigned d = 0;
      for (; d < D; ++d) {
        if (offset[d] < 2) {
          ++offset[d];
          break;
        }
        offset[d] = -2;
      }
      if (d == D)
        break;
    }
    return best != nullptr ? *best : Vec3D<T>{T(0), T(0), T(0)};
  }
};

template <class T, int D> class Oxidation {
  using IndexType = viennahrle::Index<D>;

  SmartPointer<Domain<T, D>> siInterface = nullptr;
  SmartPointer<Domain<T, D>> ambientInterface = nullptr;
  SmartPointer<Domain<T, D>> maskInterface = nullptr;

  OxidationParameters oxidationParams;
  OxidationDeformationParameters deformationParams;
  OxidationCouplingParameters couplingParams;
  OxidationMaskParameters maskParams;

  SpatialSchemeEnum spatialScheme = SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER;
  TemporalSchemeEnum temporalScheme = TemporalSchemeEnum::FORWARD_EULER;
  static constexpr int maskInteriorSign = -1;
  /// Treat the mask as a DIFFUSION BARRIER ONLY: it reduces the oxidant
  /// reaching the silicon but is mechanically transparent — the oxide expands
  /// as if at a free surface and carries the mask along on top.  Exact for a
  /// full-width mask (uniform growth = pure translation); the mode for flat
  /// Fig. 4-style validation.
  bool maskBarrierOnly_ = false;

  IndexType diffusionMinIndex{};
  IndexType diffusionMaxIndex{};
  bool diffusionBoundsSet = false;

  // Populated by apply(); available for diagnostics afterwards.
  SmartPointer<OxidationDiffusion<T, D>> diffusionField;
  SmartPointer<OxidationDeformation<T, D>> deformationField;
  T lastMaxVelocity_ = T(0);

  GpuMode gpuMode_ = GpuMode::Cpu;
  std::unordered_map<std::size_t, T> concentrationCache_;
  // Stress feedback carried between solves (see
  // OxidationDiffusion::setFeedbackCache).
  typename OxidationDiffusion<T, D>::PressureLookup pressureCache_;
  typename OxidationDiffusion<T, D>::StressLookup stressCache_;
  std::optional<OxidationMaskTracker<T, D>> maskTracker_;

public:
  const std::unordered_map<std::size_t, T> &getConcentrationCache() const {
    return concentrationCache_;
  }
  void setConcentrationCache(std::unordered_map<std::size_t, T> cache) {
    concentrationCache_ = std::move(cache);
  }

  Oxidation() = default;

  Oxidation(SmartPointer<Domain<T, D>> passedSiInterface,
            SmartPointer<Domain<T, D>> passedAmbientInterface,
            SmartPointer<Domain<T, D>> passedMaskInterface = nullptr)
      : siInterface(passedSiInterface),
        ambientInterface(passedAmbientInterface),
        maskInterface(passedMaskInterface) {}

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<Oxidation>::New(std::forward<Args>(args)...);
  }

  void setGpuMode(GpuMode mode) { gpuMode_ = mode; }

  void setSiInterface(SmartPointer<Domain<T, D>> si) { siInterface = si; }
  void setAmbientInterface(SmartPointer<Domain<T, D>> ambient) {
    ambientInterface = ambient;
  }
  void setMaskInterface(SmartPointer<Domain<T, D>> mask) {
    maskInterface = mask;
  }

  void setOxidationParameters(OxidationParameters params) {
    oxidationParams = params;
  }
  void setDeformationParameters(OxidationDeformationParameters params) {
    deformationParams = params;
  }
  void setCouplingParameters(OxidationCouplingParameters params) {
    couplingParams = params;
  }
  void setMaskParameters(OxidationMaskParameters params) {
    maskParams = params;
  }

  /// Set the spatial integration scheme for all advections.
  void setSpatialScheme(SpatialSchemeEnum scheme) { spatialScheme = scheme; }

  void setMaskBarrierOnly(bool enable) { maskBarrierOnly_ = enable; }

  /// Set the temporal integration scheme for all advections.
  void setTemporalScheme(TemporalSchemeEnum scheme) { temporalScheme = scheme; }

  /// Set the Cartesian index bounding box for the diffusion and deformation
  /// solves. If not set, bounds are auto-computed from the level-set narrow
  /// band.
  void setSolveBounds(const IndexType &minIndex, const IndexType &maxIndex) {
    diffusionMinIndex = minIndex;
    diffusionMaxIndex = maxIndex;
    diffusionBoundsSet = true;
  }

  /// Return the diffusion field populated by the most recent apply() call.
  SmartPointer<OxidationDiffusion<T, D>> getDiffusionField() const {
    return diffusionField;
  }

  /// Return the deformation field populated by the most recent apply() call.
  SmartPointer<OxidationDeformation<T, D>> getDeformationField() const {
    return deformationField;
  }

  /// Maximum interface velocity (µm/hr) from the most recent CFL-limited step.
  T getLastMaxVelocity() const { return lastMaxVelocity_; }

  /// Execute one oxidation time step of duration advectionTime.
  void apply(T advectionTime) { applyImpl(advectionTime, std::nullopt); }

  /// Execute one CFL-limited oxidation step; returns the actual time advanced.
  T applyCFLLimited(T requestedTime, T cflFactor) {
    return applyImpl(requestedTime, std::clamp(cflFactor, T(1e-3), T(0.499)));
  }

private:
  T applyImpl(T requestedTime, std::optional<T> cflFactor) {
    if (siInterface == nullptr || ambientInterface == nullptr) {
      Logger::getInstance()
          .addError("Oxidation: Si or ambient interface is null.")
          .print();
      return T(0);
    }

    if (requestedTime <= T(0))
      return T(0);

    Timer<> tStep;
    tStep.start();

    const bool hasMask = (maskInterface != nullptr);
    // The nitride is part of the oxide mechanics solve (OxidationDeformation):
    // there is no separate mask solve and no coupling loop, and each surface
    // moves with its own material's solved velocity.
    const bool monolithicMask = hasMask && !maskBarrierOnly_;
    const std::string prefix = hasMask ? "LOCOS" : "Oxidation";

    VIENNACORE_LOG_INFO(prefix + ": starting time step, requested_dt=" +
                        std::to_string(requestedTime) + " hr");

    const auto baseConcentrationCache = concentrationCache_;
    const auto basePressureCache = pressureCache_;
    const auto baseStressCache = stressCache_;
    std::string lastFieldFailureReason;

    auto solveFields = [&](T stressTimeStep) -> bool {
      auto rejectSolve = [&](const std::string &reason) {
        lastFieldFailureReason = reason;
        return false;
      };

      auto validateCoupledModel = [&](const SmartPointer<OxidationModel<T, D>>
                                          &model) {
        if (!model->hasConverged()) {
          const auto reason = model->getFailureReason();
          return rejectSolve(
              reason.empty()
                  ? "pressure-concentration coupling failed (residual=" +
                        std::to_string(model->getResidual()) + ", tolerance=" +
                        std::to_string(couplingParams.tolerance) + ")"
                  : reason);
        }
        if (!diffusionField->lastSolveConverged() ||
            !diffusionField->hasFiniteConcentrationField()) {
          return rejectSolve(
              "diffusion solve failed (residual=" +
              std::to_string(diffusionField->getNormalizedResidual()) +
              ", tolerance=" + std::to_string(oxidationParams.tolerance) + ")");
        }
        if (!deformationField->lastSolveConverged() ||
            !deformationField->hasFiniteSolution()) {
          return rejectSolve(
              "deformation solve failed (mechanics=" +
              std::to_string(deformationField->getResidual()) +
              ", continuity=" +
              std::to_string(deformationField->getLastContinuityResidual()) +
              ", stokes=" +
              std::to_string(deformationField->getLastStokesResidual()) + ")");
        }
        return true;
      };

      lastFieldFailureReason.clear();
      auto stepDeformationParams = deformationParams;
      stepDeformationParams.stressTimeStep = stressTimeStep;
      if (monolithicMask) {
        stepDeformationParams.nitrideYoungModulus = maskParams.youngModulus;
        stepDeformationParams.nitridePoissonRatio = maskParams.poissonRatio;
        stepDeformationParams.nitrideViscosity = maskParams.viscosity;
        stepDeformationParams.nitrideShearActivationVolume =
            maskParams.shearActivationVolume;
      }

      // --- Coupled diffusion + deformation solve ---

      diffusionField = OxidationDiffusion<T, D>::New(
          siInterface, ambientInterface, oxidationParams);
      diffusionField->setConcentrationCache(baseConcentrationCache);
      diffusionField->setFeedbackCache(basePressureCache, baseStressCache);
      diffusionField->setGpuMode(gpuMode_);
      if (hasMask)
        diffusionField->setMaskInterface(maskInterface, maskInteriorSign);

      deformationField = OxidationDeformation<T, D>::New(
          siInterface, ambientInterface, diffusionField, oxidationParams,
          stepDeformationParams);
      deformationField->setGpuMode(gpuMode_);
      // Barrier-only: the mask is withheld from the mechanics so the oxide
      // sees a free surface and expands outward.
      if (hasMask && !maskBarrierOnly_)
        deformationField->setMaskInterface(maskInterface, maskInteriorSign);

      auto coupledModel = OxidationModel<T, D>::New(
          diffusionField, deformationField, couplingParams);
      if (diffusionBoundsSet)
        coupledModel->setSolveBounds(diffusionMinIndex, diffusionMaxIndex);
      VIENNACORE_LOG_DEBUG(
          prefix + ": solving coupled diffusion/deformation field for dt=" +
          std::to_string(stressTimeStep) + " hr");
      Timer<> tCoupled;
      tCoupled.start();
      coupledModel->apply();
      tCoupled.finish();
      if (Logger::hasTiming())
        Logger::getInstance().addTiming("  coupled(iter=1)", tCoupled).print();
      VIENNACORE_LOG_DEBUG(prefix +
                           ": coupled diffusion/deformation solve complete");
      return validateCoupledModel(coupledModel);
    };

    auto makeAmbientVelocity = [&]() -> SmartPointer<VelocityField<T>> {
      // Barrier-only: the oxide's free surface must NOT be constrained by
      // the mask -- it expands freely and carries the mask along.
      // Monolithic: oxide and nitride share one solve; each surface moves
      // with its own material's velocity.
      if (monolithicMask)
        return OxidationMaterialVelocity<T, D>::New(deformationField, false);
      return deformationField;
    };

    T advectionTime = requestedTime;
    SmartPointer<VelocityField<T>> ambientVelocity;

    auto computeMaxVelocity =
        [&](const SmartPointer<VelocityField<T>> &passedAmbientVelocity) {
          T maxVelocity = diffusionField->getDissipationAlpha(0, -1, {});
          // The mask moves with the oxide velocity field (its own nitride
          // nodes, or carried in barrier-only mode): covered above.
          for (unsigned d = 0; d < D; ++d)
            maxVelocity =
                std::max(maxVelocity,
                         passedAmbientVelocity->getDissipationAlpha(d, -1, {}));
          return maxVelocity;
        };

    auto cflLimitedTime = [&](T trialTime, T maxVelocity) {
      if (maxVelocity <= std::numeric_limits<T>::epsilon())
        return trialTime;
      const T gridDelta = siInterface->getGrid().getGridDelta();
      return std::min(trialTime, (*cflFactor) * gridDelta / maxVelocity);
    };

    if (cflFactor) {
      T trialTime = requestedTime;
      const T minTrialTime = std::max(
          requestedTime * T(1e-10), std::numeric_limits<T>::epsilon() * T(100));
      bool accepted = false;
      for (unsigned attempt = 0; attempt < 16; ++attempt) {
        const bool predictorConverged = solveFields(trialTime);
        if (!predictorConverged) {
          const T nextTrial = trialTime * T(0.5);
          VIENNACORE_LOG_INFO(
              prefix + ": rejecting non-converged coupled predictor (" +
              lastFieldFailureReason + "), retrying with requested_dt=" +
              std::to_string(nextTrial) + " hr");
          trialTime = nextTrial;
          if (trialTime < minTrialTime)
            break;
          continue;
        }

        ambientVelocity = makeAmbientVelocity();
        T maxVelocity = computeMaxVelocity(ambientVelocity);
        if (!std::isfinite(maxVelocity))
          VIENNACORE_LOG_ERROR(prefix + ": non-finite CFL velocity estimate.");

        advectionTime = cflLimitedTime(trialTime, maxVelocity);
        VIENNACORE_LOG_INFO(
            prefix +
            ": CFL decision requested_dt=" + std::to_string(trialTime) +
            " hr, actual_dt=" + std::to_string(advectionTime) +
            " hr, max_velocity=" + std::to_string(maxVelocity) + " um/hr");

        if (advectionTime < trialTime * (T(1) - T(1e-8))) {
          const bool finalConverged = solveFields(advectionTime);
          if (!finalConverged) {
            const T nextTrial = advectionTime * T(0.5);
            VIENNACORE_LOG_INFO(
                prefix + ": rejecting non-converged CFL re-solve (" +
                lastFieldFailureReason + "), retrying with requested_dt=" +
                std::to_string(nextTrial) + " hr");
            trialTime = nextTrial;
            if (trialTime < minTrialTime)
              break;
            continue;
          }
          ambientVelocity = makeAmbientVelocity();
          maxVelocity = computeMaxVelocity(ambientVelocity);
          if (!std::isfinite(maxVelocity))
            VIENNACORE_LOG_ERROR(prefix +
                                 ": non-finite accepted CFL velocity.");

          const T verifiedTime = cflLimitedTime(advectionTime, maxVelocity);
          if (verifiedTime < advectionTime * (T(1) - T(1e-8))) {
            VIENNACORE_LOG_INFO(prefix +
                                ": rejecting CFL re-solve because accepted "
                                "velocity requires requested_dt=" +
                                std::to_string(verifiedTime) + " hr");
            trialTime = verifiedTime;
            if (trialTime < minTrialTime)
              break;
            continue;
          }
        }

        lastMaxVelocity_ = computeMaxVelocity(ambientVelocity);
        accepted = true;
        break;
      }
      if (!accepted)
        VIENNACORE_LOG_ERROR(prefix +
                             ": unable to find a converged CFL-limited step" +
                             (lastFieldFailureReason.empty()
                                  ? std::string(".")
                                  : std::string(" (last failure: ") +
                                        lastFieldFailureReason + ")."));
    } else {
      const bool fieldsConverged = solveFields(requestedTime);
      if (!fieldsConverged)
        VIENNACORE_LOG_ERROR(
            prefix + ": coupled solve failed" +
            (lastFieldFailureReason.empty()
                 ? std::string(".")
                 : std::string(" (") + lastFieldFailureReason + ")."));
      ambientVelocity = makeAmbientVelocity();
    }

    concentrationCache_ = diffusionField->getConcentrationCache();
    pressureCache_ = diffusionField->getPressureCache();
    stressCache_ = diffusionField->getStressCache();

    diffusionField->markSolved();

    // VOLTRACE: per-stage ambient-volume audit for the under-mask volume
    // leak.  Smoothed-Heaviside sum over the HRLE: band cells contribute
    // clamp(0.5 - phi, 0, 1), negative background runs contribute fully.
    // Only DIFFERENCES between stages are meaningful.
    auto ambientVolume = [&]() -> T {
      const T delta = ambientInterface->getGrid().getGridDelta();
      T cellVol = T(1);
      for (unsigned d = 0; d < D; ++d)
        cellVol *= delta;
      T vol = T(0);
      for (viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>
               it(ambientInterface->getDomain());
           !it.isFinished(); it.next()) {
        T cells = T(1);
        for (unsigned d = 0; d < D; ++d)
          cells *= T(it.getEndIndices(static_cast<int>(d)) -
                     it.getStartIndices(static_cast<int>(d)) + 1);
        if (it.isDefined())
          vol += std::min(std::max(T(0.5) - it.getValue(), T(0)), T(1)) *
                 cellVol * cells;
        else if (it.getValue() < T(0))
          vol += cellVol * cells;
      }
      return vol;
    };
    const T volBefore = ambientVolume();

    // Pre-advection clip: keep oxide outside the mask body (LOCOS only).
    if (hasMask)
      BooleanOperation<T, D>(ambientInterface, maskInterface,
                             BooleanOperationEnum::RELATIVE_COMPLEMENT)
          .apply();
    const T volAfterPreclip = ambientVolume();

    diffusionField->writePersistentFields();
    deformationField->writeFieldsToLevelSet();

    auto advect = [&](SmartPointer<Domain<T, D>> levelSet,
                      SmartPointer<VelocityField<T>> velocityField) {
      Advect<T, D> adv;
      adv.insertNextLevelSet(levelSet);
      adv.setVelocityField(velocityField);
      adv.setSpatialScheme(spatialScheme);
      adv.setTemporalScheme(temporalScheme);
      adv.setAdvectionTime(advectionTime);
      adv.apply();
    };

    Timer<> tAdvect;
    tAdvect.start();
    advect(ambientInterface, ambientVelocity);
    advect(siInterface, diffusionField);
    if (hasMask && maskBarrierOnly_)
      // Carried by the oxide: the contact face samples the oxide velocity in
      // place; OxidationMaskCarry extends that same velocity through the mask
      // body so the whole slab translates together.
      advect(maskInterface,
             OxidationMaskCarry<T, D>::New(deformationField, maskInterface));
    else if (monolithicMask) {
      if (!maskTracker_ || !maskTracker_->tracks(maskInterface))
        maskTracker_.emplace(maskInterface);
      advect(maskInterface,
             OxidationMaterialVelocity<T, D>::New(deformationField, true));
      maskTracker_->correct(
          maskInterface,
          [&](const Vec3D<T> &x) {
            return deformationField->getMaterialVelocity(x, true);
          },
          advectionTime);
    }
    tAdvect.finish();
    VIENNACORE_LOG_TIMING(std::string("  advection(") + (hasMask ? "3" : "2") +
                              " surfaces)",
                          tAdvect);
    const T volAfterAdvect = ambientVolume();

    // Post-advection clip: remove oxide that grew into the mask (LOCOS only).
    // Mask gets Interior fill first so the BooleanOp has accurate φ_mask values
    // at points inside the mask region.
    if (hasMask) {
      Interior<T, D>(maskInterface).apply();
      BooleanOperation<T, D>(ambientInterface, maskInterface,
                             BooleanOperationEnum::RELATIVE_COMPLEMENT)
          .apply();
    }
    const T volAfterPostclip = ambientVolume();
    {
      Interior<T, D> fill(ambientInterface);
      fill.setGuide(siInterface); // stop fill at Si surface
      fill.apply();
    }
    {
      const T volAfterFill = ambientVolume();
      const T toNm = std::pow(T(1000), static_cast<int>(D));
      // vol_start/vol_end are ABSOLUTE (nm^D): vol_start(N+1) != vol_end(N)
      // exposes volume changes between substeps (outside this pipeline).
      VIENNACORE_LOG_INFO(
          prefix + ": VOLTRACE dV_preclip=" +
          std::to_string((volAfterPreclip - volBefore) * toNm) + " dV_advect=" +
          std::to_string((volAfterAdvect - volAfterPreclip) * toNm) +
          " dV_postclip=" +
          std::to_string((volAfterPostclip - volAfterAdvect) * toNm) +
          " dV_fill=" +
          std::to_string((volAfterFill - volAfterPostclip) * toNm) +
          " vol_start=" + std::to_string(volBefore * toNm) +
          " vol_end=" + std::to_string(volAfterFill * toNm) + " nm^D");
    }

    // Re-write persistent fields (concentration, pressure) now that the HRLE
    // has interior points.  The first writePersistentFields() above only
    // covered the narrow band; after lsAdvect the new HRLE is again a narrow
    // band stripped of interior data.  By re-writing here we ensure that the
    // next outer step's buildNodes() can warm-start ALL oxide nodes — not just
    // the surface band — when it reads back from pointData.
    diffusionField->writePersistentFields();
    deformationField->writeFieldsToLevelSet();

    VIENNACORE_LOG_INFO(prefix + ": time step complete, actual_dt=" +
                        std::to_string(advectionTime) + " hr");

    tStep.finish();
    VIENNACORE_LOG_TIMING("── step total", tStep);

    return advectionTime;
  }
};

} // namespace viennals
