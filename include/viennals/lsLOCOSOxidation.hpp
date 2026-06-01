#pragma once

#include <lsAdvect.hpp>
#include <lsAdvectIntegrationSchemes.hpp>
#include <lsBooleanOperation.hpp>
#include <lsInterior.hpp>
#include <lsOxidationModel.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>

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
std::optional<T>
findLOCOSInterfaceY(SmartPointer<Domain<T, 2>> levelSet,
                    viennahrle::IndexType i, viennahrle::IndexType jMin,
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

/// Wrapper that executes one complete LOCOS oxidation time step.
///
/// A LOCOS step couples three physical solves and three level-set advections:
///   1. Diffusion+deformation coupled solve (OxidationModel)
///   2. Nitride mask bending solve (OxidationMaskBending)
///   3. Pre-advection boolean clip (ambientInterface \ maskInterface)
///   4. Ambient interface advection (constrained: mask-coupled under nitride)
///   5. Si/SiO2 reaction interface advection (diffusion velocity)
///   6. Mask advection (mask bending velocity)
///   7. Post-advection boolean clip (ambientInterface \ maskInterface)
///
/// The two boolean clips are mandatory: they keep the ambient and mask level
/// sets from overlapping, which would corrupt subsequent velocity evaluations.
/// This class makes both clips structural rather than hidden in user code.
///
/// LOCOS assumes the mask level set follows the usual ViennaLS solid
/// convention: negative values are inside the nitride. The wrapper applies that
/// convention internally for the diffusion, deformation, constrained ambient,
/// and mask bending fields.
///
/// Usage:
/// @code
///   auto locos = ls::LOCOSOxidation<double, 2>::New(si, ambient, mask);
///   locos->setOxidationParameters(oxParams);
///   locos->setDeformationParameters(defParams);
///   locos->setMaskParameters(maskParams);
///   locos->setSolveBounds(minIdx, maxIdx);
///   locos->setMaskBendingBounds(maskMinIdx, maskMaxIdx);
///   locos->apply(0.35);  // advection time in hours
///   locos->getDiffusionField()->getConcentration(pt);  // diagnostics
/// @endcode
template <class T, int D> class LOCOSOxidation {
  using IndexType = viennahrle::Index<D>;

  SmartPointer<Domain<T, D>> siInterface = nullptr;
  SmartPointer<Domain<T, D>> ambientInterface = nullptr;
  SmartPointer<Domain<T, D>> maskInterface = nullptr;

  OxidationParameters<T> oxidationParams;
  OxidationDeformationParameters<T> deformationParams;
  OxidationCouplingParameters<T> couplingParams;
  OxidationMaskParameters<T> maskParams;

  SpatialSchemeEnum spatialScheme =
      SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER;
  TemporalSchemeEnum temporalScheme =
      TemporalSchemeEnum::FORWARD_EULER;
  static constexpr int maskInteriorSign = -1;
  unsigned maskCouplingIterations = 8;
  T maskCouplingTolerance = 2.e-2;
  unsigned lastMaskCouplingIterations = 0;
  T lastMaskCouplingResidual = std::numeric_limits<T>::max();

  IndexType diffusionMinIndex{};
  IndexType diffusionMaxIndex{};
  bool diffusionBoundsSet = false;

  IndexType maskBendingMinIndex{};
  IndexType maskBendingMaxIndex{};
  bool maskBendingBoundsSet = false;

  // Populated by apply(); available for diagnostics afterwards.
  SmartPointer<OxidationDiffusion<T, D>> diffusionField;
  SmartPointer<OxidationDeformation<T, D>> deformationField;
  SmartPointer<OxidationMaskBending<T, D>> maskBendingField;
  T lastMaxVelocity_ = T(0);

  std::unordered_map<std::size_t, T> concentrationCache_;

public:
  LOCOSOxidation() = default;

  LOCOSOxidation(SmartPointer<Domain<T, D>> passedSiInterface,
                 SmartPointer<Domain<T, D>> passedAmbientInterface,
                 SmartPointer<Domain<T, D>> passedMaskInterface)
      : siInterface(passedSiInterface),
        ambientInterface(passedAmbientInterface),
        maskInterface(passedMaskInterface) {}

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<LOCOSOxidation>::New(std::forward<Args>(args)...);
  }

  void setSiInterface(SmartPointer<Domain<T, D>> si) { siInterface = si; }
  void setAmbientInterface(SmartPointer<Domain<T, D>> ambient) {
    ambientInterface = ambient;
  }
  void setMaskInterface(SmartPointer<Domain<T, D>> mask) {
    maskInterface = mask;
  }

  void setOxidationParameters(OxidationParameters<T> params) {
    oxidationParams = params;
  }
  void setDeformationParameters(OxidationDeformationParameters<T> params) {
    deformationParams = params;
  }
  void setCouplingParameters(OxidationCouplingParameters<T> params) {
    couplingParams = params;
  }
  void setMaskParameters(OxidationMaskParameters<T> params) {
    maskParams = params;
  }

  /// Set the spatial integration scheme for all three advections.
  /// Mirrors the lsAdvect::setSpatialScheme API.
  void setSpatialScheme(SpatialSchemeEnum scheme) { spatialScheme = scheme; }

  /// Set the temporal integration scheme for all three advections.
  /// Mirrors the lsAdvect::setTemporalScheme API.
  void setTemporalScheme(TemporalSchemeEnum scheme) {
    temporalScheme = scheme;
  }

  void setMaskCouplingIterations(unsigned iterations) {
    maskCouplingIterations = std::max(1u, iterations);
  }

  void setMaskCouplingTolerance(T tolerance) {
    maskCouplingTolerance = std::max(tolerance, T(0));
  }

  /// Set the Cartesian index bounding box for the diffusion and deformation
  /// solves. Required before the first call to apply().
  void setSolveBounds(const IndexType &minIndex, const IndexType &maxIndex) {
    diffusionMinIndex = minIndex;
    diffusionMaxIndex = maxIndex;
    diffusionBoundsSet = true;
  }

  /// Set the Cartesian index bounding box for the mask bending solve.
  /// The box should bracket the mask geometry with a one-cell margin on each
  /// side (one row below the mask bottom and one row above the mask top).
  /// Required before the first call to apply().
  void setMaskBendingBounds(const IndexType &minIndex,
                            const IndexType &maxIndex) {
    maskBendingMinIndex = minIndex;
    maskBendingMaxIndex = maxIndex;
    maskBendingBoundsSet = true;
  }

  /// Return the diffusion field populated by the most recent apply() call.
  SmartPointer<OxidationDiffusion<T, D>>
  getDiffusionField() const {
    return diffusionField;
  }

  /// Return the deformation field populated by the most recent apply() call.
  SmartPointer<OxidationDeformation<T, D>>
  getDeformationField() const {
    return deformationField;
  }

  /// Return the mask bending field populated by the most recent apply() call.
  SmartPointer<OxidationMaskBending<T, D>>
  getMaskBendingField() const {
    return maskBendingField;
  }

  unsigned getMaskCouplingIterations() const {
    return lastMaskCouplingIterations;
  }

  T getMaskCouplingResidual() const { return lastMaskCouplingResidual; }

  /// Maximum interface velocity (µm/hr) from the most recent CFL-limited step.
  /// Use this to compute the next step request rather than reusing actualDt.
  T getLastMaxVelocity() const { return lastMaxVelocity_; }

  /// Execute one complete LOCOS time step of duration advectionTime.
  void apply(T advectionTime) { applyImpl(advectionTime, std::nullopt); }

  /// Execute one LOCOS step, but limit the physical advection time by a CFL
  /// estimate computed from the freshly solved diffusion/deformation/mask
  /// velocity fields. Returns the actual advanced time.
  T applyCFLLimited(T requestedTime, T cflFactor) {
    return applyImpl(requestedTime,
                     std::clamp(cflFactor, T(1e-3), T(0.499)));
  }

private:
  static void logInfo(const std::string &message) {
    if (Logger::hasInfo())
      Logger::getInstance().addInfo(message).print();
  }

  T applyImpl(T requestedTime, std::optional<T> cflFactor) {
    if (siInterface == nullptr || ambientInterface == nullptr ||
        maskInterface == nullptr) {
      Logger::getInstance()
          .addError("LOCOSOxidation: one or more level-set interfaces are null.")
          .print();
      return T(0);
    }

    if (requestedTime <= T(0))
      return T(0);

    logInfo("LOCOS: starting time step, requested_dt=" +
            std::to_string(requestedTime) + " hr");

    auto solveFields = [&](T stressTimeStep) {
      auto stepDeformationParams = deformationParams;
      stepDeformationParams.stressTimeStep = stressTimeStep;

      // --- Coupled diffusion + deformation solve ---

      diffusionField = OxidationDiffusion<T, D>::New(
          siInterface, ambientInterface, oxidationParams);
      diffusionField->setConcentrationCache(concentrationCache_);
      diffusionField->setMaskInterface(maskInterface, maskInteriorSign);

      deformationField = OxidationDeformation<T, D>::New(
          siInterface, ambientInterface, diffusionField, oxidationParams,
          stepDeformationParams);
      deformationField->setMaskInterface(maskInterface, maskInteriorSign);
      // Velocity, pressure, and stress history are restored from
      // ambientInterface->getPointData() inside deformationField->apply()
      // via seedFromLevelSet() — no separate in-memory cache needed.

      // OxidationModel::apply() forwards solve bounds to the sub-solvers; if
      // bounds were not set, the sub-solvers auto-compute from the level-set
      // narrow band — no explicit setSolveBounds call is needed here.
      auto coupledModel = OxidationModel<T, D>::New(
          diffusionField, deformationField, couplingParams);
      if (diffusionBoundsSet)
        coupledModel->setSolveBounds(diffusionMinIndex, diffusionMaxIndex);
      logInfo("LOCOS: solving coupled diffusion/deformation field for dt=" +
              std::to_string(stressTimeStep) + " hr");
      coupledModel->apply();
      logInfo("LOCOS: coupled diffusion/deformation solve complete");

      // --- Mask bending solve ---

      // The bending solve domain is the nitride interior; oxide-side contact
      // faces drive the mask through traction.
      maskBendingField = OxidationMaskBending<T, D>::New(
          deformationField, maskInterface, maskParams, maskInteriorSign);
      maskBendingField->setAmbientInterface(ambientInterface, maskInteriorSign);
      if (maskBendingBoundsSet)
        maskBendingField->setSolveBounds(maskBendingMinIndex,
                                         maskBendingMaxIndex);
      logInfo("LOCOS: solving mask bending field");
      maskBendingField->apply();
      T initialRes = maskBendingField->getLastApplyVelocityChange();
      if (initialRes >= std::numeric_limits<T>::max() * T(0.99)) {
        logInfo("LOCOS: mask bending solve complete, residual=initial");
      } else {
        logInfo("LOCOS: mask bending solve complete, residual=" +
                std::to_string(initialRes));
      }

      lastMaskCouplingIterations = 1;
      lastMaskCouplingResidual = maskBendingField->getLastApplyVelocityChange();
      deformationField->setMaskVelocityField(maskBendingField);
      for (unsigned iteration = 1; iteration < maskCouplingIterations;
           ++iteration) {
        deformationField->setMaskVelocityField(maskBendingField);
        logInfo("LOCOS: coupling iteration " + std::to_string(iteration + 1) +
                " solving coupled field");
        coupledModel->apply();
        logInfo("LOCOS: coupling iteration " + std::to_string(iteration + 1) +
                " solving mask field");
        maskBendingField->apply();
        lastMaskCouplingIterations = iteration + 1;
        lastMaskCouplingResidual = maskBendingField->getLastApplyVelocityChange();
        logInfo("LOCOS: coupling iteration " + std::to_string(iteration + 1) +
                " residual=" + std::to_string(lastMaskCouplingResidual));
        if (lastMaskCouplingResidual <= maskCouplingTolerance)
          break;
      }
      if (lastMaskCouplingResidual <= maskCouplingTolerance) {
        logInfo("LOCOS: mask/oxide coupling converged in " +
                std::to_string(lastMaskCouplingIterations) +
                " iterations (residual=" +
                std::to_string(lastMaskCouplingResidual) + ")");
      } else {
        Logger::getInstance()
            .addWarning("LOCOSOxidation: mask/oxide coupling did not converge "
                        "after " + std::to_string(lastMaskCouplingIterations) +
                        " iterations (residual=" +
                        std::to_string(lastMaskCouplingResidual) +
                        ", tolerance=" +
                        std::to_string(maskCouplingTolerance) +
                        "). Consider increasing maskCouplingIterations.")
            .print();
      }

      concentrationCache_ = diffusionField->getConcentrationCache();
    };

    solveFields(requestedTime);

    T advectionTime = requestedTime;
    if (cflFactor) {
      T maxVelocity = diffusionField->getDissipationAlpha(0, -1, {});
      for (unsigned d = 0; d < D; ++d) {
        maxVelocity = std::max(maxVelocity,
                               deformationField->getDissipationAlpha(d, -1, {}));
        maxVelocity = std::max(maxVelocity,
                               maskBendingField->getDissipationAlpha(d, -1, {}));
      }
      lastMaxVelocity_ = maxVelocity;
      if (maxVelocity > std::numeric_limits<T>::epsilon()) {
        const T gridDelta = siInterface->getGrid().getGridDelta();
        advectionTime = std::min(advectionTime,
                                 (*cflFactor) * gridDelta / maxVelocity);
      }
      logInfo("LOCOS: CFL decision requested_dt=" +
              std::to_string(requestedTime) +
              " hr, actual_dt=" + std::to_string(advectionTime) +
              " hr, max_velocity=" + std::to_string(maxVelocity) + " um/hr");

      // The viscoelastic stress update depends on the step duration. If CFL
      // reduced the step, repeat the coupled solve with the actual duration
      // before advecting.
      if (advectionTime < requestedTime * (T(1) - T(1e-8)))
        solveFields(advectionTime);
    }

    // Ambient nodes inside the nitride are clamped to the mask bending velocity
    // instead of the free oxide velocity.
    auto constrainedAmbient =
        OxidationConstrainedAmbient<T, D>::New(
            deformationField, maskBendingField, maskInterface,
            maskInteriorSign);

    // Enforce oxide/mask non-penetration each step.
    // Sub-grid gaps (oxide surface below mask bottom by < Δx) close naturally
    // under the deformation velocity — isMaskContact is false in the gap zone,
    // so the Stokes velocity drives the oxide surface toward the mask without
    // the corner-snap artifact that a CUSTOM comparator would introduce.
    auto applyMaskContact = [&]() {
      BooleanOperation<T, D>(ambientInterface, maskInterface,
                             BooleanOperationEnum::RELATIVE_COMPLEMENT).apply();
    };

    // Prevent concurrent apply() calls inside lsAdvect's parallel velocity
    // queries: mark the diffusion solution valid so getScalarVelocity() never
    // re-enters apply() from multiple threads.
    diffusionField->markSolved();

    // --- Pre-advection clip (mandatory) ---
    // The ambient interface must not overlap the mask before advection begins.
    applyMaskContact();

    // Persist all solver fields after the clip so pointData sizes match the
    // post-clip HRLE domains that lsAdvect will remap.
    // lsInterior then fills each body's interior with all warm-start data.
    diffusionField->writePersistentFields();
    deformationField->writeFieldsToLevelSet();
    maskBendingField->writeFieldsToLevelSet();

    // --- Three level-set advections ---

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

    advect(ambientInterface, constrainedAmbient);
    advect(siInterface, diffusionField);
    advect(maskInterface, maskBendingField);

    // Fill each body's interior so the next substep can warm-start from
    // pointData across substep and outer-step boundaries.
    Interior<T, D>(ambientInterface).apply();
    Interior<T, D>(maskInterface).apply();

    // --- Post-advection clip (mandatory) ---
    // Removes any oxide that advected into the mask body during this step.
    applyMaskContact();

    logInfo("LOCOS: time step complete, actual_dt=" + std::to_string(advectionTime) +
            " hr");

    return advectionTime;
  }
};

} // namespace viennals
