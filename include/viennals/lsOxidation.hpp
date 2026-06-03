#pragma once

#include <lsAdvect.hpp>
#include <lsAdvectIntegrationSchemes.hpp>
#include <lsBooleanOperation.hpp>
#include <lsInterior.hpp>
#include <lsOxidationModel.hpp>

#include <algorithm>
#include <cmath>
#include <iostream>

#include <vcTimer.hpp>
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

/// Unified oxidation time-step orchestrator.
///
/// Couples diffusion, viscoelastic deformation, and (optionally) mask bending
/// into a single CFL-limited step. When a mask interface is provided the solver
/// operates in LOCOS mode; without one it performs standard planar or step
/// oxidation.
///
/// Standard mode (no mask):
///   1. Diffusion+deformation coupled solve (OxidationModel)
///   2. Ambient interface advection (deformation velocity)
///   3. Si/SiO2 reaction interface advection (diffusion velocity)
///   4. Interior fill for warm-start persistence
///
/// LOCOS mode (mask provided):
///   1. Diffusion+deformation coupled solve (OxidationModel)
///   2. Nitride mask bending solve (OxidationMaskBending) + coupling loop
///   3. Pre-advection boolean clip (ambientInterface \ maskInterface)
///   4. Ambient interface advection (constrained: mask-coupled under nitride)
///   5. Si/SiO2 reaction interface advection (diffusion velocity)
///   6. Mask advection (mask bending velocity)
///   7. Interior fills + post-advection boolean clip
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
///   ox->setMaskBendingBounds(maskMinIdx, maskMaxIdx);
///   ox->applyCFLLimited(0.35, 0.499);
///   ox->getDiffusionField()->getConcentration(pt);  // diagnostics
/// @endcode
template <class T, int D> class Oxidation {
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

  /// Set the spatial integration scheme for all advections.
  void setSpatialScheme(SpatialSchemeEnum scheme) { spatialScheme = scheme; }

  /// Set the temporal integration scheme for all advections.
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
  /// solves. If not set, bounds are auto-computed from the level-set narrow band.
  void setSolveBounds(const IndexType &minIndex, const IndexType &maxIndex) {
    diffusionMinIndex = minIndex;
    diffusionMaxIndex = maxIndex;
    diffusionBoundsSet = true;
  }

  /// Set the Cartesian index bounding box for the mask bending solve (LOCOS).
  void setMaskBendingBounds(const IndexType &minIndex,
                            const IndexType &maxIndex) {
    maskBendingMinIndex = minIndex;
    maskBendingMaxIndex = maxIndex;
    maskBendingBoundsSet = true;
  }

  /// Return the diffusion field populated by the most recent apply() call.
  SmartPointer<OxidationDiffusion<T, D>> getDiffusionField() const {
    return diffusionField;
  }

  /// Return the deformation field populated by the most recent apply() call.
  SmartPointer<OxidationDeformation<T, D>> getDeformationField() const {
    return deformationField;
  }

  /// Return the mask bending field (null when no mask is set).
  SmartPointer<OxidationMaskBending<T, D>> getMaskBendingField() const {
    return maskBendingField;
  }

  unsigned getMaskCouplingIterations() const {
    return lastMaskCouplingIterations;
  }

  T getMaskCouplingResidual() const { return lastMaskCouplingResidual; }

  /// Maximum interface velocity (µm/hr) from the most recent CFL-limited step.
  T getLastMaxVelocity() const { return lastMaxVelocity_; }

  /// Execute one oxidation time step of duration advectionTime.
  void apply(T advectionTime) { applyImpl(advectionTime, std::nullopt); }

  /// Execute one CFL-limited oxidation step; returns the actual time advanced.
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
    const std::string prefix = hasMask ? "LOCOS" : "Oxidation";

    logInfo(prefix + ": starting time step, requested_dt=" +
            std::to_string(requestedTime) + " hr");

    auto solveFields = [&](T stressTimeStep) {
      auto stepDeformationParams = deformationParams;
      stepDeformationParams.stressTimeStep = stressTimeStep;

      // --- Coupled diffusion + deformation solve ---

      diffusionField = OxidationDiffusion<T, D>::New(
          siInterface, ambientInterface, oxidationParams);
      diffusionField->setConcentrationCache(concentrationCache_);
      if (hasMask)
        diffusionField->setMaskInterface(maskInterface, maskInteriorSign);

      deformationField = OxidationDeformation<T, D>::New(
          siInterface, ambientInterface, diffusionField, oxidationParams,
          stepDeformationParams);
      if (hasMask)
        deformationField->setMaskInterface(maskInterface, maskInteriorSign);

      auto coupledModel = OxidationModel<T, D>::New(
          diffusionField, deformationField, couplingParams);
      if (diffusionBoundsSet)
        coupledModel->setSolveBounds(diffusionMinIndex, diffusionMaxIndex);
      logInfo(prefix + ": solving coupled diffusion/deformation field for dt=" +
              std::to_string(stressTimeStep) + " hr");
      Timer<> tCoupled;
      tCoupled.start();
      coupledModel->apply();
      tCoupled.finish();
      if (Logger::hasTiming())
        Logger::getInstance().addTiming("  coupled(iter=1)", tCoupled).print();
      logInfo(prefix + ": coupled diffusion/deformation solve complete");

      if (hasMask) {
        // --- Mask bending solve ---
        maskBendingField = OxidationMaskBending<T, D>::New(
            deformationField, maskInterface, maskParams, maskInteriorSign);
        maskBendingField->setAmbientInterface(ambientInterface, maskInteriorSign);
        if (maskBendingBoundsSet)
          maskBendingField->setSolveBounds(maskBendingMinIndex,
                                           maskBendingMaxIndex);
        logInfo(prefix + ": solving mask bending field");
        Timer<> tMask;
        tMask.start();
        maskBendingField->apply();
        tMask.finish();
        if (Logger::hasTiming())
          Logger::getInstance().addTiming("  maskBending(iter=1)", tMask).print();
        T initialRes = maskBendingField->getLastApplyVelocityChange();
        logInfo(prefix + ": mask bending solve complete, residual=" +
                (initialRes >= std::numeric_limits<T>::max() * T(0.99)
                     ? std::string("initial")
                     : std::to_string(initialRes)));

        lastMaskCouplingIterations = 1;
        lastMaskCouplingResidual = maskBendingField->getLastApplyVelocityChange();
        deformationField->setMaskVelocityField(maskBendingField);
        for (unsigned iteration = 1; iteration < maskCouplingIterations;
             ++iteration) {
          deformationField->setMaskVelocityField(maskBendingField);
          logInfo(prefix + ": coupling iteration " +
                  std::to_string(iteration + 1) + " solving coupled field");
          Timer<> tIterCoupled, tIterMask;
          tIterCoupled.start();
          coupledModel->apply();
          tIterCoupled.finish();
          logInfo(prefix + ": coupling iteration " +
                  std::to_string(iteration + 1) + " solving mask field");
          tIterMask.start();
          maskBendingField->apply();
          tIterMask.finish();
          if (Logger::hasTiming())
            Logger::getInstance()
                .addTiming("  coupled(iter=" + std::to_string(iteration + 1) + ")",
                            tIterCoupled)
                .addTiming("  maskBending(iter=" + std::to_string(iteration + 1) + ")",
                            tIterMask)
                .print();
          lastMaskCouplingIterations = iteration + 1;
          lastMaskCouplingResidual = maskBendingField->getLastApplyVelocityChange();
          logInfo(prefix + ": coupling iteration " +
                  std::to_string(iteration + 1) + " residual=" +
                  std::to_string(lastMaskCouplingResidual));
          if (lastMaskCouplingResidual <= maskCouplingTolerance)
            break;
        }
        if (lastMaskCouplingResidual <= maskCouplingTolerance) {
          logInfo(prefix + ": mask/oxide coupling converged in " +
                  std::to_string(lastMaskCouplingIterations) +
                  " iterations (residual=" +
                  std::to_string(lastMaskCouplingResidual) + ")");
        } else {
          Logger::getInstance()
              .addWarning(prefix + ": mask/oxide coupling did not converge "
                          "after " + std::to_string(lastMaskCouplingIterations) +
                          " iterations (residual=" +
                          std::to_string(lastMaskCouplingResidual) +
                          ", tolerance=" + std::to_string(maskCouplingTolerance) +
                          "). Consider increasing maskCouplingIterations.")
              .print();
        }
      } else {
        maskBendingField = nullptr;
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
        if (hasMask)
          maxVelocity = std::max(maxVelocity,
                                 maskBendingField->getDissipationAlpha(d, -1, {}));
      }
      lastMaxVelocity_ = maxVelocity;
      if (maxVelocity > std::numeric_limits<T>::epsilon()) {
        const T gridDelta = siInterface->getGrid().getGridDelta();
        advectionTime = std::min(advectionTime,
                                 (*cflFactor) * gridDelta / maxVelocity);
      }
      logInfo(prefix + ": CFL decision requested_dt=" +
              std::to_string(requestedTime) +
              " hr, actual_dt=" + std::to_string(advectionTime) +
              " hr, max_velocity=" + std::to_string(maxVelocity) + " um/hr");

      if (advectionTime < requestedTime * (T(1) - T(1e-8)))
        solveFields(advectionTime);
    }

    // When a mask is present, the ambient interface is constrained to follow
    // the mask in the contact region. Without a mask the deformation field
    // drives the ambient surface directly.
    SmartPointer<VelocityField<T>> ambientVelocity;
    if (hasMask) {
      ambientVelocity = OxidationConstrainedAmbient<T, D>::New(
          deformationField, maskBendingField, maskInterface, maskInteriorSign);
    } else {
      ambientVelocity = deformationField;
    }

    diffusionField->markSolved();

    // Pre-advection clip: keep oxide outside the mask body (LOCOS only).
    if (hasMask)
      BooleanOperation<T, D>(ambientInterface, maskInterface,
                             BooleanOperationEnum::RELATIVE_COMPLEMENT).apply();

    diffusionField->writePersistentFields();
    deformationField->writeFieldsToLevelSet();
    if (hasMask)
      maskBendingField->writeFieldsToLevelSet();

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
    if (hasMask)
      advect(maskInterface, maskBendingField);
    tAdvect.finish();
    if (Logger::hasTiming())
      Logger::getInstance()
          .addTiming(std::string("  advection(") + (hasMask ? "3" : "2") + " surfaces)",
                     tAdvect)
          .print();

    Interior<T, D>(ambientInterface).apply();
    if (hasMask)
      Interior<T, D>(maskInterface).apply();

    // Post-advection clip: remove oxide that grew into the mask (LOCOS only).
    if (hasMask)
      BooleanOperation<T, D>(ambientInterface, maskInterface,
                             BooleanOperationEnum::RELATIVE_COMPLEMENT).apply();

    logInfo(prefix + ": time step complete, actual_dt=" +
            std::to_string(advectionTime) + " hr");

    tStep.finish();
    if (Logger::hasTiming())
      Logger::getInstance().addTiming("── step total", tStep).print();

    return advectionTime;
  }
};

} // namespace viennals
