#pragma once

#include <lsAdvect.hpp>
#include <lsAdvectIntegrationSchemes.hpp>
#include <lsBooleanOperation.hpp>
#include <lsInterior.hpp>
#include <lsOxidationModel.hpp>

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

  OxidationParameters oxidationParams;
  OxidationDeformationParameters deformationParams;
  OxidationCouplingParameters couplingParams;
  OxidationMaskParameters maskParams;

  SpatialSchemeEnum spatialScheme = SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER;
  TemporalSchemeEnum temporalScheme = TemporalSchemeEnum::FORWARD_EULER;
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

  GpuMode gpuMode_ = GpuMode::Cpu;
  GpuPreconditioner gpuPreconditioner_ = GpuPreconditioner::Jacobi;
  std::unordered_map<std::size_t, T> concentrationCache_;

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
  void setGpuPreconditioner(GpuPreconditioner preconditioner) {
    gpuPreconditioner_ = preconditioner;
  }

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

  /// Set the temporal integration scheme for all advections.
  void setTemporalScheme(TemporalSchemeEnum scheme) { temporalScheme = scheme; }

  void setMaskCouplingIterations(unsigned iterations) {
    maskCouplingIterations = std::max(1u, iterations);
  }

  void setMaskCouplingTolerance(T tolerance) {
    maskCouplingTolerance = std::max(tolerance, T(0));
  }

  /// Set the Cartesian index bounding box for the diffusion and deformation
  /// solves. If not set, bounds are auto-computed from the level-set narrow
  /// band.
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
    const std::string prefix = hasMask ? "LOCOS" : "Oxidation";

    VIENNACORE_LOG_INFO(prefix + ": starting time step, requested_dt=" +
                        std::to_string(requestedTime) + " hr");

    const auto baseConcentrationCache = concentrationCache_;
    std::string lastFieldFailureReason;
    bool lastFailureWasMaskFixedPoint = false;
    T adaptiveMaskRelaxationScale = T(1);

    auto solveFields = [&](T stressTimeStep, bool logCouplingResult) -> bool {
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
              std::to_string(deformationField->getResidual()) + ", pressure=" +
              std::to_string(deformationField->getLastPressureResidual()) +
              ", stokes=" +
              std::to_string(deformationField->getLastStokesResidual()) + ")");
        }
        return true;
      };

      auto validateMaskSolve = [&]() {
        if (maskBendingField == nullptr)
          return true;
        const T maskResidual = maskBendingField->getResidual();
        if (!std::isfinite(maskResidual) ||
            maskResidual > maskParams.tolerance) {
          return rejectSolve("mask traction solve failed (residual=" +
                             std::to_string(maskResidual) + ", tolerance=" +
                             std::to_string(maskParams.tolerance) + ")");
        }
        const T couplingResidual =
            maskBendingField->getLastApplyVelocityChange();
        if (!std::isfinite(couplingResidual))
          return rejectSolve(
              "mask velocity coupling produced non-finite values");
        return true;
      };

      lastFieldFailureReason.clear();
      lastFailureWasMaskFixedPoint = false;
      auto stepDeformationParams = deformationParams;
      stepDeformationParams.stressTimeStep = stressTimeStep;

      // --- Coupled diffusion + deformation solve ---

      diffusionField = OxidationDiffusion<T, D>::New(
          siInterface, ambientInterface, oxidationParams);
      diffusionField->setConcentrationCache(baseConcentrationCache);
      diffusionField->setGpuMode(gpuMode_);
      diffusionField->setGpuPreconditioner(gpuPreconditioner_);
      if (hasMask)
        diffusionField->setMaskInterface(maskInterface, maskInteriorSign);

      deformationField = OxidationDeformation<T, D>::New(
          siInterface, ambientInterface, diffusionField, oxidationParams,
          stepDeformationParams);
      deformationField->setGpuMode(gpuMode_);
      deformationField->setGpuPreconditioner(gpuPreconditioner_);
      if (hasMask)
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
      if (!validateCoupledModel(coupledModel))
        return false;

      if (hasMask) {
        // --- Mask bending solve ---
        auto stepMaskParams = maskParams;
        stepMaskParams.stressTimeStep = stressTimeStep;
        stepMaskParams.relaxation = std::clamp(
            maskParams.relaxation * adaptiveMaskRelaxationScale, T(0.01), T(1));
        maskBendingField = OxidationMaskBending<T, D>::New(
            deformationField, maskInterface, stepMaskParams, maskInteriorSign);
        maskBendingField->setAmbientInterface(ambientInterface,
                                              maskInteriorSign);
        if (maskBendingBoundsSet)
          maskBendingField->setSolveBounds(maskBendingMinIndex,
                                           maskBendingMaxIndex);
        VIENNACORE_LOG_DEBUG(prefix + ": solving mask bending field");
        Timer<> tMask;
        tMask.start();
        try {
          maskBendingField->apply();
        } catch (const std::exception &e) {
          tMask.finish();
          return rejectSolve("mask bending solve error: " +
                             std::string(e.what()));
        }
        tMask.finish();
        if (Logger::hasTiming())
          Logger::getInstance()
              .addTiming("  maskBending(iter=1)", tMask)
              .print();
        if (!validateMaskSolve())
          return false;
        T initialRes = maskBendingField->getLastApplyVelocityChange();
        VIENNACORE_LOG_DEBUG(
            prefix + ": mask bending solve complete, residual=" +
            (initialRes >= std::numeric_limits<T>::max() * T(0.99)
                 ? std::string("initial")
                 : std::to_string(initialRes)));

        lastMaskCouplingIterations = 1;
        lastMaskCouplingResidual =
            maskBendingField->getLastApplyVelocityChange();
        deformationField->setMaskVelocityField(maskBendingField);
        for (unsigned iteration = 1; iteration < maskCouplingIterations;
             ++iteration) {
          deformationField->setMaskVelocityField(maskBendingField);
          VIENNACORE_LOG_DEBUG(prefix + ": coupling iteration " +
                               std::to_string(iteration + 1) +
                               " solving coupled field");
          Timer<> tIterCoupled, tIterMask;
          tIterCoupled.start();
          coupledModel->apply();
          tIterCoupled.finish();
          if (!validateCoupledModel(coupledModel))
            return false;
          VIENNACORE_LOG_DEBUG(prefix + ": coupling iteration " +
                               std::to_string(iteration + 1) +
                               " solving mask field");
          tIterMask.start();
          try {
            maskBendingField->apply();
          } catch (const std::exception &e) {
            tIterMask.finish();
            return rejectSolve("mask bending solve error at iteration " +
                               std::to_string(iteration + 1) + ": " + e.what());
          }
          tIterMask.finish();
          if (!validateMaskSolve())
            return false;
          if (Logger::hasTiming())
            Logger::getInstance()
                .addTiming("  coupled(iter=" + std::to_string(iteration + 1) +
                               ")",
                           tIterCoupled)
                .addTiming(
                    "  maskBending(iter=" + std::to_string(iteration + 1) + ")",
                    tIterMask)
                .print();
          lastMaskCouplingIterations = iteration + 1;
          lastMaskCouplingResidual =
              maskBendingField->getLastApplyVelocityChange();
          VIENNACORE_LOG_DEBUG(
              prefix + ": coupling iteration " + std::to_string(iteration + 1) +
              " residual=" + std::to_string(lastMaskCouplingResidual));
          if (lastMaskCouplingResidual <= maskCouplingTolerance)
            break;
        }
        const T maskAbsoluteDisplacement =
            maskBendingField->getLastApplyAbsoluteVelocityChange() *
            stressTimeStep;
        // Max physical displacement per step from the mask velocity field.
        // Independent of coupling oscillation amplitude — the oscillation check
        // (maskAbsoluteDisplacement) measures how much the velocity CHANGES
        // between iterations; this measures how far the surface actually moves.
        // When the active-set oscillates at a fixed amplitude, reducing dt
        // does not reduce maskAbsoluteDisplacement, but it does reduce this.
        T maxMaskVelocity = T(0);
        for (unsigned d = 0; d < D; ++d)
          maxMaskVelocity =
              std::max(maxMaskVelocity, maskBendingField->getDissipationAlpha(
                                            static_cast<int>(d), -1, {}));
        const T maskMaxDisplacement = maxMaskVelocity * stressTimeStep;
        const T maskDisplacementTolerance =
            maskCouplingTolerance * siInterface->getGrid().getGridDelta();
        const bool maskCouplingConverged =
            lastMaskCouplingResidual <= maskCouplingTolerance ||
            (std::isfinite(maskAbsoluteDisplacement) &&
             maskAbsoluteDisplacement <= maskDisplacementTolerance) ||
            (std::isfinite(maskMaxDisplacement) &&
             maskMaxDisplacement <= maskDisplacementTolerance);
        if (maskCouplingConverged) {
          VIENNACORE_LOG_INFO(
              prefix + ": mask/oxide coupling converged in " +
              std::to_string(lastMaskCouplingIterations) +
              " iterations (residual=" +
              std::to_string(lastMaskCouplingResidual) +
              ", displacement=" + std::to_string(maskAbsoluteDisplacement) +
              " um, maxDisplacement=" + std::to_string(maskMaxDisplacement) +
              " um)");
        } else if (logCouplingResult) {
          VIENNACORE_LOG_WARNING(
              prefix +
              ": mask/oxide coupling did not converge "
              "after " +
              std::to_string(lastMaskCouplingIterations) +
              " iterations (residual=" +
              std::to_string(lastMaskCouplingResidual) +
              ", displacement=" + std::to_string(maskAbsoluteDisplacement) +
              " um" + ", tolerance=" + std::to_string(maskCouplingTolerance) +
              "). Consider increasing maskCouplingIterations.");
        }
        if (!maskCouplingConverged) {
          lastFailureWasMaskFixedPoint = true;
          return rejectSolve(
              "mask/oxide coupling failed (residual=" +
              std::to_string(lastMaskCouplingResidual) +
              ", displacement=" + std::to_string(maskAbsoluteDisplacement) +
              " um, maxDisplacement=" + std::to_string(maskMaxDisplacement) +
              " um" + ", tolerance=" + std::to_string(maskCouplingTolerance) +
              ", relaxation=" + std::to_string(stepMaskParams.relaxation) +
              ")");
        }
        return true;
      } else {
        maskBendingField = nullptr;
      }

      return true;
    };

    auto makeAmbientVelocity = [&]() -> SmartPointer<VelocityField<T>> {
      if (hasMask) {
        return OxidationConstrainedAmbient<T, D>::New(
            deformationField, maskBendingField, maskInterface, ambientInterface,
            maskInteriorSign);
      }
      return deformationField;
    };

    T advectionTime = requestedTime;
    SmartPointer<VelocityField<T>> ambientVelocity;

    auto computeMaxVelocity =
        [&](const SmartPointer<VelocityField<T>> &passedAmbientVelocity) {
          T maxVelocity = diffusionField->getDissipationAlpha(0, -1, {});
          for (unsigned d = 0; d < D; ++d) {
            maxVelocity =
                std::max(maxVelocity,
                         passedAmbientVelocity->getDissipationAlpha(d, -1, {}));
            if (hasMask)
              maxVelocity =
                  std::max(maxVelocity,
                           maskBendingField->getDissipationAlpha(d, -1, {}));
          }
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
        const bool predictorConverged = solveFields(trialTime, false);
        if (!predictorConverged) {
          const bool dampMask = lastFailureWasMaskFixedPoint &&
                                adaptiveMaskRelaxationScale > T(0.051);
          if (dampMask)
            adaptiveMaskRelaxationScale =
                std::max(T(0.05), adaptiveMaskRelaxationScale * T(0.5));
          const T nextTrial = dampMask ? trialTime : trialTime * T(0.5);
          VIENNACORE_LOG_INFO(
              prefix +
              ": rejecting non-converged coupled predictor "
              "(" +
              (lastFieldFailureReason.empty()
                   ? "mask residual=" + std::to_string(lastMaskCouplingResidual)
                   : lastFieldFailureReason) +
              (dampMask ? ", retrying with mask relaxation scale=" +
                              std::to_string(adaptiveMaskRelaxationScale)
                        : std::string()) +
              "), retrying with requested_dt=" + std::to_string(nextTrial) +
              " hr");
          trialTime = nextTrial;
          if (trialTime < minTrialTime)
            break;
          continue;
        }

        ambientVelocity = makeAmbientVelocity();
        T maxVelocity = computeMaxVelocity(ambientVelocity);
        if (!std::isfinite(maxVelocity))
          throw std::runtime_error(prefix +
                                   ": non-finite CFL velocity estimate.");

        advectionTime = cflLimitedTime(trialTime, maxVelocity);
        VIENNACORE_LOG_INFO(
            prefix +
            ": CFL decision requested_dt=" + std::to_string(trialTime) +
            " hr, actual_dt=" + std::to_string(advectionTime) +
            " hr, max_velocity=" + std::to_string(maxVelocity) + " um/hr");

        if (advectionTime < trialTime * (T(1) - T(1e-8))) {
          const bool finalConverged = solveFields(advectionTime, false);
          if (!finalConverged) {
            const bool dampMask = lastFailureWasMaskFixedPoint &&
                                  adaptiveMaskRelaxationScale > T(0.051);
            if (dampMask)
              adaptiveMaskRelaxationScale =
                  std::max(T(0.05), adaptiveMaskRelaxationScale * T(0.5));
            const T nextTrial =
                dampMask ? advectionTime : advectionTime * T(0.5);
            VIENNACORE_LOG_INFO(
                prefix +
                ": rejecting non-converged CFL re-solve "
                "(" +
                (lastFieldFailureReason.empty()
                     ? "mask residual=" +
                           std::to_string(lastMaskCouplingResidual)
                     : lastFieldFailureReason) +
                (dampMask ? ", retrying with mask relaxation scale=" +
                                std::to_string(adaptiveMaskRelaxationScale)
                          : std::string()) +
                "), retrying with requested_dt=" + std::to_string(nextTrial) +
                " hr");
            trialTime = nextTrial;
            if (trialTime < minTrialTime)
              break;
            continue;
          }
          ambientVelocity = makeAmbientVelocity();
          maxVelocity = computeMaxVelocity(ambientVelocity);
          if (!std::isfinite(maxVelocity))
            throw std::runtime_error(prefix +
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
      if (!accepted) {
        // Last resort: if the oxide solve (diffusion + deformation) is still
        // finite — only the mask coupling diverged — freeze the mask for one
        // minimal step so the geometry can evolve past the singular contact
        // configuration.  The mask doesn't move; the oxide advances at
        // minTrialTime with no mask feedback this step.
        const bool oxideFinite =
            diffusionField && deformationField &&
            diffusionField->hasFiniteConcentrationField() &&
            deformationField->hasFiniteSolution();
        if (hasMask && oxideFinite) {
          Logger::getInstance()
              .addWarning(prefix +
                          ": all CFL attempts exhausted; freezing mask for "
                          "one step at dt=" +
                          std::to_string(minTrialTime) +
                          " hr (last failure: " + lastFieldFailureReason +
                          "). Consider increasing maskReferenceViscosity or "
                          "maskCouplingIterations.")
              .print();
          maskBendingField = nullptr; // skip mask advection this step
          ambientVelocity = makeAmbientVelocity();
          advectionTime = minTrialTime;
          lastMaxVelocity_ = computeMaxVelocity(ambientVelocity);
        } else {
          VIENNACORE_LOG_ERROR(prefix +
                               ": unable to find a converged CFL-limited step" +
                               (lastFieldFailureReason.empty()
                                    ? std::string(".")
                                    : std::string(" (last failure: ") +
                                          lastFieldFailureReason + ")."));
        }
      }
    } else {
      const bool fieldsConverged = solveFields(requestedTime, true);
      if (!fieldsConverged)
        VIENNACORE_LOG_ERROR(
            prefix + ": coupled solve failed" +
            (lastFieldFailureReason.empty()
                 ? std::string(".")
                 : std::string(" (") + lastFieldFailureReason + ")."));
      ambientVelocity = makeAmbientVelocity();
    }

    concentrationCache_ = diffusionField->getConcentrationCache();

    diffusionField->markSolved();

    // Pre-advection clip: keep oxide outside the mask body (LOCOS only).
    if (hasMask)
      BooleanOperation<T, D>(ambientInterface, maskInterface,
                             BooleanOperationEnum::RELATIVE_COMPLEMENT)
          .apply();

    diffusionField->writePersistentFields();
    deformationField->writeFieldsToLevelSet();
    if (hasMask && maskBendingField)
      maskBendingField->writeFieldsToLevelSet();

    if (hasMask && maskBendingField)
      maskBendingField->finalizeElasticAdvectionVelocity();

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
    if (hasMask && maskBendingField)
      advect(maskInterface, maskBendingField);
    tAdvect.finish();
    VIENNACORE_LOG_TIMING(std::string("  advection(") + (hasMask ? "3" : "2") +
                              " surfaces)",
                          tAdvect);

    // Post-advection clip: remove oxide that grew into the mask (LOCOS only).
    // Mask gets Interior fill first so the BooleanOp has accurate φ_mask values
    // at points inside the mask region.
    if (hasMask) {
      Interior<T, D>(maskInterface).apply();
      BooleanOperation<T, D>(ambientInterface, maskInterface,
                             BooleanOperationEnum::RELATIVE_COMPLEMENT)
          .apply();
      // Re-write mask velocity with the Interior-filled HRLE for the same
      // reason as the oxide re-write below — lsAdvect left only the narrow
      // band, so pointData is smaller than the post-Interior HRLE.
      if (maskBendingField)
        maskBendingField->writeFieldsToLevelSet();
    }
    {
      Interior<T, D> fill(ambientInterface);
      fill.setGuide(siInterface); // stop fill at Si surface
      fill.apply();
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
