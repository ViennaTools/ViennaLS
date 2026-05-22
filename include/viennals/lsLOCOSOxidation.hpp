#pragma once

#include <lsAdvect.hpp>
#include <lsAdvectIntegrationSchemes.hpp>
#include <lsBooleanOperation.hpp>
#include <lsOxidationModel.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

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
///   1. Diffusion+deformation coupled solve (OxidationCoupledModel)
///   2. Nitride mask bending solve (OxidationMaskBendingVelocityField)
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
      TemporalSchemeEnum::RUNGE_KUTTA_2ND_ORDER;
  static constexpr int maskInteriorSign = -1;
  unsigned maskCouplingIterations = 6;
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
  SmartPointer<OxidationDiffusionVelocityField<T, D>> diffusionField;
  SmartPointer<OxidationDeformationVelocityField<T, D>> deformationField;
  SmartPointer<OxidationMaskBendingVelocityField<T, D>> maskBendingField;

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
  SmartPointer<OxidationDiffusionVelocityField<T, D>>
  getDiffusionField() const {
    return diffusionField;
  }

  /// Return the deformation field populated by the most recent apply() call.
  SmartPointer<OxidationDeformationVelocityField<T, D>>
  getDeformationField() const {
    return deformationField;
  }

  /// Return the mask bending field populated by the most recent apply() call.
  SmartPointer<OxidationMaskBendingVelocityField<T, D>>
  getMaskBendingField() const {
    return maskBendingField;
  }

  unsigned getMaskCouplingIterations() const {
    return lastMaskCouplingIterations;
  }

  T getMaskCouplingResidual() const { return lastMaskCouplingResidual; }

  /// Execute one complete LOCOS time step of duration advectionTime.
  void apply(T advectionTime) {
    if (siInterface == nullptr || ambientInterface == nullptr ||
        maskInterface == nullptr) {
      Logger::getInstance()
          .addError("LOCOSOxidation: one or more level-set interfaces are null.")
          .print();
      return;
    }
    if (!diffusionBoundsSet) {
      Logger::getInstance()
          .addError(
              "LOCOSOxidation: diffusion solve bounds not set. "
              "Call setSolveBounds() before apply().")
          .print();
      return;
    }
    if (!maskBendingBoundsSet) {
      Logger::getInstance()
          .addError(
              "LOCOSOxidation: mask bending bounds not set. "
              "Call setMaskBendingBounds() before apply().")
          .print();
      return;
    }

    // --- Coupled diffusion + deformation solve ---

    diffusionField = OxidationDiffusionVelocityField<T, D>::New(
        siInterface, ambientInterface, oxidationParams);
    diffusionField->setMaskInterface(maskInterface, maskInteriorSign);
    diffusionField->setSolveBounds(diffusionMinIndex, diffusionMaxIndex);

    deformationField = OxidationDeformationVelocityField<T, D>::New(
        siInterface, ambientInterface, diffusionField, oxidationParams,
        deformationParams);
    deformationField->setMaskInterface(maskInterface, maskInteriorSign);
    deformationField->setSolveBounds(diffusionMinIndex, diffusionMaxIndex);

    auto coupledModel = OxidationCoupledModel<T, D>::New(
        diffusionField, deformationField, couplingParams);
    coupledModel->setSolveBounds(diffusionMinIndex, diffusionMaxIndex);
    coupledModel->apply();

    // --- Mask bending solve ---

    // The bending solve domain is the nitride interior; oxide-side contact
    // faces drive the mask through traction.
    maskBendingField = OxidationMaskBendingVelocityField<T, D>::New(
        deformationField, maskInterface, maskParams, maskInteriorSign);
    maskBendingField->setAmbientInterface(ambientInterface, maskInteriorSign);
    maskBendingField->setSolveBounds(maskBendingMinIndex, maskBendingMaxIndex);
    maskBendingField->apply();

    lastMaskCouplingIterations = 1;
    lastMaskCouplingResidual = maskBendingField->getLastApplyVelocityChange();
    deformationField->setMaskVelocityField(maskBendingField);
    for (unsigned iteration = 1; iteration < maskCouplingIterations; ++iteration) {
      deformationField->setMaskVelocityField(maskBendingField);
      coupledModel->apply();
      maskBendingField->apply();
      lastMaskCouplingIterations = iteration + 1;
      lastMaskCouplingResidual = maskBendingField->getLastApplyVelocityChange();
      if (lastMaskCouplingResidual <= maskCouplingTolerance)
        break;
    }

    // Ambient nodes inside the nitride are clamped to the mask bending velocity
    // instead of the free oxide velocity.
    auto constrainedAmbient =
        OxidationConstrainedAmbientVelocityField<T, D>::New(
            deformationField, maskBendingField, maskInterface,
            maskInteriorSign);

    // --- Pre-advection clip (mandatory) ---
    // The ambient interface must not overlap the mask before advection begins.
    BooleanOperation<T, D>(ambientInterface, maskInterface,
                           BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();

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

    // --- Post-advection clip (mandatory) ---
    // Corrects any penetration of the ambient interface into the mask volume
    // that accumulated during the time step.
    BooleanOperation<T, D>(ambientInterface, maskInterface,
                           BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }
};

} // namespace viennals
