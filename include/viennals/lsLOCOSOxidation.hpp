#pragma once

#include <lsAdvect.hpp>
#include <lsAdvectIntegrationSchemes.hpp>
#include <lsBooleanOperation.hpp>
#include <lsOxidationModel.hpp>

namespace viennals {

using namespace viennacore;

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
/// Level-set sign conventions used internally:
///   - diffusion / deformation / constrained-ambient: maskSign = -1
///     (the solve domain is the oxide, inside the nitride where phi_mask < 0)
///   - mask bending: maskSign = +1
///     (the solve domain is the exterior of the nitride, phi_mask > 0)
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
    diffusionField->setMaskInterface(maskInterface, -1);
    diffusionField->setSolveBounds(diffusionMinIndex, diffusionMaxIndex);

    deformationField = OxidationDeformationVelocityField<T, D>::New(
        siInterface, ambientInterface, diffusionField, oxidationParams,
        deformationParams);
    deformationField->setMaskInterface(maskInterface, -1);
    deformationField->setSolveBounds(diffusionMinIndex, diffusionMaxIndex);

    auto coupledModel = OxidationCoupledModel<T, D>::New(
        diffusionField, deformationField, couplingParams);
    coupledModel->setSolveBounds(diffusionMinIndex, diffusionMaxIndex);
    coupledModel->apply();

    // --- Mask bending solve ---

    // maskSign = +1: the bending solve domain is the exterior of the nitride
    // (phi_mask > 0); contact nodes at the oxide/nitride interface drive it.
    maskBendingField = OxidationMaskBendingVelocityField<T, D>::New(
        deformationField, maskInterface, maskParams, 1);
    maskBendingField->setSolveBounds(maskBendingMinIndex, maskBendingMaxIndex);
    maskBendingField->apply();

    // maskSign = -1: ambient nodes inside the nitride (phi_mask < 0) are
    // clamped to the mask bending velocity instead of the free oxide velocity.
    auto constrainedAmbient =
        OxidationConstrainedAmbientVelocityField<T, D>::New(
            deformationField, maskBendingField, maskInterface, -1);

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
