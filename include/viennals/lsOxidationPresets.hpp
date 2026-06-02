#pragma once

#include <lsOxidationModel.hpp>

namespace viennals {

/// Named parameter presets for common thermal oxidation processes.
///
/// Each method returns a fully populated parameter struct for a specific,
/// documented process condition. The names encode the condition deliberately —
/// a user calibrating for a different temperature or ambient should supply their
/// own struct rather than modifying these.
///
/// Sources:
///   Deal-Grove coefficients — B. E. Deal and A. S. Grove, J. Appl. Phys. 36,
///     3770 (1965); wet-oxidation values at 1000 °C from Table I.
///   Oxide mechanics — E. A. Irene, J. Electrochem. Soc. 125, 1708 (1978);
///     viscosity and elastic moduli for thermal SiO2 near 1000 °C.
///   Silicon nitride creep — S. M. Hu, J. Appl. Phys. 70, R53 (1991);
///     reference viscosity for LPCVD Si3N4 at 1000 °C.
template <class T> struct OxidationPresets {
  /// Wet oxidation at 1000 °C (Deal-Grove linear-parabolic coefficients).
  /// B/2 = D_eff = 0.157 µm²/hr, B/A = k_s = 0.74 µm/hr.
  /// Stress-coupling activation volumes from Kao et al. (1987).
  static OxidationParameters<T> wet1000CDealGrove100() {
    OxidationParameters<T> params;
    params.diffusionCoefficient = T(0.157);      // µm²/hr, B/2 = D_eff
    params.reactionRate = T(0.74);               // µm/hr, B/A = k_s
    params.transferCoefficient = T(100.);        // large → surface conc. ≈ C*
    params.equilibriumConcentration = T(1.);
    params.oxidantMoleculeDensity = T(1.);
    params.expansionCoefficient = T(2.27);       // SiO2/Si volume ratio
    params.temperature = T(1273.15);             // 1000 °C in K
    params.reactionActivationVolume = T(1.76e-35); // m³, stress activation for k_s
    params.diffusionActivationVolume = T(0.);
    return params;
  }

  /// Thermal SiO2 viscoelastic mechanics at 1000 °C.
  /// Viscosity ~1×10¹⁰ Pa·hr (Irene 1978), shear modulus ~3×10¹⁰ Pa.
  /// Pass the advection time step so the viscoelastic relaxation is consistent.
  static OxidationDeformationParameters<T> oxideMechanics1000C(T timeStep) {
    OxidationDeformationParameters<T> params;
    params.viscosity = T(1.e10);                 // Pa·hr, effective oxide viscosity
    params.bulkModulus = T(7.5e8);               // Pa
    params.shearModulus = T(3.e10);              // Pa
    params.stressTimeStep = timeStep;
    params.mechanicsIterations = 2;
    params.mechanicsTolerance = T(1.e-7);
    params.pressureIterations = 500;
    params.stokesIterations = 100;
    params.pressureTolerance = T(1.e-6);
    params.stokesTolerance = T(1.e-7);
    params.tolerance = T(1.e-7);
    return params;
  }

  /// LPCVD Si3N4 mask creep at 1000 °C (Hu 1991).
  /// referenceViscosity = 5×10¹¹ Pa·hr; creepActivationEnergy = 0 disables
  /// Arrhenius scaling — set it (J/mol) to enable temperature dependence.
  static OxidationMaskParameters<T> siliconNitrideMask1000C() {
    OxidationMaskParameters<T> params;
    params.temperature = T(1273.15);
    params.referenceTemperature = T(1273.15);
    params.referenceViscosity = T(5.e11);        // Pa·hr
    params.creepActivationEnergy = T(0.);        // set to enable Arrhenius scaling
    params.poissonRatio = T(0.27);
    params.relaxation = T(0.9);
    params.tolerance = T(5.e-6);
    return params;
  }
};

} // namespace viennals
