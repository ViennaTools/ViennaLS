#pragma once

#include <lsOxidationModel.hpp>

namespace viennals {

/// Named parameter presets for common thermal oxidation processes.
///
/// Each method returns a fully populated parameter struct for a specific,
/// documented process condition. The names encode the condition deliberately —
/// a user calibrating for a different temperature or ambient should supply
/// their own struct rather than modifying these.
///
/// Sources:
///   Deal-Grove coefficients — B. E. Deal and A. S. Grove, J. Appl. Phys. 36,
///     3770 (1965); wet-oxidation values at 1000 °C from Table I.
///   Oxide mechanics — E. A. Irene, J. Electrochem. Soc. 125, 1708 (1978);
///     viscosity and elastic moduli for thermal SiO2 near 1000 °C.
///   Silicon nitride creep — S. M. Hu, J. Appl. Phys. 70, R53 (1991);
///     reference viscosity for LPCVD Si3N4 at 1000 °C.
struct OxidationPresets {
  /// Wet oxidation at 1000 °C (Deal-Grove linear-parabolic coefficients).
  /// B/2 = D_eff = 0.157 µm²/hr, B/A = k_s = 0.74 µm/hr.
  /// Stress-coupling activation volumes from Kao et al. (1987).
  static OxidationParameters wet1000CDealGrove100() {
    OxidationParameters params;
    params.diffusionCoefficient = 0.157; // µm²/hr, B/2 = D_eff
    params.reactionRate = 0.74;          // µm/hr, B/A = k_s
    params.transferCoefficient = 100.;   // large → surface conc. ≈ C*
    params.equilibriumConcentration = 1.;
    params.oxidantMoleculeDensity = 1.;
    params.expansionCoefficient = 2.27;         // SiO2/Si volume ratio
    params.temperature = 1273.15;               // 1000 °C in K
    params.reactionActivationVolume = 1.76e-35; // m³, stress activation for k_s
    params.diffusionActivationVolume = 0.;
    return params;
  }

  /// Thermal SiO2 viscoelastic mechanics at 1000 °C.
  /// Viscosity ~1×10¹⁰ Pa·hr (Irene 1978), shear modulus ~3×10¹⁰ Pa.
  /// Pass the advection time step so the viscoelastic relaxation is consistent.
  static OxidationDeformationParameters oxideMechanics1000C(double timeStep) {
    OxidationDeformationParameters params;
    params.viscosity = 1.e10;    // Pa·hr, effective oxide viscosity
    params.bulkModulus = 7.5e8;  // Pa
    params.shearModulus = 3.e10; // Pa
    params.stressTimeStep = timeStep;
    params.mechanicsIterations = 300;
    params.mechanicsTolerance = 5.e-3;
    params.pressureIterations = 2000;
    params.stokesIterations = 1000;
    params.pressureTolerance = 1.e-3;
    params.stokesTolerance = 1.e-3;
    params.tolerance = 1.e-7;
    return params;
  }

  /// LPCVD Si3N4 mask creep at 1000 °C (Hu 1991; Senez et al. 1994).
  /// referenceViscosity = 5×10¹¹ Pa·hr at 1273.15 K with Arrhenius scaling:
  /// E_a = 3.86×10⁵ J/mol ≈ 4 eV (Senez, IEEE Trans. Electron Devices, 1994).
  /// Set temperature in OxidationMaskParameters (or via psOxidation) to get
  /// the correct viscosity at a temperature other than 1000 °C.
  static OxidationMaskParameters siliconNitrideMask1000C() {
    OxidationMaskParameters params;
    params.temperature = 1273.15;
    params.referenceTemperature = 1273.15;
    params.referenceViscosity = 5.e11;     // Pa·hr at 1000 °C
    params.creepActivationEnergy = 3.86e5; // J/mol ≈ 4 eV (Senez 1994)
    params.poissonRatio = 0.27;
    params.relaxation = 0.9;
    params.tolerance = 5.e-6;
    return params;
  }
};

} // namespace viennals
