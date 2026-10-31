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
///   Stress-dependent viscosity — P. Sutardja and W. G. Oldham, IEEE Trans.
///     Electron Devices 36, 2415 (1989); shear activation volume 300 A^3 (the
///     common TCAD default).
///   Nitride creep — V. Senez et al., IEEE Trans. Electron Devices 43, 720
///     (1996), Table II (calibrated on LOCOS shapes, 800-1100 °C).
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
    params.viscosity = 1.e10;              // Pa·hr, effective oxide viscosity
    params.shearModulus = 3.e10;           // Pa
    params.shearActivationVolume = 3.e-28; // m^3 (300 A^3)
    params.stressTimeStep = timeStep;
    return params;
  }

  /// LPCVD Si3N4 mask: E = 270 GPa, nu = 0.27, and Maxwell creep with an
  /// Eyring viscosity from Senez et al. (1996), Table II: eta_N = 6e16 /
  /// 3.1e16 / 5e15 / 1.1e15 P and V_N = 120 / 130 / 170 / 350 A^3 at 800 /
  /// 900 / 1000 / 1100 °C (relaxation ~1 h at 1000 °C, critical shear
  /// 2kT/V_N 0.21 GPa).  Between the points log(eta_N) is linear in 1/T and
  /// V_N linear in T; outside, eta_N follows the end segment's trend and V_N
  /// keeps the end value (the table covers 800-1100 °C).
  static OxidationMaskParameters siliconNitrideMask(double temperatureCelsius) {
    static constexpr double table[4][3] = {
        // °C, eta_N (Pa hr = P / 10 / 3600), V_N (m^3)
        {800., 6.0e16 / 36000., 120e-30},
        {900., 3.1e16 / 36000., 130e-30},
        {1000., 5.0e15 / 36000., 170e-30},
        {1100., 1.1e15 / 36000., 350e-30}};
    const double inverse = 1. / (temperatureCelsius + 273.15);
    int k = 0; // segment [k, k + 1]
    while (k < 2 && temperatureCelsius > table[k + 1][0])
      ++k;
    const double i0 = 1. / (table[k][0] + 273.15);
    const double i1 = 1. / (table[k + 1][0] + 273.15);
    const double w = (inverse - i0) / (i1 - i0);
    const double wT = std::clamp((temperatureCelsius - table[k][0]) /
                                     (table[k + 1][0] - table[k][0]),
                                 0., 1.);
    OxidationMaskParameters params;
    params.youngModulus = 270e9; // Pa
    params.poissonRatio = 0.27;
    params.viscosity = std::exp((1. - w) * std::log(table[k][1]) +
                                w * std::log(table[k + 1][1]));
    params.shearActivationVolume =
        (1. - wT) * table[k][2] + wT * table[k + 1][2];
    return params;
  }
  static OxidationMaskParameters siliconNitrideMask1000C() {
    return siliconNitrideMask(1000.);
  }
};

} // namespace viennals
