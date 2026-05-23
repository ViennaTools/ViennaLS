#pragma once

#include <lsOxidationModel.hpp>

namespace viennals {

/// Convenience presets for common oxidation examples.
///
/// These values are intentionally ordinary parameter objects, not hidden global
/// state. They make examples reproducible while keeping calibration data in one
/// place that users can replace with their own process/material database.
template <class T> struct OxidationProcessPresets {
  static OxidationParameters<T> wet1000CDealGrove100() {
    OxidationParameters<T> params;
    params.diffusionCoefficient = T(0.157);      // um^2/hr, B ~= 2D
    params.reactionRate = T(0.74);               // um/hr, B/A
    params.transferCoefficient = T(100.);        // gas transfer, large -> C ~= C*
    params.equilibriumConcentration = T(1.);
    params.oxidantMoleculeDensity = T(1.);
    params.expansionCoefficient = T(2.27);
    params.temperature = T(1273.15);             // 1000 C
    params.reactionActivationVolume = T(1.76e-35);
    params.diffusionActivationVolume = T(0.);
    return params;
  }

  static OxidationDeformationParameters<T> oxideMechanics1000C(T timeStep) {
    OxidationDeformationParameters<T> params;
    params.viscosity = T(1.e10);                 // Pa hr, effective oxide viscosity
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

  static OxidationMaskParameters<T> siliconNitrideMask1000C() {
    OxidationMaskParameters<T> params;
    params.temperature = T(1273.15);
    params.referenceTemperature = T(1273.15);
    params.referenceViscosity = T(5.e11);        // Pa hr
    params.creepActivationEnergy = T(0.);        // disabled unless data supplied
    params.poissonRatio = T(0.27);
    params.relaxation = T(0.9);
    params.tolerance = T(5.e-6);
    return params;
  }
};

} // namespace viennals
