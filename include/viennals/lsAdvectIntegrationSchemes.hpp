#pragma once

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>

namespace viennals {

/// Enumeration for the different spatial discretization schemes
/// used by the advection kernel
enum class SpatialSchemeEnum : unsigned {
  ENGQUIST_OSHER_1ST_ORDER = 0,
  ENGQUIST_OSHER_2ND_ORDER = 1,
  LAX_FRIEDRICHS_1ST_ORDER = 2,
  LAX_FRIEDRICHS_2ND_ORDER = 3,
  LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER = 4,
  LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER = 5,
  LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER = 6,
  LOCAL_LAX_FRIEDRICHS_1ST_ORDER = 7,
  LOCAL_LAX_FRIEDRICHS_2ND_ORDER = 8,
  STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER = 9,
  WENO_5TH_ORDER = 10
};

// Legacy naming (deprecated, will be removed in future versions)
using IntegrationSchemeEnum [[deprecated("Use SpatialSchemeEnum instead")]] =
    SpatialSchemeEnum;

/// Enumeration for the different time integration schemes
/// used to select the advection kernel
enum class TemporalSchemeEnum : unsigned {
  FORWARD_EULER = 0,
  RUNGE_KUTTA_2ND_ORDER = 1,
  RUNGE_KUTTA_3RD_ORDER = 2
};

// Forward declaration
template <class T, int D> class Advect;
} // namespace viennals

namespace lsInternal {

template <class T, int D> struct AdvectTimeIntegration {
  using AdvectType = viennals::Advect<T, D>;

  static double evolveForwardEuler(AdvectType &kernel, double maxTimeStep) {
    if (kernel.currentTimeStep < 0. || kernel.storedRates.empty())
      kernel.computeRates(maxTimeStep);

    kernel.updateLevelSet(kernel.currentTimeStep);

    kernel.rebuildLS();

    kernel.adjustLowerLayers();

    return kernel.currentTimeStep;
  }

  static double evolveRungeKutta2(AdvectType &kernel, double maxTimeStep) {
    // TVD Runge-Kutta 2nd Order (Heun's Method)
    // 1. Determine time step
    kernel.computeRates(maxTimeStep);
    const double dt = kernel.getCurrentTimeStep();

    // 2. Save u^n
    if (kernel.originalLevelSet == nullptr) {
      kernel.originalLevelSet =
          viennals::Domain<T, D>::New(kernel.levelSets.back()->getGrid());
    }
    kernel.originalLevelSet->deepCopy(kernel.levelSets.back());

    if (dt <= 0)
      return 0.;

    // Stage 1: u^(1) = u^n + dt * L(u^n)
    kernel.updateLevelSet(dt);

    // Stage 2: u^(n+1) = 1/2 u^n + 1/2 (u^(1) + dt * L(u^(1)))
    // Current level set is u^(1). Compute L(u^(1)).
    kernel.computeRates(dt);
    // Update to u* = u^(1) + dt * L(u^(1))
    kernel.updateLevelSet(dt);
    // Combine: u^(n+1) = 0.5 * u^n + 0.5 * u*
    kernel.combineLevelSets(0.5, 0.5);

    // Finalize
    kernel.rebuildLS();

    kernel.adjustLowerLayers();

    return dt;
  }

  static double evolveRungeKutta3(AdvectType &kernel, double maxTimeStep) {
    // 1. Determine the single time step 'dt' for all stages.
    kernel.computeRates(maxTimeStep);
    const double dt = kernel.getCurrentTimeStep();

    // 2. Save u^n (Deep copy to preserve topology)
    if (kernel.originalLevelSet == nullptr) {
      kernel.originalLevelSet =
          viennals::Domain<T, D>::New(kernel.levelSets.back()->getGrid());
    }
    kernel.originalLevelSet->deepCopy(kernel.levelSets.back());

    // If dt is 0 or negative, no advection is possible or needed.
    if (dt <= 0)
      return 0.;

    // Stage 1: u^(1) = u^n + dt * L(u^n)
    kernel.updateLevelSet(dt);

    // Stage 2: u^(2) = 3/4 u^n + 1/4 (u^(1) + dt * L(u^(1)))
    kernel.computeRates(dt);
    kernel.updateLevelSet(dt);
    // Combine to get u^(2) = 0.75 * u^n + 0.25 * u*.
    kernel.combineLevelSets(0.75, 0.25);

    // Stage 3: u^(n+1) = 1/3 u^n + 2/3 (u^(2) + dt * L(u^(2)))
    kernel.computeRates(dt);
    kernel.updateLevelSet(dt);
    // Combine to get u^(n+1) = 1/3 * u^n + 2/3 * u**.
    kernel.combineLevelSets(1.0 / 3.0, 2.0 / 3.0);

    // Finalize: Re-segment and renormalize the final result.
    kernel.rebuildLS();

    kernel.adjustLowerLayers();

    return dt;
  }
};

} // namespace lsInternal
