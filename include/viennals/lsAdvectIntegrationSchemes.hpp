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

  static double evolveForwardEuler(AdvectType &kernel, double maxTimeStep,
                                   bool updateLowerLayers = true) {
    if (kernel.currentTimeStep < 0. || kernel.storedRates.empty())
      kernel.computeRates(maxTimeStep);

    kernel.updateLevelSet(kernel.currentTimeStep);

    kernel.rebuildLS();

    if (updateLowerLayers)
      kernel.adjustLowerLayers();

    return kernel.currentTimeStep;
  }

  static double evolveRungeKutta2(AdvectType &kernel, double maxTimeStep) {
    // TVD Runge-Kutta 2nd Order (Heun's Method)

    // 1. Save u^n
    if (kernel.originalLevelSet == nullptr) {
      kernel.originalLevelSet =
          viennals::Domain<T, D>::New(kernel.levelSets.back()->getGrid());
    }
    kernel.originalLevelSet->deepCopy(kernel.levelSets.back());

    // Stage 1: u^(1) = u^n + dt * L(u^n)
    // Update lower layers only if we have a callback
    double dt1 = evolveForwardEuler(kernel, maxTimeStep,
                                    kernel.velocityUpdateCallback != nullptr);

    if (dt1 <= 0.)
      return 0.;

    if (kernel.velocityUpdateCallback)
      kernel.velocityUpdateCallback(kernel.levelSets.back());

    // Stage 2: u^(n+1) = 1/2 u^n + 1/2 (u^(1) + dt * L(u^(1)))
    // Current level set is u^(1). Compute L(u^(1)).
    // Update to u* = u^(1) + dt * L(u^(1))
    double dt2 = evolveForwardEuler(kernel, dt1, false);

    // Combine: u^(n+1) = 0.5 * u^n + 0.5 * u*
    kernel.combineLevelSets(0.5, 0.5);

    return 0.5 * dt1 + 0.5 * dt2;
  }

  static double evolveRungeKutta3(AdvectType &kernel, double maxTimeStep) {
    // 1. Save u^n (Deep copy to preserve topology)
    if (kernel.originalLevelSet == nullptr) {
      kernel.originalLevelSet =
          viennals::Domain<T, D>::New(kernel.levelSets.back()->getGrid());
    }
    kernel.originalLevelSet->deepCopy(kernel.levelSets.back());

    // Stage 1: u^(1) = u^n + dt * L(u^n)
    // This calculates dt based on u^n and advances to u^1.
    double dt1 = evolveForwardEuler(kernel, maxTimeStep,
                                    kernel.velocityUpdateCallback != nullptr);

    if (dt1 <= 0.)
      return 0.;

    if (kernel.velocityUpdateCallback)
      kernel.velocityUpdateCallback(kernel.levelSets.back());

    // Stage 2: u^(2) = 3/4 u^n + 1/4 (u^(1) + dt * L(u^(1)))
    // u* = u^(1) + dt * L(u^(1))
    double dt2 = evolveForwardEuler(kernel, dt1, false);
    // Combine to get u^(2) = 0.75 * u^n + 0.25 * u*.
    kernel.combineLevelSets(0.75, 0.25);

    if (kernel.velocityUpdateCallback) {
      kernel.velocityUpdateCallback(kernel.levelSets.back());
    }

    // Stage 3: u^(n+1) = 1/3 u^n + 2/3 (u^(2) + dt * L(u^(2)))
    // u** = u^(2) + dt * L(u^(2))
    double dt3 = evolveForwardEuler(kernel, dt1, false);

    // Combine to get u^(n+1) = 1/3 * u^n + 2/3 * u**.
    kernel.combineLevelSets(1.0 / 3.0, 2.0 / 3.0);

    return (dt1 + dt2 + 4.0 * dt3) / 6.0;
  }
};

} // namespace lsInternal
