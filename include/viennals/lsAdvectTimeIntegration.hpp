#pragma once

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <vcLogger.hpp>

namespace lsInternal {

template <class T, int D, class AdvectType> struct AdvectTimeIntegration {

  static double evolveForwardEuler(AdvectType &kernel, double maxTimeStep) {
    if (kernel.currentTimeStep < 0. || kernel.storedRates.empty())
      kernel.computeRates(maxTimeStep);

    kernel.updateLevelSet(kernel.currentTimeStep);

    kernel.rebuildLS();

    // Adjust all level sets below the advected one
    if (kernel.spatialScheme !=
        viennals::SpatialSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER) {
      for (unsigned i = 0; i < kernel.levelSets.size() - 1; ++i) {
        viennals::BooleanOperation<T, D>(
            kernel.levelSets[i], kernel.levelSets.back(),
            viennals::BooleanOperationEnum::INTERSECT)
            .apply();
      }
    }

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

    if (kernel.velocityUpdateCallback) {
      if (!kernel.velocityUpdateCallback(kernel.levelSets.back())) {
        VIENNACORE_LOG_WARNING(
            "Velocity update callback returned false in RK2 stage 1.");
      }
    }

    // Stage 2: u^(n+1) = 1/2 u^n + 1/2 (u^(1) + dt * L(u^(1)))
    // Current level set is u^(1). Compute L(u^(1)).
    kernel.computeRates(dt);
    // Update to u* = u^(1) + dt * L(u^(1))
    kernel.updateLevelSet(dt);
    // Combine: u^(n+1) = 0.5 * u^n + 0.5 * u*
    kernel.combineLevelSets(0.5, 0.5);

    // Finalize
    kernel.rebuildLS();

    // Adjust lower layers
    if (kernel.spatialScheme !=
        viennals::SpatialSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER) {
      for (unsigned i = 0; i < kernel.levelSets.size() - 1; ++i) {
        viennals::BooleanOperation<T, D>(
            kernel.levelSets[i], kernel.levelSets.back(),
            viennals::BooleanOperationEnum::INTERSECT)
            .apply();
      }
    }
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

    if (kernel.velocityUpdateCallback) {
      if (!kernel.velocityUpdateCallback(kernel.levelSets.back())) {
        VIENNACORE_LOG_WARNING(
            "Velocity update callback returned false in RK3 stage 1.");
      }
    }

    // Stage 2: u^(2) = 3/4 u^n + 1/4 (u^(1) + dt * L(u^(1)))
    kernel.computeRates(dt);
    kernel.updateLevelSet(dt);
    // Combine to get u^(2) = 0.75 * u^n + 0.25 * u*.
    kernel.combineLevelSets(0.75, 0.25);

    if (kernel.velocityUpdateCallback) {
      if (!kernel.velocityUpdateCallback(kernel.levelSets.back())) {
        VIENNACORE_LOG_WARNING(
            "Velocity update callback returned false in RK3 stage 2.");
      }
    }

    // Stage 3: u^(n+1) = 1/3 u^n + 2/3 (u^(2) + dt * L(u^(2)))
    kernel.computeRates(dt);
    kernel.updateLevelSet(dt);
    // Combine to get u^(n+1) = 1/3 * u^n + 2/3 * u**.
    kernel.combineLevelSets(1.0 / 3.0, 2.0 / 3.0);

    // Finalize: Re-segment and renormalize the final result.
    kernel.rebuildLS();

    // Adjust lower layers
    if (kernel.spatialScheme !=
        viennals::SpatialSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER) {
      for (unsigned i = 0; i < kernel.levelSets.size() - 1; ++i) {
        viennals::BooleanOperation<T, D>(
            kernel.levelSets[i], kernel.levelSets.back(),
            viennals::BooleanOperationEnum::INTERSECT)
            .apply();
      }
    }

    return dt;
  }
};

} // namespace lsInternal
