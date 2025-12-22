#pragma once

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>

namespace viennals {

/// This class implements the Strong Stability Preserving (SSP) Runge-Kutta
/// 3rd order time integration scheme (also known as TVD RK3).
/// It performs time integration using three stages of Euler steps and convex
/// combinations to preserve stability properties.
template <class T, int D> class AdvectRungeKutta3 : public Advect<T, D> {
  using ConstSparseIterator =
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>;
  using hrleIndexType = viennahrle::IndexType;
  using Base = Advect<T, D>;
  using Base::levelSets;
  SmartPointer<Domain<T, D>> originalLevelSet = nullptr;

public:
  using Base::Base; // inherit all constructors

private:
  // Helper function for linear combination:
  // target = wTarget * target + wSource * source
  void combineLevelSets(double wTarget,
                        const SmartPointer<Domain<T, D>> &target,
                        double wSource) {

    auto &domainDest = levelSets.back()->getDomain();
    auto &grid = levelSets.back()->getGrid();

#pragma omp parallel num_threads(domainDest.getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif
      auto &segDest = domainDest.getDomainSegment(p);

      viennahrle::Index<D> start = (p == 0)
                                       ? grid.getMinGridPoint()
                                       : domainDest.getSegmentation()[p - 1];
      viennahrle::Index<D> end =
          (p != static_cast<int>(domainDest.getNumberOfSegments()) - 1)
              ? domainDest.getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType> itDest(
          domainDest, start);
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>
          itTarget(target->getDomain(), start);

      unsigned definedValueIndex = 0;
      for (; itDest.getStartIndices() < end; ++itDest) {
        if (itDest.isDefined()) {
          itTarget.goToIndicesSequential(itDest.getStartIndices());
          T valSource = itDest.getValue();
          T valTarget = itTarget.getValue();
          segDest.definedValues[definedValueIndex++] =
              wTarget * valTarget + wSource * valSource;
        }
      }
    }
  }

  double advect(double maxTimeStep) override {
    Base::timeIntegrationScheme =
        TimeIntegrationSchemeEnum::RUNGE_KUTTA_3RD_ORDER;

    // 1. Save u^n (Deep copy to preserve topology)
    if (originalLevelSet == nullptr) {
      originalLevelSet = Domain<T, D>::New(levelSets.back()->getGrid());
    }
    originalLevelSet->deepCopy(levelSets.back());

    // 2. Determine the single time step 'dt' for all stages.
    // This is the maximum stable time step for a forward Euler step from u^n.
    Base::computeRates(maxTimeStep);
    const double dt = Base::getCurrentTimeStep();

    // If dt is 0 or negative, no advection is possible or needed.
    if (dt <= 0) {
      return 0.;
    }

    // Stage 1: u^(1) = u^n + dt * L(u^n)
    // L(u^n) is already in storedRates from the computeRates call above.
    // updateLevelSet modifies levelSets.back() from u^n to u^(1).
    Base::updateLevelSet(dt);

    // Stage 2: u^(2) = 3/4 u^n + 1/4 (u^(1) + dt * L(u^(1)))
    // The current levelSets.back() is u^(1). Compute L(u^(1)).
    Base::computeRates(dt);
    // Update levelSets.back() from u^(1) to u* = u^(1) + dt * L(u^(1)).
    Base::updateLevelSet(dt);
    // Combine to get u^(2) = 0.75 * u^n + 0.25 * u*.
    // The result is written to levelSets.back().
    combineLevelSets(0.75, originalLevelSet, 0.25);

    // Stage 3: u^(n+1) = 1/3 u^n + 2/3 (u^(2) + dt * L(u^(2)))
    // The current levelSets.back() is u^(2). Compute L(u^(2)).
    Base::computeRates(dt);
    // Update levelSets.back() from u^(2) to u** = u^(2) + dt * L(u^(2)).
    Base::updateLevelSet(dt);
    // Combine to get u^(n+1) = 1/3 * u^n + 2/3 * u**.
    // The result is written to levelSets.back().
    combineLevelSets(1.0 / 3.0, originalLevelSet, 2.0 / 3.0);

    // Finalize: Re-segment and renormalize the final result.
    Base::rebuildLS();

    // Adjust lower layers
    if (Base::integrationScheme !=
        IntegrationSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER) {
      for (unsigned i = 0; i < levelSets.size() - 1; ++i) {
        BooleanOperation<T, D>(levelSets[i], levelSets.back(),
                               BooleanOperationEnum::INTERSECT)
            .apply();
      }
    }

    return dt;
  }
};

PRECOMPILE_PRECISION_DIMENSION(AdvectRungeKutta3);

} // namespace viennals