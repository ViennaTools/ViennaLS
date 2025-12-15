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

  SmartPointer<Domain<T, D>> originalLevelSet = nullptr;

public:
  using Base::Base; // inherit all constructors

  using Base::advectedTime;
  using Base::advectionTime;
  using Base::currentTimeStep;
  using Base::levelSets;
  using Base::numberOfTimeSteps;
  using Base::performOnlySingleStep;
  using Base::storedRates;
  using Base::velocities;

private:

  // Helper function for linear combination: target = wTarget * target + wSource
  // * source
  void combineLevelSets(double wTarget, const SmartPointer<Domain<T, D>> &target,
                        double wSource, const SmartPointer<Domain<T, D>> &source) {
    // We write the result into levelSets.back() (which is passed as 'target' or
    // 'source' usually, but here we assume target is the destination)
    // Actually, to support u_new = a*u_old + b*u_new, we should write to
    // levelSets.back().
    auto &domainDest = levelSets.back()->getDomain();
    const auto &domainTarget = target->getDomain();
    const auto &domainSource = source->getDomain();

    if (domainTarget.getNumberOfSegments() !=
        domainSource.getNumberOfSegments()) {
      VIENNACORE_LOG_ERROR(
          "AdvectRungeKutta3: Topology mismatch in combineLevelSets.");
      return;
    }

#pragma omp parallel for schedule(static)
    for (int p = 0; p < static_cast<int>(domainDest.getNumberOfSegments());
         ++p) {
      auto &segDest = domainDest.getDomainSegment(p);
      const auto &segTarget = domainTarget.getDomainSegment(p);
      const auto &segSource = domainSource.getDomainSegment(p);

      if (segTarget.definedValues.size() == segSource.definedValues.size() &&
          segDest.definedValues.size() == segTarget.definedValues.size()) {
        for (size_t i = 0; i < segDest.definedValues.size(); ++i) {
          segDest.definedValues[i] = wTarget * segTarget.definedValues[i] +
                                     wSource * segSource.definedValues[i];
        }
      }
    }
  }

public:
  double advect(double maxTimeStep) override {
    if (levelSets.empty()) {
      VIENNACORE_LOG_ERROR("No level sets passed to Advect. Not advecting.");
      return std::numeric_limits<double>::max();
    }
    if (velocities == nullptr) {
      VIENNACORE_LOG_ERROR(
          "No velocity field passed to Advect. Not advecting.");
      return std::numeric_limits<double>::max();
    }

    // 1. Prepare and Expand
    Base::prepareLS();

    // 2. Save u^n (Deep copy with identical topology)
    if (originalLevelSet == nullptr) {
      originalLevelSet =
          SmartPointer<Domain<T, D>>::New(levelSets.back()->getGrid());
    }
    originalLevelSet->deepCopy(levelSets.back());

    // 3. Stage 1: u^(1) = u^n + dt * L(u^n)
    Base::computeRates(maxTimeStep);
    double dt = Base::getCurrentTimeStep();
    Base::updateLevelSet(dt);

    // 4. Stage 2: u^(2) = 3/4 u^n + 1/4 (u^(1) + dt * L(u^(1)))
    // Calculate rates based on u^(1) (current levelSets.back())
    Base::computeRates(dt);
    // Update to get u^* = u^(1) + dt * L(u^(1))
    Base::updateLevelSet(dt);
    // Combine: u^(2) = 0.75 * u^n + 0.25 * u^*
    combineLevelSets(0.75, originalLevelSet, 0.25, levelSets.back());

    // 5. Stage 3: u^(n+1) = 1/3 u^n + 2/3 (u^(2) + dt * L(u^(2)))
    // Calculate rates based on u^(2) (current levelSets.back())
    Base::computeRates(dt);
    // Update to get u^** = u^(2) + dt * L(u^(2))
    Base::updateLevelSet(dt);
    // Combine: u^(n+1) = 1/3 * u^n + 2/3 * u^**
    combineLevelSets(1.0 / 3.0, originalLevelSet, 2.0 / 3.0, levelSets.back());

    // 6. Finalize: Re-segment and renormalize only at the end
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