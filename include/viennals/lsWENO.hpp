#pragma once

#include <cmath>
#include <hrleSparseStarIterator.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsVelocityField.hpp>
#include <vcVectorType.hpp>

// Include your existing math library
#include <lsFiniteDifferences.hpp>

namespace lsInternal {

using namespace viennacore;

/// Weighted Essentially Non-Oscillatory (WENO) scheme.
/// This kernel acts as the grid-interface for the mathematical logic
/// defined in lsFiniteDifferences.hpp.
template <class T, int D, DifferentiationSchemeEnum scheme> class WENO {
  static_assert(scheme == DifferentiationSchemeEnum::WENO3 ||
                    scheme == DifferentiationSchemeEnum::WENO5 ||
                    scheme == DifferentiationSchemeEnum::WENO5_Z,
                "WENO scheme must be WENO3, WENO5, or WENO5_Z.");

  static constexpr int order = (scheme == DifferentiationSchemeEnum::WENO3) ? 3 : 5;

  SmartPointer<viennals::Domain<T, D>> levelSet;
  SmartPointer<viennals::VelocityField<T>> velocities;

  static constexpr int stencilRadius = (order + 1) / 2;

  // Iterator depth: WENO needs stencilRadius neighbors on each side.
  viennahrle::SparseStarIterator<viennahrle::Domain<T, D>, stencilRadius>
      neighborIterator;

  const bool calculateNormalVectors = true;

  // Use the existing math engine with WENO scheme
  using MathScheme = FiniteDifferences<T, scheme>;

  static T pow2(const T &value) { return value * value; }

public:
  static void prepareLS(SmartPointer<viennals::Domain<T, D>> passedlsDomain) {
    // Ensure we expand enough layers to access neighbors.
    viennals::Expand<T, D>(passedlsDomain, 2 * stencilRadius + 1).apply();
  }

  WENO(SmartPointer<viennals::Domain<T, D>> passedlsDomain,
       SmartPointer<viennals::VelocityField<T>> vel, bool calcNormal = true)
      : levelSet(passedlsDomain), velocities(vel),
        neighborIterator(levelSet->getDomain()),
        calculateNormalVectors(calcNormal) {}

  std::pair<T, T> operator()(const viennahrle::Index<D> &indices,
                             int material) {
    auto &grid = levelSet->getGrid();
    double gridDelta = grid.getGridDelta();

    VectorType<T, 3> coordinate{0., 0., 0.};
    for (unsigned i = 0; i < D; ++i) {
      coordinate[i] = indices[i] * gridDelta;
    }

    // Move iterator to the current grid point
    neighborIterator.goToIndicesSequential(indices);

    // --- Standard Normal Vector Calculation (for velocity lookup) ---
    Vec3D<T> normalVector{};
    if (calculateNormalVectors) {
      T denominator = 0;
      for (int i = 0; i < D; i++) {
        // Simple Central Difference for the normal vector is sufficient
        // and robust for velocity direction lookup.
        T pos = neighborIterator.getNeighbor(i).getValue();
        T neg = neighborIterator.getNeighbor(i + D).getValue();
        normalVector[i] = (pos - neg) * 0.5;
        denominator += normalVector[i] * normalVector[i];
      }
      if (denominator > 0) {
        denominator = 1. / std::sqrt(denominator);
        for (unsigned i = 0; i < D; ++i) {
          normalVector[i] *= denominator;
        }
      }
    }

    // --- Retrieve Velocity First ---
    double scalarVelocity = velocities->getScalarVelocity(
        coordinate, material, normalVector,
        neighborIterator.getCenter().getPointId());
    Vec3D<T> vectorVelocity = velocities->getVectorVelocity(
        coordinate, material, normalVector,
        neighborIterator.getCenter().getPointId());

    T gradPosTotal = 0;
    T gradNegTotal = 0;
    T vel_grad = 0.;

    // Array to hold the stencil values
    T stencil[2 * stencilRadius + 1];

    for (int i = 0; i < D; i++) {
      // Optimization: Only calculate what is needed based on velocity direction
      bool calcMinus = (scalarVelocity != 0.) || (vectorVelocity[i] > 0.);
      bool calcPlus = (scalarVelocity != 0.) || (vectorVelocity[i] < 0.);

      if (!calcMinus && !calcPlus)
        continue;

      // 1. GATHER STENCIL (Only if needed)
      stencil[stencilRadius] = neighborIterator.getCenter().getValue();
      for (int k = 1; k <= stencilRadius; ++k) {
        stencil[stencilRadius + k] =
            neighborIterator.getNeighbor((k - 1) * 2 * D + i).getValue();
        stencil[stencilRadius - k] =
            neighborIterator.getNeighbor((k - 1) * 2 * D + D + i).getValue();
      }

      // 2. COMPUTE DERIVATIVES & FLUX
      T wenoMinus = 0.;
      T wenoPlus = 0.;

      if (calcMinus)
        wenoMinus = MathScheme::differenceNegative(stencil, gridDelta);
      if (calcPlus)
        wenoPlus = MathScheme::differencePositive(stencil, gridDelta);

      if (scalarVelocity > 0.) {
        gradPosTotal +=
            pow2(std::max(wenoMinus, T(0))) + pow2(std::min(wenoPlus, T(0)));
      } else if (scalarVelocity < 0.) {
        gradNegTotal +=
            pow2(std::min(wenoMinus, T(0))) + pow2(std::max(wenoPlus, T(0)));
      }

      if (vectorVelocity[i] > 0.)
        vel_grad += vectorVelocity[i] * wenoMinus;
      else if (vectorVelocity[i] < 0.)
        vel_grad += vectorVelocity[i] * wenoPlus;
    }

    // Finalize Scalar term
    if (scalarVelocity > 0) {
      vel_grad += std::sqrt(gradPosTotal) * scalarVelocity;
    } else if (scalarVelocity < 0) {
      vel_grad += std::sqrt(gradNegTotal) * scalarVelocity;
    }

    // WENO is an upwind scheme, so explicit dissipation is 0.
    return {vel_grad, 0.};
  }

  void reduceTimeStepHamiltonJacobi(double &MaxTimeStep,
                                    double gridDelta) const {}
  //   // --- STABILITY IMPROVEMENT ---
  //   // High-order schemes like WENO5 combined with simple time integration
  //   (like
  //   // the likely Forward Euler used in Advect) can be less stable at
  //   CFL=0.5.
  //   // We enforce a safety factor here to ensure robustness.
  //   constexpr double wenoSafetyFactor = 0.5;
  //   MaxTimeStep *= wenoSafetyFactor;
  // }
};

} // namespace lsInternal
