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
template <class T, int D, int order> class WENO {
  static_assert(order == 3 || order == 5, "WENO order must be 3 or 5.");

  SmartPointer<viennals::Domain<T, D>> levelSet;
  SmartPointer<viennals::VelocityField<T>> velocities;

  static constexpr int stencilRadius = (order + 1) / 2;

  // Iterator depth: WENO needs stencilRadius neighbors on each side.
  viennahrle::SparseStarIterator<viennahrle::Domain<T, D>, stencilRadius>
      neighborIterator;

  const bool calculateNormalVectors = true;

  // Use the existing math engine with WENO scheme
  using MathScheme =
      FiniteDifferences<T, (order == 3) ? DifferentiationSchemeEnum::WENO3
                                        : DifferentiationSchemeEnum::WENO5>;

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

    T gradPosTotal = 0;
    T gradNegTotal = 0;

    // --- OPTIMIZATION: Store derivatives to avoid re-calculation ---
    T wenoGradMinus[D]; // Approximates derivative from left (phi_x^-)
    T wenoGradPlus[D];  // Approximates derivative from right (phi_x^+)

    // Array to hold the stencil values
    T stencil[2 * stencilRadius + 1];

    for (int i = 0; i < D; i++) {
      // 1. GATHER STENCIL
      // We map the SparseStarIterator (which uses encoded directions)
      // to the flat array expected by FiniteDifferences.

      // Center
      stencil[stencilRadius] = neighborIterator.getCenter().getValue();

      for (int k = 1; k <= stencilRadius; ++k) {
        stencil[stencilRadius + k] =
            neighborIterator.getNeighbor((k - 1) * 2 * D + i).getValue();
        stencil[stencilRadius - k] =
            neighborIterator.getNeighbor((k - 1) * 2 * D + D + i).getValue();
      }

      // 2. COMPUTE DERIVATIVES
      // Delegate the math to your existing library and store results
      wenoGradMinus[i] = MathScheme::differenceNegative(stencil, gridDelta);
      wenoGradPlus[i] = MathScheme::differencePositive(stencil, gridDelta);

      // 3. GODUNOV FLUX PREPARATION
      // We accumulate the gradient magnitudes for the scalar velocity case.
      // This is part of the standard level set equation logic (Osher-Sethian).

      // For Positive Scalar Velocity (Deposition): Use upwind selection
      gradPosTotal += pow2(std::max(wenoGradMinus[i], T(0))) +
                      pow2(std::min(wenoGradPlus[i], T(0)));

      // For Negative Scalar Velocity (Etching): Use upwind selection
      gradNegTotal += pow2(std::min(wenoGradMinus[i], T(0))) +
                      pow2(std::max(wenoGradPlus[i], T(0)));
    }

    T vel_grad = 0.;

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

    // --- Retrieve Velocity ---
    double scalarVelocity = velocities->getScalarVelocity(
        coordinate, material, normalVector,
        neighborIterator.getCenter().getPointId());
    Vec3D<T> vectorVelocity = velocities->getVectorVelocity(
        coordinate, material, normalVector,
        neighborIterator.getCenter().getPointId());

    // --- Apply Velocities ---

    // Scalar term (Etching/Deposition)
    if (scalarVelocity > 0) {
      vel_grad += std::sqrt(gradPosTotal) * scalarVelocity;
    } else {
      vel_grad += std::sqrt(gradNegTotal) * scalarVelocity;
    }

    // Vector term (Advection)
    // Here we REUSE the derivatives stored in wenoGradMinus/Plus.
    // This is the optimization compared to recalculating.
    for (int w = 0; w < D; w++) {
      if (vectorVelocity[w] > 0.) {
        vel_grad += vectorVelocity[w] * wenoGradMinus[w];
      } else {
        vel_grad += vectorVelocity[w] * wenoGradPlus[w];
      }
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