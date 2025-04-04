#pragma once

#include <hrleSparseStarIterator.hpp>

#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsVelocityField.hpp>

#include <vcVectorType.hpp>

namespace lsInternal {

using namespace viennacore;

/// Engquist-Osher integration scheme based on the
/// upwind integration scheme. Offers high performance
/// but lower accuracy for complex velocity fields.
template <class T, int D, int order> class EngquistOsher {
  SmartPointer<viennals::Domain<T, D>> levelSet;
  SmartPointer<viennals::VelocityField<T>> velocities;
  viennahrle::SparseStarIterator<viennahrle::Domain<T, D>, order>
      neighborIterator;
  bool calculateNormalVectors = true;

  static T pow2(const T &value) { return value * value; }

public:
  static void prepareLS(SmartPointer<viennals::Domain<T, D>> passedlsDomain) {
    assert(order == 1 || order == 2);
    viennals::Expand<T, D>(passedlsDomain, 2 * order + 1).apply();
  }

  EngquistOsher(SmartPointer<viennals::Domain<T, D>> passedlsDomain,
                SmartPointer<viennals::VelocityField<T>> vel,
                bool calcNormal = true)
      : levelSet(passedlsDomain), velocities(vel),
        neighborIterator(
            viennahrle::SparseStarIterator<viennahrle::Domain<T, D>, order>(
                levelSet->getDomain())),
        calculateNormalVectors(calcNormal) {}

  std::pair<T, T> operator()(const viennahrle::Index<D> &indices,
                             int material) {
    auto &grid = levelSet->getGrid();
    double gridDelta = grid.getGridDelta();

    VectorType<T, 3> coordinate{0., 0., 0.};
    for (unsigned i = 0; i < D; ++i) {
      coordinate[i] = indices[i] * gridDelta;
    }

    // move neighborIterator to current position
    neighborIterator.goToIndicesSequential(indices);

    T gradPos[D];
    T gradNeg[D];

    T gradPosTotal = 0;
    T gradNegTotal = 0;

    for (int i = 0; i < D; i++) {
      const T deltaPos = gridDelta;
      const T deltaNeg = -gridDelta;

      const T phi0 = neighborIterator.getCenter().getValue();
      const T phiPos = neighborIterator.getNeighbor(i).getValue();
      const T phiNeg = neighborIterator.getNeighbor(i + D).getValue();

      T diffPos = (phiPos - phi0) / deltaPos;
      T diffNeg = (phiNeg - phi0) / deltaNeg;

      if (order == 2) { // if second order time integration scheme is used
        const T deltaPosPos = 2 * gridDelta;
        const T deltaNegNeg = -2 * gridDelta;

        const T phiPosPos =
            neighborIterator.getNeighbor((D * order) + i).getValue();
        const T phiNegNeg =
            neighborIterator.getNeighbor((D * order) + D + i).getValue();

        const T diff00 =
            (((deltaNeg * phiPos - deltaPos * phiNeg) / (deltaPos - deltaNeg) +
              phi0)) /
            (deltaPos * deltaNeg);
        const T diffNegNeg = (((deltaNeg * phiNegNeg - deltaNegNeg * phiNeg) /
                                   (deltaNegNeg - deltaNeg) +
                               phi0)) /
                             (deltaNegNeg * deltaNeg);
        const T diffPosPos = (((deltaPos * phiPosPos - deltaPosPos * phiPos) /
                                   (deltaPosPos - deltaPos) +
                               phi0)) /
                             (deltaPosPos * deltaPos);

        if (std::signbit(diff00) == std::signbit(diffPosPos)) {
          if (std::abs(diffPosPos * deltaPos) < std::abs(diff00 * deltaNeg)) {
            diffPos -= deltaPos * diffPosPos;
          } else {
            diffPos += deltaNeg * diff00;
          }
        }

        if (std::signbit(diff00) == std::signbit(diffNegNeg)) {
          if (std::abs(diffNegNeg * deltaNeg) < std::abs(diff00 * deltaPos)) {
            diffNeg -= deltaNeg * diffNegNeg;
          } else {
            diffNeg += deltaPos * diff00;
          }
        }
      }

      gradPos[i] = diffNeg;
      gradNeg[i] = diffPos;

      gradPosTotal +=
          pow2(std::max(diffNeg, T(0))) + pow2(std::min(diffPos, T(0)));
      gradNegTotal +=
          pow2(std::min(diffNeg, T(0))) + pow2(std::max(diffPos, T(0)));
    }

    T vel_grad = 0.;

    // Calculate normal vector for velocity calculation
    // use std::array since it will be exposed to interface
    Vec3D<T> normalVector = {};
    if (calculateNormalVectors) {
      T denominator = 0;
      for (int i = 0; i < D; i++) {
        T pos = neighborIterator.getNeighbor(i).getValue() -
                neighborIterator.getCenter().getValue();
        T neg = neighborIterator.getCenter().getValue() -
                neighborIterator.getNeighbor(i + D).getValue();
        normalVector[i] = (pos + neg) * 0.5; // = 0;
        denominator += normalVector[i] * normalVector[i];
      }
      denominator = std::sqrt(denominator);
      for (unsigned i = 0; i < D; ++i) {
        normalVector[i] /= denominator;
      }
    }

    double scalarVelocity = velocities->getScalarVelocity(
        coordinate, material, normalVector,
        neighborIterator.getCenter().getPointId());
    Vec3D<T> vectorVelocity = velocities->getVectorVelocity(
        coordinate, material, normalVector,
        neighborIterator.getCenter().getPointId());

    if (scalarVelocity > 0) {
      vel_grad += std::sqrt(gradPosTotal) * scalarVelocity;
    } else {
      vel_grad += std::sqrt(gradNegTotal) * scalarVelocity;
    }

    for (int w = 0; w < D; w++) {
      if (vectorVelocity[w] > 0.) {
        vel_grad += vectorVelocity[w] * gradPos[w];
      } else {
        vel_grad += vectorVelocity[w] * gradNeg[w];
      }
    }

    return {vel_grad, 0.};
  }

  void reduceTimeStepHamiltonJacobi(double &MaxTimeStep,
                                    double gridDelta) const {}
};

} // namespace lsInternal
