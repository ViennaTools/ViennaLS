#ifndef LAX_FRIEDRICHS_SCALAR_HPP
#define LAX_FRIEDRICHS_SCALAR_HPP

#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsDomain.hpp>
#include <lsExpand.hpp>

namespace lsInternal {
template <class T, int D, int order> class lsLaxFriedrichs {
  lsDomain<T, D> &levelSet;
  hrleSparseStarIterator<hrleDomain<T, D>> neighborIterator;
  bool calculateNormalVectors = true;
  // const double alpha;

  static T pow2(const T &value) { return value * value; }

public:
  static void prepareLS(lsDomain<T, D> &passedlsDomain) {
    assert(order == 1 || order == 2);
    lsExpand<T, D>(passedlsDomain, 2 * order + 1).apply();
  }

  lsLaxFriedrichs(lsDomain<T, D> &passedlsDomain,
                  bool calcNormal = true) //, double a = 0.1)
      : levelSet(passedlsDomain),
        neighborIterator(hrleSparseStarIterator<hrleDomain<T, D>>(
            levelSet.getDomain(), order)),
        calculateNormalVectors(calcNormal) /*,
alpha(a)*/
  {
    levelSet.calculateActivePointIds();
  }

  T operator()(const hrleVectorType<hrleIndexType, D> &indices,
               lsVelocityField<T> *velocities, int material) {

    auto &grid = levelSet.getGrid();
    double gridDelta = grid.getGridDelta();

    hrleVectorType<T, 3> coordinate(0., 0., 0.);
    for (unsigned i = 0; i < D; ++i) {
      coordinate[i] = indices[i] * gridDelta;
    }

    // move neighborIterator to current position
    neighborIterator.goToIndicesSequential(indices);

    T gradPos[D];
    T gradNeg[D];

    T grad = 0.;
    T dissipation = 0.;

    hrleVectorType<T, D> normal;
    T modulus = 0;

    for (int i = 0; i < D; i++) { // iterate over dimensions

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

        const T diff00 =
            (((deltaNeg * phiPos - deltaPos * phiNeg) / (deltaPos - deltaNeg) +
              phi0)) /
            (deltaPos * deltaNeg);
        const T phiPosPos =
            neighborIterator.getNeighbor((D * order) + i).getValue();
        const T phiNegNeg =
            neighborIterator.getNeighbor((D * order) + D + i).getValue();

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

      normal[i] = (diffNeg + diffPos) * 0.5;
      modulus += normal[i] * normal[i];

      grad += pow2((diffNeg + diffPos) * 0.5);
      dissipation += (diffPos - diffNeg) * 0.5;
    }

    // Calculate normal vector for velocity calculation
    hrleVectorType<T, 3> normalVector(T(0));
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

    double scalarVelocity =
        velocities->getScalarVelocity(coordinate, material, normalVector);
    hrleVectorType<T, 3> vectorVelocity =
        velocities->getVectorVelocity(coordinate, material, normalVector);

    T alpha = 0;
    modulus = std::sqrt(modulus);
    for (unsigned i = 0; i < D; ++i) {
      alpha += std::abs(scalarVelocity * normal[i] / modulus);
    }

    T totalGrad = 0.;

    if (scalarVelocity != 0.) {
      totalGrad = scalarVelocity * std::sqrt(grad);
    }

    for (int w = 0; w < D; w++) {
      if (vectorVelocity[w] > 0.) {
        totalGrad += vectorVelocity[w] * gradPos[w];
      } else {
        totalGrad += vectorVelocity[w] * gradNeg[w];
      }
    }

    return totalGrad - ((totalGrad != 0.) ? dissipation : 0);
  }
};
} // namespace lsInternal

#endif // LAX_FRIEDRICHS_SCALAR_HPP
