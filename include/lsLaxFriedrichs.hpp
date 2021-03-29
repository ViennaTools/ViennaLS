#ifndef LS_LAX_FRIEDRICHS_HPP
#define LS_LAX_FRIEDRICHS_HPP

#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsDomain.hpp>
#include <lsExpand.hpp>

namespace lsInternal {

/// Lax Friedrichs integration scheme with constant alpha
/// value for dissipation. This alpha value should be fitted
/// based on the results of the advection and passed to the
/// advection Kernel.
template <class T, int D, int order> class lsLaxFriedrichs {
  lsSmartPointer<lsDomain<T, D>> levelSet;
  lsSmartPointer<lsVelocityField<T>> velocities;
  hrleSparseStarIterator<hrleDomain<T, D>> neighborIterator;
  bool calculateNormalVectors = true;
  const double alpha = 1.0;

  static T pow2(const T &value) { return value * value; }

public:
  static void prepareLS(lsSmartPointer<lsDomain<T, D>> passedlsDomain) {
    assert(order == 1 || order == 2);
    lsExpand<T, D>(passedlsDomain, 2 * order + 1).apply();
  }

  lsLaxFriedrichs(lsSmartPointer<lsDomain<T, D>> passedlsDomain,
                  lsSmartPointer<lsVelocityField<T>> vel,
                  bool calcNormal = true, double a = 1.0)
      : levelSet(passedlsDomain), velocities(vel),
        neighborIterator(hrleSparseStarIterator<hrleDomain<T, D>>(
            levelSet->getDomain(), order)),
        calculateNormalVectors(calcNormal), alpha(a) {}

  T operator()(const hrleVectorType<hrleIndexType, D> &indices, int material) {

    auto &grid = levelSet->getGrid();
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

    std::array<T, 3> normalVector = {};
    T normalModulus = 0;
    const bool calcNormals = calculateNormalVectors;

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

      if (calcNormals) {
        normalVector[i] = (diffNeg + diffPos) * 0.5;
        normalModulus += normalVector[i] * normalVector[i];
      }

      grad += pow2((diffNeg + diffPos) * 0.5);
      dissipation += (diffPos - diffNeg) * 0.5;
    }

    if (calcNormals) {
      normalModulus = std::sqrt(normalModulus);
      for (unsigned i = 0; i < D; ++i) {
        normalVector[i] /= normalModulus;
      }
    }

    // convert coordinate to std array for interface
    std::array<T, 3> coordArray = {coordinate[0], coordinate[1], coordinate[2]};

    double scalarVelocity = velocities->getScalarVelocity(
        coordArray, material, normalVector,
        neighborIterator.getCenter().getPointId());
    std::array<T, 3> vectorVelocity = velocities->getVectorVelocity(
        coordArray, material, normalVector,
        neighborIterator.getCenter().getPointId());

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

    return totalGrad - ((totalGrad != 0.) ? alpha * dissipation : 0);
  }
};
} // namespace lsInternal

#endif // LS_LAX_FRIEDRICHS_HPP
