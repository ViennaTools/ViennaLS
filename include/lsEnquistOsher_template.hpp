#ifndef LS_ENQUIST_OSHER_SV_TEMPLATE_HPP
#define LS_ENQUIST_OSHER_SV_TEMPLATE_HPP

#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsDomain_template.hpp>
#include <lsExpand_template.hpp>
#include <lsVelocityField_template.hpp>

namespace lsInternal {

template <class T, int D, int order> class lsEnquistOsher {
  lsDomain<T, D> &levelSet;
  hrleSparseStarIterator<hrleDomain<T, D>> neighborIterator;
  bool calculateNormalVectors = true;

  static T pow2(const T &value) { return value * value; }

public:
  static void prepareLS(lsDomain<T, D> &passedlsDomain) {
    assert(order == 1 || order == 2);
    lsExpand<T, D>(passedlsDomain).apply(2 * order + 1);
  }

  lsEnquistOsher(lsDomain<T, D> &passedlsDomain, bool calcNormal = true)
      : levelSet(passedlsDomain),
        neighborIterator(hrleSparseStarIterator<hrleDomain<T, D>>(
            levelSet.getDomain(), order)),
        calculateNormalVectors(calcNormal) {
    levelSet.calculateActivePointIds();
  }

  T operator()(const hrleVectorType<hrleIndexType, D> &indices,
               lsVelocityField<T> *velocities, int material) {
    auto &grid = levelSet.getGrid();
    double gridDelta = grid.getGridDelta();

    hrleVectorType<double, 3> coordinate(0., 0., 0.);
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

    return vel_grad;
  }
};

} // namespace lsInternal

#endif // LS_ENQUIST_OSHER_SV_TEMPLATE_HPP
