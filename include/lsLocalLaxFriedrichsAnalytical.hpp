#ifndef LS_LOCAL_LAX_FRIEDRICHS_ANALYTICAL_HPP
#define LS_LOCAL_LAX_FRIEDRICHS_ANALYTICAL_HPP

#include <hrleSparseBoxIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsDomain.hpp>
#include <lsExpand.hpp>

namespace lsInternal {

/// Lax Friedrichs integration scheme, which uses alpha values
/// provided by the user in getDissipationAlphas in lsVelocityField.
/// If it is possible to derive analytical solutions for the velocityField
/// and the alpha values, this integration scheme should be used and never
/// otherwise.
template <class T, int D, int order> class lsLocalLaxFriedrichsAnalytical {
  lsSmartPointer<lsDomain<T, D>> levelSet;
  lsSmartPointer<lsVelocityField<T>> velocities;
  hrleSparseBoxIterator<hrleDomain<T, D>> neighborIterator;

  static T pow2(const T &value) { return value * value; }

  T calculateNormalComponent(T neg, T center, T pos, T delta) {
    auto diffPos = (pos - center) / delta;
    auto diffNeg = (center - neg) / delta;
    return (diffPos + diffNeg) * 0.5;
  }

  void incrementIndices(hrleVectorType<hrleIndexType, D> &index,
                        hrleIndexType minIndex, hrleIndexType maxIndex) {
    unsigned dir = 0;
    for (; dir < D - 1; ++dir) {
      if (index[dir] < maxIndex)
        break;
      index[dir] = minIndex;
    }
    ++index[dir];
  }

public:
  static void prepareLS(lsSmartPointer<lsDomain<T, D>> passedlsDomain) {
    assert(order == 1 || order == 2);
    // at least order+1 layers since we need neighbor neighbors for
    // dissipation alpha calculation
    lsExpand<T, D>(passedlsDomain, 2 * (order + 2) + 1).apply();
  }

  // neighboriterator always needs order 2 for alpha calculation
  lsLocalLaxFriedrichsAnalytical(lsSmartPointer<lsDomain<T, D>> passedlsDomain,
                                 lsSmartPointer<lsVelocityField<T>> vel)
      : levelSet(passedlsDomain), velocities(vel),
        neighborIterator(levelSet->getDomain(), 2) {}

  T operator()(const hrleVectorType<hrleIndexType, D> &indices, int material) {

    auto &grid = levelSet->getGrid();
    double gridDelta = grid.getGridDelta();

    hrleVectorType<T, 3> coordinate(0., 0., 0.);
    for (unsigned i = 0; i < D; ++i) {
      coordinate[i] = indices[i] * gridDelta;
    }

    // move neighborIterator to current position
    neighborIterator.goToIndicesSequential(indices);

    // convert coordinate to std array for interface
    std::array<T, 3> coordArray = {coordinate[0], coordinate[1], coordinate[2]};

    T gradPos[D];
    T gradNeg[D];

    T grad = 0.;
    T dissipation = 0.;

    std::array<T, 3> normalVector = {};
    T normalModulus = 0;

    for (int i = 0; i < D; i++) { // iterate over dimensions
      hrleVectorType<hrleIndexType, D> posUnit(0);
      hrleVectorType<hrleIndexType, D> negUnit(0);
      posUnit[i] = 1;
      negUnit[i] = -1;

      const T deltaPos = gridDelta;
      const T deltaNeg = -gridDelta;

      const T phi0 = neighborIterator.getCenter().getValue();
      const T phiPos = neighborIterator.getNeighbor(posUnit).getValue();
      const T phiNeg = neighborIterator.getNeighbor(negUnit).getValue();

      T diffPos = (phiPos - phi0) / deltaPos;
      T diffNeg = (phiNeg - phi0) / deltaNeg;

      if (order == 2) { // if second order time integration scheme is used
        posUnit[i] = 2;
        negUnit[i] = -2;

        const T deltaPosPos = 2 * gridDelta;
        const T deltaNegNeg = -2 * gridDelta;

        const T diff00 =
            (((deltaNeg * phiPos - deltaPos * phiNeg) / (deltaPos - deltaNeg) +
              phi0)) /
            (deltaPos * deltaNeg);
        const T phiPosPos = neighborIterator.getNeighbor(posUnit).getValue();
        const T phiNegNeg = neighborIterator.getNeighbor(negUnit).getValue();

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

      normalVector[i] = (diffNeg + diffPos) * 0.5;
      normalModulus += normalVector[i] * normalVector[i];

      grad += pow2((diffNeg + diffPos) * 0.5);
    }

    // normalise normal vector
    normalModulus = std::sqrt(normalModulus);
    for (unsigned i = 0; i < D; ++i) {
      normalVector[i] /= normalModulus;
    }

    // Get velocities
    double scalarVelocity = velocities->getScalarVelocity(
        coordArray, material, normalVector,
        neighborIterator.getCenter().getPointId());
    std::array<T, 3> vectorVelocity = velocities->getVectorVelocity(
        coordArray, material, normalVector,
        neighborIterator.getCenter().getPointId());

    // calculate hamiltonian
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

    // calculate alphas
    T alpha[D] = {};
    {
      // alpha calculation is always on order 1 stencil
      const hrleIndexType minIndex = -1;
      const hrleIndexType maxIndex = 1;
      const unsigned numNeighbors = std::pow((maxIndex - minIndex), D);

      hrleVectorType<hrleIndexType, D> neighborIndex(minIndex);
      for (unsigned i = 0; i < numNeighbors; ++i) {

        std::array<T, 3> normal = {};
        auto center = neighborIterator.getNeighbor(neighborIndex).getValue();
        for (unsigned dir = 0; dir < D; ++dir) {
          hrleVectorType<hrleIndexType, D> unity(0);
          unity[dir] = 1;
          auto neg =
              neighborIterator.getNeighbor(neighborIndex - unity).getValue();
          auto pos =
              neighborIterator.getNeighbor(neighborIndex + unity).getValue();
          normal[dir] = calculateNormalComponent(neg, center, pos, gridDelta);
        }

        for (unsigned dir = 0; dir < D; ++dir) {
          T tempAlpha = velocities->getDissipationAlpha(dir, material, normal);
          alpha[dir] = std::max(alpha[dir], tempAlpha);
        }

        // advance to next index
        incrementIndices(neighborIndex, minIndex, maxIndex);
      }
    }

    // calculate local dissipation alphas for each direction
    // and add to dissipation term
    for (unsigned i = 0; i < D; ++i) {
      dissipation += alpha[i] * (gradNeg[i] - gradPos[i]) * 0.5;
    }

    // std::cout << neighborIterator.getCenter().getPointId() << " dissipation:
    // " << dissipation << std::endl;
    return totalGrad - ((totalGrad != 0.) ? dissipation : 0);
  }
};
} // namespace lsInternal

#endif // LS_LOCAL_LAX_FRIEDRICHS_ANALYTICAL_HPP
