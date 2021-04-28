#ifndef LS_LOCAL_LAX_FRIEDRICHS_HPP
#define LS_LOCAL_LAX_FRIEDRICHS_HPP

#include <hrleSparseBoxIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsDomain.hpp>
#include <lsExpand.hpp>

namespace lsInternal {

/// Lax Friedrichs integration scheme, which uses a first neighbour
/// stencil to calculate the alpha values for all neighbours.
/// The largest alpha value is then chosen for dissipation.
/// Slower than lsLocalLocalLaxFriedrichs or lsEngquistOsher
/// but more reliable for complex velocity fields.
template <class T, int D, int order> class lsLocalLaxFriedrichs {
  lsSmartPointer<lsDomain<T, D>> levelSet;
  lsSmartPointer<lsVelocityField<T>> velocities;
  hrleSparseBoxIterator<hrleDomain<T, D>> neighborIterator;
  const double alphaFactor;

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
  lsLocalLaxFriedrichs(lsSmartPointer<lsDomain<T, D>> passedlsDomain,
                       lsSmartPointer<lsVelocityField<T>> vel, double a = 1.0)
      : levelSet(passedlsDomain), velocities(vel),
        neighborIterator(levelSet->getDomain(), 2), alphaFactor(a) {}

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
      const unsigned numNeighbors = std::pow((maxIndex - minIndex) - 1, D);

      hrleVectorType<hrleIndexType, D> neighborIndex(minIndex);
      for (unsigned i = 0; i < numNeighbors; ++i) {
        std::array<T, 3> coords;
        for (unsigned dir = 0; dir < D; ++dir) {
          coords[dir] = coordinate[dir] + neighborIndex[dir] * gridDelta;
        }
        std::array<T, 3> normal = {};
        double normalModulus = 0.;
        auto center = neighborIterator.getNeighbor(neighborIndex).getValue();
        for (unsigned dir = 0; dir < D; ++dir) {
          hrleVectorType<hrleIndexType, D> unity(0);
          unity[dir] = 1;
          auto neg =
              neighborIterator.getNeighbor(neighborIndex - unity).getValue();
          auto pos =
              neighborIterator.getNeighbor(neighborIndex + unity).getValue();
          normal[dir] = calculateNormalComponent(neg, center, pos, gridDelta);
          normalModulus += normal[dir] * normal[dir];
        }
        normalModulus = std::sqrt(normalModulus);

        T scaVel = velocities->getScalarVelocity(
            coords, material, normal,
            neighborIterator.getCenter().getPointId());
        auto vecVel = velocities->getVectorVelocity(
            coords, material, normal,
            neighborIterator.getCenter().getPointId());

        for (unsigned dir = 0; dir < D; ++dir) {
          // normalise normal vector
          normal[i] /= normalModulus;
          T tempAlpha = std::abs((scaVel + vecVel[i]) * normal[i]);
          alpha[dir] = std::max(alpha[dir], tempAlpha);
        }

        // advance to next index
        incrementIndices(neighborIndex, minIndex, maxIndex);
      }
    }

    // calculate local dissipation alphas for each direction
    // and add to dissipation term
    for (unsigned i = 0; i < D; ++i) {
      dissipation += alphaFactor * alpha[i] * (gradNeg[i] - gradPos[i]) * 0.5;
    }

    // std::cout << neighborIterator.getCenter().getPointId() << " dissipation:
    // " << dissipation << std::endl;
    return totalGrad - ((totalGrad != 0.) ? dissipation : 0);
  }
};
} // namespace lsInternal

#endif // LS_LOCAL_LAX_FRIEDRICHS_HPP
