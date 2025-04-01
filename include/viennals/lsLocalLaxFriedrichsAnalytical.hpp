#pragma once

#include <hrleSparseBoxIterator.hpp>

#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsVelocityField.hpp>

#include <vcVectorType.hpp>

namespace lsInternal {

using namespace viennacore;

/// Lax Friedrichs integration scheme, which uses alpha values
/// provided by the user in getDissipationAlphas in lsVelocityField.
/// If it is possible to derive analytical solutions for the velocityField
/// and the alpha values, this integration scheme should be used and never
/// otherwise.
template <class T, int D, int order> class LocalLaxFriedrichsAnalytical {
  SmartPointer<viennals::Domain<T, D>> levelSet;
  SmartPointer<viennals::VelocityField<T>> velocities;
  // neighbor iterator always needs order 2 for alpha calculation
  viennahrle::SparseBoxIterator<viennahrle::Domain<T, D>, 2> neighborIterator;
  VectorType<T, 3> finalAlphas;

  static T pow2(const T &value) { return value * value; }

  T calculateNormalComponent(T neg, T center, T pos, T delta) {
    auto diffPos = (pos - center) / delta;
    auto diffNeg = (center - neg) / delta;
    return (diffPos + diffNeg) * 0.5;
  }

  static void incrementIndices(viennahrle::Index<D> &index,
                               viennahrle::IndexType minIndex,
                               viennahrle::IndexType maxIndex) {
    unsigned dir = 0;
    for (; dir < D - 1; ++dir) {
      if (index[dir] < maxIndex)
        break;
      index[dir] = minIndex;
    }
    ++index[dir];
  }

public:
  static void prepareLS(SmartPointer<viennals::Domain<T, D>> passedlsDomain) {
    assert(order == 1 || order == 2);
    // at least order+1 layers since we need neighbor neighbors for
    // dissipation alpha calculation
    viennals::Expand<T, D>(passedlsDomain, 2 * (order + 2) + 1).apply();
  }

  LocalLaxFriedrichsAnalytical(
      SmartPointer<viennals::Domain<T, D>> passedlsDomain,
      SmartPointer<viennals::VelocityField<T>> vel)
      : levelSet(passedlsDomain), velocities(vel),
        neighborIterator(levelSet->getDomain()) {
    for (int i = 0; i < 3; ++i) {
      finalAlphas[i] = 0;
    }
  }

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

    // convert coordinate to std array for interface
    Vec3D<T> coordArray{coordinate[0], coordinate[1], coordinate[2]};

    T gradPos[D];
    T gradNeg[D];

    T grad = 0.;
    T dissipation = 0.;

    Vec3D<T> normalVector = {};
    T normalModulus = 0;

    for (int i = 0; i < D; i++) { // iterate over dimensions
      viennahrle::Index<D> posUnit(0);
      viennahrle::Index<D> negUnit(0);
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
    Vec3D<T> vectorVelocity = velocities->getVectorVelocity(
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
      constexpr viennahrle::IndexType minIndex = -1;
      constexpr viennahrle::IndexType maxIndex = 1;
      constexpr int numNeighbors = hrleUtil::pow(maxIndex - minIndex + 1, D);

      viennahrle::Index<D> neighborIndex(minIndex);
      for (unsigned i = 0; i < numNeighbors; ++i) {

        Vec3D<T> normal = {};
        auto center = neighborIterator.getNeighbor(neighborIndex).getValue();
        for (unsigned dir = 0; dir < D; ++dir) {
          viennahrle::Index<D> unity(0);
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
          finalAlphas[dir] = std::max(finalAlphas[dir], tempAlpha);
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

    return {totalGrad, ((totalGrad != 0.) ? dissipation : 0)};
  }

  void reduceTimeStepHamiltonJacobi(double &MaxTimeStep, double gridDelta) {
    constexpr double alpha_maxCFL = 1.0;
    // second time step test, based on alphas

    double timeStep = 0;
    for (int i = 0; i < D; ++i) {
      timeStep += finalAlphas[i] / gridDelta;
    }

    timeStep = alpha_maxCFL / timeStep;
    MaxTimeStep = std::min(timeStep, MaxTimeStep);
  }
};
} // namespace lsInternal
