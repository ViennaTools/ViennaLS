#pragma once

#include <hrleSparseStarIterator.hpp>
#include <hrleTypes.hpp>

#include <lsDomain.hpp>
#include <lsExpand.hpp>

#include <vcVectorType.hpp>

namespace lsInternal {

using namespace viennacore;

/// Lax Friedrichs integration scheme with constant alpha
/// value for dissipation. This alpha value should be fitted
/// based on the results of the advection and passed to the
/// advection Kernel.
template <class T, int D, int order> class LaxFriedrichs {
  SmartPointer<viennals::Domain<T, D>> levelSet;
  SmartPointer<viennals::VelocityField<T>> velocities;
  viennahrle::SparseStarIterator<viennahrle::Domain<T, D>, order>
      neighborIterator;
  const double alphaFactor = 1.0;
  VectorType<T, 3> const finalAlphas;
  const bool calculateNormalVectors = true;

  static T pow2(const T &value) { return value * value; }

public:
  // static const int order_ = order;
  static void prepareLS(SmartPointer<viennals::Domain<T, D>> passedlsDomain) {
    assert(order == 1 || order == 2);
    viennals::Expand<T, D>(passedlsDomain, 2 * order + 1).apply();
  }

  LaxFriedrichs(SmartPointer<viennals::Domain<T, D>> passedlsDomain,
                SmartPointer<viennals::VelocityField<T>> vel, double alpha,
                VectorType<T, 3> &alphas, bool calcNormal)
      : levelSet(passedlsDomain), velocities(vel),
        neighborIterator(levelSet->getDomain()), alphaFactor(alpha),
        finalAlphas(alphas), calculateNormalVectors(calcNormal) {}

  std::pair<T, T> operator()(const viennahrle::Index<D> &indices,
                             int material) {

    auto &grid = levelSet->getGrid();
    double gridDelta = grid.getGridDelta();

    VectorType<T, 3> coordinate(0., 0., 0.);
    for (unsigned i = 0; i < D; ++i) {
      coordinate[i] = indices[i] * gridDelta;
    }

    // move neighborIterator to current position
    neighborIterator.goToIndicesSequential(indices);

    T gradPos[D];
    T gradNeg[D];

    T grad = 0.;
    T dissipation = 0.;

    Vec3D<T> normalVector = {};
    T normalModulus = 0;

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

      if (calculateNormalVectors) {
        normalVector[i] = (diffNeg + diffPos) * 0.5;
        normalModulus += normalVector[i] * normalVector[i];
      }

      grad += pow2((diffNeg + diffPos) * 0.5);
      dissipation += alphaFactor * finalAlphas[i] * (diffPos - diffNeg) * 0.5;
    }

    if (calculateNormalVectors) {
      normalModulus = std::sqrt(normalModulus);
      for (unsigned i = 0; i < D; ++i) {
        normalVector[i] /= normalModulus;
      }
    }

    // convert coordinate to std array for interface
    Vec3D<T> coordArray{coordinate[0], coordinate[1], coordinate[2]};

    double scalarVelocity = velocities->getScalarVelocity(
        coordArray, material, normalVector,
        neighborIterator.getCenter().getPointId());
    Vec3D<T> vectorVelocity = velocities->getVectorVelocity(
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

    return {totalGrad, ((totalGrad != 0.) ? dissipation : 0)};
  }

  void reduceTimeStepHamiltonJacobi(double &MaxTimeStep,
                                    double gridDelta) const {
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
