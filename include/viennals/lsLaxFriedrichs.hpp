#pragma once

#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsDomain.hpp>
#include <lsExpand.hpp>

#include <vcVectorUtil.hpp>

namespace lsInternal {

using namespace viennacore;

/// Lax Friedrichs integration scheme with constant alpha
/// value for dissipation. This alpha value should be fitted
/// based on the results of the advection and passed to the
/// advection Kernel.
template <class T, int D, int order> class LaxFriedrichs {
  SmartPointer<viennals::Domain<T, D>> levelSet;
  SmartPointer<viennals::VelocityField<T>> velocities;
  hrleSparseStarIterator<hrleDomain<T, D>, order> neighborIterator;
  bool calculateNormalVectors = true;
  const double alphaFactor = 1.0;
  hrleVectorType<T, 3> finalAlphas;

  static T pow2(const T &value) { return value * value; }

public:
  static const int order_ = order;
  static void prepareLS(SmartPointer<viennals::Domain<T, D>> passedlsDomain) {
    assert(order == 1 || order == 2);
    viennals::Expand<T, D>(passedlsDomain, 2 * order + 1).apply();
  }

  LaxFriedrichs(SmartPointer<viennals::Domain<T, D>> passedlsDomain,
                SmartPointer<viennals::VelocityField<T>> vel, double a = 1.0,
                bool calcNormal = true)
      : levelSet(passedlsDomain), velocities(vel),
        neighborIterator(hrleSparseStarIterator<hrleDomain<T, D>, order>(
            levelSet->getDomain())),
        calculateNormalVectors(calcNormal), alphaFactor(a) {
    for (int i = 0; i < 3; ++i) {
      finalAlphas[i] = 1.0;
    }
  }

  void setFinalAlphas(const hrleVectorType<T, 3> &alphas) {
    finalAlphas = alphas;
  }

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

    Vec3D<T> normalVector = {};
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
      dissipation += finalAlphas[i] * (diffPos - diffNeg) * 0.5;
    }

    if (calcNormals) {
      normalModulus = std::sqrt(normalModulus);
      for (unsigned i = 0; i < D; ++i) {
        normalVector[i] /= normalModulus;
      }
    }

    // convert coordinate to std array for interface
    Vec3D<T> coordArray = {coordinate[0], coordinate[1], coordinate[2]};

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

    return totalGrad - ((totalGrad != 0.) ? alphaFactor * dissipation : 0);
  }

  void reduceTimeStepHamiltonJacobi(double &MaxTimeStep,
                                    hrleCoordType gridDelta) {
    const double alpha_maxCFL = 1.0;
    // second time step test, based on alphas

    double timeStep = 0;
    for (int i = 0; i < D; ++i) {
      timeStep += finalAlphas[i] / gridDelta;
    }

    timeStep = alpha_maxCFL / timeStep;
    MaxTimeStep = std::min(timeStep, MaxTimeStep);
  }

  hrleVectorType<T, 3>
  calcAlpha(const hrleVectorType<hrleIndexType, D> &indices, int material) {
    const T gridDelta = levelSet->getGrid().getGridDelta();

    // move neighborIterator to current position
    neighborIterator.goToIndicesSequential(indices);

    Vec3D<T> coords;
    for (unsigned i = 0; i < D; ++i) {
      coords[i] = indices[i] * gridDelta;
    }

    const T deltaPos = gridDelta;
    const T deltaNeg = -gridDelta;

    Vec3D<T> normal = {};
    T normalModulus = 0.;
    const T phi0 = neighborIterator.getCenter().getValue();
    for (unsigned i = 0; i < D; ++i) {
      const T phiPos = neighborIterator.getNeighbor(i).getValue();
      const T phiNeg = neighborIterator.getNeighbor(i + D).getValue();

      T diffPos = (phiPos - phi0) / deltaPos;
      T diffNeg = (phiNeg - phi0) / deltaNeg;

      normal[i] = (diffNeg + diffPos) * 0.5;
      normalModulus += normal[i] * normal[i];
    }
    normalModulus = std::sqrt(normalModulus);
    // normalise normal vector
    for (unsigned i = 0; i < D; ++i)
      normal[i] /= normalModulus;

    T scaVel = velocities->getScalarVelocity(
        coords, material, normal, neighborIterator.getCenter().getPointId());
    auto vecVel = velocities->getVectorVelocity(
        coords, material, normal, neighborIterator.getCenter().getPointId());

    hrleVectorType<T, 3> alpha(0., 0., 0.);
    for (unsigned i = 0; i < D; ++i) {
      T tempAlpha = std::abs((scaVel + vecVel[i]) * normal[i]);
      alpha[i] = std::max(alpha[i], tempAlpha);
    }

    return alpha;
  }
};

namespace advect {
using namespace viennals;
template <class IntegrationSchemeType, class T, int D, int order>
void findGlobalAlpha(IntegrationSchemeType &integrationScheme,
                     std::vector<SmartPointer<Domain<T, D>>> &levelSets,
                     SmartPointer<VelocityField<T>> velocities) {

  auto &topDomain = levelSets.back()->getDomain();
  auto &grid = levelSets.back()->getGrid();

  const T gridDelta = grid.getGridDelta();
  const T deltaPos = gridDelta;
  const T deltaNeg = -gridDelta;

  hrleVectorType<T, 3> finalAlphas(0., 0., 0.);

#pragma omp parallel num_threads((levelSets.back())->getNumberOfSegments())
  {
    int p = 0;
#ifdef _OPENMP
    p = omp_get_thread_num();
#endif

    hrleVectorType<T, 3> localAlphas(0., 0., 0.);

    hrleVectorType<hrleIndexType, D> startVector =
        (p == 0) ? grid.getMinGridPoint() : topDomain.getSegmentation()[p - 1];

    hrleVectorType<hrleIndexType, D> endVector =
        (p != static_cast<int>(topDomain.getNumberOfSegments() - 1))
            ? topDomain.getSegmentation()[p]
            : grid.incrementIndices(grid.getMaxGridPoint());

    // an iterator for each level set
    std::vector<hrleSparseIterator<typename Domain<T, D>::DomainType>>
        iterators;
    for (auto it = levelSets.begin(); it != levelSets.end(); ++it) {
      iterators.push_back(hrleSparseIterator<typename Domain<T, D>::DomainType>(
          (*it)->getDomain()));
    }

    // neighborIterator for the top level set
    hrleSparseStarIterator<hrleDomain<T, D>, order> neighborIterator(topDomain);

    for (hrleSparseIterator<typename Domain<T, D>::DomainType> it(topDomain,
                                                                  startVector);
         it.getStartIndices() < endVector; ++it) {

      if (!it.isDefined() || std::abs(it.getValue()) > 0.5)
        continue;

      T value = it.getValue();

      // check if there is any other levelset at the same point:
      // if yes, take the velocity of the lowest levelset
      for (unsigned lowerLevelSetId = 0; lowerLevelSetId < levelSets.size();
           ++lowerLevelSetId) {
        // put iterator to same position as the top levelset
        auto indices = it.getStartIndices();
        iterators[lowerLevelSetId].goToIndicesSequential(indices);

        // if the lower surface is actually outside, i.e. its LS value
        // is lower or equal
        if (iterators[lowerLevelSetId].getValue() <= value + 1e-4) {

          // move neighborIterator to current position
          neighborIterator.goToIndicesSequential(indices);

          Vec3D<T> coords;
          for (unsigned i = 0; i < D; ++i) {
            coords[i] = indices[i] * gridDelta;
          }

          Vec3D<T> normal = {};
          T normalModulus = 0.;
          const T phi0 = neighborIterator.getCenter().getValue();
          for (unsigned i = 0; i < D; ++i) {
            const T phiPos = neighborIterator.getNeighbor(i).getValue();
            const T phiNeg = neighborIterator.getNeighbor(i + D).getValue();

            T diffPos = (phiPos - phi0) / deltaPos;
            T diffNeg = (phiNeg - phi0) / deltaNeg;

            normal[i] = (diffNeg + diffPos) * 0.5;
            normalModulus += normal[i] * normal[i];
          }
          normalModulus = std::sqrt(normalModulus);
          for (unsigned i = 0; i < D; ++i)
            normal[i] /= normalModulus;

          T scaVel = velocities->getScalarVelocity(
              coords, lowerLevelSetId, normal,
              neighborIterator.getCenter().getPointId());
          auto vecVel = velocities->getVectorVelocity(
              coords, lowerLevelSetId, normal,
              neighborIterator.getCenter().getPointId());

          for (unsigned i = 0; i < D; ++i) {
            T tempAlpha = std::abs((scaVel + vecVel[i]) * normal[i]);
            localAlphas[i] = std::max(localAlphas[i], tempAlpha);
          }

          break;
        }
      }
    }

#pragma omp critical
    {
      for (unsigned i = 0; i < D; ++i) {
        finalAlphas[i] = std::max(finalAlphas[i], localAlphas[i]);
      }
    }
  } // end of parallel section

  integrationScheme.setFinalAlphas(finalAlphas);
}
} // namespace advect
} // namespace lsInternal
