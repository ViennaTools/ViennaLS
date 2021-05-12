#ifndef LS_ADVECT_HPP
#define LS_ADVECT_HPP

#include <lsPreCompileMacros.hpp>

#include <limits>
#include <vector>

#include <hrleSparseIterator.hpp>
#include <hrleSparseStarIterator.hpp>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMarkVoidPoints.hpp>
#include <lsMessage.hpp>
#include <lsReduce.hpp>
#include <lsSmartPointer.hpp>

// Integration schemes
#include <lsConcepts.hpp>
#include <lsEnquistOsher.hpp>
#include <lsLaxFriedrichs.hpp>
#include <lsLocalLaxFriedrichs.hpp>
#include <lsLocalLaxFriedrichsAnalytical.hpp>
#include <lsLocalLocalLaxFriedrichs.hpp>
#include <lsStencilLocalLaxFriedrichsScalar.hpp>

// Velocity accessor
#include <lsVelocityField.hpp>

// #define DEBUG_LS_ADVECT_HPP
#ifdef DEBUG_LS_ADVECT_HPP
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>
#endif

/// Enumeration for the different Integration schemes
/// used by the advection kernel
enum struct lsIntegrationSchemeEnum : unsigned {
  ENGQUIST_OSHER_1ST_ORDER = 0,
  ENGQUIST_OSHER_2ND_ORDER = 1,
  LAX_FRIEDRICHS_1ST_ORDER = 2,
  LAX_FRIEDRICHS_2ND_ORDER = 3,
  LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER = 4,
  LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER = 5,
  LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER = 6,
  LOCAL_LAX_FRIEDRICHS_1ST_ORDER = 7,
  LOCAL_LAX_FRIEDRICHS_2ND_ORDER = 8,
  STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER = 9
};

/// This class is used to advance level sets over time.
/// Level sets are passed to the constructor in an std::vector, with
/// the last element being the level set to advect, or "top level set", while
/// the others are then adjusted afterwards. In order to ensure that advection
/// works correctly, the "top level set" has to include all lower level sets:
/// LS_top = LS_top U LS_i for i = {0 ... n}, where n is the number of level
/// sets. The velocities used to advect the level set are given in a concrete
/// implementation of the lsVelocityField (check Advection examples for
/// guidance)
template <class T, int D> class lsAdvect {
  std::vector<lsSmartPointer<lsDomain<T, D>>> levelSets;
  lsSmartPointer<lsVelocityField<T>> velocities = nullptr;
  lsIntegrationSchemeEnum integrationScheme =
      lsIntegrationSchemeEnum::ENGQUIST_OSHER_1ST_ORDER;
  double timeStepRatio = 0.4999;
  double dissipationAlpha = 1.0;
  bool calculateNormalVectors = true;
  bool ignoreVoids = false;
  double advectionTime = 0.;
  double advectedTime = 0.;
  unsigned numberOfTimeSteps = 0;
  bool saveAdvectionVelocities = false;
  static constexpr double wrappingLayerEpsilon = 1e-4;

  template <class IntegrationSchemeType,
            lsConcepts::IsNotSame<IntegrationSchemeType,
                                  lsInternal::lsStencilLocalLaxFriedrichsScalar<
                                      T, D, 1>> = lsConcepts::assignable>
  void reduceTimeStepHamiltonJacobi(IntegrationSchemeType &, double &) {}

  template <class IntegrationSchemeType,
            lsConcepts::IsSame<IntegrationSchemeType,
                               lsInternal::lsStencilLocalLaxFriedrichsScalar<
                                   T, D, 1>> = lsConcepts::assignable>
  void reduceTimeStepHamiltonJacobi(IntegrationSchemeType &scheme,
                                    double &MaxTimeStep) {
    const double alpha_maxCFL = 1.0;
    // second time step test, based on alphas
    hrleVectorType<T, 3> alphas = scheme.getFinalAlphas();

    auto gridDelta = levelSets.back()->getGrid().getGridDelta();

    double timeStep = 0;
    for (int i = 0; i < D; ++i) {
      timeStep += alphas[i] / gridDelta;
    }

    timeStep = alpha_maxCFL / timeStep;
    MaxTimeStep = std::min(timeStep, MaxTimeStep);
  }

  void rebuildLS() {
    // TODO: this function uses manhatten distances for renormalisation,
    // since this is the quickest. For visualisation applications, better
    // renormalisation is needed, so it might be good to implement
    // Euler distance renormalisation as an option
    auto &grid = levelSets.back()->getGrid();
    auto newlsDomain = lsSmartPointer<lsDomain<T, D>>::New(grid);
    typename lsDomain<T, D>::DomainType &newDomain = newlsDomain->getDomain();
    typename lsDomain<T, D>::DomainType &domain = levelSets.back()->getDomain();

    newDomain.initialize(domain.getNewSegmentation(),
                         domain.getAllocation() *
                             (2.0 / levelSets.back()->getLevelSetWidth()));

    // save how data should be transferred to new level set
    // list of indices into the old pointData vector
    std::vector<std::vector<unsigned>> newDataSourceIds;
    newDataSourceIds.resize(newDomain.getNumberOfSegments());

#ifdef DEBUG_LS_ADVECT_HPP
    {
      auto mesh = lsSmartPointer<lsMesh<T>>::New();
      lsToMesh<T, D>(levelSets.back(), mesh).apply();
      lsVTKWriter<T>(mesh, "lsAdvect_beforeRebuild.vtk").apply();
    }
#endif

#pragma omp parallel num_threads(newDomain.getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      auto &domainSegment = newDomain.getDomainSegment(p);

      hrleVectorType<hrleIndexType, D> startVector =
          (p == 0) ? grid.getMinGridPoint()
                   : newDomain.getSegmentation()[p - 1];

      hrleVectorType<hrleIndexType, D> endVector =
          (p != static_cast<int>(newDomain.getNumberOfSegments() - 1))
              ? newDomain.getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      // reserve a bit more to avoid reallocation
      // would expect number of points to roughly double
      newDataSourceIds[p].reserve(2.5 * domainSegment.getNumberOfPoints());

      for (hrleSparseStarIterator<typename lsDomain<T, D>::DomainType> it(
               domain, startVector);
           it.getIndices() < endVector; ++it) {

        // if the center is an active grid point
        // <1.0 since it could have been change by 0.5 max
        if (std::abs(it.getCenter().getValue()) <= 1.0) {

          int k = 0;
          for (; k < 2 * D; k++)
            if (std::signbit(it.getNeighbor(k).getValue() - 1e-7) !=
                std::signbit(it.getCenter().getValue() + 1e-7))
              break;

          // if there is at least one neighbor of opposite sign
          if (k != 2 * D) {
            if (it.getCenter().getDefinedValue() > 0.5) {
              int j = 0;
              for (; j < 2 * D; j++) {
                if (std::abs(it.getNeighbor(j).getValue()) <= 1.0)
                  if (it.getNeighbor(j).getDefinedValue() < -0.5)
                    break;
              }
              if (j == 2 * D) {
                domainSegment.insertNextDefinedPoint(
                    it.getIndices(), it.getCenter().getDefinedValue());
                newDataSourceIds[p].push_back(it.getCenter().getPointId());
                // if there is at least one active grid point, which is < -0.5
              } else {
                domainSegment.insertNextDefinedPoint(it.getIndices(), 0.5);
                newDataSourceIds[p].push_back(it.getNeighbor(j).getPointId());
              }
            } else if (it.getCenter().getDefinedValue() < -0.5) {
              int j = 0;
              for (; j < 2 * D; j++) {
                if (std::abs(it.getNeighbor(j).getValue()) <= 1.0)
                  if (it.getNeighbor(j).getDefinedValue() > 0.5)
                    break;
              }

              if (j == 2 * D) {
                domainSegment.insertNextDefinedPoint(
                    it.getIndices(), it.getCenter().getDefinedValue());
                newDataSourceIds[p].push_back(it.getCenter().getPointId());
                // if there is at least one active grid point, which is > 0.5
              } else {
                domainSegment.insertNextDefinedPoint(it.getIndices(), -0.5);
                newDataSourceIds[p].push_back(it.getNeighbor(j).getPointId());
              }
            } else {
              domainSegment.insertNextDefinedPoint(
                  it.getIndices(), it.getCenter().getDefinedValue());
              newDataSourceIds[p].push_back(it.getCenter().getPointId());
            }
          } else {
            domainSegment.insertNextUndefinedPoint(
                it.getIndices(), (it.getCenter().getDefinedValue() < 0)
                                     ? lsDomain<T, D>::NEG_VALUE
                                     : lsDomain<T, D>::POS_VALUE);
          }

        } else { // if the center is not an active grid point
          if (it.getCenter().getValue() >= 0) {
            int usedNeighbor = -1;
            T distance = lsDomain<T, D>::POS_VALUE;
            for (int i = 0; i < 2 * D; i++) {
              T value = it.getNeighbor(i).getValue();
              if (std::abs(value) <= 1.0 && (value < 0.)) {
                if (distance > value + 1.0) {
                  distance = value + 1.0;
                  usedNeighbor = i;
                }
              }
            }

            if (distance <= 1.) {
              domainSegment.insertNextDefinedPoint(it.getIndices(), distance);
              newDataSourceIds[p].push_back(
                  it.getNeighbor(usedNeighbor).getPointId());
            } else {
              domainSegment.insertNextUndefinedPoint(it.getIndices(),
                                                     lsDomain<T, D>::POS_VALUE);
            }

          } else {
            int usedNeighbor = -1;
            T distance = lsDomain<T, D>::NEG_VALUE;
            for (int i = 0; i < 2 * D; i++) {
              T value = it.getNeighbor(i).getValue();
              if (std::abs(value) <= 1.0 && (value > 0)) {
                if (distance < value - 1.0) {
                  // distance = std::max(distance, value - T(1.0));
                  distance = value - 1.0;
                  usedNeighbor = i;
                }
              }
            }

            if (distance >= -1.) {
              domainSegment.insertNextDefinedPoint(it.getIndices(), distance);
              newDataSourceIds[p].push_back(
                  it.getNeighbor(usedNeighbor).getPointId());
            } else {
              domainSegment.insertNextUndefinedPoint(it.getIndices(),
                                                     lsDomain<T, D>::NEG_VALUE);
            }
          }
        }
      }
    }

    // now copy old data into new level set
    auto &pointData = levelSets.back()->getPointData();
    if (!pointData.getScalarDataSize() || !pointData.getVectorDataSize()) {
      auto &newPointData = newlsDomain->getPointData();
      // concatenate all source ids into one vector
      newDataSourceIds.reserve(newlsDomain->getNumberOfPoints());
      for (unsigned i = 1; i < newDataSourceIds.size(); ++i) {
        for (unsigned j = 0; j < newDataSourceIds[i].size(); ++j) {
          newDataSourceIds[0].push_back(newDataSourceIds[i][j]);
        }
        // free memory of copied vectors
        std::vector<unsigned>().swap(newDataSourceIds[i]);
      }

      // scalars
      for (unsigned scalarId = 0; scalarId < pointData.getScalarDataSize();
           ++scalarId) {
        newPointData.insertNextScalarData(
            typename lsPointData<T>::ScalarDataType(),
            pointData.getScalarDataLabel(scalarId));
        auto &newScalars = *(newPointData.getScalarData(scalarId));
        auto &scalars = *(pointData.getScalarData(scalarId));
        newScalars.reserve(newlsDomain->getNumberOfPoints());

        for (unsigned i = 0; i < newDataSourceIds[0].size(); ++i) {
          newScalars.push_back(scalars[newDataSourceIds[0][i]]);
        }
      }

      // vectors
      for (unsigned vectorId = 0; vectorId < pointData.getVectorDataSize();
           ++vectorId) {
        newPointData.insertNextVectorData(
            typename lsPointData<T>::VectorDataType(),
            pointData.getVectorDataLabel(vectorId));
        auto &newVectors = *(newPointData.getVectorData(vectorId));
        auto &vectors = *(pointData.getVectorData(vectorId));
        newVectors.reserve(newlsDomain->getNumberOfPoints());

        for (unsigned i = 0; i < newDataSourceIds[0].size(); ++i) {
          newVectors.push_back(vectors[newDataSourceIds[0][i]]);
        }
      }
    }

    newDomain.finalize();
    newDomain.segment();
    levelSets.back()->deepCopy(newlsDomain);
    levelSets.back()->finalize(2);
  }

  /// internal function used as a wrapper to call specialized integrateTime
  /// with the chosen integrationScheme
  double advect(double maxTimeStep = std::numeric_limits<double>::max()) {
    // check whether a level set and velocites have been given
    if (levelSets.size() < 1) {
      lsMessage::getInstance()
          .addWarning("No level sets passed to lsAdvect. Not advecting.")
          .print();
      return std::numeric_limits<double>::max();
    }
    if (velocities == nullptr) {
      lsMessage::getInstance()
          .addWarning("No velocity field passed to lsAdvect. Not advecting.")
          .print();
      return std::numeric_limits<double>::max();
    }

    // clear all metadata as it will be invalidated
    levelSets.back()->clearMetaData();

    double currentTime = 0.;
    if (integrationScheme ==
        lsIntegrationSchemeEnum::ENGQUIST_OSHER_1ST_ORDER) {
      lsInternal::lsEnquistOsher<T, D, 1>::prepareLS(levelSets.back());
      auto is = lsInternal::lsEnquistOsher<T, D, 1>(
          levelSets.back(), velocities, calculateNormalVectors);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               lsIntegrationSchemeEnum::ENGQUIST_OSHER_2ND_ORDER) {
      lsInternal::lsEnquistOsher<T, D, 2>::prepareLS(levelSets.back());
      auto is = lsInternal::lsEnquistOsher<T, D, 2>(
          levelSets.back(), velocities, calculateNormalVectors);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               lsIntegrationSchemeEnum::LAX_FRIEDRICHS_1ST_ORDER) {
      lsInternal::lsLaxFriedrichs<T, D, 1>::prepareLS(levelSets.back());
      auto is = lsInternal::lsLaxFriedrichs<T, D, 1>(
          levelSets.back(), velocities, dissipationAlpha,
          calculateNormalVectors);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               lsIntegrationSchemeEnum::LAX_FRIEDRICHS_2ND_ORDER) {
      lsInternal::lsLaxFriedrichs<T, D, 2>::prepareLS(levelSets.back());
      auto is = lsInternal::lsLaxFriedrichs<T, D, 2>(
          levelSets.back(), velocities, dissipationAlpha,
          calculateNormalVectors);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               lsIntegrationSchemeEnum::
                   LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER) {
      lsInternal::lsLocalLaxFriedrichsAnalytical<T, D, 1>::prepareLS(
          levelSets.back());
      auto is = lsInternal::lsLocalLaxFriedrichsAnalytical<T, D, 1>(
          levelSets.back(), velocities);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               lsIntegrationSchemeEnum::LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER) {
      lsInternal::lsLocalLocalLaxFriedrichs<T, D, 1>::prepareLS(
          levelSets.back());
      auto is = lsInternal::lsLocalLocalLaxFriedrichs<T, D, 1>(
          levelSets.back(), velocities, dissipationAlpha);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               lsIntegrationSchemeEnum::LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER) {
      lsInternal::lsLocalLocalLaxFriedrichs<T, D, 2>::prepareLS(
          levelSets.back());
      auto is = lsInternal::lsLocalLocalLaxFriedrichs<T, D, 2>(
          levelSets.back(), velocities, dissipationAlpha);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               lsIntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_1ST_ORDER) {
      lsInternal::lsLocalLaxFriedrichs<T, D, 1>::prepareLS(levelSets.back());
      auto is = lsInternal::lsLocalLaxFriedrichs<T, D, 1>(
          levelSets.back(), velocities, dissipationAlpha);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               lsIntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_2ND_ORDER) {
      lsInternal::lsLocalLaxFriedrichs<T, D, 2>::prepareLS(levelSets.back());
      auto is = lsInternal::lsLocalLaxFriedrichs<T, D, 2>(
          levelSets.back(), velocities, dissipationAlpha);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               lsIntegrationSchemeEnum::
                   STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER) {
      lsInternal::lsStencilLocalLaxFriedrichsScalar<T, D, 1>::prepareLS(
          levelSets.back());
      auto is = lsInternal::lsStencilLocalLaxFriedrichsScalar<T, D, 1>(
          levelSets.back(), velocities, dissipationAlpha);
      currentTime = integrateTime(is, maxTimeStep);
    } else {
      lsMessage::getInstance()
          .addWarning("lsAdvect: Integration scheme not found. Not advecting.")
          .print();

      // if no correct scheme was found return infinity
      // to stop advection
      return std::numeric_limits<double>::max();
    }

    rebuildLS();

    // Adjust all level sets below the advected one
    // This means, that when the top levelset and one below
    // are etched, the lower one is moved with the top levelset
    // TODO: Adjust lower layers also when they have grown,
    // to allow for two different growth rates of materials
    if (integrationScheme !=
        lsIntegrationSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER) {
      for (unsigned i = 0; i < levelSets.size() - 1; ++i) {
        lsBooleanOperation<T, D>(levelSets[i], levelSets.back(),
                                 lsBooleanOperationEnum::INTERSECT)
            .apply();
      }
    }

    return currentTime;
  }

  /// Internal function used to calculate the deltas to be applied to the LS
  /// values from the given velocities and the integration scheme to be used.
  /// Level Sets below are also considered in order to adjust the advection
  /// depth accordingly if there would be a material change.
  template <class IntegrationSchemeType>
  double
  integrateTime(IntegrationSchemeType IntegrationScheme,
                double maxTimeStep = std::numeric_limits<double>::max()) {
    if (timeStepRatio >= 0.5) {
      lsMessage::getInstance()
          .addWarning("Integration time step ratio should be smaller than 0.5. "
                      "Advection might fail!")
          .print();
    }

    auto &topDomain = levelSets.back()->getDomain();
    auto &grid = levelSets.back()->getGrid();

    std::vector<std::vector<std::pair<T, T>>> totalTempRates;
    totalTempRates.resize((levelSets.back())->getNumberOfSegments());

    if (ignoreVoids) {
      lsMarkVoidPoints<T, D>(levelSets.back()).apply();
    }

#pragma omp parallel num_threads((levelSets.back())->getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      hrleVectorType<hrleIndexType, D> startVector =
          (p == 0) ? grid.getMinGridPoint()
                   : topDomain.getSegmentation()[p - 1];

      hrleVectorType<hrleIndexType, D> endVector =
          (p != static_cast<int>(topDomain.getNumberOfSegments() - 1))
              ? topDomain.getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      double tempMaxTimeStep = maxTimeStep;
      auto &tempRates = totalTempRates[p];
      tempRates.reserve(topDomain.getNumberOfPoints() /
                            double((levelSets.back())->getNumberOfSegments()) +
                        10);

      // an iterator for each level set
      std::vector<hrleSparseIterator<typename lsDomain<T, D>::DomainType>>
          iterators;
      for (auto it = levelSets.begin(); it != levelSets.end(); ++it) {
        iterators.push_back(
            hrleSparseIterator<typename lsDomain<T, D>::DomainType>(
                (*it)->getDomain()));
      }

      IntegrationSchemeType scheme(IntegrationScheme);
      auto &voidPoints =
          *(levelSets.back()->getPointData().getScalarData("VoidPointMarkers"));

      for (hrleSparseIterator<typename lsDomain<T, D>::DomainType> it(
               topDomain, startVector);
           it.getStartIndices() < endVector; ++it) {

        if (!it.isDefined() || std::abs(it.getValue()) > 0.5)
          continue;

        T value = it.getValue();
        double maxStepTime = 0;
        double cfl = timeStepRatio;

        for (int currentLevelSetId = levelSets.size() - 1;
             currentLevelSetId >= 0; --currentLevelSetId) {

          T velocity = 0;

          if (!(ignoreVoids && voidPoints[it.getPointId()])) {
            // check if there is any other levelset at the same point:
            // if yes, take the velocity of the lowest levelset
            for (unsigned lowerLevelSetId = 0;
                 lowerLevelSetId < levelSets.size(); ++lowerLevelSetId) {
              // put iterator to same position as the top levelset
              iterators[lowerLevelSetId].goToIndicesSequential(
                  it.getStartIndices());

              // if the lower surface is actually outside, i.e. its LS value is
              // lower or equal
              if (iterators[lowerLevelSetId].getValue() <=
                  value + wrappingLayerEpsilon) {
                velocity = scheme(it.getStartIndices(), lowerLevelSetId);
                break;
              }
            }
          }

          T valueBelow;
          // get value of material below (before in levelSets list)
          if (currentLevelSetId > 0) {
            iterators[currentLevelSetId - 1].goToIndicesSequential(
                it.getStartIndices());
            valueBelow = iterators[currentLevelSetId - 1].getValue();
          } else {
            valueBelow = std::numeric_limits<T>::max();
          }

          // if velocity is positive, set maximum time step possible without
          // violating the cfl condition
          if (velocity > 0.) {
            maxStepTime += cfl / velocity;
            tempRates.push_back(
                std::make_pair(velocity, -std::numeric_limits<T>::max()));
            break;
            // if velocity is 0, maximum time step is infinite
          } else if (velocity == 0.) {
            maxStepTime = std::numeric_limits<T>::max();
            tempRates.push_back(
                std::make_pair(velocity, std::numeric_limits<T>::max()));
            break;
            // if the velocity is negative apply the velocity for as long as
            // possible without infringing on material below
          } else {
            T difference = std::abs(valueBelow - value);

            if (difference >= cfl) {
              maxStepTime -= cfl / velocity;
              tempRates.push_back(
                  std::make_pair(velocity, std::numeric_limits<T>::max()));
              break;
            } else {
              maxStepTime -= difference / velocity;
              // the second part of the pair indicates how far we can move in
              // this time step until the end of the material is reached
              tempRates.push_back(std::make_pair(velocity, valueBelow));
              cfl -= difference;
              // use new LS value for next calculations
              value = valueBelow;
            }
          }
        }

        if (maxStepTime < tempMaxTimeStep)
          tempMaxTimeStep = maxStepTime;
      }

#pragma omp critical
      {
        // If scheme is STENCIL_LOCAL_LAX_FRIEDRICHS the time step is reduced
        // depending on the dissipation coefficients For all remaining schemes
        // this function is empty.
        reduceTimeStepHamiltonJacobi(scheme, tempMaxTimeStep);

        // set global timestep maximum
        if (tempMaxTimeStep < maxTimeStep)
          maxTimeStep = tempMaxTimeStep;
      }
    }

    // reduce to one layer thickness and apply new values directly to the
    // domain segments --> DO NOT CHANGE SEGMENTATION HERE (true parameter)
    lsReduce<T, D>(levelSets.back(), 1, true).apply();

    const bool saveVelocities = saveAdvectionVelocities;
    std::vector<std::vector<double>> velocityVectors(
        levelSets.back()->getNumberOfSegments());

#pragma omp parallel num_threads((levelSets.back())->getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      typename std::vector<std::pair<T, T>>::const_iterator itRS =
          totalTempRates[p].begin();
      auto &segment = topDomain.getDomainSegment(p);
      unsigned maxId = segment.getNumberOfPoints();

      if (saveVelocities) {
        velocityVectors[p].resize(maxId);
      }

      for (unsigned localId = 0; localId < maxId; ++localId) {
        T &value = segment.definedValues[localId];
        double time = maxTimeStep;

        // if there is a change in materials during one time step, deduct the
        // time taken to advect up to the end of the top material and set the LS
        // value to the one below
        while (std::abs(itRS->second - value) < std::abs(time * itRS->first)) {
          time -= std::abs((itRS->second - value) / itRS->first);
          value = itRS->second;
          ++itRS;
        }

        // now deduct the velocity times the time step we take
        value -= time * itRS->first;

        if (saveVelocities) {
          velocityVectors[p][localId] = time * itRS->first;
        }

        // this
        // is run when two materials are close but the velocity is too slow to
        // actually reach the second material, to get rid of the extra entry in
        // the TempRatesStop
        while (std::abs(itRS->second) != std::numeric_limits<T>::max())
          ++itRS;

        // advance the TempStopRates iterator by one
        ++itRS;
      }
    }

    if (saveVelocities) {
      typename lsPointData<T>::ScalarDataType pointData;
      for (unsigned i = 0; i < velocityVectors.size(); ++i) {
        pointData.insert(pointData.end(),
                         std::make_move_iterator(velocityVectors[i].begin()),
                         std::make_move_iterator(velocityVectors[i].end()));
      }
      levelSets.back()->getPointData().insertNextScalarData(
          pointData, "AdvectionVelocity");
    }

    return maxTimeStep;
  }

public:
  lsAdvect() {}

  lsAdvect(lsSmartPointer<lsDomain<T, D>> passedlsDomain) {
    levelSets.push_back(passedlsDomain);
  }

  template <class VelocityField>
  lsAdvect(lsSmartPointer<lsDomain<T, D>> passedlsDomain,
           lsSmartPointer<VelocityField> passedVelocities) {
    levelSets.push_back(passedlsDomain);
    velocities =
        std::dynamic_pointer_cast<lsVelocityField<T>>(passedVelocities);
  }

  template <class VelocityField,
            lsConcepts::IsBaseOf<lsVelocityField<T>, VelocityField> =
                lsConcepts::assignable>
  lsAdvect(lsSmartPointer<VelocityField> passedVelocities) {
    velocities =
        std::dynamic_pointer_cast<lsVelocityField<T>>(passedVelocities);
  }

  template <class VelocityField,
            lsConcepts::IsBaseOf<lsVelocityField<T>, VelocityField> =
                lsConcepts::assignable>
  lsAdvect(std::vector<lsSmartPointer<lsDomain<T, D>>> passedlsDomains,
           lsSmartPointer<VelocityField> passedVelocities)
      : levelSets(passedlsDomains) {
    velocities =
        std::dynamic_pointer_cast<lsVelocityField<T>>(passedVelocities);
  }

  /// Pushes the passed level set to the back of the list of level sets
  /// used for advection.
  void insertNextLevelSet(lsSmartPointer<lsDomain<T, D>> passedlsDomain) {
    levelSets.push_back(passedlsDomain);
  }

  /// Set the velocity field used for advection. This should be a concrete
  /// implementation of lsVelocityField
  template <class VelocityField,
            lsConcepts::IsBaseOf<lsVelocityField<T>, VelocityField> =
                lsConcepts::assignable>
  void setVelocityField(lsSmartPointer<VelocityField> passedVelocities) {
    velocities =
        std::dynamic_pointer_cast<lsVelocityField<T>>(passedVelocities);
  }

  /// Set the time until when the level set should be advected.
  /// If this takes more than one advection step, multiple will
  /// be performed. Defaults to 0, which means one advection
  /// step with the maximum time step possible according to the
  /// CFL condition(see setTimeStepRatio) will be performed.
  void setAdvectionTime(double time) { advectionTime = time; }

  /// Set the CFL condition to use during advection.
  /// The CFL condition sets the maximum distance a surface can
  /// be moved during one advection step. It MUST be below 0.5
  /// to guarantee numerical stability. Defaults to 0.4999.
  void setTimeStepRatio(const double &cfl) { timeStepRatio = cfl; }

  /// Set whether normal vectors should be calculated at each level
  /// set point. Defaults to true. If normal vectors are not required
  /// for velocity calculation, this can be set to false, in order
  /// to increase computational efficiency.
  void setCalculateNormalVectors(bool cnv) { calculateNormalVectors = cnv; }

  /// Set whether level set values, which are not part of the "top"
  /// geometrically connected part of values, should be advected.
  /// The "top" part is identified by the most positive part in the
  /// lowest dimension with INFINITE boundary conditions.
  /// Defaults to false. If set to true, only the "top" values will
  /// be advected. All others values are not changed.
  void setIgnoreVoids(bool iV) { ignoreVoids = iV; }

  /// Set whether the velocities applied to each point should be saved in
  /// the level set for debug purposes.
  void setSaveAdvectionVelocities(bool sAV) { saveAdvectionVelocities = sAV; }

  /// Get by how much the physical time was advanced during the last apply()
  /// call.
  double getAdvectedTime() { return advectedTime; }

  /// Get how many advection steps were performed during the last apply() call.
  unsigned getNumberOfTimeSteps() { return numberOfTimeSteps; }

  /// Get the value of the CFL number.
  double getTimeStepRatio() { return timeStepRatio; }

  /// Get whether normal vectors were caluclated.
  bool getCalculateNormalVectors() { return calculateNormalVectors; }

  /// Set which integration scheme should be used out of the ones specified
  /// in lsIntegrationSchemeEnum.
  void setIntegrationScheme(lsIntegrationSchemeEnum scheme) {
    integrationScheme = scheme;
  }

  /// Set the alpha dissipation coefficient.
  /// For lsLaxFriedrichs, this is used as the alpha value.
  /// For all other LaxFriedrichs schemes it is used as a
  /// scaling factor for the calculated alpha values.
  void setDissipationAlpha(const double &a) { dissipationAlpha = a; }

  /// Perform the advection.
  void apply() {
    if (advectionTime == 0.) {
      advectedTime = advect();
      numberOfTimeSteps = 1;
    } else {
      double currentTime = 0.0;
      numberOfTimeSteps = 0;
      while (currentTime < advectionTime) {
        currentTime += advect(advectionTime - currentTime);
        ++numberOfTimeSteps;
      }
      advectedTime = currentTime;
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsAdvect)

#endif // LS_ADVECT_HPP
