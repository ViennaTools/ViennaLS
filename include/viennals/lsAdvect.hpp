#pragma once

#include <lsPreCompileMacros.hpp>

#include <limits>
#include <vector>

#include <hrleSparseIterator.hpp>
#include <hrleSparseStarIterator.hpp>

#include <vcLogger.hpp>
#include <vcSmartPointer.hpp>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMarkVoidPoints.hpp>
#include <lsReduce.hpp>

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

namespace viennals {

using namespace viennacore;

/// Enumeration for the different Integration schemes
/// used by the advection kernel
enum struct IntegrationSchemeEnum : unsigned {
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
template <class T, int D> class Advect {
  std::vector<SmartPointer<Domain<T, D>>> levelSets;
  SmartPointer<VelocityField<T>> velocities = nullptr;
  IntegrationSchemeEnum integrationScheme =
      IntegrationSchemeEnum::ENGQUIST_OSHER_1ST_ORDER;
  double timeStepRatio = 0.4999;
  double dissipationAlpha = 1.0;
  bool calculateNormalVectors = true;
  bool ignoreVoids = false;
  double advectionTime = 0.;
  bool performOnlySingleStep = false;
  double advectedTime = 0.;
  unsigned numberOfTimeSteps = 0;
  bool saveAdvectionVelocities = false;
  bool updatePointData = true;
  static constexpr double wrappingLayerEpsilon = 1e-4;

  void rebuildLS() {
    // TODO: this function uses manhatten distances for renormalisation,
    // since this is the quickest. For visualisation applications, better
    // renormalisation is needed, so it might be good to implement
    // Euler distance renormalisation as an option
    auto &grid = levelSets.back()->getGrid();
    auto newlsDomain = SmartPointer<Domain<T, D>>::New(grid);
    typename Domain<T, D>::DomainType &newDomain = newlsDomain->getDomain();
    typename Domain<T, D>::DomainType &domain = levelSets.back()->getDomain();

    newDomain.initialize(domain.getNewSegmentation(),
                         domain.getAllocation() *
                             (2.0 / levelSets.back()->getLevelSetWidth()));

    const bool updateData = updatePointData;
    // save how data should be transferred to new level set
    // list of indices into the old pointData vector
    std::vector<std::vector<unsigned>> newDataSourceIds;
    if (updateData)
      newDataSourceIds.resize(newDomain.getNumberOfSegments());

#ifdef DEBUG_LS_ADVECT_HPP
    {
      auto mesh = SmartPointer<Mesh<T>>::New();
      ToMesh<T, D>(levelSets.back(), mesh).apply();
      VTKWriter<T>(mesh, "Advect_beforeRebuild.vtk").apply();
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
      if (updateData)
        newDataSourceIds[p].reserve(2.5 * domainSegment.getNumberOfPoints());

      for (hrleSparseStarIterator<typename Domain<T, D>::DomainType, 1> it(
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
                if (updateData)
                  newDataSourceIds[p].push_back(it.getCenter().getPointId());
                // if there is at least one active grid point, which is < -0.5
              } else {
                domainSegment.insertNextDefinedPoint(it.getIndices(), 0.5);
                if (updateData)
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
                if (updateData)
                  newDataSourceIds[p].push_back(it.getCenter().getPointId());
                // if there is at least one active grid point, which is > 0.5
              } else {
                domainSegment.insertNextDefinedPoint(it.getIndices(), -0.5);
                if (updateData)
                  newDataSourceIds[p].push_back(it.getNeighbor(j).getPointId());
              }
            } else {
              domainSegment.insertNextDefinedPoint(
                  it.getIndices(), it.getCenter().getDefinedValue());
              if (updateData)
                newDataSourceIds[p].push_back(it.getCenter().getPointId());
            }
          } else {
            domainSegment.insertNextUndefinedPoint(
                it.getIndices(), (it.getCenter().getDefinedValue() < 0)
                                     ? Domain<T, D>::NEG_VALUE
                                     : Domain<T, D>::POS_VALUE);
          }

        } else { // if the center is not an active grid point
          if (it.getCenter().getValue() >= 0) {
            int usedNeighbor = -1;
            T distance = Domain<T, D>::POS_VALUE;
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
              if (updateData)
                newDataSourceIds[p].push_back(
                    it.getNeighbor(usedNeighbor).getPointId());
            } else {
              domainSegment.insertNextUndefinedPoint(it.getIndices(),
                                                     Domain<T, D>::POS_VALUE);
            }

          } else {
            int usedNeighbor = -1;
            T distance = Domain<T, D>::NEG_VALUE;
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
              if (updateData)
                newDataSourceIds[p].push_back(
                    it.getNeighbor(usedNeighbor).getPointId());
            } else {
              domainSegment.insertNextUndefinedPoint(it.getIndices(),
                                                     Domain<T, D>::NEG_VALUE);
            }
          }
        }
      }
    }

    // now copy old data into new level set
    if (updateData) {
      auto &pointData = levelSets.back()->getPointData();
      newlsDomain->getPointData().translateFromMultiData(pointData,
                                                         newDataSourceIds);
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
      Logger::getInstance()
          .addWarning("No level sets passed to Advect. Not advecting.")
          .print();
      return std::numeric_limits<double>::max();
    }
    if (velocities == nullptr) {
      Logger::getInstance()
          .addWarning("No velocity field passed to Advect. Not advecting.")
          .print();
      return std::numeric_limits<double>::max();
    }

    prepareLS();

    double currentTime = 0.;
    if (integrationScheme == IntegrationSchemeEnum::ENGQUIST_OSHER_1ST_ORDER) {
      auto is = lsInternal::EnquistOsher<T, D, 1>(levelSets.back(), velocities,
                                                  calculateNormalVectors);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               IntegrationSchemeEnum::ENGQUIST_OSHER_2ND_ORDER) {
      auto is = lsInternal::EnquistOsher<T, D, 2>(levelSets.back(), velocities,
                                                  calculateNormalVectors);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               IntegrationSchemeEnum::LAX_FRIEDRICHS_1ST_ORDER) {
      auto is = lsInternal::LaxFriedrichs<T, D, 1>(levelSets.back(), velocities,
                                                   dissipationAlpha,
                                                   calculateNormalVectors);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               IntegrationSchemeEnum::LAX_FRIEDRICHS_2ND_ORDER) {
      auto is = lsInternal::LaxFriedrichs<T, D, 2>(levelSets.back(), velocities,
                                                   dissipationAlpha,
                                                   calculateNormalVectors);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               IntegrationSchemeEnum::
                   LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER) {
      auto is = lsInternal::LocalLaxFriedrichsAnalytical<T, D, 1>(
          levelSets.back(), velocities);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               IntegrationSchemeEnum::LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER) {
      auto is = lsInternal::LocalLocalLaxFriedrichs<T, D, 1>(
          levelSets.back(), velocities, dissipationAlpha);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               IntegrationSchemeEnum::LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER) {
      auto is = lsInternal::LocalLocalLaxFriedrichs<T, D, 2>(
          levelSets.back(), velocities, dissipationAlpha);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               IntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_1ST_ORDER) {
      auto is = lsInternal::LocalLaxFriedrichs<T, D, 1>(
          levelSets.back(), velocities, dissipationAlpha);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               IntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_2ND_ORDER) {
      auto is = lsInternal::LocalLaxFriedrichs<T, D, 2>(
          levelSets.back(), velocities, dissipationAlpha);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme ==
               IntegrationSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER) {
      auto is = lsInternal::StencilLocalLaxFriedrichsScalar<T, D, 1>(
          levelSets.back(), velocities, dissipationAlpha);
      currentTime = integrateTime(is, maxTimeStep);
    } else {
      Logger::getInstance()
          .addWarning("Advect: Integration scheme not found. Not advecting.")
          .print();

      // if no correct scheme was found return infinity
      // to stop advection
      return std::numeric_limits<double>::max();
    }

    rebuildLS();

    // Adjust all level sets below the advected one
    // This means, that when the top level-set and one below
    // are etched, the lower one is moved with the top level-set
    // TODO: Adjust lower layers also when they have grown,
    // to allow for two different growth rates of materials
    if (integrationScheme !=
        IntegrationSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER) {
      for (unsigned i = 0; i < levelSets.size() - 1; ++i) {
        BooleanOperation<T, D>(levelSets[i], levelSets.back(),
                               BooleanOperationEnum::INTERSECT)
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
      Logger::getInstance()
          .addWarning("Integration time step ratio should be smaller than 0.5. "
                      "Advection might fail!")
          .print();
    }

    auto &topDomain = levelSets.back()->getDomain();
    auto &grid = levelSets.back()->getGrid();

    std::vector<std::vector<std::pair<T, T>>> totalTempRates;
    totalTempRates.resize((levelSets.back())->getNumberOfSegments());

    typename PointData<T>::ScalarDataType *voidMarkerPointer;
    if (ignoreVoids) {
      MarkVoidPoints<T, D>(levelSets.back()).apply();
      auto &pointData = levelSets.back()->getPointData();
      voidMarkerPointer =
          pointData.getScalarData(MarkVoidPoints<T, D>::voidPointLabel, true);
      if (voidMarkerPointer == nullptr) {
        Logger::getInstance()
            .addWarning("Advect: Cannot find void point markers. Not "
                        "ignoring void points.")
            .print();
        ignoreVoids = false;
      }
    }

    // This function is only defined for the Lax-Friedrichs scheme.
    // The global alpha values are stored in the integration scheme.
    if (integrationScheme == IntegrationSchemeEnum::LAX_FRIEDRICHS_1ST_ORDER) {
      lsInternal::advect::findGlobalAlpha<IntegrationSchemeType, T, D, 1>(
          IntegrationScheme, levelSets, velocities);
    } else if (integrationScheme ==
               IntegrationSchemeEnum::LAX_FRIEDRICHS_2ND_ORDER) {
      lsInternal::advect::findGlobalAlpha<IntegrationSchemeType, T, D, 2>(
          IntegrationScheme, levelSets, velocities);
    }

    const bool ignoreVoidPoints = ignoreVoids;

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
      std::vector<hrleSparseIterator<typename Domain<T, D>::DomainType>>
          iterators;
      for (auto it = levelSets.begin(); it != levelSets.end(); ++it) {
        iterators.push_back(
            hrleSparseIterator<typename Domain<T, D>::DomainType>(
                (*it)->getDomain()));
      }

      IntegrationSchemeType scheme(IntegrationScheme);

      for (hrleSparseIterator<typename Domain<T, D>::DomainType> it(
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

          if (!(ignoreVoidPoints && (*voidMarkerPointer)[it.getPointId()])) {
            // check if there is any other levelset at the same point:
            // if yes, take the velocity of the lowest levelset
            for (unsigned lowerLevelSetId = 0;
                 lowerLevelSetId < levelSets.size(); ++lowerLevelSetId) {
              // put iterator to same position as the top levelset
              iterators[lowerLevelSetId].goToIndicesSequential(
                  it.getStartIndices());

              // if the lower surface is actually outside, i.e. its LS value
              // is lower or equal
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
              // the second part of the pair indicates how far we can move
              // in this time step until the end of the material is reached
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
        // If a Lax Friedrichs scheme is selected the time step is
        // reduced depending on the dissipation coefficients
        // For Engquist Osher scheme this function is empty.
        scheme.reduceTimeStepHamiltonJacobi(
            tempMaxTimeStep, levelSets.back()->getGrid().getGridDelta());

        // set global timestep maximum
        if (tempMaxTimeStep < maxTimeStep)
          maxTimeStep = tempMaxTimeStep;
      }
    }

    // reduce to one layer thickness and apply new values directly to the
    // domain segments --> DO NOT CHANGE SEGMENTATION HERE (true parameter)
    Reduce<T, D>(levelSets.back(), 1, true).apply();

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

        // if there is a change in materials during one time step, deduct
        // the time taken to advect up to the end of the top material and
        // set the LS value to the one below
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
        // is run when two materials are close but the velocity is too slow
        // to actually reach the second material, to get rid of the extra
        // entry in the TempRatesStop
        while (std::abs(itRS->second) != std::numeric_limits<T>::max())
          ++itRS;

        // advance the TempStopRates iterator by one
        ++itRS;
      }
    }

    if (saveVelocities) {
      auto &pointData = levelSets.back()->getPointData();
      // delete if already exists
      if (int i = pointData.getScalarDataIndex(velocityLabel); i != -1) {
        pointData.eraseScalarData(i);
      }

      typename PointData<T>::ScalarDataType vels;

      for (unsigned i = 0; i < velocityVectors.size(); ++i) {
        vels.insert(vels.end(),
                    std::make_move_iterator(velocityVectors[i].begin()),
                    std::make_move_iterator(velocityVectors[i].end()));
      }
      pointData.insertNextScalarData(std::move(vels), velocityLabel);
    }

    return maxTimeStep;
  }

public:
  static constexpr char velocityLabel[] = "AdvectionVelocities";

  Advect() {}

  Advect(SmartPointer<Domain<T, D>> passedlsDomain) {
    levelSets.push_back(passedlsDomain);
  }

  Advect(SmartPointer<Domain<T, D>> passedlsDomain,
         SmartPointer<VelocityField<T>> passedVelocities) {
    levelSets.push_back(passedlsDomain);
    velocities = passedVelocities;
  }

  Advect(std::vector<SmartPointer<Domain<T, D>>> passedlsDomains,
         SmartPointer<VelocityField<T>> passedVelocities)
      : levelSets(passedlsDomains) {
    velocities = passedVelocities;
  }

  /// Pushes the passed level set to the back of the list of level sets
  /// used for advection.
  void insertNextLevelSet(SmartPointer<Domain<T, D>> passedlsDomain) {
    levelSets.push_back(passedlsDomain);
  }

  /// Set the velocity field used for advection. This should be a concrete
  /// implementation of lsVelocityField
  void setVelocityField(SmartPointer<VelocityField<T>> passedVelocities) {
    velocities = passedVelocities;
  }

  /// Set the time until when the level set should be advected.
  /// If this takes more than one advection step, multiple will
  /// be performed. Defaults to 0, which means one advection
  /// step with the maximum time step possible according to the
  /// CFL condition(see setTimeStepRatio) will be performed.
  void setAdvectionTime(double time) { advectionTime = time; }

  /// If set to true, only a single advection step will be
  /// performed, even if the advection time set with
  /// setAdvectionTime(double) would require several steps to pass.
  /// Defaults to false.
  void setSingleStep(bool singleStep) { performOnlySingleStep = singleStep; }

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
  /// in IntegrationSchemeEnum.
  void setIntegrationScheme(IntegrationSchemeEnum scheme) {
    integrationScheme = scheme;
  }

  /// Set the alpha dissipation coefficient.
  /// For lsLaxFriedrichs, this is used as the alpha value.
  /// For all other LaxFriedrichs schemes it is used as a
  /// scaling factor for the calculated alpha values.
  void setDissipationAlpha(const double &a) { dissipationAlpha = a; }

  /// Set whether the point data in the old LS should
  /// be translated to the advected LS. Defaults to true.
  void setUpdatePointData(bool update) { updatePointData = update; }

  // Prepare the levelset for advection, based on the provided integration
  // scheme.
  void prepareLS() {
    // check whether a level set and velocities have been given
    if (levelSets.size() < 1) {
      Logger::getInstance()
          .addWarning("No level sets passed to Advect.")
          .print();
      return;
    }

    if (integrationScheme == IntegrationSchemeEnum::ENGQUIST_OSHER_1ST_ORDER) {
      lsInternal::EnquistOsher<T, D, 1>::prepareLS(levelSets.back());
    } else if (integrationScheme ==
               IntegrationSchemeEnum::ENGQUIST_OSHER_2ND_ORDER) {
      lsInternal::EnquistOsher<T, D, 2>::prepareLS(levelSets.back());
    } else if (integrationScheme ==
               IntegrationSchemeEnum::LAX_FRIEDRICHS_1ST_ORDER) {
      lsInternal::LaxFriedrichs<T, D, 1>::prepareLS(levelSets.back());
    } else if (integrationScheme ==
               IntegrationSchemeEnum::LAX_FRIEDRICHS_2ND_ORDER) {
      lsInternal::LaxFriedrichs<T, D, 2>::prepareLS(levelSets.back());
    } else if (integrationScheme ==
               IntegrationSchemeEnum::
                   LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER) {
      lsInternal::LocalLaxFriedrichsAnalytical<T, D, 1>::prepareLS(
          levelSets.back());
    } else if (integrationScheme ==
               IntegrationSchemeEnum::LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER) {
      lsInternal::LocalLocalLaxFriedrichs<T, D, 1>::prepareLS(levelSets.back());
    } else if (integrationScheme ==
               IntegrationSchemeEnum::LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER) {
      lsInternal::LocalLocalLaxFriedrichs<T, D, 2>::prepareLS(levelSets.back());
    } else if (integrationScheme ==
               IntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_1ST_ORDER) {
      lsInternal::LocalLaxFriedrichs<T, D, 1>::prepareLS(levelSets.back());
    } else if (integrationScheme ==
               IntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_2ND_ORDER) {
      lsInternal::LocalLaxFriedrichs<T, D, 2>::prepareLS(levelSets.back());
    } else if (integrationScheme ==
               IntegrationSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER) {
      lsInternal::StencilLocalLaxFriedrichsScalar<T, D, 1>::prepareLS(
          levelSets.back());
    } else {
      Logger::getInstance()
          .addWarning("Advect: Integration scheme not found. Not advecting.")
          .print();
    }
  }

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
        if (performOnlySingleStep)
          break;
      }
      advectedTime = currentTime;
    }
  }
};

// add all template specializations for this class
PRECOMPILE_PRECISION_DIMENSION(Advect)

} // namespace viennals
