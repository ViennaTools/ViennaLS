#ifndef LS_ADVECT_TEMPLATE_HPP
#define LS_ADVECT_TEMPLATE_HPP

#include <limits>
#include <vector>

#include <hrleCrossIterator.hpp>
#include <hrleRunsIterator.hpp>

#include <lsBooleanOperation_template.hpp>
#include <lsDomain_template.hpp>
#include <lsMessage.hpp>
#include <lsReduce_template.hpp>

// Integration schemes
#include <lsEnquistOsher_template.hpp>
#include <lsLaxFriedrichs_template.hpp>
#include <lsStencilLocalLaxFriedrichsScalar_template.hpp>

// Velocity accessor
// #include <lsMesh.hpp>
// #include <lsToMesh_template.hpp>
// #include <lsVTKWriter.hpp>
#include <lsVelocityField_template.hpp>

template <class T, int D> class lsAdvect {
public:
  enum IntegrationSchemeEnum : unsigned {
    ENGQUIST_OSHER_1ST_ORDER = 0,
    ENGQUIST_OSHER_2ND_ORDER = 1,
    LAX_FRIEDRICHS_1ST_ORDER = 2,
    LAX_FRIEDRICHS_2ND_ORDER = 3,
    STENCIL_LOCAL_LAX_FRIEDRICHS = 4
  };

private:
  std::vector<lsDomain<T, D> *> &levelSets;
  lsVelocityField<T> *velocities;
  IntegrationSchemeEnum integrationScheme = ENGQUIST_OSHER_1ST_ORDER;
  double timeStepRatio = 0.4999;
  double dissipationAlpha = 0.;
  bool calculateNormalVectors = true;

  lsAdvect();

  // SFINAE functions needed for StencilLocalLaxFriedrichs
  template <class IntegrationSchemeType,
            typename std::enable_if<!std::is_same<
                lsInternal::lsStencilLocalLaxFriedrichsScalar<T, D, 1>,
                IntegrationSchemeType>::value>::type * = nullptr>
  void reduceTimeStepHamiltonJacobi(IntegrationSchemeType &, double &) {}

  template <class IntegrationSchemeType,
            typename std::enable_if<std::is_same<
                lsInternal::lsStencilLocalLaxFriedrichsScalar<T, D, 1>,
                IntegrationSchemeType>::value>::type * = nullptr>
  void reduceTimeStepHamiltonJacobi(IntegrationSchemeType &scheme,
                                    double &MaxTimeStep) {
    // TODO Can be potentially smaller than 1 (user input???)
    const double alpha_maxCFL = 1.0;
    // second time step test, based on alphas
    hrleVectorType<T, 3> alphas = scheme.getFinalAlphas();
    hrleVectorType<T, 3> dxs = scheme.getDeltas();

    MaxTimeStep = 0;
    for (int i = 0; i < 3; ++i) {
      if (std::abs(dxs[i]) > 1e-6) {
        MaxTimeStep += alphas[i] / dxs[i];
      }
    }

    MaxTimeStep = alpha_maxCFL / MaxTimeStep;
  }

  void rebuildLS() {
    // TODO: this function uses manhatten distances for renormalisation,
    // since this is the quickest. For visualisation applications, better
    // renormalisation is needed, so it might be good to implement
    // Euler distance renormalisation as an option
    auto &grid = levelSets.back()->getGrid();
    lsDomain<T, D> newlsDomain(grid);
    typename lsDomain<T, D>::DomainType &newDomain = newlsDomain.getDomain();
    typename lsDomain<T, D>::DomainType &domain = levelSets.back()->getDomain();

    newDomain.initialize(domain.getNewSegmentation(),
                         domain.getAllocation() *
                             (2.0 / levelSets.back()->getLevelSetWidth()));

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

      for (hrleCrossIterator<typename lsDomain<T, D>::DomainType> it(
               domain, startVector);
           it.getIndices() < endVector; ++it) {

        // if the center is an active grid point
        // <1.0 since it could have been change by 0.5 max
        if (std::abs(it.getCenter().getValue()) <= 1.0) {

          int k = 0;
          for (; k < 2 * D; k++)
            if (std::signbit(it.getNeighbor(k).getValue()) !=
                std::signbit(it.getCenter().getValue()))
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
                // if there is at least one active grid point, which is < -0.5
              } else {
                domainSegment.insertNextDefinedPoint(it.getIndices(), 0.5);
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
                // if there is at least one active grid point, which is > 0.5
              } else {
                domainSegment.insertNextDefinedPoint(it.getIndices(), -0.5);
              }
            } else {
              domainSegment.insertNextDefinedPoint(
                  it.getIndices(), it.getCenter().getDefinedValue());
            }
          } else {
            domainSegment.insertNextUndefinedPoint(
                it.getIndices(), (it.getCenter().getDefinedValue() < 0)
                                     ? lsDomain<T, D>::NEG_VALUE
                                     : lsDomain<T, D>::POS_VALUE);
          }

        } else { // if the center is not an active grid point
          if (it.getCenter().getValue() >= 0) {
            T distance = lsDomain<T, D>::POS_VALUE;
            for (int i = 0; i < 2 * D; i++) {
              T value = it.getNeighbor(i).getValue();
              if (std::abs(value) <= 1.0 && (value < 0.)) {
                if (distance > value + 1.0) {
                  distance = value + 1.0;
                }
              }
            }

            if (distance <= 1.) {
              domainSegment.insertNextDefinedPoint(it.getIndices(), distance);
            } else {
              domainSegment.insertNextUndefinedPoint(it.getIndices(),
                                                     lsDomain<T, D>::POS_VALUE);
            }

          } else {
            T distance = lsDomain<T, D>::NEG_VALUE;
            for (int i = 0; i < 2 * D; i++) {
              T value = it.getNeighbor(i).getValue();
              if (std::abs(value) <= 1.0 && (value > 0)) {
                if (distance < value - 1.0) {
                  distance = std::max(distance, value - 1.0);
                }
              }
            }

            if (distance >= -1.) {
              domainSegment.insertNextDefinedPoint(it.getIndices(), distance);
            } else {
              domainSegment.insertNextUndefinedPoint(it.getIndices(),
                                                     lsDomain<T, D>::NEG_VALUE);
            }
          }
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
    double currentTime = 0.;
    if (integrationScheme == ENGQUIST_OSHER_1ST_ORDER) {
      lsInternal::lsEnquistOsher<T, D, 1>::prepareLS(*(levelSets.back()));
      auto is = lsInternal::lsEnquistOsher<T, D, 1>(*(levelSets.back()),
                                                    calculateNormalVectors);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme == ENGQUIST_OSHER_2ND_ORDER) {
      lsInternal::lsEnquistOsher<T, D, 2>::prepareLS(*(levelSets.back()));
      auto is = lsInternal::lsEnquistOsher<T, D, 2>(*(levelSets.back()),
                                                    calculateNormalVectors);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme == LAX_FRIEDRICHS_1ST_ORDER) {
      lsInternal::lsLaxFriedrichs<T, D, 1>::prepareLS(*(levelSets.back()));
      auto is = lsInternal::lsLaxFriedrichs<T, D, 1>(*(levelSets.back()),
                                                     calculateNormalVectors);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme == LAX_FRIEDRICHS_2ND_ORDER) {
      lsInternal::lsLaxFriedrichs<T, D, 2>::prepareLS(*(levelSets.back()));
      auto is = lsInternal::lsLaxFriedrichs<T, D, 2>(*(levelSets.back()),
                                                     calculateNormalVectors);
      currentTime = integrateTime(is, maxTimeStep);
    } else if (integrationScheme == STENCIL_LOCAL_LAX_FRIEDRICHS) {
      lsInternal::lsStencilLocalLaxFriedrichsScalar<T, D, 1>::prepareLS(
          *(levelSets.back()));
      auto is = lsInternal::lsStencilLocalLaxFriedrichsScalar<T, D, 1>(
          *(levelSets.back()));
      currentTime = integrateTime(is, maxTimeStep);
    } else {
      // if no correct scheme was found return infinity
      return std::numeric_limits<double>::max();
    }

    rebuildLS();

    // Adjust all level sets below the advected one
    // This means, that when the top levelset and one below
    // are etched, the lower one is moved with the top levelset
    // TODO: Adjust lower layers also when they have grown,
    // to allow for two different growth rates of materials
    if (integrationScheme != STENCIL_LOCAL_LAX_FRIEDRICHS) {
      for (unsigned i = 0; i < levelSets.size() - 1; ++i) {
        lsBooleanOperation<T, D>(*levelSets[i]).max(*(levelSets.back()));
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
                        double((levelSets.back())->getNumberOfSegments()));

      // an iterator for each level set
      std::vector<hrleRunsIterator<typename lsDomain<T, D>::DomainType>>
          iterators;
      for (auto it = levelSets.begin(); it != levelSets.end(); ++it) {
        iterators.push_back(
            hrleRunsIterator<typename lsDomain<T, D>::DomainType>(
                (*it)->getDomain()));
      }

      IntegrationSchemeType scheme(IntegrationScheme);

      for (hrleRunsIterator<typename lsDomain<T, D>::DomainType> it(
               topDomain, startVector);
           it.getStartIndices() < endVector; ++it) {

        if (!it.isDefined() || std::abs(it.getValue()) > 0.5)
          continue;

        const T value = it.getValue();
        double maxStepTime = 0;
        double cfl = timeStepRatio;

        for (int currentLevelSetId = levelSets.size() - 1;
             currentLevelSetId >= 0; --currentLevelSetId) {

          T velocity = 0;

          // check if there is any other levelset at the same point:
          // if yes, take the velocity of the lowest levelset
          for (unsigned lowerLevelSetId = 0; lowerLevelSetId < levelSets.size();
               ++lowerLevelSetId) {
            // put iterator to same position as the top levelset
            iterators[lowerLevelSetId].goToIndicesSequential(
                it.getStartIndices());

            // if the lower surface is actually outside, i.e. its LS value is
            // lower or equal
            if (iterators[lowerLevelSetId].getValue() <= value) {
              velocity =
                  scheme(it.getStartIndices(), velocities, lowerLevelSetId);
              break;
            }
          }

          T valueBelow;
          // get value of material below (before in levelSets list)
          if (currentLevelSetId > 0) {
            // do not need to go to indices, since we should already be there
            // iterators[currentLevelSetId - 1].goToIndicesSequential(
            //     it.getStartIndices());
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
    lsReduce<T, D>(*levelSets.back()).apply(1, true);

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

    return maxTimeStep;
  }

public:
  lsAdvect(std::vector<lsDomain<T, D> *> &passedlsDomains,
           lsVelocityField<T> &passedVelocities)
      : levelSets(passedlsDomains), velocities(&passedVelocities) {}

  void setTimeStepRatio(const double &cfl) { timeStepRatio = cfl; }

  void setCalculateNormalVectors(bool cnv) { calculateNormalVectors = cnv; }

  double getTimeStepRatio() { return timeStepRatio; }

  bool getCalculateNormalVectors() { return calculateNormalVectors; }

  void setIntegrationScheme(unsigned scheme) {
    integrationScheme = static_cast<IntegrationSchemeEnum>(scheme);
  }

  void setIntegrationScheme(IntegrationSchemeEnum scheme) {
    integrationScheme = scheme;
  }

  void setDissipationAlpha(const double &a) { dissipationAlpha = a; }

  double apply() { return advect(); }

  unsigned apply(const double timeDelta) {
    double currentTime = 0.;
    unsigned counter = 0;
    while (currentTime < timeDelta) {
      currentTime += advect(timeDelta - currentTime);
      ++counter;
    }
    return counter;
  }
};

#endif // LS_ADVECT_TEMPLATE_HPP
