#ifndef LS_FAST_ADVECT_HPP
#define LS_FAST_ADVECT_HPP

#include <unordered_map>

#include <hrleDenseIterator.hpp>
#include <hrleSparseIterator.hpp>
#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsMessage.hpp>
#include <lsPreCompileMacros.hpp>

#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFastAdvectDistributions.hpp>

/// This class advects the level set according to a given distribution.
/// This distribution is overlayed at every cell. All cells within
/// this distribution are then filled, with cells at the edge marked
/// with the correct level set values. Therefore, the surface can
/// be shifted long distances in one step. This algorithm is therefore
/// preferable to normal advection if there is growth/reduction by a geometrical
/// directional distribution.
template <class T, int D> class lsFastAdvect {
  lsDomain<T, D> *levelSet = nullptr;
  const lsFastAdvectDistribution<hrleCoordType, D> *dist = nullptr;

  static void incrementIndices(hrleVectorType<hrleIndexType, D> &indices,
                               const hrleVectorType<hrleIndexType, D> &min,
                               const hrleVectorType<hrleIndexType, D> &max) {
    int dim = 0;
    for (; dim < D - 1; ++dim) {
      if (indices[dim] < max[dim])
        break;
      indices[dim] = min[dim];
    }
    ++indices[dim];
  }

public:
  lsFastAdvect() {}

  template <class DistType>
  lsFastAdvect(lsDomain<T, D> &passedLevelSet, DistType &passedDist)
      : levelSet(&passedLevelSet), dist(&passedDist) {}

  void setLevelSet(lsDomain<T, D> &passedLevelSet) {
    levelSet = &passedLevelSet;
  }

  void setAdvectionDistribution(
      const lsFastAdvectDistribution<hrleCoordType, D> &distribution) {
    dist = &distribution;
  }

  // iterate through all points of new cell set and check whether distributions
  // on the old cell set will set the point
  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set passed to lsFastAdvect. Not Advecting.")
          .print();
      return;
    }
    if (dist == nullptr) {
      lsMessage::getInstance()
          .addWarning("No lsFastAdvectDistribution passed to lsFastAdvect. Not "
                      "Advecting.")
          .print();
      return;
    }

    typedef typename lsDomain<T, D>::DomainType DomainType;

    // expand for normalvector calculation
    lsExpand<T, D>(*levelSet, 3).apply();
    lsCalculateNormalVectors<T, D>(*levelSet).apply();

    auto &normalVectors = levelSet->getNormalVectors();

    auto &domain = levelSet->getDomain();

    auto &grid = levelSet->getGrid();
    auto gridDelta = grid.getGridDelta();
    lsDomain<T, D> newLevelSet(grid);
    auto &newDomain = newLevelSet.getDomain();

    // initialize as single threaded
    newDomain.initialize();

    // find bounds of distribution
    std::array<hrleCoordType, 2 * D> distBounds;
    dist->getBounds(distBounds);
    // check bounds of distribution since they must be > gridDelta
    for (unsigned i = 0; i < 2 * D; ++i) {
      if (std::abs(distBounds[i]) < gridDelta / 2) {
        lsMessage::getInstance()
            .addWarning(
                "Distribution passed to lsFastAdvect is too small. It must be "
                "> gridDelta / 2 in every direction. Not performing Advection.")
            .print();
        return;
      }
    }

    hrleVectorType<hrleIndexType, D> distMin, distMax;

    // find bounding box of old domain
    hrleIndexType bounds[2 * D];
    domain.getDomainBounds(bounds);
    hrleVectorType<hrleIndexType, D> min, max;
    for (unsigned i = 0; i < D; ++i) {
      // translate from coords to indices
      distMin[i] = distBounds[2 * i] / gridDelta - 1;
      distMax[i] = distBounds[2 * i + 1] / gridDelta + 1;

      min[i] =
          bounds[2 * i] + ((grid.isNegBoundaryInfinite(i)) ? distMin[i] : 0);
      max[i] = bounds[2 * i + 1] +
               ((grid.isPosBoundaryInfinite(i)) ? distMax[i] : 0);
    }

    bool currentFullRun = false;
    bool currentEmptyRun = false;

    hrleVectorType<hrleIndexType, D> lastIndex = min;

    // set first undefined run
    {
      // get iterator to min
      hrleConstSparseIterator<DomainType> checkIt(levelSet->getDomain());
      newDomain.insertNextUndefinedPoint(0, grid.getMinGridPoint(),
                                         checkIt.getValue());
    }

    // Iterate through the bounds of new lsDomain lexicographically
    for (hrleVectorType<hrleIndexType, D> currentIndex = min;
         currentIndex <= max; incrementIndices(currentIndex, min, max)) {
      // if point is already full in old level set, skip it
      // TODO: do not always initialize new iterator
      // but use the same and just increment it when necessary
      {
        hrleConstSparseIterator<DomainType> checkIt(levelSet->getDomain(),
                                                    currentIndex);
        // if run is negative undefined
        if (checkIt.getValue() < -0.5) {
          if (!currentFullRun) {
            newDomain.insertNextUndefinedPoint(0, currentIndex,
                                               levelSet->NEG_VALUE);
            currentFullRun = true;
            currentEmptyRun = false;
          }
          continue;
        }
      }

      hrleVectorType<hrleCoordType, D> currentCoords;
      hrleVectorType<hrleIndexType, D> currentDistMin;
      hrleVectorType<hrleIndexType, D> currentDistMax;

      for (unsigned i = 0; i < D; ++i) {
        currentCoords[i] = currentIndex[i] * gridDelta;

        currentDistMin[i] = currentIndex[i] + distMin[i];
        if (currentDistMin[i] < grid.getMinGridPoint(i)) {
          currentDistMin[i] = grid.getMinGridPoint(i);
        }
        currentDistMax[i] = currentIndex[i] + distMax[i];
        if (currentDistMin[i] > grid.getMaxGridPoint(i)) {
          currentDistMin[i] = grid.getMaxGridPoint(i);
        }
      }

      double distance = 0.5;

      // Now start iterator over space around current index
      for (hrleConstSparseIterator<DomainType> distIt(levelSet->getDomain(),
                                                      currentDistMin);
           distIt.getStartIndices() <= currentDistMax; ++distIt) {
        if (!distIt.isDefined() || std::abs(distIt.getValue()) > 0.5) {
          continue;
        }

        hrleVectorType<hrleIndexType, D> distIndex = distIt.getStartIndices();

        // if we are outside min/max go to next index inside
        {
          bool outside = false;
          for (unsigned i = 0; i < D; ++i) {
            if (distIndex[i] < currentDistMin[i] ||
                distIndex[i] > currentDistMax[i]) {
              outside = true;
            }
          }
          if (outside) {
            incrementIndices(distIndex, currentDistMin, currentDistMax);
            --distIndex[0];
            distIt.goToIndices(distIndex);
            continue;
          }
        }

        hrleVectorType<hrleCoordType, D> distCoords;

        auto distNormal = normalVectors[distIt.getPointId()];
        double vectorMax = 0.;
        for (unsigned i = 0; i < D; ++i) {
          distCoords[i] = distIndex[i] * gridDelta;
          if (std::abs(distNormal[i]) > vectorMax) {
            vectorMax = std::abs(distNormal[i]);
          }
        }
        for (unsigned i = 0; i < D; ++i) {
          // shift distcoords to surface from grid point
          distCoords[i] -=
              distIt.getValue() * gridDelta * distNormal[i] * vectorMax;
        }

        std::array<hrleCoordType, D> localCoords;
        for (unsigned i = 0; i < D; ++i) {
          localCoords[i] = currentCoords[i] - distCoords[i];
        }

        if (!dist->isInside(localCoords, 2 * gridDelta)) {
          continue;
        }

        // get filling fraction from distance to dist surface
        double tmpDistance = dist->getSignedDistance(localCoords) / gridDelta;

        // if cell is far within a distribution, set it filled
        if (tmpDistance < -0.5) {
          distance = -1.;
          break;
        }

        // if distance is smaller, set the new one
        if (tmpDistance < distance) {
          distance = tmpDistance;
        }
      }

      // if there is a lexicographical jump, add undefined runs
      if (currentIndex[1] > lastIndex[1] ||
          (D == 3 && currentIndex[2] > lastIndex[2])) {
        auto minRunIndex = currentIndex;

        for (unsigned i = 1; i < D; ++i) {
          if (currentIndex[i] > lastIndex[i]) {
            for (unsigned j = 0; j < i; ++j) {
              minRunIndex[j] = grid.getMinGridPoint(j);
            }
          }
        }

        if (currentEmptyRun) {
          newDomain.insertNextUndefinedPoint(0, minRunIndex,
                                             levelSet->POS_VALUE);
        } else {
          newDomain.insertNextUndefinedPoint(0, minRunIndex,
                                             levelSet->NEG_VALUE);
        }
        lastIndex = currentIndex;
      }

      // now set the correct runs for the calculated value
      if (distance < -0.5) {
        if (!currentFullRun) {
          newDomain.insertNextUndefinedPoint(0, currentIndex,
                                             levelSet->NEG_VALUE);
          currentFullRun = true;
          currentEmptyRun = false;
        }
      } else if (distance >= 0.5) {
        if (!currentEmptyRun) {
          newDomain.insertNextUndefinedPoint(0, currentIndex,
                                             levelSet->POS_VALUE);
          currentEmptyRun = true;
          currentFullRun = false;
        }
      } else {
        newDomain.insertNextDefinedPoint(0, currentIndex, distance);
        currentEmptyRun = false;
        currentFullRun = false;
      }

    } // domainBounds for

    // insert final undefined run
    hrleVectorType<hrleIndexType, D> finalRun = grid.getMaxGridPoint();
    ++finalRun[D - 1];
    if (currentEmptyRun) {
      newDomain.insertNextUndefinedPoint(0, finalRun, levelSet->POS_VALUE);
    } else {
      newDomain.insertNextUndefinedPoint(0, finalRun, levelSet->NEG_VALUE);
    }

    newDomain.finalize();
    newDomain.segment();
    levelSet->deepCopy(newLevelSet);
    lsExpand<T, D>(*levelSet, 2).apply();
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsFastAdvect)

#endif // LS_FAST_ADVECT_HPP
