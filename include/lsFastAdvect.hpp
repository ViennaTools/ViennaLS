#ifndef LS_FAST_ADVECT_HPP
#define LS_FAST_ADVECT_HPP

#include <unordered_map>

#include <hrleDenseIterator.hpp>
#include <hrleSparseIterator.hpp>
#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsMessage.hpp>

#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFastAdvectDistributions.hpp>
#include <lsToDiskMesh.hpp>
#include <lsVTKWriter.hpp>

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

    // calculate normal vectors and save in a hash
    std::unordered_map<hrleVectorType<hrleIndexType, D>,
                       hrleVectorType<hrleCoordType, D>,
                       typename hrleVectorType<hrleIndexType, D>::hash>
        normalVectors;
    for (hrleConstSparseStarIterator<DomainType> neighborIt(
             levelSet->getDomain());
         !neighborIt.isFinished(); ++neighborIt) {
      auto &center = neighborIt.getCenter();
      if (!center.isDefined() || std::abs(center.getValue()) > 0.5)
        continue;

      hrleVectorType<hrleCoordType, D> n;
      T denominator = 0;
      for (int i = 0; i < D; i++) {
        n[i] = (neighborIt.getNeighbor(i).getValue() -
                neighborIt.getNeighbor(i + D).getValue()) *
               0.5;
        denominator += n[i] * n[i];
      }

      denominator = std::sqrt(denominator);
      if (std::abs(denominator) < 1e-12) {
        for (unsigned i = 0; i < D; ++i)
          n[i] = 0.;
      } else {
        for (unsigned i = 0; i < D; ++i) {
          n[i] /= denominator;
        }
      }

      auto pair =
          normalVectors.insert(std::make_pair(center.getStartIndices(), n));
      if (!pair.second) {
        std::cout << "Could not insert " << center.getStartIndices()
                  << " because of " << (pair.first)->first << std::endl;
      }
    }

    auto &domain = levelSet->getDomain();

    auto &grid = levelSet->getGrid();
    auto gridDelta = grid.getGridDelta();
    lsDomain<T, D> newLevelSet(grid);
    auto &newDomain = newLevelSet.getDomain();

    newDomain.initialize(domain.getNewSegmentation(), domain.getAllocation());

    // find bounds of distribution
    hrleCoordType distBounds[2 * D];
    dist->getBounds(distBounds);
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
      if(domain.getNumberOfSegments() != 0) {
        auto segment = domain.getDomainSegment(0);
        if (segment.runTypes[0].size() != 0) {
          // set the first runtype of old domain
          auto firstRunValue = segment.undefinedValues[hrleRunTypeValues::UNDEF_PT - segment.runTypes[0][0]];
          newDomain.insertNextUndefinedPoint(0, grid.getMinGridPoint(),
                                             firstRunValue);
        }
      } else {
        newDomain.insertNextUndefinedPoint(0, grid.getMinGridPoint(),
                                           levelSet->POS_VALUE);
      }
    }

    // Iterate through the bounds of new lsDomain lexicographically
    for (hrleVectorType<hrleIndexType, D> currentIndex = min;
         currentIndex <= max; incrementIndices(currentIndex, min, max)) {
      // if point is already full in old level set, skip it
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
            // std::cout << "N ";
          } else {
            // std::cout << "# ";
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
          for(unsigned i = 0; i < D; ++i) {
            if(distIndex[i] < currentDistMin[i] || distIndex[i] > currentDistMax[i]) {
              outside = true;
            }
          }
          if(outside) {
            incrementIndices(distIndex, currentDistMin, currentDistMax);
            --distIndex[0];
            distIt.goToIndices(distIndex);
            continue;
          }
        }

        hrleVectorType<hrleCoordType, D> distCoords;

        auto normalsIt = normalVectors.find(distIndex);
        auto distNormal = normalsIt->second;
        double vectorMax = 0.;
        for (unsigned i = 0; i < D; ++i) {
          distCoords[i] = distIndex[i] * gridDelta;
          if (std::abs(distNormal[i]) > vectorMax) {
            vectorMax = std::abs(distNormal[i]);
          }
        }
        for (unsigned i = 0; i < D; ++i) {
          // shift distcoords to surface from grid point
          distCoords[i] -= distIt.getValue() * gridDelta * distNormal[i] * vectorMax;
        }

        hrleVectorType<hrleCoordType, D> localCoords =
            currentCoords - distCoords;

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
          // std::cout << "F ";
        } else {
          // std::cout << "_ ";
        }
      } else if (distance >= 0.5) {
        if (!currentEmptyRun) {
          newDomain.insertNextUndefinedPoint(0, currentIndex,
                                             levelSet->POS_VALUE);
          currentEmptyRun = true;
          currentFullRun = false;
          // std::cout << "E ";
        } else {
          // std::cout << ". ";
        }
      } else {
        newDomain.insertNextDefinedPoint(0, currentIndex, distance);
        currentEmptyRun = false;
        currentFullRun = false;
        // std::cout << "D ";
      }

    } // domainBounds for
    // std::cout << std::endl;

    // insert final undefined run
    hrleVectorType<hrleIndexType, D> finalRun = grid.getMaxGridPoint();
    ++finalRun[D - 1];
    if(currentEmptyRun) {
      newDomain.insertNextUndefinedPoint(0, finalRun, levelSet->POS_VALUE);
    } else {
      newDomain.insertNextUndefinedPoint(0, finalRun, levelSet->NEG_VALUE);
    }

    // newDomain.print();

    newDomain.finalize();
    levelSet->deepCopy(newLevelSet);
    lsExpand<T, D>(*levelSet, 2).apply();
  }
};

#endif // LS_FAST_ADVECT_HPP
