#ifndef LS_FAST_ADVECT_HPP
#define LS_FAST_ADVECT_HPP

#include <unordered_map>

#include <hrleDenseIterator.hpp>
#include <hrleSparseIterator.hpp>
#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsMessage.hpp>
#include <lsPreCompileMacros.hpp>

#include <lsToDiskMesh.hpp>
// #include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFastAdvectDistributions.hpp>

#include <lsToMesh.hpp>
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

    // Extract the original surface as a point cloud of grid
    // points shifted to the surface (disk mesh)
    lsMesh surfaceMesh;
    lsToDiskMesh<T, D>(*levelSet, surfaceMesh).apply();
    // lsVTKWriter(surfaceMesh, "surfaceMesh.vtk").apply();
    typedef std::vector<std::array<double, 3>> SurfaceNodesType;
    const SurfaceNodesType &surfaceNodes = surfaceMesh.getNodes();

    auto &domain = levelSet->getDomain();

    auto &grid = levelSet->getGrid();
    auto gridDelta = grid.getGridDelta();

    // find bounds of distribution
    std::array<hrleCoordType, 2 * D> distBounds;
    dist->getBounds(distBounds);
    // check bounds of distribution since they must be > gridDelta
    for (unsigned i = 0; i < 2 * D; ++i) {
      if (std::abs(distBounds[i]) < gridDelta) {
        lsMessage::getInstance()
            .addWarning(
                "Distribution passed to lsFastAdvect is too small. It must be "
                "> gridDelta in every direction. Not performing Advection.")
            .print();
        return;
      }
    }

    // TODO: need to add support for periodic boundary conditions!
    hrleVectorType<hrleIndexType, D> distMin, distMax;

    std::cout << "minPoint: " << domain.getDomainSegment(0).definedValues[0] << std::endl;
    std::cout << "maxPoint: " << domain.getDomainSegment(domain.getNumberOfSegments()-1).definedValues.back() << std::endl;

    bool minPointNegative = domain.getDomainSegment(0).definedValues[0] < 0.;
    bool maxPointNegative = domain.getDomainSegment(domain.getNumberOfSegments()-1).definedValues.back() < 0.;

    // find bounding box of old domain
    hrleIndexType bounds[2 * D];
    domain.getDomainBounds(bounds);
    hrleVectorType<hrleIndexType, D> min, max;
    for (unsigned i = 0; i < D; ++i) {
      // translate from coords to indices
      distMin[i] = distBounds[2 * i] / gridDelta - 2;
      distMax[i] = distBounds[2 * i + 1] / gridDelta + 2;

      // std::cout << "dim " << i << ": " << std::endl;
      // std::cout << "Runbreaks: " << domain.getMinRunBreak(i) << ", " << domain.getMaxRunBreak(i) << std::endl;
      // std::cout << "MinGrid: " << grid.getMinGridPoint(i) << ", " << grid.getMaxGridPoint(i) << std::endl;

      // use the extent of the diskMesh to identify bounding box of new
      // level set
      // TODO: respect periodic boundary condition
      min[i] = surfaceMesh.minimumExtent[i] / gridDelta;
      if(grid.isNegBoundaryInfinite(i) && minPointNegative) {
        min[i] += distMax[i] - 2;
      } else {
        min[i] += distMin[i];
      }
      // if calculated index is out of bounds, set the extent
      // TODO: need to add periodic BNC handling here
      if(min[i] < grid.getMinGridPoint(i)) {
        min[i] = grid.getMinGridPoint(i);
      }

      max[i] = surfaceMesh.maximumExtent[i] / gridDelta;
      if(grid.isPosBoundaryInfinite(i) && maxPointNegative) {
        max[i] += distMin[i] + 2;
      } else {
        max[i] += distMax[i];
      }
      if(max[i] > grid.getMaxGridPoint(i)) {
        max[i] = grid.getMaxGridPoint(i);
      }

      std::cout << i << " min: " << min[i] << ", max: " << max[i] << std::endl;
    }
    std::cout << "segments: " << domain.getNumberOfSegments() << std::endl;
    // initialize with segmentation for whole range
    typename hrleDomain<T, D>::hrleIndexPoints segmentation;

    {
      unsigned long long numPoints = 1;
      unsigned long long pointsPerDimension[D];
      for (unsigned i = 0; i < D; ++i) {
        pointsPerDimension[i] = numPoints;
        numPoints *= max[i] - min[i];
        std::cout << "pointsPerDimension " << i << ": " << pointsPerDimension[i]
                  << std::endl;
      }
      unsigned long numberOfSegments = domain.getNumberOfSegments();
      unsigned long long pointsPerSegment = numPoints / numberOfSegments;
      std::cout << "numPoints: " << numPoints << std::endl;
      std::cout << "pointsPerSegment: " << pointsPerSegment
                << " * N: " << pointsPerSegment * numberOfSegments << std::endl;
      unsigned long long pointId = 0;
      for (unsigned i = 0; i < numberOfSegments - 1; ++i) {
        pointId = pointsPerSegment * (i + 1);
        hrleVectorType<hrleIndexType, D> segmentPoint;
        for (int j = D - 1; j >= 0; --j) {
          segmentPoint[j] = pointId / (pointsPerDimension[j]) + min[j];
          pointId %= pointsPerDimension[j];
        }
        std::cout << i << ": " << segmentPoint << std::endl;
        segmentation.push_back(segmentPoint);
      }
    }

    typedef std::vector<std::pair<hrleVectorType<hrleIndexType, D>, T>>
        PointValueVector;
    std::vector<PointValueVector> newPoints;
    newPoints.resize(domain.getNumberOfSegments());

    constexpr double cutoffValue = 1.0;

// set up multithreading
#pragma omp parallel num_threads(domain.getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      hrleVectorType<hrleIndexType, D> startVector =
          (p == 0) ? min : segmentation[p - 1];

      hrleVectorType<hrleIndexType, D> endVector =
          (p != static_cast<int>(domain.getNumberOfSegments() - 1))
              ? segmentation[p]
              : grid.incrementIndices(max);

      hrleConstSparseIterator<DomainType> checkIt(levelSet->getDomain(),
                                                  startVector);

      // Iterate through the bounds of new lsDomain lexicographically
      for (hrleVectorType<hrleIndexType, D> currentIndex = startVector;
           currentIndex <= endVector;
           incrementIndices(currentIndex, min, max)) {
        // if point is already full in old level set, skip it
        {
          checkIt.goToIndicesSequential(currentIndex);
          // if run is already negative undefined, just ignore the point
          if (checkIt.getValue() < -0.5) {
            continue;
          }
        }

        hrleVectorType<hrleCoordType, D> currentCoords;
        hrleVectorType<hrleCoordType, D> currentDistMin;
        hrleVectorType<hrleCoordType, D> currentDistMax;

        for (unsigned i = 0; i < D; ++i) {
          currentCoords[i] = currentIndex[i] * gridDelta;

          currentDistMin[i] = currentIndex[i] + distMin[i];
          if (currentDistMin[i] < grid.getMinGridPoint(i)) {
            currentDistMin[i] = grid.getMinGridPoint(i);
          }
          currentDistMin[i] *= gridDelta;

          currentDistMax[i] = currentIndex[i] + distMax[i];
          if (currentDistMin[i] > grid.getMaxGridPoint(i)) {
            currentDistMin[i] = grid.getMaxGridPoint(i);
          }
          currentDistMax[i] *= gridDelta;
        }

        double distance = std::numeric_limits<double>::max();

        // now check which surface points contribute to currentIndex
        for (typename SurfaceNodesType::const_iterator surfIt =
                 surfaceNodes.begin();
             surfIt != surfaceNodes.end(); ++surfIt) {

          auto &currentNode = *surfIt;

          // if we are outside min/max go to next index inside
          {
            bool outside = false;
            for (unsigned i = 0; i < D; ++i) {
              if ((currentNode[i] < currentDistMin[i]) ||
                  (currentNode[i] > currentDistMax[i])) {
                outside = true;
                break;
              }
            }
            if (outside) {
              continue;
            }
          }

          std::array<hrleCoordType, D> localCoords;
          for (unsigned i = 0; i < D; ++i) {
            localCoords[i] = currentCoords[i] - currentNode[i];
          }

          // TODO: does this really save time? Try without it.
          if (!dist->isInside(localCoords, 2 * gridDelta)) {
            continue;
          }

          // get filling fraction from distance to dist surface
          double tmpDistance = dist->getSignedDistance(localCoords) / gridDelta;

          // if cell is far within a distribution, set it filled
          if (tmpDistance < -cutoffValue) {
            distance = std::numeric_limits<T>::lowest();
            break;
          }

          // if distance is smaller, set the new one
          if (tmpDistance < distance) {
            distance = tmpDistance;
          }
        }

        if (std::abs(distance) <= cutoffValue) {
          newPoints[p].push_back(std::make_pair(currentIndex, distance));
        }

      } // domainBounds for
    }   // parallel region

    // copy all points into the first vector
    {
      std::cout << "Balancing: " << std::endl;
      std::cout << "0: " << newPoints[0].size() << std::endl;
      unsigned long long numberOfPoints = newPoints[0].size();
      for (unsigned i = 1; i < domain.getNumberOfSegments(); ++i) {
        numberOfPoints += newPoints[i].size();
      }
      newPoints[0].reserve(numberOfPoints);
      for (unsigned i = 1; i < domain.getNumberOfSegments(); ++i) {
        std::cout << i << ": " << newPoints[i].size() << std::endl;
        // newPoints[0].reserve(newPoints[0].size() + newPoints[i].size());
        std::move(std::begin(newPoints[i]), std::end(newPoints[i]),
                  std::back_inserter(newPoints[0]));
        // newPoints[i].clear();
      }
    }

    lsMesh mesh;
    // for(auto it = newPoints[0].begin(); it != newPoints[0].end(); ++it) {

    // }
    // lsDomain<T, D> newLevelSet(grid);
    // auto &newDomain = newLevelSet.getDomain();
    // lsToMesh<T, D>(*levelSet, mesh, false).apply();
    // lsVTKWriter(mesh, "beforeInsert.vtk").apply();

    levelSet->insertPoints(newPoints[0]);
    lsToMesh<T, D>(*levelSet, mesh, false).apply();
    lsVTKWriter(mesh, "afterInsert.vtk").apply();
    levelSet->setLevelSetWidth(1);
    lsExpand<T, D>(*levelSet, 2).apply();
    lsToMesh<T, D>(*levelSet, mesh, false).apply();
    lsVTKWriter(mesh, "afterExpand.vtk").apply();

    levelSet->getDomain().segment();
    levelSet->finalize(1);

    // levelSet->deepCopy(newLevelSet);
    lsExpand<T, D>(*levelSet, 2).apply();
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsFastAdvect)

#endif // LS_FAST_ADVECT_HPP
