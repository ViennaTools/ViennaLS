#ifndef LS_FAST_ADVECT_HPP
#define LS_FAST_ADVECT_HPP

#include <unordered_map>

#include <hrleDenseIterator.hpp>
#include <hrleSparseIterator.hpp>
#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsMessage.hpp>
#include <lsPreCompileMacros.hpp>

#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsGeometricAdvectDistributions.hpp>
#include <lsFromMesh.hpp>
#include <lsToDiskMesh.hpp>

// TODO: remove
#include <lsVTKWriter.hpp>

/// This class advects the level set according to a given distribution.
/// This distribution is overlayed at every cell. All cells within
/// this distribution are then filled, with cells at the edge marked
/// with the correct level set values. Therefore, the surface can
/// be shifted long distances in one step. This algorithm is therefore
/// preferable to normal advection if there is growth/reduction by a geometric
/// directional distribution.
template <class T, int D> class lsGeometricAdvect {
  lsDomain<T, D> *levelSet = nullptr;
  const lsGeometricAdvectDistribution<hrleCoordType, D> *dist = nullptr;
  static constexpr T numericEps = 10 * std::numeric_limits<T>::epsilon();

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
  lsGeometricAdvect() {}

  template <class DistType>
  lsGeometricAdvect(lsDomain<T, D> &passedLevelSet, DistType &passedDist)
      : levelSet(&passedLevelSet), dist(&passedDist) {}

  void setLevelSet(lsDomain<T, D> &passedLevelSet) {
    levelSet = &passedLevelSet;
  }

  void setAdvectionDistribution(
      const lsGeometricAdvectDistribution<hrleCoordType, D> &distribution) {
    dist = &distribution;
  }

  // iterate through all points of new cell set and check whether distributions
  // on the old cell set will set the point
  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set passed to lsGeometricAdvect. Not Advecting.")
          .print();
      return;
    }
    if (dist == nullptr) {
      lsMessage::getInstance()
          .addWarning("No lsGeometricAdvectDistribution passed to lsGeometricAdvect. Not "
                      "Advecting.")
          .print();
      return;
    }

    typedef typename lsDomain<T, D>::DomainType DomainType;

    // Extract the original surface as a point cloud of grid
    // points shifted to the surface (disk mesh)
    lsMesh surfaceMesh;
    lsToDiskMesh<T, D>(*levelSet, surfaceMesh).apply();
    typedef std::vector<std::array<double, 3>> SurfaceNodesType;
    const SurfaceNodesType &surfaceNodes = surfaceMesh.getNodes();

    auto &domain = levelSet->getDomain();

    auto &grid = levelSet->getGrid();
    auto gridDelta = grid.getGridDelta();

    // find bounds of distribution
    std::array<hrleCoordType, 6> distBounds;
    dist->getBounds(distBounds);

    // TODO: need to add support for periodic boundary conditions!
    hrleVectorType<hrleIndexType, D> distMin, distMax;

    bool minPointNegative = domain.getDomainSegment(0).definedValues[0] < 0.;
    bool maxPointNegative =
        domain.getDomainSegment(domain.getNumberOfSegments() - 1)
            .definedValues.back() < 0.;

    // find bounding box of old domain
    hrleIndexType bounds[6];
    domain.getDomainBounds(bounds);
    hrleVectorType<hrleIndexType, D> min, max;
    for (unsigned i = 0; i < D; ++i) {
      // translate from coords to indices
      distMin[i] = distBounds[2 * i] / gridDelta - 2;
      distMax[i] = distBounds[2 * i + 1] / gridDelta + 2;

      // use the extent of the diskMesh to identify bounding box of new
      // level set
      // TODO: respect periodic boundary condition
      min[i] = surfaceMesh.minimumExtent[i] / gridDelta;
      if (grid.isNegBoundaryInfinite(i) && minPointNegative) {
        min[i] -=  2;
      } else {
        min[i] += distMin[i];
      }
      // if calculated index is out of bounds, set the extent
      // TODO: need to add periodic BNC handling here
      if (min[i] < grid.getMinGridPoint(i)) {
        min[i] = grid.getMinGridPoint(i);
      }

      max[i] = surfaceMesh.maximumExtent[i] / gridDelta;
      if (grid.isPosBoundaryInfinite(i) && maxPointNegative) {
        max[i] += 2;
      } else {
        max[i] += distMax[i];
      }
      if (max[i] > grid.getMaxGridPoint(i)) {
        max[i] = grid.getMaxGridPoint(i);
      }
    }
    // initialize with segmentation for whole range
    typename hrleDomain<T, D>::hrleIndexPoints segmentation;

    {
      unsigned long long numPoints = 1;
      unsigned long long pointsPerDimension[D];
      for (unsigned i = 0; i < D; ++i) {
        pointsPerDimension[i] = numPoints;
        numPoints *= max[i] - min[i];
      }
      unsigned long numberOfSegments = domain.getNumberOfSegments();
      unsigned long long pointsPerSegment = numPoints / numberOfSegments;
      unsigned long long pointId = 0;
      for (unsigned i = 0; i < numberOfSegments - 1; ++i) {
        pointId = pointsPerSegment * (i + 1);
        hrleVectorType<hrleIndexType, D> segmentPoint;
        for (int j = D - 1; j >= 0; --j) {
          segmentPoint[j] = pointId / (pointsPerDimension[j]) + min[j];
          pointId %= pointsPerDimension[j];
        }
        segmentation.push_back(segmentPoint);
      }
    }

    typedef std::vector<std::pair<hrleVectorType<hrleIndexType, D>, T>>
        PointValueVector;
    std::vector<PointValueVector> newPoints;
    newPoints.resize(domain.getNumberOfSegments());

    constexpr T cutoffValue = 1.0 + numericEps;

// set up multithreading
#pragma omp parallel num_threads(domain.getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      hrleVectorType<hrleIndexType, D> startVector;
      if (p == 0) {
        startVector = min;
      } else {
        startVector = segmentation[p - 1];
        incrementIndices(startVector, min, max);
      }

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

        std::array<hrleCoordType, 3> currentCoords;
        std::array<hrleCoordType, 3> currentDistMin;
        std::array<hrleCoordType, 3> currentDistMax;

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

        T distance = std::numeric_limits<double>::max();

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

          // TODO: does this really save time? Try without it.
          if (!dist->isInside(currentNode, currentCoords, 2 * gridDelta)) {
            continue;
          }

          // get filling fraction from distance to dist surface
          T tmpDistance = dist->getSignedDistance(currentNode, currentCoords) / gridDelta;

          // if cell is far within a distribution, set it filled
          if (tmpDistance <= -cutoffValue) {
            distance = std::numeric_limits<T>::lowest();
            break;
          }

          // if distance is smaller, set the new one
          if (tmpDistance < distance) {
            distance = tmpDistance;
          }
        }

        if (std::abs(distance) <= cutoffValue) {
          newPoints[p].push_back(
              std::make_pair(currentIndex, distance - numericEps));
        }

      } // domainBounds for
    }   // parallel region

    // copy all points into the first vector
    {
      unsigned long long numberOfPoints = newPoints[0].size();
      for (unsigned i = 1; i < domain.getNumberOfSegments(); ++i) {
        numberOfPoints += newPoints[i].size();
      }
      newPoints[0].reserve(numberOfPoints);
      for (unsigned i = 1; i < domain.getNumberOfSegments(); ++i) {
        std::move(std::begin(newPoints[i]), std::end(newPoints[i]),
                  std::back_inserter(newPoints[0]));
      }
    }

    lsMesh mesh;
    // output all points directly to mesh
    {
      std::vector<double> scalarData;
      for (auto it = newPoints[0].begin(); it != newPoints[0].end(); ++it) {
        std::array<double, 3> node = {};
        for (unsigned i = 0; i < D; ++i) {
          node[i] = T((it->first)[i]) * gridDelta;
        }

        mesh.insertNextNode(node);
        std::array<unsigned, 1> vertex;
        vertex[0] = mesh.vertices.size();
        mesh.insertNextVertex(vertex);
        scalarData.push_back(it->second);
      }
      mesh.insertNextScalarData(scalarData, "LSValues");
    }

    // TODO: remove
    lsVTKWriter(mesh, "afterAdvect.vtk").apply();

    lsFromMesh<T, D>(*levelSet, mesh).apply();
    lsPrune<T, D>(*levelSet).apply();

    levelSet->getDomain().segment();
    levelSet->finalize(1);

    lsExpand<T, D>(*levelSet, 2).apply();
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsGeometricAdvect)

#endif // LS_FAST_ADVECT_HPP
