#pragma once

#include <hrleSparseIterator.hpp>

#include <lsBooleanOperation.hpp>
#include <lsConcepts.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromMesh.hpp>
#include <lsGeometricAdvectDistributions.hpp>
#include <lsPreCompileMacros.hpp>
#include <lsToDiskMesh.hpp>

#include <vcLogger.hpp>
#include <vcSmartPointer.hpp>
#include <vcVectorType.hpp>

#ifndef NDEBUG // if in debug build
#include <lsCheck.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>
#endif

namespace viennals {

using namespace viennacore;

/// This class advects the level set according to a given distribution.
/// This distribution is overlayed at every grid point of the old surface. All
/// cells within this distribution are then filled, with cells at the edge
/// marked with the correct level set values. Therefore, the surface can be
/// shifted long distances in one step. This algorithm is therefore preferable
/// to normal advection if there is growth/reduction by a purely geometric
/// directional distribution.
template <class T, int D> class GeometricAdvect {
  using hrleIndexType = viennahrle::IndexType;
  using hrleCoordType = viennahrle::CoordType;

  SmartPointer<Domain<T, D>> levelSet = nullptr;
  SmartPointer<Domain<T, D>> maskLevelSet = nullptr;
  SmartPointer<const GeometricAdvectDistribution<hrleCoordType, D>> dist =
      nullptr;
  static constexpr T cutoffValue =
      T(1.) + std::numeric_limits<T>::epsilon() * T(100);

  static void incrementIndices(viennahrle::Index<D> &indices,
                               const viennahrle::Index<D> &min,
                               const viennahrle::Index<D> &max) {
    int dim = 0;
    for (; dim < D - 1; ++dim) {
      if (indices[dim] < max[dim])
        break;
      indices[dim] = min[dim];
    }
    ++indices[dim];
  }

  template <class K, class V, template <class...> class MapType, class... Ts>
  MapType<V, K> inverseTranslator(MapType<K, V, Ts...> &map) {
    MapType<V, K> inv;
    std::for_each(map.begin(), map.end(), [&inv](const std::pair<K, V> &p) {
      inv.insert(std::make_pair(p.second, p.first));
    });
    return inv;
  }

public:
  GeometricAdvect() = default;

  template <class DistType,
            lsConcepts::IsBaseOf<GeometricAdvectDistribution<hrleCoordType, D>,
                                 DistType> = lsConcepts::assignable>
  GeometricAdvect(SmartPointer<Domain<T, D>> passedLevelSet,
                  SmartPointer<DistType> passedDist,
                  SmartPointer<Domain<T, D>> passedMaskLevelSet = nullptr)
      : levelSet(passedLevelSet), maskLevelSet(passedMaskLevelSet) {
    dist = std::dynamic_pointer_cast<
        GeometricAdvectDistribution<hrleCoordType, D>>(passedDist);
  }

  /// Set the levelset which should be advected.
  void setLevelSet(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  /// Set which advection distribution to use. Must be derived from
  /// GeometricAdvectDistribution.
  template <class DistType,
            lsConcepts::IsBaseOf<GeometricAdvectDistribution<hrleCoordType, D>,
                                 DistType> = lsConcepts::assignable>
  void setAdvectionDistribution(SmartPointer<DistType> passedDist) {
    dist = std::dynamic_pointer_cast<
        GeometricAdvectDistribution<hrleCoordType, D>>(passedDist);
  }

  /// Set the levelset, which should be used as a mask. This level set
  /// has to be wrapped by the levelset set by setLevelSet, so the mask
  /// is entirely inside the advected level set.
  void setMaskLevelSet(SmartPointer<Domain<T, D>> passedMaskLevelSet) {
    maskLevelSet = passedMaskLevelSet;
  }

  /// Perform geometrical advection.
  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addWarning("No level set passed to GeometricAdvect. Not Advecting.")
          .print();
      return;
    }
    if (dist == nullptr) {
      Logger::getInstance()
          .addWarning("No GeometricAdvectDistribution passed to "
                      "GeometricAdvect. Not "
                      "Advecting.")
          .print();
      return;
    }

    // levelSet must have at least a width of 3
    Expand<T, D>(levelSet, 3).apply();

    if (maskLevelSet != nullptr) {
      Expand<T, D>(maskLevelSet, 3).apply();
    }

    typedef typename Domain<T, D>::DomainType DomainType;

    auto &domain = levelSet->getDomain();

    auto &grid = levelSet->getGrid();
    auto gridDelta = grid.getGridDelta();

    // Extract the original surface as a point cloud of grid
    // points shifted to the surface (disk mesh)
    auto surfaceMesh = SmartPointer<Mesh<hrleCoordType>>::New();
    auto pointIdTranslator =
        SmartPointer<typename ToDiskMesh<T, D>::TranslatorType>::New();
    ToDiskMesh<T, D, hrleCoordType>(levelSet, surfaceMesh, pointIdTranslator)
        .apply();
    *pointIdTranslator = inverseTranslator(*pointIdTranslator);

    // find bounds of distribution
    auto distBounds = dist->getBounds();

    // TODO: need to add support for periodic boundary conditions!
    viennahrle::Index<D> distMin, distMax;

    bool minPointNegative = domain.getDomainSegment(0).definedValues[0] < 0.;
    bool maxPointNegative =
        domain.getDomainSegment(domain.getNumberOfSegments() - 1)
            .definedValues.back() < 0.;
    bool distIsPositive = true;

    // find bounding box of old domain
    hrleIndexType bounds[6];
    domain.getDomainBounds(bounds);
    viennahrle::Index<D> min, max;
    for (unsigned i = 0; i < D; ++i) {
      // translate from coords to indices
      distMin[i] =
          distBounds[2 * i] / gridDelta + ((distBounds[2 * i] < 0) ? -2 : 2);
      distMax[i] = distBounds[2 * i + 1] / gridDelta +
                   ((distBounds[2 * i + 1] < 0) ? -2 : 2);
      if (distBounds[2 * i] >= 0) {
        distIsPositive = false;
      }

      // use the extent of the diskMesh to identify bounding box of new
      // level set
      // TODO: respect periodic boundary condition
      min[i] = surfaceMesh->minimumExtent[i] / gridDelta;
      // TODO also do the same thing for positive point and etching
      if (grid.isNegBoundaryInfinite(i) && minPointNegative && distMin[i] < 0) {
        min[i] -= 2;
      } else {
        if (distIsPositive) {
          min[i] += distMin[i];
        } else {
          min[i] -= distMin[i];
        }
      }
      // if calculated index is out of bounds, set the extent
      // TODO: need to add periodic BNC handling here
      if (min[i] < grid.getMinGridPoint(i)) {
        min[i] = grid.getMinGridPoint(i);
      }

      max[i] = surfaceMesh->maximumExtent[i] / gridDelta;
      if (grid.isPosBoundaryInfinite(i) && maxPointNegative && distMax[i] > 0) {
        max[i] += 2;
      } else {
        if (distIsPositive) {
          max[i] += distMax[i];
        } else {
          max[i] -= distMax[i];
        }
      }
      if (max[i] > grid.getMaxGridPoint(i)) {
        max[i] = grid.getMaxGridPoint(i);
      }
    }

    // Remove contribute points if they are part of the mask
    // If a mask is supplied, remove all contribute points which
    // lie on (or inside) the mask
    if (maskLevelSet != nullptr) {
      // Go over all contribute points and see if they are on the mask surface
      auto &maskDomain = maskLevelSet->getDomain();
      auto values = surfaceMesh->cellData.getScalarData("LSValues");
      auto valueIt = values->begin();

      auto newSurfaceMesh = SmartPointer<Mesh<hrleCoordType>>::New();
      PointData<hrleCoordType>::ScalarDataType newValues;
      viennahrle::ConstSparseIterator<DomainType> maskIt(maskDomain);
      for (auto &node : surfaceMesh->getNodes()) {
        viennahrle::Index<D> index;
        for (unsigned i = 0; i < D; ++i) {
          index[i] = std::round(node[i] / gridDelta);
        }
        // can do sequential, because surfaceNodes are lexicographically sorted
        // from lsToDiskMesh
        maskIt.goToIndicesSequential(index);
        // if it is a mask point, mark it to maybe use it in new level set
        if (!maskIt.isDefined() || !(maskIt.getValue() < *valueIt + 1e-5)) {
          newSurfaceMesh->insertNextNode(node);
          newValues.push_back(*valueIt);
          // insert vertex
          std::array<unsigned, 1> vertex{};
          vertex[0] = newSurfaceMesh->nodes.size();
          newSurfaceMesh->insertNextVertex(vertex);
        }
        ++valueIt;
      }
      newSurfaceMesh->cellData.insertNextScalarData(newValues, "LSValues");
      // use new mesh as surfaceMesh
      newSurfaceMesh->minimumExtent = surfaceMesh->minimumExtent;
      newSurfaceMesh->maximumExtent = surfaceMesh->maximumExtent;
      surfaceMesh = newSurfaceMesh;
    }

#ifndef NDEBUG // if in debug build
    {
      Logger::getInstance()
          .addDebug("GeomAdvect: Writing debug meshes")
          .print();
      VTKWriter<hrleCoordType>(surfaceMesh, FileFormatEnum::VTP,
                               "DEBUG_lsGeomAdvectMesh_contributewoMask.vtp")
          .apply();
      auto mesh = SmartPointer<Mesh<T>>::New();
      if (maskLevelSet != nullptr) {
        ToMesh<T, D>(maskLevelSet, mesh).apply();
        VTKWriter<T>(mesh, FileFormatEnum::VTP,
                     "DEBUG_lsGeomAdvectMesh_mask.vtp")
            .apply();
      }
      ToMesh<T, D>(levelSet, mesh).apply();
      VTKWriter<T>(mesh, FileFormatEnum::VTP,
                   "DEBUG_lsGeomAdvectMesh_initial.vtp")
          .apply();
    }

#endif

    const auto &surfaceNodes = surfaceMesh->getNodes();

    // initialize with segmentation for whole range
    typename viennahrle::Domain<T, D>::IndexPoints segmentation;

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
        viennahrle::Index<D> segmentPoint;
        for (int j = D - 1; j >= 0; --j) {
          segmentPoint[j] = pointId / (pointsPerDimension[j]) + min[j];
          pointId %= pointsPerDimension[j];
        }
        segmentation.push_back(segmentPoint);
      }
    }

    typedef std::vector<std::pair<viennahrle::Index<D>, T>> PointValueVector;
    std::vector<PointValueVector> newPoints;
    newPoints.resize(domain.getNumberOfSegments());

    const T initialDistance = (distIsPositive)
                                  ? std::numeric_limits<double>::max()
                                  : std::numeric_limits<double>::lowest();

#ifndef NDEBUG
    {
      std::ostringstream oss;
      oss << "GeomAdvect: Min: " << min << ", Max: " << max << std::endl;
      Logger::getInstance().addDebug(oss.str()).print();
    }
#endif
// set up multithreading
#pragma omp parallel num_threads(domain.getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      viennahrle::Index<D> startVector;
      if (p == 0) {
        startVector = min;
      } else {
        startVector = segmentation[p - 1];
        incrementIndices(startVector, min, max);
      }

      viennahrle::Index<D> endVector =
          (p != static_cast<int>(domain.getNumberOfSegments() - 1))
              ? segmentation[p]
              : grid.incrementIndices(max);

      viennahrle::ConstSparseIterator<DomainType> checkIt(levelSet->getDomain(),
                                                          startVector);

      // Mask iterator for checking whether inside mask or not
      SmartPointer<viennahrle::ConstSparseIterator<DomainType>> maskIt =
          nullptr;
      if (maskLevelSet != nullptr) {
        maskIt = SmartPointer<viennahrle::ConstSparseIterator<DomainType>>::New(
            maskLevelSet->getDomain(), startVector);
      }

      // Iterate through the bounds of new lsDomain lexicographically
      for (viennahrle::Index<D> currentIndex = startVector;
           currentIndex <= endVector;
           incrementIndices(currentIndex, min, max)) {
        // if point is already full in old level set, skip it
        checkIt.goToIndicesSequential(currentIndex);
        T oldValue = checkIt.getValue();
        // if run is already negative undefined, just ignore the point
        if (distIsPositive) {
          if (oldValue < -cutoffValue) {
            continue;
          }
        } else if (oldValue > cutoffValue) {
          continue;
        }

        VectorType<hrleCoordType, 3> currentCoords{};
        VectorType<hrleCoordType, 3> currentDistMin{};
        VectorType<hrleCoordType, 3> currentDistMax{};

        for (unsigned i = 0; i < D; ++i) {
          currentCoords[i] = currentIndex[i] * gridDelta;

          currentDistMin[i] = currentIndex[i] - std::abs(distMin[i]);
          if (currentDistMin[i] < grid.getMinGridPoint(i)) {
            currentDistMin[i] = grid.getMinGridPoint(i);
          }
          currentDistMin[i] *= gridDelta;

          currentDistMax[i] = currentIndex[i] + std::abs(distMax[i]);
          if (currentDistMin[i] > grid.getMaxGridPoint(i)) {
            currentDistMin[i] = grid.getMaxGridPoint(i);
          }
          currentDistMax[i] *= gridDelta;
        }

        T distance = initialDistance;

        unsigned long currentPointId = 0;
        // now check which surface points contribute to currentIndex
        for (auto surfIt = surfaceNodes.begin(); surfIt != surfaceNodes.end();
             ++surfIt, ++currentPointId) {

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
          T tmpDistance = dist->getSignedDistance(
                              currentNode, currentCoords,
                              pointIdTranslator->find(currentPointId)->second) /
                          gridDelta;

          // if cell is far within a distribution, set it filled
          if (distIsPositive) {
            if (tmpDistance <= -cutoffValue) {
              distance = std::numeric_limits<T>::lowest();
              break;
            }

            if (tmpDistance < distance) {
              distance = tmpDistance;
            }
          } else {
            if (tmpDistance >= cutoffValue) {
              distance = std::numeric_limits<T>::max();
              break;
            }

            if (tmpDistance > distance) {
              distance = tmpDistance;
            }
          }
        }

        // TODO: There are still issues with positive box distributions
        // if there is a mask used!
        // if point is part of the mask, keep smaller value
        if (maskLevelSet != nullptr) {
          maskIt->goToIndicesSequential(currentIndex);

          // if dist is positive, flip logic of comparison
          if (distIsPositive ^
              (std::abs(oldValue - maskIt->getValue()) < 1e-6)) {
            if (!distIsPositive && std::abs(oldValue) <= cutoffValue) {
              newPoints[p].push_back(std::make_pair(currentIndex, oldValue));
              continue;
            }
          } else {
            if (distance != initialDistance) {
              distance = std::min(maskIt->getValue(), distance);
            } else if (distIsPositive || oldValue >= 0.) {
              newPoints[p].push_back(std::make_pair(currentIndex, oldValue));
              continue;
            }
          }
        }

        if (std::abs(distance) <= cutoffValue) {
          // avoid using distribution in wrong direction
          if (distIsPositive && oldValue >= 0.) {
            newPoints[p].push_back(std::make_pair(currentIndex, distance));
          } else if (!distIsPositive && oldValue <= 0.) {
            // if we are etching, need to make sure, we are not inside mask
            if (maskIt == nullptr || maskIt->getValue() > -cutoffValue) {
              newPoints[p].push_back(std::make_pair(currentIndex, distance));
            }
          } else {
            // this only happens if distribution is very small, < 2 * gridDelta
            newPoints[p].push_back(std::make_pair(currentIndex, oldValue));
          }
        }
      }
    }

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

    auto mesh = SmartPointer<Mesh<T>>::New();
    // output all points directly to mesh
    {
      std::vector<T> scalarData;
      for (auto it = newPoints[0].begin(); it != newPoints[0].end(); ++it) {
        Vec3D<T> node{0., 0., 0.};
        for (unsigned i = 0; i < D; ++i) {
          node[i] = T((it->first)[i]) * gridDelta;
        }

        mesh->insertNextNode(node);
        std::array<unsigned, 1> vertex{};
        vertex[0] = mesh->vertices.size();
        mesh->insertNextVertex(vertex);
        scalarData.push_back(it->second);
      }
      mesh->cellData.insertNextScalarData(scalarData, "LSValues");
    }

#ifndef NDEBUG // if in debug build
    Logger::getInstance().addDebug("GeomAdvect: Writing final mesh...").print();
    VTKWriter<T>(mesh, FileFormatEnum::VTP, "DEBUG_lsGeomAdvectMesh_final.vtp")
        .apply();
#endif

    FromMesh<T, D>(levelSet, mesh).apply();

#ifndef NDEBUG // if in debug build
    Logger::getInstance().addDebug("GeomAdvect: Writing final LS...").print();
    ToMesh<T, D>(levelSet, mesh).apply();
    VTKWriter<T>(mesh, FileFormatEnum::VTP, "DEBUG_lsGeomAdvectLS_final.vtp")
        .apply();
#endif

    Prune<T, D>(levelSet).apply();

    levelSet->getDomain().segment();
    levelSet->finalize(1);

    Expand<T, D>(levelSet, 2).apply();
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(GeometricAdvect)

} // namespace viennals
