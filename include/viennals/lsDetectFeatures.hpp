#pragma once

#include <hrleCartesianPlaneIterator.hpp>
#include <hrleSparseBoxIterator.hpp>
#include <lsCalculateNormalVectors.hpp>
#include <lsCurvatureFormulas.hpp>
#include <lsDomain.hpp>

#include <vcSmartPointer.hpp>
#include <vcVectorType.hpp>

namespace viennals {

using namespace viennacore;

enum struct FeatureDetectionEnum : unsigned {
  CURVATURE = 0,
  NORMALS_ANGLE = 1,
};

/// This class detects features of the level set function. This class offers two
/// methods to determine features of the surface: based on the mean curvature,
/// and based on the angle between surface normals. The curvature-based
/// algorithm is the default as it leads to more accurate results and should be
/// preferred in general.
template <class T, int D> class DetectFeatures {
  typedef typename Domain<T, D>::DomainType hrleDomainType;
  SmartPointer<Domain<T, D>> levelSet = nullptr;
  FeatureDetectionEnum method = FeatureDetectionEnum::CURVATURE;
  T flatLimit = 1.;
  T flatLimit2 = 1.;
  std::vector<T> flaggedCells;

public:
  static constexpr char featureMarkersLabel[] = "FeatureMarkers";

  DetectFeatures() = default;

  DetectFeatures(SmartPointer<Domain<T, D>> passedLevelSet)
      : levelSet(passedLevelSet) {}

  DetectFeatures(SmartPointer<Domain<T, D>> passedLevelSet, T passedLimit)
      : levelSet(passedLevelSet), flatLimit(passedLimit),
        flatLimit2(flatLimit * flatLimit) {}

  DetectFeatures(SmartPointer<Domain<T, D>> passedLevelSet, T passedLimit,
                 FeatureDetectionEnum passedMethod)
      : levelSet(passedLevelSet), method(passedMethod), flatLimit(passedLimit),
        flatLimit2(flatLimit * flatLimit) {}

  void setDetectionThreshold(T threshold) {
    flatLimit = threshold;
    flatLimit2 = flatLimit * flatLimit;
  }

  /// Set which algorithm to use to detect features. The curvature-based
  /// algorithm should always be preferred, while the normals-based algorithm is
  /// just provided for experimental use.
  void setDetectionMethod(FeatureDetectionEnum passedMethod) {
    method = passedMethod;
  }

  /// Execute the algorithm.
  void apply() {
    if (method == FeatureDetectionEnum::CURVATURE) {
      FeatureDetectionCurvature();
    } else {
      FeatureDetectionNormals();
    }

    // insert into pointData of levelSet
    auto &pointData = levelSet->getPointData();
    auto vectorDataPointer = pointData.getScalarData(featureMarkersLabel, true);
    // if it does not exist, insert new feature vector
    if (vectorDataPointer == nullptr) {
      pointData.insertNextScalarData(flaggedCells, featureMarkersLabel);
    } else {
      // if it does exist, just swap the old with the new values
      *vectorDataPointer = std::move(flaggedCells);
    }
  }

private:
  // Detects Features of the level set by calculating the absolute Curvature of
  // each active grid point (levelset value <= 0.5). In 3D the Gaussian
  // Curvature is also calculated to detect minimal surfaces. The minimal
  // curvature value that should be considered a feature is passed to the
  // constructor 0.0 Curvature describes a flat plane, the bigger the passed
  // parameter gets the more grid points will be detected as features.
  void FeatureDetectionCurvature() {
    flaggedCells.clear();

    auto grid = levelSet->getGrid();
    typename Domain<T, D>::DomainType &domain = levelSet->getDomain();
    std::vector<std::vector<T>> flagsReserve(levelSet->getNumberOfSegments());

#pragma omp parallel num_threads((levelSet)->getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      auto &flagsSegment = flagsReserve[p];
      flagsSegment.reserve(
          levelSet->getDomain().getDomainSegment(p).getNumberOfPoints());

      viennahrle::Index<D> const startVector =
          (p == 0) ? grid.getMinGridPoint() : domain.getSegmentation()[p - 1];

      viennahrle::Index<D> const endVector =
          (p != static_cast<int>(domain.getNumberOfSegments() - 1))
              ? domain.getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      for (viennahrle::CartesianPlaneIterator<typename Domain<T, D>::DomainType,
                                              1>
               neighborIt(levelSet->getDomain(), startVector);
           neighborIt.getIndices() < endVector; neighborIt.next()) {

        auto &center = neighborIt.getCenter();
        if (!center.isDefined()) {
          continue;
        } else if (std::abs(center.getValue()) > 0.5) {
          flagsSegment.push_back(0);
          continue;
        }

        T curve = lsInternal::meanCurvature(neighborIt);
        if (std::abs(curve) > flatLimit) {
          flagsSegment.push_back(1);
        } else {
          if constexpr (D == 2) {
            flagsSegment.push_back(0);
          } else {
            curve = lsInternal::gaussianCurvature(neighborIt);
            if (std::abs(curve) > flatLimit2)
              flagsSegment.push_back(1);
            else
              flagsSegment.push_back(0);
          }
        }
      }
    }

    flaggedCells.reserve(levelSet->getNumberOfPoints());
    for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i)
      flaggedCells.insert(flaggedCells.end(), flagsReserve[i].begin(),
                          flagsReserve[i].end());
  }

  // Detects Features of the level set by comparing the angle of each normal
  // vector on the surface to its adjacent normal vectors. The minimal angle
  // that should be considered a feature is passed to the class constructor.
  void FeatureDetectionNormals() {
    // Clear results from previous run
    flaggedCells.clear();

    auto &grid = levelSet->getGrid();
    auto &domain = levelSet->getDomain();
    T cosAngleTreshold = std::cos(flatLimit);

    // CALCULATE NORMALS
    Expand<T, D>(levelSet, 3).apply();
    CalculateNormalVectors<T, D>(levelSet).apply();
    const auto &normals = *(levelSet->getPointData().getVectorData(
        CalculateNormalVectors<T, D>::normalVectorsLabel));

    std::vector<std::vector<T>> flagsReserve(levelSet->getNumberOfSegments());

    // Compare angles between normal vectors
#pragma omp parallel num_threads((levelSet)->getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      Vec3D<T> zeroVector{};

      std::vector<T> &flagsSegment = flagsReserve[p];
      flagsSegment.reserve(
          levelSet->getDomain().getDomainSegment(p).getNumberOfPoints());

      viennahrle::Index<D> const startVector =
          (p == 0) ? grid.getMinGridPoint() : domain.getSegmentation()[p - 1];

      viennahrle::Index<D> const endVector =
          (p != static_cast<int>(domain.getNumberOfSegments() - 1))
              ? domain.getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      for (viennahrle::SparseBoxIterator<typename Domain<T, D>::DomainType, 1>
               neighborIt(levelSet->getDomain(), startVector);
           neighborIt.getIndices() < endVector; neighborIt.next()) {
        if (!neighborIt.getCenter().isDefined()) {
          continue;
        } else if (std::abs(neighborIt.getCenter().getValue()) >= 0.5) {
          flagsSegment.push_back(0);
          continue;
        }

        Vec3D<T> centerNormal = normals[neighborIt.getCenter().getPointId()];

        bool flag = false;

        constexpr unsigned numNeighbors = (D == 3) ? 27 : 9;
        for (unsigned dir = 0; dir < numNeighbors; dir++) {
          auto neighbor = neighborIt.getNeighbor(dir);
          if (!neighbor.isDefined())
            continue;
          Vec3D<T> currentNormal = normals[neighbor.getPointId()];

          if (currentNormal != zeroVector) {
            T skp = 0.;
            // Calculate scalar product
            for (int j = 0; j < D; j++) {
              skp += currentNormal[j] * centerNormal[j];
            }
            // Vectors are normlized so skp = cos(alpha)
            if ((cosAngleTreshold - skp) >= 0.) {
              flag = true;
              break;
            }
          }
        }

        if (flag) {
          flagsSegment.push_back(1);
        } else {
          flagsSegment.push_back(0);
        }
      }
    }

    for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i)
      flaggedCells.insert(flaggedCells.end(), flagsReserve[i].begin(),
                          flagsReserve[i].end());
  }
};

} // namespace viennals
