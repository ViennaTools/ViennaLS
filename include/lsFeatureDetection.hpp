#ifndef LS_FEATURE_DETECTION_HPP
#define LS_FEATURE_DETECTION_HPP

#include <lsDomain.hpp>

#include <hrleSparseStarIterator.hpp>

#include <hrleSparseBoxIterator.hpp>

#include <hrleCartesianPlaneIterator.hpp>

#include <lsCurvatureFormulas.hpp>

enum struct FeatureDetectionMethod : unsigned {
  CURVATURE = 0,
  NORMALS_ANGLE = 1,
};

// Detects Features of the level set function for more information see private
// functions. This class offers two methods, curvature based feature detection
// should always be the preferred choice.

template <class T, int D> class lsFeatureDetection {

private:
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;

  lsSmartPointer<lsDomain<T, D>> levelSet = nullptr;

  FeatureDetectionMethod method = FeatureDetectionMethod::CURVATURE;

  T flatBoundary = 0.;

  std::string outputName;

  std::vector<T> flaggedCells;

public:
  lsFeatureDetection(
      lsSmartPointer<lsDomain<T, D>> passedLevelSet, T passedBoundary,
      FeatureDetectionMethod passedMethod = FeatureDetectionMethod::CURVATURE,
      std::string passedOutputName = "Features")
      : levelSet(passedLevelSet), flatBoundary(passedBoundary),
        method(passedMethod), outputName(passedOutputName) {}

  void setOutputName(std::string passedOutputName) {
    outputName = passedOutputName;
  }

  void setFeatureDetectionMethod(FeatureDetectionMethod passedMethod) {
    method = passedMethod;
  }

  void apply() {

    if (method == FeatureDetectionMethod::CURVATURE) {
      FeatureDetectionCurvature();
    } else {
      FeatureDetectionNormals();
    }

    // insert into pointData of levelSet
    auto &pointData = levelSet->getPointData();
    auto vectorDataPointer = pointData.getScalarData(outputName);
    // if it does not exist, insert new normals vector
    if (vectorDataPointer == nullptr) {
      pointData.insertNextScalarData(flaggedCells, outputName);
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
    flaggedCells.reserve(levelSet->getNumberOfPoints());

    auto grid = levelSet->getGrid();

    typename lsDomain<T, D>::DomainType &domain = levelSet->getDomain();

    std::vector<std::vector<T>> flagsReserve(levelSet->getNumberOfSegments());

    double pointsPerSegment =
        double(2 * levelSet->getDomain().getNumberOfPoints()) /
        double(levelSet->getLevelSetWidth());

    // In 2 Dimensions there are no minimal Surfaces -> Do not need to calculate
    // Gaussian Curvature
    if (D == 2) {

#pragma omp parallel num_threads((levelSet)->getNumberOfSegments())
      {
        int p = 0;
#ifdef _OPENMP
        p = omp_get_thread_num();
#endif

        lsInternal::curvaturGeneralFormula<T, D> curvatureCalculator(
            levelSet->getGrid().getGridDelta());

        auto &flagsSegment = flagsReserve[p];
        flagsSegment.reserve(pointsPerSegment);

        hrleVectorType<hrleIndexType, D> startVector =
            (p == 0) ? grid.getMinGridPoint() : domain.getSegmentation()[p - 1];

        hrleVectorType<hrleIndexType, D> endVector =
            (p != static_cast<int>(domain.getNumberOfSegments() - 1))
                ? domain.getSegmentation()[p]
                : grid.incrementIndices(grid.getMaxGridPoint());

        for (hrleCartesianPlaneIterator<typename lsDomain<T, D>::DomainType>
                 neighborIt(levelSet->getDomain(), startVector, 1);
             neighborIt.getIndices() < endVector; neighborIt.next()) {

          if (!neighborIt.getCenter().isDefined() ||
              std::abs(neighborIt.getCenter().getValue()) > 0.5) {
            continue;
          }

          T curve = curvatureCalculator.meanCurvature(neighborIt);

          if (std::abs(curve) > flatBoundary) {
            flagsSegment.push_back(1);
          } else {
            flagsSegment.push_back(0);
          }
        }
      }

    } else {

#pragma omp parallel num_threads((levelSet)->getNumberOfSegments())
      {
        int p = 0;
#ifdef _OPENMP
        p = omp_get_thread_num();
#endif

        lsInternal::curvaturGeneralFormula<T, D> curvatureCalculator(
            levelSet->getGrid().getGridDelta());

        auto &flagsSegment = flagsReserve[p];
        flagsSegment.reserve(pointsPerSegment);

        hrleVectorType<hrleIndexType, D> startVector =
            (p == 0) ? grid.getMinGridPoint() : domain.getSegmentation()[p - 1];

        hrleVectorType<hrleIndexType, D> endVector =
            (p != static_cast<int>(domain.getNumberOfSegments() - 1))
                ? domain.getSegmentation()[p]
                : grid.incrementIndices(grid.getMaxGridPoint());

        for (hrleCartesianPlaneIterator<typename lsDomain<T, D>::DomainType>
                 neighborIt(levelSet->getDomain(), startVector, 1);
             neighborIt.getIndices() < endVector; neighborIt.next()) {

          if (!neighborIt.getCenter().isDefined() ||
              std::abs(neighborIt.getCenter().getValue()) > 0.5) {
            continue;
          }

          std::array<T, 2> curves =
              curvatureCalculator.meanGaussianCurvature(neighborIt);

          // Minimal surfaces can have Mean Curvature equal to 0 at non-flat
          // points on the surface additionally check Gaussian Curvature if the
          // point is flat
          if (std::abs(curves[0]) > flatBoundary ||
              std::abs(curves[1]) > flatBoundary) {
            flagsSegment.push_back(1);
          } else {
            flagsSegment.push_back(0);
          }
        }
      }
    }

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

    T cosAngleTreshold = std::cos(flatBoundary);
    ;

    std::vector<hrleVectorType<hrleIndexType, D>> combinations;

    if (D == 3) {

      for (int i = 0; i < D; i++) {

        hrleVectorType<hrleIndexType, D> posUnit(0);
        hrleVectorType<hrleIndexType, D> negUnit(0);

        int first_pos = i;
        int second_pos = (i + 1) % D;

        posUnit[first_pos] = 1;
        negUnit[first_pos] = -1;

        combinations.push_back(posUnit);
        combinations.push_back(negUnit);

        posUnit[second_pos] = 1;
        negUnit[second_pos] = 1;

        combinations.push_back(posUnit);
        combinations.push_back(negUnit);

        posUnit[second_pos] = -1;
        negUnit[second_pos] = -1;

        combinations.push_back(posUnit);
        combinations.push_back(negUnit);
      }

      combinations.push_back(hrleVectorType<hrleIndexType, D>(1, 1, 1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(-1, -1, -1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(-1, -1, 1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(1, -1, -1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(-1, 1, -1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(-1, 1, 1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(1, -1, 1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(1, 1, -1));
    } else {

      combinations.push_back(hrleVectorType<hrleIndexType, D>(1, 1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(-1, -1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(-1, 1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(1, -1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(0, 1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(1, 0));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(0, -1));
      combinations.push_back(hrleVectorType<hrleIndexType, D>(-1, 0));
    }

    std::vector<std::array<T, D>> normals;

    std::vector<std::vector<std::array<T, D>>> normalsVector(
        levelSet->getNumberOfSegments());

    normals.reserve(levelSet->getNumberOfPoints());
    flaggedCells.reserve(levelSet->getNumberOfPoints());

    auto grid = levelSet->getGrid();

    typename lsDomain<T, D>::DomainType &domain = levelSet->getDomain();

    std::vector<std::vector<T>> flagsReserve(levelSet->getNumberOfSegments());

    double pointsPerSegment =
        double(2 * levelSet->getDomain().getNumberOfPoints()) /
        double(levelSet->getLevelSetWidth());

    // Calculate all Surface normals
#pragma omp parallel num_threads((levelSet)->getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      std::array<T, D> zeroVector;

      for (int i = 0; i < D; i++) {
        zeroVector[i] = 0.;
      }

      auto &normalsSegment = normalsVector[p];
      normalsSegment.reserve(pointsPerSegment);

      hrleVectorType<hrleIndexType, D> startVector =
          (p == 0) ? grid.getMinGridPoint() : domain.getSegmentation()[p - 1];

      hrleVectorType<hrleIndexType, D> endVector =
          (p != static_cast<int>(domain.getNumberOfSegments() - 1))
              ? domain.getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      for (hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType>
               neighborIt(levelSet->getDomain(), startVector);
           neighborIt.getIndices() < endVector; neighborIt.next()) {

        if (!neighborIt.getCenter().isDefined()) {
          continue;
        } else if (std::abs(neighborIt.getCenter().getValue()) >= 0.5) {
          normalsSegment.push_back(zeroVector);
          continue;
        }

        std::array<T, D> n;

        T norm = 0.;

        for (int i = 0; i < D; i++) {

          T pos1 = neighborIt.getNeighbor(i).getValue();

          T neg1 = neighborIt.getNeighbor(i + D).getValue();

          n[i] = (pos1 - neg1) * 0.5;

          norm += n[i] * n[i];
        }

        norm = std::sqrt(norm);

        for (int j = 0; j < D; j++)
          n[j] /= norm;

        normalsSegment.push_back(n);
      }
    }

    for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i)
      normals.insert(normals.end(), normalsVector[i].begin(),
                     normalsVector[i].end());

      // Compare angles between normal vectors
#pragma omp parallel num_threads((levelSet)->getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      std::array<T, D> zeroVector;

      for (int i = 0; i < D; i++) {
        zeroVector[i] = 0.;
      }

      std::vector<T> &flagsSegment = flagsReserve[p];
      flagsSegment.reserve(pointsPerSegment);

      hrleVectorType<hrleIndexType, D> startVector =
          (p == 0) ? grid.getMinGridPoint() : domain.getSegmentation()[p - 1];

      hrleVectorType<hrleIndexType, D> endVector =
          (p != static_cast<int>(domain.getNumberOfSegments() - 1))
              ? domain.getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      for (hrleSparseBoxIterator<typename lsDomain<T, D>::DomainType>
               neighborIt(levelSet->getDomain(), startVector, 1);
           neighborIt.getIndices() < endVector; neighborIt.next()) {

        if (!neighborIt.getCenter().isDefined() ||
            std::abs(neighborIt.getCenter().getValue()) >= 0.5) {
          continue;
        }

        std::array<T, D> centerNormal =
            normals[neighborIt.getCenter().getPointId()];

        bool flag = false;

        for (unsigned dir = 0; dir < (D * D * D); dir++) {

          std::array<T, D> currentNormal =
              normals[neighborIt.getNeighbor(dir).getPointId()];

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

#endif // LS_FEATURE_DETECTION_HPP