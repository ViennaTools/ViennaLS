#ifndef LS_CALCULATE_CURVATURES_HPP
#define LS_CALCULATE_CURVATURES_HPP

#include <lsDomain.hpp>

#include <hrleCartesianPlaneIterator.hpp>

#include <lsCurvatureFormulas.hpp>

enum struct lsCurvatureType : unsigned {
  MEAN_CURVATURE = 0,
  GAUSSIAN_CURVATURE = 1,
  MEAN_AND_GAUSSIAN_CURVATURE = 2
};

// Calculates the Curvature 2D, Mean Curvature and/or Gaussian Curvature (3D)
// for the passed lsDomain for all points with level set values <= 0.5. The
// result is saved in the lsDomain. requires a levelset width >=4.

template <class T, int D> class lsCalculateCurvatures {

  lsSmartPointer<lsDomain<T, D>> levelSet = nullptr;
  T maxValue = 0.5;
  lsCurvatureType type = lsCurvatureType::MEAN_CURVATURE;

public:
  lsCalculateCurvatures() {}

  lsCalculateCurvatures(lsSmartPointer<lsDomain<T, D>> passedLevelSet,
                        T passedMaxValue = 0.5)
      : levelSet(passedLevelSet), maxValue(passedMaxValue) {}

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  void setCurvatureType(lsCurvatureType passedType) {
    if (D == 2) {
      return;
    }
    type = passedType;
  }

  void setMaxValue(const T passedMaxValue) { maxValue = passedMaxValue; }

  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsCalculateCurvatures.")
          .print();
    }

    if (levelSet->getLevelSetWidth() < (maxValue * 5) + 1) {
      lsMessage::getInstance()
          .addWarning("lsCalculateCurvatures: Level set width must be "
                      "greater than " +
                      std::to_string((maxValue * 5) + 1) + " !")
          .print();
    }

    std::vector<std::vector<T>> meanCurvaturesVector(
        levelSet->getNumberOfSegments());
    std::vector<std::vector<T>> gaussCurvaturesVector(
        levelSet->getNumberOfSegments());
    double pointsPerSegment =
        double(2 * levelSet->getDomain().getNumberOfPoints()) /
        double(levelSet->getLevelSetWidth());

    auto grid = levelSet->getGrid();

    //! Calculate Curvatures
#pragma omp parallel num_threads(levelSet->getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      lsInternal::curvaturGeneralFormula<T, D> curvatureCalculator(
          levelSet->getGrid().getGridDelta());

      auto &meanCurvatures = meanCurvaturesVector[p];
      auto &gaussCurvatures = gaussCurvaturesVector[p];

      if (type == lsCurvatureType::MEAN_CURVATURE) {
        meanCurvatures.reserve(pointsPerSegment);
      } else if (type == lsCurvatureType::GAUSSIAN_CURVATURE) {
        gaussCurvatures.reserve(pointsPerSegment);
      } else {
        meanCurvatures.reserve(pointsPerSegment);
        gaussCurvatures.reserve(pointsPerSegment);
      }

      hrleVectorType<hrleIndexType, D> startVector =
          (p == 0) ? grid.getMinGridPoint()
                   : levelSet->getDomain().getSegmentation()[p - 1];

      hrleVectorType<hrleIndexType, D> endVector =
          (p != static_cast<int>(levelSet->getNumberOfSegments() - 1))
              ? levelSet->getDomain().getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      for (hrleCartesianPlaneIterator<typename lsDomain<T, D>::DomainType>
               neighborIt(levelSet->getDomain(), startVector);
           neighborIt.getIndices() < endVector; neighborIt.next()) {

        auto &center = neighborIt.getCenter();
        if (!center.isDefined()) {
          continue;
        } else if (std::abs(center.getValue()) > maxValue) {
          continue;
        }

        // calculate curvatures
        if (type == lsCurvatureType::MEAN_CURVATURE) {
          meanCurvatures.push_back(
              curvatureCalculator.meanCurvature(neighborIt));
        } else if (type == lsCurvatureType::GAUSSIAN_CURVATURE) {
          gaussCurvatures.push_back(
              curvatureCalculator.gaussianCurvature(neighborIt));
        } else {
          std::array<T, 2> curves =
              curvatureCalculator.meanGaussianCurvature(neighborIt);
          meanCurvatures.push_back(curves[0]);
          gaussCurvatures.push_back(curves[1]);
        }
      }
    }

    // put all curvature values in the correct ordering

    if (type == lsCurvatureType::MEAN_CURVATURE ||
        type == lsCurvatureType::MEAN_AND_GAUSSIAN_CURVATURE) {
      unsigned numberOfCurvatures = 0;
      for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i) {
        numberOfCurvatures += meanCurvaturesVector[i].size();
      }
      meanCurvaturesVector[0].reserve(numberOfCurvatures);

      for (unsigned i = 1; i < levelSet->getNumberOfSegments(); ++i) {
        meanCurvaturesVector[0].insert(meanCurvaturesVector[0].end(),
                                       meanCurvaturesVector[i].begin(),
                                       meanCurvaturesVector[i].end());
      }
    }

    if (type == lsCurvatureType::GAUSSIAN_CURVATURE ||
        type == lsCurvatureType::MEAN_AND_GAUSSIAN_CURVATURE) {
      unsigned numberOfCurvatures = 0;
      for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i) {
        numberOfCurvatures += gaussCurvaturesVector[i].size();
      }
      gaussCurvaturesVector[0].reserve(numberOfCurvatures);

      for (unsigned i = 1; i < levelSet->getNumberOfSegments(); ++i) {
        gaussCurvaturesVector[0].insert(gaussCurvaturesVector[0].end(),
                                        gaussCurvaturesVector[i].begin(),
                                        gaussCurvaturesVector[i].end());
      }
    }

    // insert into pointData of levelSet
    if (type == lsCurvatureType::MEAN_CURVATURE ||
        type == lsCurvatureType::MEAN_AND_GAUSSIAN_CURVATURE) {
      auto &pointData = levelSet->getPointData();
      auto scalarDataPointer = pointData.getScalarData("MeanCurvatures");
      // if it does not exist, insert new normals vector
      if (scalarDataPointer == nullptr) {
        pointData.insertNextScalarData(meanCurvaturesVector[0],
                                       "MeanCurvatures");
      } else {
        // if it does exist, just swap the old with the new values
        *scalarDataPointer = std::move(meanCurvaturesVector[0]);
      }
    }

    if (type == lsCurvatureType::GAUSSIAN_CURVATURE ||
        type == lsCurvatureType::MEAN_AND_GAUSSIAN_CURVATURE) {

      auto &pointData = levelSet->getPointData();
      auto scalarDataPointer = pointData.getScalarData("GaussianCurvatures");
      // if it does not exist, insert new normals vector
      if (scalarDataPointer == nullptr) {
        pointData.insertNextScalarData(gaussCurvaturesVector[0],
                                       "GaussianCurvatures");
      } else {
        // if it does exist, just swap the old with the new values
        *scalarDataPointer = std::move(gaussCurvaturesVector[0]);
      }
    }
  }
};

#endif // LS_CALCULATE_CURVATURES_HPP