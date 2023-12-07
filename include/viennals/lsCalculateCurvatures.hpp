#ifndef LS_CALCULATE_CURVATURES_HPP
#define LS_CALCULATE_CURVATURES_HPP

#include <hrleCartesianPlaneIterator.hpp>
#include <lsCurvatureFormulas.hpp>
#include <lsDomain.hpp>

enum struct lsCurvatureEnum : unsigned {
  MEAN_CURVATURE = 0,
  GAUSSIAN_CURVATURE = 1,
  MEAN_AND_GAUSSIAN_CURVATURE = 2
};

// Calculates the Mean Curvature and/or Gaussian Curvature (3D)
// for the passed lsDomain for all points with level set values <= 0.5. The
// result is saved in the lsDomain.
template <class T, int D> class lsCalculateCurvatures {
  lsSmartPointer<lsDomain<T, D>> levelSet = nullptr;
  T maxValue = 0.5;
  lsCurvatureEnum type = lsCurvatureEnum::MEAN_CURVATURE;

public:
  static constexpr char meanCurvatureLabel[] = "MeanCurvatures";
  static constexpr char gaussianCurvatureLabel[] = "GaussianCurvatures";

  lsCalculateCurvatures() {}

  lsCalculateCurvatures(lsSmartPointer<lsDomain<T, D>> passedLevelSet)
      : levelSet(passedLevelSet) {}

  lsCalculateCurvatures(lsSmartPointer<lsDomain<T, D>> passedLevelSet,
                        lsCurvatureEnum method)
      : levelSet(passedLevelSet), type(method) {}

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  void setCurvatureType(lsCurvatureEnum passedType) {
    // in 2D there is only one option so ignore
    if constexpr (D == 3) {
      type = passedType;
    } else {
      if (passedType != type) {
        lsMessage::getInstance().addWarning(
            "lsCalculateCurvatures: Could not set curvature type because 2D "
            "only supports mean curvature.");
      }
    }
  }

  void setMaxValue(const T passedMaxValue) { maxValue = passedMaxValue; }

  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsCalculateCurvatures.")
          .print();
    }

    // need second neighbours
    if (unsigned minWidth = std::ceil((maxValue * 8) + 1);
        levelSet->getLevelSetWidth() < minWidth) {
      lsMessage::getInstance()
          .addWarning("lsCalculateCurvatures: Level set width must be "
                      "at least " +
                      std::to_string(minWidth) + " !")
          .print();
    }

    std::vector<std::vector<T>> meanCurvaturesVector(
        levelSet->getNumberOfSegments());
    std::vector<std::vector<T>> gaussCurvaturesVector(
        levelSet->getNumberOfSegments());

    auto grid = levelSet->getGrid();
    const bool calculateMean =
        (type == lsCurvatureEnum::MEAN_CURVATURE) ||
        (type == lsCurvatureEnum::MEAN_AND_GAUSSIAN_CURVATURE);
    const bool calculateGauss =
        (type == lsCurvatureEnum::GAUSSIAN_CURVATURE) ||
        (type == lsCurvatureEnum::MEAN_AND_GAUSSIAN_CURVATURE);

    //! Calculate Curvatures
#pragma omp parallel num_threads(levelSet->getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      auto &meanCurvatures = meanCurvaturesVector[p];
      auto &gaussCurvatures = gaussCurvaturesVector[p];

      if (calculateMean) {
        meanCurvatures.reserve(
            levelSet->getDomain().getDomainSegment(p).getNumberOfPoints());
      }
      if (calculateGauss) {
        gaussCurvatures.reserve(
            levelSet->getDomain().getDomainSegment(p).getNumberOfPoints());
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
          if (calculateMean)
            meanCurvatures.push_back(0.);
          if (calculateGauss)
            gaussCurvatures.push_back(0.);
          continue;
        }

        // calculate curvatures
        if (calculateMean) {
          meanCurvatures.push_back(lsInternal::meanCurvature(neighborIt));
        }
        if (calculateGauss) {
          gaussCurvatures.push_back(lsInternal::gaussianCurvature(neighborIt));
        }
      }
    }

    // put all curvature values in the correct ordering
    if (calculateMean) {
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

    if (calculateGauss) {
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
    if (calculateMean) {
      auto &pointData = levelSet->getPointData();
      auto scalarDataPointer = pointData.getScalarData(meanCurvatureLabel);
      // if it does not exist, insert new normals vector
      if (scalarDataPointer == nullptr) {
        pointData.insertNextScalarData(meanCurvaturesVector[0],
                                       meanCurvatureLabel);
      } else {
        // if it does exist, just swap the old with the new values
        *scalarDataPointer = std::move(meanCurvaturesVector[0]);
      }
    }

    if (calculateGauss) {
      auto &pointData = levelSet->getPointData();
      auto scalarDataPointer = pointData.getScalarData(gaussianCurvatureLabel);
      // if it does not exist, insert new normals vector
      if (scalarDataPointer == nullptr) {
        pointData.insertNextScalarData(gaussCurvaturesVector[0],
                                       gaussianCurvatureLabel);
      } else {
        // if it does exist, just swap the old with the new values
        *scalarDataPointer = std::move(gaussCurvaturesVector[0]);
      }
    }
  }
};

#endif // LS_CALCULATE_CURVATURES_HPP