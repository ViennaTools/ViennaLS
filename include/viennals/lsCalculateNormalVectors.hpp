#pragma once

#include <lsPreCompileMacros.hpp>

#include <algorithm>

#include <hrleSparseStarIterator.hpp>

#include <lsDomain.hpp>

#include <vcLogger.hpp>
#include <vcSmartPointer.hpp>
#include <vcVectorType.hpp>

namespace viennals {

using namespace viennacore;

/// This algorithm is used to compute the normal vectors for all points
/// with level set values <= 0.5. The result is saved in the lsPointData of the
/// lsDomain and can be retrieved with
/// lsDomain.getPointData().getVectorData("Normals"). Since neighbors in each
/// cartesian direction are necessary for the calculation, the level set width
/// must be >=3.
template <class T, int D> class CalculateNormalVectors {
  SmartPointer<Domain<T, D>> levelSet = nullptr;
  T maxValue = 0.5;

public:
  static constexpr char normalVectorsLabel[] = "Normals";

  CalculateNormalVectors() = default;

  CalculateNormalVectors(SmartPointer<Domain<T, D>> passedLevelSet,
                         T passedMaxValue = 0.5)
      : levelSet(passedLevelSet), maxValue(passedMaxValue) {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  void setMaxValue(const T passedMaxValue) { maxValue = passedMaxValue; }

  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addWarning("No level set was passed to CalculateNormalVectors.")
          .print();
    }

    if (levelSet->getLevelSetWidth() < (maxValue * 4) + 1) {
      Logger::getInstance()
          .addWarning("CalculateNormalVectors: Level set width must be "
                      "greater than " +
                      std::to_string((maxValue * 4) + 1) + "!")
          .print();
    }

    std::vector<std::vector<Vec3D<T>>> normalVectorsVector(
        levelSet->getNumberOfSegments());
    double pointsPerSegment =
        double(2 * levelSet->getDomain().getNumberOfPoints()) /
        double(levelSet->getLevelSetWidth());

    auto grid = levelSet->getGrid();

    //! Calculate Normalvectors
#pragma omp parallel num_threads(levelSet->getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      auto &normalVectors = normalVectorsVector[p];
      normalVectors.reserve(pointsPerSegment);

      viennahrle::Index<D> const startVector =
          (p == 0) ? grid.getMinGridPoint()
                   : levelSet->getDomain().getSegmentation()[p - 1];

      viennahrle::Index<D> const endVector =
          (p != static_cast<int>(levelSet->getNumberOfSegments() - 1))
              ? levelSet->getDomain().getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      for (viennahrle::ConstSparseStarIterator<
               typename Domain<T, D>::DomainType, 1>
               neighborIt(levelSet->getDomain(), startVector);
           neighborIt.getIndices() < endVector; neighborIt.next()) {

        auto &center = neighborIt.getCenter();
        if (!center.isDefined()) {
          continue;
        } else if (std::abs(center.getValue()) > maxValue) {
          // push an empty vector to keep ordering correct
          Vec3D<T> tmp = {};
          normalVectors.push_back(tmp);
          continue;
        }

        Vec3D<T> n;

        T denominator = 0;
        for (int i = 0; i < D; i++) {
          T pos = neighborIt.getNeighbor(i).getValue() - center.getValue();
          T neg = center.getValue() - neighborIt.getNeighbor(i + D).getValue();
          n[i] = (pos + neg) * 0.5;
          denominator += n[i] * n[i];
        }

        denominator = std::sqrt(denominator);
        if (std::abs(denominator) < 1e-12) {
          std::ostringstream oss;
          oss << "CalculateNormalVectors: Vector of length 0 at "
              << neighborIt.getIndices();
          Logger::getInstance().addWarning(oss.str()).print();
          for (unsigned i = 0; i < D; ++i)
            n[i] = 0.;
        } else {
          for (unsigned i = 0; i < D; ++i) {
            n[i] /= denominator;
          }
        }

        normalVectors.push_back(n);
      }
    }

    // copy all normals
    unsigned numberOfNormals = 0;
    for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i) {
      numberOfNormals += normalVectorsVector[i].size();
    }
    normalVectorsVector[0].reserve(numberOfNormals);

    for (unsigned i = 1; i < levelSet->getNumberOfSegments(); ++i) {
      normalVectorsVector[0].insert(normalVectorsVector[0].end(),
                                    normalVectorsVector[i].begin(),
                                    normalVectorsVector[i].end());
    }

    // insert into pointData of levelSet
    auto &pointData = levelSet->getPointData();
    auto vectorDataPointer = pointData.getVectorData(normalVectorsLabel, true);
    // if it does not exist, insert new normals vector
    if (vectorDataPointer == nullptr) {
      pointData.insertNextVectorData(normalVectorsVector[0],
                                     normalVectorsLabel);
    } else {
      // if it does exist, just swap the old with the new values
      *vectorDataPointer = std::move(normalVectorsVector[0]);
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(CalculateNormalVectors)

} // namespace viennals
