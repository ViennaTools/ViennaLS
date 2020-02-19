#ifndef LS_CALCULATE_NORMAL_VECTORS_HPP
#define LS_CALCULATE_NORMAL_VECTORS_HPP

#include <lsPreCompileMacros.hpp>

#include <algorithm>

#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsDomain.hpp>
#include <lsMessage.hpp>

/// This algorithm is used to compute the normal vectors for all points
/// with level set values <= 0.5. The result is saved in the lsDomain and
/// can be retrieved with lsDomain.getNormalVectors().
/// Since neighbors in each cartesian direction are necessary for
/// the calculation, the levelset width must be >=3.
template <class T, int D> class lsCalculateNormalVectors {
  lsDomain<T, D> *levelSet = nullptr;
  T maxValue = 0.5;

public:
  lsCalculateNormalVectors() {}

  lsCalculateNormalVectors(lsDomain<T, D> &passedLevelSet,
                           T passedMaxValue = 0.5)
      : levelSet(&passedLevelSet), maxValue(passedMaxValue) {}

  void setLevelSet(lsDomain<T, D> &passedLevelSet) {
    levelSet = &passedLevelSet;
  }

  void setMaxValue(const T passedMaxValue) { maxValue = passedMaxValue; }

  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsCalculateNormalVectors.")
          .print();
    }

    if (levelSet->getLevelSetWidth() < (maxValue * 4) + 1) {
      lsMessage::getInstance()
          .addWarning("lsCalculateNormalVectors: Level set width must be "
                      "greater than " +
                      std::to_string((maxValue * 4) + 1) + " 2!")
          .print();
    }

    std::vector<std::vector<std::array<T, D>>> normalVectorsVector(
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

      std::vector<std::array<T, D>> &normalVectors = normalVectorsVector[p];
      normalVectors.reserve(pointsPerSegment);

      hrleVectorType<hrleIndexType, D> startVector =
          (p == 0) ? grid.getMinGridPoint()
                   : levelSet->getDomain().getSegmentation()[p - 1];

      hrleVectorType<hrleIndexType, D> endVector =
          (p != static_cast<int>(levelSet->getNumberOfSegments() - 1))
              ? levelSet->getDomain().getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      for (hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType>
               neighborIt(levelSet->getDomain(), startVector);
           neighborIt.getIndices() < endVector; neighborIt.next()) {

        auto &center = neighborIt.getCenter();
        if (!center.isDefined()) {
          continue;
        } else if (std::abs(center.getValue()) > maxValue) {
          // push an empty vector to keep ordering correct
          std::array<T, D> tmp = {};
          normalVectors.push_back(tmp);
          continue;
        }

        std::array<T, D> n;

        T denominator = 0;
        for (int i = 0; i < D; i++) {
          T pos = neighborIt.getNeighbor(i).getValue() - center.getValue();
          T neg = center.getValue() - neighborIt.getNeighbor(i + D).getValue();
          n[i] = (pos + neg) * 0.5;
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

        normalVectors.push_back(n);
      }
    }

    // copy all normals
    auto &normals = levelSet->getNormalVectors();
    normals.clear();
    unsigned numberOfNormals = 0;
    for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i) {
      numberOfNormals += normalVectorsVector[i].size();
    }
    normals.reserve(numberOfNormals);

    for (unsigned i = 0; i < levelSet->getNumberOfSegments(); ++i) {
      normals.insert(normals.end(), normalVectorsVector[i].begin(),
                     normalVectorsVector[i].end());
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsCalculateNormalVectors)

#endif // LS_CALCULATE_NORMAL_VECTORS_HPP
