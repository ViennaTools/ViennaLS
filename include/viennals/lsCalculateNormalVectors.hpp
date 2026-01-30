#pragma once

#include <lsPreCompileMacros.hpp>

#include <algorithm>

#include <hrleSparseStarIterator.hpp>

#include <lsDomain.hpp>
#include <lsExpand.hpp>

#include <vcLogger.hpp>
#include <vcSmartPointer.hpp>
#include <vcVectorType.hpp>

namespace viennals {

using namespace viennacore;

enum class lsNormalCalculationMethodEnum {
  CENTRAL_DIFFERENCES,
  ONE_SIDED_MIN_MOD
};

/// This algorithm is used to compute the normal vectors for all points
/// with level set values <= maxValue (default 0.5). The result is saved in
/// the lsPointData of the lsDomain and can be retrieved with
/// lsDomain.getPointData().getVectorData("Normals").
///
/// The algorithm uses central differences to compute gradients and normalizes
/// them to unit vectors. Since neighbors in each cartesian direction are
/// necessary for the calculation, the level set width must be >= (maxValue * 4)
/// + 1. If the level set width is insufficient, it will be automatically
/// expanded.
///
/// The calculation is parallelized using OpenMP across level set segments.
/// Normal vectors are computed using finite differences and normalized to unit
/// length. Points with zero gradient magnitude are assigned zero normal
/// vectors.
template <class T, int D> class CalculateNormalVectors {
  SmartPointer<Domain<T, D>> levelSet = nullptr;
  T maxValue = 0.5;
  lsNormalCalculationMethodEnum method =
      lsNormalCalculationMethodEnum::CENTRAL_DIFFERENCES;

  // Constants for calculation
  static constexpr T DEFAULT_MAX_VALUE = 0.5;
  static constexpr T EPSILON = 1e-12;
  static constexpr T FINITE_DIFF_FACTOR = 0.5;

public:
  static constexpr char normalVectorsLabel[] = "Normals";

  CalculateNormalVectors() = default;

  CalculateNormalVectors(
      SmartPointer<Domain<T, D>> passedLevelSet,
      T passedMaxValue = DEFAULT_MAX_VALUE,
      lsNormalCalculationMethodEnum passedMethod =
          lsNormalCalculationMethodEnum::CENTRAL_DIFFERENCES)
      : levelSet(passedLevelSet), maxValue(passedMaxValue), method(passedMethod) {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  void setMaxValue(const T passedMaxValue) {
    if (passedMaxValue <= 0) {
      VIENNACORE_LOG_WARNING(
          "CalculateNormalVectors: maxValue should be positive. "
          "Using default value " +
          std::to_string(DEFAULT_MAX_VALUE) + ".");
      maxValue = DEFAULT_MAX_VALUE;
    } else {
      maxValue = passedMaxValue;
    }
  }

  void setMethod(lsNormalCalculationMethodEnum passedMethod) {
    method = passedMethod;
  }

  SmartPointer<Domain<T, D>> getLevelSet() const { return levelSet; }

  T getMaxValue() const { return maxValue; }

  /// Check if normal vectors are already calculated for the level set
  bool hasNormalVectors() const {
    if (levelSet == nullptr)
      return false;
    auto &pointData = levelSet->getPointData();
    return pointData.getVectorData(normalVectorsLabel) != nullptr;
  }

  void apply() {
    if (levelSet == nullptr) {
      VIENNACORE_LOG_ERROR(
          "No level set was passed to CalculateNormalVectors.");
      return;
    }

    switch (method) {
    case lsNormalCalculationMethodEnum::CENTRAL_DIFFERENCES:
      calculateCentralDifferences();
      break;
    case lsNormalCalculationMethodEnum::ONE_SIDED_MIN_MOD:
      calculateOneSidedMinMod();
      break;
    }
  }

private:
  void calculateCentralDifferences() {
    if (levelSet->getLevelSetWidth() < (maxValue * 4) + 1) {
      VIENNACORE_LOG_WARNING("CalculateNormalVectors: Level set width must be "
                             "greater than " +
                             std::to_string((maxValue * 4) + 1) +
                             ". Expanding level set to " +
                             std::to_string((maxValue * 4) + 1) + ".");
      Expand<T, D>(levelSet, (maxValue * 4) + 1).apply();
    }

    std::vector<std::vector<Vec3D<T>>> normalVectorsVector(
        levelSet->getNumberOfSegments());

    // Estimate memory requirements per thread to improve cache performance
    double pointsPerSegment =
        double(2 * levelSet->getDomain().getNumberOfPoints()) /
        double(levelSet->getLevelSetWidth());

    auto grid = levelSet->getGrid();

    // Calculate Normalvectors
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
          Vec3D<T> tmp{};
          normalVectors.push_back(tmp);
          continue;
        }

        Vec3D<T> n{};

        T denominator = 0;
        for (int i = 0; i < D; i++) {
          viennahrle::Index<D> posIdx(0);
          posIdx[i] = 1;
          viennahrle::Index<D> negIdx(0);
          negIdx[i] = -1;
          T pos = neighborIt.getNeighbor(posIdx).getValue() - center.getValue();
          T neg = center.getValue() - neighborIt.getNeighbor(negIdx).getValue();
          n[i] = (pos + neg) * FINITE_DIFF_FACTOR;
          denominator += n[i] * n[i];
        }

        denominator = std::sqrt(denominator);
        if (std::abs(denominator) < EPSILON) {
          VIENNACORE_LOG_WARNING(
              "CalculateNormalVectors: Vector of length 0 at " +
              neighborIt.getIndices().to_string());
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
    insertIntoPointData(normalVectorsVector);
  }

  void calculateOneSidedMinMod() {
    // This method does not require expansion, since it is robust to missing
    // neighbors.
    auto &domain = levelSet->getDomain();
    auto &grid = levelSet->getGrid();

    // Directly write to a single vector indexed by point ID.
    // This is thread-safe since each thread works on a distinct range of
    // points.
    std::vector<Vec3D<T>> normalVectors(levelSet->getNumberOfPoints());

#pragma omp parallel num_threads(levelSet->getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      viennahrle::Index<D> startVector =
          (p == 0) ? grid.getMinGridPoint()
                   : levelSet->getDomain().getSegmentation()[p - 1];
      viennahrle::Index<D> endVector =
          (p != static_cast<int>(levelSet->getNumberOfSegments() - 1))
              ? levelSet->getDomain().getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      viennahrle::SparseStarIterator<typename Domain<T, D>::DomainType, 1>
          neighborIt(domain, startVector);

      for (; neighborIt.getIndices() < endVector; neighborIt.next()) {
        if (neighborIt.getCenter().isDefined()) {
          Vec3D<T> grad{};
          for (int i = 0; i < D; ++i) {
            viennahrle::Index<D> posIdx(0);
            posIdx[i] = 1;
            viennahrle::Index<D> negIdx(0);
            negIdx[i] = -1;
            bool negDefined = neighborIt.getNeighbor(negIdx).isDefined();
            bool posDefined = neighborIt.getNeighbor(posIdx).isDefined();

            if (negDefined && posDefined) {
              T valNeg = neighborIt.getNeighbor(negIdx).getValue();
              T valCenter = neighborIt.getCenter().getValue();
              T valPos = neighborIt.getNeighbor(posIdx).getValue();

              const bool centerSign = valCenter >= 0;
              const bool negSign = valNeg > 0;
              const bool posSign = valPos > 0;

              const T d_neg = valCenter - valNeg;
              const T d_pos = valPos - valCenter;

              if (centerSign != negSign && centerSign != posSign) {
                // Center is an extremum, use minmod to be safe
                grad[i] = 0.;
                    // (std::abs(d_pos) < std::abs(d_neg)) ? d_pos : d_neg;
              } else if (centerSign != negSign) {
                // Interface is on the negative side, use backward difference
                grad[i] = d_neg;
              } else if (centerSign != posSign) {
                // Interface is on the positive side, use forward difference
                grad[i] = d_pos;
              } else {
                // No sign change, use minmod to handle sharp features smoothly
                grad[i] = 
                    (std::abs(d_pos) < std::abs(d_neg)) ? d_pos : d_neg;
              }
            } else if (negDefined) {
              grad[i] = (neighborIt.getCenter().getValue() -
                         neighborIt.getNeighbor(negIdx).getValue());
            } else if (posDefined) {
              grad[i] = (neighborIt.getNeighbor(posIdx).getValue() -
                         neighborIt.getCenter().getValue());
            } else {
              grad[i] = 0;
            }
          }
          Normalize(grad);
          normalVectors[neighborIt.getCenter().getPointId()] = grad;
        }
      }
    }

    // insert into pointData of levelSet
    auto &pointData = levelSet->getPointData();
    auto vectorDataPointer = pointData.getVectorData(normalVectorsLabel, true);
    // if it does not exist, insert new normals vector
    if (vectorDataPointer == nullptr) {
      pointData.insertNextVectorData(normalVectors, normalVectorsLabel);
    } else {
      // if it does exist, just swap the old with the new values
      *vectorDataPointer = std::move(normalVectors);
    }
  }

  void insertIntoPointData(
      std::vector<std::vector<Vec3D<T>>> &normalVectorsVector) {
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
