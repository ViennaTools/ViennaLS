#ifndef LS_CALCULATE_NORMAL_VECTORS_HPP
#define LS_CALCULATE_NORMAL_VECTORS_HPP

#include <lsPreCompileMacros.hpp>

#include <algorithm>

#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsDomain.hpp>
#include <lsMessage.hpp>

/// This algorithm is used to compute the normal vectors for all points
/// defined in the level set. The result is saved in the lsDomain and
/// can be retrieved with lsDomain.getNormalVectors().
/// Since neighbors in each cartesian direction are necessary for
/// the calculation, the levelset width must be >=3.
template <class T, int D> class lsCalculateNormalVectors {
  lsDomain<T, D> *domain = nullptr;
  bool onlyActivePoints = false;

public:
  lsCalculateNormalVectors() {}

  lsCalculateNormalVectors(lsDomain<T, D> &passedDomain,
                           bool passedOnlyActivePoints = false)
      : domain(&passedDomain), onlyActivePoints(passedOnlyActivePoints) {}

  void setLevelSet(lsDomain<T, D> &passedDomain) { domain = &passedDomain; }

  /// Set whether normal vectors should only be calculated for level set
  /// points <=0.5. Defaults to false.
  void setOnlyActivePoints(bool passedOnlyActivePoints) {
    onlyActivePoints = passedOnlyActivePoints;
  }

  void apply() {
    if (domain == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsCalculateNormalVectors.")
          .print();
    }

    if (domain->getLevelSetWidth() < 3) {
      lsMessage::getInstance()
          .addWarning("lsCalculateNormalVectors: Level set width must be "
                      "greater than 2!")
          .print();
    }

    std::vector<std::vector<hrleVectorType<T, D>>> normalVectorsVector(
        domain->getNumberOfSegments());
    double pointsPerSegment =
        double(2 * domain->getDomain().getNumberOfPoints()) /
        double(domain->getLevelSetWidth() *
               domain->getDomain().getNumberOfSegments());

    auto grid = domain->getGrid();

    //! Calculate Normalvectors
#pragma omp parallel num_threads(domain->getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      std::vector<hrleVectorType<T, D>> &normalVectors = normalVectorsVector[p];
      normalVectors.reserve(pointsPerSegment);

      hrleVectorType<hrleIndexType, D> startVector =
          (p == 0) ? grid.getMinGridPoint()
                   : domain->getDomain().getSegmentation()[p - 1];

      hrleVectorType<hrleIndexType, D> endVector =
          (p != static_cast<int>(domain->getNumberOfSegments() - 1))
              ? domain->getDomain().getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      for (hrleConstSparseStarIterator<typename lsDomain<T, D>::DomainType>
               neighborIt(domain->getDomain(), startVector);
           neighborIt.getIndices() < endVector; neighborIt.next()) {

        auto &center = neighborIt.getCenter();
        if (!center.isDefined() ||
            (onlyActivePoints && std::abs(center.getValue()) > 0.5))
          continue;

        hrleVectorType<T, D> n;

        T denominator = 0;
        for (int i = 0; i < D; i++) {
          T pos = neighborIt.getNeighbor(i).getValue() - center.getValue();
          T neg = center.getValue() - neighborIt.getNeighbor(i + D).getValue();
          n[i] = (pos + neg) * 0.5;
          denominator += n[i] * n[i];
        }

        denominator = std::sqrt(denominator);
        for (unsigned i = 0; i < D; ++i) {
          n[i] /= denominator;
        }

        normalVectors.push_back(n);
      }
    }

    // copy all normals
    auto &normals = domain->getNormalVectors();
    normals.clear();
    unsigned numberOfNormals = 0;
    for (unsigned i = 0; i < domain->getNumberOfSegments(); ++i) {
      numberOfNormals += normalVectorsVector[i].size();
    }
    normals.reserve(numberOfNormals);

    for (unsigned i = 0; i < domain->getNumberOfSegments(); ++i) {
      normals.insert(normals.end(), normalVectorsVector[i].begin(),
                     normalVectorsVector[i].end());
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsCalculateNormalVectors)

#endif // LS_CALCULATE_NORMAL_VECTORS_HPP
