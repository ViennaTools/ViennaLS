#ifndef LS_EXPAND_HPP
#define LS_EXPAND_HPP

#include <lsPreCompileMacros.hpp>

#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>
#include <lsDomain.hpp>
#include <lsMessage.hpp>

/// Expands the leveleSet to the specified number of layers.
/// The largest value in the levelset is thus width*0.5
/// Returns the number of added points
template <class T, int D> class lsExpand {
  typedef typename lsDomain<T, D>::GridType GridType;
  lsSmartPointer<lsDomain<T, D>> levelSet = nullptr;
  int width = 0;

public:
  lsExpand() {}

  lsExpand(lsSmartPointer<lsDomain<T, D>> passedlsDomain)
      : levelSet(passedlsDomain) {}

  lsExpand(lsSmartPointer<lsDomain<T, D>> passedlsDomain, int passedWidth)
      : levelSet(passedlsDomain), width(passedWidth) {}

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  /// Set how far the level set should be extended. Points
  /// with value width*0.5 will be added by this algorithm.
  void setWidth(int passedWidth) { width = passedWidth; }

  /// Apply the expansion to the specified width
  void apply() {
    if (width <= levelSet->getLevelSetWidth())
      return;

    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set passed to lsExpand. Not expanding.")
          .print();
    }

    if (levelSet->getNumberOfPoints() == 0) {
      return;
    }

    const T totalLimit = width * 0.5;
    const int startWidth = levelSet->getLevelSetWidth();
    const int numberOfRequiredCycles = width - startWidth;

    for (int currentCycle = 0; currentCycle < numberOfRequiredCycles;
         ++currentCycle) {

      const int allocationFactor =
          1 + 1.0 / static_cast<double>(startWidth + currentCycle);
      const T limit = (startWidth + currentCycle + 1) * T(0.5);

      auto &grid = levelSet->getGrid();
      auto newlsDomain = lsSmartPointer<lsDomain<T, D>>::New(grid);
      typename lsDomain<T, D>::DomainType &newDomain = newlsDomain->getDomain();
      typename lsDomain<T, D>::DomainType &domain = levelSet->getDomain();

      newDomain.initialize(domain.getNewSegmentation(),
                           domain.getAllocation() * allocationFactor);

#pragma omp parallel num_threads(newDomain.getNumberOfSegments())
      {
        int p = 0;
#ifdef _OPENMP
        p = omp_get_thread_num();
#endif

        auto &domainSegment = newDomain.getDomainSegment(p);

        hrleVectorType<hrleIndexType, D> startVector =
            (p == 0) ? grid.getMinGridPoint()
                     : newDomain.getSegmentation()[p - 1];

        hrleVectorType<hrleIndexType, D> endVector =
            (p != static_cast<int>(newDomain.getNumberOfSegments() - 1))
                ? newDomain.getSegmentation()[p]
                : grid.incrementIndices(grid.getMaxGridPoint());

        for (hrleSparseStarIterator<typename lsDomain<T, D>::DomainType>
                 neighborIt(domain, startVector);
             neighborIt.getIndices() < endVector; neighborIt.next()) {

          auto &centerIt = neighborIt.getCenter();
          if (std::abs(centerIt.getValue()) <= totalLimit) {
            domainSegment.insertNextDefinedPoint(neighborIt.getIndices(),
                                                 centerIt.getValue());
          } else {
            if (centerIt.getValue() > 0.) {
              T distance = lsDomain<T, D>::POS_VALUE;
              for (int i = 0; i < 2 * D; i++) {
                distance = std::min(
                    distance, neighborIt.getNeighbor(i).getValue() + T(1));
              }
              if (distance <= limit) {
                domainSegment.insertNextDefinedPoint(neighborIt.getIndices(),
                                                     distance);
              } else {
                // TODO: use insertNextUndefinedRunType
                domainSegment.insertNextUndefinedPoint(
                    neighborIt.getIndices(), lsDomain<T, D>::POS_VALUE);
              }
            } else {
              T distance = lsDomain<T, D>::NEG_VALUE;
              for (int i = 0; i < 2 * D; i++) {
                distance = std::max(
                    distance, neighborIt.getNeighbor(i).getValue() - T(1));
              }
              if (distance >= -limit) {
                domainSegment.insertNextDefinedPoint(neighborIt.getIndices(),
                                                     distance);
              } else {
                // TODO: use insertNextUndefinedRunType
                domainSegment.insertNextUndefinedPoint(
                    neighborIt.getIndices(), lsDomain<T, D>::NEG_VALUE);
              }
            }
          }
        }
      }
      newDomain.finalize();
      levelSet->deepCopy(newlsDomain);
    }
    levelSet->getDomain().segment();
    levelSet->finalize(width);
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsExpand)

#endif // LS_EXPAND_HPP
