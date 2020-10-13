#ifndef LS_PRUNE_HPP
#define LS_PRUNE_HPP

#include <lsPreCompileMacros.hpp>

#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsDomain.hpp>

/// Removes all level set points, which do not have
/// at least one oppositely signed neighbour (Meaning
/// they do not lie directly at the interface).
/// Afterwards the level set will occupy the least memory
/// possible.
template <class T, int D> class lsPrune {
  lsSmartPointer<lsDomain<T, D>> levelSet = nullptr;

public:
  lsPrune() {}

  lsPrune(lsSmartPointer<lsDomain<T, D>> passedlsDomain)
      : levelSet(passedlsDomain){};

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  /// removes all grid points, which do not have at least one opposite signed
  /// neighbour
  /// returns the number of removed points
  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsPrune.")
          .print();
      return;
    }
    if(levelSet->getNumberOfPoints() == 0) {
      return;
    }

    auto &grid = levelSet->getGrid();
    auto newlsDomain = lsSmartPointer<lsDomain<T, D>>::New(grid);
    typename lsDomain<T, D>::DomainType &newDomain = newlsDomain->getDomain();
    typename lsDomain<T, D>::DomainType &domain = levelSet->getDomain();

    newDomain.initialize(domain.getNewSegmentation(), domain.getAllocation());

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
        bool centerSign = std::signbit(centerIt.getValue());
        if (centerIt.isDefined()) {
          int i = 0;
          for (; i < 2 * D; i++) {
            // Use signbit here instead of numericEps because it makes a clearer cut
            // between negative and positive numbers. Eps resulted in problems with
            // exact 0.0 LS values.
            if(std::signbit(neighborIt.getNeighbor(i).getValue()) != centerSign)
              break;
          }
          if (i != 2 * D) {
            domainSegment.insertNextDefinedPoint(neighborIt.getIndices(),
                                                 centerIt.getValue());
          } else {
            // TODO: it is more efficient to insertNextUndefinedRunType, since
            // we know it already exists
            domainSegment.insertNextUndefinedPoint(
                neighborIt.getIndices(), centerSign
                                             ? lsDomain<T, D>::NEG_VALUE
                                             : lsDomain<T, D>::POS_VALUE);
          }
        } else {
          domainSegment.insertNextUndefinedPoint(
              neighborIt.getIndices(), centerSign ? lsDomain<T, D>::NEG_VALUE
                                                  : lsDomain<T, D>::POS_VALUE);
        }
      }
    }
    // distribute evenly across segments and copy
    newDomain.finalize();
    newDomain.segment();
    levelSet->deepCopy(newlsDomain);
    levelSet->finalize(2);
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsPrune)

#endif // LS_PRUNE_HPP
