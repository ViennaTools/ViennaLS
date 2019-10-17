#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>
#include <lsDomain_template.hpp>

template <class T, int D> class lsReduce {
  typedef typename lsDomain<T, D>::GridType GridType;
  lsDomain<T, D> &levelSet;

public:
  lsReduce(lsDomain<T, D> &passedlsDomain) : levelSet(passedlsDomain){};

  /// Reduces the leveleSet to the specified number of layers.
  /// The largest value in the levelset is thus width*0.5
  /// Returns the number of added points
  void apply(int width, bool noNewSegment = false) {
    if (width >= levelSet.getLevelSetWidth())
      return;

    const T valueLimit = width * 0.5;

    auto &grid = levelSet.getGrid();
    lsDomain<T, D> newlsDomain(levelSet.getGrid());
    typename lsDomain<T, D>::DomainType &newDomain = newlsDomain.getDomain();
    typename lsDomain<T, D>::DomainType &domain = levelSet.getDomain();

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

      for (hrleSparseIterator<typename lsDomain<T, D>::DomainType> it(
               domain, startVector);
           it.getStartIndices() < endVector; ++it) {
        T currentValue = it.getValue();
        if (it.isDefined() && std::abs(currentValue) <= valueLimit) {
          domainSegment.insertNextDefinedPoint(it.getStartIndices(),
                                               currentValue);
        } else {
          // TODO: use insertNextUndefinedRunType
          domainSegment.insertNextUndefinedPoint(
              it.getStartIndices(), (currentValue < 0)
                                        ? lsDomain<T, D>::NEG_VALUE
                                        : lsDomain<T, D>::POS_VALUE);
        }
      }
    }

    // distribute evenly across segments and copy
    newDomain.finalize();
    if (!noNewSegment)
      newDomain.segment();
    levelSet.deepCopy(newlsDomain);
    levelSet.finalize(width);
  }
};
