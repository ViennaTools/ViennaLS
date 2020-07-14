#ifndef LS_REDUCE_HPP
#define LS_REDUCE_HPP

#include <lsPreCompileMacros.hpp>

#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>
#include <lsDomain.hpp>

/// Reduce the level set size to the specified width.
/// This means all level set points with value <= 0.5*width
/// are removed, reducing the memory footprint of the lsDomain.
template <class T, int D> class lsReduce {
  typedef typename lsDomain<T, D>::GridType GridType;
  lsSmartPointer<lsDomain<T, D>> levelSet = nullptr;
  int width = 0;
  bool noNewSegment = false;

public:
  lsReduce() {}

  lsReduce(lsSmartPointer<lsDomain<T, D>> passedlsDomain)
      : levelSet(passedlsDomain){};

  lsReduce(lsSmartPointer<lsDomain<T, D>> passedlsDomain, int passedWidth,
           bool passedNoNewSegment = false)
      : levelSet(passedlsDomain), width(passedWidth),
        noNewSegment(passedNoNewSegment){};

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  /// Set which level set points should be kept.
  /// All points with a level set value >0.5*width will be
  /// removed by this algorithm.
  void setWidth(int passedWidth) { width = passedWidth; }

  /// Set whether to segment the level set after algorithm
  /// is finished. This means points will be evenly distributed
  /// across points. Defaults to true.
  void setNoNewSegment(bool passedNoNewSegment) {
    noNewSegment = passedNoNewSegment;
  }

  /// Reduces the leveleSet to the specified number of layers.
  /// The largest value in the levelset is thus width*0.5
  /// Returns the number of added points
  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsReduce.")
          .print();
      return;
    }

    if (width >= levelSet->getLevelSetWidth())
      return;

    const T valueLimit = width * 0.5;

    auto &grid = levelSet->getGrid();
    auto newlsDomain = lsSmartPointer<lsDomain<T, D>>::New(levelSet->getGrid());
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
    levelSet->deepCopy(newlsDomain);
    levelSet->finalize(width);
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsReduce)

#endif // LS_REDUCE_HPP
