#pragma once

#include <lsDomain.hpp>
#include <lsPreCompileMacros.hpp>
#include <vcVectorType.hpp>

namespace viennals {

using namespace viennacore;

/// Reduce the level set size to the specified width.
/// This means all level set points with value <= 0.5*width
/// are removed, reducing the memory footprint of the lsDomain.
template <class T, int D> class Reduce {
  typedef typename Domain<T, D>::GridType GridType;
  SmartPointer<Domain<T, D>> levelSet = nullptr;
  int width = 0;
  bool noNewSegment = false;
  bool updatePointData = true;

public:
  Reduce() = default;

  Reduce(SmartPointer<Domain<T, D>> passedlsDomain)
      : levelSet(passedlsDomain) {};

  Reduce(SmartPointer<Domain<T, D>> passedlsDomain, int passedWidth,
         bool passedNoNewSegment = false)
      : levelSet(passedlsDomain), width(passedWidth),
        noNewSegment(passedNoNewSegment) {};

  void setLevelSet(SmartPointer<Domain<T, D>> passedlsDomain) {
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

  /// Set whether to update the point data stored in the LS
  /// during this algorithm. Defaults to true.
  void setUpdatePointData(bool update) { updatePointData = update; }

  /// Reduces the leveleSet to the specified number of layers.
  /// The largest value in the levelset is thus width*0.5
  /// Returns the number of added points
  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addWarning("No level set was passed to Reduce.")
          .print();
      return;
    }

    if (width >= levelSet->getLevelSetWidth())
      return;

    const T valueLimit = width * 0.5;

    auto &grid = levelSet->getGrid();
    auto newlsDomain = SmartPointer<Domain<T, D>>::New(levelSet->getGrid());
    typename Domain<T, D>::DomainType &newDomain = newlsDomain->getDomain();
    typename Domain<T, D>::DomainType &domain = levelSet->getDomain();

    if (noNewSegment)
      newDomain.initialize(domain.getSegmentation(), domain.getAllocation());
    else
      newDomain.initialize(domain.getNewSegmentation(), domain.getAllocation());

    const bool updateData = updatePointData;
    // save how data should be transferred to new level set
    // list of indices into the old pointData vector
    std::vector<std::vector<unsigned>> newDataSourceIds;
    if (updateData)
      newDataSourceIds.resize(newDomain.getNumberOfSegments());

#pragma omp parallel num_threads(newDomain.getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      auto &domainSegment = newDomain.getDomainSegment(p);

      viennahrle::Index<D> const startVector =
          (p == 0) ? grid.getMinGridPoint()
                   : newDomain.getSegmentation()[p - 1];

      viennahrle::Index<D> const endVector =
          (p != static_cast<int>(newDomain.getNumberOfSegments() - 1))
              ? newDomain.getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      for (viennahrle::SparseIterator<typename Domain<T, D>::DomainType> it(
               domain, startVector);
           it.getStartIndices() < endVector; ++it) {
        T currentValue = it.getValue();
        if (it.isDefined() && std::abs(currentValue) <= valueLimit) {
          domainSegment.insertNextDefinedPoint(it.getStartIndices(),
                                               currentValue);
          if (updateData)
            newDataSourceIds[p].push_back(it.getPointId());
        } else {
          // TODO: use insertNextUndefinedRunType
          domainSegment.insertNextUndefinedPoint(it.getStartIndices(),
                                                 (currentValue < 0)
                                                     ? Domain<T, D>::NEG_VALUE
                                                     : Domain<T, D>::POS_VALUE);
        }
      }
    }

    // now copy old data into new level set
    if (updateData) {
      newlsDomain->getPointData().translateFromMultiData(
          levelSet->getPointData(), newDataSourceIds);
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
PRECOMPILE_PRECISION_DIMENSION(Reduce)

} // namespace viennals
