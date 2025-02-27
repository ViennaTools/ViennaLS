#pragma once

#include <lsPreCompileMacros.hpp>

#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>
#include <lsDomain.hpp>

namespace viennals {

using namespace viennacore;

/// Expands the levelSet to the specified number of layers.
/// The largest value in the levelset is thus width*0.5
/// Returns the number of added points
template <class T, int D> class Expand {
  typedef typename Domain<T, D>::GridType GridType;
  SmartPointer<Domain<T, D>> levelSet = nullptr;
  int width = 0;
  bool updatePointData = true;

public:
  Expand() = default;

  Expand(SmartPointer<Domain<T, D>> passedlsDomain)
      : levelSet(passedlsDomain) {}

  Expand(SmartPointer<Domain<T, D>> passedlsDomain, int passedWidth)
      : levelSet(passedlsDomain), width(passedWidth) {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  /// Set how far the level set should be extended. Points
  /// with value width*0.5 will be added by this algorithm.
  void setWidth(int passedWidth) { width = passedWidth; }

  /// Set whether to update the point data stored in the LS
  /// during this algorithm. Defaults to true.
  void setUpdatePointData(bool update) { updatePointData = update; }

  /// Apply the expansion to the specified width
  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addWarning("No level set passed to Expand. Not expanding.")
          .print();
      return;
    }

    if (width <= levelSet->getLevelSetWidth())
      return;

    if (levelSet->getNumberOfPoints() == 0)
      return;

    const T totalLimit = width * 0.5;
    const int startWidth = levelSet->getLevelSetWidth();
    const int numberOfRequiredCycles = width - startWidth;

    for (int currentCycle = 0; currentCycle < numberOfRequiredCycles;
         ++currentCycle) {

      const int allocationFactor =
          1 + 1.0 / static_cast<double>(startWidth + currentCycle);
      const T limit = (startWidth + currentCycle + 1) * T(0.5);

      auto &grid = levelSet->getGrid();
      auto newlsDomain = SmartPointer<Domain<T, D>>::New(grid);
      typename Domain<T, D>::DomainType &newDomain = newlsDomain->getDomain();
      typename Domain<T, D>::DomainType &domain = levelSet->getDomain();

      newDomain.initialize(domain.getNewSegmentation(),
                           domain.getAllocation() * allocationFactor);

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

        hrleVectorType<hrleIndexType, D> startVector =
            (p == 0) ? grid.getMinGridPoint()
                     : newDomain.getSegmentation()[p - 1];

        hrleVectorType<hrleIndexType, D> endVector =
            (p != static_cast<int>(newDomain.getNumberOfSegments() - 1))
                ? newDomain.getSegmentation()[p]
                : grid.incrementIndices(grid.getMaxGridPoint());

        for (hrleSparseStarIterator<typename Domain<T, D>::DomainType, 1>
                 neighborIt(domain, startVector);
             neighborIt.getIndices() < endVector; neighborIt.next()) {

          auto &centerIt = neighborIt.getCenter();
          if (std::abs(centerIt.getValue()) <= totalLimit) {
            domainSegment.insertNextDefinedPoint(neighborIt.getIndices(),
                                                 centerIt.getValue());
            if (updateData)
              newDataSourceIds[p].push_back(centerIt.getPointId());
          } else {
            if (centerIt.getValue() > -std::numeric_limits<T>::epsilon()) {
              T distance = Domain<T, D>::POS_VALUE;
              int neighbor = -1;
              for (int i = 0; i < 2 * D; i++) {
                T newValue = neighborIt.getNeighbor(i).getValue() + T(1);
                if (distance > newValue) {
                  distance = newValue;
                  neighbor = i;
                }
              }
              if (distance <= limit) {
                domainSegment.insertNextDefinedPoint(neighborIt.getIndices(),
                                                     distance);
                if (updateData)
                  newDataSourceIds[p].push_back(
                      neighborIt.getNeighbor(neighbor).getPointId());
              } else {
                // TODO: use insertNextUndefinedRunType
                domainSegment.insertNextUndefinedPoint(neighborIt.getIndices(),
                                                       Domain<T, D>::POS_VALUE);
              }
            } else {
              T distance = Domain<T, D>::NEG_VALUE;
              int neighbor = -1;
              for (int i = 0; i < 2 * D; i++) {
                T newValue = neighborIt.getNeighbor(i).getValue() - T(1);
                if (distance < newValue) {
                  distance = newValue;
                  neighbor = i;
                }
              }
              if (distance >= -limit) {
                domainSegment.insertNextDefinedPoint(neighborIt.getIndices(),
                                                     distance);
                if (updateData)
                  newDataSourceIds[p].push_back(
                      neighborIt.getNeighbor(neighbor).getPointId());
              } else {
                // TODO: use insertNextUndefinedRunType
                domainSegment.insertNextUndefinedPoint(neighborIt.getIndices(),
                                                       Domain<T, D>::NEG_VALUE);
              }
            }
          }
        }
      }

      // now copy old data into new level set
      if (updateData) {
        newlsDomain->getPointData().translateFromMultiData(
            levelSet->getPointData(), newDataSourceIds);
      }

      newDomain.finalize();
      levelSet->deepCopy(newlsDomain);
    }
    levelSet->getDomain().segment();
    levelSet->finalize(width);
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(Expand)

} // namespace viennals
