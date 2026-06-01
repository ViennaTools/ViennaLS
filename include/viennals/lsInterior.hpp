#pragma once

#include <lsPreCompileMacros.hpp>

#include <hrleSparseStarIterator.hpp>
#include <lsDomain.hpp>

#include <vcVectorType.hpp>

namespace viennals {

using namespace viennacore;

/// Activates interior points in a level set.
/// The largest value in the levelset is thus width*0.5
/// Returns the number of added points
template <class T, int D> class Interior {
  SmartPointer<Domain<T, D>> levelSet = nullptr;
  int width = 0;
  bool updatePointData = true;

public:
  Interior() = default;

  Interior(SmartPointer<Domain<T, D>> passedlsDomain)
      : levelSet(passedlsDomain) {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  /// Set whether to update the point data stored in the LS
  /// during this algorithm. Defaults to true.
  void setUpdatePointData(bool update) { updatePointData = update; }

  /// Apply the interior definition until no more points can be added
  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addError("No level set passed to Interior. Not activating.")
          .print();
      return;
    }

    if (levelSet->getNumberOfPoints() == 0)
      return;

    // Store initial domain extents by iterating through all domain points
    auto &domain = levelSet->getDomain();
    auto &grid = levelSet->getGrid();

    viennahrle::Index<D> initialMinExtent = grid.getMaxGridPoint();
    viennahrle::Index<D> initialMaxExtent = grid.getMinGridPoint();

    // Use ConstSparseStarIterator to find all points and track extents
    viennahrle::ConstSparseStarIterator<typename Domain<T, D>::DomainType, 1>
        extentIt(domain, grid.getMinGridPoint());
    viennahrle::Index<D> extentEnd = grid.incrementIndices(grid.getMaxGridPoint());

    bool foundPoints = false;
    while (extentIt.getIndices() < extentEnd) {
      if (extentIt.getCenter().isDefined()) {
        viennahrle::Index<D> idx = extentIt.getIndices();
        foundPoints = true;
        for (int d = 0; d < D; ++d) {
          if (idx[d] < initialMinExtent[d])
            initialMinExtent[d] = idx[d];
          if (idx[d] > initialMaxExtent[d])
            initialMaxExtent[d] = idx[d];
        }
      }
      extentIt.next();
    }

    // If no points found, use grid bounds as fallback
    if (!foundPoints) {
      initialMinExtent = grid.getMinGridPoint();
      initialMaxExtent = grid.getMaxGridPoint();
    }

    const int startWidth = levelSet->getLevelSetWidth();
    const T startLimit = startWidth * 0.5;

    for (int currentCycle = 0; currentCycle < 100;
         ++currentCycle) {
      int addedPoints = 0;

      const T limit = (startWidth + currentCycle + 1);

      auto &grid = levelSet->getGrid();
      auto newlsDomain = SmartPointer<Domain<T, D>>::New(grid);
      auto &newDomain = newlsDomain->getDomain();
      auto &domain = levelSet->getDomain();

      newDomain.initialize(domain.getNewSegmentation(),
                           domain.getAllocation());

      const bool updateData = updatePointData;
      // save how data should be transferred to new level set
      // list of indices into the old pointData vector
      std::vector<std::vector<unsigned>> newDataSourceIds;
      if (updateData)
        newDataSourceIds.resize(newDomain.getNumberOfSegments());

#pragma omp parallel num_threads(newDomain.getNumberOfSegments()) reduction(+:addedPoints)
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

        for (viennahrle::ConstSparseStarIterator<
                 typename Domain<T, D>::DomainType, 1>
                 neighborIt(domain, startVector);
             neighborIt.getIndices() < endVector; neighborIt.next()) {

          auto &centerIt = neighborIt.getCenter();
          viennahrle::Index<D> idx = neighborIt.getIndices();

          // Check if point is within initial domain extents (all axes, both directions)
          bool withinInitialExtents = true;
          for (int d = 0; d < D; ++d) {
            // Check extents
            if (idx[d] < initialMinExtent[d] || idx[d] > initialMaxExtent[d]) {
              withinInitialExtents = false;
              break;
            }
          }

          if (centerIt.getValue() == Domain<T, D>::NEG_VALUE) {
              // Interior/negative side. Track the best *defined* neighbor
              // separately: undefined sentinels (POS_VALUE / NEG_VALUE) must
              // not be used as pointData sources — getPointId() on an
              // undefined run returns an invalid index.
              T distance = Domain<T, D>::NEG_VALUE;
              T definedDistance = Domain<T, D>::NEG_VALUE;
              int definedNeighbor = -1;
              for (int i = 0; i < 2 * D; i++) {
                auto &nb = neighborIt.getNeighbor(i);
                T newValue = nb.getValue() - T(1);
                if (distance < newValue)
                  distance = newValue;
                if (nb.isDefined() && definedDistance < newValue) {
                  definedDistance = newValue;
                  definedNeighbor = i;
                }
              }
              // Only activate when a proper defined neighbor drives the fill.
              if (definedDistance >= -limit && withinInitialExtents &&
                  definedNeighbor != -1) {
                addedPoints++;
                domainSegment.insertNextDefinedPoint(neighborIt.getIndices(),
                                                     definedDistance);
                if (updateData)
                  newDataSourceIds[p].push_back(
                      neighborIt.getNeighbor(definedNeighbor).getPointId());
              } else {
                // insertNextUndefinedRunType
                domainSegment.insertNextUndefinedPoint(neighborIt.getIndices(),
                                                       Domain<T, D>::NEG_VALUE);
              }
          } else if (centerIt.getValue() == Domain<T, D>::POS_VALUE) {
              // Ignore positive phi values (exterior region)
              domainSegment.insertNextUndefinedPoint(neighborIt.getIndices(),
                                                     Domain<T, D>::POS_VALUE);
          } else {
            domainSegment.insertNextDefinedPoint(neighborIt.getIndices(),
                                                 centerIt.getValue());
            if (updateData)
              newDataSourceIds[p].push_back(centerIt.getPointId());
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
      if (addedPoints == 0)
        break;
    }
    levelSet->getDomain().segment();
    // levelSet->finalize(width);
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(Interior)

} // namespace viennals
