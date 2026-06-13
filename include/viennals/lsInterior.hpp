#pragma once

#include <lsPreCompileMacros.hpp>

#include <hrleSparseIterator.hpp>
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
  SmartPointer<Domain<T, D>> guide = nullptr;
  int width = 0;
  bool updatePointData = true;

public:
  Interior() = default;

  Interior(SmartPointer<Domain<T, D>> passedlsDomain)
      : levelSet(passedlsDomain) {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  /// Set a guide level set: the fill will not activate points that are inside
  /// the guide (φ_guide < 0).  Pass the Si surface to stop the oxide Interior
  /// fill at the Si-SiO2 interface.
  void setGuide(SmartPointer<Domain<T, D>> g) { guide = g; }

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

    const int startWidth = levelSet->getLevelSetWidth();

    for (int currentCycle = 0; currentCycle < 100; ++currentCycle) {
      int addedPoints = 0;

      const T limit = (startWidth + currentCycle + 1);

      auto &grid = levelSet->getGrid();
      auto newlsDomain = SmartPointer<Domain<T, D>>::New(grid);
      auto &newDomain = newlsDomain->getDomain();
      auto &domain = levelSet->getDomain();

      newDomain.initialize(domain.getNewSegmentation(), domain.getAllocation());

      const bool updateData = updatePointData;
      // save how data should be transferred to new level set
      // list of indices into the old pointData vector
      std::vector<std::vector<unsigned>> newDataSourceIds;
      if (updateData)
        newDataSourceIds.resize(newDomain.getNumberOfSegments());

#pragma omp parallel num_threads(newDomain.getNumberOfSegments())              \
    reduction(+ : addedPoints)
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

        // Per-thread guide iterator: advances in lock-step with main iterator.
        using GuideIt =
            viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>;
        std::unique_ptr<GuideIt> guideIt;
        if (guide != nullptr)
          guideIt = std::make_unique<GuideIt>(guide->getDomain(), startVector);

        for (viennahrle::ConstSparseStarIterator<
                 typename Domain<T, D>::DomainType, 1>
                 neighborIt(domain, startVector);
             neighborIt.getIndices() < endVector; neighborIt.next()) {

          auto &centerIt = neighborIt.getCenter();

          // If a guide is set, block activation for points inside it (φ < 0).
          bool insideGuide = false;
          if (guideIt) {
            guideIt->goToIndicesSequential(neighborIt.getIndices());
            insideGuide = (guideIt->getValue() < T(0));
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
            // Only activate when a proper defined neighbor drives the fill
            // and the point is not inside the guide (e.g. Si substrate).
            if (!insideGuide && definedDistance >= -limit &&
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
    // segment() is intentionally skipped here.  When the Interior fill
    // produces a dense HRLE (many small runs from a deep fill), the
    // balanced-partition boundaries computed by getNewSegmentation() land
    // in the middle of existing runs.  The SparseIterator inside segment()
    // then starts at the run's true beginning (before the boundary), causing
    // adjacent threads to overlap and producing a corrupt multi-segment domain
    // that goToIndices later walks off.  The single-segment domain left by
    // the final finalize()+deepCopy() is valid for all downstream consumers:
    // ConstSparseIterator in buildNodes(), writePersistentFields(), and
    // lsAdvect (which re-segments its output internally).
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(Interior)

} // namespace viennals
