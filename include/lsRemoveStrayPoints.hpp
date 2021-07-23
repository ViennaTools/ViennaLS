#pragma once

#include <lsDomain.hpp>
#include <lsMarkVoidPoints.hpp>
#include <lsSmartPointer.hpp>

/// This algorithm can be used to remove all LS values
/// which are not part of a so-called top surface.
/// This surface is detected using the lsMarkVoidPoints
/// algorithm, according to the method chosen by the user.
/// This method is set using setVoidTopSurface, which
/// is equivalent to the corresponding member function
/// of lsMarkVoidPoints.
template <class T, int D> class lsRemoveStrayPoints {
  lsSmartPointer<lsDomain<T, D>> levelSet = nullptr;
  lsVoidTopSurfaceEnum voidTopSurface = lsVoidTopSurfaceEnum::LARGEST;

public:
  lsRemoveStrayPoints() {}

  lsRemoveStrayPoints(lsSmartPointer<lsDomain<T, D>> passedLevelSet)
      : levelSet(passedLevelSet) {}

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  /// Set how the algorithm should pick the surface which will not
  /// be removed. Defaults to the surface with the most LS points.
  void setVoidTopSurface(lsVoidTopSurfaceEnum topSurface) {
    voidTopSurface = topSurface;
  }

  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsPrune.")
          .print();
      return;
    }
    if (levelSet->getNumberOfPoints() == 0) {
      return;
    }

    // Mark which points are voids
    {
      lsMarkVoidPoints<T, D> marker;
      marker.setLevelSet(levelSet);
      marker.setVoidTopSurface(voidTopSurface);
      marker.apply();
    }

    auto voidMarkers = levelSet->getPointData().getScalarData("VoidPointMarkers");
    if(voidMarkers == nullptr) {
      lsMessage::getInstance().addWarning("lsRemoveStrayPoints: No scalar data for void point markers found. Cannot remove stray points.").print();
    }

    // now iterate through the domain and remove points which are void points
    auto &grid = levelSet->getGrid();
    auto newlsDomain = lsSmartPointer<lsDomain<T, D>>::New(grid);
    typename lsDomain<T, D>::DomainType &newDomain = newlsDomain->getDomain();
    typename lsDomain<T, D>::DomainType &domain = levelSet->getDomain();

    newDomain.initialize(domain.getNewSegmentation(), domain.getAllocation());

    std::vector<typename lsDomain<T, D>::PointValueVectorType> newPoints;
    newPoints.resize(newDomain.getNumberOfSegments());

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

      for (hrleConstSparseIterator<typename lsDomain<T, D>::DomainType>
               it(domain, startVector);
           it.getStartIndices() < endVector; it.next()) {
        if(it.isDefined() && !voidMarkers->at(it.getPointId())) {
          newPoints[p].push_back(std::make_pair(it.getStartIndices(), it.getValue()));
        }
      }
    }

    // merge points
    for(unsigned i = 1; i < newDomain.getNumberOfSegments(); ++i) {
      newPoints[0].insert(newPoints[0].end(),
                          newPoints[i].begin(),
                          newPoints[i].end());
    }

    newlsDomain->insertPoints(newPoints[0]);

    // distribute evenly across segments and copy
    newDomain.finalize();
    newDomain.segment();
    levelSet->deepCopy(newlsDomain);
    levelSet->finalize(2);
  }
};