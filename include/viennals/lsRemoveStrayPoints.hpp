#pragma once

#include <lsDomain.hpp>
#include <lsMarkVoidPoints.hpp>
#include <lsReduce.hpp>
#include <lsReinitialize.hpp>

namespace viennals {

using namespace viennacore;

/// This algorithm can be used to remove all LS values
/// which are not part of a so-called top surface.
/// This surface is detected using the MarkVoidPoints
/// algorithm, according to the method chosen by the user.
/// This method is set using setVoidTopSurface, which
/// is equivalent to the corresponding member function
/// of lsMarkVoidPoints.
/// A Euclidean level set is reduced to the narrow band the void detection
/// expects, where the bands of separate surfaces do not touch, keeps the point
/// data of the remaining points and is reinitialized from their planes to its
/// width (Reinitialize).
template <class T, int D> class RemoveStrayPoints {
  SmartPointer<Domain<T, D>> levelSet = nullptr;
  VoidTopSurfaceEnum voidTopSurface = VoidTopSurfaceEnum::LARGEST;

public:
  RemoveStrayPoints() = default;

  RemoveStrayPoints(SmartPointer<Domain<T, D>> passedLevelSet)
      : levelSet(passedLevelSet) {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  /// Set how the algorithm should pick the surface which will not
  /// be removed. Defaults to the surface with the most LS points.
  void setVoidTopSurface(VoidTopSurfaceEnum topSurface) {
    voidTopSurface = topSurface;
  }

  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addError("No level set was passed to lsPrune.")
          .print();
      return;
    }
    if (levelSet->getNumberOfPoints() == 0) {
      return;
    }

    const bool euclidean =
        levelSet->getNormalization() == NormalizationEnum::EUCLIDEAN;
    const int width = levelSet->getLevelSetWidth();
    if (euclidean)
      Reduce<T, D>(levelSet, 3).apply();

    // Mark which points are voids
    {
      MarkVoidPoints<T, D> marker;
      marker.setLevelSet(levelSet);
      marker.setVoidTopSurface(voidTopSurface);
      marker.apply();
    }

    auto voidMarkers =
        levelSet->getPointData().getScalarData("VoidPointMarkers", true);
    if (voidMarkers == nullptr) {
      Logger::getInstance()
          .addError("RemoveStrayPoints: No scalar data for void point "
                    "markers found. Cannot remove stray points.")
          .print();
    }

    // now iterate through the domain and remove points which are void points
    auto &grid = levelSet->getGrid();
    auto newlsDomain = SmartPointer<Domain<T, D>>::New(grid);
    typename Domain<T, D>::DomainType &newDomain = newlsDomain->getDomain();
    typename Domain<T, D>::DomainType &domain = levelSet->getDomain();

    newDomain.initialize(domain.getNewSegmentation(), domain.getAllocation());

    std::vector<typename Domain<T, D>::PointValueVectorType> newPoints;
    newPoints.resize(newDomain.getNumberOfSegments());
    // point ids of the remaining points, in the same order
    std::vector<std::vector<unsigned>> sourceIds(newPoints.size());

#pragma omp parallel for
    for (unsigned p = 0; p < newDomain.getNumberOfSegments(); ++p) {

      auto &domainSegment = newDomain.getDomainSegment(p);

      viennahrle::Index<D> const startVector =
          (p == 0) ? grid.getMinGridPoint()
                   : newDomain.getSegmentation()[p - 1];

      viennahrle::Index<D> const endVector =
          (p != static_cast<int>(newDomain.getNumberOfSegments() - 1))
              ? newDomain.getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      for (viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>
               it(domain, startVector);
           it.getStartIndices() < endVector; it.next()) {
        if (it.isDefined() && !voidMarkers->at(it.getPointId())) {
          newPoints[p].push_back(
              std::make_pair(it.getStartIndices(), it.getValue()));
          sourceIds[p].push_back(it.getPointId());
        }
      }
    }

    // merge points
    for (unsigned i = 1; i < newDomain.getNumberOfSegments(); ++i) {
      newPoints[0].insert(newPoints[0].end(), newPoints[i].begin(),
                          newPoints[i].end());
    }

    newlsDomain->insertPoints(newPoints[0]);
    if (euclidean) {
      for (unsigned i = 1; i < sourceIds.size(); ++i)
        sourceIds[0].insert(sourceIds[0].end(), sourceIds[i].begin(),
                            sourceIds[i].end());
      newlsDomain->getPointData().translateFromData(levelSet->getPointData(),
                                                    sourceIds[0]);
      newlsDomain->setNormalization(levelSet->getNormalization());
    }

    // distribute evenly across segments and copy
    newDomain.finalize();
    newDomain.segment();
    levelSet->deepCopy(newlsDomain);
    levelSet->finalize(2);
    if (euclidean)
      Reinitialize<T, D>(levelSet, width).apply();
  }
};

} // namespace viennals
