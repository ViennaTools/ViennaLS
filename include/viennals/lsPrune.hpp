#pragma once

#include <lsPreCompileMacros.hpp>

#include <hrleSparseStarIterator.hpp>
#include <hrleVectorType.hpp>

#include <lsDomain.hpp>

namespace viennals {

using namespace viennacore;

/// Removes all level set points, which do not have
/// at least one oppositely signed neighbour (Meaning
/// they do not lie directly at the interface).
/// Afterwards the level set will occupy the least memory
/// possible.
template <class T, int D> class Prune {
  SmartPointer<Domain<T, D>> levelSet = nullptr;
  bool updatePointData = true;
  bool removeStrayZeros = false;

  template <class Numeric> static bool isNegative(const Numeric a) {
    return a <= -std::numeric_limits<Numeric>::epsilon();
  }

  template <class Numeric>
  static bool isSignDifferent(const Numeric a, const Numeric b) {
    return (isNegative(a) ^ isNegative(b));
  }

  template <class Numeric>
  static bool isSignDifferentOrZero(const Numeric a, const Numeric b) {
    if (a == 0. || b == 0.)
      return true;
    return isSignDifferent(a, b);
  }

  template <class hrleIterator, class Compare>
  static bool checkNeighbourSigns(const hrleIterator &it, Compare comp) {
    for (int i = 0; i < 2 * D; i++) {
      if (comp(it.getCenter().getValue(), it.getNeighbor(i).getValue())) {
        return true;
      }
    }
    return false;
  }

  // small helper to check whether LS function is monotone
  // around a zero value
  static bool isMonotone(const T a, const T c) {
    return a == 0. || c == 0. || isSignDifferent(a, c);
  }

public:
  Prune() {}

  Prune(SmartPointer<Domain<T, D>> passedlsDomain)
      : levelSet(passedlsDomain) {};

  void setLevelSet(SmartPointer<Domain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  /// Set whether to update the point data stored in the LS
  /// during this algorithm. Defaults to true.
  void setUpdatePointData(bool update) { updatePointData = update; }

  /// Set whether to remove exact zero values between grid
  /// points with the same sign
  void setRemoveStrayZeros(bool rsz) { removeStrayZeros = rsz; }

  /// removes all grid points, which do not have at least one opposite signed
  /// neighbour
  /// returns the number of removed points
  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addWarning("No level set was passed to Prune.")
          .print();
      return;
    }
    if (levelSet->getNumberOfPoints() == 0) {
      return;
    }

    auto &grid = levelSet->getGrid();
    auto newlsDomain = SmartPointer<Domain<T, D>>::New(grid);
    typename Domain<T, D>::DomainType &newDomain = newlsDomain->getDomain();
    typename Domain<T, D>::DomainType &domain = levelSet->getDomain();

    newDomain.initialize(domain.getNewSegmentation(), domain.getAllocation());

    const bool updateData = updatePointData;
    const bool removeZeros = removeStrayZeros;
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
        bool centerSign = isNegative(centerIt.getValue());
        if (centerIt.isDefined()) {
          bool keepPoint = true;
          // if center is exact zero, always treat as unprunable
          if (std::abs(centerIt.getValue()) != 0.) {
            if (removeZeros) {
              keepPoint =
                  checkNeighbourSigns(neighborIt, isSignDifferentOrZero<T>);
            } else {
              keepPoint = checkNeighbourSigns(neighborIt, isSignDifferent<T>);
            }
          }

          if (removeZeros) {
            // if the centre point is 0.0 and the level set values
            // along each grid dimension are not monotone, it is
            // a numerical glitch and should be removed
            if (std::abs(centerIt.getValue()) == 0.) {
              bool overWritePoint = false;
              T undefVal = 0.;
              for (int i = 0; i < D; i++) {
                const auto &negVal = neighborIt.getNeighbor(i).getValue();
                const auto &posVal = neighborIt.getNeighbor(D + i).getValue();

                // if LS function is not monotone around the zero value,
                // set the points value to that of the lower neighbour
                if (!isMonotone(negVal, posVal)) {
                  overWritePoint = true;
                  undefVal =
                      std::abs(negVal) < std::abs(posVal) ? negVal : posVal;
                  break;
                }
              }
              if (overWritePoint) {
                domainSegment.insertNextDefinedPoint(neighborIt.getIndices(),
                                                     undefVal);
                continue;
              }
            }
          }

          if (keepPoint) {
            domainSegment.insertNextDefinedPoint(neighborIt.getIndices(),
                                                 centerIt.getValue());
            if (updateData)
              newDataSourceIds[p].push_back(centerIt.getPointId());
          } else {
            // TODO: it is more efficient to insertNextUndefinedRunType, since
            // we know it already exists
            domainSegment.insertNextUndefinedPoint(
                neighborIt.getIndices(),
                centerSign ? Domain<T, D>::NEG_VALUE : Domain<T, D>::POS_VALUE);
          }
        } else {
          domainSegment.insertNextUndefinedPoint(
              neighborIt.getIndices(),
              centerSign ? Domain<T, D>::NEG_VALUE : Domain<T, D>::POS_VALUE);
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
    newDomain.segment();
    levelSet->deepCopy(newlsDomain);
    levelSet->finalize(2);
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(Prune)

} // namespace viennals
