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
  bool updatePointData = true;
  bool removeStrayZeros = false;

  template <class Numeric> bool isNegative(const Numeric a) {
    return a <= -std::numeric_limits<Numeric>::epsilon();
  }

  template <class Numeric>
  bool isSignDifferent(const Numeric a, const Numeric b) {
    return (isNegative(a) ^ isNegative(b));
  }

  // small helper to check whether LS function is monotone
  // around a zero value
  bool isMonotone(const T a, const T b, const T c) {
    const auto diff1 = a - b;
    const auto diff2 = b - c;
    return !(isSignDifferent(diff1, diff2) && a != 0. && c != 0.);
  }

public:
  lsPrune() {}

  lsPrune(lsSmartPointer<lsDomain<T, D>> passedlsDomain)
      : levelSet(passedlsDomain){};

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedlsDomain) {
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
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsPrune.")
          .print();
      return;
    }
    if (levelSet->getNumberOfPoints() == 0) {
      return;
    }

    auto &grid = levelSet->getGrid();
    auto newlsDomain = lsSmartPointer<lsDomain<T, D>>::New(grid);
    typename lsDomain<T, D>::DomainType &newDomain = newlsDomain->getDomain();
    typename lsDomain<T, D>::DomainType &domain = levelSet->getDomain();

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
        bool centerSign = isNegative(centerIt.getValue());
        if (centerIt.isDefined()) {
          int i = 0;
          for (; i < 2 * D; i++) {
            if (isSignDifferent(neighborIt.getNeighbor(i).getValue(),
                                centerIt.getValue())) {
              break;
            }
          }

          if(removeStrayZeros) {
            // if the centre point is 0.0 and the level set values
            // along each grid dimension are not monotone, it is
            // a numerical glitch and should be removed
            if(const auto &midVal = centerIt.getValue(); midVal == std::abs(0.)) {
              int undefVal = 0;
              for (int i = 0; i < D; i++) {
                const auto &negVal = neighborIt.getNeighbor(i).getValue();
                const auto &posVal = neighborIt.getNeighbor(D + i).getValue();

                if(!isMonotone(negVal, midVal, posVal)) {
                  undefVal = isNegative(negVal) ? -1 : 1;
                  break;
                }
              }
              if(undefVal != 0) {
                domainSegment.insertNextUndefinedPoint(
                  neighborIt.getIndices(), undefVal == -1
                                              ? lsDomain<T, D>::NEG_VALUE
                                              : lsDomain<T, D>::POS_VALUE);
                continue;
              }
            }
          }

          if (i != 2 * D) {
            domainSegment.insertNextDefinedPoint(neighborIt.getIndices(),
                                                 centerIt.getValue());
            if (updateData)
              newDataSourceIds[p].push_back(centerIt.getPointId());
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
PRECOMPILE_PRECISION_DIMENSION(lsPrune)

#endif // LS_PRUNE_HPP
