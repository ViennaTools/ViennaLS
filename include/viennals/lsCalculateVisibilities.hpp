#pragma once

#include <lsDomain.hpp>

namespace viennals {
using namespace viennacore;

template <class NumericType, int D> class CalculateVisibilities {
  SmartPointer<Domain<NumericType, D>> levelSet;
  int openBoundaryDirection;
  bool isOpenBoundaryNegative;

public:
  static constexpr char visibilitiesLabel[] = "Visibilities";

  CalculateVisibilities(
      const SmartPointer<Domain<NumericType, D>> &passedLevelSet,
      int passedOpenBoundaryDirection, bool passedIsOpenBoundaryNegative)
      : levelSet(passedLevelSet),
        openBoundaryDirection(passedOpenBoundaryDirection),
        isOpenBoundaryNegative(passedIsOpenBoundaryNegative) {}

  void apply() {

    auto &domain = levelSet->getDomain();
    auto &grid = levelSet->getGrid();
    const NumericType max = std::numeric_limits<NumericType>::max();

    auto numDefinedPoints = domain.getNumberOfPoints();
    std::vector<NumericType> visibilities(numDefinedPoints);

    std::vector<hrleIndexType> oldIndices(
        D - 1 - openBoundaryDirection,
        std::numeric_limits<hrleIndexType>::max());

    unsigned int size = 1;
    for (int i = 0; i < openBoundaryDirection; ++i) {
      assert(!grid.isPosBoundaryInfinite(i));
      assert(!grid.isNegBoundaryInfinite(i));

      size *= (grid.getMaxIndex(i) - grid.getMinIndex(i) + 1);
    }

    std::vector<NumericType> minValues(size, max);

    hrleSizeType id = 0;

    hrleSparseIterator<typename Domain<NumericType, D>::DomainType> it(
        domain, !isOpenBoundaryNegative);
    while (!it.isFinished()) {

      for (int i = 0; i < D - 1 - openBoundaryDirection; ++i) {
        bool b = false;
        if (oldIndices[i] !=
            it.getStartIndices(i + openBoundaryDirection + 1)) {
          oldIndices[i] = it.getStartIndices(i + openBoundaryDirection + 1);
          b = true;
        }
        if (b)
          minValues.assign(size, max);
      }

      unsigned int pos_begin = 0;
      unsigned int pos_end = 0;

      for (int i = openBoundaryDirection - 1; i >= 0; --i) {
        pos_begin *= (grid.getMaxIndex(i) - grid.getMinIndex(i) + 1);
        pos_end *= (grid.getMaxIndex(i) - grid.getMinIndex(i) + 1);
        pos_begin += (it.getStartIndices(i) - grid.getMinIndex(i));
        pos_end += (it.getEndIndices(i) - grid.getMinIndex(i));
      }

      if (it.isDefined()) {
        visibilities[isOpenBoundaryNegative ? id
                                            : (numDefinedPoints - 1 - id)] =
            (it.getValue() < minValues.at(pos_begin)) ? 1. : 0.;
        ++id;
      }

      for (unsigned int i = pos_begin; i <= pos_end; ++i)
        minValues.at(i) = std::min(minValues.at(i), it.getValue());

      if (isOpenBoundaryNegative) {
        it.next();
      } else {
        it.previous();
      }
    }

    levelSet->getPointData().insertNextScalarData(visibilities,
                                                  visibilitiesLabel);

    assert(id == numDefinedPoints);
  }
};

} // namespace viennals