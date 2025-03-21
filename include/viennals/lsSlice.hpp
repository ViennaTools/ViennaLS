#pragma once

#include <lsDomain.hpp>
#include <lsPreCompileMacros.hpp>

namespace viennals {

using namespace viennacore;

/// Extract a 2D slice from a 3D level set domain at a fixed position along one
/// axis. The resulting 2D domain contains points where the 3D domain intersects
/// the slice plane and is always inserted into the x-y plane. Use with caution,
/// as the result might be an empty domain if the slice does not intersect any
/// defined points in the source domain.
template <class T> class Slice {
  SmartPointer<Domain<T, 3>> sourceLevelSet = nullptr;
  SmartPointer<Domain<T, 2>> sliceLevelSet = nullptr;

  int sliceDimension = 0; // (0=x, 1=y, 2=z)
  T slicePosition = 0;    // Position of the slice

public:
  Slice() = default;

  Slice(SmartPointer<Domain<T, 3>> passedDomain,
        SmartPointer<Domain<T, 2>> passedSliceDomain,
        int passedSliceDimension = 0, T passedSlicePosition = 0)
      : sourceLevelSet(passedDomain), sliceLevelSet(passedSliceDomain),
        sliceDimension(passedSliceDimension),
        slicePosition(passedSlicePosition) {}

  void setSourceLevelSet(SmartPointer<Domain<T, 3>> passedDomain) {
    sourceLevelSet = passedDomain;
  }

  void setSliceLevelSet(SmartPointer<Domain<T, 2>> passedDomain) {
    sliceLevelSet = passedDomain;
  }

  void setSliceDimension(const int dimension) {
    if (dimension >= 0 && dimension < 3) {
      sliceDimension = dimension;
    } else {
      Logger::getInstance()
          .addError("Invalid slice dimension. Must be 0 (x), 1 (y), or 2 (z)")
          .print();
    }
  }

  void setSlicePosition(T position) { slicePosition = position; }

  void apply() {
    if (sourceLevelSet == nullptr) {
      Logger::getInstance()
          .addError("No source level-set passed to Slice")
          .print();
      return;
    }

    if (sliceLevelSet == nullptr) {
      Logger::getInstance()
          .addError("No slice level-set passed to Slice")
          .print();
      return;
    }

    const auto &sourceGrid = sourceLevelSet->getGrid();
    auto const gridDelta = sourceGrid.getGridDelta();

    // Create bounds for the slice domain
    double sliceBounds[4];
    BoundaryConditionEnum sliceBoundaryConds[2];

    for (int i = 0, d = 0; d < 3; d++) {
      if (d != sliceDimension) {
        sliceBounds[2 * i] = sourceGrid.getMinGridPoint(d) * gridDelta;
        sliceBounds[2 * i + 1] = sourceGrid.getMaxGridPoint(d) * gridDelta;
        sliceBoundaryConds[i] = sourceGrid.getBoundaryConditions(d);
        i++;
      }
    }

    // Container for the extracted points
    std::vector<std::pair<hrleVectorType<hrleIndexType, 2>, T>> pointData;

    // slice index
    const int sliceIndex = static_cast<int>(slicePosition / gridDelta);

    // Iterate through the source domain
    hrleConstSparseIterator<typename Domain<T, 3>::DomainType> it(
        sourceLevelSet->getDomain());

    while (!it.isFinished()) {
      if (!it.isDefined()) {
        it.next();
        continue;
      }

      auto indices = it.getStartIndices();
      if (indices[sliceDimension] == sliceIndex) {
        // Extract the value
        T value = it.getValue();

        // Create a new 2D index
        hrleVectorType<hrleIndexType, 2> sliceIndices;
        for (int d = 0, j = 0; d < 3; d++) {
          if (d != sliceDimension) {
            sliceIndices[j++] = indices[d];
          }
        }

        // Add to our collection
        pointData.emplace_back(sliceIndices, value);
      }

      it.next();
    }

    // Insert the extracted points into the slice domain, issue a warning if
    // there are no points to insert
    if (pointData.empty()) {
      Logger::getInstance().addWarning("No points extracted in Slice").print();
    } else {
      auto slice = SmartPointer<Domain<T, 2>>::New(
          pointData, sliceBounds, sliceBoundaryConds, gridDelta);
      sliceLevelSet->deepCopy(slice);
    }
  }
};

// Add template specializations for this class
PRECOMPILE_PRECISION(Slice)

} // namespace viennals
