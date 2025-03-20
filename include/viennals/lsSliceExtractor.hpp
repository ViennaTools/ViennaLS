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
template <class T> class SliceExtractor {
private:
  SmartPointer<Domain<T, 3>> sourceDomain = nullptr;
  SmartPointer<Domain<T, 2>> sliceDomain = nullptr;

  int sliceDimension = 0; // (0=x, 1=y, 2=z)
  T slicePosition = 0;    // Position of the slice
  T tolerance = 1e-6;     // Tolerance for matching slice position

public:
  SliceExtractor() {}

  SliceExtractor(SmartPointer<Domain<T, 3>> passedDomain,
                 SmartPointer<Domain<T, 2>> passedSliceDomain,
                 int passedSliceDimension = 0, T passedSlicePosition = 0)
      : sourceDomain(passedDomain), sliceDomain(passedSliceDomain),
        sliceDimension(passedSliceDimension),
        slicePosition(passedSlicePosition) {}

  void setSourceDomain(SmartPointer<Domain<T, 3>> passedDomain) {
    sourceDomain = passedDomain;
  }

  void setSliceDomain(SmartPointer<Domain<T, 2>> passedDomain) {
    sliceDomain = passedDomain;
  }

  void setSliceDimension(int dimension) {
    if (dimension >= 0 && dimension < 3) {
      sliceDimension = dimension;
    } else {
      Logger::getInstance()
          .addError("Invalid slice dimension. Must be 0 (x), 1 (y), or 2 (z)")
          .print();
    }
  }

  void setSlicePosition(T position) { slicePosition = position; }

  void setTolerance(T passedTolerance) { tolerance = passedTolerance; }

  void apply() {
    if (sourceDomain == nullptr) {
      Logger::getInstance()
          .addError("Source domain is null in SliceExtractor")
          .print();
      return;
    }

    if (sliceDomain == nullptr) {
      Logger::getInstance()
          .addError("Slice domain is null in SliceExtractor")
          .print();
      return;
    }

    const auto &sourceGrid = sourceDomain->getGrid();
    auto const gridDelta = sourceGrid.getGridDelta();

    if (tolerance >= gridDelta) {
      Logger::getInstance()
          .addWarning("Tolerance is greater equal grid delta in "
                      "SliceExtractor. This might lead to unexpected results.")
          .print();
    }

    // Create bounds for the slice domain
    double sliceBounds[4];
    int sliceIdx = 0;
    for (int d = 0; d < 3; d++) {
      if (d != sliceDimension) {
        sliceBounds[2 * sliceIdx] = sourceGrid.getMinGridPoint(d) * gridDelta;
        sliceBounds[2 * sliceIdx + 1] =
            sourceGrid.getMaxGridPoint(d) * gridDelta;
        sliceIdx++;
      }
    }

    // Container for the extracted points
    std::vector<std::pair<hrleVectorType<hrleIndexType, 2>, T>> pointData;

    // Iterate through the source domain
    hrleConstSparseIterator<typename Domain<T, 3>::DomainType> it(
        sourceDomain->getDomain());

    while (!it.isFinished()) {
      if (!it.isDefined()) {
        it.next();
        continue;
      }

      auto indices = it.getStartIndices();

      // Check if this point is in requested slice (within tolerance)
      T coord = indices[sliceDimension] * gridDelta;
      if (std::abs(coord - slicePosition) <= tolerance) {
        // Extract the value
        T value = it.getValue();

        // Create a new 2D index
        hrleVectorType<hrleIndexType, 2> sliceIndices;
        int sliceIdx = 0;
        for (int d = 0; d < 3; d++) {
          if (d != sliceDimension) {
            sliceIndices[sliceIdx++] = indices[d];
          }
        }

        // Add to our collection
        pointData.push_back(std::make_pair(sliceIndices, value));
      }

      it.next();
    }

    // Insert the extracted points into the slice domain, issue a warning if
    // there are no points to insert
    if (pointData.empty()) {
      Logger::getInstance()
          .addWarning("No points extracted in SliceExtractor")
          .print();
    } else {
      sliceDomain->insertPoints(pointData);
    }
  }
};

// Add template specializations for this class
PRECOMPILE_PRECISION(SliceExtractor)

} // namespace viennals
