#pragma once

#include <unordered_map>
#include <limits>

#include <lsPreCompileMacros.hpp>
#include <hrleDenseCellIterator.hpp>
#include <lsDomain.hpp>

namespace viennals {

using namespace viennacore;

/// Computes an estimate of the area where two level sets differ.
/// Custom increment values can be set for specific x and y ranges.
template <class T, int D> class CompareArea {
  typedef typename Domain<T, D>::DomainType hrleDomainType;

  SmartPointer<Domain<T, D>> levelSetTarget = nullptr;
  SmartPointer<Domain<T, D>> levelSetSample = nullptr;
  hrleVectorType<hrleIndexType, D> minIndex, maxIndex;

  unsigned long int cellCount = 0;
  unsigned long int customCellCount = 0;
  double mismatchedArea = 0.;
  double customMismatchedArea = 0.;

  hrleIndexType xRangeMin = std::numeric_limits<hrleIndexType>::lowest();
  hrleIndexType xRangeMax = std::numeric_limits<hrleIndexType>::max();
  hrleIndexType yRangeMin = std::numeric_limits<hrleIndexType>::lowest();
  hrleIndexType yRangeMax = std::numeric_limits<hrleIndexType>::max();
  bool useCustomXIncrement = false;
  bool useCustomYIncrement = false;
  
  unsigned short int customXIncrement = 0;
  unsigned short int customYIncrement = 0;
  unsigned short int defaultIncrement = 1;
  
  double gridDelta = 0.0;

  bool checkAndCalculateBounds() {
    if (levelSetTarget == nullptr || levelSetSample == nullptr) {
      Logger::getInstance()
          .addError("Missing level set in CompareArea.")
          .print();
      return false;
    }

    // Check if the grids are compatible
    const auto &gridTarget = levelSetTarget->getGrid();
    const auto &gridSample = levelSetSample->getGrid();
    
    if (gridTarget.getGridDelta() != gridSample.getGridDelta()) {
      Logger::getInstance()
          .addError("Grid delta mismatch in CompareArea. The grid deltas of the two level sets must be equal.")
          .print();
      return false;
    } else {
      gridDelta = gridTarget.getGridDelta();
    }

    // Initialize min and max indices
    for (unsigned i = 0; i < D; ++i) {
      minIndex[i] = std::numeric_limits<hrleIndexType>::max();
      maxIndex[i] = std::numeric_limits<hrleIndexType>::lowest();
    }

    // Get the boundaries of the two level sets
    auto &domainTarget = levelSetTarget->getDomain();
    auto &domainSample = levelSetSample->getDomain();

    // Calculate actual bounds
    for (unsigned i = 0; i < D; ++i) {
      minIndex[i] = std::min({minIndex[i], 
                             (gridTarget.isNegBoundaryInfinite(i)) ? domainTarget.getMinRunBreak(i) : gridTarget.getMinBounds(i),
                             (gridSample.isNegBoundaryInfinite(i)) ? domainSample.getMinRunBreak(i) : gridSample.getMinBounds(i)});

      maxIndex[i] = std::max({maxIndex[i], 
                             (gridTarget.isPosBoundaryInfinite(i)) ? domainTarget.getMaxRunBreak(i) : gridTarget.getMaxBounds(i),
                             (gridSample.isPosBoundaryInfinite(i)) ? domainSample.getMaxRunBreak(i) : gridSample.getMaxBounds(i)});
    }
    
    return true;
  }

public:
  CompareArea() {}

  CompareArea(SmartPointer<Domain<T, D>> passedLevelSetTarget,
             SmartPointer<Domain<T, D>> passedLevelSetSample)
      : levelSetTarget(passedLevelSetTarget), levelSetSample(passedLevelSetSample) {}

  /// Sets the target level set.
  void setLevelSetTarget(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetTarget = passedLevelSet;
  }

  /// Sets the sample level set.
  void setLevelSetSample(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetSample = passedLevelSet;
  }

  /// Set default increment value
  void setDefaultIncrement(unsigned short int increment) { defaultIncrement = increment; }

  /// Sets the x-range and custom increment value
  void setXRangeAndIncrement(hrleIndexType minXRange, hrleIndexType maxXRange, unsigned short int Xincrement) {
    xRangeMin = minXRange;
    xRangeMax = maxXRange;
    customXIncrement = Xincrement;
    useCustomXIncrement = true;
  }

  /// Sets the y-range and custom increment value
  void setYRangeAndIncrement(hrleIndexType minYRange, hrleIndexType maxYRange, unsigned short int Yincrement) {
    yRangeMin = minYRange;
    yRangeMax = maxYRange;
    customYIncrement = Yincrement;
    useCustomYIncrement = true;
  }

  /// Returns the computed area mismatch.
  double getAreaMismatch() const { return static_cast<double>(cellCount) * gridDelta * gridDelta; }

  /// Returns the computed area mismatch, with custom increments applied.
  double getCustomAreaMismatch() const { return static_cast<double>(customCellCount) * gridDelta * gridDelta; }

  /// Returns the number of cells where the level sets differ.
  unsigned long int getCellCount() const { return cellCount; }

  /// Returns the number of cells where the level sets differ, with custom increments applied.
  unsigned long int getCustomCellCount() const { return customCellCount; }

  /// Computes the area difference between the two level sets.
  void apply() {
    // Calculate the bounds for iteration
    if (!checkAndCalculateBounds()) {
      cellCount = std::numeric_limits<unsigned long int>::max(); // Error indicator
      customCellCount = std::numeric_limits<unsigned long int>::max(); // Error indicator
      return;
    }
    
    // Set up dense cell iterators for both level sets
    hrleConstDenseCellIterator<typename Domain<T, D>::DomainType> itTarget(levelSetTarget->getDomain(), minIndex);
    hrleConstDenseCellIterator<typename Domain<T, D>::DomainType> itSample(levelSetSample->getDomain(), minIndex);

    cellCount = 0;
    customCellCount = 0;

    // Iterate through the domain defined by the bounding box
    for (; itTarget.getIndices() < maxIndex; itTarget.next()) {
      itSample.goToIndicesSequential(itTarget.getIndices());

      // Compare cell values at the current position
      T centerValueTarget = 0.;
      T centerValueSample = 0.;

      // Calculate the average value of the corners for both level sets
      for (int i = 0; i < (1 << D); ++i) {
        centerValueTarget += itTarget.getCorner(i).getValue();
        centerValueSample += itSample.getCorner(i).getValue();
      }
      centerValueTarget /= (1 << D);
      centerValueSample /= (1 << D);

      // Determine if cell is inside and outside for level sets
      bool isDifferent = (centerValueTarget <= 0.) != (centerValueSample <= 0.);

      // If one cell is inside and the other is outside, apply the increment
      if (isDifferent) {
        // Always increment simple cell count by 1
        cellCount += 1;
        
        // For custom cell count, apply special rules
        unsigned short int incrementToAdd = defaultIncrement;
        
        bool inXRange = useCustomXIncrement && 
                (itTarget.getIndices()[0] * gridDelta >= xRangeMin && itTarget.getIndices()[0] * gridDelta <= xRangeMax);
        bool inYRange = useCustomYIncrement && 
                (itTarget.getIndices()[1] * gridDelta >= yRangeMin && itTarget.getIndices()[1] * gridDelta <= yRangeMax);

        if (inXRange && inYRange) {
          // Apply both increments
          incrementToAdd = customXIncrement + customYIncrement;
        } else if (inXRange) {
          incrementToAdd = customXIncrement;
        } else if (inYRange) {
          incrementToAdd = customYIncrement;
        }
        
        customCellCount += incrementToAdd;
      }
    }
  }
};

// Precompile for common precision and dimensions
PRECOMPILE_PRECISION_DIMENSION(CompareArea)

} // namespace viennals
