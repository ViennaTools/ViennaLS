#pragma once

#include <hrleSparseIterator.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMesh.hpp>
#include <lsPreCompileMacros.hpp>
#include <lsReduce.hpp>

#include <cmath>
#include <limits>

namespace viennals {

using namespace viennacore;

/// Calculate distance measure between two level sets by comparing their SDF
/// values on a sparse field. This class iterates over the points in the sparse
/// field of the iterated level set and calculates differences with the
/// corresponding values of the expanded level set.
///
/// The expanded level set is expected to be expanded in such a way that the
/// sparse field of the iterated level set always overlaps with defined values
/// in the expanded level set. If a level set with a width < 50 is passed as the
/// expaned one, the function will automatically expand it to a width of 50, a
/// value which might be sufficient for many purposes. However, making sure that
/// the level set is sufficiently expanded is the responsibility of the user.
///
/// The iterated level set is expected to be sparse. The reduction is performed
/// automatically if this is not the case.
///
/// The code is currently intended for 2D level sets only.

template <class T, int D = 2> class CompareSparseField {
  using hrleIndexType = viennahrle::IndexType;

  SmartPointer<Domain<T, D>> levelSetExpanded = nullptr;
  SmartPointer<Domain<T, D>> levelSetIterated = nullptr;

  // Variables for x and y range restrictions
  T xRangeMin = std::numeric_limits<T>::lowest();
  T xRangeMax = std::numeric_limits<T>::max();
  T yRangeMin = std::numeric_limits<T>::lowest();
  T yRangeMax = std::numeric_limits<T>::max();
  bool useXRange = false;
  bool useYRange = false;

  // Fields to store the calculation results
  T sumSquaredDifferences = 0.0;
  T sumDifferences = 0.0;
  unsigned numPoints = 0;
  unsigned numSkippedPoints = 0;

  // Add mesh output capability
  SmartPointer<Mesh<T>> outputMesh = nullptr;

  // Add bool for filling the iterated level set with distances
  bool fillIteratedWithDistances = false;

  // Expansion width for the expanded level set
  int expandedLevelSetWidth = 50;

  bool checkAndCalculateBounds() {
    if (levelSetExpanded == nullptr || levelSetIterated == nullptr) {
      Logger::getInstance()
          .addError("Missing level set in CompareSparseField.")
          .print();
      return false;
    }

    // Check if the grids are compatible
    const auto &gridExpanded = levelSetExpanded->getGrid();
    const auto &gridIterated = levelSetIterated->getGrid();

    if (gridExpanded.getGridDelta() != gridIterated.getGridDelta()) {
      Logger::getInstance()
          .addError("Grid delta mismatch in CompareSparseField. The grid "
                    "deltas of the two level sets must be equal.")
          .print();
      return false;
    }

    // // Check if the x extents of both level sets are equal
    // const auto &domainExpanded = levelSetExpanded->getDomain();
    // const auto &domainIterated = levelSetIterated->getDomain();

    // hrleIndexType expandedMinX = gridExpanded.isNegBoundaryInfinite(0)
    //                                ? domainExpanded.getMinRunBreak(0)
    //                                : gridExpanded.getMinIndex(0);
    // hrleIndexType expandedMaxX = gridExpanded.isPosBoundaryInfinite(0)
    //                                ? domainExpanded.getMaxRunBreak(0)
    //                                : gridExpanded.getMaxIndex(0);
    // hrleIndexType iteratedMinX = gridIterated.isNegBoundaryInfinite(0)
    //                                ? domainIterated.getMinRunBreak(0)
    //                                : gridIterated.getMinIndex(0);
    // hrleIndexType iteratedMaxX = gridIterated.isPosBoundaryInfinite(0)
    //                                ? domainIterated.getMaxRunBreak(0)
    //                                : gridIterated.getMaxIndex(0);

    // if (expandedMinX != iteratedMinX || expandedMaxX != iteratedMaxX) {
    //   Logger::getInstance()
    //       .addWarning("X extent mismatch in CompareSparseField. The x extents
    //       "
    //                   "of both level sets must be equal.")
    //       .print();
    //   return false;
    // }

    // Check if expanded level set width is sufficient
    if (levelSetExpanded->getLevelSetWidth() < expandedLevelSetWidth) {
      VIENNACORE_LOG_WARNING(
          "Expanded level set width is insufficient. It must have a width of "
          "at least " +
          std::to_string(expandedLevelSetWidth) + ". \n" +
          " CORRECTION: The expansion was performed. \n"
          "ALTERNATIVE: Alternatively, please expand the expanded yourself "
          "using lsExpand before passing it to this function. \n");
      Expand<T, D>(levelSetExpanded, expandedLevelSetWidth).apply();
    }

    // Reduce the iterated level set to a sparse field if necessary
    if (levelSetIterated->getLevelSetWidth() > 1) {
      VIENNACORE_LOG_WARNING(
          "Iterated level set width is too large. It must be reduced to a "
          "sparse field. \n"
          " CORRECTION: The reduction was performed. \n"
          "ALTERNATIVE: Alternatively, please reduce the iterated yourself "
          "using lsReduce before passing it to this function. \n");
      Reduce<T, D>(levelSetIterated, 1).apply();
    }

    return true;
  }

public:
  CompareSparseField() {}

  CompareSparseField(SmartPointer<Domain<T, D>> passedLevelSetExpanded,
                     SmartPointer<Domain<T, D>> passedLevelSetIterated)
      : levelSetExpanded(passedLevelSetExpanded),
        levelSetIterated(passedLevelSetIterated) {}

  void setLevelSetExpanded(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetExpanded = passedLevelSet;
  }

  void setLevelSetIterated(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetIterated = passedLevelSet;
  }

  /// Set the x-coordinate range to restrict the comparison area
  void setXRange(T minXRange, T maxXRange) {
    xRangeMin = minXRange;
    xRangeMax = maxXRange;
    useXRange = true;
  }

  /// Set the y-coordinate range to restrict the comparison area
  void setYRange(T minYRange, T maxYRange) {
    yRangeMin = minYRange;
    yRangeMax = maxYRange;
    useYRange = true;
  }

  /// Clear the x-range restriction
  void clearXRange() {
    useXRange = false;
    xRangeMin = std::numeric_limits<T>::lowest();
    xRangeMax = std::numeric_limits<T>::max();
  }

  /// Clear the y-range restriction
  void clearYRange() {
    useYRange = false;
    yRangeMin = std::numeric_limits<T>::lowest();
    yRangeMax = std::numeric_limits<T>::max();
  }

  /// Set the output mesh where difference values will be stored
  void setOutputMesh(SmartPointer<Mesh<T>> passedMesh) {
    outputMesh = passedMesh;
  }

  /// Set whether to fill the iterated level set with distances
  void setFillIteratedWithDistances(bool fill) {
    fillIteratedWithDistances = fill;
  }

  /// Set the expansion width for the expanded level set
  /// This value will be used if the expanded level set needs to be expanded
  /// automatically during the apply() call
  void setExpandedLevelSetWidth(int width) {
    if (width <= 0) {
      VIENNACORE_LOG_WARNING(
          "Expansion width must be positive. Using default value of 50.");
      expandedLevelSetWidth = 50;
    } else {
      expandedLevelSetWidth = width;
    }
  }

  /// Apply the comparison and calculate the sum of squared differences.
  void apply() {
    // Perform compatibility checks
    if (!checkAndCalculateBounds()) {
      // If checks fail, return NaN
      sumSquaredDifferences = std::numeric_limits<T>::quiet_NaN();
      sumDifferences = std::numeric_limits<T>::quiet_NaN();
      numPoints = 0;
      numSkippedPoints = 0;
      return;
    }

    const auto &gridExpanded = levelSetExpanded->getGrid();
    const auto gridDelta = gridExpanded.getGridDelta();

    sumSquaredDifferences = 0.0;
    sumDifferences = 0.0;
    numPoints = 0;
    numSkippedPoints = 0;

    // Direct storage for points and differences
    std::vector<Vec3D<T>> nodeCoordinates; // 3D necessary for lsMesh
    std::vector<std::array<unsigned, 1>> vertexIndices;
    std::vector<T> differenceValues;
    std::vector<T> squaredDifferenceValues;

    // Prepare mesh output if needed
    const bool generateMesh = outputMesh != nullptr;
    if (generateMesh) {
      outputMesh->clear();

      // Initialize mesh extent
      for (unsigned i = 0; i < D; ++i) {
        outputMesh->minimumExtent[i] = std::numeric_limits<T>::max();
        outputMesh->maximumExtent[i] = std::numeric_limits<T>::lowest();
      }

      // Reserve space for mesh data
      nodeCoordinates.reserve(levelSetIterated->getNumberOfPoints());
      vertexIndices.reserve(levelSetIterated->getNumberOfPoints());
      differenceValues.reserve(levelSetIterated->getNumberOfPoints());
      squaredDifferenceValues.reserve(levelSetIterated->getNumberOfPoints());
    }

    // Prepare for point data filling if needed
    auto &iteratedPointData = levelSetIterated->getPointData();
    std::vector<T> pointDataDistances;
    if (fillIteratedWithDistances) {
      pointDataDistances.reserve(levelSetIterated->getNumberOfPoints());
    }

    // Create sparse iterators for the level sets
    viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>
        itIterated(levelSetIterated->getDomain());
    viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>
        itExpanded(levelSetExpanded->getDomain());

    // Iterate over all defined points in the iterated level set
    while (!itIterated.isFinished()) {
      if (!itIterated.isDefined()) {
        // this block is necessary to skip undefined points, I tested it:
        // std::cout << "Skipping undefined point" << std::endl;
        itIterated.next();
        continue;
      }

      auto indices = itIterated.getStartIndices();

      // Calculate coordinates
      T xCoord = indices[0] * gridDelta;
      T yCoord = indices[1] * gridDelta;
      T zCoord = (D == 3) ? indices[2] * gridDelta
                          : 0.0; // Always use 0 for z-coordinate in 2D

      // Skip if outside the specified x-range
      if (useXRange && (xCoord < xRangeMin || xCoord > xRangeMax)) {
        itIterated.next();
        continue;
      }

      // Skip if outside the specified y-range
      if (useYRange && (yCoord < yRangeMin || yCoord > yRangeMax)) {
        itIterated.next();
        continue;
      }

      // Get iterated value
      T valueIterated = itIterated.getValue();

      itExpanded.goToIndicesSequential(indices);
      T valueExpanded = itExpanded.getValue();

      // Check for infinite or extreme values that might cause numerical
      // issues
      if (!itExpanded.isDefined() || std::isinf(valueExpanded) ||
          std::isinf(valueIterated)) {
        numSkippedPoints++;
        itIterated.next();
        continue;
      }

      // Calculate difference and add to sum
      T diff = std::abs(valueExpanded - valueIterated) * gridDelta;
      T diffSquared = diff * diff;
      sumDifferences += diff;
      sumSquaredDifferences += diffSquared;
      numPoints++;

      // Store difference in mesh if required
      if (generateMesh) {
        // Create a new point with the coordinates of the iterated level set
        // point
        Vec3D<T> coords{xCoord, yCoord, zCoord}; // lsMesh needs 3D

        // Store the coordinates
        nodeCoordinates.push_back(coords);

        // Store the difference value (squared and absolute)
        differenceValues.push_back(diff);
        squaredDifferenceValues.push_back(diffSquared);

        // Create a vertex for this point
        std::array<unsigned, 1> vertex = {
            static_cast<unsigned>(nodeCoordinates.size() - 1)};
        vertexIndices.push_back(vertex);

        // Update the mesh extent
        for (unsigned i = 0; i < D; ++i) {
          outputMesh->minimumExtent[i] =
              std::min(outputMesh->minimumExtent[i], coords[i]);
          outputMesh->maximumExtent[i] =
              std::max(outputMesh->maximumExtent[i], coords[i]);
        }
      }

      if (fillIteratedWithDistances) {
        // Store the distance value in the point data
        pointDataDistances.push_back(diff);
      }
      // Move to next point
      itIterated.next();
    }

    // Finalize mesh output
    if (generateMesh && !nodeCoordinates.empty()) {
      // Store the node coordinates directly
      outputMesh->nodes = std::move(nodeCoordinates);

      // Store the vertices for point visualization
      outputMesh->vertices = std::move(vertexIndices);

      // Store the difference values as point data
      outputMesh->pointData.insertNextScalarData(std::move(differenceValues),
                                                 "Absolute differences");
      outputMesh->pointData.insertNextScalarData(
          std::move(squaredDifferenceValues), "Squared differences");
    }

    if (fillIteratedWithDistances) {
      iteratedPointData.insertNextScalarData(std::move(pointDataDistances),
                                             "DistanceToExpanded");
    }
  }

  /// Return the sum of squared differences calculated by apply().
  T getSumSquaredDifferences() const { return sumSquaredDifferences; }

  /// Return the sum of differences calculated by apply().
  T getSumDifferences() const { return sumDifferences; }

  /// Return the number of points used in the comparison.
  unsigned getNumPoints() const { return numPoints; }

  /// Return the number of skipped points during the comparison.
  unsigned getNumSkippedPoints() const { return numSkippedPoints; }

  /// Calculate the root mean square error from previously computed values.
  T getRMSE() const {
    return (numPoints > 0) ? std::sqrt(sumSquaredDifferences / numPoints)
                           : std::numeric_limits<T>::infinity();
  }
};

// Add template specializations for this class
PRECOMPILE_PRECISION_DIMENSION(CompareSparseField)

} // namespace viennals
