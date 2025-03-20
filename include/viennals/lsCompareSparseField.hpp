#pragma once

#include <hrleDenseIterator.hpp>
#include <hrleSparseIterator.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMesh.hpp>
#include <lsPreCompileMacros.hpp>
#include <lsReduce.hpp>

#include <cmath>
#include <limits>
#include <unordered_map>

namespace viennals {

using namespace viennacore;

/// Calculate distance measure between two level sets by comparing their SDF
/// values on a sparse field. This class iterates over the points in the sparse
/// field of the sample level set and calculates differences with the target
/// level set. The target level set needs to be expanded whereas the sample
/// level set reduced to a sparse field if necessary. The code is currently
/// intended for 2D level sets only.
template <class T, int D = 2> class CompareSparseField {
  SmartPointer<Domain<T, D>> levelSetTarget = nullptr;
  SmartPointer<Domain<T, D>> levelSetSample = nullptr;

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

  // Add mesh output capability
  SmartPointer<Mesh<T>> outputMesh = nullptr;

  bool checkAndCalculateBounds() {
    if (levelSetTarget == nullptr || levelSetSample == nullptr) {
      Logger::getInstance()
          .addWarning("Missing level set in CompareSparseField.")
          .print();
      return false;
    }

    // Check if the grids are compatible
    const auto &gridTarget = levelSetTarget->getGrid();
    const auto &gridSample = levelSetSample->getGrid();

    if (gridTarget.getGridDelta() != gridSample.getGridDelta()) {
      Logger::getInstance()
          .addWarning("Grid delta mismatch in CompareSparseField. The grid "
                      "deltas of the two level sets must be equal.")
          .print();
      return false;
    }

    // Check if the x extents of both level sets are equal
    const auto &domainTarget = levelSetTarget->getDomain();
    const auto &domainSample = levelSetSample->getDomain();

    hrleIndexType targetMinX = gridTarget.isNegBoundaryInfinite(0)
                                   ? domainTarget.getMinRunBreak(0)
                                   : gridTarget.getMinIndex(0);
    hrleIndexType targetMaxX = gridTarget.isPosBoundaryInfinite(0)
                                   ? domainTarget.getMaxRunBreak(0)
                                   : gridTarget.getMaxIndex(0);
    hrleIndexType sampleMinX = gridSample.isNegBoundaryInfinite(0)
                                   ? domainSample.getMinRunBreak(0)
                                   : gridSample.getMinIndex(0);
    hrleIndexType sampleMaxX = gridSample.isPosBoundaryInfinite(0)
                                   ? domainSample.getMaxRunBreak(0)
                                   : gridSample.getMaxIndex(0);

    if (targetMinX != sampleMinX || targetMaxX != sampleMaxX) {
      Logger::getInstance()
          .addWarning("X extent mismatch in CompareNarrowBand. The x extents "
                      "of both level sets must be equal.")
          .print();
      return false;
    }

    // Check if target level set width is sufficient
    if (levelSetTarget->getLevelSetWidth() < 50) {
      Logger::getInstance()
          .addWarning(
              "Target level set width is insufficient. It must have a width of "
              "at least 50. \n"
              " CORRECTION: The expansion was performed. \n"
              "ALTERNATIVE: Alternatively, please expand the target yourself "
              "using lsExpand before passing it to this function. \n")
          .print();
      Expand<T, D>(levelSetTarget, 50).apply();
    }

    // Reduce the sample level set to a sparse field if necessary
    if (levelSetSample->getLevelSetWidth() > 1) {
      Logger::getInstance()
          .addWarning(
              "Sample level set width is too large. It must be reduced to a "
              "sparse field. \n"
              " CORRECTION: The reduction was performed. \n"
              "ALTERNATIVE: Alternatively, please reduce the sample yourself "
              "using lsReduce before passing it to this function. \n")
          .print();
      Reduce<T, D>(levelSetSample, 1).apply();
    }

    return true;
  }

public:
  CompareSparseField() {
    static_assert(
        D == 2 &&
        "CompareSparseField is currently only implemented for 2D level sets.");
  }

  CompareSparseField(SmartPointer<Domain<T, D>> passedLevelSetTarget,
                     SmartPointer<Domain<T, D>> passedLevelSetSample)
      : levelSetTarget(passedLevelSetTarget),
        levelSetSample(passedLevelSetSample) {
    static_assert(
        D == 2 &&
        "CompareSparseField is currently only implemented for 2D level sets.");
  }

  void setLevelSetTarget(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetTarget = passedLevelSet;
  }

  void setLevelSetSample(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetSample = passedLevelSet;
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

  /// Apply the comparison and calculate the sum of squared differences.
  void apply() {
    // Perform compatibility checks
    if (!checkAndCalculateBounds()) {
      // If checks fail, return NaN
      sumSquaredDifferences = std::numeric_limits<T>::quiet_NaN();
      sumDifferences = std::numeric_limits<T>::quiet_NaN();
      numPoints = 0;
      return;
    }

    const auto &gridTarget = levelSetTarget->getGrid();
    const auto gridDelta = gridTarget.getGridDelta();

    sumSquaredDifferences = 0.0;
    sumDifferences = 0.0;
    numPoints = 0;

    // Direct storage for points and differences
    std::vector<std::array<T, 3>> nodeCoordinates; // 3D necessary for lsMesh
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
      nodeCoordinates.reserve(levelSetSample->getNumberOfPoints());
      vertexIndices.reserve(levelSetSample->getNumberOfPoints());
      differenceValues.reserve(levelSetSample->getNumberOfPoints());
      squaredDifferenceValues.reserve(levelSetSample->getNumberOfPoints());
    }

    // Create sparse iterators for the level sets
    hrleConstSparseIterator<typename Domain<T, D>::DomainType> itSample(
        levelSetSample->getDomain());
    hrleConstSparseIterator<typename Domain<T, D>::DomainType> itTarget(
        levelSetTarget->getDomain());

    // Iterate over all defined points in the sample level set
    while (!itSample.isFinished()) {
      if (!itSample.isDefined()) {
        // this block is necessary to skip undefined points, I tested it:
        // std::cout << "Skipping undefined point" << std::endl;
        itSample.next();
        continue;
      }

      auto indices = itSample.getStartIndices();

      // Calculate coordinates
      T xCoord = indices[0] * gridDelta;
      T yCoord = indices[1] * gridDelta;
      T zCoord = 0.0; // Always use 0 for z-coordinate in 2D

      // Skip if outside the specified x-range
      if (useXRange && (xCoord < xRangeMin || xCoord > xRangeMax)) {
        itSample.next();
        continue;
      }

      // Skip if outside the specified y-range
      if (useYRange && (yCoord < yRangeMin || yCoord > yRangeMax)) {
        itSample.next();
        continue;
      }

      // Get sample value
      T valueSample = itSample.getValue();

      itTarget.goToIndicesSequential(indices);
      T valueTarget = itTarget.getValue();

      // Check for infinite or extreme values that might cause numerical issues
      if (!itTarget.isDefined() || std::isinf(valueTarget) ||
          std::isinf(valueSample)) {
        itSample.next();
        continue;
      }

      // Calculate difference and add to sum
      T diff = std::abs(valueTarget - valueSample);
      T diffSquared = diff * diff;
      sumDifferences += diff;
      sumSquaredDifferences += diffSquared;
      numPoints++;

      // Store difference in mesh if required
      if (generateMesh) {
        // Create a new point with the coordinates of the sample level set point
        std::array<T, 3> coords = {xCoord, yCoord, zCoord}; // lsMesh needs 3D

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

      // Move to next point
      itSample.next();
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
  }

  /// Return the sum of squared differences calculated by apply().
  T getSumSquaredDifferences() const { return sumSquaredDifferences; }

  /// Return the sum of differences calculated by apply().
  T getSumDifferences() const { return sumDifferences; }

  /// Return the number of points used in the comparison.
  unsigned getNumPoints() const { return numPoints; }

  /// Calculate the root mean square error from previously computed values.
  T getRMSE() const {
    return (numPoints > 0) ? std::sqrt(sumSquaredDifferences / numPoints)
                           : std::numeric_limits<T>::infinity();
  }
};

// Add template specializations for this class
PRECOMPILE_PRECISION_DIMENSION(CompareSparseField)

} // namespace viennals
