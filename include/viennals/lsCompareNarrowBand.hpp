#pragma once

#include <cmath>
#include <hrleDenseCellIterator.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsPreCompileMacros.hpp>

namespace viennals {

using namespace viennacore;

/// Calculate distance measure between two level sets by comparing their SDF
/// values on a narrow band. Returns the sum of squared differences between
/// corresponding grid points.
template <class T, int D> class CompareNarrowBand {
  SmartPointer<Domain<T, D>> levelSetTarget = nullptr;
  SmartPointer<Domain<T, D>> levelSetSample = nullptr;
  hrleVectorType<hrleIndexType, D> minIndex, maxIndex;

  // Variables for x and y range restrictions
  T xRangeMin = std::numeric_limits<T>::lowest();
  T xRangeMax = std::numeric_limits<T>::max();
  T yRangeMin = std::numeric_limits<T>::lowest();
  T yRangeMax = std::numeric_limits<T>::max();
  bool useXRange = false;
  bool useYRange = false;

  bool checkAndCalculateBounds() {
    if (levelSetTarget == nullptr || levelSetSample == nullptr) {
      Logger::getInstance()
          .addWarning("Missing level set in CompareNarrowBand.")
          .print();
      return false;
    }

    // Check if the grids are compatible
    const auto &gridTarget = levelSetTarget->getGrid();
    const auto &gridSample = levelSetSample->getGrid();

    if (gridTarget.getGridDelta() != gridSample.getGridDelta()) {
      Logger::getInstance()
          .addWarning("Grid delta mismatch in CompareNarrowBand. The grid "
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

    // Expand the sample level set using lsExpand to a default width of 5
    if (levelSetSample->getLevelSetWidth() < 5) {
      Logger::getInstance()
          .addWarning("Sample level set width is insufficient. Expanding it to "
                      "a width of 5.")
          .print();
      Expand<T, D> expander;
      expander.setLevelSet(levelSetSample);
      expander.setWidth(5);
      expander.apply();
    }

    // Check if target level set width is sufficient
    if (levelSetTarget->getLevelSetWidth() <
        levelSetSample->getLevelSetWidth() + 50) {
      Logger::getInstance()
          .addWarning(
              "Target level set width is insufficient. It must exceed sample "
              "width by least 50. \n"
              " CORRECTION: The expansion was performed. \n"
              "ALTERNATIVE: Alternatively, please expand the target yourself "
              "using lsExpand before passing it to this function. \n")
          .print();
      Expand<T, D> expander;
      expander.setLevelSet(levelSetTarget);
      expander.setWidth(levelSetSample->getLevelSetWidth() + 50);
      expander.apply();
    }

    // Initialize min and max indices
    for (unsigned i = 0; i < D; ++i) {
      minIndex[i] = std::numeric_limits<hrleIndexType>::max();
      maxIndex[i] = std::numeric_limits<hrleIndexType>::lowest();
    }

    // Calculate actual bounds
    for (unsigned i = 0; i < D; ++i) {
      minIndex[i] = std::min({minIndex[i],
                              (gridTarget.isNegBoundaryInfinite(i))
                                  ? domainTarget.getMinRunBreak(i)
                                  : gridTarget.getMinIndex(i),
                              (gridSample.isNegBoundaryInfinite(i))
                                  ? domainSample.getMinRunBreak(i)
                                  : gridSample.getMinIndex(i)});

      maxIndex[i] = std::max({maxIndex[i],
                              (gridTarget.isPosBoundaryInfinite(i))
                                  ? domainTarget.getMaxRunBreak(i)
                                  : gridTarget.getMaxIndex(i),
                              (gridSample.isPosBoundaryInfinite(i))
                                  ? domainSample.getMaxRunBreak(i)
                                  : gridSample.getMaxIndex(i)});
    }

    return true;
  }

public:
  CompareNarrowBand() {}

  CompareNarrowBand(SmartPointer<Domain<T, D>> passedLevelSetTarget,
                    SmartPointer<Domain<T, D>> passedlevelSetSample)
      : levelSetTarget(passedLevelSetTarget),
        levelSetSample(passedlevelSetSample) {}

  void setLevelSetTarget(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetTarget = passedLevelSet;
  }

  void setlevelSetSample(SmartPointer<Domain<T, D>> passedLevelSet) {
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

  // Fields to store the calculation results
  T sumSquaredDifferences = 0.0;
  unsigned numPoints = 0;

  // Add mesh output capability
  SmartPointer<Mesh<T>> outputMesh = nullptr;
  bool outputSquaredDifferences = false;

  /// Set the output mesh where difference values will be stored
  void setOutputMesh(SmartPointer<Mesh<T>> passedMesh) {
    outputMesh = passedMesh;
  }

  /// Set whether to output squared differences (true) or raw differences
  /// (false)
  void setOutputSquaredDifferences(bool value) {
    outputSquaredDifferences = value;
  }

  /// Apply the comparison and calculate the sum of squared differences.
  void apply() {
    // Perform compatibility checks and calculate bounds
    if (!checkAndCalculateBounds()) {
      // If checks fail, return NaN
      sumSquaredDifferences = std::numeric_limits<T>::quiet_NaN();
      numPoints = 0;
      return;
    }

    const auto &gridTarget = levelSetTarget->getGrid();
    double gridDelta = gridTarget.getGridDelta();

    // Set up iterators for both level sets
    hrleConstDenseCellIterator<typename Domain<T, D>::DomainType> itSample(
        levelSetSample->getDomain(), minIndex);
    hrleConstDenseCellIterator<typename Domain<T, D>::DomainType> itTarget(
        levelSetTarget->getDomain(), minIndex);

    sumSquaredDifferences = 0.0;
    numPoints = 0;

    // Prepare mesh output if needed
    std::unordered_map<hrleVectorType<hrleIndexType, D>, size_t,
                       typename hrleVectorType<hrleIndexType, D>::hash>
        pointIdMapping;
    std::vector<T> differenceValues;
    size_t currentPointId = 0;

    if (outputMesh != nullptr) {
      outputMesh->clear();

      // Initialize mesh extent
      for (unsigned i = 0; i < D; ++i) {
        outputMesh->minimumExtent[i] = std::numeric_limits<T>::max();
        outputMesh->maximumExtent[i] = std::numeric_limits<T>::lowest();
      }

      // Prepare scalar data field for differences
      outputMesh->cellData.insertNextScalarData(
          typename PointData<T>::ScalarDataType(),
          outputSquaredDifferences ? "SquaredDifferences" : "Differences");
    }

    // Iterate through the domain defined by the bounding box
    for (; itSample.getIndices() < maxIndex; itSample.next()) {
      // Check if current point is within specified x and y ranges
      T xCoord = itSample.getIndices()[0] * gridDelta;
      T yCoord = (D > 1) ? itSample.getIndices()[1] * gridDelta : 0;

      // Skip if outside the specified x-range
      if (useXRange && (xCoord < xRangeMin || xCoord > xRangeMax)) {
        continue;
      }

      // Skip if outside the specified y-range (only check in 2D and 3D)
      if (D > 1 && useYRange && (yCoord < yRangeMin || yCoord > yRangeMax)) {
        continue;
      }

      // Move the second iterator to the same position
      itTarget.goToIndicesSequential(itSample.getIndices());

      // Get values at current position
      T valueTarget = 0.0;
      T valueSample = 0.0;

      // Calculate average value at cell center
      for (int i = 0; i < (1 << D); ++i) {
        valueSample += itSample.getCorner(i).getValue();
        valueTarget += itTarget.getCorner(i).getValue();
      }
      valueTarget /= (1 << D);
      valueSample /= (1 << D);

      // Check for infinite or extreme values that might cause numerical issues
      if (std::isinf(valueTarget) || std::isinf(valueSample) ||
          std::abs(valueTarget) > 1000 || std::abs(valueSample) > 1000) {
        continue;
      }

      // Calculate difference and add to sum
      T diff = (valueTarget - valueSample) * gridDelta;
      T diffSquared = diff * diff;
      sumSquaredDifferences += diffSquared;
      numPoints++;

      // Store difference in mesh if required
      if (outputMesh != nullptr) {
        std::array<unsigned, 1 << D> voxel;
        bool addVoxel = true;

        // Insert all points of voxel into pointList
        for (unsigned i = 0; i < (1 << D); ++i) {
          hrleVectorType<hrleIndexType, D> index;
          for (unsigned j = 0; j < D; ++j) {
            index[j] =
                itSample.getIndices(j) + itSample.getCorner(i).getOffset()[j];
            if (index[j] > maxIndex[j]) {
              addVoxel = false;
              break;
            }
          }
          if (addVoxel) {
            auto pointIdValue = std::make_pair(index, currentPointId);
            auto pointIdPair = pointIdMapping.insert(pointIdValue);
            voxel[i] = pointIdPair.first->second;
            if (pointIdPair.second) {
              ++currentPointId;
            }
          } else {
            break;
          }
        }

        if (addVoxel) {
          if constexpr (D == 3) {
            std::array<unsigned, 8> hexa{voxel[0], voxel[1], voxel[3],
                                         voxel[2], voxel[4], voxel[5],
                                         voxel[7], voxel[6]};
            outputMesh->hexas.push_back(hexa);
          } else if constexpr (D == 2) {
            std::array<unsigned, 4> quad{voxel[0], voxel[1], voxel[3],
                                         voxel[2]};
            outputMesh->tetras.push_back(quad);
          } else if constexpr (D == 1) {
            std::array<unsigned, 2> line{voxel[0], voxel[1]};
            outputMesh->lines.push_back(line);
          }

          // Add difference value to cell data
          outputMesh->cellData.getScalarData(0)->push_back(
              outputSquaredDifferences ? diffSquared : diff);
        }
      }
    }

    // Insert points into mesh if needed
    if (outputMesh != nullptr && !pointIdMapping.empty()) {
      outputMesh->nodes.resize(pointIdMapping.size());
      for (auto it = pointIdMapping.begin(); it != pointIdMapping.end(); ++it) {
        std::array<T, 3> coords{};
        for (unsigned i = 0; i < D; ++i) {
          coords[i] = gridDelta * it->first[i];

          if (coords[i] < outputMesh->minimumExtent[i]) {
            outputMesh->minimumExtent[i] = coords[i];
          } else if (coords[i] > outputMesh->maximumExtent[i]) {
            outputMesh->maximumExtent[i] = coords[i];
          }
        }
        outputMesh->nodes[it->second] = coords;
      }
    }
  }

  /// Apply with mesh output - convenience method that ensures the mesh is
  /// generated
  void applyWithMeshOutput(SmartPointer<Mesh<T>> outputMesh,
                           bool outputSquared = true) {
    setOutputMesh(outputMesh);
    setOutputSquaredDifferences(outputSquared);
    apply();
  }

  /// Return the sum of squared differences calculated by apply().
  T getSumSquaredDifferences() const { return sumSquaredDifferences; }

  /// Return the number of points used in the comparison.
  unsigned getNumPoints() const { return numPoints; }

  /// Calculate the root mean square error from previously computed values.
  T getRMSE() const {
    return (numPoints > 0) ? std::sqrt(sumSquaredDifferences / numPoints)
                           : std::numeric_limits<T>::infinity();
  }
};

// Add template specializations for this class
PRECOMPILE_PRECISION_DIMENSION(CompareNarrowBand)

} // namespace viennals
