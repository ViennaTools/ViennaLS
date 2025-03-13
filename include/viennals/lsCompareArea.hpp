#pragma once

#include <limits>
#include <unordered_map>

#include <hrleDenseCellIterator.hpp>
#include <lsDomain.hpp>
#include <lsMesh.hpp>
#include <lsPreCompileMacros.hpp>

namespace viennals {

using namespace viennacore;

/// Computes an estimate of the area where two level sets differ.
/// The area is calculated by iterating through the bounding box of the two
/// level sets and comparing the cell values. The grid delta is used as the unit
/// of area. Custom increment values can be set for specific x and y ranges,
/// allowing to count certain areas multiple times or skip them. Optionally, a
/// passed mesh can be filled with the area information, allowing for
/// visualization of the differences.
/// The code is currently itended for 2D level sets only.
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

  // Mesh output related members
  SmartPointer<Mesh<T>> outputMesh = nullptr;
  bool generateMesh = false;
  std::unordered_map<hrleVectorType<hrleIndexType, D>, size_t,
                     typename hrleVectorType<hrleIndexType, D>::hash>
      pointIdMapping;
  size_t currentPointId = 0;

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
          .addError("Grid delta mismatch in CompareArea. The grid deltas of "
                    "the two level sets must be equal.")
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
                              (gridTarget.isNegBoundaryInfinite(i))
                                  ? domainTarget.getMinRunBreak(i)
                                  : gridTarget.getMinBounds(i),
                              (gridSample.isNegBoundaryInfinite(i))
                                  ? domainSample.getMinRunBreak(i)
                                  : gridSample.getMinBounds(i)});

      maxIndex[i] = std::max({maxIndex[i],
                              (gridTarget.isPosBoundaryInfinite(i))
                                  ? domainTarget.getMaxRunBreak(i)
                                  : gridTarget.getMaxBounds(i),
                              (gridSample.isPosBoundaryInfinite(i))
                                  ? domainSample.getMaxRunBreak(i)
                                  : gridSample.getMaxBounds(i)});
    }

    return true;
  }

public:
  CompareArea() {}

  CompareArea(SmartPointer<Domain<T, D>> passedLevelSetTarget,
              SmartPointer<Domain<T, D>> passedLevelSetSample)
      : levelSetTarget(passedLevelSetTarget),
        levelSetSample(passedLevelSetSample) {}

  /// Sets the target level set.
  void setLevelSetTarget(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetTarget = passedLevelSet;
  }

  /// Sets the sample level set.
  void setLevelSetSample(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetSample = passedLevelSet;
  }

  /// Set default increment value
  void setDefaultIncrement(unsigned short int increment) {
    defaultIncrement = increment;
  }

  /// Sets the x-range and custom increment value
  void setXRangeAndIncrement(hrleIndexType minXRange, hrleIndexType maxXRange,
                             unsigned short int Xincrement) {
    xRangeMin = minXRange;
    xRangeMax = maxXRange;
    customXIncrement = Xincrement;
    useCustomXIncrement = true;
  }

  /// Sets the y-range and custom increment value
  void setYRangeAndIncrement(hrleIndexType minYRange, hrleIndexType maxYRange,
                             unsigned short int Yincrement) {
    yRangeMin = minYRange;
    yRangeMax = maxYRange;
    customYIncrement = Yincrement;
    useCustomYIncrement = true;
  }

  /// Set the output mesh where difference areas will be stored
  void setOutputMesh(SmartPointer<Mesh<T>> passedMesh) {
    outputMesh = passedMesh;
  }

  /// Apply the comparison and generate a mesh visualization of the differences.
  /// Each cell in the mesh will have a material ID:
  ///   0: Areas where both level sets are inside
  ///   1: Areas where only one level set is inside (mismatched areas)
  void applyWithMeshOutput(SmartPointer<Mesh<T>> passedMesh) {
    if (passedMesh == nullptr) {
      Logger::getInstance()
          .addWarning("No mesh was passed to CompareArea::applyWithMeshOutput.")
          .print();
      return;
    }

    outputMesh = passedMesh;
    outputMesh->clear();
    generateMesh = true;
    apply();
    generateMesh = false;
  }

  /// Returns the computed area mismatch.
  double getAreaMismatch() const {
    return static_cast<double>(cellCount) * gridDelta * gridDelta;
  }

  /// Returns the computed area mismatch, with custom increments applied.
  double getCustomAreaMismatch() const {
    return static_cast<double>(customCellCount) * gridDelta * gridDelta;
  }

  /// Returns the number of cells where the level sets differ.
  unsigned long int getCellCount() const { return cellCount; }

  /// Returns the number of cells where the level sets differ, with custom
  /// increments applied.
  unsigned long int getCustomCellCount() const { return customCellCount; }

  /// Computes the area difference between the two level sets.
  void apply() {
    // Calculate the bounds for iteration
    if (!checkAndCalculateBounds()) {
      cellCount =
          std::numeric_limits<unsigned long int>::max(); // Error indicator
      customCellCount =
          std::numeric_limits<unsigned long int>::max(); // Error indicator
      return;
    }

    // Set up dense cell iterators for both level sets
    hrleConstDenseCellIterator<typename Domain<T, D>::DomainType> itTarget(
        levelSetTarget->getDomain(), minIndex);
    hrleConstDenseCellIterator<typename Domain<T, D>::DomainType> itSample(
        levelSetSample->getDomain(), minIndex);

    cellCount = 0;
    customCellCount = 0;

    // Initialize mesh-related variables if generating a mesh
    if (generateMesh && outputMesh != nullptr) {
      // Save the extent of the resulting mesh
      for (unsigned i = 0; i < D; ++i) {
        outputMesh->minimumExtent[i] = std::numeric_limits<T>::max();
        outputMesh->maximumExtent[i] = std::numeric_limits<T>::lowest();
      }

      pointIdMapping.clear();
      currentPointId = 0;

      // Prepare mesh for material IDs and custom values
      outputMesh->cellData.insertNextScalarData(
          typename PointData<T>::ScalarDataType(), "Material");
      auto &materialIds = *(outputMesh->cellData.getScalarData(0));

      // Add a second scalar field for custom increment values
      outputMesh->cellData.insertNextScalarData(
          typename PointData<T>::ScalarDataType(), "CustomIncrement");
      auto &incrementValues = *(outputMesh->cellData.getScalarData(1));
    }

    // Reference to the material IDs and increment values
    auto materialIdPtr =
        generateMesh ? &(*(outputMesh->cellData.getScalarData(0))) : nullptr;
    auto incrementValuesPtr =
        generateMesh ? &(*(outputMesh->cellData.getScalarData(1))) : nullptr;

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

      // Determine if cell is inside for each level set
      bool insideTarget = (centerValueTarget <= 0.);
      bool insideSample = (centerValueSample <= 0.);
      bool isDifferent = insideTarget != insideSample;

      // Calculate range checks once per cell
      bool inXRange = useCustomXIncrement &&
                      (itTarget.getIndices()[0] * gridDelta >= xRangeMin &&
                       itTarget.getIndices()[0] * gridDelta <= xRangeMax);
      bool inYRange = useCustomYIncrement &&
                      (itTarget.getIndices()[1] * gridDelta >= yRangeMin &&
                       itTarget.getIndices()[1] * gridDelta <= yRangeMax);

      // Calculate increment to add based on ranges
      unsigned short int incrementToAdd = defaultIncrement;
      if (inXRange && inYRange) {
        // Apply both increments
        incrementToAdd = customXIncrement + customYIncrement;
      } else if (inXRange) {
        incrementToAdd = customXIncrement;
      } else if (inYRange) {
        incrementToAdd = customYIncrement;
      }

      // If cells differ, update the counters
      if (isDifferent) {
        // Always increment simple cell count by 1
        cellCount += 1;

        // For custom cell count, apply the calculated increment
        customCellCount += incrementToAdd;
      }

      // For mesh generation, process all cells where at least one level set is
      // inside
      if (generateMesh && (insideTarget || insideSample)) {
        // Material ID: 0 for matching inside, 1 for mismatched
        int materialId = isDifferent ? 1 : 0;

        std::array<unsigned, 1 << D> voxel;
        bool addVoxel = true;

        // Insert all points of voxel into pointList
        for (unsigned i = 0; i < (1 << D); ++i) {
          hrleVectorType<hrleIndexType, D> index;
          for (unsigned j = 0; j < D; ++j) {
            index[j] =
                itTarget.getIndices(j) + itTarget.getCorner(i).getOffset()[j];
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

          materialIdPtr->push_back(materialId);

          // Store increment value (0 if not different)
          incrementValuesPtr->push_back(
              isDifferent ? static_cast<T>(incrementToAdd) : 0);
        }
      }
    }

    // Finalize mesh if generating one
    if (generateMesh && outputMesh != nullptr && !pointIdMapping.empty()) {
      // Insert points into the mesh
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

    // Calculate the mismatched area based on the cell count
    mismatchedArea = static_cast<double>(cellCount) * gridDelta * gridDelta;
    customMismatchedArea =
        static_cast<double>(customCellCount) * gridDelta * gridDelta;
  }
};

// Precompile for common precision and dimensions
PRECOMPILE_PRECISION_DIMENSION(CompareArea)

} // namespace viennals
