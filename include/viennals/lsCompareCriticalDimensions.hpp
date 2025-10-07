#pragma once

#include <hrleSparseIterator.hpp>
#include <lsDomain.hpp>
#include <lsMesh.hpp>
#include <lsPreCompileMacros.hpp>
#include <lsToSurfaceMesh.hpp>

#include <cmath>
#include <limits>
#include <vector>

namespace viennals {

using namespace viennacore;

/// Compares critical dimensions (surface positions) between two level sets.
/// Critical dimensions are defined as the maximum or minimum positions where
/// the surface (SDF = 0) exists within a specified range.
///
/// - If X range is specified: finds Y coordinates where surface exists, then
///   identifies max/min Y positions
/// - If Y range is specified: finds X coordinates where surface exists, then
///   identifies max/min X positions
///
/// The surface position is interpolated from grid points where the SDF crosses
/// zero. Multiple ranges can be specified to compare different critical
/// dimensions.
///
/// The code is currently intended for 2D level sets only.

template <class T, int D = 2> class CompareCriticalDimensions {
  using hrleIndexType = viennahrle::IndexType;

  SmartPointer<Domain<T, D>> levelSetReference = nullptr;
  SmartPointer<Domain<T, D>> levelSetCompare = nullptr;

  // Structure to hold a range specification
  struct RangeSpec {
    bool isXRange; // true if X range, false if Y range
    T rangeMin;
    T rangeMax;
    bool findMaximum; // true for maximum, false for minimum
  };

  std::vector<RangeSpec> rangeSpecs;

  // Structure to hold critical dimension results
  struct CriticalDimensionResult {
    bool isXRange;
    T rangeMin;
    T rangeMax;
    bool findMaximum;
    T positionReference; // Critical dimension position in reference LS
    T positionCompare;   // Critical dimension position in compare LS
    T difference;        // Absolute difference
    bool valid;          // Whether both critical dimensions were found
  };

  std::vector<CriticalDimensionResult> results;

  // Optional mesh output
  SmartPointer<Mesh<T>> outputMesh = nullptr;

  bool checkInputs() {
    if (levelSetReference == nullptr || levelSetCompare == nullptr) {
      Logger::getInstance()
          .addWarning("Missing level set in CompareCriticalDimensions.")
          .print();
      return false;
    }

    // Check if the grids are compatible
    const auto &gridReference = levelSetReference->getGrid();
    const auto &gridCompare = levelSetCompare->getGrid();

    if (gridReference.getGridDelta() != gridCompare.getGridDelta()) {
      Logger::getInstance()
          .addWarning("Grid delta mismatch in CompareCriticalDimensions. The "
                      "grid deltas of the two level sets must be equal.")
          .print();
      return false;
    }

    if (rangeSpecs.empty()) {
      Logger::getInstance()
          .addWarning("No ranges specified in CompareCriticalDimensions.")
          .print();
      return false;
    }

    return true;
  }

  // Extract surface positions from mesh nodes within the specified range
  // Returns the position coordinates (Y if isXRange=true, X if isXRange=false)
  std::vector<T> findSurfaceCrossings(SmartPointer<Mesh<T>> surfaceMesh,
                                      bool isXRange, T scanMin, T scanMax,
                                      T perpMin, T perpMax) {
    std::vector<T> crossings;

    // Iterate through all surface mesh nodes
    for (const auto &node : surfaceMesh->nodes) {
      T xCoord = node[0];
      T yCoord = node[1];

      // Check if point is in our scan range
      bool inRange = false;
      if (isXRange) {
        // X range specified - check if X is in scan range
        inRange = (xCoord >= scanMin && xCoord <= scanMax);
      } else {
        // Y range specified - check if Y is in scan range
        inRange = (yCoord >= scanMin && yCoord <= scanMax);
      }

      if (inRange) {
        // Extract perpendicular coordinate
        T perpCoord = isXRange ? yCoord : xCoord;

        // Check if perpendicular coordinate is in range
        if (perpCoord >= perpMin && perpCoord <= perpMax) {
          crossings.push_back(perpCoord);
        }
      }
    }

    return crossings;
  }

  // Find critical dimension (max or min) from a list of surface crossings
  std::pair<bool, T> findCriticalDimension(const std::vector<T> &crossings,
                                           bool findMaximum) {
    if (crossings.empty()) {
      return {false, 0.0};
    }

    if (findMaximum) {
      return {true, *std::max_element(crossings.begin(), crossings.end())};
    } else {
      return {true, *std::min_element(crossings.begin(), crossings.end())};
    }
  }

public:
  CompareCriticalDimensions() {
    static_assert(D == 2 && "CompareCriticalDimensions is currently only "
                            "implemented for 2D level sets.");
  }

  CompareCriticalDimensions(SmartPointer<Domain<T, D>> passedLevelSetReference,
                            SmartPointer<Domain<T, D>> passedLevelSetCompare)
      : levelSetReference(passedLevelSetReference),
        levelSetCompare(passedLevelSetCompare) {
    static_assert(D == 2 && "CompareCriticalDimensions is currently only "
                            "implemented for 2D level sets.");
  }

  void setLevelSetReference(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetReference = passedLevelSet;
  }

  void setLevelSetCompare(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetCompare = passedLevelSet;
  }

  /// Add an X range to find maximum or minimum Y position
  void addXRange(T minX, T maxX, bool findMaximum = true) {
    RangeSpec spec;
    spec.isXRange = true;
    spec.rangeMin = minX;
    spec.rangeMax = maxX;
    spec.findMaximum = findMaximum;
    rangeSpecs.push_back(spec);
  }

  /// Add a Y range to find maximum or minimum X position
  void addYRange(T minY, T maxY, bool findMaximum = true) {
    RangeSpec spec;
    spec.isXRange = false;
    spec.rangeMin = minY;
    spec.rangeMax = maxY;
    spec.findMaximum = findMaximum;
    rangeSpecs.push_back(spec);
  }

  /// Clear all range specifications
  void clearRanges() { rangeSpecs.clear(); }

  /// Set the output mesh where critical dimension locations will be stored
  void setOutputMesh(SmartPointer<Mesh<T>> passedMesh) {
    outputMesh = passedMesh;
  }

  /// Apply the comparison
  void apply() {
    results.clear();

    if (!checkInputs()) {
      return;
    }

    const auto &grid = levelSetReference->getGrid();
    const T gridDelta = grid.getGridDelta();

    // Convert both level sets to surface meshes once
    auto surfaceMeshRef = SmartPointer<Mesh<T>>::New();
    auto surfaceMeshCmp = SmartPointer<Mesh<T>>::New();

    ToSurfaceMesh<T, D>(levelSetReference, surfaceMeshRef).apply();
    ToSurfaceMesh<T, D>(levelSetCompare, surfaceMeshCmp).apply();

    // Get the domain extent to determine perpendicular range
    // Use the union of both domains' bounds to ensure all points are considered
    const auto &gridRef = levelSetReference->getGrid();
    const auto &gridCmp = levelSetCompare->getGrid();
    T xMin =
        std::min(gridRef.getMinBounds(0), gridCmp.getMinBounds(0)) * gridDelta;
    T xMax =
        std::max(gridRef.getMaxBounds(0), gridCmp.getMaxBounds(0)) * gridDelta;
    T yMin =
        std::min(gridRef.getMinBounds(1), gridCmp.getMinBounds(1)) * gridDelta;
    T yMax =
        std::max(gridRef.getMaxBounds(1), gridCmp.getMaxBounds(1)) * gridDelta;

    // Process each range specification
    for (const auto &spec : rangeSpecs) {
      CriticalDimensionResult result;
      result.isXRange = spec.isXRange;
      result.rangeMin = spec.rangeMin;
      result.rangeMax = spec.rangeMax;
      result.findMaximum = spec.findMaximum;
      result.valid = false;

      // Determine scan and perpendicular ranges
      T scanMin, scanMax, perpMin, perpMax;
      if (spec.isXRange) {
        // X range specified - scan along X, find Y positions
        scanMin = spec.rangeMin;
        scanMax = spec.rangeMax;
        perpMin = yMin;
        perpMax = yMax;
      } else {
        // Y range specified - scan along Y, find X positions
        scanMin = spec.rangeMin;
        scanMax = spec.rangeMax;
        perpMin = xMin;
        perpMax = xMax;
      }

      // Find all surface crossings from the mesh nodes
      auto crossingsRef = findSurfaceCrossings(
          surfaceMeshRef, spec.isXRange, scanMin, scanMax, perpMin, perpMax);
      auto crossingsCmp = findSurfaceCrossings(
          surfaceMeshCmp, spec.isXRange, scanMin, scanMax, perpMin, perpMax);

      // Find critical dimensions
      auto [validRef, cdRef] =
          findCriticalDimension(crossingsRef, spec.findMaximum);
      auto [validCmp, cdCmp] =
          findCriticalDimension(crossingsCmp, spec.findMaximum);

      if (validRef && validCmp) {
        result.valid = true;
        result.positionReference = cdRef;
        result.positionCompare = cdCmp;
        result.difference = std::abs(cdRef - cdCmp);
      }

      results.push_back(result);
    }

    // Generate mesh if requested
    if (outputMesh != nullptr) {
      generateMesh();
    }
  }

  /// Get the number of critical dimensions compared
  size_t getNumCriticalDimensions() const { return results.size(); }

  /// Get a specific critical dimension result
  bool getCriticalDimensionResult(size_t index, T &positionReference,
                                  T &positionCompare, T &difference) const {
    if (index >= results.size() || !results[index].valid) {
      return false;
    }
    positionReference = results[index].positionReference;
    positionCompare = results[index].positionCompare;
    difference = results[index].difference;
    return true;
  }

  /// Get mean absolute difference across all valid critical dimensions
  T getMeanDifference() const {
    T sum = 0.0;
    size_t count = 0;
    for (const auto &result : results) {
      if (result.valid) {
        sum += result.difference;
        count++;
      }
    }
    return count > 0 ? sum / count : std::numeric_limits<T>::quiet_NaN();
  }

  /// Get maximum difference across all valid critical dimensions
  T getMaxDifference() const {
    T maxDiff = 0.0;
    for (const auto &result : results) {
      if (result.valid) {
        maxDiff = std::max(maxDiff, result.difference);
      }
    }
    return maxDiff;
  }

  /// Get RMSE across all valid critical dimensions
  T getRMSE() const {
    T sumSquared = 0.0;
    size_t count = 0;
    for (const auto &result : results) {
      if (result.valid) {
        sumSquared += result.difference * result.difference;
        count++;
      }
    }
    return count > 0 ? std::sqrt(sumSquared / count)
                     : std::numeric_limits<T>::quiet_NaN();
  }

  /// Get all valid results
  std::vector<T> getAllDifferences() const {
    std::vector<T> differences;
    for (const auto &result : results) {
      if (result.valid) {
        differences.push_back(result.difference);
      }
    }
    return differences;
  }

private:
  void generateMesh() {
    outputMesh->clear();

    std::vector<Vec3D<T>> nodeCoordinates;
    std::vector<std::array<unsigned, 1>> vertexIndices;
    std::vector<T> differenceValues;
    std::vector<T> referenceValues;
    std::vector<T> compareValues;

    for (unsigned i = 0; i < D; ++i) {
      outputMesh->minimumExtent[i] = std::numeric_limits<T>::max();
      outputMesh->maximumExtent[i] = std::numeric_limits<T>::lowest();
    }

    unsigned pointId = 0;
    for (const auto &result : results) {
      if (!result.valid)
        continue;

      // Create points for reference and compare positions
      Vec3D<T> coordRef, coordCmp;

      if (result.isXRange) {
        // Critical dimension is in Y, position is along X range
        T xMid = (result.rangeMin + result.rangeMax) / 2.0;
        coordRef = {xMid, result.positionReference, 0.0};
        coordCmp = {xMid, result.positionCompare, 0.0};
      } else {
        // Critical dimension is in X, position is along Y range
        T yMid = (result.rangeMin + result.rangeMax) / 2.0;
        coordRef = {result.positionReference, yMid, 0.0};
        coordCmp = {result.positionCompare, yMid, 0.0};
      }

      // Add reference point
      nodeCoordinates.push_back(coordRef);
      vertexIndices.push_back({pointId++});
      differenceValues.push_back(result.difference);
      referenceValues.push_back(result.positionReference);
      compareValues.push_back(result.positionCompare);

      // Add compare point
      nodeCoordinates.push_back(coordCmp);
      vertexIndices.push_back({pointId++});
      differenceValues.push_back(result.difference);
      referenceValues.push_back(result.positionReference);
      compareValues.push_back(result.positionCompare);

      // Update extent
      for (unsigned i = 0; i < D; ++i) {
        outputMesh->minimumExtent[i] =
            std::min({outputMesh->minimumExtent[i], coordRef[i], coordCmp[i]});
        outputMesh->maximumExtent[i] =
            std::max({outputMesh->maximumExtent[i], coordRef[i], coordCmp[i]});
      }
    }

    if (!nodeCoordinates.empty()) {
      outputMesh->nodes = std::move(nodeCoordinates);
      outputMesh->vertices = std::move(vertexIndices);
      outputMesh->pointData.insertNextScalarData(std::move(differenceValues),
                                                 "Difference");
      outputMesh->pointData.insertNextScalarData(std::move(referenceValues),
                                                 "ReferencePosition");
      outputMesh->pointData.insertNextScalarData(std::move(compareValues),
                                                 "ComparePosition");
    }
  }
};

// Add template specializations for this class
PRECOMPILE_PRECISION_DIMENSION(CompareCriticalDimensions)

} // namespace viennals
