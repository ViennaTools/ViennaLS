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
/// Note for the future: lsToDiskMesh could be used instead of lsToSurfaceMesh,
/// which is probably more efficient but slightly less accurate.

template <class T, int D = 2> class CompareCriticalDimensions {
  using hrleIndexType = viennahrle::IndexType;

  SmartPointer<Domain<T, D>> levelSetTarget = nullptr;
  SmartPointer<Domain<T, D>> levelSetSample = nullptr;

  // Structure to hold a range specification
  struct RangeSpec {
    int measureDimension;
    std::array<T, D> minBounds;
    std::array<T, D> maxBounds;
    int measureDimension;
    std::array<T, D> minBounds;
    std::array<T, D> maxBounds;
    bool findMaximum; // true for maximum, false for minimum
  };

  std::vector<RangeSpec> rangeSpecs;

  // Structure to hold critical dimension results
  struct CriticalDimensionResult {
    int measureDimension;
    std::array<T, D> minBounds;
    std::array<T, D> maxBounds;
    int measureDimension;
    std::array<T, D> minBounds;
    std::array<T, D> maxBounds;
    bool findMaximum;
    T positionTarget; // Critical dimension position in target LS
    T positionSample; // Critical dimension position in sample LS
    T difference;     // Absolute difference
    bool valid;       // Whether both critical dimensions were found
  };

  std::vector<CriticalDimensionResult> results;

  // Optional mesh output
  SmartPointer<Mesh<T>> outputMesh = nullptr;

  bool checkInputs() {
    if (levelSetTarget == nullptr || levelSetSample == nullptr) {
      VIENNACORE_LOG_WARNING("Missing level set in CompareCriticalDimensions.");
      return false;
    }

    // Check if the grids are compatible
    const auto &gridTarget = levelSetTarget->getGrid();
    const auto &gridSample = levelSetSample->getGrid();

    if (gridTarget.getGridDelta() != gridSample.getGridDelta()) {
      VIENNACORE_LOG_WARNING(
          "Grid delta mismatch in CompareCriticalDimensions. The "
          "grid deltas of the two level sets must be equal.");
      return false;
    }

    if (rangeSpecs.empty()) {
      VIENNACORE_LOG_WARNING(
          "No ranges specified in CompareCriticalDimensions.");
      return false;
    }

    return true;
  }

  // Extract surface positions from mesh nodes within the specified range
  // Returns the position coordinates (Y if isXRange=true, X if isXRange=false)
  std::vector<T> findSurfaceCrossings(SmartPointer<Mesh<T>> surfaceMesh,
                                      const RangeSpec &spec) {
                                      const RangeSpec &spec) {
    std::vector<T> crossings;

    // Iterate through all surface mesh nodes
    for (const auto &node : surfaceMesh->nodes) {
      // Check if point is in our scan range
      bool inRange = true;
      for (int i = 0; i < D; ++i) {
        if (i == spec.measureDimension)
          continue;
        if (node[i] < spec.minBounds[i] || node[i] > spec.maxBounds[i]) {
          inRange = false;
          break;
        }
      bool inRange = true;
      for (int i = 0; i < D; ++i) {
        if (i == spec.measureDimension)
          continue;
        if (node[i] < spec.minBounds[i] || node[i] > spec.maxBounds[i]) {
          inRange = false;
          break;
        }
      }

      if (inRange) {
        T val = node[spec.measureDimension];
        if (val >= spec.minBounds[spec.measureDimension] &&
            val <= spec.maxBounds[spec.measureDimension]) {
          crossings.push_back(val);
        T val = node[spec.measureDimension];
        if (val >= spec.minBounds[spec.measureDimension] &&
            val <= spec.maxBounds[spec.measureDimension]) {
          crossings.push_back(val);
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
  CompareCriticalDimensions() {}
  CompareCriticalDimensions() {}

  CompareCriticalDimensions(SmartPointer<Domain<T, D>> passedLevelSetTarget,
                            SmartPointer<Domain<T, D>> passedLevelSetSample)
      : levelSetTarget(passedLevelSetTarget),
        levelSetSample(passedLevelSetSample) {}
        levelSetSample(passedLevelSetSample) {}

  void setLevelSetTarget(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetTarget = passedLevelSet;
  }

  void setLevelSetSample(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetSample = passedLevelSet;
  }

  /// Add a generic range specification
  void addRange(int measureDimension, const std::array<T, D> &minBounds,
                const std::array<T, D> &maxBounds, bool findMaximum = true) {
  /// Add a generic range specification
  void addRange(int measureDimension, const std::array<T, D> &minBounds,
                const std::array<T, D> &maxBounds, bool findMaximum = true) {
    RangeSpec spec;
    spec.measureDimension = measureDimension;
    spec.minBounds = minBounds;
    spec.maxBounds = maxBounds;
    spec.measureDimension = measureDimension;
    spec.minBounds = minBounds;
    spec.maxBounds = maxBounds;
    spec.findMaximum = findMaximum;
    rangeSpecs.push_back(spec);
  }

  /// Add an X range to find maximum or minimum Y position
  void addXRange(T minX, T maxX, bool findMaximum = true) {
    if constexpr (D == 2) {
      std::array<T, D> minBounds = {minX, std::numeric_limits<T>::lowest()};
      std::array<T, D> maxBounds = {maxX, std::numeric_limits<T>::max()};
      addRange(1, minBounds, maxBounds, findMaximum);
    } else {
      VIENNACORE_LOG_WARNING(
          "addXRange is only supported for 2D. Use addRange for 3D.");
    }
  }

  /// Add an X range to find maximum or minimum Y position
  void addXRange(T minX, T maxX, bool findMaximum = true) {
    if constexpr (D == 2) {
      std::array<T, D> minBounds = {minX, std::numeric_limits<T>::lowest()};
      std::array<T, D> maxBounds = {maxX, std::numeric_limits<T>::max()};
      addRange(1, minBounds, maxBounds, findMaximum);
    } else {
      VIENNACORE_LOG_WARNING(
          "addXRange is only supported for 2D. Use addRange for 3D.");
    }
  }

  /// Add a Y range to find maximum or minimum X position
  void addYRange(T minY, T maxY, bool findMaximum = true) {
    if constexpr (D == 2) {
      std::array<T, D> minBounds = {std::numeric_limits<T>::lowest(), minY};
      std::array<T, D> maxBounds = {std::numeric_limits<T>::max(), maxY};
      addRange(0, minBounds, maxBounds, findMaximum);
    } else {
      VIENNACORE_LOG_WARNING(
          "addYRange is only supported for 2D. Use addRange for 3D.");
    }
    if constexpr (D == 2) {
      std::array<T, D> minBounds = {std::numeric_limits<T>::lowest(), minY};
      std::array<T, D> maxBounds = {std::numeric_limits<T>::max(), maxY};
      addRange(0, minBounds, maxBounds, findMaximum);
    } else {
      VIENNACORE_LOG_WARNING(
          "addYRange is only supported for 2D. Use addRange for 3D.");
    }
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

    const auto &grid = levelSetTarget->getGrid();
    const T gridDelta = grid.getGridDelta();

    // Convert both level sets to surface meshes once
    auto surfaceMeshRef = SmartPointer<Mesh<T>>::New();
    auto surfaceMeshCmp = SmartPointer<Mesh<T>>::New();

    ToSurfaceMesh<T, D>(levelSetTarget, surfaceMeshRef).apply();
    ToSurfaceMesh<T, D>(levelSetSample, surfaceMeshCmp).apply();

    // Pre-allocate results vector to avoid race conditions
    results.resize(rangeSpecs.size());

    // Process each range specification in parallel
#pragma omp parallel for
    for (size_t specIdx = 0; specIdx < rangeSpecs.size(); ++specIdx) {
      const auto &spec = rangeSpecs[specIdx];
      CriticalDimensionResult result;
      result.measureDimension = spec.measureDimension;
      result.minBounds = spec.minBounds;
      result.maxBounds = spec.maxBounds;
      result.measureDimension = spec.measureDimension;
      result.minBounds = spec.minBounds;
      result.maxBounds = spec.maxBounds;
      result.findMaximum = spec.findMaximum;
      result.valid = false;

      // Find all surface crossings from the mesh nodes
      auto crossingsRef = findSurfaceCrossings(surfaceMeshRef, spec);
      auto crossingsCmp = findSurfaceCrossings(surfaceMeshCmp, spec);
      auto crossingsRef = findSurfaceCrossings(surfaceMeshRef, spec);
      auto crossingsCmp = findSurfaceCrossings(surfaceMeshCmp, spec);

      // Find critical dimensions
      auto [validRef, cdRef] =
          findCriticalDimension(crossingsRef, spec.findMaximum);
      auto [validCmp, cdCmp] =
          findCriticalDimension(crossingsCmp, spec.findMaximum);

      if (validRef && validCmp) {
        result.valid = true;
        result.positionTarget = cdRef;
        result.positionSample = cdCmp;
        result.difference = std::abs(cdRef - cdCmp);
      }

      // Direct assignment instead of push_back to avoid race condition
      results[specIdx] = result;
    }

    // Generate mesh if requested
    if (outputMesh != nullptr) {
      generateMesh();
    }
  }

  /// Get the number of critical dimensions compared
  size_t getNumCriticalDimensions() const { return results.size(); }

  /// Get a specific critical dimension result
  bool getCriticalDimensionResult(size_t index, T &positionTarget,
                                  T &positionSample, T &difference) const {
    if (index >= results.size() || !results[index].valid) {
      return false;
    }
    positionTarget = results[index].positionTarget;
    positionSample = results[index].positionSample;
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
    std::vector<T> targetValues;
    std::vector<T> sampleValues;

    for (unsigned i = 0; i < 3; ++i) {
      outputMesh->minimumExtent[i] =
          (i < D) ? std::numeric_limits<T>::max() : 0.0;
      outputMesh->maximumExtent[i] =
          (i < D) ? std::numeric_limits<T>::lowest() : 0.0;
    }

    unsigned pointId = 0;
    for (const auto &result : results) {
      if (!result.valid)
        continue;

      // Create points for target and sample positions
      Vec3D<T> coordTarget{};
      Vec3D<T> coordSample{};

      for (int i = 0; i < D; ++i) {
        if (i == result.measureDimension) {
          coordTarget[i] = result.positionTarget;
          coordSample[i] = result.positionSample;
        } else {
          // Use midpoint of range for other dimensions
          T mid = (result.minBounds[i] + result.maxBounds[i]) / 2.0;
          // If bounds are infinite, use 0
          if (std::abs(mid) > std::numeric_limits<T>::max() / 2.0)
            mid = 0.0;
          coordTarget[i] = mid;
          coordSample[i] = mid;
        }
      }

      // Add target point
      nodeCoordinates.push_back(coordTarget);
      vertexIndices.push_back({pointId++});
      differenceValues.push_back(result.difference);
      targetValues.push_back(result.positionTarget);
      sampleValues.push_back(result.positionSample);

      // Add sample point
      nodeCoordinates.push_back(coordSample);
      vertexIndices.push_back({pointId++});
      differenceValues.push_back(result.difference);
      targetValues.push_back(result.positionTarget);
      sampleValues.push_back(result.positionSample);

      // Update extent
      for (unsigned i = 0; i < D; ++i) {
        outputMesh->minimumExtent[i] = std::min(
            {outputMesh->minimumExtent[i], coordTarget[i], coordSample[i]});
        outputMesh->maximumExtent[i] = std::max(
            {outputMesh->maximumExtent[i], coordTarget[i], coordSample[i]});
      }
    }

    if (!nodeCoordinates.empty()) {
      outputMesh->nodes = std::move(nodeCoordinates);
      outputMesh->vertices = std::move(vertexIndices);
      outputMesh->pointData.insertNextScalarData(std::move(differenceValues),
                                                 "Difference");
      outputMesh->pointData.insertNextScalarData(std::move(targetValues),
                                                 "TargetPosition");
      outputMesh->pointData.insertNextScalarData(std::move(sampleValues),
                                                 "SamplePosition");
    }
  }
};

// Add template specializations for this class
PRECOMPILE_PRECISION_DIMENSION(CompareCriticalDimensions)

} // namespace viennals
