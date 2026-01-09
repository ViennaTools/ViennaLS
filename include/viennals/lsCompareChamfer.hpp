#pragma once

#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMesh.hpp>
#include <lsPreCompileMacros.hpp>
#include <lsToSurfaceMesh.hpp>

#include <vcKDTree.hpp>

#include <cmath>
#include <limits>

namespace viennals {

using namespace viennacore;

/// Calculate Chamfer distance between two level sets by comparing their
/// zero-level-set surfaces. The Chamfer distance is a bidirectional metric
/// that measures the average nearest-neighbor distance between two point sets.
///
/// This class extracts the surface representations of both level sets and
/// computes:
/// - Forward distance: average distance from target surface to sample surface
/// - Backward distance: average distance from sample surface to target surface
/// - Chamfer distance: average of forward and backward distances
/// - RMS Chamfer distance: root mean square of nearest-neighbor distances
/// - Maximum distance: maximum nearest-neighbor distance across both directions
///
/// The code works for 2D and 3D level sets. Surfaces are represented as line
/// segments in 2D and triangles in 3D.
///
/// Both level sets must have a width of at least 2 to extract surfaces. If
/// not, they will be automatically expanded.
///
/// Note for the future: lsToDiskMesh could be used instead of lsToSurfaceMesh,
/// which is probably more efficient but slightly less accurate.

template <class T, int D = 2> class CompareChamfer {
  SmartPointer<Domain<T, D>> levelSetTarget = nullptr;
  SmartPointer<Domain<T, D>> levelSetSample = nullptr;

  // Results
  T forwardDistance = 0.0;    // Target → Sample
  T backwardDistance = 0.0;   // Sample → Target
  T chamferDistance = 0.0;    // Average of forward and backward
  T rmsChamferDistance = 0.0; // RMS of all nearest-neighbor distances
  T maxDistance = 0.0;        // Maximum nearest-neighbor distance
  unsigned numTargetPoints = 0;
  unsigned numSamplePoints = 0;

  // Optional mesh output
  SmartPointer<Mesh<T>> outputMeshTarget = nullptr;
  SmartPointer<Mesh<T>> outputMeshSample = nullptr;

  bool checkCompatibility() {
    if (levelSetTarget == nullptr || levelSetSample == nullptr) {
      VIENNACORE_LOG_WARNING("Missing level set in CompareChamfer.");
      return false;
    }

    // Check if the grids are compatible
    const auto &gridTarget = levelSetTarget->getGrid();
    const auto &gridSample = levelSetSample->getGrid();

    if (gridTarget.getGridDelta() != gridSample.getGridDelta()) {
      VIENNACORE_LOG_WARNING("Grid delta mismatch in CompareChamfer. The grid "
                             "deltas of the two level sets must be equal.");
      return false;
    }

    return true;
  }

public:
  CompareChamfer() {}

  CompareChamfer(SmartPointer<Domain<T, D>> passedLevelSetTarget,
                 SmartPointer<Domain<T, D>> passedLevelSetSample)
      : levelSetTarget(passedLevelSetTarget),
        levelSetSample(passedLevelSetSample) {}

  /// Set the target level set
  void setLevelSetTarget(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetTarget = passedLevelSet;
  }

  /// Set the sample level set
  void setLevelSetSample(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetSample = passedLevelSet;
  }

  /// Set output mesh for target surface points with distance data
  void setOutputMeshTarget(SmartPointer<Mesh<T>> passedMesh) {
    outputMeshTarget = passedMesh;
  }

  /// Set output mesh for sample surface points with distance data
  void setOutputMeshSample(SmartPointer<Mesh<T>> passedMesh) {
    outputMeshSample = passedMesh;
  }

  /// Apply the Chamfer distance calculation
  void apply() {
    // Check compatibility
    if (!checkCompatibility()) {
      // Set NaN results to indicate failure
      forwardDistance = std::numeric_limits<T>::quiet_NaN();
      backwardDistance = std::numeric_limits<T>::quiet_NaN();
      chamferDistance = std::numeric_limits<T>::quiet_NaN();
      rmsChamferDistance = std::numeric_limits<T>::quiet_NaN();
      maxDistance = std::numeric_limits<T>::quiet_NaN();
      numTargetPoints = 0;
      numSamplePoints = 0;
      return;
    }

    // Ensure both level sets have sufficient width for surface extraction
    constexpr int minimumWidth = 2;

    auto workingTarget = levelSetTarget;
    auto workingSample = levelSetSample;

    if (levelSetTarget->getLevelSetWidth() < minimumWidth) {
      workingTarget = SmartPointer<Domain<T, D>>::New(levelSetTarget);
      Expand<T, D>(workingTarget, minimumWidth).apply();
      VIENNACORE_LOG_INFO(
          "CompareChamfer: Expanded target level set to width " +
          std::to_string(minimumWidth) + " for surface extraction.");
    }

    if (levelSetSample->getLevelSetWidth() < minimumWidth) {
      workingSample = SmartPointer<Domain<T, D>>::New(levelSetSample);
      Expand<T, D>(workingSample, minimumWidth).apply();
      VIENNACORE_LOG_INFO(
          "CompareChamfer: Expanded sample level set to width " +
          std::to_string(minimumWidth) + " for surface extraction.");
    }

    // Extract surface meshes
    auto targetSurfaceMesh = SmartPointer<Mesh<T>>::New();
    auto sampleSurfaceMesh = SmartPointer<Mesh<T>>::New();

    ToSurfaceMesh<T, D>(workingTarget, targetSurfaceMesh).apply();
    ToSurfaceMesh<T, D>(workingSample, sampleSurfaceMesh).apply();

    // Get surface points
    const auto &targetNodes = targetSurfaceMesh->getNodes();
    const auto &sampleNodes = sampleSurfaceMesh->getNodes();

    numTargetPoints = targetNodes.size();
    numSamplePoints = sampleNodes.size();

    if (numTargetPoints == 0 || numSamplePoints == 0) {
      VIENNACORE_LOG_WARNING(
          "CompareChamfer: One or both surfaces have no points. "
          "Cannot compute Chamfer distance.");
      forwardDistance = std::numeric_limits<T>::infinity();
      backwardDistance = std::numeric_limits<T>::infinity();
      chamferDistance = std::numeric_limits<T>::infinity();
      rmsChamferDistance = std::numeric_limits<T>::infinity();
      maxDistance = std::numeric_limits<T>::infinity();
      return;
    }

    // Convert nodes to format suitable for KDTree
    std::vector<std::array<T, D>> targetPoints(numTargetPoints);
    std::vector<std::array<T, D>> samplePoints(numSamplePoints);

    for (unsigned i = 0; i < numTargetPoints; ++i) {
      for (unsigned d = 0; d < D; ++d) {
        targetPoints[i][d] = targetNodes[i][d];
      }
    }

    for (unsigned i = 0; i < numSamplePoints; ++i) {
      for (unsigned d = 0; d < D; ++d) {
        samplePoints[i][d] = sampleNodes[i][d];
      }
    }

    // Build KD-trees
    KDTree<T, std::array<T, D>> targetTree(targetPoints);
    targetTree.build();

    KDTree<T, std::array<T, D>> sampleTree(samplePoints);
    sampleTree.build();

    // Compute forward distance (target → sample)
    T sumForward = 0.0;
    T sumSquaredForward = 0.0;
    T maxForward = 0.0;
    std::vector<T> targetDistances;
    if (outputMeshTarget != nullptr) {
      targetDistances.reserve(numTargetPoints);
    }

    for (unsigned i = 0; i < numTargetPoints; ++i) {
      auto nearest = sampleTree.findNearest(targetPoints[i]);
      if (nearest.has_value()) {
        T dist = nearest.value().second;
        sumForward += dist;
        sumSquaredForward += dist * dist;
        maxForward = std::max(maxForward, dist);
        if (outputMeshTarget != nullptr) {
          targetDistances.push_back(dist);
        }
      }
    }

    forwardDistance = sumForward / numTargetPoints;

    // Compute backward distance (sample → target)
    T sumBackward = 0.0;
    T sumSquaredBackward = 0.0;
    T maxBackward = 0.0;
    std::vector<T> sampleDistances;
    if (outputMeshSample != nullptr) {
      sampleDistances.reserve(numSamplePoints);
    }

    for (unsigned i = 0; i < numSamplePoints; ++i) {
      auto nearest = targetTree.findNearest(samplePoints[i]);
      if (nearest.has_value()) {
        T dist = nearest.value().second;
        sumBackward += dist;
        sumSquaredBackward += dist * dist;
        maxBackward = std::max(maxBackward, dist);
        if (outputMeshSample != nullptr) {
          sampleDistances.push_back(dist);
        }
      }
    }

    backwardDistance = sumBackward / numSamplePoints;

    // Compute overall metrics
    chamferDistance =
        (sumForward + sumBackward) / (numTargetPoints + numSamplePoints);
    rmsChamferDistance = std::sqrt((sumSquaredForward + sumSquaredBackward) /
                                   (numTargetPoints + numSamplePoints));
    maxDistance = std::max(maxForward, maxBackward);

    // Generate output meshes if requested
    if (outputMeshTarget != nullptr && !targetDistances.empty()) {
      outputMeshTarget->clear();
      outputMeshTarget->nodes = targetNodes;

      // Copy topology for visualization
      if constexpr (D == 2) {
        outputMeshTarget->lines = targetSurfaceMesh->lines;
      } else {
        outputMeshTarget->triangles = targetSurfaceMesh->triangles;
      }

      // Add distance data
      outputMeshTarget->pointData.insertNextScalarData(
          std::move(targetDistances), "DistanceToSample");

      // Set mesh extent
      for (unsigned d = 0; d < 3; ++d) {
        outputMeshTarget->minimumExtent[d] =
            (d < D) ? std::numeric_limits<T>::max() : 0.0;
        outputMeshTarget->maximumExtent[d] =
            (d < D) ? std::numeric_limits<T>::lowest() : 0.0;
      }
      for (unsigned d = 0; d < D; ++d) {
        for (const auto &node : targetNodes) {
          outputMeshTarget->minimumExtent[d] =
              std::min(outputMeshTarget->minimumExtent[d], node[d]);
          outputMeshTarget->maximumExtent[d] =
              std::max(outputMeshTarget->maximumExtent[d], node[d]);
        }
      }
    }

    if (outputMeshSample != nullptr && !sampleDistances.empty()) {
      outputMeshSample->clear();
      outputMeshSample->nodes = sampleNodes;

      // Copy topology for visualization
      if constexpr (D == 2) {
        outputMeshSample->lines = sampleSurfaceMesh->lines;
      } else {
        outputMeshSample->triangles = sampleSurfaceMesh->triangles;
      }

      // Add distance data
      outputMeshSample->pointData.insertNextScalarData(
          std::move(sampleDistances), "DistanceToTarget");

      // Set mesh extent
      for (unsigned d = 0; d < 3; ++d) {
        outputMeshSample->minimumExtent[d] =
            (d < D) ? std::numeric_limits<T>::max() : 0.0;
        outputMeshSample->maximumExtent[d] =
            (d < D) ? std::numeric_limits<T>::lowest() : 0.0;
      }
      for (unsigned d = 0; d < D; ++d) {
        for (const auto &node : sampleNodes) {
          outputMeshSample->minimumExtent[d] =
              std::min(outputMeshSample->minimumExtent[d], node[d]);
          outputMeshSample->maximumExtent[d] =
              std::max(outputMeshSample->maximumExtent[d], node[d]);
        }
      }
    }
  }

  /// Get the forward distance (average distance from target to sample)
  T getForwardDistance() const { return forwardDistance; }

  /// Get the backward distance (average distance from sample to target)
  T getBackwardDistance() const { return backwardDistance; }

  /// Get the Chamfer distance (average of forward and backward)
  T getChamferDistance() const { return chamferDistance; }

  /// Get the RMS Chamfer distance
  T getRMSChamferDistance() const { return rmsChamferDistance; }

  /// Get the maximum nearest-neighbor distance
  T getMaxDistance() const { return maxDistance; }

  /// Get the number of target surface points
  unsigned getNumTargetPoints() const { return numTargetPoints; }

  /// Get the number of sample surface points
  unsigned getNumSamplePoints() const { return numSamplePoints; }
};

// Add template specializations for this class
PRECOMPILE_PRECISION_DIMENSION(CompareChamfer)

} // namespace viennals
