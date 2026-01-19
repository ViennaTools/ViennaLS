#pragma once

#include <lsDomain.hpp>
#include <lsPreCompileMacros.hpp>
#include <lsExpand.hpp>
#include <hrleSparseStarIterator.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace viennals {

using namespace viennacore;

/// Resamples a level set to a new grid resolution (gridDelta).
/// This is essential for Multi-Resolution Advection strategies.
template <class T, int D> class Resample {
  SmartPointer<Domain<T, D>> inputLevelSet = nullptr;
  SmartPointer<Domain<T, D>> outputLevelSet = nullptr;
  double targetGridDelta = 1.0;
  bool refine = false;
  bool subGridCorrection = true;

public:
  Resample() = default;

  Resample(SmartPointer<Domain<T, D>> input,
           SmartPointer<Domain<T, D>> output,
           double newGridDelta)
      : inputLevelSet(input), outputLevelSet(output), targetGridDelta(newGridDelta) {}

  void setInputLevelSet(SmartPointer<Domain<T, D>> input) {
    inputLevelSet = input;
  }

  void setOutputLevelSet(SmartPointer<Domain<T, D>> output) {
    outputLevelSet = output;
  }

  void setTargetGridDelta(double newGridDelta) {
    targetGridDelta = newGridDelta;
  }

  void setSubGridCorrection(bool active) {
    subGridCorrection = active;
  }

  void apply() {
    if (inputLevelSet == nullptr || outputLevelSet == nullptr) {
      Logger::getInstance()
          .addError("Resample: Input or Output LevelSet not set.")
          .print();
      return;
    }

    auto &inGrid = inputLevelSet->getGrid();
    double sourceDelta = inGrid.getGridDelta();
    double ratio = targetGridDelta / sourceDelta;

    // Ensure input has enough width to cover the coarse grid points
    // We need the fine LS to cover at least ~3 layers of the coarse grid
    int requiredWidth = std::ceil(3.0 * ratio) + 2;
    Expand<T, D>(inputLevelSet, requiredWidth).apply();

    // Initialize output with new grid delta but same physical boundaries
    double bounds[2 * D];
    auto bcs = inGrid.getBoundaryConditions();

    // We expanded the level set, so we should also expand the bounds
    // to ensure the new coarse grid covers the expanded data.
    double padding = requiredWidth * sourceDelta * 2.;
    
    for(int i=0; i<D; ++i) {
        double minCoord = inGrid.getMinIndex(i) * sourceDelta;
        double maxCoord = inGrid.getMaxIndex(i) * sourceDelta;
        if (bcs[i] == viennahrle::BoundaryType::INFINITE_BOUNDARY) {
            minCoord -= padding;
            maxCoord += padding;
        }
        bounds[2 * i] = minCoord;
        bounds[2 * i + 1] = maxCoord;
    }

    auto tempLS = SmartPointer<Domain<T, D>>::New(bounds, bcs.data(), targetGridDelta);

    // Strategy:
    // 1. Iterate over the INPUT (Sparse).
    // 2. Map input coordinate to OUTPUT index.
    // 3. Insert value.
    // Note: This is a "Point Injection" strategy. For high-quality upsampling,
    // a Neural Network should be applied AFTER this step on the output domain.
    
    std::vector<std::pair<viennahrle::Index<D>, T>> newPoints;
    
    // Reserve estimation
    double invRatio = sourceDelta / targetGridDelta;
    double volRatio = std::pow(invRatio, D);
    // newPoints.reserve(inputLevelSet->getNumberOfPoints() * volRatio);

    viennahrle::ConstSparseStarIterator<typename Domain<T, D>::DomainType, 1> it(inputLevelSet->getDomain());
    for (; !it.isFinished(); it.next()) {
      if (!it.getCenter().isDefined()) continue;

      // Get Physical Coordinate
      // Map to New Index
      viennahrle::Index<D> newIdx;
      auto indices = it.getIndices();
      for(int i=0; i<D; ++i) {
          double coord = indices[i] * sourceDelta;
          newIdx[i] = std::round(coord / targetGridDelta);
      }

      // Pre-calculate gradients for sub-grid correction
      T gradients[D];
      bool hasGradient[D];
      if (subGridCorrection) {
          for(int i=0; i<D; ++i) {
              if (it.getNeighbor(i).isDefined() && it.getNeighbor(i+D).isDefined()) {
                  T posVal = it.getNeighbor(i).getValue();
                  T negVal = it.getNeighbor(i+D).getValue();
                  gradients[i] = (posVal - negVal) * 0.5;
                  hasGradient[i] = true;
              } else {
                  gradients[i] = 0;
                  hasGradient[i] = false;
              }
          }
      }

      // Determine if we are upsampling (target finer than source)
      // If so, we need to fill multiple fine points per coarse point
      int splatRadius = 0;
      if (targetGridDelta < sourceDelta) {
          splatRadius = std::ceil((sourceDelta / targetGridDelta) * 0.5 + 1e-6);
      }

      // Iterate over the splat footprint (or just 0 for downsampling)
      // Using a flat loop to handle arbitrary dimension D
      long long numSplat = std::pow(2*splatRadius + 1, D);
      for(long long k=0; k<numSplat; ++k) {
          viennahrle::Index<D> offset;
          long long temp = k;
          bool inBounds = true;
          
          for(int d=0; d<D; ++d) {
              int dimLen = 2*splatRadius + 1;
              int val = (temp % dimLen) - splatRadius;
              temp /= dimLen;
              offset[d] = val;
              
              // Check if this fine point is within the coarse cell's influence
              if (std::abs(offset[d]) * targetGridDelta > sourceDelta * 0.5 + 1e-6) {
                  inBounds = false;
              }
          }
          
          if (!inBounds) continue;

          viennahrle::Index<D> fineIdx = newIdx + offset;
          T correction = 0;
          for(int d=0; d<D; ++d) {
              if (subGridCorrection && hasGradient[d]) {
                  double dist = fineIdx[d] * targetGridDelta - indices[d] * sourceDelta;
                  correction += gradients[d] * (dist / sourceDelta);
              }
          }
          T val = (it.getCenter().getValue() + correction) * (sourceDelta / targetGridDelta);
          newPoints.push_back(std::make_pair(fineIdx, val));
      }
    }

    // Sort by index to group duplicates
    std::sort(newPoints.begin(), newPoints.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });

    // Filter duplicates: keep the one with smallest absolute value (closest to
    // interface)
    std::vector<std::pair<viennahrle::Index<D>, T>> uniquePoints;
    if (!newPoints.empty()) {
      auto current = newPoints[0];
      for (size_t i = 1; i < newPoints.size(); ++i) {
        if (newPoints[i].first == current.first) {
          // Duplicate index, keep the one closer to 0
          if (std::abs(newPoints[i].second) < std::abs(current.second)) {
            current = newPoints[i];
          }
        } else {
          uniquePoints.push_back(current);
          current = newPoints[i];
        }
      }
      uniquePoints.push_back(current);
    }

    tempLS->insertPoints(uniquePoints);
    tempLS->finalize(inputLevelSet->getLevelSetWidth());
    
    outputLevelSet->deepCopy(tempLS);
  }
};

PRECOMPILE_PRECISION_DIMENSION(Resample)

} // namespace viennals