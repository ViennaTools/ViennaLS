#pragma once

#include <lsDomain.hpp>
#include <vcVectorUtil.hpp>

namespace viennals {
using namespace viennacore;

template <class NumericType, int D> class CalculateVisibilities {
  SmartPointer<Domain<NumericType, D>> levelSet;
  Vec3D<NumericType> direction;
  const NumericType epsilon = static_cast<NumericType>(1e-6);

public:
  static constexpr char visibilitiesLabel[] = "Visibilities";

  CalculateVisibilities(
      const SmartPointer<Domain<NumericType, D>> &passedLevelSet,
      const Vec3D<NumericType> passedDirection)
      : levelSet(passedLevelSet), direction(passedDirection) {}

  void apply() {

    auto &domain = levelSet->getDomain();
    auto &grid = levelSet->getGrid();

    // *** Determine extents of domain ***
    Vec3D<NumericType> minDefinedPoint;
    Vec3D<NumericType> maxDefinedPoint;
    // Initialize with extreme values
    for (int i = 0; i < D; ++i) {
        minDefinedPoint[i] = std::numeric_limits<NumericType>::max();
        maxDefinedPoint[i] = std::numeric_limits<NumericType>::lowest();
    }
    // Iterate through all defined points in the domain
    for (hrleSparseIterator<typename Domain<NumericType, D>::DomainType> it(domain); !it.isFinished(); it.next()) {
      if (!it.isDefined()) continue; // Skip undefined points

      // Get the coordinate of the current point
      auto point = it.getStartIndices();
      for (int i = 0; i < D; ++i) {
        // Compare to update min and max defined points
        NumericType coord = point[i]; // * grid.getGridDelta();
        minDefinedPoint[i] = std::min(minDefinedPoint[i], coord);
        maxDefinedPoint[i] = std::max(maxDefinedPoint[i], coord);
      }
    }
    //****************************

    // Invert the vector
    auto dir = Normalize(Inv(direction));

    auto numDefinedPoints = domain.getNumberOfPoints();
    std::vector<NumericType> visibilities(numDefinedPoints);

      // std::vector<Vec3D<hrleIndexType>> visitedCells;
    
      hrleSizeType id = 0;
      hrleSparseIterator<typename Domain<NumericType, D>::DomainType> it(domain);

      while (!it.isFinished()) {
          if (!it.isDefined()) {
              it.next();
              continue;
          }

          // Starting position of the point
          Vec3D<NumericType> currentPos;
          for (int i = 0; i < D; ++i) {
              currentPos[i] = it.getStartIndices(i);
          }

          // Start tracing the ray
          NumericType minLevelSetValue = it.getValue(); // Starting level set value
          Vec3D<NumericType> rayPos = currentPos;
          bool visibility = true;

          while(1) {
            // Update the ray position
            for (int i = 0; i < D; ++i) {
              rayPos[i] += dir[i];
            }

            // Determine the nearest grid cell (round to nearest index)
            Vec3D<hrleIndexType> nearestCell;
            for (int i = 0; i < D; ++i) {
              nearestCell[i] = static_cast<hrleIndexType>(rayPos[i]);
            }

            // // Before adding a cell, check if it's already visited
            // if (std::find(visitedCells.begin(), visitedCells.end(), nearestCell) == visitedCells.end()) {
            //     visitedCells.push_back(nearestCell);
            // }

            // Check if the nearest cell is within bounds
            bool outOfBounds = false;
            for (int i = 0; i < D; ++i) {
              if (nearestCell[i] < minDefinedPoint[i] || nearestCell[i] > maxDefinedPoint[i]) {
                outOfBounds = true;
                break;
              }
            }

            if (outOfBounds) {
              break; // Ray is outside the grid
            }

            // Access the level set value at the nearest cell
            NumericType neighborValue = std::numeric_limits<NumericType>::max();
            hrleSparseIterator<typename Domain<NumericType, D>::DomainType> neighborIt(domain);
            neighborIt.goToIndices(nearestCell);
            if (neighborIt.isDefined()) {
              neighborValue = neighborIt.getValue();
            }

            // Update the minimum value encountered
            if (neighborValue < minLevelSetValue) {
              visibility = false;
              break;
            }
        }

        // Update visibility for this point
        visibilities[id++] = visibility ? 1.0 : 0.0;
        it.next();
      }

    auto &pointData = levelSet->getPointData();
    // delete if already exists
    if (int i = pointData.getScalarDataIndex("Visibilities"); i != -1) {
      pointData.eraseScalarData(i);
    }
    pointData.insertNextScalarData(visibilities, "Visibilities");

    assert(id == numDefinedPoints);
  }
};

} // namespace viennals