#pragma once

#include <hrleSparseIterator.hpp>
#include <lsDomain.hpp>
#include <vcVectorType.hpp>

namespace viennals {
using namespace viennacore;

template <class T, int D> class CalculateVisibilities {
  using hrleDomainType = typename Domain<T, D>::DomainType;

  SmartPointer<Domain<T, D>> levelSet;
  Vec3D<T> direction;
  const T epsilon = static_cast<T>(1e-6);
  const std::string visibilitiesLabel;

public:
  CalculateVisibilities(const SmartPointer<Domain<T, D>> &passedLevelSet,
                        const Vec3D<T> &passedDirection,
                        std::string label = "Visibilities")
      : levelSet(passedLevelSet), direction(passedDirection),
        visibilitiesLabel(std::move(label)) {}

  void apply() {
    auto &domain = levelSet->getDomain();
    auto &grid = levelSet->getGrid();

    // *** Determine extents of domain ***
    Vec3D<T> minDefinedPoint;
    Vec3D<T> maxDefinedPoint;
    // Initialize with extreme values
    for (int i = 0; i < D; ++i) {
      minDefinedPoint[i] = std::numeric_limits<T>::max();
      maxDefinedPoint[i] = std::numeric_limits<T>::lowest();
    }
    // Iterate through all defined points in the domain
    for (viennahrle::SparseIterator<hrleDomainType> it(domain);
         !it.isFinished(); it.next()) {
      if (!it.isDefined())
        continue; // Skip undefined points

      // Get the coordinate of the current point
      auto point = it.getStartIndices();
      for (int i = 0; i < D; ++i) {
        // Compare to update min and max defined points
        T coord = point[i]; // * grid.getGridDelta();
        minDefinedPoint[i] = std::min(minDefinedPoint[i], coord);
        maxDefinedPoint[i] = std::max(maxDefinedPoint[i], coord);
      }
    }
    //****************************

    // Invert the vector
    auto dir = Normalize(Inv(direction));
    std::vector<T> visibilities(domain.getNumberOfPoints(), static_cast<T>(-1));

#pragma omp parallel num_threads(levelSet->getNumberOfSegments())
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      const viennahrle::Index<D> startVector =
          (p == 0) ? grid.getMinGridPoint() : domain.getSegmentation()[p - 1];

      const viennahrle::Index<D> endVector =
          (p != static_cast<int>(domain.getNumberOfSegments() - 1))
              ? domain.getSegmentation()[p]
              : grid.incrementIndices(grid.getMaxGridPoint());

      for (viennahrle::SparseIterator<hrleDomainType> it(domain, startVector);
           it.getStartIndices() < endVector; ++it) {

        if (!it.isDefined())
          continue;

        // Starting position of the point
        Vec3D<T> currentPos;
        for (int i = 0; i < D; ++i) {
          currentPos[i] = it.getStartIndices(i);
        }

        // Start tracing the ray
        T minLevelSetValue = it.getValue(); // Starting level set value
        Vec3D<T> rayPos = currentPos;

        // Step once to skip immediate neighbor
        for (int i = 0; i < D; ++i)
          rayPos[i] += dir[i];

        bool visibility = true;

        while (true) {
          // Update the ray position
          for (int i = 0; i < D; ++i)
            rayPos[i] += dir[i];

          // Determine the nearest grid cell (round to nearest index)
          viennahrle::Index<D> nearestCell;
          for (int i = 0; i < D; ++i)
            nearestCell[i] = static_cast<viennahrle::IndexType>(rayPos[i]);

          // Check if the nearest cell is within bounds
          bool outOfBounds = false;
          for (int i = 0; i < D; ++i) {
            if (nearestCell[i] < minDefinedPoint[i] ||
                nearestCell[i] > maxDefinedPoint[i]) {
              outOfBounds = true;
              break;
            }
          }

          if (outOfBounds)
            break; // Ray is outside the grid

          // Access the level set value at the nearest cell
          T value =
              viennahrle::SparseIterator<hrleDomainType>(domain, nearestCell)
                  .getValue();

          // Update the minimum value encountered
          if (value < minLevelSetValue) {
            visibility = false;
            break;
          }
        }

        // Update visibility for this point
        visibilities[it.getPointId()] = visibility ? 1.0 : 0.0;
      }
    }

    int unassignedCount = 0;
    for (size_t i = 0; i < visibilities.size(); ++i) {
      if (visibilities[i] < 0) {
        std::cerr << "Unassigned visibility at point ID: " << i << std::endl;
        ++unassignedCount;
      }
    }
    if (unassignedCount > 0) {
      std::cerr << "[Error] Total unassigned points: " << unassignedCount
                << std::endl;
    }

    auto &pointData = levelSet->getPointData();
    // delete if already exists
    if (int i = pointData.getScalarDataIndex(visibilitiesLabel); i != -1) {
      pointData.eraseScalarData(i);
    }
    pointData.insertNextScalarData(visibilities, visibilitiesLabel);
  }
};

} // namespace viennals