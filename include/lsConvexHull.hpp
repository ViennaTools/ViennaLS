#ifndef LS_CONVEX_HULL_HPP
#define LS_CONVEX_HULL_HPP

#include <algorithm>
#include <lsPreCompileMacros.hpp>

#include <hrleVectorType.hpp>

#include <lsGeometries.hpp>
#include <lsMesh.hpp>
#include <lsMessage.hpp>

/// This algorithm creates a convex hull mesh from a
/// point cloud. This is done using the gift wrapping approach.
template <class T, int D> class lsConvexHull {
  lsMesh *mesh = nullptr;
  lsPointCloud<T, D> *pointCloud = nullptr;

public:
  lsConvexHull(lsMesh &passedMesh, lsPointCloud<T, D> &passedPointCloud)
      : mesh(&passedMesh), pointCloud(&passedPointCloud) {}

  void setMesh(lsMesh &passedMesh) { mesh = &passedMesh; }

  void setPointCloud(lsPointCloud<T, D> &passedPointCloud) {
    pointCloud = &passedPointCloud;
  }

  void apply() {
    if (mesh == nullptr) {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsConvexHull.")
          .print();
      return;
    }
    if (pointCloud == nullptr) {
      lsMessage::getInstance()
          .addWarning("No point cloud was passed to lsConvexHull.")
          .print();
      return;
    }

    mesh->clear();
    auto &points = pointCloud->points;

    // find the outmost point in highest cartesian direction
    auto it = std::min_element(points.begin(), points.end());
    unsigned currentIndex = std::distance(points.begin(), it);

    std::vector<unsigned> hullIndices(1, currentIndex);

    // until we reach start again
    while (true) {
      // for all points in the cloud
      unsigned nextIndex = (currentIndex + 1) % (points.size());
      for (unsigned i = 0; i < points.size(); ++i) {
        if (i == currentIndex)
          continue;
        if (i == nextIndex)
          continue;

        // create surface element and check if all points are to its right
        auto distance = points[i] - points[currentIndex];
        // in 2D normal is a 90 degree rotation
        auto normal = points[nextIndex] - points[currentIndex];
        if (D == 2) {
          std::swap(normal[0], normal[1]);
          normal[0] = -normal[0];
        }

        // check if point is to the left of current line
        if (DotProduct(distance, normal) < 0) {
          nextIndex = i;
        }
      }

      if (nextIndex == hullIndices.front()) {
        break;
      } else {
        hullIndices.push_back(nextIndex);
        currentIndex = nextIndex;
      }
    }

    // now make the 2D lines
    auto &elements = mesh->getElements<D>();
    for (unsigned i = 0; i < hullIndices.size(); ++i) {
      hrleVectorType<double, 3> node(points[hullIndices[i]][0],
                                     points[hullIndices[i]][1],
                                     (D == 2) ? 0. : points[hullIndices[i]][2]);
      mesh->nodes.push_back(node);

      elements.push_back(
          hrleVectorType<unsigned, 2>(i, (i + 1) % hullIndices.size()));
    }
  }
};

#endif // LS_CONVEX_HULL_HPP
