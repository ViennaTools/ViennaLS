#pragma once

#include <lsPreCompileMacros.hpp>

#include <algorithm>
#include <list>
#include <unordered_map>

#include <lsGeometries.hpp>
#include <lsMesh.hpp>

#include <vcLogger.hpp>
#include <vcSmartPointer.hpp>
#include <vcVectorType.hpp>

namespace viennals {

using namespace viennacore;

/// This algorithm creates a convex hull mesh from a
/// point cloud. This is done using the gift wrapping approach.
/// The points in the point cloud MUST be unique, otherwise this will fail.
template <class T, int D> class ConvexHull {
  typedef VectorType<unsigned, D - 1> EdgeType;

  SmartPointer<Mesh<T>> mesh = nullptr;
  SmartPointer<PointCloud<T, D>> pointCloud = nullptr;
  std::vector<EdgeType> visitedEdges;
  std::vector<VectorType<unsigned, D>> hullElements;
  std::list<EdgeType> remainingEdges;

  // go through all points in the set and find the correct element
  // to do this find the normal vector and dot product to see which side
  // a point is on
  unsigned pivotEdge(EdgeType currentEdge) const {
    auto &points = pointCloud->points;
    unsigned nextIndex = 0;

    // start from index which is neither one of the edges
    while (nextIndex == currentEdge[0] ||
           (D == 3 && nextIndex == currentEdge[1])) {
      ++nextIndex;
    }

    for (unsigned i = 0; i < points.size(); ++i) {
      if (i == currentEdge[0] || i == nextIndex)
        continue;
      if (D == 3 && i == currentEdge[1])
        continue;

      // surface element normal and distance to point i
      VectorType<T, D> normal;
      VectorType<T, D> distance = points[i] - points[currentEdge[0]];
      // create surface element and check if all points are to its right
      if constexpr (D == 2) {
        // in 2D normal is a 90 degree rotation
        normal = points[nextIndex] - points[currentEdge[0]];
        if (D == 2) {
          std::swap(normal[0], normal[1]);
          normal[0] = -normal[0];
        }
      } else if constexpr (D == 3) {
        auto v1 = points[currentEdge[1]] - points[currentEdge[0]];
        auto v2 = points[nextIndex] - points[currentEdge[0]];
        normal = CrossProduct(v1, v2);
        Normalize(normal);
      }

      auto product = DotProduct(distance, normal);
      // if dot product is very small, point is very close to plane
      // we need to check if we already have the correct point
      // or if it is the next correct one
      if (std::abs(product) < 1e-9) {
        // check if suggested triangle intersects with any other
        // in the same plane
        VectorType<unsigned, D> triangle;
        triangle[0] = currentEdge[0];
        triangle[1] = currentEdge[1];
        triangle[2] = i;
        if (doesTriangleClip(triangle)) {
          continue;
        }

        EdgeType edges[2];
        edges[0][0] = currentEdge[0];
        edges[0][1] = i;
        edges[1][0] = currentEdge[1];
        edges[1][1] = i;

        if (!wasEdgeVisited(edges[0]) && !wasEdgeVisited(edges[1])) {
          nextIndex = i;
        }
      }
      // check if point is to the right of current element
      else if (product > 0) {
        nextIndex = i;
      }
    }
    return nextIndex;
  }

  // check if edge was already visited
  bool wasEdgeVisited(const EdgeType &edge) const {
    for (unsigned i = 0; i < visitedEdges.size(); ++i) {
      // in 2D one match is enough, in 3D both points must be the same
      if (edge[0] == visitedEdges[i][0]) {
        if constexpr (D == 3) {
          if (edge[1] == visitedEdges[i][1]) {
            return true;
          }
        } else {
          return true;
        }
      }

      // in 3D, an edge is also the same if they contain the same nodes
      if constexpr (D == 3) {
        if (edge[0] == visitedEdges[i][1] && edge[1] == visitedEdges[i][0]) {
          return true;
        }
      }
    }
    return false;
  }

  // find whether an edge has already been used for a triangle once
  typename std::list<EdgeType>::iterator findEdge(EdgeType edge) {
    for (typename std::list<EdgeType>::iterator it = remainingEdges.begin();
         it != remainingEdges.end(); ++it) {
      bool containsElem[2] = {false, false};
      for (unsigned i = 0; i < D - 1; ++i) {
        // check all possible permutations of elements
        for (unsigned k = 0; k < D - 1; ++k) {
          if ((*it)[i] == edge[k]) {
            containsElem[i] = true;
          }
        }
      }
      if (containsElem[0] && containsElem[1]) {
        return it;
      }
    }
    return remainingEdges.end();
  }

  // return whether two triangles with shared nodes intersect
  bool intersectSharedNode(const VectorType<unsigned, D> &triangle1,
                           const VectorType<unsigned, D> &triangle2) const {
    // find which nodes are shared and which nodes are others
    std::vector<unsigned> sharedNodes;
    // otherNodes2 contains the nodes of triangle1 which are not shared
    std::vector<unsigned> otherNodes;

    for (unsigned i = 0; i < D; ++i) {
      unsigned shared = sharedNodes.size();
      for (unsigned j = 0; j < D; ++j) {
        if (triangle1[i] == triangle2[j]) {
          // node numbers of t2
          sharedNodes.push_back(j);
        }
      }
      // if triangle1[i] was not shared, push to other nodes
      if (shared == sharedNodes.size()) {
        otherNodes.push_back(i);
      }
    }

    assert(!sharedNodes.empty());

    auto &points = pointCloud->points;

    // for each shared node, calculate the edges of t2 leading away from it.
    // then generate the average vector and its angle to the edges.
    // then calculate the dot product of each other node of t1 to the average
    // vector. If it is smaller than the dot product between the average vector
    // and an edge, the node must be outside of triangle. If it is inside, one
    // triangle must clip the other.
    for (unsigned int &sharedNode : sharedNodes) {
      VectorType<T, D> averageVector;
      // save the dot product between average and one edge
      double centerEdgeDot;
      {
        auto shared = sharedNode;
        VectorType<T, D> edge1 = Normalize(points[triangle2[(shared + 1) % D]] -
                                           points[triangle2[shared]]);
        VectorType<T, D> edge2 = Normalize(points[triangle2[(shared + 2) % D]] -
                                           points[triangle2[shared]]);
        averageVector = edge1 + edge2;
        averageVector /= 2;
        Normalize(averageVector);
        centerEdgeDot = DotProduct(averageVector, edge2);
      }

      // go over other nodes of triangle1
      for (unsigned int &otherNode : otherNodes) {
        VectorType<T, D> vector = Normalize(points[triangle1[otherNode]] -
                                            points[triangle2[sharedNode]]);

        if (DotProduct(vector, averageVector) > centerEdgeDot + 1e-9) {
          return true;
        }
      }
    }
    // if no other node was ever within two edges, the triangles do not clip
    return false;
  }

  // check if triangle defined by two edges clips any other triangle
  bool doesTriangleClip(const VectorType<unsigned, D> &triangle) const {
    auto &points = pointCloud->points;

    auto triangleNormal =
        calculateNormal(points[triangle[1]] - points[triangle[0]],
                        points[triangle[2]] - points[triangle[0]]);

    // unsigned shareEdges = 0;
    for (unsigned i = 0; i < hullElements.size(); ++i) {
      unsigned inOneTriangle = 0;
      for (unsigned j = 0; j < D; ++j) {
        // check all possible permutations of elements
        for (unsigned k = 0; k < D; ++k) {
          if (hullElements[i][j] == triangle[k]) {
            ++inOneTriangle;
          }
        }
      }

      // if they share at least one node, they might clip, so check
      if (inOneTriangle > 0) {
        // check if they are in the same plane
        auto normal2 = calculateNormal(
            points[hullElements[i][1]] - points[hullElements[i][0]],
            points[hullElements[i][2]] - points[hullElements[i][0]]);

        bool skip = false;
        for (unsigned d = 0; d < D; ++d) {
          double diff =
              std::abs(std::abs(triangleNormal[d]) - std::abs(normal2[d]));
          if (diff > 1e-6)
            skip = true;
        }

        if (skip)
          continue;

        if (intersectSharedNode(triangle, hullElements[i]))
          return true;
        if (intersectSharedNode(hullElements[i], triangle))
          return true;
      }
    }
    return false;
  }

  // calculate the normal vector of two vectors
  VectorType<T, D> calculateNormal(const VectorType<T, D> &v1,
                                   const VectorType<T, D> &v2) const {
    VectorType<T, D> newNormal;
    newNormal[0] = v1[1] * v2[2] - v1[2] * v2[1];
    newNormal[1] = v1[2] * v2[0] - v1[0] * v2[2];
    newNormal[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return newNormal;
  }

  // calculate Area of triangle defined by 2 vectors
  T calculateArea(const VectorType<T, D> &v1,
                  const VectorType<T, D> &v2) const {
    VectorType<T, D> newNormal = CrossProduct(v1, v2);
    return std::sqrt(DotProduct(newNormal, newNormal));
  }

public:
  ConvexHull() = default;

  ConvexHull(SmartPointer<Mesh<T>> passedMesh,
             SmartPointer<PointCloud<T, D>> passedPointCloud)
      : mesh(passedMesh), pointCloud(passedPointCloud) {}

  void setMesh(SmartPointer<Mesh<T>> passedMesh) { mesh = passedMesh; }

  void setPointCloud(SmartPointer<PointCloud<T, D>> passedPointCloud) {
    pointCloud = passedPointCloud;
  }

  void apply() {
    if (mesh == nullptr) {
      Logger::getInstance()
          .addWarning("No mesh was passed to ConvexHull.")
          .print();
      return;
    }
    if (pointCloud == nullptr) {
      Logger::getInstance()
          .addWarning("No point cloud was passed to ConvexHull.")
          .print();
      return;
    }

    mesh->clear();
    auto &points = pointCloud->points;
    if (points.size() == 0)
      return;

    {
      // find first hull point
      auto it = std::min_element(points.begin(), points.end());
      unsigned currentIndex = std::distance(points.begin(), it);

      if (D == 2) {
        remainingEdges.push_back(EdgeType(currentIndex));
      } else {
        // need to check if there is second point at same z coord
        int currentPointIndex = -1;
        double hullMetric = std::numeric_limits<double>::max();
        for (int i = 0; i < points.size(); ++i) {
          if (i == currentIndex)
            continue;
          if (std::abs(points[i][2] - points[currentIndex][2]) < 1e-7) {
            auto diff = points[currentIndex] - points[i];
            // choose closest point if points are in z plane
            if (const double currentHullMetric = DotProduct(diff, diff);
                currentHullMetric < hullMetric) {
              hullMetric = currentHullMetric;
              currentPointIndex = i;
            }
          }
        }

        EdgeType edge;
        edge[0] = currentIndex;
        // if there was no point at same z, find through pivot
        if (currentPointIndex == -1) {
          VectorType<T, D> newPoint;
          newPoint = points[currentIndex];
          newPoint[0] += 1.0;
          // take newPoint for fake edge to find correct first point
          points.push_back(newPoint);
          edge[1] = points.size() - 1;
          edge[1] = pivotEdge(edge);
          points.pop_back();
        } else {
          edge[1] = currentPointIndex;
        }

        remainingEdges.push_back(edge);
      }
    }

    // until there are no more edges to check
    while (!remainingEdges.empty()) {
      // for all points in the cloud
      EdgeType currentEdge = remainingEdges.back();
      remainingEdges.pop_back();

      if (wasEdgeVisited(currentEdge))
        continue;

      // find next point in hull
      unsigned nextIndex = pivotEdge(currentEdge);

      // set new hull element and save
      VectorType<unsigned, D> element;
      for (unsigned i = 0; i < D - 1; ++i) {
        element[i] = currentEdge[i];
      }
      element[D - 1] = nextIndex;
      hullElements.push_back(element);

      // mark edge as visited
      visitedEdges.push_back(currentEdge);
      if (D == 2) {
        remainingEdges.push_back(EdgeType(nextIndex));
      } else {
        EdgeType edge;
        edge[0] = nextIndex;
        edge[1] = currentEdge[1];
        // if edge is already part of another element
        // and we just added it, we need to remove it
        // from the list
        auto it = findEdge(edge);
        if (it != remainingEdges.end()) {
          remainingEdges.erase(it);
        } else {
          remainingEdges.push_back(edge);
        }

        edge[0] = currentEdge[0];
        edge[1] = nextIndex;
        it = findEdge(edge);
        if (it != remainingEdges.end()) {
          remainingEdges.erase(it);
        } else {
          remainingEdges.push_back(edge);
        }
      }
    }

    // now make the mesh
    auto &elements = mesh->template getElements<D>();
    std::unordered_map<unsigned, unsigned> oldToNewNodes;
    for (unsigned i = 0; i < hullElements.size(); ++i) {
      // here add translation of old index to new index for new points
      std::array<unsigned, D> newElement;
      for (unsigned j = 0; j < D; ++j) {
        auto insertion = oldToNewNodes.insert(
            std::make_pair(hullElements[i][j], oldToNewNodes.size()));

        // if insertion took place, add the point to the nodes
        if (insertion.second) {
          Vec3D<T> node{points[hullElements[i][j]][0],
                        points[hullElements[i][j]][1],
                        (D == 2) ? T(0.) : points[hullElements[i][j]][2]};
          mesh->nodes.push_back(node);
        }

        newElement[j] = insertion.first->second;
      }

      elements.push_back(newElement);
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(ConvexHull)

} // namespace viennals
