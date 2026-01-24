#pragma once

#include <lsPreCompileMacros.hpp>

#include <iostream>
#include <map>

#include <hrleSparseCellIterator.hpp>
#include <hrleSparseStarIterator.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFiniteDifferences.hpp>
#include <lsMarchingCubes.hpp>
#include <lsMesh.hpp>

namespace viennals {

using namespace viennacore;

/// Extract an explicit Mesh<> instance from an lsDomain.
/// The interface is then described by explicit surface elements:
/// Lines in 2D, Triangles in 3D.
template <class T, int D> class ToSurfaceMesh {
  typedef typename Domain<T, D>::DomainType hrleDomainType;

  SmartPointer<Domain<T, D>> levelSet = nullptr;
  SmartPointer<Mesh<T>> mesh = nullptr;
  // std::vector<hrleIndexType> meshNodeToPointIdMapping;
  const T epsilon;
  bool updatePointData = true;
  bool checkCornersAndEdges = true;

public:
  explicit ToSurfaceMesh(double eps = 1e-12) : epsilon(eps) {}

  ToSurfaceMesh(const SmartPointer<Domain<T, D>> passedLevelSet,
                SmartPointer<Mesh<T>> passedMesh, double eps = 1e-12)
      : levelSet(passedLevelSet), mesh(passedMesh), epsilon(eps) {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  void setMesh(SmartPointer<Mesh<T>> passedMesh) { mesh = passedMesh; }

  void setUpdatePointData(bool update) { updatePointData = update; }

  void setCheckCornersAndEdges(bool check) { checkCornersAndEdges = check; }

  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addError("No level set was passed to ToSurfaceMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      Logger::getInstance()
          .addError("No mesh was passed to ToSurfaceMesh.")
          .print();
      return;
    }

    if (levelSet->getNumberOfPoints() == 0) {
      VIENNACORE_LOG_WARNING(
          "ToSurfaceMesh: Level set is empty. No mesh will be created.");
      return;
    }

    mesh->clear();
    const unsigned int corner0[12] = {0, 1, 2, 0, 4, 5, 6, 4, 0, 1, 3, 2};
    const unsigned int corner1[12] = {1, 3, 3, 2, 5, 7, 7, 6, 4, 5, 7, 6};
    const unsigned int direction[12] = {0, 1, 0, 1, 0, 1, 0, 1, 2, 2, 2, 2};

    // test if level set function consists of at least 2 layers of
    // defined grid points
    if (levelSet->getLevelSetWidth() < 2) {
      VIENNACORE_LOG_WARNING("Levelset is less than 2 layers wide. Expanding "
                             "levelset to 2 layers.");
      Expand<T, D>(levelSet, 2).apply();
    }

    typedef std::map<viennahrle::Index<D>, unsigned> nodeContainerType;

    nodeContainerType nodes[D];
    nodeContainerType faceNodes[D];
    typename nodeContainerType::iterator nodeIt;
    const bool updateData = updatePointData;

    // save how data should be transferred to new level set
    // list of indices into the old pointData vector
    std::vector<std::vector<unsigned>> newDataSourceIds;
    // there is no multithreading here, so just use 1
    if (updateData)
      newDataSourceIds.resize(1);

    // Cache for gradients to avoid re-calculation
    std::vector<Vec3D<T>> nodeNormals;
    if (checkCornersAndEdges) {
      nodeNormals.resize(levelSet->getNumberOfPoints());

      auto &domain = levelSet->getDomain();
      auto &grid = levelSet->getGrid();

      // Parallel pre-calculation of gradients
#pragma omp parallel num_threads(domain.getNumberOfSegments())
      {
        int p = 0;
#ifdef _OPENMP
        p = omp_get_thread_num();
#endif

        viennahrle::Index<D> startVector =
            (p == 0) ? grid.getMinGridPoint() : domain.getSegmentation()[p - 1];
        viennahrle::Index<D> endVector =
            (p != static_cast<int>(domain.getNumberOfSegments() - 1))
                ? domain.getSegmentation()[p]
                : grid.incrementIndices(grid.getMaxGridPoint());

        viennahrle::SparseStarIterator<hrleDomainType, 1> neighborIt(
            domain, startVector);

        for (; neighborIt.getIndices() < endVector; neighborIt.next()) {
          if (neighborIt.getCenter().isDefined()) {
            Vec3D<T> grad{};
            for (int i = 0; i < D; ++i) {
              bool negDefined = neighborIt.getNeighbor(i + D).isDefined();
              bool posDefined = neighborIt.getNeighbor(i).isDefined();

              if (negDefined && posDefined) {
                T valNeg = neighborIt.getNeighbor(i + D).getValue();
                T valCenter = neighborIt.getCenter().getValue();
                T valPos = neighborIt.getNeighbor(i).getValue();

                bool centerSign = valCenter > 0;
                bool negSign = valNeg > 0;
                bool posSign = valPos > 0;

                if (centerSign != negSign && centerSign == posSign) {
                  grad[i] = (valCenter - valNeg) / grid.getGridDelta();
                } else if (centerSign != posSign && centerSign == negSign) {
                  grad[i] = (valPos - valCenter) / grid.getGridDelta();
                } else {
                  grad[i] = (valPos - valNeg) / (2 * grid.getGridDelta());
                }
              } else if (negDefined) {
                grad[i] = (neighborIt.getCenter().getValue() - neighborIt.getNeighbor(i + D).getValue()) / grid.getGridDelta();
              } else if (posDefined) {
                grad[i] = (neighborIt.getNeighbor(i).getValue() - neighborIt.getCenter().getValue()) / grid.getGridDelta();
              } else {
                grad[i] = 0;
              }
            }
            Normalize(grad);
            nodeNormals[neighborIt.getCenter().getPointId()] = grad;
          }
        }
      }
    }

    // iterate over all active points
    for (viennahrle::ConstSparseCellIterator<hrleDomainType> cellIt(
             levelSet->getDomain());
         !cellIt.isFinished(); cellIt.next()) {

      for (int u = 0; u < D; u++) {
        while (!nodes[u].empty() &&
               nodes[u].begin()->first <
                   viennahrle::Index<D>(cellIt.getIndices()))
          nodes[u].erase(nodes[u].begin());

        while (!faceNodes[u].empty() &&
               faceNodes[u].begin()->first <
                   viennahrle::Index<D>(cellIt.getIndices()))
          faceNodes[u].erase(faceNodes[u].begin());
      }

      unsigned signs = 0;
      for (int i = 0; i < (1 << D); i++) {
        if (cellIt.getCorner(i).getValue() >= T(0))
          signs |= (1 << i);
      }

      // all corners have the same sign, so no surface here
      if (signs == 0)
        continue;
      if (signs == (1 << (1 << D)) - 1)
        continue;

      bool perfectCornerFound = false;
      if (checkCornersAndEdges) {
        auto getGradient = [&](int cornerID) {
          auto corner = cellIt.getCorner(cornerID);

          if (corner.isDefined()) {
            return nodeNormals[corner.getPointId()];
          }
          return Vec3D<T>{};
        };

        // Check for perfect corner (2D only)
        if constexpr (D == 2) {
          int insideCount = 0;
          int cornerIdx = -1;
          for (int i = 0; i < 4; ++i) {
            if (signs & (1 << i))
              insideCount++;
          }

          // Identify the "corner" vertex of the cell
          if (insideCount == 1) {
            // Convex corner: the inside vertex is the corner
            for (int i = 0; i < 4; ++i)
              if (signs & (1 << i))
                cornerIdx = i;
          } else if (insideCount == 3) {
            // Concave corner: the outside vertex is the corner
            for (int i = 0; i < 4; ++i)
              if (!(signs & (1 << i)))
                cornerIdx = i;
          }

          if (cornerIdx != -1) {
            int n1 = cornerIdx ^ 1;
            int n2 = cornerIdx ^ 2;

            Vec3D<T> norm1 = getGradient(n1);
            Vec3D<T> norm2 = getGradient(n2);

            // if ((std::abs(DotProduct(norm1, norm2)) < 0.1e-4) && (DotProduct(norm1, norm2) > -1e-8)) {
            // if (DotProduct(norm1, norm2) < 0.9) {
              int edgeX = -1, edgeY = -1;
              if (cornerIdx == 0) {
                edgeX = 0;
                edgeY = 3;
              } else if (cornerIdx == 1) {
                edgeX = 0;
                edgeY = 1;
              } else if (cornerIdx == 2) {
                edgeX = 2;
                edgeY = 3;
              } else if (cornerIdx == 3) {
                edgeX = 2;
                edgeY = 1;
              }

              // Helper to get/create node
              auto getNode = [&](int edge) {
                unsigned p0 = corner0[edge];
                unsigned p1 = corner1[edge];
                auto dir = direction[edge];
                viennahrle::Index<D> d(cellIt.getIndices());
                d += viennahrle::BitMaskToIndex<D>(p0);

                if (nodes[dir].find(d) != nodes[dir].end()) {
                  return nodes[dir][d];
                }

                // Create node
                Vec3D<T> cc{};
                std::size_t currentPointId = 0;
                for (int z = 0; z < D; z++) {
                  if (z != dir) {
                    cc[z] = static_cast<double>(
                        cellIt.getIndices(z) +
                        viennahrle::BitMaskToIndex<D>(p0)[z]);
                  } else {
                    T d0 = cellIt.getCorner(p0).getValue();
                    T d1 = cellIt.getCorner(p1).getValue();
                    if (d0 == -d1) {
                      currentPointId = cellIt.getCorner(p0).getPointId();
                      cc[z] = static_cast<T>(cellIt.getIndices(z)) + 0.5;
                    } else {
                      if (std::abs(d0) <= std::abs(d1)) {
                        currentPointId = cellIt.getCorner(p0).getPointId();
                        cc[z] = static_cast<T>(cellIt.getIndices(z)) +
                                (d0 / (d0 - d1));
                      } else {
                        currentPointId = cellIt.getCorner(p1).getPointId();
                        cc[z] = static_cast<T>(cellIt.getIndices(z) + 1) -
                                (d1 / (d1 - d0));
                      }
                    }
                    cc[z] = std::max(cc[z], cellIt.getIndices(z) + epsilon);
                    cc[z] =
                        std::min(cc[z], (cellIt.getIndices(z) + 1) - epsilon);
                  }
                  cc[z] = levelSet->getGrid().getGridDelta() * cc[z];
                }
                unsigned nodeId = mesh->insertNextNode(cc);
                nodes[dir][d] = nodeId;
                if (updateData)
                  newDataSourceIds[0].push_back(currentPointId);
                return nodeId;
              };

              unsigned nX = getNode(edgeX);
              unsigned nY = getNode(edgeY);

              auto pX = mesh->getNodes()[nX];
              auto pY = mesh->getNodes()[nY];

              double d1 = DotProduct(norm1, pX);
              double d2 = DotProduct(norm2, pY);
              double det = norm1[0] * norm2[1] - norm1[1] * norm2[0];

              Vec3D<T> cornerPos{};
              if (std::abs(det) > 1e-6) {
                cornerPos[0] = (d1 * norm2[1] - d2 * norm1[1]) / det;
                cornerPos[1] = (d2 * norm1[0] - d1 * norm2[0]) / det;
              } else {
                for (int i = 0; i < D; ++i)
                  cornerPos[i] = (pX[i] + pY[i]) * 0.5;
              }
              // cornerPos[2] = 0.0;

              // Check if point is valid (does not cross the linear segment in a wrong way)
              // The intersection point should be on the "void" side of the linear segment (dot > 0)
              Vec3D<T> midPoint{};
              for (int i = 0; i < D; ++i)
                  midPoint[i] = (pX[i] + pY[i]) * 0.5;
              Vec3D<T> avgNorm = norm1 + norm2;
              T dot = DotProduct(cornerPos - midPoint, avgNorm);

              if (dot > 0) {
                perfectCornerFound = true;
                unsigned nCorner = mesh->insertNextNode(cornerPos);
                if (updateData)
                  newDataSourceIds[0].push_back(newDataSourceIds[0].back());

                // Add lines
                mesh->insertNextElement(std::array<unsigned, 2>{nX, nCorner});
                mesh->insertNextElement(std::array<unsigned, 2>{nCorner, nY});
              }
            // }
          }
        } else if constexpr (D == 3) {
          int insideCount = 0;
          int cornerIdx = -1;
          for (int i = 0; i < 8; ++i) {
            if (signs & (1 << i))
              insideCount++;
          }

          bool isConvex = false;
          if (insideCount == 1) {
            isConvex = true;
            for (int i = 0; i < 8; ++i)
              if (signs & (1 << i))
                cornerIdx = i;
          } else if (insideCount == 7) {
            isConvex = false;
            for (int i = 0; i < 8; ++i)
              if (!(signs & (1 << i)))
                cornerIdx = i;
          }

          if (cornerIdx != -1) {
            // Neighbors of cornerIdx in the cube
            int n1 = cornerIdx ^ 1;
            int n2 = cornerIdx ^ 2;
            int n4 = cornerIdx ^ 4;

            Vec3D<T> normX = getGradient(n1);
            Vec3D<T> normY = getGradient(n2);
            Vec3D<T> normZ = getGradient(n4);

            // Helper to get/create node on edge
            auto getNode = [&](unsigned pA, unsigned pB) {
              int edge = -1;
              for (int e = 0; e < 12; ++e) {
                if ((corner0[e] == pA && corner1[e] == pB) ||
                    (corner0[e] == pB && corner1[e] == pA)) {
                  edge = e;
                  break;
                }
              }
              unsigned p0 = corner0[edge];
              unsigned p1 = corner1[edge];
              auto dir = direction[edge];
              viennahrle::Index<D> d(cellIt.getIndices());
              d += viennahrle::BitMaskToIndex<D>(p0);

              if (nodes[dir].find(d) != nodes[dir].end()) {
                return nodes[dir][d];
              }

              Vec3D<T> cc{};
              std::size_t currentPointId = 0;
              for (int z = 0; z < D; z++) {
                if (z != dir) {
                  cc[z] = static_cast<double>(
                      cellIt.getIndices(z) +
                      viennahrle::BitMaskToIndex<D>(p0)[z]);
                } else {
                  T d0 = cellIt.getCorner(p0).getValue();
                  T d1 = cellIt.getCorner(p1).getValue();
                  if (d0 == -d1) {
                    currentPointId = cellIt.getCorner(p0).getPointId();
                    cc[z] = static_cast<T>(cellIt.getIndices(z)) + 0.5;
                  } else {
                    if (std::abs(d0) <= std::abs(d1)) {
                      currentPointId = cellIt.getCorner(p0).getPointId();
                      cc[z] = static_cast<T>(cellIt.getIndices(z)) +
                              (d0 / (d0 - d1));
                    } else {
                      currentPointId = cellIt.getCorner(p1).getPointId();
                      cc[z] = static_cast<T>(cellIt.getIndices(z) + 1) -
                              (d1 / (d1 - d0));
                    }
                  }
                  cc[z] = std::max(cc[z], cellIt.getIndices(z) + epsilon);
                  cc[z] = std::min(cc[z], (cellIt.getIndices(z) + 1) - epsilon);
                }
                cc[z] = levelSet->getGrid().getGridDelta() * cc[z];
              }
              unsigned nodeId = mesh->insertNextNode(cc);
              nodes[dir][d] = nodeId;
              if (updateData)
                newDataSourceIds[0].push_back(currentPointId);
              return nodeId;
            };

            unsigned nX = getNode(cornerIdx, n1);
            unsigned nY = getNode(cornerIdx, n2);
            unsigned nZ = getNode(cornerIdx, n4);

            auto pX = mesh->getNodes()[nX];
            auto pY = mesh->getNodes()[nY];
            auto pZ = mesh->getNodes()[nZ];

            T d1 = DotProduct(normX, pX);
            T d2 = DotProduct(normY, pY);
            T d3 = DotProduct(normZ, pZ);

            auto intersectPlanes =
                [](Vec3D<T> n1, T D1, Vec3D<T> n2, T D2, Vec3D<T> n3,
                   T D3) -> std::pair<bool, Vec3D<T>> {
              Vec3D<T> cross23 = CrossProduct(n2, n3);
              T det = DotProduct(n1, cross23);
              if (std::abs(det) < 1e-6)
                return {false, {}};
              Vec3D<T> p = (D1 * cross23 + D2 * CrossProduct(n3, n1) +
                            D3 * CrossProduct(n1, n2)) /
                           det;
              return {true, p};
            };

            // Face Plane X: x = cellPos[0] + (corner is 1 ? gridDelta : 0)
            Vec3D<T> faceNormX = {1, 0, 0};
            T faceDistX = (cellIt.getIndices(0) + ((cornerIdx & 1) ? 1 : 0)) *
                          levelSet->getGrid().getGridDelta();

            Vec3D<T> faceNormY = {0, 1, 0};
            T faceDistY = (cellIt.getIndices(1) + ((cornerIdx & 2) ? 1 : 0)) *
                          levelSet->getGrid().getGridDelta();

            Vec3D<T> faceNormZ = {0, 0, 1};
            T faceDistZ = (cellIt.getIndices(2) + ((cornerIdx & 4) ? 1 : 0)) *
                          levelSet->getGrid().getGridDelta();

            // F_yz: On Face X. Intersection of Plane Y and Plane Z.
            auto resF_yz =
                intersectPlanes(normY, d2, normZ, d3, faceNormX, faceDistX);
            // F_xz: On Face Y. Intersection of Plane X and Plane Z.
            auto resF_xz =
                intersectPlanes(normX, d1, normZ, d3, faceNormY, faceDistY);
            // F_xy: On Face Z. Intersection of Plane X and Plane Y.
            auto resF_xy =
                intersectPlanes(normX, d1, normY, d2, faceNormZ, faceDistZ);

            Vec3D<T> cross23 = CrossProduct(normY, normZ);
            T det = DotProduct(normX, cross23);

            bool isOrthogonal = (std::abs(DotProduct(normX, normY)) < 0.1) &&
                                (std::abs(DotProduct(normX, normZ)) < 0.1) &&
                                (std::abs(DotProduct(normY, normZ)) < 0.1);

            if (isOrthogonal && std::abs(det) > 1e-3 && resF_yz.first &&
                resF_xz.first &&
                resF_xy.first) {
              Vec3D<T> cornerPos =
                  (d1 * cross23 + d2 * CrossProduct(normZ, normX) +
                   d3 * CrossProduct(normX, normY)) /
                  det;

              perfectCornerFound = true;
              unsigned nCorner = mesh->insertNextNode(cornerPos);
              if (updateData)
                newDataSourceIds[0].push_back(newDataSourceIds[0].back());

              auto handleFaceNode = [&](int axis, Vec3D<T> pos) {
                viennahrle::Index<D> fIdx(cellIt.getIndices());
                if (cornerIdx & (1 << axis))
                  fIdx[axis]++;

                if (faceNodes[axis].find(fIdx) != faceNodes[axis].end()) {
                  return faceNodes[axis][fIdx];
                }
                unsigned nid = mesh->insertNextNode(pos);
                faceNodes[axis][fIdx] = nid;
                if (updateData)
                  newDataSourceIds[0].push_back(newDataSourceIds[0].back());
                return nid;
              };

              unsigned nF_yz = handleFaceNode(0, resF_yz.second);
              unsigned nF_xz = handleFaceNode(1, resF_xz.second);
              unsigned nF_xy = handleFaceNode(2, resF_xy.second);

              auto addQuad = [&](unsigned a, unsigned b, unsigned c,
                                 unsigned d) {
                if (isConvex) {
                  mesh->insertNextElement(std::array<unsigned, 3>{a, b, c});
                  mesh->insertNextElement(std::array<unsigned, 3>{a, c, d});
                } else {
                  mesh->insertNextElement(std::array<unsigned, 3>{a, c, b});
                  mesh->insertNextElement(std::array<unsigned, 3>{a, d, c});
                }
              };

              addQuad(nCorner, nF_xy, nX, nF_xz); // Plane X
              addQuad(nCorner, nF_yz, nY, nF_xy); // Plane Y
              addQuad(nCorner, nF_xz, nZ, nF_yz); // Plane Z
            }
          }

          int edgeAxis = -1;
          bool isConvexEdge = false;

          if (insideCount == 2) {
            int c1 = -1, c2 = -1;
            for (int i = 0; i < 8; ++i) {
              if (signs & (1 << i)) {
                if (c1 == -1)
                  c1 = i;
                else
                  c2 = i;
              }
            }
            int diff = c1 ^ c2;
            if (diff == 1)
              edgeAxis = 0;
            else if (diff == 2)
              edgeAxis = 1;
            else if (diff == 4)
              edgeAxis = 2;

            if (edgeAxis != -1)
              isConvexEdge = true;
          } else if (insideCount == 6) {
            int c1 = -1, c2 = -1;
            for (int i = 0; i < 8; ++i) {
              if (!(signs & (1 << i))) {
                if (c1 == -1)
                  c1 = i;
                else
                  c2 = i;
              }
            }
            int diff = c1 ^ c2;
            if (diff == 1)
              edgeAxis = 0;
            else if (diff == 2)
              edgeAxis = 1;
            else if (diff == 4)
              edgeAxis = 2;

            if (edgeAxis != -1)
              isConvexEdge = false;
          }

          if (edgeAxis != -1) {
            auto getNode = [&](unsigned pA, unsigned pB) {
              int edge = -1;
              for (int e = 0; e < 12; ++e) {
                if ((corner0[e] == pA && corner1[e] == pB) ||
                    (corner0[e] == pB && corner1[e] == pA)) {
                  edge = e;
                  break;
                }
              }
              unsigned p0 = corner0[edge];
              unsigned p1 = corner1[edge];
              auto dir = direction[edge];
              viennahrle::Index<D> d(cellIt.getIndices());
              d += viennahrle::BitMaskToIndex<D>(p0);

              if (nodes[dir].find(d) != nodes[dir].end()) {
                return nodes[dir][d];
              }

              Vec3D<T> cc{};
              std::size_t currentPointId = 0;
              for (int z = 0; z < D; z++) {
                if (z != dir) {
                  cc[z] = static_cast<double>(
                      cellIt.getIndices(z) +
                      viennahrle::BitMaskToIndex<D>(p0)[z]);
                } else {
                  T d0 = cellIt.getCorner(p0).getValue();
                  T d1 = cellIt.getCorner(p1).getValue();
                  if (d0 == -d1) {
                    currentPointId = cellIt.getCorner(p0).getPointId();
                    cc[z] = static_cast<T>(cellIt.getIndices(z)) + 0.5;
                  } else {
                    if (std::abs(d0) <= std::abs(d1)) {
                      currentPointId = cellIt.getCorner(p0).getPointId();
                      cc[z] = static_cast<T>(cellIt.getIndices(z)) +
                              (d0 / (d0 - d1));
                    } else {
                      currentPointId = cellIt.getCorner(p1).getPointId();
                      cc[z] = static_cast<T>(cellIt.getIndices(z) + 1) -
                              (d1 / (d1 - d0));
                    }
                  }
                  cc[z] = std::max(cc[z], cellIt.getIndices(z) + epsilon);
                  cc[z] =
                      std::min(cc[z], (cellIt.getIndices(z) + 1) - epsilon);
                }
                cc[z] = levelSet->getGrid().getGridDelta() * cc[z];
              }
              unsigned nodeId = mesh->insertNextNode(cc);
              nodes[dir][d] = nodeId;
              if (updateData)
                newDataSourceIds[0].push_back(currentPointId);
              return nodeId;
            };

            std::map<int, unsigned> nodes0, nodes1;
            auto checkSlice =
                [&](int sliceIdx,
                    std::map<int, unsigned> &nodesMap) -> std::pair<bool, Vec3D<T>> {
              int corner = -1;
              for (int i = 0; i < 8; ++i) {
                if (((i >> edgeAxis) & 1) != sliceIdx)
                  continue;
                bool isInside = (signs & (1 << i));
                if (isConvexEdge == isInside)
                  corner = i;
              }
              if (corner == -1)
                return {false, {}};

              int axes[2];
              int idx = 0;
              for (int axis = 0; axis < 3; ++axis) {
                if (axis != edgeAxis)
                  axes[idx++] = axis;
              }

              int n1 = corner ^ (1 << axes[0]);
              int n2 = corner ^ (1 << axes[1]);

              Vec3D<T> g1 = getGradient(n1);
              Vec3D<T> g2 = getGradient(n2);

              if (std::abs(DotProduct(g1, g2)) < 0.1) {
                unsigned node1 = getNode(corner, n1);
                unsigned node2 = getNode(corner, n2);
                nodesMap[axes[0]] = node1;
                nodesMap[axes[1]] = node2;

                Vec3D<T> p1 = mesh->getNodes()[node1];
                Vec3D<T> p2 = mesh->getNodes()[node2];

                double det =
                    g1[axes[0]] * g2[axes[1]] - g1[axes[1]] * g2[axes[0]];
                if (std::abs(det) < 1e-6)
                  return {false, {}};

                double val1 =
                    g1[axes[0]] * p1[axes[0]] + g1[axes[1]] * p1[axes[1]];
                double val2 =
                    g2[axes[0]] * p2[axes[0]] + g2[axes[1]] * p2[axes[1]];

                Vec3D<T> res;
                res[edgeAxis] = p1[edgeAxis];
                res[axes[0]] = (val1 * g2[axes[1]] - val2 * g1[axes[1]]) / det;
                res[axes[1]] = (val2 * g1[axes[0]] - val1 * g2[axes[0]]) / det;
                return {true, res};
              }
              return {false, {}};
            };

            auto res0 = checkSlice(0, nodes0);
            auto res1 = checkSlice(1, nodes1);

            if (res0.first && res1.first) {
              perfectCornerFound = true;

              auto handleFaceNode = [&](int sliceIdx, Vec3D<T> pos) {
                viennahrle::Index<D> fIdx(cellIt.getIndices());
                if (sliceIdx == 1)
                  fIdx[edgeAxis]++;

                if (faceNodes[edgeAxis].find(fIdx) != faceNodes[edgeAxis].end()) {
                  return faceNodes[edgeAxis][fIdx];
                }
                unsigned nid = mesh->insertNextNode(pos);
                faceNodes[edgeAxis][fIdx] = nid;
                if (updateData)
                  newDataSourceIds[0].push_back(newDataSourceIds[0].back());
                return nid;
              };

              unsigned s0 = handleFaceNode(0, res0.second);
              unsigned s1 = handleFaceNode(1, res1.second);

              for (auto const &[axis, n0] : nodes0) {
                unsigned n1 = nodes1[axis];
                Vec3D<T> pA = mesh->getNodes()[s0];
                Vec3D<T> pB = mesh->getNodes()[n0];
                Vec3D<T> pC = mesh->getNodes()[n1];
                Vec3D<T> center;
                for (int k = 0; k < 3; ++k)
                  center[k] = (cellIt.getIndices(k) + 0.5) *
                              levelSet->getGrid().getGridDelta();

                Vec3D<T> n = CrossProduct(pB - pA, pC - pA);
                Vec3D<T> dir = pA - center;

                if (DotProduct(n, dir) < 0) {
                  mesh->insertNextElement(std::array<unsigned, 3>{s0, n1, n0});
                  mesh->insertNextElement(std::array<unsigned, 3>{s0, s1, n1});
                } else {
                  mesh->insertNextElement(std::array<unsigned, 3>{s0, n0, n1});
                  mesh->insertNextElement(std::array<unsigned, 3>{s0, n1, s1});
                }
              }
            }
          }
        }
      }

      if (perfectCornerFound)
        continue;


      // for each element
      const int *Triangles =
          (D == 2) ? lsInternal::MarchingCubes::polygonize2d(signs)
                   : lsInternal::MarchingCubes::polygonize3d(signs);

      for (; Triangles[0] != -1; Triangles += D) {
        std::array<unsigned, D> nod_numbers;

        // for each node
        for (int n = 0; n < D; n++) {
          const int edge = Triangles[n];

          unsigned p0 = corner0[edge];
          unsigned p1 = corner1[edge];

          // determine direction of edge
          auto dir = direction[edge];

          // look for existing surface node
          viennahrle::Index<D> d(cellIt.getIndices());
          d += viennahrle::BitMaskToIndex<D>(p0);

          nodeIt = nodes[dir].find(d);
          if (nodeIt != nodes[dir].end()) {
            nod_numbers[n] = nodeIt->second;
          } else { // if node does not exist yet

            // calculate coordinate of new node
            Vec3D<T> cc{}; // initialise with zeros
            std::size_t currentPointId = 0;
            for (int z = 0; z < D; z++) {
              if (z != dir) {
                // TODO might not need BitMaskToVector here, just check if z bit
                // is set
                cc[z] =
                    static_cast<double>(cellIt.getIndices(z) +
                                        viennahrle::BitMaskToIndex<D>(p0)[z]);
              } else {
                T d0, d1;

                d0 = cellIt.getCorner(p0).getValue();
                d1 = cellIt.getCorner(p1).getValue();

                // calculate the surface-grid intersection point
                if (d0 == -d1) { // includes case where d0=d1=0
                  currentPointId = cellIt.getCorner(p0).getPointId();
                  cc[z] = static_cast<T>(cellIt.getIndices(z)) + 0.5;
                } else {
                  if (std::abs(d0) <= std::abs(d1)) {
                    currentPointId = cellIt.getCorner(p0).getPointId();
                    cc[z] =
                        static_cast<T>(cellIt.getIndices(z)) + (d0 / (d0 - d1));
                  } else {
                    currentPointId = cellIt.getCorner(p1).getPointId();
                    cc[z] = static_cast<T>(cellIt.getIndices(z) + 1) -
                            (d1 / (d1 - d0));
                  }
                }
                cc[z] = std::max(cc[z], cellIt.getIndices(z) + epsilon);
                cc[z] = std::min(cc[z], (cellIt.getIndices(z) + 1) - epsilon);
              }
              cc[z] = levelSet->getGrid().getGridDelta() * cc[z];
            }

            // insert new node
            nod_numbers[n] =
                mesh->insertNextNode(cc); // insert new surface node
            nodes[dir][d] = nod_numbers[n];

            if (updateData)
              newDataSourceIds[0].push_back(currentPointId);
          }
        }

        if (!triangleMisformed(nod_numbers))
          mesh->insertNextElement(nod_numbers); // insert new surface element
      }
    }

    // now copy old data into new level set
    if (updateData) {
      mesh->getPointData().translateFromMultiData(levelSet->getPointData(),
                                                  newDataSourceIds);
    }
  }

private:
  static bool triangleMisformed(const std::array<unsigned, D> &nodeNumbers) {
    if constexpr (D == 3) {
      return nodeNumbers[0] == nodeNumbers[1] ||
             nodeNumbers[0] == nodeNumbers[2] ||
             nodeNumbers[1] == nodeNumbers[2];
    } else {
      return nodeNumbers[0] == nodeNumbers[1];
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(ToSurfaceMesh)

} // namespace viennals
