#pragma once

#include <lsPreCompileMacros.hpp>

#include <map>

#include <hrleSparseCellIterator.hpp>
#include <hrleSparseStarIterator.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFiniteDifferences.hpp>
#include <lsMarchingCubes.hpp>
#include <lsMesh.hpp>
#include <lsCalculateNormalVectors.hpp>

namespace viennals {

using namespace viennacore;

/// Extract an explicit Mesh<> instance from an lsDomain.
/// The interface is then described by explicit surface elements:
/// Lines in 2D, Triangles in 3D.
template <class T, int D> class ToSurfaceMesh {
  typedef typename Domain<T, D>::DomainType hrleDomainType;

  SmartPointer<Domain<T, D>> levelSet = nullptr;
  SmartPointer<Mesh<T>> mesh = nullptr;
  const T epsilon;
  bool updatePointData = true;
  bool sharpCorners = true;

public:
  explicit ToSurfaceMesh(double eps = 1e-12) : epsilon(eps) {}

  ToSurfaceMesh(const SmartPointer<Domain<T, D>> passedLevelSet,
                SmartPointer<Mesh<T>> passedMesh, double eps = 1e-12)
      : levelSet(passedLevelSet), mesh(passedMesh), epsilon(eps) {
        // if constexpr (D == 3) {
        //   sharpCorners = false;
        // }
      }

  void setLevelSet(SmartPointer<Domain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  void setMesh(SmartPointer<Mesh<T>> passedMesh) { mesh = passedMesh; }

  void setUpdatePointData(bool update) { updatePointData = update; }

  void setSharpCorners(bool check) { sharpCorners = check; }

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

    // test if level set function consists of at least 2 layers of
    // defined grid points
    if (levelSet->getLevelSetWidth() < 2) {
      VIENNACORE_LOG_WARNING("Levelset is less than 2 layers wide. Expanding "
                             "levelset to 2 layers.");
      Expand<T, D>(levelSet, 2).apply();
    }

    viennahrle::ConstSparseIterator<hrleDomainType> valueIt(levelSet->getDomain());

    typedef std::map<viennahrle::Index<D>, unsigned> nodeContainerType;

    nodeContainerType nodes[D];
    nodeContainerType faceNodes[D];
    nodeContainerType cornerNodes;
    typename nodeContainerType::iterator nodeIt;

    // save how data should be transferred to new level set
    // list of indices into the old pointData vector
    std::vector<std::vector<unsigned>> newDataSourceIds;
    // there is no multithreading here, so just use 1
    if (updatePointData)
      newDataSourceIds.resize(1);

    typename PointData<T>::VectorDataType *normalVectorData = nullptr;
    if (sharpCorners) {
      viennals::CalculateNormalVectors<T, D> normalCalculator(levelSet);
      normalCalculator.setMethod(
          viennals::lsNormalCalculationMethodEnum::ONE_SIDED_MIN_MOD);
      normalCalculator.setMaxValue(std::numeric_limits<T>::max());
      normalCalculator.apply();
      normalVectorData = levelSet->getPointData().getVectorData(
          viennals::CalculateNormalVectors<T, D>::normalVectorsLabel);
    }

    // iterate over all cells with active points
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

        if (u == 0) {
          while (!cornerNodes.empty() && cornerNodes.begin()->first < viennahrle::Index<D>(cellIt.getIndices()))
            cornerNodes.erase(cornerNodes.begin());
        }
      }

      unsigned signs = 0;
      bool hasZero = false;
      for (int i = 0; i < (1 << D); i++) {
        T val = cellIt.getCorner(i).getValue();
        if (val >= T(0))
          signs |= (1 << i);
        if (std::abs(val) <= epsilon)
          hasZero = true;
      }

      // all corners have the same sign, so no surface here
      if (signs == 0)
        continue;
      if (signs == (1 << (1 << D)) - 1 && !hasZero)
        continue;
      
      bool perfectCornerFound = false;
      if (sharpCorners) {
        auto getGradient = [&](int cornerID) {
          auto corner = cellIt.getCorner(cornerID);

          if (corner.isDefined()) {
            return (*normalVectorData)[corner.getPointId()];
          }
          return Vec3D<T>{};
        };

        // Check for perfect corner (2D only)
        if constexpr (D == 2) {
          perfectCornerFound = generateSharpCorner2D(
              cellIt, nodes, newDataSourceIds, valueIt);
        } else if constexpr (D == 3) {
          perfectCornerFound = generateSharpEdge3D(
              cellIt, nodes, faceNodes, newDataSourceIds, valueIt);
          if (!perfectCornerFound) {
            perfectCornerFound = generateSharpCorner3D(
                cellIt, nodes, cornerNodes, faceNodes, newDataSourceIds, valueIt);
          }
          if (!perfectCornerFound) {
            perfectCornerFound = generateSharpL3D(
                cellIt, nodes, faceNodes, newDataSourceIds, valueIt);
          }
        }
      }

      if (perfectCornerFound)
        continue;



      if constexpr (D == 3) {
        // Stitch to perfect corners/edges
        for (int axis = 0; axis < 3; ++axis) {
          for (int d = 0; d < 2; ++d) {
            viennahrle::Index<D> faceKey(cellIt.getIndices());
            if (d == 1)
              faceKey[axis]++;

            auto it = faceNodes[axis].find(faceKey);
            if (it != faceNodes[axis].end()) {
              const unsigned faceNodeId = it->second;
              const int *Triangles = lsInternal::MarchingCubes::polygonize3d(signs);

              auto getNode = [&](int edge) -> unsigned {
                return this->getNode(cellIt, edge, nodes,
                                     newDataSourceIds);
              };

              for (; Triangles[0] != -1; Triangles += 3) {
                std::vector<unsigned> face_edge_nodes;
                for (int i = 0; i < 3; ++i) {
                  int edge = Triangles[i];
                  int c0 = corner0[edge];
                  int c1 = corner1[edge];
                  bool onFace = (((c0 >> axis) & 1) == d) && (((c1 >> axis) & 1) == d);
                  if (onFace) {
                    face_edge_nodes.push_back(getNode(edge));
                  }
                }
                if (face_edge_nodes.size() == 2) {
                  mesh->insertNextElement(std::array<unsigned, 3>{
                      face_edge_nodes[0], face_edge_nodes[1], faceNodeId});
                }
              }
            }
          }
        }
      }

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
            nod_numbers[n] =
                this->getNode(cellIt, edge, nodes, newDataSourceIds);
          }
        }

        if (!triangleMisformed(nod_numbers))
          mesh->insertNextElement(nod_numbers); // insert new surface element
      }
    }

    // now copy old data into new level set
    if (updatePointData) {
      mesh->getPointData().translateFromMultiData(levelSet->getPointData(),
                                                  newDataSourceIds);
    }
  }

private:
  unsigned getNode(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt, int edge,
      std::map<viennahrle::Index<D>, unsigned> *nodes,
      std::vector<std::vector<unsigned>> &newDataSourceIds) {

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
        cc[z] = static_cast<double>(cellIt.getIndices(z) +
                                    viennahrle::BitMaskToIndex<D>(p0)[z]);
      } else {
        T d0 = cellIt.getCorner(p0).getValue();
        T d1 = cellIt.getCorner(p1).getValue();
        if (d0 == -d1) {
          currentPointId = cellIt.getCorner(p0).getPointId();
          cc[z] = static_cast<T>(cellIt.getIndices(z)) + T(0.5);
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
        cc[z] =
            std::min(cc[z], (cellIt.getIndices(z) + 1) - epsilon);
      }
      cc[z] = levelSet->getGrid().getGridDelta() * cc[z];
    }
    unsigned nodeId = mesh->insertNextNode(cc);
    nodes[dir][d] = nodeId;
    if (updatePointData)
      newDataSourceIds[0].push_back(currentPointId);
    return nodeId;
  }

  void stitchToNeighbor(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt,
      int axis, bool isHighFace, unsigned faceNodeId,
      std::map<viennahrle::Index<D>, unsigned> *nodes,
      viennahrle::ConstSparseIterator<hrleDomainType> &valueIt) {
    // Backward stitching: Check if neighbor on this face is "past" and needs stitching
    viennahrle::Index<D> neighborIdx = cellIt.getIndices();
    if (isHighFace) neighborIdx[axis]++; else neighborIdx[axis]--;

    if (neighborIdx < cellIt.getIndices()) {
      unsigned nSigns = 0;
      auto &grid = levelSet->getGrid();
      for(int i=0; i<8; ++i) {
        viennahrle::Index<D> cIdx = neighborIdx + viennahrle::BitMaskToIndex<D>(i);
        if (!grid.isOutsideOfDomain(cIdx)) {
          valueIt.goToIndices(cIdx);
          if (valueIt.getValue() >= 0) nSigns |= (1<<i);
        }
      }

      if (nSigns != 0 && nSigns != 255) {
        const int *nTriangles = lsInternal::MarchingCubes::polygonize3d(nSigns);
        int nFaceD = isHighFace ? 0 : 1;

        for (; nTriangles[0] != -1; nTriangles += 3) {
          std::vector<unsigned> face_edge_nodes;
          for (int i = 0; i < 3; ++i) {
            int edge = nTriangles[i];
            int c0 = corner0[edge];
            int c1 = corner1[edge];
            bool onFace = (((c0 >> axis) & 1) == nFaceD) && (((c1 >> axis) & 1) == nFaceD);
            if (onFace) {
              unsigned p0 = corner0[edge];
              auto dir = direction[edge];
              viennahrle::Index<D> d = neighborIdx + viennahrle::BitMaskToIndex<D>(p0);
              auto itN = nodes[dir].find(d);
              if (itN != nodes[dir].end()) {
                face_edge_nodes.push_back(itN->second);
              }
            }
          }
          if (face_edge_nodes.size() == 2) {
            mesh->insertNextElement(std::array<unsigned, 3>{
                face_edge_nodes[0], face_edge_nodes[1], faceNodeId});
          }
        }
      }
    }
  }

  const unsigned int corner0[12] = {0, 1, 2, 0, 4, 5, 6, 4, 0, 1, 3, 2};
  const unsigned int corner1[12] = {1, 3, 3, 2, 5, 7, 7, 6, 4, 5, 7, 6};
  const unsigned int direction[12] = {0, 1, 0, 1, 0, 1, 0, 1, 2, 2, 2, 2};

  bool generateCanonicalSharpCorner2D(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt,
      int transform, Vec3D<T> &cornerPos) {

    const auto *normalVectorData = levelSet->getPointData().getVectorData(
        viennals::CalculateNormalVectors<T, D>::normalVectorsLabel);

    if (normalVectorData == nullptr)
      return false;

    auto getTransformedGradient = [&](int cornerID) {
      auto corner = cellIt.getCorner(cornerID ^ transform);
      if (corner.isDefined()) {
        auto normal = (*normalVectorData)[corner.getPointId()];
        if ((transform & 1) != 0)
          normal[0] = -normal[0];
        if ((transform & 2) != 0)
          normal[1] = -normal[1];
        return normal;
      }
      return Vec3D<T>{};
    };

    Vec3D<T> norm1 = getTransformedGradient(1); // neighbor in x
    Vec3D<T> norm2 = getTransformedGradient(2); // neighbor in y

    if (std::abs(DotProduct(norm1, norm2)) >= 0.8) {
      return false;
    }

    auto getValue = [&](int cornerIdx) {
      return cellIt.getCorner(cornerIdx ^ transform).getValue();
    };

    auto calculateNodePos = [&](int edge, Vec3D<T> &pos) {
      unsigned p0 = corner0[edge];
      unsigned p1 = corner1[edge];
      auto dir = direction[edge];

      for (int z = 0; z < D; z++) {
        if (z != dir) {
          pos[z] = static_cast<T>(viennahrle::BitMaskToIndex<D>(p0)[z]);
        } else {
          T d0 = getValue(p0);
          T d1 = getValue(p1);
          if (d0 == -d1) {
            pos[z] = T(0.5);
          } else {
            if (std::abs(d0) <= std::abs(d1)) {
              pos[z] = (d0 / (d0 - d1));
            } else {
              pos[z] = T(1.0) - (d1 / (d1 - d0));
            }
          }
          pos[z] = std::max(pos[z], epsilon);
          pos[z] = std::min(pos[z], T(1.0) - epsilon);
        }
      }
    };

    Vec3D<T> pX, pY;
    calculateNodePos(0, pX); // Edge along x-axis from corner 0
    calculateNodePos(3, pY); // Edge along y-axis from corner 0

    double d1 = DotProduct(norm1, pX);
    double d2 = DotProduct(norm2, pY);
    double det = norm1[0] * norm2[1] - norm1[1] * norm2[0];

    if (std::abs(det) > 1e-6) {
      cornerPos[0] = (d1 * norm2[1] - d2 * norm1[1]) / det;
      cornerPos[1] = (d2 * norm1[0] - d1 * norm2[0]) / det;
    } else {
      for (int i = 0; i < D; ++i)
        cornerPos[i] = (pX[i] + pY[i]) * T(0.5);
    }

    return true;
  }

  bool generateSharpCorner2D(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt,
      std::map<viennahrle::Index<D>, unsigned> *nodes,
      std::vector<std::vector<unsigned>> &newDataSourceIds,
      viennahrle::ConstSparseIterator<hrleDomainType> &valueIt) {

    int countNeg = 0;
    int countPos = 0;
    int countZero = 0;
    int negIdx = -1;
    int posIdx = -1;

    for (int i = 0; i < 4; ++i) {
      T val = cellIt.getCorner(i).getValue();
      if (val < -epsilon) {
        countNeg++;
        negIdx = i;
      } else if (val > epsilon) {
        countPos++;
        posIdx = i;
      } else {
        countZero++;
      }
    }

    int cornerIdx = -1;

    if (countNeg == 1 && (countPos + countZero == 3)) {
      cornerIdx = negIdx;
    } else if (countPos == 1 && (countNeg + countZero == 3)) {
      cornerIdx = posIdx;
    }

    if (cornerIdx != -1) {
      // Check if this corner is also a corner in any neighboring cell
      bool isSharedCorner = false;
      viennahrle::Index<D> pIdx =
          cellIt.getIndices() + viennahrle::BitMaskToIndex<D>(cornerIdx);
      T pVal = cellIt.getCorner(cornerIdx).getValue();
      auto &grid = levelSet->getGrid();

      for (int i = 0; i < (1 << D); ++i) {
        if (i == cornerIdx)
          continue;

        viennahrle::Index<D> neighborIndices;
        for (int k = 0; k < D; ++k)
          neighborIndices[k] = pIdx[k] - ((i >> k) & 1);

        bool neighborIsCorner = true;

        // Check edge-connected neighbors in the neighbor cell
        for (int k = 0; k < D; ++k) {
          int neighborLocal = i ^ (1 << k);
          viennahrle::Index<D> checkIdx =
              neighborIndices +
              viennahrle::BitMaskToIndex<D>(neighborLocal);
          if (grid.isOutsideOfDomain(checkIdx)) {
            checkIdx = grid.globalIndices2LocalIndices(checkIdx);
          }
          valueIt.goToIndices(checkIdx);
          T nVal = valueIt.getValue();
          if (((pVal > epsilon) && (nVal > epsilon)) ||
              ((pVal < -epsilon) && (nVal < -epsilon))) {
            neighborIsCorner = false;
            break;
          }
        }

        if (neighborIsCorner) {
          isSharedCorner = true;
          break;
        }
      }

      if (!isSharedCorner) {
        Vec3D<T> cornerPos;
        if (generateCanonicalSharpCorner2D(cellIt, cornerIdx,
                                           cornerPos)) {

          // inverse transform cornerPos from canonical local to this cell's
          // local
          if ((cornerIdx & 1) != 0)
            cornerPos[0] = T(1.0) - cornerPos[0];
          if ((cornerIdx & 2) != 0)
            cornerPos[1] = T(1.0) - cornerPos[1];

          // convert to global coordinates
          for (int i = 0; i < D; ++i) {
            cornerPos[i] = (cornerPos[i] + cellIt.getIndices(i)) *
                           levelSet->getGrid().getGridDelta();
          }

          // Determine edges for this corner
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

          unsigned nX =
              this->getNode(cellIt, edgeX, nodes, newDataSourceIds);
          unsigned nY =
              this->getNode(cellIt, edgeY, nodes, newDataSourceIds);

          unsigned nCorner = mesh->insertNextNode(cornerPos);
          if (updatePointData)
            newDataSourceIds[0].push_back(
                newDataSourceIds[0].back()); // TODO: improve point data source

          mesh->insertNextElement(std::array<unsigned, 2>{nX, nCorner});
          mesh->insertNextElement(std::array<unsigned, 2>{nCorner, nY});
          return true;
        }
      }
    }
    return false;
  }


  bool generateSharpL3D(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt,
      std::map<viennahrle::Index<D>, unsigned> *nodes,
      std::map<viennahrle::Index<D>, unsigned> *faceNodes,
      std::vector<std::vector<unsigned>> &newDataSourceIds,
      viennahrle::ConstSparseIterator<hrleDomainType> &valueIt) {

    int countNeg = 0;
    int countPos = 0;
    int negMask = 0;
    int posMask = 0;

    for (int i = 0; i < 8; ++i) {
      T val = cellIt.getCorner(i).getValue();
      if (val < -epsilon) {
        countNeg++;
        negMask |= (1 << i);
      } else if (val > epsilon) {
        countPos++;
        posMask |= (1 << i);
      }
    }

    bool inverted = false;
    int mask = 0;
    if (countNeg == 3) {
      mask = negMask;
      inverted = false;
    } else if (countPos == 3) {
      mask = posMask;
      inverted = true;
    } else {
      return false;
    }

    // Check for L-shape: 3 points on a face
    int C = -1, A = -1, B = -1;
    for (int i = 0; i < 8; ++i) {
      if ((mask >> i) & 1) {
        int neighbors = 0;
        for (int k = 0; k < 3; ++k) {
          if ((mask >> (i ^ (1 << k))) & 1) neighbors++;
        }
        if (neighbors == 2) {
          C = i;
          break;
        }
      }
    }

    if (C == -1) return false;

    // Identify A and B
    for (int k = 0; k < 3; ++k) {
      int n = C ^ (1 << k);
      if ((mask >> n) & 1) {
        if (A == -1) A = n;
        else B = n;
      }
    }

    const auto *normalVectorData = levelSet->getPointData().getVectorData(
        viennals::CalculateNormalVectors<T, D>::normalVectorsLabel);
    if (!normalVectorData) return false;

    auto getNormal = [&](int idx) {
      auto corner = cellIt.getCorner(idx);
      if (corner.isDefined()) {
        Vec3D<T> n = (*normalVectorData)[corner.getPointId()];
        if (inverted) for (int i = 0; i < 3; ++i) n[i] = -n[i];
        return n;
      }
      return Vec3D<T>{};
    };

    auto getValue = [&](int idx) {
      T val = cellIt.getCorner(idx).getValue();
      return inverted ? -val : val;
    };

    auto getInterp = [&](int p_a, int p_b) {
      const T v_a = getValue(p_a);
      const T v_b = getValue(p_b);
      if (std::abs(v_a - v_b) < epsilon) return T(0.5);
      return v_a / (v_a - v_b);
    };

    // Determine axes
    int axisA = 0; while ((C ^ A) != (1 << axisA)) axisA++;
    int axisB = 0; while ((C ^ B) != (1 << axisB)) axisB++;
    int axisZ = 3 - axisA - axisB;

    // Solve 2D face problem to find sharp point on a face
    auto solveFace = [&](int axis, int corner, int n1, int n2) -> std::optional<Vec3D<T>> {
        viennahrle::Index<D> faceIdx = cellIt.getIndices();
        if ((corner >> axis) & 1) faceIdx[axis]++;
        
        if (faceNodes[axis].find(faceIdx) != faceNodes[axis].end()) {
            unsigned nodeId = faceNodes[axis][faceIdx];
            Vec3D<T> pos = mesh->nodes[nodeId];
            for(int d=0; d<3; ++d) pos[d] = pos[d]/levelSet->getGrid().getGridDelta() - cellIt.getIndices(d);
            return pos;
        }

        Vec3D<T> N1 = getNormal(n1);
        Vec3D<T> N2 = getNormal(n2);
        
        if (std::abs(DotProduct(N1, N2)) >= 0.1) return std::nullopt;

        int u = (axis + 1) % 3;
        int v = (axis + 2) % 3;
        
        int axis1 = 0; while((corner^n1) != (1<<axis1)) axis1++;
        T t1 = getInterp(corner, n1);
        
        int axis2 = 0; while((corner^n2) != (1<<axis2)) axis2++;
        T t2 = getInterp(corner, n2);
        
        Vec3D<T> P1; P1[0] = ((corner>>0)&1); P1[1] = ((corner>>1)&1); P1[2] = ((corner>>2)&1);
        P1[axis1] = ((corner>>axis1)&1) ? (1.0-t1) : t1;
        
        Vec3D<T> P2; P2[0] = ((corner>>0)&1); P2[1] = ((corner>>1)&1); P2[2] = ((corner>>2)&1);
        P2[axis2] = ((corner>>axis2)&1) ? (1.0-t2) : t2;
        
        double det = N1[u]*N2[v] - N1[v]*N2[u];
        if (std::abs(det) < 1e-6) return std::nullopt;
        
        double c1 = N1[u]*P1[u] + N1[v]*P1[v];
        double c2 = N2[u]*P2[u] + N2[v]*P2[v];
        
        Vec3D<T> P;
        P[axis] = P1[axis];
        P[u] = (c1*N2[v] - c2*N1[v]) / det;
        P[v] = (c2*N1[u] - c1*N2[u]) / det;
        
        auto distSq = [](const Vec3D<T> &a, const Vec3D<T> &b) {
          return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                 (a[2] - b[2]) * (a[2] - b[2]);
        };

        if (distSq(P, P1) > 9.0 || distSq(P, P2) > 9.0) return std::nullopt;

        Vec3D<T> globalP = P;
        for(int d=0; d<3; ++d) globalP[d] = (globalP[d] + cellIt.getIndices(d)) * levelSet->getGrid().getGridDelta();
        unsigned nodeId = mesh->insertNextNode(globalP);
        faceNodes[axis][faceIdx] = nodeId;
        if (updatePointData) newDataSourceIds[0].push_back(cellIt.getCorner(corner).getPointId());

        stitchToNeighbor(cellIt, axis, (corner >> axis) & 1, nodeId, nodes, valueIt);
        
        return P;
    };

    int D_corner = A ^ (1 << axisB);
    int A_z = A ^ (1 << axisZ);
    int B_z = B ^ (1 << axisZ);
    int C_z = C ^ (1 << axisZ);

    auto P_A_opt = solveFace(axisA, A, D_corner, A_z);
    auto P_B_opt = solveFace(axisB, B, D_corner, B_z);
    
    if (!P_A_opt || !P_B_opt) return false;
    
    Vec3D<T> P_A = *P_A_opt;
    Vec3D<T> P_B = *P_B_opt;
    
    // Construct S
    Vec3D<T> S;
    S[axisA] = P_B[axisA];
    S[axisB] = P_A[axisB];
    S[axisZ] = (P_A[axisZ] + P_B[axisZ]) * 0.5;

    unsigned nS = mesh->insertNextNode([&](){
        Vec3D<T> gS = S;
        for(int i=0; i<3; ++i) gS[i] = (gS[i] + cellIt.getIndices(i)) * levelSet->getGrid().getGridDelta();
        return gS;
    }());
    if (updatePointData) newDataSourceIds[0].push_back(cellIt.getCorner(C).getPointId());

    // Calculate S_Face
    Vec3D<T> S_Face = S;
    S_Face[axisZ] = static_cast<T>((C >> axisZ) & 1);

    unsigned nS_Face = mesh->insertNextNode([&](){
        Vec3D<T> gS = S_Face;
        for(int i=0; i<3; ++i) gS[i] = (gS[i] + cellIt.getIndices(i)) * levelSet->getGrid().getGridDelta();
        return gS;
    }());
    if (updatePointData) newDataSourceIds[0].push_back(cellIt.getCorner(C).getPointId());

    stitchToNeighbor(cellIt, axisZ, (C >> axisZ) & 1, nS_Face, nodes, valueIt);

    // Get face nodes
    viennahrle::Index<D> faceIdxA = cellIt.getIndices(); if ((A >> axisA) & 1) faceIdxA[axisA]++;
    unsigned nNA = faceNodes[axisA][faceIdxA];
    
    viennahrle::Index<D> faceIdxB = cellIt.getIndices(); if ((B >> axisB) & 1) faceIdxB[axisB]++;
    unsigned nNB = faceNodes[axisB][faceIdxB];

    // Get boundary intersection nodes
    auto getEdgeNode = [&](int v1, int v2) {
        int edgeIdx = -1;
        for (int e = 0; e < 12; ++e) {
            if ((corner0[e] == static_cast<unsigned>(v1) && corner1[e] == static_cast<unsigned>(v2)) ||
                (corner0[e] == static_cast<unsigned>(v2) && corner1[e] == static_cast<unsigned>(v1))) {
                edgeIdx = e;
                break;
            }
        }
        return this->getNode(cellIt, edgeIdx, nodes, newDataSourceIds);
    };

    unsigned nI_AD = getEdgeNode(A, D_corner);
    unsigned nI_AZ = getEdgeNode(A, A_z);
    unsigned nI_BD = getEdgeNode(B, D_corner);
    unsigned nI_BZ = getEdgeNode(B, B_z);
    unsigned nI_CZ = getEdgeNode(C, C_z);

    // Winding order
    int parity = (C&1) + ((C>>1)&1) + ((C>>2)&1);
    bool flip = (parity % 2 != 0) ^ inverted;
    
    int vA = 0; while((A^C) != (1<<vA)) vA++;
    int vB = 0; while((B^C) != (1<<vB)) vB++;
    bool is_cyclic = ((vA+1)%3 == vB);
    if (!is_cyclic) {
        std::swap(nNA, nNB);
        std::swap(nI_AD, nI_BD);
        std::swap(nI_AZ, nI_BZ);
    }

    if (!flip) {
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nNA, nI_AD});
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nI_AD, nS_Face});
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nS_Face, nI_BD});
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nI_BD, nNB});
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nNB, nI_BZ});
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nI_BZ, nI_CZ});
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nI_CZ, nI_AZ});
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nI_AZ, nNA});
    } else {
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nI_AD, nNA});
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nS_Face, nI_AD});
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nI_BD, nS_Face});
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nNB, nI_BD});
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nI_BZ, nNB});
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nI_CZ, nI_BZ});
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nI_AZ, nI_CZ});
        mesh->insertNextElement(std::array<unsigned, 3>{nS, nNA, nI_AZ});
    }

    return true;
  }

  bool generateCanonicalSharpEdge3D(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt,
      int transform, int axis, bool inverted,
      std::map<viennahrle::Index<D>, unsigned> *nodes,
      std::map<viennahrle::Index<D>, unsigned> *faceNodes,
      std::vector<std::vector<unsigned>> &newDataSourceIds,
      viennahrle::ConstSparseIterator<hrleDomainType> &valueIt) {

    const auto *normalVectorData = levelSet->getPointData().getVectorData(
        viennals::CalculateNormalVectors<T, D>::normalVectorsLabel);

    if (normalVectorData == nullptr)
      return false;

    // Helper to map global vector to canonical frame
    auto toCanonical = [&](Vec3D<T> v) {
      Vec3D<T> res;
      // 1. Apply reflection (transform)
      if (transform & 1) v[0] = -v[0];
      if (transform & 2) v[1] = -v[1];
      if (transform & 4) v[2] = -v[2];
      // 2. Apply rotation (axis permutation)
      // axis=0: XYZ -> XYZ (0,1,2)
      // axis=1: XYZ -> YZX (1,2,0) - puts Y in X place
      // axis=2: XYZ -> ZXY (2,0,1) - puts Z in X place
      if (axis == 0) res = v;
      else if (axis == 1) res = Vec3D<T>{v[1], v[2], v[0]};
      else res = Vec3D<T>{v[2], v[0], v[1]};
      return res;
    };

    // Helper to map canonical point back to global frame
    auto fromCanonical = [&](Vec3D<T> v) {
      Vec3D<T> res;
      // 1. Inverse rotation
      if (axis == 0) res = v;
      else if (axis == 1) res = Vec3D<T>{v[2], v[0], v[1]};
      else res = Vec3D<T>{v[1], v[2], v[0]};
      // 2. Inverse reflection (same as forward for 0/1 swap)
      if (transform & 1) res[0] = T(1.0) - res[0];
      if (transform & 2) res[1] = T(1.0) - res[1];
      if (transform & 4) res[2] = T(1.0) - res[2];
      return res;
    };

    // Helper to get value/normal with inversion handling
    auto getValue = [&](int idx) {
      T val = cellIt.getCorner(idx).getValue();
      return inverted ? -val : val;
    };

    auto getNormal = [&](int idx) {
      auto corner = cellIt.getCorner(idx);
      if (corner.isDefined()) {
        Vec3D<T> n = (*normalVectorData)[corner.getPointId()];
        if (inverted) {
          for (int i = 0; i < 3; ++i)
            n[i] = -n[i];
        }
        return n;
      }
      return Vec3D<T>{};
    };

    // Calculate face nodes P0 (at x=0) and P1 (at x=1) in canonical frame
    Vec3D<T> P[2];
    unsigned faceNodeIds[2];

    for (int k = 0; k < 2; ++k) {
      // Check if face node already exists
      // k=0 -> x=0 face (left). In global, this corresponds to the face perpendicular to 'axis'
      // at the 'transform' side.
      // If transform has bit 'axis' set, then x=0 in canonical is x=1 in global (relative to cell origin).
      // We need to be careful with the face index key.
      // Let's compute the global index of the face.
      viennahrle::Index<D> faceIdx = cellIt.getIndices();
      // If k=0 (canonical x=0), and transform has bit 'axis' set (flipped), this is global x=1 face.
      // If k=1 (canonical x=1), and transform has bit 'axis' set, this is global x=0 face.
      bool isHighFace = (k == 1) ^ ((transform >> axis) & 1);
      if (isHighFace) faceIdx[axis]++;

      auto it = faceNodes[axis].find(faceIdx);
      if (it != faceNodes[axis].end()) {
        faceNodeIds[k] = it->second;
        Vec3D<T> relPos = mesh->nodes[faceNodeIds[k]];
        T gridDelta = levelSet->getGrid().getGridDelta();
        for (int d = 0; d < 3; ++d) {
          relPos[d] /= gridDelta;
          relPos[d] -= static_cast<T>(cellIt.getIndices(d));
        }
        P[k] = toCanonical(relPos);
      } else {
        // Calculate P[k]
        // For k=0 (x=0): Active corner is 0 (000). Neighbors 2 (010) and 4 (001).
        // For k=1 (x=1): Active corner is 1 (100). Neighbors 3 (110) and 5 (101).
        // We need to map these canonical indices back to global to get normals/values.
        
        // Inverse map indices
        auto mapIdx = [&](int c) {
          // c is canonical index.
          // 1. Inverse rotation:
          int r = 0;
          if (axis == 0) r = c;
          else if (axis == 1) r = ((c&1)<<2) | ((c&2)>>1) | ((c&4)>>1); // x->z, y->x, z->y? No.
          // axis=1: x_can = y_glob. y_can = z_glob. z_can = x_glob.
          // bit 0 (x_can) -> bit 1 (y_glob). bit 1 (y_can) -> bit 2 (z_glob). bit 2 (z_can) -> bit 0 (x_glob).
          else r = ((c&1)<<1) | ((c&2)<<1) | ((c&4)>>2); // axis=2: x_can=z_glob...
          
          // Actually, simpler: construct vector, permute, reconstruct index
          int x = (c&1), y = (c>>1)&1, z = (c>>2)&1;
          int gx, gy, gz;
          if (axis == 0) { gx=x; gy=y; gz=z; }
          else if (axis == 1) { gx=z; gy=x; gz=y; } // x_can comes from y_glob? No, toCanonical: res[0]=v[1]. So x_can = y_glob.
          else { gx=y; gy=z; gz=x; }
          
          // 2. Inverse reflection
          if (transform & 1) gx = 1-gx;
          if (transform & 2) gy = 1-gy;
          if (transform & 4) gz = 1-gz;
          
          return gx | (gy<<1) | (gz<<2);
        };

        int c0 = k; // 0 or 1
        int c_y = c0 | 2; // Neighbor in Y (canonical)
        int c_z = c0 | 4; // Neighbor in Z (canonical)

        Vec3D<T> n_y = toCanonical(getNormal(mapIdx(c_y)));
        Vec3D<T> n_z = toCanonical(getNormal(mapIdx(c_z)));

        if (std::abs(DotProduct(n_y, n_z)) >= 0.1) {
          return false;
        }

        // Solve 2D problem in YZ plane
        // Line 1: passes through intersection on edge c0-c_y. Normal n_y (projected).
        // Line 2: passes through intersection on edge c0-c_z. Normal n_z (projected).
        
        auto getInterp = [&](int p_a, int p_b) {
           T v_a = getValue(mapIdx(p_a));
           T v_b = getValue(mapIdx(p_b));
           return v_a / (v_a - v_b);
        };

        T t_y = getInterp(c0, c_y); // Fraction along Y axis
        T t_z = getInterp(c0, c_z); // Fraction along Z axis
        
        // In canonical local frame relative to c0:
        // Point on Y-edge: (0, t_y, 0)
        // Point on Z-edge: (0, 0, t_z)
        // Normals n_y and n_z.
        // Solve for (y, z):
        // n_y.y * (y - t_y) + n_y.z * (z - 0) = 0
        // n_z.y * (y - 0) + n_z.z * (z - t_z) = 0
        
        double det = n_y[1] * n_z[2] - n_y[2] * n_z[1];
        if (std::abs(det) < 1e-6) return false;
        
        double d1 = n_y[1] * t_y;
        double d2 = n_z[2] * t_z;
        
        double y = (d1 * n_z[2] - d2 * n_y[2]) / det;
        double z = (d2 * n_y[1] - d1 * n_z[1]) / det;
        
        P[k] = Vec3D<T>{(T)k, (T)y, (T)z};
        
        // Check if the point is too far away from the intersection points
        // 2 grid deltas = 2.0 in canonical coordinates
        auto distSq = [](const Vec3D<T> &a, const Vec3D<T> &b) {
          return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                 (a[2] - b[2]) * (a[2] - b[2]);
        };

        if (distSq(P[k], Vec3D<T>{(T)k, t_y, 0}) > 9.0 ||
            distSq(P[k], Vec3D<T>{(T)k, 0, t_z}) > 9.0) {
          return false;
        }

        // Store
        Vec3D<T> globalP = fromCanonical(P[k]);
        for(int d=0; d<3; ++d) globalP[d] = (globalP[d] + cellIt.getIndices(d)) * levelSet->getGrid().getGridDelta();
        
        faceNodeIds[k] = mesh->insertNextNode(globalP);
        faceNodes[axis][faceIdx] = faceNodeIds[k];
        if (updatePointData) newDataSourceIds[0].push_back(cellIt.getCorner(mapIdx(c0)).getPointId());

        stitchToNeighbor(cellIt, axis, isHighFace, faceNodeIds[k], nodes, valueIt);
      }
    }

    // Calculate edge intersection points in canonical frame
    // E02: on edge 0-2 (Y axis). (0, t_02, 0)
    // E04: on edge 0-4 (Z axis). (0, 0, t_04)
    // E13: on edge 1-3 (Y axis). (1, t_13, 0)
    // E15: on edge 1-5 (Z axis). (1, 0, t_15)
    
    // We need to get the node IDs for these. We can use the existing getNode helper, 
    // but we need to map canonical edges to global edges.
    auto getEdgeNode = [&](int c1, int c2) {
       // Map canonical corners to global
       auto mapIdx = [&](int c) {
          int x = (c&1), y = (c>>1)&1, z = (c>>2)&1;
          int gx, gy, gz;
          if (axis == 0) { gx=x; gy=y; gz=z; }
          else if (axis == 1) { gx=z; gy=x; gz=y; }
          else { gx=y; gy=z; gz=x; }
          if (transform & 1) gx = 1-gx;
          if (transform & 2) gy = 1-gy;
          if (transform & 4) gz = 1-gz;
          return gx | (gy<<1) | (gz<<2);
       };
       int g1 = mapIdx(c1);
       int g2 = mapIdx(c2);
       
       // Find edge index
       int edgeIdx = -1;
       for (int e = 0; e < 12; ++e) {
        if ((corner0[e] == static_cast<unsigned>(g1) && corner1[e] == static_cast<unsigned>(g2)) ||
            (corner0[e] == static_cast<unsigned>(g2) && corner1[e] == static_cast<unsigned>(g1))) {
          edgeIdx = e;
          break;
        }
       }
       return this->getNode(cellIt, edgeIdx, nodes, newDataSourceIds);
    };

    unsigned n02 = getEdgeNode(0, 2);
    unsigned n04 = getEdgeNode(0, 4);
    unsigned n13 = getEdgeNode(1, 3);
    unsigned n15 = getEdgeNode(1, 5);

    // Generate Quads (split into triangles)
    // Quad 1 (Y-interface): E02, P0, P1, E13. Normal +Y.
    // Winding for +Y normal: E02 -> P0 -> P1 -> E13.
    mesh->insertNextElement(std::array<unsigned, 3>{n02, faceNodeIds[0], faceNodeIds[1]});
    mesh->insertNextElement(std::array<unsigned, 3>{n02, faceNodeIds[1], n13});

    // Quad 2 (Z-interface): E04, E15, P1, P0. Normal +Z.
    // Winding for +Z normal: E04 -> E15 -> P1 -> P0.
    mesh->insertNextElement(std::array<unsigned, 3>{n04, n15, faceNodeIds[1]});
    mesh->insertNextElement(std::array<unsigned, 3>{n04, faceNodeIds[1], faceNodeIds[0]});

    return true;
  }

  bool generateSharpEdge3D(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt,
      std::map<viennahrle::Index<D>, unsigned> *nodes,
      std::map<viennahrle::Index<D>, unsigned> *faceNodes,
      std::vector<std::vector<unsigned>> &newDataSourceIds,
      viennahrle::ConstSparseIterator<hrleDomainType> &valueIt) {

    int countNeg = 0;
    int countPos = 0;
    int negMask = 0;
    int posMask = 0;

    for (int i = 0; i < 8; ++i) {
      T val = cellIt.getCorner(i).getValue();
      if (val < -epsilon) {
        countNeg++;
        negMask |= (1 << i);
      } else if (val > epsilon) {
        countPos++;
        posMask |= (1 << i);
      }
    }

    bool inverted = false;
    int mask = 0;
    if (countNeg == 2) {
      mask = negMask;
      inverted = false;
    } else if (countPos == 2) {
      mask = posMask;
      inverted = true;
    } else {
      return false;
    }

    // Check if the two vertices share an edge
    int v1 = -1, v2 = -1;
    for (int i = 0; i < 8; ++i) {
      if ((mask >> i) & 1) {
        if (v1 == -1) v1 = i;
        else v2 = i;
      }
    }

    int diff = v1 ^ v2;
    if ((diff & (diff - 1)) != 0) return false; // Not connected by a single edge

    int axis = 0;
    while ((diff >> (axis + 1)) > 0) axis++;

    // Transform maps v1 to 0.
    int transform = v1;

    return generateCanonicalSharpEdge3D(cellIt, transform, axis, inverted, nodes, faceNodes, newDataSourceIds, valueIt);
  }

  bool generateCanonicalSharpCorner3D(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt,
      int transform, bool inverted,
      std::map<viennahrle::Index<D>, unsigned> *nodes,
      std::map<viennahrle::Index<D>, unsigned> *faceNodes,
      std::vector<std::vector<unsigned>> &newDataSourceIds,
      viennahrle::ConstSparseIterator<hrleDomainType> &valueIt,
      Vec3D<T> &cornerPos) {

    const auto *normalVectorData = levelSet->getPointData().getVectorData(
        viennals::CalculateNormalVectors<T, D>::normalVectorsLabel);

    if (normalVectorData == nullptr)
      return false;

    // Helper to map global vector to canonical frame (just reflection)
    auto toCanonical = [&](Vec3D<T> v) {
      if (transform & 1)
        v[0] = -v[0];
      if (transform & 2)
        v[1] = -v[1];
      if (transform & 4)
        v[2] = -v[2];
      return v;
    };

    // Helper to map canonical point back to global frame
    auto fromCanonical = [&](Vec3D<T> v) {
      if (transform & 1)
        v[0] = T(1.0) - v[0];
      if (transform & 2)
        v[1] = T(1.0) - v[1];
      if (transform & 4)
        v[2] = T(1.0) - v[2];
      return v;
    };

    auto getValue = [&](int idx) {
      T val = cellIt.getCorner(idx).getValue();
      return inverted ? -val : val;
    };

    auto getNormal = [&](int idx) {
      auto corner = cellIt.getCorner(idx);
      if (corner.isDefined()) {
        Vec3D<T> n = (*normalVectorData)[corner.getPointId()];
        if (inverted) {
          for (int i = 0; i < 3; ++i)
            n[i] = -n[i];
        }
        return n;
      }
      return Vec3D<T>{};
    };

    auto getInterp = [&](int p_a, int p_b) {
      const T v_a = getValue(p_a);
      const T v_b = getValue(p_b);
      if (std::abs(v_a - v_b) < epsilon)
        return T(0.5);
      return v_a / (v_a - v_b);
    };

    Vec3D<T> norm1, norm2, norm3;
    double d1, d2, d3;
    Vec3D<T> P1, P2, P4;

    const int g1 = transform ^ 1;
    const int g2 = transform ^ 2;
    const int g4 = transform ^ 4;

    norm1 = toCanonical(getNormal(g1));
    norm2 = toCanonical(getNormal(g2));
    norm3 = toCanonical(getNormal(g4));

    const T t1 = getInterp(transform, g1);
    const T t2 = getInterp(transform, g2);
    const T t4 = getInterp(transform, g4);

    P1 = {t1, 0, 0};
    P2 = {0, t2, 0};
    P4 = {0, 0, t4};

    d1 = DotProduct(norm1, P1);
    d2 = DotProduct(norm2, P2);
    d3 = DotProduct(norm3, P4);

    if (std::abs(DotProduct(norm1, norm2)) >= 0.1 ||
        std::abs(DotProduct(norm1, norm3)) >= 0.1 ||
        std::abs(DotProduct(norm2, norm3)) >= 0.1) {
      return false;
    }

    // We solve the system of 3 plane equations:
    // norm1 . S = d1
    // norm2 . S = d2
    // norm3 . S = d3
    // This is a linear system Ax = b, where A is the matrix of normals, x is
    // the corner position S, and b is the vector of dot products.

    // Using Cramer's rule to solve the system

    const Vec3D<T> c23 = CrossProduct(norm2, norm3);
    const Vec3D<T> c31 = CrossProduct(norm3, norm1);
    const Vec3D<T> c12 = CrossProduct(norm1, norm2);

    const double det = DotProduct(norm1, c23);

    if (std::abs(det) < epsilon)
      return false; // Planes are parallel or linearly dependent

    const T invDet = 1.0 / det;
    Vec3D<T> S;
    for (int i = 0; i < 3; ++i) {
      S[i] = (d1 * c23[i] + d2 * c31[i] + d3 * c12[i]) * invDet;
    }

    auto distSq = [](const Vec3D<T> &a, const Vec3D<T> &b) {
      return (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
             (a[2] - b[2]) * (a[2] - b[2]);
    };

    if (distSq(S, P1) > 9.0 || distSq(S, P2) > 9.0 || distSq(S, P4) > 9.0) {
      return false;
    }

    // Transform corner position back from canonical to global coordinates
    Vec3D<T> globalS = fromCanonical(S);
    for (int d = 0; d < 3; ++d) {
      cornerPos[d] = (globalS[d] + cellIt.getIndices(d)) *
                     levelSet->getGrid().getGridDelta();
    }

    return true;
  }

  bool generateSharpCorner3D(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt,
      std::map<viennahrle::Index<D>, unsigned> *nodes,
      std::map<viennahrle::Index<D>, unsigned> &cornerNodes,
      std::map<viennahrle::Index<D>, unsigned> *faceNodes,
      std::vector<std::vector<unsigned>> &newDataSourceIds,
      viennahrle::ConstSparseIterator<hrleDomainType> &valueIt) {

#if defined(__GNUC__) || defined(__clang__)
#define POPCOUNT __builtin_popcount
#else
    auto popcount = [](unsigned i) {
      i = i - ((i >> 1) & 0x55555555);
      i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
      return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
    };
#define POPCOUNT popcount
#endif

    int countNeg = 0;
    int countPos = 0;
    int negMask = 0;
    int posMask = 0;

    for (int i = 0; i < 8; ++i) {
      T val = cellIt.getCorner(i).getValue();
      if (val < -epsilon) {
        countNeg++;
        negMask |= (1 << i);
      } else if (val > epsilon) {
        countPos++;
        posMask |= (1 << i);
      }
    }

    bool inverted = false;
    int transform = -1;
    // bool is3vs5 = false;

    if (countNeg == 1) {
        transform = 0;
        while (!((negMask >> transform) & 1)) transform++;
        inverted = false;
    } else if (countPos == 1) {
        transform = 0;
        while (!((posMask >> transform) & 1)) transform++;
        inverted = true;
    }

#undef POPCOUNT
    
    if (transform == -1)
      return false;

    const auto *normalVectorData = levelSet->getPointData().getVectorData(
        viennals::CalculateNormalVectors<T, D>::normalVectorsLabel);

    Vec3D<T> cornerPos;
    if (generateCanonicalSharpCorner3D(cellIt, transform, inverted, nodes, faceNodes, newDataSourceIds, valueIt, cornerPos)) {
        unsigned nCorner;
        viennahrle::Index<D> cornerIdx = cellIt.getIndices();
        cornerIdx += viennahrle::BitMaskToIndex<D>(transform);

        auto it = cornerNodes.find(cornerIdx);
        if (it != cornerNodes.end()) {
          nCorner = it->second;
        } else {
          nCorner = mesh->insertNextNode(cornerPos);
          cornerNodes[cornerIdx] = nCorner;
          if (updatePointData)
            newDataSourceIds[0].push_back(cellIt.getCorner(transform).getPointId());
        }

        // Get edge nodes
        auto getEdgeNode = [&](int neighbor) {
            int v1 = transform;
            int v2 = neighbor;
            int edgeIdx = -1;
            for (int e = 0; e < 12; ++e) {
                if ((corner0[e] == static_cast<unsigned>(v1) && corner1[e] == static_cast<unsigned>(v2)) ||
                    (corner0[e] == static_cast<unsigned>(v2) && corner1[e] == static_cast<unsigned>(v1))) {
                    edgeIdx = e;
                    break;
                }
            }
            return this->getNode(cellIt, edgeIdx, nodes, newDataSourceIds);
        };

        auto getEdgeNode2 = [&](int v1, int v2) {
            int edgeIdx = -1;
            for (int e = 0; e < 12; ++e) {
                if ((corner0[e] == static_cast<unsigned>(v1) && corner1[e] == static_cast<unsigned>(v2)) ||
                    (corner0[e] == static_cast<unsigned>(v2) && corner1[e] == static_cast<unsigned>(v1))) {
                    edgeIdx = e;
                    break;
                }
            }
            return this->getNode(cellIt, edgeIdx, nodes, newDataSourceIds);
        };

        unsigned n1 = getEdgeNode(transform ^ 1);
        unsigned n2 = getEdgeNode(transform ^ 2);
        unsigned n4 = getEdgeNode(transform ^ 4);

        const auto *normalVectorData = levelSet->getPointData().getVectorData(
            viennals::CalculateNormalVectors<T, D>::normalVectorsLabel);

        // Try to generate cube corner topology
        if (normalVectorData) {
            // Re-calculate canonical normals and plane constants
            auto toCanonical = [&](Vec3D<T> v) {
                if (transform & 1) v[0] = -v[0];
                if (transform & 2) v[1] = -v[1];
                if (transform & 4) v[2] = -v[2];
                return v;
            };
            
            auto fromCanonical = [&](Vec3D<T> v) {
                if (transform & 1) v[0] = T(1.0) - v[0];
                if (transform & 2) v[1] = T(1.0) - v[1];
                if (transform & 4) v[2] = T(1.0) - v[2];
                return v;
            };

            auto getNormal = [&](int idx) {
                auto corner = cellIt.getCorner(idx);
                if (corner.isDefined()) {
                    Vec3D<T> n = (*normalVectorData)[corner.getPointId()];
                    if (inverted) for (int i = 0; i < 3; ++i) n[i] = -n[i];
                    return n;
                }
                return Vec3D<T>{};
            };
            
            auto getValue = [&](int idx) {
                T val = cellIt.getCorner(idx).getValue();
                return inverted ? -val : val;
            };

            auto getInterp = [&](int p_a, int p_b) {
                const T v_a = getValue(p_a);
                const T v_b = getValue(p_b);
                if (std::abs(v_a - v_b) < epsilon) return T(0.5);
                return v_a / (v_a - v_b);
            };

            Vec3D<T> norm1, norm2, norm3;
            double d1, d2, d3;

            // Identify neighbors
            int n1_idx = -1, n2_idx = -1, n3_idx = -1;
            int neighbors[] = {transform ^ 1, transform ^ 2, transform ^ 4};
            
            for (int n_idx : neighbors) {
                T val = cellIt.getCorner(n_idx).getValue();
                // If signs are opposite (product < 0), this is a face neighbor
                if (val * cellIt.getCorner(transform).getValue() < 0) {
                    if (n1_idx == -1) n1_idx = n_idx;
                    else n2_idx = n_idx;
                } else {
                    n3_idx = n_idx;
                }
            }
            
            // Standard 1-vs-7 case
            norm1 = toCanonical(getNormal(transform ^ 1));
            norm2 = toCanonical(getNormal(transform ^ 2));
            norm3 = toCanonical(getNormal(transform ^ 4));
            
            d1 = norm1[0] * getInterp(transform, transform ^ 1);
            d2 = norm2[1] * getInterp(transform, transform ^ 2);
            d3 = norm3[2] * getInterp(transform, transform ^ 4);

            auto solve2x2 = [&](double a1, double b1, double c1,
                                double a2, double b2, double c2,
                                T &res1, T &res2) {
                double det = a1 * b2 - a2 * b1;
                if (std::abs(det) < 1e-6) return false;
                res1 = (c1 * b2 - c2 * b1) / det;
                res2 = (a1 * c2 - a2 * c1) / det;
                return true;
            };

            Vec3D<T> Fx, Fy, Fz;
            Fx[0] = 0.0; Fy[1] = 0.0; Fz[2] = 0.0;
            
            bool valid = true;
            valid &= solve2x2(norm2[1], norm2[2], d2 - norm2[0] * Fx[0], norm3[1], norm3[2], d3 - norm3[0] * Fx[0], Fx[1], Fx[2]);
            valid &= solve2x2(norm1[0], norm1[2], d1 - norm1[1] * Fy[1], norm3[0], norm3[2], d3 - norm3[1] * Fy[1], Fy[0], Fy[2]);
            valid &= solve2x2(norm1[0], norm1[1], d1 - norm1[2] * Fz[2], norm2[0], norm2[1], d2 - norm2[2] * Fz[2], Fz[0], Fz[1]);

            auto checkBounds = [&](Vec3D<T> p) {
                return p[0] >= -1e-4 && p[0] <= 1.0+1e-4 &&
                       p[1] >= -1e-4 && p[1] <= 1.0+1e-4 &&
                       p[2] >= -1e-4 && p[2] <= 1.0+1e-4;
            };
            
            if (valid && checkBounds(Fx) && checkBounds(Fy) && checkBounds(Fz)) {
                auto getFaceNode = [&](int axis, Vec3D<T> pos) {
                    viennahrle::Index<D> faceIdx = cellIt.getIndices();
                    if ((transform >> axis) & 1) faceIdx[axis]++;
                    
                    if (faceNodes[axis].find(faceIdx) != faceNodes[axis].end()) {
                        return faceNodes[axis][faceIdx];
                    }
                    Vec3D<T> globalPos = fromCanonical(pos);
                    for(int d=0; d<3; ++d) globalPos[d] = (globalPos[d] + cellIt.getIndices(d)) * levelSet->getGrid().getGridDelta();
                    unsigned nodeId = mesh->insertNextNode(globalPos);
                    faceNodes[axis][faceIdx] = nodeId;
                    if (updatePointData) newDataSourceIds[0].push_back(cellIt.getCorner(transform).getPointId());
                    stitchToNeighbor(cellIt, axis, (transform >> axis) & 1, nodeId, nodes, valueIt);
                    return nodeId;
                };

                unsigned nFx = getFaceNode(0, Fx);
                unsigned nFy = getFaceNode(1, Fy);
                unsigned nFz = getFaceNode(2, Fz);

                // Calculate parity of the transform (number of reflections)
                int parity = (transform & 1) + ((transform >> 1) & 1) + ((transform >> 2) & 1);

                // Flip winding if parity is odd XOR inverted is true.
                bool flip = (parity % 2 != 0) ^ inverted;

                if (!flip) {
                    mesh->insertNextElement(std::array<unsigned, 3>{nCorner, nFy, n1});
                    mesh->insertNextElement(std::array<unsigned, 3>{nCorner, n1, nFz});
                    mesh->insertNextElement(std::array<unsigned, 3>{nCorner, nFz, n2});
                    mesh->insertNextElement(std::array<unsigned, 3>{nCorner, n2, nFx});
                    mesh->insertNextElement(std::array<unsigned, 3>{nCorner, nFx, n4});
                    mesh->insertNextElement(std::array<unsigned, 3>{nCorner, n4, nFy});
                } else {
                    mesh->insertNextElement(std::array<unsigned, 3>{nCorner, n1, nFy});
                    mesh->insertNextElement(std::array<unsigned, 3>{nCorner, nFz, n1});
                    mesh->insertNextElement(std::array<unsigned, 3>{nCorner, n2, nFz});
                    mesh->insertNextElement(std::array<unsigned, 3>{nCorner, nFx, n2});
                    mesh->insertNextElement(std::array<unsigned, 3>{nCorner, n4, nFx});
                    mesh->insertNextElement(std::array<unsigned, 3>{nCorner, nFy, n4});
                }
                return true;
            }
        }

        // Triangles
        // Cycle around normal (1,1,1) in canonical is 1->2->4->1 (x->y->z->x)
        // Triangles: (n1, S, n4), (n4, S, n2), (n2, S, n1)
        // If inverted, flip winding.
        
        // Calculate parity of the transform (number of reflections)
        int parity = (transform & 1) + ((transform >> 1) & 1) + ((transform >> 2) & 1);

        // Flip winding if parity is odd XOR inverted is true.
        bool flip = (parity % 2 != 0) ^ inverted;

        if (!flip) {
            mesh->insertNextElement(std::array<unsigned, 3>{n1, nCorner, n4});
            mesh->insertNextElement(std::array<unsigned, 3>{n4, nCorner, n2});
            mesh->insertNextElement(std::array<unsigned, 3>{n2, nCorner, n1});
        } else {
            mesh->insertNextElement(std::array<unsigned, 3>{n1, n4, nCorner});
            mesh->insertNextElement(std::array<unsigned, 3>{n4, n2, nCorner});
            mesh->insertNextElement(std::array<unsigned, 3>{n2, n1, nCorner});
        }
        
        return true;
    }

    return false;
  }

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
