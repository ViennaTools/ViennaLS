#pragma once

#include <lsPreCompileMacros.hpp>

#include <map>
#include <unordered_map>
#include <unordered_set>

#include <hrleSparseCellIterator.hpp>
#include <hrleSparseStarIterator.hpp>
#include <lsCalculateNormalVectors.hpp>
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
protected:
  using lsDomainType = viennals::Domain<T, D>;
  using hrleDomainType = typename lsDomainType::DomainType;
  using hrleIndex = viennahrle::Index<D>;
  using ConstSparseIterator = viennahrle::ConstSparseIterator<hrleDomainType>;

  std::vector<SmartPointer<lsDomainType>> levelSets;
  SmartPointer<Mesh<T>> mesh = nullptr;

  const T epsilon;
  const double minNodeDistanceFactor;
  bool updatePointData = true;
  bool generateSharpCorners = false;

  struct I3 {
    int x, y, z;
    bool operator==(const I3 &o) const {
      return x == o.x && y == o.y && z == o.z;
    }
  };

  struct I3Hash {
    size_t operator()(const I3 &k) const {
      // 64-bit mix
      uint64_t a = (uint64_t)(uint32_t)k.x;
      uint64_t b = (uint64_t)(uint32_t)k.y;
      uint64_t c = (uint64_t)(uint32_t)k.z;
      uint64_t h = a * 0x9E3779B185EBCA87ULL;
      h ^= b + 0xC2B2AE3D27D4EB4FULL + (h << 6) + (h >> 2);
      h ^= c + 0x165667B19E3779F9ULL + (h << 6) + (h >> 2);
      return (size_t)h;
    }
  };

  // Context for mesh generation to share between helpers
  std::unordered_map<I3, unsigned, I3Hash> nodeIdByBin;
  std::unordered_set<I3, I3Hash> uniqueElements;
  std::vector<Vec3D<T>> currentNormals;
  std::vector<T> currentMaterials;
  double currentGridDelta;
  T currentMaterialId = 0;
  SmartPointer<lsDomainType> currentLevelSet = nullptr;
  typename PointData<T>::VectorDataType *normalVectorData = nullptr;

  // Store sharp corner nodes created during this cell iteration
  std::vector<std::pair<unsigned, Vec3D<T>>> matSharpCornerNodes;

  static constexpr unsigned int corner0[12] = {0, 1, 2, 0, 4, 5,
                                               6, 4, 0, 1, 3, 2};
  static constexpr unsigned int corner1[12] = {1, 3, 3, 2, 5, 7,
                                               7, 6, 4, 5, 7, 6};
  static constexpr unsigned int direction[12] = {0, 1, 0, 1, 0, 1,
                                                 0, 1, 2, 2, 2, 2};

public:
  explicit ToSurfaceMesh(double mnd = 0.05, double eps = 1e-12)
      : minNodeDistanceFactor(mnd), epsilon(eps) {}

  ToSurfaceMesh(const SmartPointer<lsDomainType> passedLevelSet,
                SmartPointer<Mesh<T>> passedMesh, double mnd = 0.05,
                double eps = 1e-12)
      : levelSets({passedLevelSet}), mesh(passedMesh),
        minNodeDistanceFactor(mnd), epsilon(eps) {}

  void setLevelSet(SmartPointer<lsDomainType> passedlsDomain) {
    levelSets = {passedlsDomain};
    normalVectorData = nullptr;
  }

  void setMesh(SmartPointer<Mesh<T>> passedMesh) { mesh = passedMesh; }

  void setUpdatePointData(bool update) { updatePointData = update; }

  void setSharpCorners(bool check) { generateSharpCorners = check; }

  virtual void apply() {
    currentLevelSet = levelSets[0];
    currentGridDelta = currentLevelSet->getGrid().getGridDelta();
    if (currentLevelSet == nullptr) {
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

    if (currentLevelSet->getNumberOfPoints() == 0) {
      VIENNACORE_LOG_WARNING(
          "ToSurfaceMesh: Level set is empty. No mesh will be created.");
      return;
    }

    mesh->clear();
    if (currentLevelSet != nullptr) {
      if (currentLevelSet->getLevelSetWidth() < 2) {
        VIENNACORE_LOG_WARNING("Levelset is less than 2 layers wide. Expanding "
                               "levelset to 2 layers.");
        Expand<T, D>(currentLevelSet, 2).apply();
      }
    }

    ConstSparseIterator valueIt(currentLevelSet->getDomain());

    typedef std::map<hrleIndex, unsigned> nodeContainerType;

    nodeContainerType nodes[D];
    nodeContainerType faceNodes[D];
    nodeContainerType cornerNodes;
    const bool sharpCorners = generateSharpCorners;
    typename nodeContainerType::iterator nodeIt;

    // save how data should be transferred to new level set
    // list of indices into the old pointData vector
    std::vector<std::vector<unsigned>> newDataSourceIds;
    // there is no multithreading here, so just use 1
    if (updatePointData)
      newDataSourceIds.resize(1);

    // If sharp corners are enabled, calculate normals first as they are needed
    // for feature reconstruction
    if (sharpCorners) {
      CalculateNormalVectors<T, D> normalCalculator(currentLevelSet);
      normalCalculator.setMethod(
          NormalCalculationMethodEnum::ONE_SIDED_MIN_MOD);
      normalCalculator.setMaxValue(std::numeric_limits<T>::max());
      normalCalculator.apply();
      normalVectorData = currentLevelSet->getPointData().getVectorData(
          CalculateNormalVectors<T, D>::normalVectorsLabel);
    }

    // iterate over all cells with active points
    for (viennahrle::ConstSparseCellIterator<hrleDomainType> cellIt(
             currentLevelSet->getDomain());
         !cellIt.isFinished(); cellIt.next()) {

      // Clear node caches for dimensions that have moved out of scope
      if (!sharpCorners) {
        for (int u = 0; u < D; u++) {
          while (!nodes[u].empty() &&
                 nodes[u].begin()->first < hrleIndex(cellIt.getIndices()))
            nodes[u].erase(nodes[u].begin());

          while (!faceNodes[u].empty() &&
                 faceNodes[u].begin()->first < hrleIndex(cellIt.getIndices()))
            faceNodes[u].erase(faceNodes[u].begin());

          if (u == 0) {
            while (!cornerNodes.empty() &&
                   cornerNodes.begin()->first < hrleIndex(cellIt.getIndices()))
              cornerNodes.erase(cornerNodes.begin());
          }
        }
      }

      // Calculate signs of all corners to determine the marching cubes case
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

      // Attempt to generate sharp features if enabled
      bool perfectCornerFound = false;
      if (sharpCorners) {
        int countNeg = 0;
        int countPos = 0;
        int negMask = 0;
        int posMask = 0;
        for (int i = 0; i < (1 << D); ++i) {
          T val = cellIt.getCorner(i).getValue();
          if (val < -epsilon) {
            countNeg++;
            negMask |= (1 << i);
          } else if (val > epsilon) {
            countPos++;
            posMask |= (1 << i);
          } else {
            if (val >= 0) {
              countPos++;
              posMask |= (1 << i);
            } else {
              countNeg++;
              negMask |= (1 << i);
            }
          }
        }

        // Check for perfect corner (2D only)
        if constexpr (D == 2) {
          perfectCornerFound = generateSharpCorner2D(
              cellIt, countNeg, countPos, negMask, posMask, nodes,
              &newDataSourceIds[0], valueIt);
        } else if constexpr (D == 3) {
          if (countNeg == 2 || countPos == 2) {
            // Try to generate a sharp edge (2 corners active)
            perfectCornerFound = generateSharpEdge3D(
                cellIt, countNeg, countPos, negMask, posMask, nodes, faceNodes,
                &newDataSourceIds[0], valueIt);
          } else if (countNeg == 1 || countPos == 1) {
            // Try to generate a sharp corner (1 corner active)
            perfectCornerFound = generateSharpCorner3D(
                cellIt, countNeg, countPos, negMask, posMask, nodes,
                cornerNodes, faceNodes, &newDataSourceIds[0], valueIt);
          } else if (countNeg == 3 || countPos == 3) {
            // Try to generate an L-shape corner (3 corners active)
            perfectCornerFound = generateSharpL3D(
                cellIt, countNeg, countPos, negMask, posMask, nodes, faceNodes,
                &newDataSourceIds[0], valueIt);
          }
        }
      }

      // If a sharp feature was successfully generated, skip standard Marching
      // Cubes for this cell
      if (perfectCornerFound)
        continue;

      if constexpr (D == 3) {
        // Stitch to perfect corners/edges
        // Check if neighbors have generated nodes on shared faces that we need
        // to connect to
        for (int axis = 0; axis < 3; ++axis) {
          for (int d = 0; d < 2; ++d) {
            hrleIndex faceKey(cellIt.getIndices());
            if (d == 1)
              faceKey[axis]++;

            auto it = faceNodes[axis].find(faceKey);
            if (it != faceNodes[axis].end()) {
              const unsigned faceNodeId = it->second;
              const int *Triangles =
                  lsInternal::MarchingCubes::polygonize3d(signs);

              for (; Triangles[0] != -1; Triangles += 3) {
                std::vector<unsigned> face_edge_nodes;
                for (int i = 0; i < 3; ++i) {
                  int edge = Triangles[i];
                  int c0 = corner0[edge];
                  int c1 = corner1[edge];
                  bool onFace =
                      (((c0 >> axis) & 1) == d) && (((c1 >> axis) & 1) == d);
                  if (onFace) {
                    face_edge_nodes.push_back(
                        getNode(cellIt, edge, nodes, &newDataSourceIds[0]));
                  }
                }
                if (face_edge_nodes.size() == 2) {
                  insertElement(
                      {face_edge_nodes[0], face_edge_nodes[1], faceNodeId});
                }
              }
            }
          }
        }
      }

      // Standard Marching Cubes / Squares algorithm
      // for each element
      const int *Triangles =
          (D == 2) ? lsInternal::MarchingCubes::polygonize2d(signs)
                   : lsInternal::MarchingCubes::polygonize3d(signs);

      for (; Triangles[0] != -1; Triangles += D) {
        std::array<unsigned, D> nodeNumbers;

        // for each node
        for (int n = 0; n < D; n++) {
          const int edge = Triangles[n];

          unsigned p0 = corner0[edge];
          unsigned p1 = corner1[edge];

          // determine direction of edge
          auto dir = direction[edge];

          // look for existing surface node
          hrleIndex d(cellIt.getIndices());
          d += viennahrle::BitMaskToIndex<D>(p0);

          nodeIt = nodes[dir].find(d);
          if (nodeIt != nodes[dir].end()) {
            nodeNumbers[n] = nodeIt->second;
          } else { // if node does not exist yet
            nodeNumbers[n] = getNode(cellIt, edge, nodes, &newDataSourceIds[0]);
          }
        }

        insertElement(nodeNumbers); // insert new surface element
      }
    }

    scaleMesh();

    // now copy old data into new level set
    if (updatePointData && !newDataSourceIds[0].empty()) {
      mesh->getPointData().translateFromMultiData(
          currentLevelSet->getPointData(), newDataSourceIds);
    }
  }

protected:
  void scaleMesh() {
    // Scale the mesh to global coordinates
    for (auto &node : mesh->nodes) {
      for (int i = 0; i < D; ++i)
        node[i] *= currentGridDelta;
    }
    for (int i = 0; i < D; ++i) {
      mesh->minimumExtent[i] *= currentGridDelta;
      mesh->maximumExtent[i] *= currentGridDelta;
    }
  }

  // Helper to get or create a node on an edge using linear interpolation
  // Compute the interpolated position of a node on an edge without inserting it
  Vec3D<T> computeNodePosition(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt, int edge) {
    unsigned p0 = corner0[edge];
    unsigned p1 = corner1[edge];
    auto dir = direction[edge];

    Vec3D<T> cc{};
    for (int z = 0; z < D; z++) {
      if (z != dir) {
        cc[z] = static_cast<double>(cellIt.getIndices(z) +
                                    viennahrle::BitMaskToIndex<D>(p0)[z]);
      } else {
        T d0 = cellIt.getCorner(p0).getValue();
        T d1 = cellIt.getCorner(p1).getValue();
        if (d0 == -d1) {
          cc[z] = static_cast<T>(cellIt.getIndices(z)) + T(0.5);
        } else {
          if (std::abs(d0) <= std::abs(d1)) {
            cc[z] = static_cast<T>(cellIt.getIndices(z)) + (d0 / (d0 - d1));
          } else {
            cc[z] = static_cast<T>(cellIt.getIndices(z) + 1) - (d1 / (d1 - d0));
          }
        }
        cc[z] = std::max(cc[z], cellIt.getIndices(z) + epsilon);
        cc[z] = std::min(cc[z], (cellIt.getIndices(z) + 1) - epsilon);
      }
    }
    return cc;
  }

  unsigned getNode(viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt,
                   int edge, std::map<hrleIndex, unsigned> *nodes,
                   std::vector<unsigned> *newDataSourceIds) {

    unsigned p0 = corner0[edge];
    unsigned p1 = corner1[edge];
    auto dir = direction[edge];
    hrleIndex d(cellIt.getIndices());
    d += viennahrle::BitMaskToIndex<D>(p0);

    if (nodes[dir].find(d) != nodes[dir].end()) {
      return nodes[dir][d];
    }

    // Create node
    Vec3D<T> cc = computeNodePosition(cellIt, edge);

    std::size_t currentPointId = 0;
    T d0 = cellIt.getCorner(p0).getValue();
    T d1 = cellIt.getCorner(p1).getValue();
    if (std::abs(d0) <= std::abs(d1)) {
      currentPointId = cellIt.getCorner(p0).getPointId();
    } else {
      currentPointId = cellIt.getCorner(p1).getPointId();
    }
    if (updatePointData && newDataSourceIds)
      newDataSourceIds->push_back(currentPointId);

    unsigned nodeId = insertNode(cc);
    nodes[dir][d] = nodeId;
    return nodeId;
  }

  // Helper to stitch a sharp feature node on a face to the standard mesh in a
  // neighboring cell
  void
  stitchToNeighbor(viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt,
                   int axis, bool isHighFace, unsigned faceNodeId,
                   std::map<hrleIndex, unsigned> *nodes,
                   ConstSparseIterator &valueIt) {
    // Backward stitching: Check if neighbor on this face is "past" and needs
    // stitching
    hrleIndex neighborIdx = cellIt.getIndices();
    if (isHighFace)
      neighborIdx[axis]++;
    else
      neighborIdx[axis]--;

    if (neighborIdx < cellIt.getIndices()) {
      unsigned nSigns = 0;
      auto &grid = currentLevelSet->getGrid();
      for (int i = 0; i < 8; ++i) {
        hrleIndex cIdx = neighborIdx + viennahrle::BitMaskToIndex<D>(i);
        if (!grid.isOutsideOfDomain(cIdx)) {
          valueIt.goToIndices(cIdx);
          if (valueIt.getValue() >= 0)
            nSigns |= (1 << i);
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
            bool onFace = (((c0 >> axis) & 1) == nFaceD) &&
                          (((c1 >> axis) & 1) == nFaceD);
            if (onFace) {
              unsigned p0 = corner0[edge];
              auto dir = direction[edge];
              hrleIndex d = neighborIdx + viennahrle::BitMaskToIndex<D>(p0);
              auto itN = nodes[dir].find(d);
              if (itN != nodes[dir].end()) {
                face_edge_nodes.push_back(itN->second);
              }
            }
          }
          if (face_edge_nodes.size() == 2) {
            insertElement(std::array<unsigned, 3>{
                face_edge_nodes[0], face_edge_nodes[1], faceNodeId});
          }
        }
      }
    }
  }

  // Solves for a sharp corner position in 2D by intersecting normals
  bool generateCanonicalSharpCorner2D(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt,
      int transform, Vec3D<T> &cornerPos) const {
    assert(normalVectorData &&
           "Normal vector data must be computed for sharp corner generation.");

    auto getTransformedGradient = [&](int cornerID) {
      auto corner = cellIt.getCorner(cornerID ^ transform);
      if (corner.isDefined()) {
        auto normal = (*normalVectorData)[corner.getPointId()];
        normal[2] = 0;
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

    auto calculateNodePos = [&](int edge, Vec3D<T> &pos) {
      unsigned p0 = corner0[edge];
      unsigned p1 = corner1[edge];
      auto dir = direction[edge];

      for (int z = 0; z < D; z++) {
        if (z != dir) {
          pos[z] = static_cast<T>(viennahrle::BitMaskToIndex<D>(p0)[z]);
        } else {
          const T d0 = cellIt.getCorner(p0 ^ transform).getValue();
          const T d1 = cellIt.getCorner(p1 ^ transform).getValue();
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

    Vec3D<T> pX{}, pY{};
    calculateNodePos(0, pX); // Edge along x-axis from corner 0
    calculateNodePos(3, pY); // Edge along y-axis from corner 0

    double d1 = DotProduct(norm1, pX);
    double d2 = DotProduct(norm2, pY);
    double det = norm1[0] * norm2[1] - norm1[1] * norm2[0];

    if (std::abs(det) > 1e-6 && std::isfinite(det)) {
      cornerPos[0] = (d1 * norm2[1] - d2 * norm1[1]) / det;
      cornerPos[1] = (d2 * norm1[0] - d1 * norm2[0]) / det;
    } else {
      for (int i = 0; i < D; ++i)
        cornerPos[i] = (pX[i] + pY[i]) * T(0.5);
    }

    return true;
  }

  // Wrapper for 2D sharp corner generation handling different
  // rotations/reflections
  bool generateSharpCorner2D(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt, int countNeg,
      int countPos, int negMask, int posMask,
      std::map<hrleIndex, unsigned> *nodes,
      std::vector<unsigned> *newDataSourceIds, ConstSparseIterator &valueIt) {

    int cornerIdx = -1;
    if (countNeg == 1) {
      for (int i = 0; i < 4; ++i)
        if ((negMask >> i) & 1) {
          cornerIdx = i;
          break;
        }
    } else if (countPos == 1) {
      for (int i = 0; i < 4; ++i)
        if ((posMask >> i) & 1) {
          cornerIdx = i;
          break;
        }
    }

    if (cornerIdx != -1) {
      // Check if this corner is also a corner in any neighboring cell
      bool isSharedCorner = false;
      auto pIdx =
          cellIt.getIndices() + viennahrle::BitMaskToIndex<D>(cornerIdx);
      T pVal = cellIt.getCorner(cornerIdx).getValue();
      auto &grid = currentLevelSet->getGrid();

      for (int i = 0; i < (1 << D); ++i) {
        if (i == cornerIdx)
          continue;

        hrleIndex neighborIndices;
        for (int k = 0; k < D; ++k)
          neighborIndices[k] = pIdx[k] - ((i >> k) & 1);

        bool neighborIsCorner = true;

        // Check edge-connected neighbors in the neighbor cell
        for (int k = 0; k < D; ++k) {
          int neighborLocal = i ^ (1 << k);
          auto checkIdx =
              neighborIndices + viennahrle::BitMaskToIndex<D>(neighborLocal);
          if (grid.isOutsideOfDomain(checkIdx)) {
            checkIdx = grid.globalIndices2LocalIndices(checkIdx);
          }
          valueIt.goToIndices(checkIdx);
          T nVal = valueIt.getValue();
          bool pSign = (pVal >= 0);
          bool nSign = (nVal >= 0);
          if (pSign == nSign) {
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
        Vec3D<T> cornerPos{};
        if (generateCanonicalSharpCorner2D(cellIt, cornerIdx, cornerPos)) {

          // inverse transform cornerPos from canonical local to this cell's
          // local
          if ((cornerIdx & 1) != 0)
            cornerPos[0] = T(1.0) - cornerPos[0];
          if ((cornerIdx & 2) != 0)
            cornerPos[1] = T(1.0) - cornerPos[1];

          // convert to global coordinates
          for (int i = 0; i < D; ++i) {
            cornerPos[i] = (cornerPos[i] + cellIt.getIndices(i));
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

          unsigned nX = getNode(cellIt, edgeX, nodes, newDataSourceIds);
          unsigned nY = getNode(cellIt, edgeY, nodes, newDataSourceIds);
          unsigned nCorner = insertNode(cornerPos);

          // Store this sharp corner node for potential use by derived classes
          matSharpCornerNodes.push_back({nCorner, cornerPos});

          if (updatePointData && newDataSourceIds) {
            assert(!newDataSourceIds->empty());
            newDataSourceIds->push_back(
                newDataSourceIds->back()); // TODO: improve point data source
          }

          insertElement({nX, nCorner});
          insertElement({nCorner, nY});
          return true;
        }
      }
    }
    return false;
  }

  // Generates geometry for an "L-shaped" configuration (3 active corners) in 3D
  bool
  generateSharpL3D(viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt,
                   int countNeg, int countPos, int negMask, int posMask,
                   std::map<hrleIndex, unsigned> *nodes,
                   std::map<hrleIndex, unsigned> *faceNodes,
                   std::vector<unsigned> *newDataSourceIds,
                   ConstSparseIterator &valueIt) {

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
          if ((mask >> (i ^ (1 << k))) & 1)
            neighbors++;
        }
        if (neighbors == 2) {
          C = i;
          break;
        }
      }
    }

    if (C == -1)
      return false;

    // Identify A and B
    for (int k = 0; k < 3; ++k) {
      int n = C ^ (1 << k);
      if ((mask >> n) & 1) {
        if (A == -1)
          A = n;
        else
          B = n;
      }
    }

    assert(normalVectorData &&
           "Normal vector data must be computed for sharp feature generation.");

    auto getNormal = [&](int idx) {
      auto corner = cellIt.getCorner(idx);
      if (corner.isDefined()) {
        Vec3D<T> n = (*normalVectorData)[corner.getPointId()];
        if (inverted)
          for (int i = 0; i < 3; ++i)
            n[i] = -n[i];
        return n;
      }
      return Vec3D<T>{};
    };

    // Determine axes
    int axisA = 0;
    while ((C ^ A) != (1 << axisA))
      axisA++;
    int axisB = 0;
    while ((C ^ B) != (1 << axisB))
      axisB++;
    int axisZ = 3 - axisA - axisB;

    // Solve 2D face problem to find sharp point on a face
    auto calculateFace = [&](int axis, int corner, int n1,
                             int n2) -> std::optional<Vec3D<T>> {
      hrleIndex faceIdx = cellIt.getIndices();
      if ((corner >> axis) & 1)
        faceIdx[axis]++;

      if (faceNodes[axis].find(faceIdx) != faceNodes[axis].end()) {
        unsigned nodeId = faceNodes[axis][faceIdx];
        Vec3D<T> pos = mesh->nodes[nodeId];
        for (int d = 0; d < 3; ++d)
          pos[d] = pos[d] - cellIt.getIndices(d);
        return pos;
      }

      Vec3D<T> N1 = getNormal(n1);
      Vec3D<T> N2 = getNormal(n2);

      if (std::abs(DotProduct(N1, N2)) >= 0.8)
        return std::nullopt;

      int u = (axis + 1) % 3;
      int v = (axis + 2) % 3;

      int axis1 = 0;
      while ((corner ^ n1) != (1 << axis1))
        axis1++;
      T t1 = getInterp(corner, n1, cellIt, inverted);

      int axis2 = 0;
      while ((corner ^ n2) != (1 << axis2))
        axis2++;
      T t2 = getInterp(corner, n2, cellIt, inverted);

      Vec3D<T> P1;
      P1[0] = ((corner >> 0) & 1);
      P1[1] = ((corner >> 1) & 1);
      P1[2] = ((corner >> 2) & 1);
      P1[axis1] = ((corner >> axis1) & 1) ? (1.0 - t1) : t1;

      Vec3D<T> P2;
      P2[0] = ((corner >> 0) & 1);
      P2[1] = ((corner >> 1) & 1);
      P2[2] = ((corner >> 2) & 1);
      P2[axis2] = ((corner >> axis2) & 1) ? (1.0 - t2) : t2;

      double det = N1[u] * N2[v] - N1[v] * N2[u];
      if (std::abs(det) < 1e-6)
        return std::nullopt;

      double c1 = N1[u] * P1[u] + N1[v] * P1[v];
      double c2 = N2[u] * P2[u] + N2[v] * P2[v];

      Vec3D<T> P;
      P[axis] = P1[axis];
      P[u] = (c1 * N2[v] - c2 * N1[v]) / det;
      P[v] = (c2 * N1[u] - c1 * N2[u]) / det;

      if (Norm2(P - P1) > 9.0 || Norm2(P - P2) > 9.0)
        return std::nullopt;

      return P;
    };

    auto commitFace = [&](int axis, int corner, Vec3D<T> P) -> unsigned {
      hrleIndex faceIdx = cellIt.getIndices();
      if ((corner >> axis) & 1)
        faceIdx[axis]++;

      if (faceNodes[axis].find(faceIdx) != faceNodes[axis].end()) {
        return faceNodes[axis][faceIdx];
      }

      Vec3D<T> globalP = P;
      for (int d = 0; d < 3; ++d)
        globalP[d] = (globalP[d] + cellIt.getIndices(d));
      unsigned nodeId = insertNode(globalP);
      faceNodes[axis][faceIdx] = nodeId;
      if (updatePointData && newDataSourceIds)
        newDataSourceIds->push_back(cellIt.getCorner(corner).getPointId());

      stitchToNeighbor(cellIt, axis, (corner >> axis) & 1, nodeId, nodes,
                       valueIt);

      return nodeId;
    };

    int D_corner = A ^ (1 << axisB);
    int A_z = A ^ (1 << axisZ);
    int B_z = B ^ (1 << axisZ);
    int C_z = C ^ (1 << axisZ);

    auto P_A_opt = calculateFace(axisA, A, D_corner, A_z);
    auto P_B_opt = calculateFace(axisB, B, D_corner, B_z);

    if (!P_A_opt.has_value() || !P_B_opt.has_value())
      return false;

    auto const &P_A = P_A_opt.value();
    auto const &P_B = P_B_opt.value();

    // Construct S
    Vec3D<T> S;
    S[axisA] = P_B[axisA];
    S[axisB] = P_A[axisB];
    S[axisZ] = (P_A[axisZ] + P_B[axisZ]) * 0.5;

    unsigned nS = insertNode([&]() {
      Vec3D<T> gS = S;
      for (int i = 0; i < 3; ++i)
        gS[i] = (gS[i] + cellIt.getIndices(i));
      return gS;
    }());
    if (updatePointData && newDataSourceIds)
      newDataSourceIds->push_back(cellIt.getCorner(C).getPointId());

    // Calculate S_Face
    Vec3D<T> S_Face = S;
    S_Face[axisZ] = static_cast<T>((C >> axisZ) & 1);

    unsigned nS_Face;
    auto faceIdxZ = cellIt.getIndices();
    if ((C >> axisZ) & 1)
      faceIdxZ[axisZ]++;

    if (faceNodes[axisZ].find(faceIdxZ) != faceNodes[axisZ].end()) {
      nS_Face = faceNodes[axisZ][faceIdxZ];
    } else {
      nS_Face = insertNode([&]() {
        Vec3D<T> gS = S_Face;
        for (int i = 0; i < 3; ++i)
          gS[i] = (gS[i] + cellIt.getIndices(i));
        return gS;
      }());
      faceNodes[axisZ][faceIdxZ] = nS_Face;
      if (updatePointData)
        newDataSourceIds[0].push_back(cellIt.getCorner(C).getPointId());

      stitchToNeighbor(cellIt, axisZ, (C >> axisZ) & 1, nS_Face, nodes,
                       valueIt);
    }

    // Get face nodes
    unsigned nNA = commitFace(axisA, A, P_A);
    unsigned nNB = commitFace(axisB, B, P_B);

    // Get boundary intersection nodes
    auto getEdgeNode = [&](int v1, int v2) {
      int edgeIdx = -1;
      for (int e = 0; e < 12; ++e) {
        if ((corner0[e] == static_cast<unsigned>(v1) &&
             corner1[e] == static_cast<unsigned>(v2)) ||
            (corner0[e] == static_cast<unsigned>(v2) &&
             corner1[e] == static_cast<unsigned>(v1))) {
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
    int parity = (C & 1) + ((C >> 1) & 1) + ((C >> 2) & 1);
    bool flip = (parity % 2 != 0) ^ inverted;

    int vA = 0;
    while ((A ^ C) != (1 << vA))
      vA++;
    int vB = 0;
    while ((B ^ C) != (1 << vB))
      vB++;
    bool is_cyclic = ((vA + 1) % 3 == vB);
    if (!is_cyclic) {
      std::swap(nNA, nNB);
      std::swap(nI_AD, nI_BD);
      std::swap(nI_AZ, nI_BZ);
    }

    if (!flip) {
      insertElement({nS, nNA, nI_AD});
      insertElement({nS, nI_AD, nS_Face});
      insertElement({nS, nS_Face, nI_BD});
      insertElement({nS, nI_BD, nNB});
      insertElement({nS, nNB, nI_BZ});
      insertElement({nS, nI_BZ, nI_CZ});
      insertElement({nS, nI_CZ, nI_AZ});
      insertElement({nS, nI_AZ, nNA});
    } else {
      insertElement({nS, nI_AD, nNA});
      insertElement({nS, nS_Face, nI_AD});
      insertElement({nS, nI_BD, nS_Face});
      insertElement({nS, nNB, nI_BD});
      insertElement({nS, nI_BZ, nNB});
      insertElement({nS, nI_CZ, nI_BZ});
      insertElement({nS, nI_AZ, nI_CZ});
      insertElement({nS, nNA, nI_AZ});
    }

    return true;
  }

  // Generates geometry for a sharp edge (2 active corners) in 3D (canonical
  // orientation)
  bool generateCanonicalSharpEdge3D(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt,
      int transform, int axis, bool inverted,
      std::map<hrleIndex, unsigned> *nodes,
      std::map<hrleIndex, unsigned> *faceNodes,
      std::vector<unsigned> *newDataSourceIds, ConstSparseIterator &valueIt) {
    assert(normalVectorData &&
           "Normal vector data must be computed for sharp edge generation.");

    // Helper to map global vector to canonical frame
    auto toCanonical = [&](Vec3D<T> v) {
      Vec3D<T> res;
      // 1. Apply reflection (transform)
      if (transform & 1)
        v[0] = -v[0];
      if (transform & 2)
        v[1] = -v[1];
      if (transform & 4)
        v[2] = -v[2];
      // 2. Apply rotation (axis permutation)
      // axis=0: XYZ -> XYZ (0,1,2)
      // axis=1: XYZ -> YZX (1,2,0) - puts Y in X place
      // axis=2: XYZ -> ZXY (2,0,1) - puts Z in X place
      if (axis == 0)
        res = v;
      else if (axis == 1)
        res = Vec3D<T>{v[1], v[2], v[0]};
      else
        res = Vec3D<T>{v[2], v[0], v[1]};
      return res;
    };

    // Helper to map canonical point back to global frame
    auto fromCanonical = [&](Vec3D<T> v) {
      Vec3D<T> res;
      // 1. Inverse rotation
      if (axis == 0)
        res = v;
      else if (axis == 1)
        res = Vec3D<T>{v[2], v[0], v[1]};
      else
        res = Vec3D<T>{v[1], v[2], v[0]};
      // 2. Inverse reflection (same as forward for 0/1 swap)
      if (transform & 1)
        res[0] = T(1.0) - res[0];
      if (transform & 2)
        res[1] = T(1.0) - res[1];
      if (transform & 4)
        res[2] = T(1.0) - res[2];
      return res;
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
    bool nodeExists[2] = {false, false};

    // Inverse map indices helper
    auto mapIdx = [&](int c) {
      // c is canonical index.
      // 1. Inverse rotation:
      int x = (c & 1), y = (c >> 1) & 1, z = (c >> 2) & 1;
      int gx, gy, gz;
      if (axis == 0) {
        gx = x;
        gy = y;
        gz = z;
      } else if (axis == 1) {
        gx = z;
        gy = x;
        gz = y;
      } else {
        gx = y;
        gy = z;
        gz = x;
      }

      // 2. Inverse reflection
      if (transform & 1)
        gx = 1 - gx;
      if (transform & 2)
        gy = 1 - gy;
      if (transform & 4)
        gz = 1 - gz;

      return gx | (gy << 1) | (gz << 2);
    };

    // First pass: Calculate and check validity without modifying mesh
    for (int k = 0; k < 2; ++k) {
      // Check if face node already exists
      // k=0 -> x=0 face (left). In global, this corresponds to the face
      // perpendicular to 'axis' at the 'transform' side. If transform has bit
      // 'axis' set, then x=0 in canonical is x=1 in global (relative to cell
      // origin). We need to be careful with the face index key. Let's compute
      // the global index of the face.
      hrleIndex faceIdx = cellIt.getIndices();
      // If k=0 (canonical x=0), and transform has bit 'axis' set (flipped),
      // this is global x=1 face. If k=1 (canonical x=1), and transform has bit
      // 'axis' set, this is global x=0 face.
      bool isHighFace = (k == 1) ^ ((transform >> axis) & 1);
      if (isHighFace)
        faceIdx[axis]++;

      auto it = faceNodes[axis].find(faceIdx);
      if (it != faceNodes[axis].end()) {
        nodeExists[k] = true;
        faceNodeIds[k] = it->second;
        Vec3D<T> relPos = mesh->nodes[faceNodeIds[k]];
        for (int d = 0; d < 3; ++d) {
          relPos[d] -= static_cast<T>(cellIt.getIndices(d));
        }
        P[k] = toCanonical(relPos);
      } else {
        // Calculate P[k]
        // For k=0 (x=0): Active corner is 0 (000). Neighbors 2 (010) and 4
        // (001). For k=1 (x=1): Active corner is 1 (100). Neighbors 3 (110) and
        // 5 (101). We need to map these canonical indices back to global to get
        // normals/values.

        int c0 = k;       // 0 or 1
        int c_y = c0 | 2; // Neighbor in Y (canonical)
        int c_z = c0 | 4; // Neighbor in Z (canonical)

        c0 = mapIdx(c0);
        c_y = mapIdx(c_y);
        c_z = mapIdx(c_z);

        Vec3D<T> n_y = toCanonical(getNormal(c_y));
        Vec3D<T> n_z = toCanonical(getNormal(c_z));

        if (std::abs(DotProduct(n_y, n_z)) >= 0.8) {
          return false;
        }

        // Solve 2D problem in YZ plane
        // Line 1: passes through intersection on edge c0-c_y. Normal n_y
        // (projected). Line 2: passes through intersection on edge c0-c_z.
        // Normal n_z (projected).
        T t_y = getInterp(c0, c_y, cellIt, inverted); // Fraction along Y axis
        T t_z = getInterp(c0, c_z, cellIt, inverted); // Fraction along Z axis

        // In canonical local frame relative to c0:
        // Point on Y-edge: (0, t_y, 0)
        // Point on Z-edge: (0, 0, t_z)
        // Normals n_y and n_z.
        // Solve for (y, z):
        // n_y.y * (y - t_y) + n_y.z * (z - 0) = 0
        // n_z.y * (y - 0) + n_z.z * (z - t_z) = 0

        double det = n_y[1] * n_z[2] - n_y[2] * n_z[1];
        if (std::abs(det) < 1e-6)
          return false;

        double d1 = n_y[1] * t_y;
        double d2 = n_z[2] * t_z;

        double y = (d1 * n_z[2] - d2 * n_y[2]) / det;
        double z = (d2 * n_y[1] - d1 * n_z[1]) / det;

        P[k] = Vec3D<T>{(T)k, (T)y, (T)z};

        // Check if the point is too far away from the intersection points
        // 2 grid deltas = 2.0 in canonical coordinates
        if (Norm2(P[k] - Vec3D<T>{(T)k, t_y, 0}) > 9.0 ||
            Norm2(P[k] - Vec3D<T>{(T)k, 0, t_z}) > 9.0) {
          return false;
        }
      }
    }

    // Second pass: Commit nodes and stitch
    for (int k = 0; k < 2; ++k) {
      if (!nodeExists[k]) {
        int c0 = k;
        hrleIndex faceIdx = cellIt.getIndices();
        bool isHighFace = (k == 1) ^ ((transform >> axis) & 1);
        if (isHighFace)
          faceIdx[axis]++;

        Vec3D<T> globalP = fromCanonical(P[k]);
        for (int d = 0; d < 3; ++d)
          globalP[d] = (globalP[d] + cellIt.getIndices(d));

        faceNodeIds[k] = insertNode(globalP);
        faceNodes[axis][faceIdx] = faceNodeIds[k];
        if (updatePointData && newDataSourceIds)
          newDataSourceIds->push_back(
              cellIt.getCorner(mapIdx(c0)).getPointId());

        stitchToNeighbor(cellIt, axis, isHighFace, faceNodeIds[k], nodes,
                         valueIt);
      }
    }

    // Calculate edge intersection points in canonical frame
    // E02: on edge 0-2 (Y axis). (0, t_02, 0)
    // E04: on edge 0-4 (Z axis). (0, 0, t_04)
    // E13: on edge 1-3 (Y axis). (1, t_13, 0)
    // E15: on edge 1-5 (Z axis). (1, 0, t_15)

    // We need to get the node IDs for these. We can use the existing getNode
    // helper, but we need to map canonical edges to global edges.
    auto getEdgeNode = [&](int c1, int c2) {
      // Map canonical corners to global
      int g1 = mapIdx(c1);
      int g2 = mapIdx(c2);

      // Find edge index
      int edgeIdx = -1;
      for (int e = 0; e < 12; ++e) {
        if ((corner0[e] == static_cast<unsigned>(g1) &&
             corner1[e] == static_cast<unsigned>(g2)) ||
            (corner0[e] == static_cast<unsigned>(g2) &&
             corner1[e] == static_cast<unsigned>(g1))) {
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
    insertElement({n02, faceNodeIds[0], faceNodeIds[1]});
    insertElement({n02, faceNodeIds[1], n13});

    // Quad 2 (Z-interface): E04, E15, P1, P0. Normal +Z.
    // Winding for +Z normal: E04 -> E15 -> P1 -> P0.
    insertElement({n04, n15, faceNodeIds[1]});
    insertElement({n04, faceNodeIds[1], faceNodeIds[0]});

    return true;
  }

  // Wrapper for 3D sharp edge generation handling different
  // rotations/reflections
  bool generateSharpEdge3D(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt, int countNeg,
      int countPos, int negMask, int posMask,
      std::map<hrleIndex, unsigned> *nodes,
      std::map<hrleIndex, unsigned> *faceNodes,
      std::vector<unsigned> *newDataSourceIds, ConstSparseIterator &valueIt) {

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
        if (v1 == -1)
          v1 = i;
        else
          v2 = i;
      }
    }

    int diff = v1 ^ v2;
    if ((diff & (diff - 1)) != 0)
      return false; // Not connected by a single edge

    int axis = 0;
    while ((diff >> (axis + 1)) > 0)
      axis++;

    // Transform maps v1 to 0.
    int transform = v1;

    return generateCanonicalSharpEdge3D(cellIt, transform, axis, inverted,
                                        nodes, faceNodes, newDataSourceIds,
                                        valueIt);
  }

  // Generates geometry for a sharp corner (1 active corner) in 3D (canonical
  // orientation)
  bool generateCanonicalSharpCorner3D(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt,
      int transform, bool inverted, std::map<hrleIndex, unsigned> *nodes,
      std::map<hrleIndex, unsigned> *faceNodes,
      std::vector<unsigned> *newDataSourceIds, ConstSparseIterator &valueIt,
      Vec3D<T> &cornerPos) {

    const auto *normalVectorData =
        currentLevelSet->getPointData().getVectorData(
            viennals::CalculateNormalVectors<T, D>::normalVectorsLabel);
    assert(normalVectorData &&
           "Normal vector data must be available for sharp corner generation");

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

    Vec3D<T> norm1, norm2, norm3;
    double d1, d2, d3;
    Vec3D<T> P1, P2, P4;

    const int g1 = transform ^ 1;
    const int g2 = transform ^ 2;
    const int g4 = transform ^ 4;

    norm1 = toCanonical(getNormal(g1));
    norm2 = toCanonical(getNormal(g2));
    norm3 = toCanonical(getNormal(g4));

    const T t1 = getInterp(transform, g1, cellIt, inverted);
    const T t2 = getInterp(transform, g2, cellIt, inverted);
    const T t4 = getInterp(transform, g4, cellIt, inverted);

    P1 = {t1, 0, 0};
    P2 = {0, t2, 0};
    P4 = {0, 0, t4};

    d1 = DotProduct(norm1, P1);
    d2 = DotProduct(norm2, P2);
    d3 = DotProduct(norm3, P4);

    if (std::abs(DotProduct(norm1, norm2)) >= 0.8 ||
        std::abs(DotProduct(norm1, norm3)) >= 0.8 ||
        std::abs(DotProduct(norm2, norm3)) >= 0.8) {
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

    if (Norm2(S - P1) > 9.0 || Norm2(S - P2) > 9.0 || Norm2(S - P4) > 9.0) {
      return false;
    }

    // Transform corner position back from canonical to global coordinates
    Vec3D<T> globalS = fromCanonical(S);
    for (int d = 0; d < 3; ++d) {
      cornerPos[d] = (globalS[d] + cellIt.getIndices(d));
    }

    return true;
  }

  // Wrapper for 3D sharp corner generation handling different
  // rotations/reflections
  bool generateSharpCorner3D(
      viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt, int countNeg,
      int countPos, int negMask, int posMask,
      std::map<hrleIndex, unsigned> *nodes,
      std::map<hrleIndex, unsigned> &cornerNodes,
      std::map<hrleIndex, unsigned> *faceNodes,
      std::vector<unsigned> *newDataSourceIds, ConstSparseIterator &valueIt) {

    bool inverted = false;
    int transform = -1;
    // bool is3vs5 = false;

    if (countNeg == 1) {
      transform = 0;
      while (!((negMask >> transform) & 1))
        transform++;
      inverted = false;
    } else if (countPos == 1) {
      transform = 0;
      while (!((posMask >> transform) & 1))
        transform++;
      inverted = true;
    }

    if (transform == -1)
      return false;

    Vec3D<T> cornerPos;
    if (generateCanonicalSharpCorner3D(cellIt, transform, inverted, nodes,
                                       faceNodes, newDataSourceIds, valueIt,
                                       cornerPos)) {
      unsigned nCorner;
      hrleIndex cornerIdx = cellIt.getIndices();
      cornerIdx += viennahrle::BitMaskToIndex<D>(transform);

      auto it = cornerNodes.find(cornerIdx);
      if (it != cornerNodes.end()) {
        nCorner = it->second;
      } else {
        nCorner = insertNode(cornerPos);
        cornerNodes[cornerIdx] = nCorner;

        // Store this sharp corner node for potential use by derived classes
        matSharpCornerNodes.push_back({nCorner, cornerPos});

        if (updatePointData && newDataSourceIds)
          newDataSourceIds->push_back(cellIt.getCorner(transform).getPointId());
      }

      // Get edge nodes
      auto getEdgeNode = [&](int neighbor) {
        int v1 = transform;
        int v2 = neighbor;
        int edgeIdx = -1;
        for (int e = 0; e < 12; ++e) {
          if ((corner0[e] == static_cast<unsigned>(v1) &&
               corner1[e] == static_cast<unsigned>(v2)) ||
              (corner0[e] == static_cast<unsigned>(v2) &&
               corner1[e] == static_cast<unsigned>(v1))) {
            edgeIdx = e;
            break;
          }
        }
        return this->getNode(cellIt, edgeIdx, nodes, newDataSourceIds);
      };

      unsigned n1 = getEdgeNode(transform ^ 1);
      unsigned n2 = getEdgeNode(transform ^ 2);
      unsigned n4 = getEdgeNode(transform ^ 4);

      // Try to generate cube corner topology
      if (normalVectorData) {
        // Re-calculate canonical normals and plane constants
        auto toCanonical = [&](Vec3D<T> v) {
          if (transform & 1)
            v[0] = -v[0];
          if (transform & 2)
            v[1] = -v[1];
          if (transform & 4)
            v[2] = -v[2];
          return v;
        };

        auto fromCanonical = [&](Vec3D<T> v) {
          if (transform & 1)
            v[0] = T(1.0) - v[0];
          if (transform & 2)
            v[1] = T(1.0) - v[1];
          if (transform & 4)
            v[2] = T(1.0) - v[2];
          return v;
        };

        auto getNormal = [&](int idx) {
          auto corner = cellIt.getCorner(idx);
          if (corner.isDefined()) {
            Vec3D<T> n = (*normalVectorData)[corner.getPointId()];
            if (inverted)
              for (int i = 0; i < 3; ++i)
                n[i] = -n[i];
            return n;
          }
          return Vec3D<T>{};
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
            if (n1_idx == -1)
              n1_idx = n_idx;
            else
              n2_idx = n_idx;
          } else {
            n3_idx = n_idx;
          }
        }

        // Standard 1-vs-7 case
        norm1 = toCanonical(getNormal(transform ^ 1));
        norm2 = toCanonical(getNormal(transform ^ 2));
        norm3 = toCanonical(getNormal(transform ^ 4));

        d1 = norm1[0] * getInterp(transform, transform ^ 1, cellIt, inverted);
        d2 = norm2[1] * getInterp(transform, transform ^ 2, cellIt, inverted);
        d3 = norm3[2] * getInterp(transform, transform ^ 4, cellIt, inverted);

        auto solve2x2 = [&](double a1, double b1, double c1, double a2,
                            double b2, double c2, T &res1, T &res2) {
          double det = a1 * b2 - a2 * b1;
          if (std::abs(det) < 1e-6)
            return false;
          res1 = (c1 * b2 - c2 * b1) / det;
          res2 = (a1 * c2 - a2 * c1) / det;
          return true;
        };

        Vec3D<T> Fx, Fy, Fz;
        Fx[0] = 0.0;
        Fy[1] = 0.0;
        Fz[2] = 0.0;

        bool valid = true;
        valid &= solve2x2(norm2[1], norm2[2], d2 - norm2[0] * Fx[0], norm3[1],
                          norm3[2], d3 - norm3[0] * Fx[0], Fx[1], Fx[2]);
        valid &= solve2x2(norm1[0], norm1[2], d1 - norm1[1] * Fy[1], norm3[0],
                          norm3[2], d3 - norm3[1] * Fy[1], Fy[0], Fy[2]);
        valid &= solve2x2(norm1[0], norm1[1], d1 - norm1[2] * Fz[2], norm2[0],
                          norm2[1], d2 - norm2[2] * Fz[2], Fz[0], Fz[1]);

        auto checkBounds = [&](Vec3D<T> p) {
          return p[0] >= -1e-4 && p[0] <= 1.0 + 1e-4 && p[1] >= -1e-4 &&
                 p[1] <= 1.0 + 1e-4 && p[2] >= -1e-4 && p[2] <= 1.0 + 1e-4;
        };

        if (valid && checkBounds(Fx) && checkBounds(Fy) && checkBounds(Fz)) {
          auto getFaceNode = [&](int axis, Vec3D<T> pos) {
            hrleIndex faceIdx = cellIt.getIndices();
            if ((transform >> axis) & 1)
              faceIdx[axis]++;

            if (faceNodes[axis].find(faceIdx) != faceNodes[axis].end()) {
              return faceNodes[axis][faceIdx];
            }
            Vec3D<T> globalPos = fromCanonical(pos);
            for (int d = 0; d < 3; ++d)
              globalPos[d] = (globalPos[d] + cellIt.getIndices(d));
            unsigned nodeId = insertNode(globalPos);
            faceNodes[axis][faceIdx] = nodeId;
            if (updatePointData)
              newDataSourceIds->push_back(
                  cellIt.getCorner(transform).getPointId());
            stitchToNeighbor(cellIt, axis, (transform >> axis) & 1, nodeId,
                             nodes, valueIt);
            return nodeId;
          };

          unsigned nFx = getFaceNode(0, Fx);
          unsigned nFy = getFaceNode(1, Fy);
          unsigned nFz = getFaceNode(2, Fz);

          // Calculate parity of the transform (number of reflections)
          int parity =
              (transform & 1) + ((transform >> 1) & 1) + ((transform >> 2) & 1);

          // Flip winding if parity is odd XOR inverted is true.
          bool flip = (parity % 2 != 0) ^ inverted;

          if (!flip) {
            insertElement({nCorner, nFy, n1});
            insertElement({nCorner, n1, nFz});
            insertElement({nCorner, nFz, n2});
            insertElement({nCorner, n2, nFx});
            insertElement({nCorner, nFx, n4});
            insertElement({nCorner, n4, nFy});
          } else {
            insertElement({nCorner, n1, nFy});
            insertElement({nCorner, nFz, n1});
            insertElement({nCorner, n2, nFz});
            insertElement({nCorner, nFx, n2});
            insertElement({nCorner, n4, nFx});
            insertElement({nCorner, nFy, n4});
          }
          return true;
        }
      }

      // Triangles
      // Cycle around normal (1,1,1) in canonical is 1->2->4->1 (x->y->z->x)
      // Triangles: (n1, S, n4), (n4, S, n2), (n2, S, n1)
      // If inverted, flip winding.

      // Calculate parity of the transform (number of reflections)
      int parity =
          (transform & 1) + ((transform >> 1) & 1) + ((transform >> 2) & 1);

      // Flip winding if parity is odd XOR inverted is true.
      bool flip = (parity % 2 != 0) ^ inverted;

      if (!flip) {
        insertElement({n1, nCorner, n4});
        insertElement({n4, nCorner, n2});
        insertElement({n2, nCorner, n1});
      } else {
        insertElement({n1, n4, nCorner});
        insertElement({n4, n2, nCorner});
        insertElement({n2, n1, nCorner});
      }

      return true;
    }

    return false;
  }

  static inline bool
  triangleMisformed(const std::array<unsigned, D> &nodeNumbers) noexcept {
    if constexpr (D == 3) {
      return nodeNumbers[0] == nodeNumbers[1] ||
             nodeNumbers[0] == nodeNumbers[2] ||
             nodeNumbers[1] == nodeNumbers[2];
    } else {
      return nodeNumbers[0] == nodeNumbers[1];
    }
  }

  static inline Vec3D<T> calculateNormal(const Vec3D<T> &nodeA,
                                         const Vec3D<T> &nodeB,
                                         const Vec3D<T> &nodeC) noexcept {
    return CrossProduct(nodeB - nodeA, nodeC - nodeA);
  }

  T getInterp(int p_a, int p_b,
              const viennahrle::ConstSparseCellIterator<hrleDomainType> &cellIt,
              bool inverted) const {
    auto getValue = [&](int idx) {
      T val = cellIt.getCorner(idx).getValue();
      return inverted ? -val : val;
    };
    const T v_a = getValue(p_a);
    const T v_b = getValue(p_b);
    if (std::abs(v_a - v_b) < epsilon)
      return T(0.5);
    return v_a / (v_a - v_b);
  }

  void insertElement(const std::array<unsigned, D> &nodeNumbers) {
    if (triangleMisformed(nodeNumbers))
      return;

    auto elementI3 = [&](const std::array<unsigned, D> &element) -> I3 {
      if constexpr (D == 2)
        return {(int)element[0], (int)element[1], 0};
      else
        return {(int)element[0], (int)element[1], (int)element[2]};
    };

    if (uniqueElements.insert(elementI3(nodeNumbers)).second) {
      Vec3D<T> normal;
      if constexpr (D == 2) {
        normal = Vec3D<T>{
            -(mesh->nodes[nodeNumbers[1]][1] - mesh->nodes[nodeNumbers[0]][1]),
            mesh->nodes[nodeNumbers[1]][0] - mesh->nodes[nodeNumbers[0]][0],
            T(0)};
      } else {
        normal = calculateNormal(mesh->nodes[nodeNumbers[0]],
                                 mesh->nodes[nodeNumbers[1]],
                                 mesh->nodes[nodeNumbers[2]]);
      }

      double n2 =
          normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
      if (n2 > epsilon) {
        mesh->insertNextElement(nodeNumbers);
        T invn = static_cast<T>(1.) / std::sqrt(static_cast<T>(n2));
        for (int d = 0; d < D; d++) {
          normal[d] *= invn;
        }
        currentNormals.push_back(normal);
        currentMaterials.push_back(static_cast<T>(currentMaterialId));
      }
    }
  }

  unsigned insertNode(Vec3D<T> const &pos) {
    auto quantize = [&](const Vec3D<T> &p) -> I3 {
      const T inv = T(1) / minNodeDistanceFactor;
      return {(int)std::llround(p[0] * inv), (int)std::llround(p[1] * inv),
              (int)std::llround(p[2] * inv)};
    };

    int nodeIdx = -1;
    if (minNodeDistanceFactor > 0) {
      auto q = quantize(pos);
      auto it = nodeIdByBin.find(q);
      if (it != nodeIdByBin.end())
        nodeIdx = it->second;
    }

    if (nodeIdx >= 0) {
      for (int i = 0; i < 3; ++i) {
        mesh->nodes[nodeIdx][i] = (mesh->nodes[nodeIdx][i] + pos[i]) * T(0.5);
      }
      return nodeIdx;
    } else {
      unsigned newNodeId = mesh->insertNextNode(pos);
      if (minNodeDistanceFactor > 0) {
        nodeIdByBin.emplace(quantize(pos), newNodeId);
      }
      for (int a = 0; a < D; a++) {
        if (pos[a] < mesh->minimumExtent[a])
          mesh->minimumExtent[a] = pos[a];
        if (pos[a] > mesh->maximumExtent[a])
          mesh->maximumExtent[a] = pos[a];
      }
      return newNodeId;
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(ToSurfaceMesh)

} // namespace viennals
