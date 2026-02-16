#pragma once

#include <lsMaterialMap.hpp>
#include <lsToSurfaceMesh.hpp>

namespace viennals {

using namespace viennacore;

template <class NumericType, int D>
class ToMultiSurfaceMesh : public ToSurfaceMesh<NumericType, D> {
  using lsDomainType = viennals::Domain<NumericType, D>;
  using hrleDomainType = typename lsDomainType::DomainType;
  using hrleIndex = viennahrle::Index<D>;

  using ToSurfaceMesh<NumericType, D>::mesh;
  using ToSurfaceMesh<NumericType, D>::levelSets;
  using ToSurfaceMesh<NumericType, D>::currentLevelSet;
  using ToSurfaceMesh<NumericType, D>::currentGridDelta;
  using ToSurfaceMesh<NumericType, D>::currentMaterialId;
  using ToSurfaceMesh<NumericType, D>::epsilon;
  using ToSurfaceMesh<NumericType, D>::nodeIdByBin;
  using ToSurfaceMesh<NumericType, D>::uniqueElements;
  using ToSurfaceMesh<NumericType, D>::currentNormals;
  using ToSurfaceMesh<NumericType, D>::currentMaterials;
  using ToSurfaceMesh<NumericType, D>::normalVectorData;
  using ToSurfaceMesh<NumericType, D>::generateSharpCorners;
  using ToSurfaceMesh<NumericType, D>::matSharpCornerNodes;

  SmartPointer<MaterialMap> materialMap = nullptr;

  // Structure to track sharp corner nodes from each material
  struct SharpCornerNode {
    unsigned nodeId;
    Vec3D<NumericType> position;
  };

  // Find the closest sharp corner node from materials below l.
  // Returns the node ID or -1 if none found within threshold.
  static int findNearestLowerCorner(
      unsigned l, const Vec3D<NumericType> &pos,
      const std::unordered_map<unsigned, std::vector<SharpCornerNode>>
          &sharpCornerNodes) {
    int snapId = -1;
    NumericType minDist2 = 0.1;
    for (unsigned m = 0; m < l; ++m) {
      auto it = sharpCornerNodes.find(m);
      if (it == sharpCornerNodes.end())
        continue;
      for (const auto &scn : it->second) {
        NumericType dist2 = 0;
        for (int i = 0; i < D; ++i) {
          NumericType d = scn.position[i] - pos[i];
          dist2 += d * d;
        }
        if (dist2 < minDist2) {
          minDist2 = dist2;
          snapId = scn.nodeId;
        }
      }
    }
    return snapId;
  }

  // In 2D, check if cornerPos lies on any edge of a previous material's mesh.
  // If so, split that edge to include cornerNodeId.
  void splitPreviousMaterialEdge(unsigned l, unsigned cornerNodeId,
                                 const Vec3D<NumericType> &cornerPos,
                                 size_t meshSizeBefore) {
    if constexpr (D != 2)
      return;
    if (l == 0)
      return;

    for (size_t lineIdx = 0; lineIdx < mesh->lines.size(); ++lineIdx) {
      unsigned mat = currentMaterials[lineIdx];
      if (mat >= l)
        continue;

      auto &line = mesh->lines[lineIdx];
      unsigned node0 = line[0];
      unsigned node1 = line[1];

      auto const &pos0 = mesh->nodes[node0];
      auto const &pos1 = mesh->nodes[node1];
      auto lineVec = pos1 - pos0;

      NumericType lineLen2 = 0;
      NumericType dot = 0;
      for (int i = 0; i < D; ++i) {
        lineLen2 += lineVec[i] * lineVec[i];
        dot += (cornerPos[i] - pos0[i]) * lineVec[i];
      }

      if (lineLen2 < epsilon)
        continue;

      NumericType t = dot / lineLen2;
      if (t <= 0 || t >= 1)
        continue;

      // Distance from corner to closest point on line
      NumericType dist2 = 0;
      for (int i = 0; i < D; ++i) {
        NumericType d = cornerPos[i] - (pos0[i] + t * lineVec[i]);
        dist2 += d * d;
      }
      if (dist2 >= 0.01)
        continue;

      // Corner is on this edge - find which recent lines already cover
      // the split segments
      auto reassignSegment = [&](int idx, unsigned a, unsigned b) {
        currentMaterials[idx] = mat;
        mesh->lines[idx] = {a, b};
        currentNormals[idx] = currentNormals[lineIdx];
      };

      int seg1Idx = -1;
      int seg2Idx = -1;
      for (size_t k = meshSizeBefore; k < mesh->lines.size(); ++k) {
        const auto &rl = mesh->lines[k];
        if ((rl[0] == node0 && rl[1] == cornerNodeId) ||
            (rl[0] == cornerNodeId && rl[1] == node0))
          seg1Idx = k;
        if ((rl[0] == cornerNodeId && rl[1] == node1) ||
            (rl[0] == node1 && rl[1] == cornerNodeId))
          seg2Idx = k;
      }

      if (seg1Idx != -1 && seg2Idx != -1) {
        reassignSegment(seg1Idx, node0, cornerNodeId);
        reassignSegment(seg2Idx, cornerNodeId, node1);
        mesh->lines.erase(mesh->lines.begin() + lineIdx);
        currentMaterials.erase(currentMaterials.begin() + lineIdx);
        currentNormals.erase(currentNormals.begin() + lineIdx);
      } else if (seg1Idx != -1) {
        reassignSegment(seg1Idx, node0, cornerNodeId);
        line[0] = cornerNodeId;
        line[1] = node1;
      } else if (seg2Idx != -1) {
        reassignSegment(seg2Idx, cornerNodeId, node1);
        line[1] = cornerNodeId;
      } else {
        line[1] = cornerNodeId;
        mesh->lines.push_back({cornerNodeId, node1});
        currentMaterials.push_back(mat);
        currentNormals.push_back(currentNormals[lineIdx]);
      }

      break; // Corner can only be on one edge
    }
  }

  // Handle snapping/splitting of sharp corners between materials.
  void
  snapSharpCorners(unsigned l, size_t meshSizeBefore,
                   std::unordered_map<unsigned, std::vector<SharpCornerNode>>
                       &sharpCornerNodes) {
    for (const auto &cornerPair : matSharpCornerNodes) {
      unsigned cornerNodeId = cornerPair.first;
      Vec3D<NumericType> cornerPos = cornerPair.second;

      int snapToNodeId =
          (l > 0) ? findNearestLowerCorner(l, cornerPos, sharpCornerNodes) : -1;

      if (snapToNodeId >= 0) {
        // Snap: replace cornerNodeId with snapToNodeId in recent elements
        if constexpr (D == 2) {
          using I3 = typename ToSurfaceMesh<NumericType, D>::I3;
          for (size_t i = mesh->lines.size(); i-- > meshSizeBefore;) {
            for (int j = 0; j < 2; ++j) {
              if (mesh->lines[i][j] == cornerNodeId)
                mesh->lines[i][j] = snapToNodeId;
            }

            unsigned n0 = mesh->lines[i][0];
            unsigned n1 = mesh->lines[i][1];
            bool isDuplicate = uniqueElements.count({(int)n0, (int)n1, 0}) ||
                               uniqueElements.count({(int)n1, (int)n0, 0});

            if (n0 == n1 || isDuplicate) {
              mesh->lines.erase(mesh->lines.begin() + i);
              currentMaterials.erase(currentMaterials.begin() + i);
              currentNormals.erase(currentNormals.begin() + i);
            } else {
              uniqueElements.insert({(int)n0, (int)n1, 0});
            }
          }
        }
        sharpCornerNodes[l].push_back(
            {static_cast<unsigned>(snapToNodeId), cornerPos});
      } else {
        // No snap - keep new corner and split any overlapping edge
        sharpCornerNodes[l].push_back({cornerNodeId, cornerPos});
        splitPreviousMaterialEdge(l, cornerNodeId, cornerPos, meshSizeBefore);
      }
    }
  }

public:
  ToMultiSurfaceMesh(double minNodeDistFactor = 0.05, double eps = 1e-12)
      : ToSurfaceMesh<NumericType, D>(minNodeDistFactor, eps) {}

  ToMultiSurfaceMesh(SmartPointer<lsDomainType> passedLevelSet,
                     SmartPointer<viennals::Mesh<NumericType>> passedMesh,
                     double minNodeDistFactor = 0.05, double eps = 1e-12)
      : ToSurfaceMesh<NumericType, D>(passedLevelSet, passedMesh,
                                      minNodeDistFactor, eps) {}

  ToMultiSurfaceMesh(
      std::vector<SmartPointer<lsDomainType>> const &passedLevelSets,
      SmartPointer<viennals::Mesh<NumericType>> passedMesh,
      double minNodeDistFactor = 0.05, double eps = 1e-12)
      : ToSurfaceMesh<NumericType, D>(minNodeDistFactor, eps) {
    levelSets = passedLevelSets;
    mesh = passedMesh;
  }

  ToMultiSurfaceMesh(SmartPointer<viennals::Mesh<NumericType>> passedMesh,
                     double minNodeDistFactor = 0.05, double eps = 1e-12)
      : ToSurfaceMesh<NumericType, D>(minNodeDistFactor, eps) {
    mesh = passedMesh;
  }

  void insertNextLevelSet(SmartPointer<lsDomainType> passedLevelSet) {
    levelSets.push_back(passedLevelSet);
  }

  void clearLevelSets() { levelSets.clear(); }

  void setMaterialMap(SmartPointer<MaterialMap> passedMaterialMap) {
    materialMap = passedMaterialMap;
  }

  void apply() override {
    if (levelSets.empty()) {
      Logger::getInstance()
          .addError("No level set was passed to ToMultiSurfaceMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      Logger::getInstance()
          .addError("No mesh was passed to ToMultiSurfaceMesh.")
          .print();
      return;
    }

    mesh->clear();
    currentGridDelta = levelSets.front()->getGrid().getGridDelta();
    for (unsigned i = 0; i < D; ++i) {
      mesh->minimumExtent[i] = std::numeric_limits<NumericType>::max();
      mesh->maximumExtent[i] = std::numeric_limits<NumericType>::lowest();
    }

    using nodeContainerType = std::map<viennahrle::Index<D>, unsigned>;

    nodeContainerType nodes[D];
    nodeContainerType faceNodes[D];
    nodeContainerType cornerNodes;

    nodeIdByBin.clear();
    uniqueElements.clear();
    currentNormals.clear();
    currentMaterials.clear();

    const bool useMaterialMap = materialMap != nullptr;
    const bool sharpCorners = generateSharpCorners;

    if (sharpCorners && D == 3) {
      Logger::getInstance()
          .addWarning("Sharp corner generation in 3D is experimental and may "
                      "produce suboptimal meshes. Use with caution.")
          .print();
    }

    // an iterator for each level set
    std::vector<viennahrle::ConstSparseCellIterator<hrleDomainType>> cellIts;
    for (const auto &ls : levelSets)
      cellIts.emplace_back(ls->getDomain());

    // Explicit storage for sharp corner nodes, keyed by material index
    std::unordered_map<unsigned, std::vector<SharpCornerNode>> sharpCornerNodes;

    for (unsigned l = 0; l < levelSets.size(); l++) {
      currentLevelSet = levelSets[l];
      if (useMaterialMap) {
        currentMaterialId = materialMap->getMaterialId(l);
      } else {
        currentMaterialId = static_cast<unsigned>(l);
      }

      normalVectorData = nullptr;
      if (sharpCorners) {
        CalculateNormalVectors<NumericType, D> normalCalculator(
            currentLevelSet);
        normalCalculator.setMethod(
            NormalCalculationMethodEnum::ONE_SIDED_MIN_MOD);
        normalCalculator.setMaxValue(std::numeric_limits<NumericType>::max());
        normalCalculator.apply();
        normalVectorData = currentLevelSet->getPointData().getVectorData(
            CalculateNormalVectors<NumericType, D>::normalVectorsLabel);
      }

      viennahrle::ConstSparseIterator<hrleDomainType> valueIt(
          currentLevelSet->getDomain());

      // iterate over all active surface points
      for (auto cellIt = cellIts[l]; !cellIt.isFinished(); cellIt.next()) {
        for (int u = 0; u < D; u++) {
          while (!nodes[u].empty() &&
                 nodes[u].begin()->first <
                     viennahrle::Index<D>(cellIt.getIndices()))
            nodes[u].erase(nodes[u].begin());

          if constexpr (D == 3) {
            while (!faceNodes[u].empty() &&
                   faceNodes[u].begin()->first <
                       viennahrle::Index<D>(cellIt.getIndices()))
              faceNodes[u].erase(faceNodes[u].begin());

            if (u == 0) {
              while (!cornerNodes.empty() &&
                     cornerNodes.begin()->first <
                         viennahrle::Index<D>(cellIt.getIndices()))
                cornerNodes.erase(cornerNodes.begin());
            }
          }
        }

        unsigned signs = 0;
        for (int i = 0; i < (1 << D); i++) {
          if (cellIt.getCorner(i).getValue() >= NumericType(0))
            signs |= (1 << i);
        }

        // all corners have the same sign, so no surface here
        if (signs == 0)
          continue;
        if (signs == (1 << (1 << D)) - 1)
          continue;

        // Check if this cell is at a material boundary with the immediate
        // material below
        bool atMaterialBoundary = false;
        unsigned touchingMaterial = 0; // Which material we're touching
        if (sharpCorners && l > 0) {
          touchingMaterial = l - 1;
          viennahrle::ConstSparseIterator<hrleDomainType> prevIt(
              levelSets[touchingMaterial]->getDomain());
          for (int i = 0; i < (1 << D); ++i) {
            hrleIndex cornerIdx =
                cellIt.getIndices() + viennahrle::BitMaskToIndex<D>(i);
            prevIt.goToIndices(cornerIdx);
            if (prevIt.getValue() <= 0) { // Inside previous material
              atMaterialBoundary = true;
              break;
            }
          }
        }

        // Attempt to generate sharp features if enabled
        bool perfectCornerFound = false;
        if (sharpCorners) {
          int countNeg = 0;
          int countPos = 0;
          int negMask = 0;
          int posMask = 0;
          for (int i = 0; i < (1 << D); ++i) {
            NumericType val = cellIt.getCorner(i).getValue();
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

          // Clear the base class's recent corner nodes before generation
          matSharpCornerNodes.clear();

          // Track mesh size before generating sharp corners
          size_t meshSizeBefore = 0;
          if constexpr (D == 2) {
            meshSizeBefore = mesh->lines.size();
            perfectCornerFound =
                this->generateSharpCorner2D(cellIt, countNeg, countPos, negMask,
                                            posMask, nodes, nullptr, valueIt);
          } else if constexpr (D == 3) {
            if (countNeg == 2 || countPos == 2) {
              perfectCornerFound = this->generateSharpEdge3D(
                  cellIt, countNeg, countPos, negMask, posMask, nodes,
                  faceNodes, nullptr, valueIt);
            } else if (countNeg == 1 || countPos == 1) {
              perfectCornerFound = this->generateSharpCorner3D(
                  cellIt, countNeg, countPos, negMask, posMask, nodes,
                  cornerNodes, faceNodes, nullptr, valueIt);
            } else if (countNeg == 3 || countPos == 3) {
              perfectCornerFound = this->generateSharpL3D(
                  cellIt, countNeg, countPos, negMask, posMask, nodes,
                  faceNodes, nullptr, valueIt);
            }
          }

          // If sharp corners were generated, check if they should snap to
          // existing corners or split edges of previous materials
          if (D == 2 && perfectCornerFound && !matSharpCornerNodes.empty()) {
            snapSharpCorners(l, meshSizeBefore, sharpCornerNodes);
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
                const int *Triangles =
                    lsInternal::MarchingCubes::polygonize3d(signs);

                for (; Triangles[0] != -1; Triangles += 3) {
                  std::vector<unsigned> face_edge_nodes;
                  for (int i = 0; i < 3; ++i) {
                    int edge = Triangles[i];
                    int c0 = this->corner0[edge];
                    int c1 = this->corner1[edge];
                    bool onFace =
                        (((c0 >> axis) & 1) == d) && (((c1 >> axis) & 1) == d);
                    if (onFace) {
                      face_edge_nodes.push_back(
                          this->getNode(cellIt, edge, nodes, nullptr));
                    }
                  }
                  if (face_edge_nodes.size() == 2) {
                    this->insertElement(std::array<unsigned, 3>{
                        face_edge_nodes[0], face_edge_nodes[1], faceNodeId});
                  }
                }
              }
            }
          }
        }

        // for each element
        const int *Triangles;
        if constexpr (D == 2) {
          Triangles = lsInternal::MarchingCubes::polygonize2d(signs);
        } else {
          Triangles = lsInternal::MarchingCubes::polygonize3d(signs);
        }

        for (; Triangles[0] != -1; Triangles += D) {
          std::array<unsigned, D> nodeNumbers;

          // for each node
          for (int n = 0; n < D; n++) {
            const int edge = Triangles[n];

            unsigned p0 = this->corner0[edge];

            // determine direction of edge
            unsigned dir = this->direction[edge];

            // look for existing surface node
            viennahrle::Index<D> d(cellIt.getIndices());
            auto p0B = viennahrle::BitMaskToIndex<D>(p0);
            d += p0B;

            auto nodeIt = nodes[dir].find(d);
            if (nodeIt != nodes[dir].end()) {
              nodeNumbers[n] = nodeIt->second;

              // Check if this cached node is a sharp corner that should be
              // inherited
              if (sharpCorners && l > 0 && atMaterialBoundary) {
                unsigned cachedNodeId = nodeIt->second;
                // Check if this is a corner from the material below
                auto it = sharpCornerNodes.find(touchingMaterial);
                if (it != sharpCornerNodes.end()) {
                  for (const auto &scn : it->second) {
                    if (scn.nodeId == cachedNodeId) {
                      // Only inherit if not already inherited
                      auto &lCorners = sharpCornerNodes[l];
                      bool alreadyInherited = false;
                      for (const auto &existing : lCorners) {
                        if (existing.nodeId == cachedNodeId) {
                          alreadyInherited = true;
                          break;
                        }
                      }
                      if (!alreadyInherited)
                        lCorners.push_back({cachedNodeId, scn.position});
                      break;
                    }
                  }
                }
              }
            } else {
              // Node does not exist yet - try snapping to a sharp corner
              // from the touching material at material boundaries
              int snapNodeId = -1;
              Vec3D<NumericType> snapPos{};
              if (sharpCorners && l > 0 && atMaterialBoundary) {
                auto tmIt = sharpCornerNodes.find(touchingMaterial);
                if (tmIt != sharpCornerNodes.end() && !tmIt->second.empty()) {
                  Vec3D<NumericType> nodePos =
                      this->computeNodePosition(cellIt, edge);

                  NumericType minDist2 = 0.1; // Threshold: 0.5 grid spacings
                  for (const auto &scn : tmIt->second) {
                    NumericType dist2 = 0;
                    for (int i = 0; i < D; ++i) {
                      NumericType dd = scn.position[i] - nodePos[i];
                      dist2 += dd * dd;
                    }
                    if (dist2 < minDist2) {
                      minDist2 = dist2;
                      snapNodeId = scn.nodeId;
                      snapPos = scn.position;
                    }
                  }
                }
              }

              if (snapNodeId >= 0) {
                nodeNumbers[n] = snapNodeId;
                nodes[dir][d] = snapNodeId; // Cache it
                // Inherit this corner so upper materials can find it
                sharpCornerNodes[l].push_back(
                    {static_cast<unsigned>(snapNodeId), snapPos});
              } else {
                nodeNumbers[n] = this->getNode(cellIt, edge, nodes, nullptr);
              }
            }
          }

          // Check for duplicates before inserting (standard Marching Cubes)
          bool isDuplicate = false;
          if constexpr (D == 2) {
            using I3 = typename ToSurfaceMesh<NumericType, D>::I3;
            if (uniqueElements.find({(int)nodeNumbers[1], (int)nodeNumbers[0],
                                     0}) != uniqueElements.end())
              isDuplicate = true;
          }

          if (!isDuplicate) {
            this->insertElement(nodeNumbers);
          }
        }
      }
    }

    this->scaleMesh();

    mesh->cellData.insertNextScalarData(currentMaterials, "MaterialIds");
    mesh->cellData.insertNextVectorData(currentNormals, "Normals");
    mesh->triangles.shrink_to_fit();
    mesh->nodes.shrink_to_fit();
  }
};

PRECOMPILE_PRECISION_DIMENSION(ToMultiSurfaceMesh)

} // namespace viennals