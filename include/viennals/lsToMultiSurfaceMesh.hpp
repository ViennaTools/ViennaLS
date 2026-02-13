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
  using ConstSparseIterator = viennahrle::ConstSparseIterator<hrleDomainType>;

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
    unsigned materialId;
    Vec3D<NumericType> position;
  };

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

    typedef std::map<viennahrle::Index<D>, unsigned> nodeContainerType;

    nodeContainerType nodes[D];
    nodeContainerType faceNodes[D];
    nodeContainerType cornerNodes;

    nodeIdByBin.clear();
    uniqueElements.clear();
    currentNormals.clear();
    currentMaterials.clear();

    typename nodeContainerType::iterator nodeIt;

    const bool useMaterialMap = materialMap != nullptr;
    const bool sharpCorners = generateSharpCorners;

    auto elementI3 = [&](const std::array<unsigned, D> &element) ->
        typename ToSurfaceMesh<NumericType, D>::I3 {
          if constexpr (D == 2)
            return {(int)element[0], (int)element[1], 0};
          else
            return {(int)element[0], (int)element[1], (int)element[2]};
        };

    // an iterator for each level set
    std::vector<viennahrle::ConstSparseCellIterator<hrleDomainType>> cellIts;
    for (const auto &ls : levelSets)
      cellIts.emplace_back(ls->getDomain());

    // Explicit storage for sharp corner nodes from each material
    std::vector<SharpCornerNode> sharpCornerNodes;

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
        // static unsigned boundaryCount = 0;
        if (sharpCorners && l > 0) {
          // Only check the material immediately below (l-1)
          int m = l - 1;
          viennahrle::ConstSparseIterator<hrleDomainType> prevIt(
              levelSets[m]->getDomain());
          for (int i = 0; i < (1 << D); ++i) {
            hrleIndex cornerIdx =
                cellIt.getIndices() + viennahrle::BitMaskToIndex<D>(i);
            prevIt.goToIndices(cornerIdx);
            if (prevIt.getValue() <= 0) { // Inside previous material
              atMaterialBoundary = true;
              touchingMaterial = m;
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
          // existing corners
          if (perfectCornerFound && !matSharpCornerNodes.empty()) {

            // Iterate over all generated sharp corners and check for snapping
            for (const auto &cornerPair : matSharpCornerNodes) {
              unsigned cornerNodeId = cornerPair.first;
              Vec3D<NumericType> cornerPos = cornerPair.second;
              int snapToNodeId = -1;

              if (l > 0) {
                NumericType minDist2 = 0.1;
                for (const auto &scn : sharpCornerNodes) {
                  if (scn.materialId < l) {
                    NumericType dist2 = 0;
                    for (int i = 0; i < D; ++i) {
                      NumericType d = scn.position[i] - cornerPos[i];
                      dist2 += d * d;
                    }
                    if (dist2 < minDist2) {
                      minDist2 = dist2;
                      snapToNodeId = scn.nodeId;
                    }
                  }
                }
              }

              if (snapToNodeId >= 0) {
                // Snap to existing sharp corner
                if constexpr (D == 2) {
                  for (size_t i = mesh->lines.size(); i-- > meshSizeBefore;) {
                    for (int j = 0; j < 2; ++j) {
                      if (mesh->lines[i][j] == cornerNodeId)
                        mesh->lines[i][j] = snapToNodeId;
                    }

                    // Check for duplicates or degenerate lines
                    unsigned n0 = mesh->lines[i][0];
                    unsigned n1 = mesh->lines[i][1];

                    bool isDuplicate = false;
                    using I3 = typename ToSurfaceMesh<NumericType, D>::I3;
                    if (uniqueElements.find({(int)n0, (int)n1, 0}) !=
                        uniqueElements.end())
                      isDuplicate = true;
                    if (uniqueElements.find({(int)n1, (int)n0, 0}) !=
                        uniqueElements.end())
                      isDuplicate = true;

                    if (n0 == n1 || isDuplicate) {
                      mesh->lines.erase(mesh->lines.begin() + i);
                      currentMaterials.erase(currentMaterials.begin() + i);
                      currentNormals.erase(currentNormals.begin() + i);
                    } else {
                      // Add to uniqueElements to prevent future duplicates
                      uniqueElements.insert({(int)n0, (int)n1, 0});
                    }
                  }
                }

                // Add to sharpCornerNodes for current material
                SharpCornerNode scn;
                scn.nodeId = snapToNodeId;
                scn.materialId = l;
                scn.position = cornerPos;
                sharpCornerNodes.push_back(scn);

              } else {
                // No snap - keep new corner
                SharpCornerNode scn;
                scn.nodeId = cornerNodeId;
                scn.materialId = l;
                scn.position = cornerPos;
                sharpCornerNodes.push_back(scn);

                // Check if this corner lies on an edge of a previous material's
                // mesh If so, split that edge to include this corner node
                if constexpr (D == 2) {
                  if (l > 0) {
                    // Find which material owns each line segment
                    for (size_t lineIdx = 0; lineIdx < mesh->lines.size();
                         ++lineIdx) {
                      unsigned mat = currentMaterials[lineIdx];

                      // Only check lines from previous materials
                      if (mat >= l)
                        continue;

                      auto &line = mesh->lines[lineIdx];
                      unsigned node0 = line[0];
                      unsigned node1 = line[1];

                      // Get positions of the line endpoints
                      Vec3D<NumericType> pos0 = mesh->nodes[node0];
                      Vec3D<NumericType> pos1 = mesh->nodes[node1];

                      // Check if corner lies on this line segment
                      // Calculate distance from corner to line
                      Vec3D<NumericType> lineVec = pos1 - pos0;
                      Vec3D<NumericType> cornerVec = cornerPos - pos0;

                      NumericType lineLen2 = 0;
                      NumericType dot = 0;
                      for (int i = 0; i < D; ++i) {
                        lineLen2 += lineVec[i] * lineVec[i];
                        dot += cornerVec[i] * lineVec[i];
                      }

                      if (lineLen2 < epsilon)
                        continue; // Degenerate line

                      NumericType t = dot / lineLen2;

                      // Check if corner projects onto the line segment (0 < t <
                      // 1)
                      if (t <= 0 || t >= 1)
                        continue;

                      // Calculate closest point on line
                      Vec3D<NumericType> closestPt;
                      for (int i = 0; i < D; ++i) {
                        closestPt[i] = pos0[i] + t * lineVec[i];
                      }

                      // Check distance from corner to line
                      NumericType dist2 = 0;
                      for (int i = 0; i < D; ++i) {
                        NumericType d = cornerPos[i] - closestPt[i];
                        dist2 += d * d;
                      }

                      // If corner is very close to the line, split it
                      if (dist2 < 0.01) { // Same threshold as corner snapping
                        // Check for duplicates in Material 3's recent lines
                        int seg1Idx = -1;
                        int seg2Idx = -1;
                        for (size_t k = meshSizeBefore; k < mesh->lines.size();
                             ++k) {
                          const auto &l3 = mesh->lines[k];
                          if ((l3[0] == node0 && l3[1] == cornerNodeId) ||
                              (l3[0] == cornerNodeId && l3[1] == node0)) {
                            seg1Idx = k;
                          }
                          if ((l3[0] == cornerNodeId && l3[1] == node1) ||
                              (l3[0] == node1 && l3[1] == cornerNodeId)) {
                            seg2Idx = k;
                          }
                        }

                        if (seg1Idx != -1 && seg2Idx != -1) {
                          currentMaterials[seg1Idx] = mat;
                          mesh->lines[seg1Idx] = {node0, cornerNodeId};
                          currentNormals[seg1Idx] = currentNormals[lineIdx];

                          currentMaterials[seg2Idx] = mat;
                          mesh->lines[seg2Idx] = {cornerNodeId, node1};
                          currentNormals[seg2Idx] = currentNormals[lineIdx];

                          line[1] = line[0]; // Degenerate original line
                        } else if (seg1Idx != -1) {
                          currentMaterials[seg1Idx] = mat;
                          mesh->lines[seg1Idx] = {node0, cornerNodeId};
                          currentNormals[seg1Idx] = currentNormals[lineIdx];
                          line[0] = cornerNodeId;
                          line[1] = node1;
                        } else if (seg2Idx != -1) {
                          currentMaterials[seg2Idx] = mat;
                          mesh->lines[seg2Idx] = {cornerNodeId, node1};
                          currentNormals[seg2Idx] = currentNormals[lineIdx];
                          line[1] = cornerNodeId;
                        } else {
                          // Replace this line with two lines that include the
                          // corner
                          line[1] =
                              cornerNodeId; // First segment: node0 -> corner

                          // Add second segment: corner -> node1
                          std::array<unsigned, 2> newLine = {cornerNodeId,
                                                             node1};
                          mesh->lines.push_back(newLine);
                          currentMaterials.push_back(
                              mat); // Same material as original
                          currentNormals.push_back(
                              currentNormals[lineIdx]); // Same normal
                        }

                        break; // Corner can only be on one edge
                      }
                    }
                  }
                }
              }
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
            unsigned p1 = this->corner1[edge];

            // determine direction of edge
            unsigned dir = this->direction[edge];

            // look for existing surface node
            viennahrle::Index<D> d(cellIt.getIndices());
            auto p0B = viennahrle::BitMaskToIndex<D>(p0);
            d += p0B;

            nodeIt = nodes[dir].find(d);
            if (nodeIt != nodes[dir].end()) {
              nodeNumbers[n] = nodeIt->second;

              // Check if this cached node is a sharp corner that should be
              // inherited
              if (sharpCorners && l > 0 && atMaterialBoundary) {
                unsigned cachedNodeId = nodeIt->second;
                // Check if this is a corner from the material below
                bool foundCorner = false;
                for (const auto &scn : sharpCornerNodes) {
                  if (scn.nodeId == cachedNodeId &&
                      scn.materialId == touchingMaterial) {
                    // This is a corner from the lower material - inherit it
                    SharpCornerNode inheritedCorner;
                    inheritedCorner.nodeId = cachedNodeId;
                    inheritedCorner.materialId = l;
                    inheritedCorner.position = scn.position;
                    sharpCornerNodes.push_back(inheritedCorner);
                    foundCorner = true;
                    break;
                  }
                }
              }
            } else {
              // if node does not exist yet
              // For materials after the first with sharp corners, check for
              // nearby sharp corner nodes from previous materials at material
              // boundaries
              if (sharpCorners && l > 0 && atMaterialBoundary &&
                  !sharpCornerNodes.empty()) {
                // Calculate where this node would be
                unsigned p0 = this->corner0[edge];
                unsigned p1 = this->corner1[edge];
                Vec3D<NumericType> nodePos{};
                for (int z = 0; z < D; z++) {
                  if (z != dir) {
                    nodePos[z] = static_cast<NumericType>(
                        cellIt.getIndices(z) +
                        viennahrle::BitMaskToIndex<D>(p0)[z]);
                  } else {
                    NumericType d0 = cellIt.getCorner(p0).getValue();
                    NumericType d1 = cellIt.getCorner(p1).getValue();
                    if (d0 == -d1) {
                      nodePos[z] =
                          static_cast<NumericType>(cellIt.getIndices(z)) + 0.5;
                    } else {
                      nodePos[z] = static_cast<NumericType>(
                          cellIt.getIndices(z) + (d0 / (d0 - d1)));
                    }
                  }
                }

                // Search for nearby sharp corner from the specific material
                // we're touching
                int snapNodeId = -1;
                NumericType minDist2 = 0.1; // Threshold: 0.5 grid spacings
                for (const auto &scn : sharpCornerNodes) {
                  // Only snap to corners from the material we're directly
                  // touching
                  if (scn.materialId != touchingMaterial)
                    continue;

                  // Check distance in all dimensions
                  NumericType dist2 = 0;
                  for (int i = 0; i < D; ++i) {
                    NumericType d = scn.position[i] - nodePos[i];
                    dist2 += d * d;
                  }

                  if (dist2 < minDist2) {
                    minDist2 = dist2;
                    snapNodeId = scn.nodeId;
                  }
                }

                if (snapNodeId >= 0) {
                  // Snap to existing sharp corner
                  nodeNumbers[n] = snapNodeId;
                  nodes[dir][d] = snapNodeId; // Cache it

                  // Find the corner we're snapping to and add it to current
                  // material's corners This allows upper materials to find it
                  // when searching for this material's corners
                  for (const auto &scn : sharpCornerNodes) {
                    if (scn.nodeId == snapNodeId) {
                      SharpCornerNode inheritedCorner;
                      inheritedCorner.nodeId = snapNodeId;
                      inheritedCorner.materialId =
                          l; // Current material inherits this corner
                      inheritedCorner.position = scn.position;
                      sharpCornerNodes.push_back(inheritedCorner);
                      break;
                    }
                  }

                } else {
                  // Create new node normally
                  nodeNumbers[n] = this->getNode(cellIt, edge, nodes, nullptr);
                }
              } else {
                // First material or no sharp corners - create normally
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

    // // DEBUG: Check specific nodes
    // std::vector<unsigned> debugNodes = {179, 210};
    // for (unsigned nId : debugNodes) {
    //   if (nId < mesh->nodes.size()) {
    //     std::cout << "DEBUG: Node " << nId << " at " << mesh->nodes[nId]
    //               << " is connected to elements:" << std::endl;
    //     if constexpr (D == 2) {
    //       for (size_t i = 0; i < mesh->lines.size(); ++i) {
    //         if (mesh->lines[i][0] == nId || mesh->lines[i][1] == nId) {
    //           std::cout << "  Line " << i << " [" << mesh->lines[i][0] << ",
    //           "
    //                     << mesh->lines[i][1]
    //                     << "] Material: " << currentMaterials[i] <<
    //                     std::endl;
    //         }
    //       }
    //     }
    //   }
    // }

    mesh->cellData.insertNextScalarData(currentMaterials, "MaterialIds");
    mesh->cellData.insertNextVectorData(currentNormals, "Normals");
    mesh->triangles.shrink_to_fit();
    mesh->nodes.shrink_to_fit();
  }
};

PRECOMPILE_PRECISION_DIMENSION(ToMultiSurfaceMesh)

} // namespace viennals