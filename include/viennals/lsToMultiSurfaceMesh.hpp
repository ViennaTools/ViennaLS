#pragma once

#include <lsMaterialMap.hpp>
#include <lsToSurfaceMesh.hpp>
#include <iostream>

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

  // Find the closest sharp corner node from materials below l, subject to a
  // grid-level layer thickness check.  Returns the node ID or -1 if none found.
  //
  // activeGridIdx : global grid index of the "active" (single-minority sign)
  //   corner of the cell in which material l's sharp corner was generated.
  // lsLAtActive   : LS value of material l at that grid corner.
  //
  // For conformal deposition of material l (thickness t) on material m:
  //   LS_l = LS_m - t  everywhere, so
  //   delta = LS_m(activeGridIdx) - LS_l(activeGridIdx) ≈ t.
  // Skipping snap when delta > layerThreshold prevents merging corners that
  // belong to physically distinct surfaces separated by a real layer.
  int findNearestLowerCorner(
      unsigned l, const Vec3D<NumericType> &pos,
      const std::unordered_map<unsigned, std::vector<SharpCornerNode>>
          &sharpCornerNodes,
      const hrleIndex &activeGridIdx, NumericType lsLAtActive) {
    int snapId = -1;
    NumericType minDist2 = NumericType(this->minNodeDistanceFactor);
    // Minimum layer thickness (in grid units) that suppresses snapping.
    const NumericType layerThreshold = NumericType(this->minNodeDistanceFactor);
    for (unsigned m = 0; m < l; ++m) {
      auto it = sharpCornerNodes.find(m);
      if (it == sharpCornerNodes.end())
        continue;
      // Evaluate LS_m at the active grid corner of material l's corner cell.
      // delta approximates the local layer thickness between material m and l.
      viennahrle::ConstSparseIterator<hrleDomainType> mIt(
          levelSets[m]->getDomain());
      mIt.goToIndices(activeGridIdx);
      const NumericType lsMAtActive = mIt.getValue();
      // Only apply the layer-thickness check when the active corner is inside
      // BOTH materials (LS_m < 0 and LS_l < 0). If LS_m > 0, the active corner
      // lies outside material m (e.g. inside a trench cavity), meaning the two
      // surfaces meet at a genuine junction — not a thin-layer overlap — and the
      // snap should be allowed regardless of delta.
      if (lsMAtActive < NumericType(0) &&
          lsMAtActive - lsLAtActive > layerThreshold)
        continue; // thin conformal layer present — do not snap
      for (const auto &scn : it->second) {
        NumericType dist2 = NumericType(0);
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

      NumericType lineLen2 = NumericType(0);
      NumericType dot = NumericType(0);
      for (int i = 0; i < D; ++i) {
        lineLen2 += lineVec[i] * lineVec[i];
        dot += (cornerPos[i] - pos0[i]) * lineVec[i];
      }

      if (lineLen2 < epsilon)
        continue;

      NumericType t = dot / lineLen2;
      if (t <= NumericType(0) || t >= NumericType(1))
        continue;

      // Distance from corner to closest point on line
      NumericType dist2 = NumericType(0);
      for (int i = 0; i < D; ++i) {
        NumericType d = cornerPos[i] - (pos0[i] + t * lineVec[i]);
        dist2 += d * d;
      }
      if (dist2 >= NumericType(this->minNodeDistanceFactor))
        continue;

      // Split the edge unconditionally to ensure the bottom material remains closed
      uniqueElements.erase({(int)node0, (int)node1, 0});
      uniqueElements.insert({(int)node0, (int)cornerNodeId, 0});
      uniqueElements.insert({(int)cornerNodeId, (int)node1, 0});

      line[1] = cornerNodeId;
      mesh->lines.push_back({cornerNodeId, node1});
      currentMaterials.push_back(mat);
      Vec3D<NumericType> normal = currentNormals[lineIdx];
      currentNormals.push_back(normal);

      // Remove any duplicate segments from the current material (l)
      for (size_t k = meshSizeBefore; k < mesh->lines.size(); ) {
        if (currentMaterials[k] != static_cast<NumericType>(l)) {
          ++k;
          continue;
        }

        const auto &rl = mesh->lines[k];
        bool isSeg1 = (rl[0] == node0 && rl[1] == cornerNodeId) || (rl[0] == cornerNodeId && rl[1] == node0);
        bool isSeg2 = (rl[0] == cornerNodeId && rl[1] == node1) || (rl[0] == node1 && rl[1] == cornerNodeId);

        if (isSeg1 || isSeg2) {
          mesh->lines.erase(mesh->lines.begin() + k);
          currentMaterials.erase(currentMaterials.begin() + k);
          currentNormals.erase(currentNormals.begin() + k);
        } else {
          ++k;
        }
      }

      break; // Corner can only be on one edge
    }
  }

  bool getSegmentIntersection(const Vec3D<NumericType> &p1,
                              const Vec3D<NumericType> &p2,
                              const Vec3D<NumericType> &p3,
                              const Vec3D<NumericType> &p4,
                              Vec3D<NumericType> &intersection,
                              NumericType &t) const {
    NumericType x1 = p1[0], y1 = p1[1];
    NumericType x2 = p2[0], y2 = p2[1];
    NumericType x3 = p3[0], y3 = p3[1];
    NumericType x4 = p4[0], y4 = p4[1];

    NumericType denom = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
    if (std::abs(denom) < epsilon)
      return false;

    NumericType ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denom;
    NumericType ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denom;

    if (ua > epsilon && ua < 1.0 - epsilon && ub > epsilon &&
        ub < 1.0 - epsilon) {
      intersection[0] = x1 + ua * (x2 - x1);
      intersection[1] = y1 + ua * (y2 - y1);
      if constexpr (D == 3)
        intersection[2] = 0;
      t = ua;
      return true;
    }
    return false;
  }

  // In 2D, check if the polymer edge at pIdx crosses BOTH corner edges of any
  // lower material corner node that lies in the same cell (cellIndices).
  // When both crossings are found simultaneously:
  //   - Split each mask corner edge at its intersection node.
  //   - Keep Part A (pNode0→intNode1) and Part C (intNode2→pNode1) of the
  //     polymer edge.
  //   - DISCARD Part B (intNode1→intNode2), which lies fictitiously inside the
  //     lower material.
  void handleCellCornerCrossings(
      size_t pIdx,
      const std::unordered_map<unsigned, std::vector<SharpCornerNode>>
          &sharpCornerNodes,
      unsigned l, const hrleIndex &cellIndices) {
    if constexpr (D != 2)
      return;

    unsigned pNode0 = mesh->lines[pIdx][0];
    unsigned pNode1 = mesh->lines[pIdx][1];
    Vec3D<NumericType> pp0 = mesh->nodes[pNode0];
    Vec3D<NumericType> pp1 = mesh->nodes[pNode1];

    for (unsigned m = 0; m < l; ++m) {
      auto cit = sharpCornerNodes.find(m);
      if (cit == sharpCornerNodes.end())
        continue;

      for (const auto &mscn : cit->second) {
        // Cell-local check: the lower material's corner must be in this cell
        bool inCell = true;
        for (int i = 0; i < D; ++i) {
          if (mscn.position[i] < NumericType(cellIndices[i]) - epsilon ||
              mscn.position[i] > NumericType(cellIndices[i] + 1) + epsilon) {
            inCell = false;
            break;
          }
        }
        if (!inCell)
          continue;

        // Find the two corner edges of material m incident to this corner node
        size_t edge1Idx = SIZE_MAX, edge2Idx = SIZE_MAX;
        for (size_t eIdx = 0; eIdx < pIdx; ++eIdx) {
          if (static_cast<unsigned>(currentMaterials[eIdx]) != m)
            continue;
          if (mesh->lines[eIdx][0] != mscn.nodeId &&
              mesh->lines[eIdx][1] != mscn.nodeId)
            continue;
          if (edge1Idx == SIZE_MAX)
            edge1Idx = eIdx;
          else if (edge2Idx == SIZE_MAX) {
            edge2Idx = eIdx;
            break;
          }
        }
        if (edge1Idx == SIZE_MAX || edge2Idx == SIZE_MAX)
          continue;

        // Check if the polymer edge crosses BOTH corner edges
        Vec3D<NumericType> int1{}, int2{};
        NumericType t1 = 0, t2 = 0;
        bool cross1 = getSegmentIntersection(pp0, pp1,
            mesh->nodes[mesh->lines[edge1Idx][0]],
            mesh->nodes[mesh->lines[edge1Idx][1]], int1, t1);
        bool cross2 = getSegmentIntersection(pp0, pp1,
            mesh->nodes[mesh->lines[edge2Idx][0]],
            mesh->nodes[mesh->lines[edge2Idx][1]], int2, t2);

        if (!cross1 || !cross2)
          continue;

        // Sort by parameter along polymer edge so t1 < t2
        if (t1 > t2) {
          std::swap(t1, t2);
          std::swap(int1, int2);
          std::swap(edge1Idx, edge2Idx);
        }

        // Insert intersection nodes
        unsigned intNode1 = this->insertNode(int1);
        unsigned intNode2 = this->insertNode(int2);

        // Save metadata before any array modifications
        NumericType matM1 = currentMaterials[edge1Idx];
        Vec3D<NumericType> normalM1 = currentNormals[edge1Idx];
        NumericType matM2 = currentMaterials[edge2Idx];
        Vec3D<NumericType> normalM2 = currentNormals[edge2Idx];
        NumericType matL = currentMaterials[pIdx];
        Vec3D<NumericType> normalL = currentNormals[pIdx];
        unsigned e1n0 = mesh->lines[edge1Idx][0];
        unsigned e1n1 = mesh->lines[edge1Idx][1];
        unsigned e2n0 = mesh->lines[edge2Idx][0];
        unsigned e2n1 = mesh->lines[edge2Idx][1];

        // Split mask edge1 at intNode1
        uniqueElements.erase({(int)e1n0, (int)e1n1, 0});
        mesh->lines[edge1Idx][1] = intNode1;
        uniqueElements.insert({(int)e1n0, (int)intNode1, 0});
        mesh->lines.push_back({intNode1, e1n1});
        currentMaterials.push_back(matM1);
        currentNormals.push_back(normalM1);
        uniqueElements.insert({(int)intNode1, (int)e1n1, 0});

        // Split mask edge2 at intNode2
        uniqueElements.erase({(int)e2n0, (int)e2n1, 0});
        mesh->lines[edge2Idx][1] = intNode2;
        uniqueElements.insert({(int)e2n0, (int)intNode2, 0});
        mesh->lines.push_back({intNode2, e2n1});
        currentMaterials.push_back(matM2);
        currentNormals.push_back(normalM2);
        uniqueElements.insert({(int)intNode2, (int)e2n1, 0});

        // Polymer edge:
        //   Part A (pNode0→intNode1) — kept
        //   Part B (intNode1→intNode2) — DISCARDED (fictitiously inside mask)
        //   Part C (intNode2→pNode1) — kept
        uniqueElements.erase({(int)pNode0, (int)pNode1, 0});
        mesh->lines[pIdx][1] = intNode1;
        uniqueElements.insert({(int)pNode0, (int)intNode1, 0});
        mesh->lines.push_back({intNode2, pNode1});
        currentMaterials.push_back(matL);
        currentNormals.push_back(normalL);
        uniqueElements.insert({(int)intNode2, (int)pNode1, 0});

        return; // Only one lower-material corner can be in this cell
      }
    }
  }

  // In 2D, check if the newly added edge (at the back of mesh->lines)
  // passes close to any sharp corner of a previous material.
  // If so, split the edge.
  void splitCurrentMaterialEdge(
      unsigned l, const std::unordered_map<unsigned, std::vector<SharpCornerNode>>
                      &sharpCornerNodes) {
    if constexpr (D != 2)
      return;
    if (l == 0)
      return;
    if (mesh->lines.empty())
      return;

    size_t lineIdx = mesh->lines.size() - 1;
    auto &line = mesh->lines[lineIdx];
    unsigned node0 = line[0];
    unsigned node1 = line[1];

    auto const &pos0 = mesh->nodes[node0];
    auto const &pos1 = mesh->nodes[node1];
    auto lineVec = pos1 - pos0;

    NumericType lineLen2 = NumericType(0);
    for (int i = 0; i < D; ++i)
      lineLen2 += lineVec[i] * lineVec[i];

    if (lineLen2 < epsilon)
      return;

    for (unsigned m = 0; m < l; ++m) {
      auto it = sharpCornerNodes.find(m);
      if (it == sharpCornerNodes.end())
        continue;

      for (const auto &scn : it->second) {
        unsigned cornerNodeId = scn.nodeId;
        const auto &cornerPos = scn.position;

        if (cornerNodeId == node0 || cornerNodeId == node1)
          continue;

        NumericType dot = NumericType(0);
        for (int i = 0; i < D; ++i)
          dot += (cornerPos[i] - pos0[i]) * lineVec[i];

        NumericType t = dot / lineLen2;
        if (t <= NumericType(0) || t >= NumericType(1))
          continue;

        NumericType dist2 = NumericType(0);
        for (int i = 0; i < D; ++i) {
          NumericType d = cornerPos[i] - (pos0[i] + t * lineVec[i]);
          dist2 += d * d;
        }

        if (dist2 < NumericType(this->minNodeDistanceFactor)) {
          bool keep0 = true;
          bool keep1 = true;

          // Check reverse edges (shared boundaries)
          if (uniqueElements.count({(int)cornerNodeId, (int)node0, 0})) keep0 = false;
          if (uniqueElements.count({(int)node1, (int)cornerNodeId, 0})) keep1 = false;

          // Check direct duplicates
          if (uniqueElements.count({(int)node0, (int)cornerNodeId, 0})) keep0 = false;
          if (uniqueElements.count({(int)cornerNodeId, (int)node1, 0})) keep1 = false;

          uniqueElements.erase({(int)node0, (int)node1, 0});

          if (keep0 && keep1) {
            line[1] = cornerNodeId;
            mesh->lines.push_back({cornerNodeId, node1});
            currentMaterials.push_back(currentMaterials[lineIdx]);
            currentNormals.push_back(currentNormals[lineIdx]);
            uniqueElements.insert({(int)node0, (int)cornerNodeId, 0});
            uniqueElements.insert({(int)cornerNodeId, (int)node1, 0});
          } else if (keep0) {
            line[1] = cornerNodeId;
            uniqueElements.insert({(int)node0, (int)cornerNodeId, 0});
          } else if (keep1) {
            line[0] = cornerNodeId;
            uniqueElements.insert({(int)cornerNodeId, (int)node1, 0});
          } else {
            mesh->lines.erase(mesh->lines.begin() + lineIdx);
            currentMaterials.erase(currentMaterials.begin() + lineIdx);
            currentNormals.erase(currentNormals.begin() + lineIdx);
          }
          return;
        }
      }
    }
  }

  // Handle snapping/splitting of sharp corners between materials.
  // activeGridIdx and lsLAtActive identify the active grid corner of the
  // current cell (see findNearestLowerCorner) for the layer thickness check.
  // cellIndices is the integer grid index of the current cell, used to restrict
  // the double-crossing check to lower material corners in the same cell.
  void
  snapSharpCorners(unsigned l, size_t meshSizeBefore,
                   std::unordered_map<unsigned, std::vector<SharpCornerNode>>
                       &sharpCornerNodes,
                   const hrleIndex &activeGridIdx, NumericType lsLAtActive,
                   const hrleIndex &cellIndices) {
    for (const auto &cornerPair : matSharpCornerNodes) {
      unsigned cornerNodeId = cornerPair.first;
      Vec3D<NumericType> cornerPos = cornerPair.second;

      int snapToNodeId =
          (l > 0) ? findNearestLowerCorner(l, cornerPos, sharpCornerNodes,
                                            activeGridIdx, lsLAtActive)
                  : -1;

      if (snapToNodeId >= 0) {
        // Snap: replace cornerNodeId with snapToNodeId in recent elements
        if (snapToNodeId != (int)cornerNodeId) {
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
        }
        sharpCornerNodes[l].push_back(
            {static_cast<unsigned>(snapToNodeId), cornerPos});
      } else {
        // No snap - keep new corner and split any overlapping edge
        sharpCornerNodes[l].push_back({cornerNodeId, cornerPos});
        splitPreviousMaterialEdge(l, cornerNodeId, cornerPos, meshSizeBefore);
        // For thin-layer case: the endpoints (nX, nY) of the polymer's sharp
        // corner edges may lie on previous material surface edges.
        // Collect all endpoints first, then apply all splits before crossing
        // checks. Interleaving the two would risk index-shifting: erases inside
        // splitPreviousMaterialEdge's duplicate-removal loop can invalidate
        // pIdx before handleCellCornerCrossings is called.
        size_t sizeNow = mesh->lines.size();
        std::vector<std::pair<unsigned, Vec3D<NumericType>>> endpointsToSplit;
        for (size_t pIdx = meshSizeBefore; pIdx < sizeNow; ++pIdx) {
          if (static_cast<unsigned>(currentMaterials[pIdx]) != l)
            continue;
          unsigned n0 = mesh->lines[pIdx][0];
          unsigned n1 = mesh->lines[pIdx][1];
          if (n0 != cornerNodeId)
            endpointsToSplit.push_back({n0, mesh->nodes[n0]});
          if (n1 != cornerNodeId)
            endpointsToSplit.push_back({n1, mesh->nodes[n1]});
        }
        for (auto &[nodeId, pos] : endpointsToSplit)
          splitPreviousMaterialEdge(l, nodeId, pos, meshSizeBefore);
        // Run crossing checks after all lower-material splits are complete.
        sizeNow = mesh->lines.size();
        for (size_t pIdx = meshSizeBefore; pIdx < sizeNow; ++pIdx) {
          if (static_cast<unsigned>(currentMaterials[pIdx]) != l)
            continue;
          handleCellCornerCrossings(pIdx, sharpCornerNodes, l, cellIndices);
        }
      }
    }
  }

public:
  ToMultiSurfaceMesh(double minNodeDistFactor = 0.01, double eps = 1e-12)
      : ToSurfaceMesh<NumericType, D>(minNodeDistFactor, eps) {}

  ToMultiSurfaceMesh(SmartPointer<lsDomainType> passedLevelSet,
                     SmartPointer<viennals::Mesh<NumericType>> passedMesh,
                     double minNodeDistFactor = 0.01, double eps = 1e-12)
      : ToSurfaceMesh<NumericType, D>(passedLevelSet, passedMesh,
                                      minNodeDistFactor, eps) {}

  ToMultiSurfaceMesh(
      std::vector<SmartPointer<lsDomainType>> const &passedLevelSets,
      SmartPointer<viennals::Mesh<NumericType>> passedMesh,
      double minNodeDistFactor = 0.01, double eps = 1e-12)
      : ToSurfaceMesh<NumericType, D>(minNodeDistFactor, eps) {
    levelSets = passedLevelSets;
    mesh = passedMesh;
  }

  ToMultiSurfaceMesh(SmartPointer<viennals::Mesh<NumericType>> passedMesh,
                     double minNodeDistFactor = 0.01, double eps = 1e-12)
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
            if (prevIt.getValue() <= NumericType(0)) { // Inside previous material
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
              if (val >= NumericType(0)) {
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
          // existing corners or split edges of previous materials.
          // Identify the "active" (single-minority sign) grid corner so the
          // layer thickness check inside snapSharpCorners can use grid-level LS
          // values instead of unreliable sub-grid interpolation.
          if (D == 2 && perfectCornerFound && !matSharpCornerNodes.empty()) {
            int activeBit = -1;
            if (countNeg == 1) {
              for (int i = 0; i < (1 << D); ++i)
                if ((negMask >> i) & 1) { activeBit = i; break; }
            } else if (countPos == 1) {
              for (int i = 0; i < (1 << D); ++i)
                if ((posMask >> i) & 1) { activeBit = i; break; }
            }
            hrleIndex activeGridIdx = cellIt.getIndices();
            NumericType lsLAtActive = NumericType(0);
            if (activeBit >= 0) {
              activeGridIdx += viennahrle::BitMaskToIndex<D>(activeBit);
              lsLAtActive = cellIt.getCorner(activeBit).getValue();
            }
            snapSharpCorners(l, meshSizeBefore, sharpCornerNodes,
                             activeGridIdx, lsLAtActive, cellIt.getIndices());
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

                  NumericType minDist2 = NumericType(this->minNodeDistanceFactor);
                  for (const auto &scn : tmIt->second) {
                    NumericType dist2 = NumericType(0);
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
            if (sharpCorners && l > 0) {
              size_t mcIdx = mesh->lines.size() - 1;
              splitCurrentMaterialEdge(l, sharpCornerNodes);
              // Only call if the edge still exists at mcIdx; splitCurrentMaterialEdge
              // may have erased it (keep-neither case) or appended a new edge (keep-both),
              // but mcIdx always refers to the original inserted edge.
              if (mesh->lines.size() > mcIdx &&
                  static_cast<unsigned>(currentMaterials[mcIdx]) == l)
                handleCellCornerCrossings(mcIdx, sharpCornerNodes,
                                          l, cellIt.getIndices());
            }
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
