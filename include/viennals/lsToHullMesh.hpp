#pragma once

#include <hrleDenseIterator.hpp>
#include <map>
#include <tuple>

#include <lsDomain.hpp>
#include <lsMaterialMap.hpp>
#include <lsPreCompileMacros.hpp>
#include <lsToMultiSurfaceMesh.hpp>
#include <unordered_map>
#include <utility>

namespace viennals {

using namespace viennacore;

/// Extracts a hull (surface outline) mesh from a stack of level sets using
/// ToMultiSurfaceMesh with closed boundary caps. Material IDs are assigned per
/// cell. setSharpCorners(true/false) controls whether sharp corner generation
/// is enabled during surface extraction.
template <class T, int D> class ToHullMesh {
  typedef typename Domain<T, D>::DomainType hrleDomainType;
  using LevelSetsType = std::vector<SmartPointer<Domain<T, D>>>;
  LevelSetsType levelSets;
  SmartPointer<Mesh<T>> multiMesh = nullptr;
  SmartPointer<MaterialMap> materialMap = nullptr;
  bool generateSharpCorners = false;

  /// Generate hull mesh using ToMultiSurfaceMesh and boundary caps.
  /// generateSharpCorners is forwarded to the converter. Works for both
  /// 2D (line segments) and 3D (triangles).
  void generateHull() {
    if (levelSets.empty()) {
      VIENNACORE_LOG_WARNING(
          "ToHullMesh: No level sets provided! Hull mesh will be empty.");
      return;
    }

    if (multiMesh == nullptr) {
      VIENNACORE_LOG_WARNING("ToHullMesh: No mesh provided!");
      return;
    }

    auto &grid = levelSets[0]->getGrid();
    double gridDelta = grid.getGridDelta();

    int totalMinimum[D];
    int totalMaximum[D];
    for (int i = 0; i < D; ++i) {
      totalMinimum[i] = std::numeric_limits<int>::max();
      totalMaximum[i] = std::numeric_limits<int>::lowest();
    }

    for (auto &ls : levelSets) {
      if (ls->getNumberOfPoints() == 0)
        continue;
      auto &dom = ls->getDomain();
      auto &grd = ls->getGrid();
      for (int i = 0; i < D; ++i) {
        int minP = (grd.isNegBoundaryInfinite(i)) ? dom.getMinRunBreak(i) - 1
                                                  : grd.getMinBounds(i);
        int maxP = (grd.isPosBoundaryInfinite(i)) ? dom.getMaxRunBreak(i) + 1
                                                  : grd.getMaxBounds(i);
        if (minP < totalMinimum[i])
          totalMinimum[i] = minP;
        if (maxP > totalMaximum[i])
          totalMaximum[i] = maxP;
      }
    }

    // Generate multi-surface mesh with sharp corners
    ToMultiSurfaceMesh<T, D> converter;
    for (auto &ls : levelSets)
      converter.insertNextLevelSet(ls);
    if (materialMap)
      converter.setMaterialMap(materialMap);
    converter.setMesh(multiMesh);
    converter.setSharpCorners(generateSharpCorners);
    converter.apply();

    // remove normals from cell data
    multiMesh->cellData.eraseVectorData(0);

    // Add closed boundary caps
    if constexpr (D == 2) {
      auto capBoundary = [&](int boundaryAxis, double fixedValue,
                             double rangeMin, double rangeMax) {
        int varyingAxis = 1 - boundaryAxis;
        std::vector<std::pair<double, unsigned>> boundaryNodes;
        double tolerance = 1e-6 * gridDelta;

        for (unsigned i = 0; i < multiMesh->nodes.size(); ++i) {
          if (std::abs(multiMesh->nodes[i][boundaryAxis] - fixedValue) <
              tolerance) {
            if (multiMesh->nodes[i][varyingAxis] >= rangeMin - tolerance &&
                multiMesh->nodes[i][varyingAxis] <= rangeMax + tolerance) {
              boundaryNodes.push_back({multiMesh->nodes[i][varyingAxis], i});
            }
          }
        }

        auto addCorner = [&](double val) {
          for (const auto &bn : boundaryNodes)
            if (std::abs(bn.first - val) < tolerance)
              return;

          if (levelSets.back()->getNumberOfPoints() == 0)
            return;

          viennahrle::Index<D> idx;
          idx[boundaryAxis] = std::round(fixedValue / gridDelta);
          idx[varyingAxis] = std::round(val / gridDelta);

          auto &grid = levelSets.back()->getGrid();
          if (grid.isOutsideOfDomain(idx))
            idx = grid.globalIndices2LocalIndices(idx);

          // Clamp to grid bounds to ensure iterator safety
          auto minGrid = grid.getMinGridPoint();
          auto maxGrid = grid.getMaxGridPoint();
          for (int i = 0; i < D; ++i) {
            if (idx[i] < minGrid[i])
              idx[i] = minGrid[i];
            if (idx[i] > maxGrid[i])
              idx[i] = maxGrid[i];
          }

          viennahrle::ConstDenseIterator<hrleDomainType> it(
              levelSets.back()->getDomain());
          it.goToIndices(idx);
          if (it.getValue() > 0.)
            return;

          Vec3D<T> pos;
          pos[boundaryAxis] = fixedValue;
          pos[varyingAxis] = val;
          pos[2] = 0;
          unsigned id = multiMesh->insertNextNode(pos);
          boundaryNodes.push_back({val, id});
        };
        addCorner(rangeMin);
        addCorner(rangeMax);

        std::sort(boundaryNodes.begin(), boundaryNodes.end());

        if (boundaryNodes.size() < 2)
          return;

        auto matIds = multiMesh->cellData.getScalarData("MaterialIds");
        for (size_t i = 0; i < boundaryNodes.size() - 1; ++i) {
          double mid =
              (boundaryNodes[i].first + boundaryNodes[i + 1].first) * 0.5;

          int matId = -1;
          int boundaryCandidate = -1;
          viennahrle::Index<D> queryIdx;
          queryIdx[boundaryAxis] = std::round(fixedValue / gridDelta);
          queryIdx[varyingAxis] = std::round(mid / gridDelta);

          for (unsigned l = 0; l < levelSets.size(); ++l) {
            viennahrle::ConstDenseIterator<hrleDomainType> it(
                levelSets[l]->getDomain());
            it.goToIndices(queryIdx);
            double val = it.getValue();
            if (val < -1e-6 * gridDelta) {
              matId = (materialMap) ? materialMap->getMaterialId(l) : (int)l;
              break;
            }
            if (val <= 1e-6 * gridDelta && boundaryCandidate == -1) {
              boundaryCandidate =
                  (materialMap) ? materialMap->getMaterialId(l) : (int)l;
            }
          }
          if (matId == -1)
            matId = boundaryCandidate;

          if (matId != -1) {
            multiMesh->insertNextLine(
                {static_cast<unsigned>(boundaryNodes[i].second),
                 static_cast<unsigned>(boundaryNodes[i + 1].second)});
            if (matIds)
              matIds->push_back(matId);
          }
        }
      };

      double minX = totalMinimum[0] * gridDelta;
      double maxX = totalMaximum[0] * gridDelta;
      double minY = totalMinimum[1] * gridDelta;
      double maxY = totalMaximum[1] * gridDelta;

      capBoundary(0, minX, minY, maxY);
      capBoundary(0, maxX, minY, maxY);
      capBoundary(1, minY, minX, maxX);
      capBoundary(1, maxY, minX, maxX);

    } else {
      // 3D implementation: Marching Squares on boundary faces

      // Map to find existing nodes on grid edges/corners
      // EdgeKey: (dir, i, j, k) where dir is edge direction (0=x, 1=y, 2=z)
      using EdgeKey = std::tuple<int, int, int, int>;
      std::map<EdgeKey, unsigned> edgeToNode;
      std::map<std::tuple<int, int, int>, unsigned> existingGridNodes;
      double invGridDelta = 1.0 / gridDelta;

      for (unsigned i = 0; i < multiMesh->nodes.size(); ++i) {
        const auto &pos = multiMesh->nodes[i];
        int intCoords[3];
        bool isInt[3];
        int intCount = 0;
        for (int d = 0; d < 3; ++d) {
          double val = pos[d] * invGridDelta;
          intCoords[d] = std::round(val);
          isInt[d] = std::abs(val - intCoords[d]) < 1e-4;
          if (isInt[d])
            intCount++;
        }

        if (intCount == 3) {
          existingGridNodes[{intCoords[0], intCoords[1], intCoords[2]}] = i;
        } else if (intCount == 2) {
          int dir = -1;
          for (int d = 0; d < 3; ++d)
            if (!isInt[d])
              dir = d;
          int start[3];
          for (int d = 0; d < 3; ++d) {
            if (d == dir)
              start[d] = std::floor(pos[d] * invGridDelta);
            else
              start[d] = intCoords[d];
          }
          edgeToNode[{dir, start[0], start[1], start[2]}] = i;
        }
      }

      std::map<std::tuple<int, int, int>, unsigned> gridNodeMap;
      auto getGridNode = [&](int i, int j, int k) {
        std::tuple<int, int, int> key = {i, j, k};
        auto it = gridNodeMap.find(key);
        if (it != gridNodeMap.end())
          return it->second;
        auto it2 = existingGridNodes.find(key);
        if (it2 != existingGridNodes.end())
          return it2->second;

        Vec3D<T> pos;
        pos[0] = i * gridDelta;
        pos[1] = j * gridDelta;
        pos[2] = k * gridDelta;
        unsigned id = multiMesh->insertNextNode(pos);
        gridNodeMap[key] = id;
        return id;
      };

      // Marching Squares Cases (Triangles)
      // 0-3: Corners (BL, BR, TR, TL), 4-7: Edges (B, R, T, L)
      const std::vector<std::vector<int>> msTris = {{},
                                                    {0, 4, 7},
                                                    {1, 5, 4},
                                                    {0, 1, 5, 0, 5, 7},
                                                    {2, 6, 5},
                                                    {0, 4, 7, 2, 6, 5},
                                                    {1, 2, 6, 1, 6, 4},
                                                    {0, 1, 2, 0, 2, 6, 0, 6, 7},
                                                    {3, 7, 6},
                                                    {0, 4, 6, 0, 6, 3},
                                                    {1, 5, 4, 3, 7, 6},
                                                    {0, 1, 5, 0, 5, 6, 0, 6, 3},
                                                    {3, 2, 5, 3, 5, 7},
                                                    {0, 4, 5, 0, 5, 2, 0, 2, 3},
                                                    {1, 2, 3, 1, 3, 7, 1, 7, 4},
                                                    {0, 1, 2, 0, 2, 3}};

      for (int d = 0; d < D; ++d) {
        for (int side = 0; side < 2; ++side) {
          int boundaryCoord = (side == 0) ? totalMinimum[d] : totalMaximum[d];
          int u = (d + 1) % D;
          int v = (d + 2) % D;
          if (D == 3 && u > v)
            std::swap(u, v);

          int uMin = totalMinimum[u];
          int uMax = totalMaximum[u];
          int vMin = (D == 3) ? totalMinimum[v] : 0;
          int vMax = (D == 3) ? totalMaximum[v] : 0;

          bool flip = (side == 0);

          for (int i = uMin; i < uMax; ++i) {
            for (int j = vMin; j < ((D == 3) ? vMax : 1); ++j) {
              int cI[4] = {i, i + 1, i + 1, i};
              int cJ[4] = {j, j, j + 1, j + 1};

              int mask = 0;
              int matId = -1;

              for (int k = 0; k < 4; ++k) {
                viennahrle::Index<D> idx;
                idx[d] = boundaryCoord;
                idx[u] = cI[k];
                idx[v] = cJ[k];

                for (unsigned l = 0; l < levelSets.size(); ++l) {
                  viennahrle::ConstDenseIterator<hrleDomainType> it(
                      levelSets[l]->getDomain());
                  it.goToIndices(idx);
                  if (it.getValue() <= 0) {
                    mask |= (1 << k);
                    if (matId == -1)
                      matId = (materialMap) ? materialMap->getMaterialId(l)
                                            : (int)l;
                    break;
                  }
                }
              }

              if (mask == 0)
                continue;
              if (matId == -1)
                matId = 0;

              auto &tris = msTris[mask];
              for (size_t t = 0; t < tris.size(); t += 3) {
                unsigned nodes[3];
                for (int k = 0; k < 3; ++k) {
                  int pt = tris[t + k];
                  if (pt < 4) {
                    int coords[3];
                    coords[d] = boundaryCoord;
                    coords[u] = cI[pt];
                    coords[v] = cJ[pt];
                    nodes[k] = getGridNode(coords[0], coords[1], coords[2]);
                  } else {
                    int edgeDir, ei;
                    ei = boundaryCoord;
                    int s[3];
                    s[d] = ei;

                    if (pt == 4) { // Bottom
                      edgeDir = u;
                      s[u] = i;
                      s[v] = j;
                    } else if (pt == 5) { // Right
                      edgeDir = v;
                      s[u] = i + 1;
                      s[v] = j;
                    } else if (pt == 6) { // Top
                      edgeDir = u;
                      s[u] = i;
                      s[v] = j + 1;
                    } else { // Left
                      edgeDir = v;
                      s[u] = i;
                      s[v] = j;
                    }
                    nodes[k] = edgeToNode[{edgeDir, s[0], s[1], s[2]}];

                    if (nodes[k] == 0 &&
                        edgeToNode.find({0, 0, 0, 0}) == edgeToNode.end()) {
                      // Fallback: create midpoint if node missing
                      Vec3D<T> pos;
                      pos[d] = ei * gridDelta;
                      if (pt == 4) {
                        pos[u] = (i + 0.5) * gridDelta;
                        pos[v] = j * gridDelta;
                      } else if (pt == 5) {
                        pos[u] = (i + 1) * gridDelta;
                        pos[v] = (j + 0.5) * gridDelta;
                      } else if (pt == 6) {
                        pos[u] = (i + 0.5) * gridDelta;
                        pos[v] = (j + 1) * gridDelta;
                      } else {
                        pos[u] = i * gridDelta;
                        pos[v] = (j + 0.5) * gridDelta;
                      }
                      nodes[k] = multiMesh->insertNextNode(pos);
                    }
                  }
                }

                if (flip)
                  std::swap(nodes[1], nodes[2]);

                multiMesh->insertNextTriangle({nodes[0], nodes[1], nodes[2]});
                auto matIds = multiMesh->cellData.getScalarData("MaterialIds");
                if (matIds)
                  matIds->push_back(matId);
              }
            }
          }
        }
      }
    }
  }

public:
  ToHullMesh() = default;

  ToHullMesh(SmartPointer<Mesh<T>> passedMesh) { multiMesh = passedMesh; }

  ToHullMesh(SmartPointer<Mesh<T>> passedMesh,
             SmartPointer<Domain<T, D>> levelSet) {
    multiMesh = passedMesh;
    levelSets.push_back(levelSet);
  }

  ToHullMesh(SmartPointer<Mesh<T>> passedMesh,
             std::vector<SmartPointer<Domain<T, D>>> const &levelSetsVector) {
    multiMesh = passedMesh;
    levelSets = levelSetsVector;
  }

  void setMesh(SmartPointer<Mesh<T>> passedMesh) { multiMesh = passedMesh; }

  /// Level sets wrapping other level sets have to be inserted last.
  void insertNextLevelSet(SmartPointer<Domain<T, D>> levelSet) {
    levelSets.push_back(levelSet);
  }

  void clearLevelSets() { levelSets.clear(); }

  void setSharpCorners(bool passedSharpCorners) {
    generateSharpCorners = passedSharpCorners;
  }

  void setMaterialMap(SmartPointer<MaterialMap> passedMaterialMap) {
    materialMap = passedMaterialMap;
  }

  void apply() { generateHull(); }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(ToHullMesh)

} // namespace viennals
