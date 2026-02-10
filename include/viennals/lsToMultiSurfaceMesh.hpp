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

  SmartPointer<MaterialMap> materialMap = nullptr;

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

          if constexpr (D == 2) {
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
            } else {
              // if node does not exist yet
              nodeNumbers[n] = this->getNode(cellIt, edge, nodes, nullptr);
            }
          }

          this->insertElement(nodeNumbers);
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