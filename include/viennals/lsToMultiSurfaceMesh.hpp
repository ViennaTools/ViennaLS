#pragma once

#include <lsPreCompileMacros.hpp>

#include <iostream>
#include <map>
#include <unordered_map>
#include <unordered_set>

#include <hrleSparseCellIterator.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMarchingCubes.hpp>
#include <lsMaterialMap.hpp>
#include <lsMesh.hpp>

namespace viennals {

using namespace viennacore;

template <class NumericType, int D, class T = NumericType>
class ToMultiSurfaceMesh {

  typedef Domain<NumericType, D> lsDomainType;
  typedef typename Domain<NumericType, D>::DomainType hrleDomainType;

  std::vector<SmartPointer<lsDomainType>> levelSets;
  SmartPointer<Mesh<T>> mesh = nullptr;
  SmartPointer<MaterialMap> materialMap = nullptr;

  const double epsilon;
  const double minNodeDistanceFactor;

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

public:
  ToMultiSurfaceMesh(double eps = 1e-12, double minNodeDistFactor = 0.05)
      : epsilon(eps), minNodeDistanceFactor(minNodeDistFactor) {}

  ToMultiSurfaceMesh(SmartPointer<lsDomainType> passedLevelSet,
                     SmartPointer<Mesh<T>> passedMesh, double eps = 1e-12,
                     double minNodeDistFactor = 0.05)
      : mesh(passedMesh), epsilon(eps),
        minNodeDistanceFactor(minNodeDistFactor) {
    levelSets.push_back(passedLevelSet);
  }

  ToMultiSurfaceMesh(
      std::vector<SmartPointer<lsDomainType>> const &passedLevelSets,
      SmartPointer<Mesh<T>> passedMesh, double eps = 1e-12,
      double minNodeDistFactor = 0.05)
      : levelSets(passedLevelSets), mesh(passedMesh), epsilon(eps),
        minNodeDistanceFactor(minNodeDistFactor) {}

  ToMultiSurfaceMesh(SmartPointer<Mesh<T>> passedMesh, double eps = 1e-12,
                     double minNodeDistFactor = 0.05)
      : mesh(passedMesh), epsilon(eps),
        minNodeDistanceFactor(minNodeDistFactor) {}

  void setMesh(SmartPointer<Mesh<T>> passedMesh) { mesh = passedMesh; }

  void insertNextLevelSet(SmartPointer<lsDomainType> passedLevelSet) {
    levelSets.push_back(passedLevelSet);
  }

  void clearLevelSets() { levelSets.clear(); }

  void setMaterialMap(SmartPointer<MaterialMap> passedMaterialMap) {
    materialMap = passedMaterialMap;
  }

  void apply() {
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
    const auto gridDelta = levelSets.front()->getGrid().getGridDelta();
    for (unsigned i = 0; i < D; ++i) {
      mesh->minimumExtent[i] = std::numeric_limits<NumericType>::max();
      mesh->maximumExtent[i] = std::numeric_limits<NumericType>::lowest();
    }

    constexpr unsigned int corner0[12] = {0, 1, 2, 0, 4, 5, 6, 4, 0, 1, 3, 2};
    constexpr unsigned int corner1[12] = {1, 3, 3, 2, 5, 7, 7, 6, 4, 5, 7, 6};
    constexpr unsigned int direction[12] = {0, 1, 0, 1, 0, 1, 0, 1, 2, 2, 2, 2};

    typedef std::map<viennahrle::Index<D>, unsigned> nodeContainerType;

    nodeContainerType nodes[D];
    const double minNodeDistance = gridDelta * minNodeDistanceFactor;
    std::unordered_map<I3, unsigned, I3Hash> nodeIdByBin;
    std::unordered_set<I3, I3Hash> uniqueElements;

    typename nodeContainerType::iterator nodeIt;

    std::vector<Vec3D<T>> normals;
    std::vector<T> materials;

    const bool checkNodeFlag = minNodeDistanceFactor > 0;
    const bool useMaterialMap = materialMap != nullptr;

    auto quantize = [&minNodeDistance](const Vec3D<T> &p) -> I3 {
      const T inv = T(1) / minNodeDistance;
      return {(int)std::llround(p[0] * inv), (int)std::llround(p[1] * inv),
              (int)std::llround(p[2] * inv)};
    };

    auto elementI3 = [&](const std::array<unsigned, D> &element) -> I3 {
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
      // iterate over all active surface points
      for (auto cellIt = cellIts[l]; !cellIt.isFinished(); cellIt.next()) {
        for (int u = 0; u < D; u++) {
          while (!nodes[u].empty() &&
                 nodes[u].begin()->first <
                     viennahrle::Index<D>(cellIt.getIndices()))
            nodes[u].erase(nodes[u].begin());
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

        // for each element
        const int *Triangles;
        if constexpr (D == 2) {
          Triangles = lsInternal::MarchingCubes::polygonize2d(signs);
        } else {
          Triangles = lsInternal::MarchingCubes::polygonize3d(signs);
        }

        for (; Triangles[0] != -1; Triangles += D) {
          std::array<unsigned, D> nod_numbers;

          // for each node
          for (int n = 0; n < D; n++) {
            const int edge = Triangles[n];

            unsigned p0 = corner0[edge];
            unsigned p1 = corner1[edge];

            // determine direction of edge
            unsigned dir = direction[edge];

            // look for existing surface node
            viennahrle::Index<D> d(cellIt.getIndices());
            auto p0B = viennahrle::BitMaskToIndex<D>(p0);
            d += p0B;

            nodeIt = nodes[dir].find(d);
            if (nodeIt != nodes[dir].end()) {
              nod_numbers[n] = nodeIt->second;
            } else {
              // if node does not exist yet
              // calculate coordinate of new node
              Vec3D<T> cc{}; // initialise with zeros
              for (int z = 0; z < D; z++) {
                if (z != dir) {
                  // TODO might not need BitMaskToVector here, just check if z
                  // bit is set
                  cc[z] = static_cast<T>(cellIt.getIndices(z) + p0B[z]);
                } else {
                  auto d0 = static_cast<T>(cellIt.getCorner(p0).getValue());
                  auto d1 = static_cast<T>(cellIt.getCorner(p1).getValue());

                  // calculate the surface-grid intersection point
                  if (d0 == -d1) { // includes case where d0=d1=0
                    cc[z] = static_cast<T>(cellIt.getIndices(z)) + 0.5;
                  } else {
                    if (std::abs(d0) <= std::abs(d1)) {
                      cc[z] = static_cast<T>(cellIt.getIndices(z)) +
                              (d0 / (d0 - d1));
                    } else {
                      cc[z] = static_cast<T>(cellIt.getIndices(z) + 1) -
                              (d1 / (d1 - d0));
                    }
                  }
                  cc[z] = std::max(
                      cc[z], static_cast<T>(cellIt.getIndices(z) + epsilon));
                  cc[z] = std::min(
                      cc[z],
                      static_cast<T>((cellIt.getIndices(z) + 1) - epsilon));
                }
                cc[z] *= gridDelta;
              }

              int nodeIdx = -1;
              if (checkNodeFlag) {
                auto q = quantize(cc);
                auto it = nodeIdByBin.find(q);
                if (it != nodeIdByBin.end())
                  nodeIdx = it->second;
              }
              if (nodeIdx >= 0) {
                nod_numbers[n] = nodeIdx;
              } else {
                // insert new node
                nod_numbers[n] =
                    mesh->insertNextNode(cc); // insert new surface node
                nodes[dir][d] = nod_numbers[n];
                if (checkNodeFlag)
                  nodeIdByBin.emplace(quantize(cc), nod_numbers[n]);

                for (int a = 0; a < D; a++) {
                  if (cc[a] < mesh->minimumExtent[a])
                    mesh->minimumExtent[a] = cc[a];
                  if (cc[a] > mesh->maximumExtent[a])
                    mesh->maximumExtent[a] = cc[a];
                }
              }
            }
          }

          // check for degenerate triangles
          if (!triangleMisformed(nod_numbers)) {

            // only insert if not already present
            if (uniqueElements.insert(elementI3(nod_numbers)).second) {
              Vec3D<T> normal;
              if constexpr (D == 2) {
                normal = Vec3D<T>{-(mesh->nodes[nod_numbers[1]][1] -
                                    mesh->nodes[nod_numbers[0]][1]),
                                  mesh->nodes[nod_numbers[1]][0] -
                                      mesh->nodes[nod_numbers[0]][0],
                                  T(0)};
              } else {
                normal = calculateNormal(mesh->nodes[nod_numbers[0]],
                                         mesh->nodes[nod_numbers[1]],
                                         mesh->nodes[nod_numbers[2]]);
              }

              double n2 = normal[0] * normal[0] + normal[1] * normal[1] +
                          normal[2] * normal[2];
              if (n2 > epsilon) {
                // insert new element
                mesh->insertNextElement(nod_numbers);

                // insert normal
                T invn = static_cast<T>(1.) / std::sqrt(static_cast<T>(n2));
                for (int d = 0; d < D; d++) {
                  normal[d] *= invn;
                }
                normals.push_back(normal);

                // insert material
                if (useMaterialMap) {
                  materials.push_back(materialMap->getMaterialId(l));
                } else {
                  materials.push_back(static_cast<T>(l));
                }
              }
            }
          }
        }
      }
    }

    mesh->cellData.insertNextScalarData(materials, "MaterialIds");
    mesh->cellData.insertNextVectorData(normals, "Normals");
    mesh->triangles.shrink_to_fit();
    mesh->nodes.shrink_to_fit();
  }

private:
  static inline bool
  triangleMisformed(const std::array<unsigned, D> &nod_numbers) noexcept {
    if constexpr (D == 3) {
      return nod_numbers[0] == nod_numbers[1] ||
             nod_numbers[0] == nod_numbers[2] ||
             nod_numbers[1] == nod_numbers[2];
    } else {
      return nod_numbers[0] == nod_numbers[1];
    }
  }

  static inline Vec3D<NumericType>
  calculateNormal(const Vec3D<NumericType> &nodeA,
                  const Vec3D<NumericType> &nodeB,
                  const Vec3D<NumericType> &nodeC) noexcept {
    return CrossProduct(nodeB - nodeA, nodeC - nodeA);
  }
};

PRECOMPILE_PRECISION_DIMENSION(ToMultiSurfaceMesh)

} // namespace viennals