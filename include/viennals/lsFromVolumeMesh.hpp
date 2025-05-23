#pragma once

#include <lsPreCompileMacros.hpp>

#include <map>

#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMesh.hpp>

#include <vcLogger.hpp>
#include <vcSmartPointer.hpp>

namespace viennals {

using namespace viennacore;

/// This class creates a level set from a tetrahedral mesh.
/// If the mesh contains a scalar data array called "Material",
/// one level set for each material will be created and stored
/// in the supplied std::vector<Domain<T,D>> object.
template <class T, int D> class FromVolumeMesh {
public:
  using LevelSetType = SmartPointer<Domain<T, D>>;
  using LevelSetsType = std::vector<LevelSetType>;
  using GridType = typename Domain<T, D>::GridType;

private:
  LevelSetsType levelSets;
  SmartPointer<Mesh<T>> mesh = nullptr;
  GridType grid;
  bool removeBoundaryTriangles = true;
  bool gridSet = false;

public:
  FromVolumeMesh() = default;

  FromVolumeMesh(const GridType &passedGrid, SmartPointer<Mesh<T>> passedMesh,
                 bool passedRemoveBoundaryTriangles = true)
      : mesh(passedMesh), grid(passedGrid),
        removeBoundaryTriangles(passedRemoveBoundaryTriangles), gridSet(true) {}

  void setGrid(const GridType &passedGrid) {
    grid = passedGrid;
    gridSet = true;
  }

  void setMesh(SmartPointer<Mesh<T>> passedMesh) { mesh = passedMesh; }

  void setRemoveBoundaryTriangles(bool passedRemoveBoundaryTriangles) {
    removeBoundaryTriangles = passedRemoveBoundaryTriangles;
  }

  std::vector<LevelSetType> getLevelSets() const { return levelSets; }

  void apply() {
    if (!gridSet) {
      Logger::getInstance()
          .addWarning("No grid has been set in FromVolumeMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      Logger::getInstance()
          .addWarning("No mesh was passed to FromVolumeMesh.")
          .print();
      return;
    }

    // get the unique material numbers for explicit booling
    std::vector<int> materialInts;
    typename PointData<T>::ScalarDataType *materialData =
        mesh->cellData.getScalarData("Material", true);
    if (materialData != nullptr) {
      // make unique list of materialIds
      materialInts =
          std::vector<int>(materialData->begin(), materialData->end());
      std::sort(materialInts.begin(), materialInts.end());
      auto it = std::unique(materialInts.begin(), materialInts.end());
      materialInts.erase(it, materialInts.end());
    } else {
      // no materials are defined
      materialInts.push_back(0);
    }

    // Map for all surfaceElements and their corresponding material
    typedef std::map<VectorType<unsigned int, D>, std::pair<int, int>>
        triangleMapType;
    triangleMapType surfaceElements;

    unsigned numberOfElements =
        (D == 3) ? mesh->tetras.size() : mesh->triangles.size();
    for (unsigned int i = 0; i < numberOfElements; ++i) {
      for (int j = 0; j < D + 1; j++) {
        VectorType<unsigned int, D> currentSurfaceElement;
        for (int k = 0; k < D; k++) {
          currentSurfaceElement[k] =
              mesh->template getElements<D + 1>()[i][(j + k) % (D + 1)];
        }

        // std::bitset<2 * D> flags;
        // flags.set();
        //
        // // if triangle at border skip
        // for (int k = 0; k < D; k++) {
        //   for (int l = 0; l < D; l++) {
        //     if (mesh->nodes[currentSurfaceElement[k]][l] < Geometry.Max[l] -
        //     eps) {
        //       flags.reset(l + D);
        //     }
        //     if (mesh->nodes[currentSurfaceElement[k]][l] > Geometry.Min[l] +
        //     eps) {
        //       flags.reset(l);
        //     }
        //   }
        // }
        //
        // flags &= remove_flags;

        // if (is_open_boundary_negative) flags.reset(open_boundary_direction);
        // else flags.reset(open_boundary_direction+D);
        //
        // if (flags.any())
        //   continue;

        // currentSurfaceElement.sort();
        std::sort(currentSurfaceElement.begin(), currentSurfaceElement.end());

        Vec3D<T> currentElementPoints[D + 1];
        for (int k = 0; k < D; k++) {
          currentElementPoints[k] = mesh->nodes[currentSurfaceElement[k]];
        }

        // get the other point of the element as well
        currentElementPoints[D] =
            mesh->nodes[mesh->template getElements<D + 1>()[i]
                                                           [(j + D) % (D + 1)]];

        typename triangleMapType::iterator it =
            surfaceElements.lower_bound(currentSurfaceElement);
        if ((it != surfaceElements.end()) &&
            (it->first == currentSurfaceElement)) {
          if (Orientation(currentElementPoints)) {
            if (it->second.second != materialInts.back() + 1) {
              Logger::getInstance()
                  .addWarning(
                      "Coinciding surface elements with same orientation in "
                      "Element: " +
                      std::to_string(i))
                  .print();
            }
            it->second.second =
                (materialData == nullptr) ? 0 : (*materialData)[i];
          } else {
            if (it->second.first != materialInts.back() + 1) {
              Logger::getInstance()
                  .addWarning(
                      "Coinciding surface elements with same orientation in "
                      "Element: " +
                      std::to_string(i))
                  .print();
            }
            it->second.first =
                (materialData == nullptr) ? 0 : (*materialData)[i];
          }

          if (it->second.first == it->second.second)
            surfaceElements.erase(it);

        } else {
          if (Orientation(currentElementPoints)) {
            surfaceElements.insert(
                it, std::make_pair(currentSurfaceElement,
                                   std::make_pair(materialInts.back() + 1,
                                                  (materialData == nullptr)
                                                      ? 0
                                                      : (*materialData)[i])));
          } else {
            surfaceElements.insert(
                it, std::make_pair(currentSurfaceElement,
                                   std::make_pair((materialData == nullptr)
                                                      ? 0
                                                      : (*materialData)[i],
                                                  materialInts.back() + 1)));
          }
        }
      }
    }

    // for all materials/for each surface
    // resize to empty levelsets, but they need grid information
    {
      levelSets.clear();
      for (unsigned i = 0; i < materialInts.size(); ++i) {
        auto ls = LevelSetType::New(grid);
        levelSets.push_back(ls);
      }
    }

    auto levelSetIterator = levelSets.begin();
    for (int &materialInt : materialInts) {
      auto currentSurface = SmartPointer<Mesh<T>>::New();
      auto &meshElements = currentSurface->template getElements<D>();
      for (auto it = surfaceElements.begin(); it != surfaceElements.end();
           ++it) {
        if ((materialInt >= it->second.first) &&
            (materialInt < it->second.second)) {
          std::array<unsigned, D> element{it->first[0], it->first[1]};
          if constexpr (D == 3)
            element[2] = it->first[2];
          meshElements.push_back(element);
        } else if ((materialInt >= it->second.second) &&
                   (materialInt < it->second.first)) {
          // swap first two elements since triangle has different orientation
          std::array<unsigned, D> element{it->first[1], it->first[0]};
          if constexpr (D == 3)
            element[2] = it->first[2];
          meshElements.push_back(element);
        }
      }

      // replace Nodes of Geometry by Nodes of individual surface
      constexpr unsigned int undefined_node =
          std::numeric_limits<unsigned int>::max();
      std::vector<unsigned int> nodeReplacements(mesh->nodes.size(),
                                                 undefined_node);
      unsigned int NodeCounter = 0;

      for (unsigned int k = 0; k < meshElements.size(); ++k) {
        for (int h = 0; h < D; h++) {
          unsigned int origin_node = meshElements[k][h];
          if (nodeReplacements[origin_node] == undefined_node) {
            nodeReplacements[origin_node] = NodeCounter++;
            currentSurface->nodes.push_back(mesh->nodes[origin_node]);
          }
          meshElements[k][h] = nodeReplacements[origin_node];
        }
      }

      // create level set from surface
      FromSurfaceMesh<T, D>(*levelSetIterator, currentSurface,
                            removeBoundaryTriangles)
          .apply();

      ++levelSetIterator;
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(FromVolumeMesh)

} // namespace viennals
