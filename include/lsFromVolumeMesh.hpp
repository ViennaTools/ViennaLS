#ifndef LS_FROM_VOLUME_MESH_HPP
#define LS_FROM_VOLUME_MESH_HPP

#include <lsPreCompileMacros.hpp>

#include <map>

#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMesh.hpp>
#include <lsMessage.hpp>

/// This class creates a level set from a tetrahedral mesh.
/// If the mesh contains a scalar data array called "Material",
/// one level set for each material will be created and stored
/// in the supplied std::vector<lsDomain<T,D>> object.
template <class T, int D> class lsFromVolumeMesh {
  std::vector<lsSmartPointer<lsDomain<T, D>>> levelSets;
  lsSmartPointer<lsMesh> mesh = nullptr;
  bool removeBoundaryTriangles = true;

public:
  lsFromVolumeMesh() {}

  lsFromVolumeMesh(std::vector<lsSmartPointer<lsDomain<T, D>>> passedLevelSets,
                   lsSmartPointer<lsMesh> passedMesh,
                   bool passedRemoveBoundaryTriangles = true)
      : levelSets(passedLevelSets), mesh(passedMesh),
        removeBoundaryTriangles(passedRemoveBoundaryTriangles) {}

  void
  setLevelSets(std::vector<lsSmartPointer<lsDomain<T, D>>> passedLevelSets) {
    levelSets = passedLevelSets;
  }

  void setMesh(lsSmartPointer<lsMesh> passedMesh) { mesh = passedMesh; }

  void setRemoveBoundaryTriangles(bool passedRemoveBoundaryTriangles) {
    removeBoundaryTriangles = passedRemoveBoundaryTriangles;
  }

  void apply() {
    if (levelSets.empty()) {
      lsMessage::getInstance()
          .addWarning("No level set vector was passed to lsFromVolumeMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsFromVolumeMesh.")
          .print();
      return;
    }

    // get the unique material numbers for explicit booling
    std::vector<int> materialInts;
    typename lsPointData::ScalarDataType *materialData =
        mesh->getScalarData("Material");
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
    typedef std::map<hrleVectorType<unsigned int, D>, std::pair<int, int>>
        triangleMapType;
    triangleMapType surfaceElements;

    unsigned numberOfElements =
        (D == 3) ? mesh->tetras.size() : mesh->triangles.size();
    for (unsigned int i = 0; i < numberOfElements; ++i) {
      for (int j = 0; j < D + 1; j++) {
        hrleVectorType<unsigned int, D> currentSurfaceElement;
        for (int k = 0; k < D; k++) {
          currentSurfaceElement[k] =
              mesh->getElements<D + 1>()[i][(j + k) % (D + 1)];
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

        currentSurfaceElement.sort();

        hrleVectorType<double, D> currentElementPoints[D + 1];
        for (int k = 0; k < D; k++) {
          currentElementPoints[k] = mesh->nodes[currentSurfaceElement[k]];
        }

        // get the other point of the element as well
        currentElementPoints[D] =
            mesh->nodes[mesh->getElements<D + 1>()[i][(j + D) % (D + 1)]];

        typename triangleMapType::iterator it =
            surfaceElements.lower_bound(currentSurfaceElement);
        if ((it != surfaceElements.end()) &&
            (it->first == currentSurfaceElement)) {
          if (Orientation(currentElementPoints)) {
            if (it->second.second != materialInts.back() + 1) {
              lsMessage::getInstance()
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
              lsMessage::getInstance()
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
    levelSets.resize(materialInts.size());
    auto levelSetIterator = levelSets.begin();
    for (auto matIt = materialInts.begin(); matIt != materialInts.end();
         ++matIt) {
      auto currentSurface = lsSmartPointer<lsMesh>::New();
      auto &meshElements = currentSurface->getElements<D>();
      for (auto it = surfaceElements.begin(); it != surfaceElements.end();
           ++it) {
        if (((*matIt) >= it->second.first) && ((*matIt) < it->second.second)) {
          std::array<unsigned, D> element{it->first[0], it->first[1]};
          if (D == 3)
            element[2] = it->first[2];
          meshElements.push_back(element);
        } else if (((*matIt) >= it->second.second) &&
                   ((*matIt) < it->second.first)) {
          // swap first two elements since triangle has different orientation
          std::array<unsigned, D> element{it->first[1], it->first[0]};
          if (D == 3)
            element[2] = it->first[2];
          meshElements.push_back(element);
        }
      }

      // replace Nodes of Geometry by Nodes of individual surface
      const unsigned int undefined_node =
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
      lsFromSurfaceMesh<T, D>(*levelSetIterator, currentSurface,
                              removeBoundaryTriangles)
          .apply();

      ++levelSetIterator;
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsFromVolumeMesh)

#endif // LS_FROM_VOLUME_MESH_HPP
