#ifndef LS_FROM_VOLUME_MESH_HPP
#define LS_FROM_VOLUME_MESH_HPP

#include <lsPreCompileMacros.hpp>

#include <map>

#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMesh.hpp>
#include <lsMessage.hpp>
#include <lsVTKWriter.hpp> // TODO remove

/**
  This class creates a level set from a tetrahedral mesh.
  If the mesh contains a scalar data array called "Material"
  One level set for each material will be created and stored
  in the supplied std::vector<lsDomain<T,D>> object.
*/
template <class T, int D> class lsFromVolumeMesh {
  std::vector<lsDomain<T, D>> *levelSets = nullptr;
  lsMesh *mesh = nullptr;
  bool removeBoundaryTriangles = true;

public:
  lsFromVolumeMesh(std::vector<lsDomain<T, D>> &passedLevelSets,
                   lsMesh &passedMesh,
                   bool passedRemoveBoundaryTriangles = true)
      : levelSets(&passedLevelSets), mesh(&passedMesh),
        removeBoundaryTriangles(passedRemoveBoundaryTriangles) {}

  void setLevelSets(std::vector<lsDomain<T, D>> &passedLevelSets) {
    levelSets = &passedLevelSets;
  }

  void setMesh(lsMesh &passedMesh) { mesh = &passedMesh; }

  void setRemoveBoundaryTriangles(bool passedRemoveBoundaryTriangles) {
    removeBoundaryTriangles = passedRemoveBoundaryTriangles;
  }

  void apply() {
    if (levelSets == nullptr) {
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
    int scalarDataIndex = -1;
    // see if there is a scalar data array specifying "Material"
    {
      auto it = std::find(mesh->scalarDataLabels.begin(),
                          mesh->scalarDataLabels.end(), "Material");
      if (it != mesh->scalarDataLabels.end()) {
        scalarDataIndex = std::distance(mesh->scalarDataLabels.begin(), it);

        for (auto materialIt = mesh->scalarData[scalarDataIndex].begin();
             materialIt != mesh->scalarData[scalarDataIndex].end();
             ++materialIt) {
          if (std::find(materialInts.begin(), materialInts.end(),
                        *materialIt) == materialInts.end()) {
            materialInts.push_back(static_cast<int>(*materialIt));
          }
        }
        std::sort(materialInts.begin(), materialInts.end());
      } else {
        materialInts.push_back(0);
      }
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
            it->second.second = (scalarDataIndex == -1)
                                    ? 0
                                    : mesh->scalarData[scalarDataIndex][i];
          } else {
            if (it->second.first != materialInts.back() + 1) {
              lsMessage::getInstance()
                  .addWarning(
                      "Coinciding surface elements with same orientation in "
                      "Element: " +
                      std::to_string(i))
                  .print();
            }
            it->second.first = (scalarDataIndex == -1)
                                   ? 0
                                   : mesh->scalarData[scalarDataIndex][i];
          }

          if (it->second.first == it->second.second)
            surfaceElements.erase(it);

        } else {
          if (Orientation(currentElementPoints)) {
            surfaceElements.insert(
                it, std::make_pair(
                        currentSurfaceElement,
                        std::make_pair(
                            materialInts.back() + 1,
                            (scalarDataIndex == -1)
                                ? 0
                                : mesh->scalarData[scalarDataIndex][i])));
          } else {
            surfaceElements.insert(
                it,
                std::make_pair(
                    currentSurfaceElement,
                    std::make_pair((scalarDataIndex == -1)
                                       ? 0
                                       : mesh->scalarData[scalarDataIndex][i],
                                   materialInts.back() + 1)));
          }
        }
      }
    }

    // for all materials/for each surface
    levelSets->resize(materialInts.size());
    auto levelSetIterator = levelSets->begin();
    for (auto matIt = materialInts.begin(); matIt != materialInts.end();
         ++matIt) {
      lsMesh currentSurface;
      auto &meshElements = currentSurface.getElements<D>();
      for (auto it = surfaceElements.begin(); it != surfaceElements.end();
           ++it) {
        if (((*matIt) >= it->second.first) && ((*matIt) < it->second.second)) {
          meshElements.push_back(it->first);
        } else if (((*matIt) >= it->second.second) &&
                   ((*matIt) < it->second.first)) {
          meshElements.push_back(it->first);
          std::swap(meshElements.back()[0], meshElements.back()[1]);
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
            currentSurface.nodes.push_back(mesh->nodes[origin_node]);
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
