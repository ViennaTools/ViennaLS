#pragma once

#include <unordered_map>

#include <lsPreCompileMacros.hpp>

#include <hrleDenseCellIterator.hpp>
#include <lsDomain.hpp>
#include <lsMaterialMap.hpp>
#include <lsMesh.hpp>

namespace viennals {

using namespace viennacore;

/// Creates a mesh, which consists only of quads/hexas for completely
/// filled grid cells in the level set. Interfaces will not be smooth
/// but stepped. (This can be used to create meshes for finite difference
/// algorithms)
template <class T, int D> class ToVoxelMesh {
  typedef typename Domain<T, D>::DomainType hrleDomainType;

  std::vector<SmartPointer<Domain<T, D>>> levelSets;
  SmartPointer<Mesh<T>> mesh = nullptr;
  SmartPointer<MaterialMap> materialMap = nullptr;
  viennahrle::Index<D> minIndex, maxIndex;

  void calculateBounds() {
    // set to zero
    for (unsigned i = 0; i < D; ++i) {
      minIndex[i] = std::numeric_limits<viennahrle::IndexType>::max();
      maxIndex[i] = std::numeric_limits<viennahrle::IndexType>::lowest();
    }
    for (unsigned l = 0; l < levelSets.size(); ++l) {
      auto &grid = levelSets[l]->getGrid();
      auto &domain = levelSets[l]->getDomain();
      for (unsigned i = 0; i < D; ++i) {
        minIndex[i] = std::min(minIndex[i], (grid.isNegBoundaryInfinite(i))
                                                ? domain.getMinRunBreak(i)
                                                : grid.getMinBounds(i));

        maxIndex[i] = std::max(maxIndex[i], (grid.isPosBoundaryInfinite(i))
                                                ? domain.getMaxRunBreak(i)
                                                : grid.getMaxBounds(i));
      }
    }
  }

public:
  ToVoxelMesh() = default;

  ToVoxelMesh(SmartPointer<Mesh<T>> passedMesh) : mesh(passedMesh) {}

  ToVoxelMesh(SmartPointer<Domain<T, D>> passedLevelSet,
              SmartPointer<Mesh<T>> passedMesh)
      : mesh(passedMesh) {
    levelSets.push_back(passedLevelSet);
  }

  ToVoxelMesh(const std::vector<SmartPointer<Domain<T, D>>> passedLevelSets,
              SmartPointer<Mesh<T>> passedMesh)
      : mesh(passedMesh) {
    levelSets = passedLevelSets;
  }

  /// Push level set to the list of level sets used for output.
  /// If more than one are specified, the voxels will be marked
  /// using a material number for each level set and output into
  /// a single mesh.
  void insertNextLevelSet(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSets.push_back(passedLevelSet);
  }

  void clearLevelSets() { levelSets.clear(); }

  void setMesh(SmartPointer<Mesh<T>> passedMesh) { mesh = passedMesh; }

  void setMaterialMap(SmartPointer<MaterialMap> passedMaterialMap) {
    materialMap = passedMaterialMap;
  }

  void apply() {
    if (levelSets.size() < 1) {
      Logger::getInstance()
          .addError("No Level Sets supplied to ToVoxelMesh! Not converting.")
          .print();
    }
    if (mesh == nullptr) {
      Logger::getInstance()
          .addError("No mesh was passed to ToVoxelMesh.")
          .print();
      return;
    }

    mesh->clear();
    auto &levelSet = *(levelSets.back());
    auto &grid = levelSet.getGrid();

    calculateBounds();

    // save the extent of the resulting mesh
    for (unsigned i = 0; i < D; ++i) {
      mesh->minimumExtent[i] = std::numeric_limits<T>::max();
      mesh->maximumExtent[i] = std::numeric_limits<T>::lowest();
    }

    std::unordered_map<viennahrle::Index<D>, size_t,
                       typename viennahrle::Index<D>::hash>
        pointIdMapping;
    size_t currentPointId = 0;

    // prepare mesh for material ids
    mesh->cellData.insertNextScalarData(typename PointData<T>::ScalarDataType(),
                                        "Material");
    auto &materialIds = *(mesh->cellData.getScalarData(0));
    const bool useMaterialMap = materialMap != nullptr;

    // set up iterators for all materials
    std::vector<
        viennahrle::ConstDenseCellIterator<typename Domain<T, D>::DomainType>>
        iterators;
    for (auto it = levelSets.begin(); it != levelSets.end(); ++it) {
      iterators.emplace_back((*it)->getDomain(), minIndex);
    }

    // move iterator for lowest material id and then adjust others if they are
    // needed
    for (; iterators.front().getIndices() < maxIndex;
         iterators.front().next()) {
      // go over all materials
      for (unsigned materialId = 0; materialId < levelSets.size();
           ++materialId) {

        auto &cellIt = iterators[materialId];

        cellIt.goToIndicesSequential(iterators.front().getIndices());

        // find out whether the centre of the box is inside
        T centerValue = 0.;
        for (int i = 0; i < (1 << D); ++i) {
          centerValue += cellIt.getCorner(i).getValue();
        }

        if (centerValue <= 0.) {
          std::array<unsigned, 1 << D> voxel;
          bool addVoxel = false;
          // now insert all points of voxel into pointList
          for (unsigned i = 0; i < (1 << D); ++i) {
            viennahrle::Index<D> index;
            addVoxel = true;
            for (unsigned j = 0; j < D; ++j) {
              index[j] =
                  cellIt.getIndices(j) + cellIt.getCorner(i).getOffset()[j];
              if (index[j] > maxIndex[j]) {
                addVoxel = false;
                break;
              }
            }
            if (addVoxel) {
              auto pointIdValue = std::make_pair(index, currentPointId);
              auto pointIdPair = pointIdMapping.insert(pointIdValue);
              voxel[i] = pointIdPair.first->second;
              if (pointIdPair.second) {
                ++currentPointId;
              }
            } else {
              break;
            }
          }

          // create element if inside domain bounds
          if (addVoxel) {
            int material = materialId;
            if (useMaterialMap)
              material = materialMap->getMaterialId(materialId);

            if constexpr (D == 3) {
              // reorder elements for hexas to be ordered correctly
              std::array<unsigned, 8> hexa{voxel[0], voxel[1], voxel[3],
                                           voxel[2], voxel[4], voxel[5],
                                           voxel[7], voxel[6]};
              mesh->hexas.push_back(hexa);
              materialIds.push_back(material);
            } else {
              std::array<unsigned, 4> tetra{voxel[0], voxel[2], voxel[3],
                                            voxel[1]};
              mesh->tetras.push_back(tetra);
              materialIds.push_back(material);
            }
          }
          // jump out of material for loop
          break;
        }
      }
    }

    // now insert points
    double gridDelta = grid.getGridDelta();
    mesh->nodes.resize(pointIdMapping.size());
    for (auto it = pointIdMapping.begin(); it != pointIdMapping.end(); ++it) {
      Vec3D<T> coords;
      for (unsigned i = 0; i < D; ++i) {
        coords[i] = gridDelta * it->first[i];

        // save extent
        if (coords[i] < mesh->minimumExtent[i]) {
          mesh->minimumExtent[i] = coords[i];
        } else if (coords[i] > mesh->maximumExtent[i]) {
          mesh->maximumExtent[i] = coords[i];
        }
      }
      mesh->nodes[it->second] = coords;
    }
  }
};

// add all template specializations for this class
PRECOMPILE_PRECISION_DIMENSION(ToVoxelMesh)

} // namespace viennals
