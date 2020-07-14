#ifndef LS_TO_VOXEL_MESH_HPP
#define LS_TO_VOXEL_MESH_HPP

#include <unordered_map>

#include <lsPreCompileMacros.hpp>

#include <hrleDenseCellIterator.hpp>
#include <lsDomain.hpp>
#include <lsMesh.hpp>

/// Creates a mesh, which consists only of quads/hexas for completely
/// filled grid cells in the level set. Interfaces will not be smooth
/// but stepped. (This can be used to create meshes for finite difference
/// algorithms)
template <class T, int D> class lsToVoxelMesh {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;

  std::vector<lsSmartPointer<lsDomain<T, D>>> levelSets;
  lsSmartPointer<lsMesh> mesh = nullptr;
  hrleVectorType<hrleIndexType, D> minIndex, maxIndex;

  void calculateBounds() {
    // set to zero
    for (unsigned i = 0; i < D; ++i) {
      minIndex[i] = std::numeric_limits<hrleIndexType>::max();
      maxIndex[i] = std::numeric_limits<hrleIndexType>::lowest();
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
  lsToVoxelMesh() {}

  lsToVoxelMesh(lsSmartPointer<lsMesh> passedMesh) : mesh(passedMesh) {}

  lsToVoxelMesh(lsSmartPointer<lsDomain<T, D>> passedLevelSet,
                lsSmartPointer<lsMesh> passedMesh)
      : mesh(passedMesh) {
    levelSets.push_back(passedLevelSet);
  }

  lsToVoxelMesh(
      const std::vector<lsSmartPointer<lsDomain<T, D>>> passedLevelSets,
      lsSmartPointer<lsMesh> passedMesh)
      : mesh(passedMesh) {
    levelSets = passedLevelSets;
  }

  /// Push level set to the list of level sets used for output.
  /// If more than one are specified, the voxels will be marked
  /// using a material number for each level set and output into
  /// a single mesh.
  void insertNextLevelSet(lsSmartPointer<lsDomain<T, D>> passedLevelSet) {
    levelSets.push_back(passedLevelSet);
  }

  void setMesh(lsSmartPointer<lsMesh> passedMesh) { mesh = passedMesh; }

  void apply() {
    if (levelSets.size() < 1) {
      lsMessage::getInstance()
          .addWarning(
              "No Level Sets supplied to lsToVoxelMesh! Not converting.")
          .print();
    }
    if (mesh == nullptr) {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsToVoxelMesh.")
          .print();
      return;
    }

    mesh->clear();
    auto &levelSet = *(levelSets.back());
    auto &grid = levelSet.getGrid();

    calculateBounds();

    std::unordered_map<hrleVectorType<hrleIndexType, D>, size_t,
                       typename hrleVectorType<hrleIndexType, D>::hash>
        pointIdMapping;
    size_t currentPointId = 0;

    // prepare mesh for material ids
    mesh->insertNextScalarData(lsPointData::ScalarDataType(), "Material");
    auto &materialIds = *(mesh->getScalarData(0));

    // set up iterators for all materials
    std::vector<hrleConstDenseCellIterator<typename lsDomain<T, D>::DomainType>>
        iterators;
    for (auto it = levelSets.begin(); it != levelSets.end(); ++it) {
      iterators.push_back(
          hrleConstDenseCellIterator<typename lsDomain<T, D>::DomainType>(
              (*it)->getDomain(), minIndex));
    }

    // move iterator for lowest material id and then adjust others if they are
    // needed
    for (; iterators.back().getIndices() < maxIndex; iterators.back().next()) {
      // go over all materials
      for (unsigned materialId = 0; materialId < levelSets.size();
           ++materialId) {

        auto &cellIt = iterators[materialId];

        cellIt.goToIndicesSequential(iterators.back().getIndices());

        // find out whether the centre of the box is inside
        T centerValue = 0.;
        for (int i = 0; i < (1 << D); ++i) {
          centerValue += cellIt.getCorner(i).getValue();
        }

        if (centerValue <= 0.) {
          std::array<unsigned, 1 << D> voxel;
          // now insert all points of voxel into pointList
          for (unsigned i = 0; i < (1 << D); ++i) {
            hrleVectorType<hrleIndexType, D> index;
            for (unsigned j = 0; j < D; ++j) {
              index[j] =
                  cellIt.getIndices(j) + cellIt.getCorner(i).getOffset()[j];
            }
            auto pointIdValue = std::make_pair(index, currentPointId);

            auto pointIdPair = pointIdMapping.insert(pointIdValue);
            voxel[i] = pointIdPair.first->second;
            if (pointIdPair.second) {
              ++currentPointId;
            }
          }

          // create element
          if (D == 3) {
            // reorder elements for hexas to be ordered correctly
            std::array<unsigned, 8> hexa{voxel[0], voxel[1], voxel[3],
                                         voxel[2], voxel[4], voxel[5],
                                         voxel[7], voxel[6]};
            mesh->hexas.push_back(hexa);
            materialIds.push_back(materialId);
          } else {
            std::array<unsigned, 3> triangle{voxel[0], voxel[1], voxel[2]};
            mesh->triangles.push_back(triangle);
            materialIds.push_back(materialId);
            triangle[0] = voxel[3];
            mesh->triangles.push_back(triangle);
            materialIds.push_back(materialId);
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
      std::array<double, 3> coords{};
      for (unsigned i = 0; i < D; ++i) {
        coords[i] = gridDelta * it->first[i];
      }
      mesh->nodes[it->second] = coords;
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsToVoxelMesh)

#endif // LS_TO_VOXEL_MESH_HPP
