#ifndef LS_TO_VOXEL_MESH_HPP
#define LS_TO_VOXEL_MESH_HPP

#include <unordered_map>

#include <lsPreCompileMacros.hpp>

#include <hrleDenseCellIterator.hpp>
#include <lsDomain.hpp>
#include <lsMesh.hpp>

// std::hash specialisation in order to use unordered_map
namespace std {
template <int D> struct hash<hrleVectorType<hrleIndexType, D>> {
  std::size_t operator()(const hrleVectorType<hrleIndexType, D> &v) const {
    using std::hash;
    using std::size_t;
    using std::string;

    /* Better hash combine functions might be:
      hash(a)<<1 + hash(a) + hash(b)
      or from boost:
      size_t hash_combine( size_t lhs, size_t rhs ) {
        lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
        return lhs;
      }
      https://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes
    */

    // Compute individual hash values for first,
    // second and third and combine them using XOR
    // and bit shifting:

    size_t result = hash<hrleIndexType>()(v[0]);
    result ^= hash<hrleIndexType>()(v[1]) << 1;
    if (D == 3) {
      result = (result >> 1) ^ (hash<hrleIndexType>()(v[3]) << 1);
    }
    return result;
  }
};

} // namespace std

/// Creates a mesh, which consists only of quads/hexas for completely
/// filled grid cells in the level set. Interfaces will not be smooth
/// but stepped. (This can be used to create meshes for finite difference
/// algorithms)
template <class T, int D> class lsToVoxelMesh {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;

  std::vector<const lsDomain<T, D> *> levelSets;
  lsMesh *mesh = nullptr;
  hrleVectorType<hrleIndexType, D> minIndex, maxIndex;

  lsToVoxelMesh();

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
  lsToVoxelMesh(lsMesh &passedMesh) : mesh(&passedMesh) {}

  lsToVoxelMesh(const lsDomain<T, D> &passedLevelSet, lsMesh &passedMesh)
      : mesh(&passedMesh) {
    levelSets.push_back(&passedLevelSet);
  }

  lsToVoxelMesh(const std::vector<const lsDomain<T, D> *> &passedLevelSets,
                lsMesh &passedMesh)
      : mesh(&passedMesh) {
    levelSets = passedLevelSets;
  }

  void insertNextLevelSet(const lsDomain<T, D> &passedLevelSet) {
    levelSets.push_back(&passedLevelSet);
  }

  void setMesh(lsMesh &passedMesh) { mesh = &passedMesh; }

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

    std::unordered_map<hrleVectorType<hrleIndexType, D>, size_t> pointIdMapping;
    size_t currentPointId = 0;

    // prepare mesh for material ids
    mesh->scalarDataLabels.push_back("Material ID");
    mesh->scalarData.push_back(std::vector<double>());

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
          centerValue = cellIt.getCorner(i).getValue();
        }

        if (centerValue <= 0.) {
          unsigned voxel[1 << D];
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
            hrleVectorType<unsigned, 8> hexa(voxel);
            std::swap(hexa[2], hexa[3]);
            std::swap(hexa[6], hexa[7]);
            mesh->hexas.push_back(hexa);
            mesh->scalarData[0].push_back(materialId);
          } else {
            hrleVectorType<unsigned, 3> triangle(voxel[0], voxel[1], voxel[2]);
            mesh->triangles.push_back(triangle);
            mesh->scalarData[0].push_back(materialId);
            triangle[0] = voxel[3];
            mesh->triangles.push_back(triangle);
            mesh->scalarData[0].push_back(materialId);
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
      hrleVectorType<double, D> coords;
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
