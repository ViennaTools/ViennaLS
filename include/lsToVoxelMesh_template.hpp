#ifndef LS_TO_VOXEL_MESH_TEMPLATE_HPP
#define LS_TO_VOXEL_MESH_TEMPLATE_HPP

#include <hrleDenseCellIterator.hpp>
#include <lsDomain_template.hpp>
// #include <hrleSparseCellIterator.hpp>
#include <lsMesh.hpp>

/// Creates a mesh, which consists only of quads/hexas for completely
/// filled grid cells in the level set
template <class T, int D> class lsToVoxelMesh {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;

  const lsDomain<T, D> &levelSet;
  lsMesh &mesh;

  lsToVoxelMesh();

public:
  lsToVoxelMesh(const lsDomain<T, D> &passedLevelSet, lsMesh &passedMesh)
      : levelSet(passedLevelSet), mesh(passedMesh) {}

  void apply() {
    mesh.clear();
    auto &grid = levelSet.getGrid();
    auto &domain = levelSet.getDomain();

    hrleVectorType<hrleIndexType, D> minIndex, maxIndex;
    for (unsigned i = 0; i < D; ++i) {
      minIndex[i] = (grid.isNegBoundaryInfinite(i)) ? domain.getMinRunBreak(i)
                                                    : grid.getMinBounds(i);
      maxIndex[i] = (grid.isPosBoundaryInfinite(i)) ? domain.getMaxRunBreak(i)
                                                    : grid.getMaxBounds(i);
    }

    std::vector<std::vector<hrleVectorType<double, 3>>> allNodes;
    std::vector<std::vector<hrleVectorType<unsigned, 3>>> allTriangles;
    std::vector<std::vector<hrleVectorType<unsigned, 8>>> allHexas;
    allNodes.resize(domain.getNumberOfSegments());
    allTriangles.resize(domain.getNumberOfSegments());
    allHexas.resize(domain.getNumberOfSegments());

    // TODO: parallel section does not really help here, since
    // we are memory bound --> might force single core for simplicity
    #pragma omp parallel num_threads(domain.getNumberOfSegments())
    {
      int p = 0;
    #ifdef _OPENMP
      p = omp_get_thread_num();
    #endif

      hrleVectorType<hrleIndexType, D> startVector =
          (p == 0) ? minIndex
                   : domain.getSegmentation()[p - 1];

      hrleVectorType<hrleIndexType, D> endVector =
          (p != static_cast<int>(domain.getNumberOfSegments() - 1))
              ? domain.getSegmentation()[p]
              : maxIndex;

      auto &nodes = allNodes[p];
      auto &triangles = allTriangles[p];
      auto &hexas = allHexas[p];

      // iterate over all active points
      for (hrleConstDenseCellIterator<hrleDomainType> cellIt(
               levelSet.getDomain(), startVector);
           cellIt.getIndices() < endVector; cellIt.next()) {

        unsigned signs = 0;
        for (int i = 0; i < (1 << D); i++) {
          if (cellIt.getCorner(i).getValue() >= T(0))
            signs |= (1 << i);
        }

        // all corners are negative, make a hexa
        if (signs == 0) {
          unsigned voxel[1 << D];
          for (unsigned i = 0; i < (1 << D); i++) {
            hrleVectorType<double, 3> coord(0.);
            auto index = cellIt.getIndices() + cellIt.getCorner(i).getOffset();
            for (unsigned j = 0; j < D; ++j) {
              coord[j] = levelSet.getGrid().getGridDelta() * index[j];
            }
            voxel[i] = std::distance(
                nodes.begin(),
                std::find(nodes.begin(), nodes.end(), coord));
            if (voxel[i] == nodes.size()) {
              nodes.push_back(coord);
            }
          }
          if (D == 3) {
            // reorder elements for hexas to be ordered correctly
            hrleVectorType<unsigned, 8> hexa(voxel);
            std::swap(hexa[2], hexa[3]);
            std::swap(hexa[6], hexa[7]);
            hexas.push_back(hexa);
          } else {
            hrleVectorType<unsigned, 3> triangle(voxel[0], voxel[1], voxel[2]);
            triangles.push_back(triangle);
            triangle[0] = voxel[3];
            triangles.push_back(triangle);
          }
        }
      }
    }
    // now put all nodes / triangles / hexas together and push to mesh
    unsigned offset = 0;
    unsigned numberOfTriangles = 0;
    unsigned numberOfHexas = 0;
    for(unsigned i = 0; i < domain.getNumberOfSegments(); ++i){
      numberOfTriangles += allTriangles[i].size();
      numberOfHexas += allHexas[i].size();
    }
    mesh.triangles.reserve(numberOfTriangles);
    mesh.hexas.reserve(numberOfHexas);

    for(unsigned i=0; i<domain.getNumberOfSegments(); ++i){
      mesh.nodes.insert(mesh.nodes.end(), allNodes[i].begin(), allNodes[i].end());

      hrleVectorType<unsigned, 3> tOffset(offset);
      hrleVectorType<unsigned, 8> hOffset(offset);
      for(unsigned j = 0; j < allTriangles[i].size(); ++j){
        mesh.triangles.push_back(allTriangles[i][j] + tOffset);
      }
      for(unsigned j = 0; j < allHexas[i].size(); ++j){
        mesh.hexas.push_back(allHexas[i][j] + hOffset);
      }
      // mesh.triangles.insert(mesh.triangles.end(), allTriangles[i].begin(), allTriangles[i].end());
      // mesh.hexas.insert(mesh.hexas.end(), allHexas[i].begin(), allHexas[i].end());
      offset += allNodes[i].size();
    }

    // now need to get rid of duplicate nodes
    if(domain.getNumberOfSegments()>1)
      mesh.removeDuplicateNodes();
  }
};

#endif // LS_TO_VOXEL_MESH_TEMPLATE_HPP
