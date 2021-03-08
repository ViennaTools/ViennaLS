#ifndef LS_TO_SURFACE_MESH_HPP
#define LS_TO_SURFACE_MESH_HPP

#include <lsPreCompileMacros.hpp>

#include <iostream>
#include <map>

#include <hrleSparseCellIterator.hpp>
#include <lsDomain.hpp>
#include <lsMarchingCubes.hpp>
#include <lsMesh.hpp>
#include <lsMessage.hpp>

/// Extract an explicit lsMesh<> instance from an lsDomain.
/// The interface is then described by explciit surface elements:
/// Lines in 2D, Triangles in 3D.
template <class T, int D> class lsToSurfaceMesh {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;

  lsSmartPointer<lsDomain<T, D>> levelSet = nullptr;
  lsSmartPointer<lsMesh<T>> mesh = nullptr;
  // std::vector<hrleIndexType> meshNodeToPointIdMapping;
  const T epsilon;

public:
  lsToSurfaceMesh(double eps = 1e-12) : epsilon(eps) {}

  lsToSurfaceMesh(const lsSmartPointer<lsDomain<T, D>> passedLevelSet,
                  lsSmartPointer<lsMesh<T>> passedMesh, double eps = 1e-12)
      : levelSet(passedLevelSet), mesh(passedMesh), epsilon(eps) {}

  void setLevelSet(lsSmartPointer<lsDomain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  void setMesh(lsSmartPointer<lsMesh<T>> passedMesh) { mesh = passedMesh; }

  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsToSurfaceMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsToSurfaceMesh.")
          .print();
      return;
    }

    if (levelSet->getNumberOfPoints() == 0) {
      return;
    }

    mesh->clear();
    const unsigned int corner0[12] = {0, 1, 2, 0, 4, 5, 6, 4, 0, 1, 3, 2};
    const unsigned int corner1[12] = {1, 3, 3, 2, 5, 7, 7, 6, 4, 5, 7, 6};
    const unsigned int direction[12] = {0, 1, 0, 1, 0, 1, 0, 1, 2, 2, 2, 2};

    // test if level set function consists of at least 2 layers of
    // defined grid points
    if (levelSet->getLevelSetWidth() < 2) {
      lsMessage::getInstance()
          .addWarning("Levelset is less than 2 layers wide. Export might fail!")
          .print();
    }

    typedef typename std::map<hrleVectorType<hrleIndexType, D>, unsigned>
        nodeContainerType;

    nodeContainerType nodes[D];

    typename nodeContainerType::iterator nodeIt;

    lsInternal::lsMarchingCubes marchingCubes;

    // iterate over all active points
    for (hrleConstSparseCellIterator<hrleDomainType> cellIt(
             levelSet->getDomain());
         !cellIt.isFinished(); cellIt.next()) {

      for (int u = 0; u < D; u++) {
        while (!nodes[u].empty() &&
               nodes[u].begin()->first <
                   hrleVectorType<hrleIndexType, D>(cellIt.getIndices()))
          nodes[u].erase(nodes[u].begin());
      }

      unsigned signs = 0;
      for (int i = 0; i < (1 << D); i++) {
        if (cellIt.getCorner(i).getValue() >= T(0))
          signs |= (1 << i);
      }

      // all corners have the same sign, so no surface here
      if (signs == 0)
        continue;
      if (signs == (1 << (1 << D)) - 1)
        continue;

      // for each element
      const int *Triangles = (D == 2) ? marchingCubes.polygonize2d(signs)
                                      : marchingCubes.polygonize3d(signs);

      for (; Triangles[0] != -1; Triangles += D) {
        std::array<unsigned, D> nod_numbers;

        // for each node
        for (int n = 0; n < D; n++) {
          const int edge = Triangles[n];

          unsigned p0 = corner0[edge];
          unsigned p1 = corner1[edge];

          // determine direction of edge
          int dir = direction[edge];

          // look for existing surface node
          hrleVectorType<hrleIndexType, D> d(cellIt.getIndices());
          d += BitMaskToVector<D, hrleIndexType>(p0);

          nodeIt = nodes[dir].find(d);
          if (nodeIt != nodes[dir].end()) {
            nod_numbers[n] = nodeIt->second;
          } else { // if node does not exist yet

            // calculate coordinate of new node
            std::array<T, 3> cc{}; // initialise with zeros
            for (int z = 0; z < D; z++) {
              if (z != dir) {
                // TODO might not need BitMaskToVector here, just check if z bit
                // is set
                cc[z] = double(cellIt.getIndices(z) +
                               BitMaskToVector<D, hrleIndexType>(p0)[z]);
              } else {
                T d0, d1;

                d0 = cellIt.getCorner(p0).getValue();
                d1 = cellIt.getCorner(p1).getValue();

                // calculate the surface-grid intersection point
                if (d0 == -d1) { // includes case where d0=d1=0
                  // meshNodeToPointIdMapping.push_back(
                  //     cellIt.getCorner(p0).getPointId());
                  cc[z] = static_cast<T>(cellIt.getIndices(z)) + 0.5;
                } else {
                  if (std::abs(d0) <= std::abs(d1)) {
                    // meshNodeToPointIdMapping.push_back(
                    //     cellIt.getCorner(p0).getPointId());
                    cc[z] =
                        static_cast<T>(cellIt.getIndices(z)) + (d0 / (d0 - d1));
                  } else {
                    // meshNodeToPointIdMapping.push_back(
                    //     cellIt.getCorner(p1).getPointId());
                    cc[z] = static_cast<T>(cellIt.getIndices(z) + 1) -
                            (d1 / (d1 - d0));
                  }
                }
                cc[z] = std::max(cc[z], cellIt.getIndices(z) + epsilon);
                cc[z] = std::min(cc[z], (cellIt.getIndices(z) + 1) - epsilon);
              }
              cc[z] = levelSet->getGrid().getGridDelta() * cc[z];
            }

            // insert new node
            nod_numbers[n] =
                mesh->insertNextNode(cc); // insert new surface node
            nodes[dir][d] = nod_numbers[n];
          }
        }

        mesh->insertNextElement(nod_numbers); // insert new surface element
      }
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsToSurfaceMesh)

#endif // LS_TO_SURFACE_MESH_HPP
