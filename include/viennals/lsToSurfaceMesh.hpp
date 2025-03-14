#pragma once

#include <lsPreCompileMacros.hpp>

#include <iostream>
#include <map>

#include <hrleSparseCellIterator.hpp>
#include <lsDomain.hpp>
#include <lsMarchingCubes.hpp>
#include <lsMesh.hpp>

namespace viennals {

using namespace viennacore;

/// Extract an explicit Mesh<> instance from an lsDomain.
/// The interface is then described by explicit surface elements:
/// Lines in 2D, Triangles in 3D.
template <class T, int D> class ToSurfaceMesh {
  typedef typename Domain<T, D>::DomainType hrleDomainType;

  SmartPointer<Domain<T, D>> levelSet = nullptr;
  SmartPointer<Mesh<T>> mesh = nullptr;
  // std::vector<hrleIndexType> meshNodeToPointIdMapping;
  const T epsilon;
  bool updatePointData = true;

public:
  explicit ToSurfaceMesh(double eps = 1e-12) : epsilon(eps) {}

  ToSurfaceMesh(const SmartPointer<Domain<T, D>> passedLevelSet,
                SmartPointer<Mesh<T>> passedMesh, double eps = 1e-12)
      : levelSet(passedLevelSet), mesh(passedMesh), epsilon(eps) {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  void setMesh(SmartPointer<Mesh<T>> passedMesh) { mesh = passedMesh; }

  void setUpdatePointData(bool update) { updatePointData = update; }

  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addWarning("No level set was passed to ToSurfaceMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      Logger::getInstance()
          .addWarning("No mesh was passed to ToSurfaceMesh.")
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
      Logger::getInstance()
          .addWarning("Levelset is less than 2 layers wide. Export might fail!")
          .print();
    }

    typedef std::map<hrleVectorType<hrleIndexType, D>, unsigned>
        nodeContainerType;

    nodeContainerType nodes[D];
    typename nodeContainerType::iterator nodeIt;
    const bool updateData = updatePointData;

    // save how data should be transferred to new level set
    // list of indices into the old pointData vector
    std::vector<std::vector<unsigned>> newDataSourceIds;
    // there is no multithreading here, so just use 1
    if (updateData)
      newDataSourceIds.resize(1);

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
      const int *Triangles =
          (D == 2) ? lsInternal::MarchingCubes::polygonize2d(signs)
                   : lsInternal::MarchingCubes::polygonize3d(signs);

      for (; Triangles[0] != -1; Triangles += D) {
        std::array<unsigned, D> nod_numbers;

        // for each node
        for (int n = 0; n < D; n++) {
          const int edge = Triangles[n];

          unsigned p0 = corner0[edge];
          unsigned p1 = corner1[edge];

          // determine direction of edge
          auto dir = direction[edge];

          // look for existing surface node
          hrleVectorType<hrleIndexType, D> d(cellIt.getIndices());
          d += hrleUtil::BitMaskToVector<D, hrleIndexType>(p0);

          nodeIt = nodes[dir].find(d);
          if (nodeIt != nodes[dir].end()) {
            nod_numbers[n] = nodeIt->second;
          } else { // if node does not exist yet

            // calculate coordinate of new node
            std::array<T, 3> cc{}; // initialise with zeros
            std::size_t currentPointId = 0;
            for (int z = 0; z < D; z++) {
              if (z != dir) {
                // TODO might not need BitMaskToVector here, just check if z bit
                // is set
                cc[z] = static_cast<double>(
                    cellIt.getIndices(z) +
                    hrleUtil::BitMaskToVector<D, hrleIndexType>(p0)[z]);
              } else {
                T d0, d1;

                d0 = cellIt.getCorner(p0).getValue();
                d1 = cellIt.getCorner(p1).getValue();

                // calculate the surface-grid intersection point
                if (d0 == -d1) { // includes case where d0=d1=0
                  currentPointId = cellIt.getCorner(p0).getPointId();
                  cc[z] = static_cast<T>(cellIt.getIndices(z)) + 0.5;
                } else {
                  if (std::abs(d0) <= std::abs(d1)) {
                    currentPointId = cellIt.getCorner(p0).getPointId();
                    cc[z] =
                        static_cast<T>(cellIt.getIndices(z)) + (d0 / (d0 - d1));
                  } else {
                    currentPointId = cellIt.getCorner(p1).getPointId();
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

            if (updateData)
              newDataSourceIds[0].push_back(currentPointId);
          }
        }

        if (!triangleMisformed(nod_numbers))
          mesh->insertNextElement(nod_numbers); // insert new surface element
      }
    }

    // now copy old data into new level set
    if (updateData) {
      mesh->getPointData().translateFromMultiData(levelSet->getPointData(),
                                                  newDataSourceIds);
    }
  }

private:
  static bool triangleMisformed(const std::array<unsigned, D> &nodeNumbers) {
    if constexpr (D == 3) {
      return nodeNumbers[0] == nodeNumbers[1] ||
             nodeNumbers[0] == nodeNumbers[2] ||
             nodeNumbers[1] == nodeNumbers[2];
    } else {
      return nodeNumbers[0] == nodeNumbers[1];
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(ToSurfaceMesh)

} // namespace viennals
