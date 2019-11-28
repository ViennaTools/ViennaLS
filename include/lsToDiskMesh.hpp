#ifndef LS_TO_DISK_MESH_HPP
#define LS_TO_DISK_MESH_HPP

#include <lsPreCompileMacros.hpp>

#include <hrleSparseIterator.hpp>

#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>

/// This class creates a mesh from the level set
/// with all grid points with a level set value <= 0.5.
/// These grid points are shifted in space towards the
/// direction of their normal vector by grid delta * LS value.
/// Grid delta and the origin grid point are saved for each point.
/// This allows for a simple setup of disks for ray tracing.
template <class T, int D> class lsToDiskMesh {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;

  lsDomain<T, D> *levelSet = nullptr;
  lsMesh *mesh = nullptr;

public:
  lsToDiskMesh() {}

  lsToDiskMesh(lsDomain<T, D> &passedLevelSet, lsMesh &passedMesh)
      : levelSet(&passedLevelSet), mesh(&passedMesh) {}

  void setLevelSet(lsDomain<T, D> &passedLevelSet) {
    levelSet = &passedLevelSet;
  }

  void setMesh(lsMesh &passedMesh) { mesh = &passedMesh; }

  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsToDiskMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsToDiskMesh.")
          .print();
      return;
    }

    mesh->clear();

    lsExpand<T, D>(*levelSet, 3).apply();
    lsCalculateNormalVectors<T, D>(*levelSet, true).apply();

    const T gridDelta = levelSet->getGrid().getGridDelta();
    const auto &normalVectors = levelSet->getNormalVectors();

    // set up data arrays
    std::vector<double> values(normalVectors.size());
    std::vector<double> gridSpacing(normalVectors.size());
    std::vector<std::array<double, 3>> normals(normalVectors.size());

    unsigned pointId = 0;

    for (hrleConstSparseIterator<hrleDomainType> it(levelSet->getDomain());
         !it.isFinished(); ++it) {
      if (!it.isDefined() || (std::abs(it.getValue()) > 0.5))
        continue;

      // insert vertex
      std::array<unsigned, 1> vertex;
      vertex[0] = mesh->nodes.size();
      mesh->insertNextVertex(vertex);

      // insert corresponding node shifted by ls value in direction of the
      // normal vector
      std::array<double, 3> node;
      node[2] = 0.;
      for (unsigned i = 0; i < D; ++i) {
        // original position
        node[i] = double(it.getStartIndices(i)) * gridDelta;
        // shift position
        node[i] -= it.getValue() * gridDelta * normalVectors[pointId][i];
      }
      mesh->insertNextNode(node);

      // add data into mesh
      values[pointId] = it.getValue();
      gridSpacing[pointId] = gridDelta;
      // copy normal
      if (D == 2)
        normals[pointId][2] = 0.;
      for (unsigned i = 0; i < D; ++i) {
        normals[pointId][i] = normalVectors[pointId][i];
      }

      ++pointId;
    }

    mesh->insertNextScalarData(values, "LSValues");
    mesh->insertNextScalarData(gridSpacing, "gridSpacing");
    mesh->insertNextVectorData(normals, "Normals");
  }
};

#endif // LS_TO_DISK_MESH_HPP
