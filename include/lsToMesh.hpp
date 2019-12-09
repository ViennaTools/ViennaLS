#ifndef LS_TO_MESH_HPP
#define LS_TO_MESH_HPP

#include <lsPreCompileMacros.hpp>

#include <iostream>

#include <hrleSparseIterator.hpp>
#include <lsDomain.hpp>
#include <lsMesh.hpp>

/// Extract the regular grid, on which the level set values are
/// defined, to an explicit lsMesh. The Vertices will contain
/// the level set value stored at its location. (This is very useful
/// for debugging)
template <class T, int D> class lsToMesh {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;

  const lsDomain<T, D> *levelSet = nullptr;
  lsMesh *mesh = nullptr;
  bool onlyDefined;
  bool onlyActive;

public:
  lsToMesh(){};

  lsToMesh(const lsDomain<T, D> &passedLevelSet, lsMesh &passedMesh,
           bool passedOnlyDefined = true, bool passedOnlyActive = false)
      : levelSet(&passedLevelSet), mesh(&passedMesh),
        onlyDefined(passedOnlyDefined), onlyActive(passedOnlyActive) {}

  void setLevelSet(lsDomain<T, D> &passedlsDomain) {
    levelSet = &passedlsDomain;
  }

  void setMesh(lsMesh &passedMesh) { mesh = &passedMesh; }

  void setOnlyDefined(bool passedOnlyDefined) {
    onlyDefined = passedOnlyDefined;
  }

  void setOnlyActive(bool passedOnlyActive) { onlyActive = passedOnlyActive; }

  void apply() {
    if (levelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning("No level set was passed to lsToMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      lsMessage::getInstance()
          .addWarning("No mesh was passed to lsToMesh.")
          .print();
      return;
    }

    mesh->clear();

    std::vector<double> scalarData;
    std::vector<double> subLS;

    const T gridDelta = levelSet->getGrid().getGridDelta();

    for (hrleConstSparseIterator<hrleDomainType> it(levelSet->getDomain());
         !it.isFinished(); ++it) {
      if ((onlyDefined && !it.isDefined()) ||
          (onlyActive && std::abs(it.getValue()) > 0.5))
        continue;

      // insert vertex
      std::array<unsigned, 1> vertex;
      vertex[0] = mesh->nodes.size();
      mesh->insertNextVertex(vertex);

      // insert corresponding node
      std::array<double, 3> node;
      if (D == 2)
        node[2] = 0.;
      for (unsigned i = 0; i < D; ++i) {
        node[i] = double(it.getStartIndices(i)) * gridDelta;
      }
      mesh->insertNextNode(node);

      // insert LS value
      scalarData.push_back(it.getValue());
      subLS.push_back(it.getSegmentId());
    }

    mesh->insertNextScalarData(scalarData, "LSValues");
    mesh->insertNextScalarData(subLS, "SegmentID");
    if (levelSet->getPointData().getScalarDataSize() > 0 ||
        levelSet->getPointData().getVectorDataSize() > 0) {
      mesh->append(levelSet->getPointData());
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsToMesh)

#endif // LS_TO_MESH_HPP
