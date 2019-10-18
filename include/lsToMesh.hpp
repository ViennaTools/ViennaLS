#ifndef LS_TO_MESH_HPP
#define LS_TO_MESH_HPP

#include <lsPreCompileMacros.hpp>

#include <iostream>
#include <map>

#include <hrleSparseIterator.hpp>
#include <lsDomain.hpp>
#include <lsMesh.hpp>

template <class T, int D> class lsToMesh {
  typedef typename lsDomain<T, D>::DomainType hrleDomainType;

  const lsDomain<T, D> &levelSet;
  lsMesh &mesh;
  const bool onlyDefined;
  const bool onlyActive;

  lsToMesh();

public:
  lsToMesh(const lsDomain<T, D> &passedLevelSet, lsMesh &passedMesh,
           bool passedOnlyDefined = true, bool passedOnlyActive = false)
      : levelSet(passedLevelSet), mesh(passedMesh),
        onlyDefined(passedOnlyDefined), onlyActive(passedOnlyActive) {}

  void apply() {
    mesh.clear();

    std::vector<double> scalarData;
    std::vector<double> subLS;

    for (hrleConstSparseIterator<hrleDomainType> it(levelSet.getDomain());
         !it.isFinished(); ++it) {
      if ((onlyDefined && !it.isDefined()) ||
          (onlyActive && std::abs(it.getValue()) > 0.5))
        continue;

      // insert vertex
      hrleVectorType<unsigned, 1> vertex;
      vertex[0] = mesh.nodes.size();
      mesh.insertNextVertex(vertex);

      // insert corresponding node
      hrleVectorType<double, 3> node;
      node[2] = 0.;
      T gridDelta = levelSet.getGrid().getGridDelta();
      for (unsigned i = 0; i < D; ++i) {
        node[i] = double(it.getStartIndices(i)) * gridDelta;
      }
      mesh.insertNextNode(node);

      // insert LS value
      scalarData.push_back(it.getValue());
      subLS.push_back(it.getSegmentId());
    }

    mesh.insertNextScalarData(scalarData, "LSValues");
    mesh.insertNextScalarData(subLS, "SegmentID");
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(lsToMesh)

#endif // LS_TO_MESH_HPP
