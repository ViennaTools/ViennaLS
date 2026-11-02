#pragma once

#include <lsPreCompileMacros.hpp>

#include <lsDomain.hpp>
#include <lsMesh.hpp>

namespace viennals {

using namespace viennacore;

/// Import the regular grid, on which the level set values are
/// defined, from an explicit Mesh<>. The Vertices must be defined,
/// as well as a scalar data field "LSValues" holding the signed
/// distance in grid spacings. Points further than 2.5 grid spacings
/// from the surface are ignored. The points are sorted before they are
/// inserted, unless setSortPointList(false) is used, in which case they
/// must already be sorted lexicographically, last coordinate first.
template <class T, int D> class FromMesh {
  SmartPointer<Domain<T, D>> levelSet = nullptr;
  SmartPointer<Mesh<T>> mesh = nullptr;
  bool sortPointList = true;

public:
  FromMesh() = default;

  FromMesh(SmartPointer<Domain<T, D>> passedLevelSet,
           SmartPointer<Mesh<T>> passedMesh)
      : levelSet(passedLevelSet), mesh(passedMesh) {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedlsDomain) {
    levelSet = passedlsDomain;
  }

  void setMesh(const SmartPointer<Mesh<T>> passedMesh) { mesh = passedMesh; }

  void setSortPointList(bool passedSortPointList) {
    sortPointList = passedSortPointList;
  }

  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addError("No level set was passed to FromMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      Logger::getInstance()
          .addError("No Mesh was supplied to FromMesh.")
          .print();
      return;
    }

    auto &nodes = mesh->getNodes();
    auto values = mesh->getPointData().getScalarData("LSValues", true);
    if (values == nullptr) {
      values = mesh->getCellData().getScalarData("LSValues", true);
    }

    if (values == nullptr) {
      Logger::getInstance()
          .addError("Mesh does not contain level set values (\"LSValues\").")
          .print();
      return;
    }

    const T gridDelta = levelSet->getGrid().getGridDelta();
    const std::size_t numberOfPoints = std::min(nodes.size(), values->size());

    typename Domain<T, D>::PointValueVectorType points;
    points.reserve(numberOfPoints);
    for (std::size_t i = 0; i < numberOfPoints; ++i) {
      // only read in points within the first 5 layers, to ignore
      // undefined points
      if (std::abs((*values)[i]) > 2.5)
        continue;

      viennahrle::Index<D> index;
      for (unsigned j = 0; j < D; ++j) {
        index[j] = std::round(nodes[i][j] / gridDelta);
      }
      points.emplace_back(index, (*values)[i]);
    }

    levelSet->insertPoints(std::move(points), sortPointList);
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(FromMesh)

} // namespace viennals
