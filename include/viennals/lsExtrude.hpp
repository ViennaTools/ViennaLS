#pragma once

#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsToSurfaceMesh.hpp>

namespace viennals {

using namespace viennacore;

/// Extrudes a 2D Level Set into a 3D domain. The axis in which should be
/// extruded can be set and boundary conditions in the 3D domain must be
/// specified.
template <class T> class Extrude {
  SmartPointer<Domain<T, 2>> inputLevelSet = nullptr;
  SmartPointer<Domain<T, 3>> outputLevelSet = nullptr;
  Vec2D<T> extent = {0., 0.};
  int extrudeDim = 0;
  std::array<BoundaryConditionEnum<3>, 3> boundaryConds;

public:
  Extrude() {}
  Extrude(SmartPointer<Domain<T, 2>> passedInputLS,
          SmartPointer<Domain<T, 3>> passedOutputLS, Vec2D<T> passedExtent,
          const int passedExtrudeDim,
          std::array<BoundaryConditionEnum<3>, 3> passedBoundaryConds)
      : inputLevelSet(passedInputLS), outputLevelSet(passedOutputLS),
        extent(passedExtent), extrudeDim(passedExtrudeDim),
        boundaryConds(passedBoundaryConds) {}

  void setInputLevelSet(SmartPointer<Domain<T, 2>> passedInputLS) {
    inputLevelSet = passedInputLS;
  }

  // The 3D output LS will be overwritten by the extruded LS
  void setOutputLevelSet(SmartPointer<Domain<T, 3>> &passedOutputLS) {
    outputLevelSet = passedOutputLS;
  }

  // Set the min and max extent in the extruded dimension
  void setExtent(Vec2D<T> passedExtent) { extent = passedExtent; }

  // Set which index of the added dimension (x: 0, y: 1, z: 2)
  void setExtrudeDimension(const int passedExtrudeDim) {
    extrudeDim = passedExtrudeDim;
  }

  void setBoundaryConditions(
      std::array<BoundaryConditionEnum<3>, 3> passedBoundaryConds) {
    boundaryConds = passedBoundaryConds;
  }

  void setBoundaryConditions(BoundaryConditionEnum<3> passedBoundaryConds[3]) {
    for (int i = 0; i < 3; i++)
      boundaryConds[i] = passedBoundaryConds[i];
  }

  void apply() {
    if (inputLevelSet == nullptr) {
      Logger::getInstance()
          .addWarning("No input Level Set supplied to Extrude! Not converting.")
          .print();
    }
    if (outputLevelSet == nullptr) {
      Logger::getInstance()
          .addWarning(
              "No output Level Set supplied to Extrude! Not converting.")
          .print();
      return;
    }

    // x and y of the input LS get transformed to these indices
    const auto extrudeDims = getExtrudeDims();

    // create new domain based on 2D extent
    {
      const T gridDelta = inputLevelSet->getGrid().getGridDelta();
      auto minBounds = inputLevelSet->getGrid().getMinBounds();
      auto maxBounds = inputLevelSet->getGrid().getMaxBounds();

      double domainBounds[2 * 3];
      domainBounds[2 * extrudeDim] = extent[0];
      domainBounds[2 * extrudeDim + 1] = extent[1];
      domainBounds[2 * extrudeDims[0]] = gridDelta * minBounds[0];
      domainBounds[2 * extrudeDims[0] + 1] = gridDelta * maxBounds[0];
      domainBounds[2 * extrudeDims[1]] = gridDelta * minBounds[1];
      domainBounds[2 * extrudeDims[1] + 1] = gridDelta * maxBounds[1];

      auto tmpLevelSet = SmartPointer<Domain<T, 3>>::New(
          domainBounds, boundaryConds.data(), gridDelta);
      outputLevelSet->deepCopy(tmpLevelSet);
    }

    auto surface = SmartPointer<Mesh<T>>::New();
    ToSurfaceMesh<T, 2>(inputLevelSet, surface).apply();

    auto &lines = surface->template getElements<2>();
    auto &nodes = surface->getNodes();
    const unsigned numNodes = nodes.size();

    // add new nodes shifted by the extent
    for (unsigned i = 0; i < numNodes; i++) {
      nodes[i][extrudeDims[1]] = nodes[i][1];
      nodes[i][extrudeDims[0]] = nodes[i][0];
      nodes[i][extrudeDim] = extent[1];

      nodes.push_back(nodes[i]);

      nodes[i][extrudeDim] = extent[0];
    }

    // add triangles in places of lines
    for (unsigned i = 0; i < lines.size(); i++) {
      std::array<unsigned, 3> triangle = {lines[i][1], lines[i][0],
                                          lines[i][0] + numNodes};
      if (extrudeDim == 1) {
        std::swap(triangle[0], triangle[2]);
      }
      surface->insertNextTriangle(triangle);
      triangle[0] = lines[i][0] + numNodes;
      triangle[1] = lines[i][1] + numNodes;
      triangle[2] = lines[i][1];
      if (extrudeDim == 1) {
        std::swap(triangle[0], triangle[2]);
      }
      surface->insertNextTriangle(triangle);
    }
    surface->template getElements<2>().clear(); // remove lines

    FromSurfaceMesh<T, 3>(outputLevelSet, surface).apply();
  }

private:
  inline std::array<unsigned, 2> getExtrudeDims() const {
    assert(extrudeDim < 3);
    if (extrudeDim == 0) {
      return {1, 2};
    } else if (extrudeDim == 1) {
      return {0, 2};
    } else {
      return {0, 1};
    }
  }
};

} // namespace viennals
