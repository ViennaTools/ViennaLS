#ifndef LS_EXTRUDE_HPP
#define LS_EXTRUDE_HPP

#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsToSurfaceMesh.hpp>

// Extrude a 2D Level Set into a 3D domain.
template <class T> class lsExtrude {
  lsSmartPointer<lsDomain<T, 2>> inputLevelSet = nullptr;
  lsSmartPointer<lsDomain<T, 3>> outputLevelSet = nullptr;
  std::array<T, 2> extent = {0., 0.};
  int extrudeDim = 0;

public:
  lsExtrude() {}
  lsExtrude(lsSmartPointer<lsDomain<T, 2>> passedInputLS,
            lsSmartPointer<lsDomain<T, 3>> passedOutputLS,
            std::array<T, 2> passedExtent, const int passedExtrudeDim = 0)
      : inputLevelSet(passedInputLS), outputLevelSet(passedOutputLS),
        extent(passedExtent), extrudeDim(passedExtrudeDim) {}

  void setInputLevelSet(lsSmartPointer<lsDomain<T, 2>> passedInputLS) {
    inputLevelSet = passedInputLS;
  }

  // The 3D output LS will be overwritten by the extruded LS
  void setOutputLevelSet(lsSmartPointer<lsDomain<T, 3>> &passedOutputLS) {
    outputLevelSet = passedOutputLS;
  }

  // Set the min and max extent in the extruded dimension
  void setExtent(std::array<T, 2> passedExtent) { extent = passedExtent; }

  // Set which index of the added dimension (x: 0, y: 1, z: 2)
  void setExtrudeDimension(const int passedExtrudeDim) {
    extrudeDim = passedExtrudeDim;
  }

  void apply() {
    if (inputLevelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning(
              "No input Level Set supplied to lsExtrude! Not converting.")
          .print();
    }
    if (outputLevelSet == nullptr) {
      lsMessage::getInstance()
          .addWarning(
              "No output Level Set supplied to lsExtrude! Not converting.")
          .print();
      return;
    }

    // x and y of the input LS get transformed to these indices
    const auto extrudeDims = getExtrudeDims();

    // create new domain based on 2D extent and boundary conditions
    {
      const T gridDelta = inputLevelSet->getGrid().getGridDelta();
      const auto inputBoundaryConds =
          inputLevelSet->getGrid().getBoundaryConditions();
      auto minBounds = inputLevelSet->getGrid().getMinBounds();
      auto maxBounds = inputLevelSet->getGrid().getMaxBounds();

      double domainBounds[2 * 3];
      domainBounds[2 * extrudeDim] = extent[0];
      domainBounds[2 * extrudeDim + 1] = extent[1];
      domainBounds[2 * extrudeDims[0]] = gridDelta * minBounds[0];
      domainBounds[2 * extrudeDims[0] + 1] = gridDelta * maxBounds[0];
      domainBounds[2 * extrudeDims[1]] = gridDelta * minBounds[1];
      domainBounds[2 * extrudeDims[1] + 1] = gridDelta * maxBounds[1];

      lsBoundaryConditionEnum<3> boundaryConds[3];
      boundaryConds[extrudeDim] =
          convertBoundaryCondition(inputBoundaryConds[extrudeDims[0]]);
      boundaryConds[extrudeDims[0]] =
          convertBoundaryCondition(inputBoundaryConds[extrudeDims[0]]);
      boundaryConds[extrudeDims[1]] =
          convertBoundaryCondition(inputBoundaryConds[extrudeDims[1]]);

      auto tmpLevelSet = lsSmartPointer<lsDomain<T, 3>>::New(
          domainBounds, boundaryConds, gridDelta);
      outputLevelSet->deepCopy(tmpLevelSet);
    }

    auto surface = lsSmartPointer<lsMesh<T>>::New();
    lsToSurfaceMesh<T, 2>(inputLevelSet, surface).apply();

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
      std::array<unsigned, 3> triangle = {lines[i][0], lines[i][1],
                                          lines[i][0] + numNodes};
      surface->insertNextTriangle(triangle);
      triangle[0] = lines[i][1];
      triangle[1] = lines[i][1] + numNodes;
      triangle[2] = lines[i][0] + numNodes;
      surface->insertNextTriangle(triangle);
    }
    surface->template getElements<2>().clear(); // remove lines

    lsFromSurfaceMesh<T, 3>(outputLevelSet, surface).apply();
  }

private:
  lsBoundaryConditionEnum<3>
  convertBoundaryCondition(lsBoundaryConditionEnum<2> boundaryCond) {
    switch (boundaryCond) {
    case lsBoundaryConditionEnum<2>::PERIODIC_BOUNDARY:
      return lsBoundaryConditionEnum<3>::PERIODIC_BOUNDARY;
    case lsBoundaryConditionEnum<2>::INFINITE_BOUNDARY:
      return lsBoundaryConditionEnum<3>::INFINITE_BOUNDARY;
    case lsBoundaryConditionEnum<2>::POS_INFINITE_BOUNDARY:
      return lsBoundaryConditionEnum<3>::POS_INFINITE_BOUNDARY;
    case lsBoundaryConditionEnum<2>::NEG_INFINITE_BOUNDARY:
      return lsBoundaryConditionEnum<3>::NEG_INFINITE_BOUNDARY;
    case lsBoundaryConditionEnum<2>::REFLECTIVE_BOUNDARY:
      return lsBoundaryConditionEnum<3>::REFLECTIVE_BOUNDARY;
    default:
      return lsBoundaryConditionEnum<3>::INFINITE_BOUNDARY;
    }
  }

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

#endif // LS_EXTRUDE_HPP