#pragma once

#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsToSurfaceMesh.hpp>

#include <hrleSparseIterator.hpp>

namespace viennals {

using namespace viennacore;

/// Extrudes a 2D Level Set into a 3D domain. The axis in which should be
/// extruded can be set and boundary conditions in the 3D domain must be
/// specified.
template <class T> class Extrude {
  using hrleIndexType = viennahrle::IndexType;
  SmartPointer<Domain<T, 2>> inputLevelSet = nullptr;
  SmartPointer<Domain<T, 3>> outputLevelSet = nullptr;
  Vec2D<T> extent = {0., 0.};
  int extrudeDim = 0;
  std::array<BoundaryConditionEnum, 3> boundaryConds = {};

public:
  Extrude() = default;

  Extrude(SmartPointer<Domain<T, 2>> passedInputLS,
          SmartPointer<Domain<T, 3>> passedOutputLS, Vec2D<T> passedExtent,
          const int passedExtrudeDim = 2,
          BoundaryConditionEnum passedBoundaryCond =
              BoundaryConditionEnum::INFINITE_BOUNDARY)
      : inputLevelSet(passedInputLS), outputLevelSet(passedOutputLS),
        extent(passedExtent), extrudeDim(passedExtrudeDim),
        boundaryConds{{passedInputLS->getGrid().getBoundaryConditions(0),
                       passedInputLS->getGrid().getBoundaryConditions(1),
                       passedBoundaryCond}} {}

  Extrude(SmartPointer<Domain<T, 2>> passedInputLS,
          SmartPointer<Domain<T, 3>> passedOutputLS, Vec2D<T> passedExtent,
          const int passedExtrudeDim,
          std::array<BoundaryConditionEnum, 3> passedBoundaryConds)
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
      std::array<BoundaryConditionEnum, 3> passedBoundaryConds) {
    boundaryConds = passedBoundaryConds;
  }

  void setBoundaryConditions(BoundaryConditionEnum passedBoundaryConds[3]) {
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

    {
      std::vector<std::pair<viennahrle::Index<3>, T>> points3D;

      auto &domain2D = inputLevelSet->getDomain();
      auto &grid2D = inputLevelSet->getGrid();
      const T gridDelta = grid2D.getGridDelta();
      auto minBounds = grid2D.getMinBounds();
      auto maxBounds = grid2D.getMaxBounds();

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

      for (viennahrle::SparseIterator<typename Domain<T, 2>::DomainType> it(
               domain2D);
           !it.isFinished(); ++it) {
        if (!it.isDefined())
          continue;

        const auto index2D = it.getStartIndices();
        const T value = it.getValue();

        const hrleIndexType zStart = std::floor(extent[0] / gridDelta) - 1;
        const hrleIndexType zEnd = std::ceil(extent[1] / gridDelta) + 1;
        for (hrleIndexType z = zStart; z <= zEnd; ++z) {
          viennahrle::Index<3> index3D;
          index3D[0] = index2D[0];
          index3D[1] = index2D[1];
          index3D[2] = z;

          points3D.emplace_back(index3D, value);
        }
      }

      outputLevelSet->insertPoints(points3D);
      outputLevelSet->finalize(2);
    }
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
