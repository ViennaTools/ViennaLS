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
  int extrusionAxis = 0;
  std::array<BoundaryConditionEnum, 3> boundaryConds{};

public:
  Extrude() = default;

  Extrude(SmartPointer<Domain<T, 2>> passedInputLS,
          SmartPointer<Domain<T, 3>> passedOutputLS, Vec2D<T> passedExtent,
          int passedExtrusionAxis,
          BoundaryConditionEnum passedBoundaryCond =
              BoundaryConditionEnum::INFINITE_BOUNDARY)
      : inputLevelSet(passedInputLS), outputLevelSet(passedOutputLS),
        extent(passedExtent), extrusionAxis(passedExtrusionAxis),
        boundaryConds{{passedInputLS->getGrid().getBoundaryConditions(0),
                       passedInputLS->getGrid().getBoundaryConditions(1),
                       passedBoundaryCond}} {}

  Extrude(SmartPointer<Domain<T, 2>> passedInputLS,
          SmartPointer<Domain<T, 3>> passedOutputLS, Vec2D<T> passedExtent,
          int passedExtrusionAxis,
          std::array<BoundaryConditionEnum, 3> passedBoundaryConds)
      : inputLevelSet(passedInputLS), outputLevelSet(passedOutputLS),
        extent(passedExtent), extrusionAxis(passedExtrusionAxis),
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

  void setExtrusionAxis(int passedExtrusionAxis) {
    extrusionAxis = passedExtrusionAxis;
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
          .addError("No input Level Set supplied to Extrude.")
          .print();
    }
    if (outputLevelSet == nullptr) {
      Logger::getInstance()
          .addError("No output Level Set supplied to Extrude.")
          .print();
      return;
    }
    if (extrusionAxis < 0 || extrusionAxis > 2) {
      Logger::getInstance()
          .addError("Extrusion axis must be between 0 and 2.")
          .print();
      return;
    }

    std::vector<std::pair<viennahrle::Index<3>, T>> points3D;

    auto &domain2D = inputLevelSet->getDomain();
    auto &grid2D = inputLevelSet->getGrid();
    const T gridDelta = grid2D.getGridDelta();
    auto minBounds = grid2D.getMinBounds();
    auto maxBounds = grid2D.getMaxBounds();

    double domainBounds[2 * 3];
    unsigned xx, yy;
    if (extrusionAxis == 0) {
      domainBounds[0] = extent[0];
      domainBounds[1] = extent[1];
      domainBounds[2] = minBounds[0] * gridDelta;
      domainBounds[3] = maxBounds[0] * gridDelta;
      domainBounds[4] = minBounds[1] * gridDelta;
      domainBounds[5] = maxBounds[1] * gridDelta;
      xx = 1;
      yy = 2;
    } else if (extrusionAxis == 1) {
      domainBounds[0] = minBounds[0] * gridDelta;
      domainBounds[1] = maxBounds[0] * gridDelta;
      domainBounds[2] = extent[0];
      domainBounds[3] = extent[1];
      domainBounds[4] = minBounds[1] * gridDelta;
      domainBounds[5] = maxBounds[1] * gridDelta;
      xx = 0;
      yy = 2;
    } else if (extrusionAxis == 2) {
      domainBounds[0] = minBounds[0] * gridDelta;
      domainBounds[1] = maxBounds[0] * gridDelta;
      domainBounds[2] = minBounds[1] * gridDelta;
      domainBounds[3] = maxBounds[1] * gridDelta;
      domainBounds[4] = extent[0];
      domainBounds[5] = extent[1];
      xx = 0;
      yy = 1;
    }

    auto tmpLevelSet = SmartPointer<Domain<T, 3>>::New(
        domainBounds, boundaryConds.data(), gridDelta);
    outputLevelSet->deepCopy(tmpLevelSet);

    const hrleIndexType extStart = std::floor(extent[0] / gridDelta);
    const hrleIndexType extEnd = std::ceil(extent[1] / gridDelta);

    for (viennahrle::SparseIterator<typename Domain<T, 2>::DomainType> it(
             domain2D);
         !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;

      const auto index2D = it.getStartIndices();
      const T value = it.getValue();

      for (hrleIndexType ext = extStart; ext <= extEnd; ++ext) {
        viennahrle::Index<3> index3D;
        index3D[extrusionAxis] = ext;
        index3D[xx] = index2D[0];
        index3D[yy] = index2D[1];

        points3D.emplace_back(index3D, value);
      }
    }

    outputLevelSet->insertPoints(points3D);
    outputLevelSet->finalize(2);
  }
};

PRECOMPILE_PRECISION(Extrude);

} // namespace viennals
