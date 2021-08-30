#include <iostream>
#include <numeric>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsTestAsserts.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Check if Boolean Operation works with exact zero LS values
  \example BooleanOperationExactZero.cpp
*/

using NumericType = double;
constexpr int D = 2;

using LSType = lsSmartPointer<lsDomain<NumericType, D>>;

void writeLS(LSType levelSet, std::string fileName) {
  auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();
  lsToMesh(levelSet, mesh).apply();
  lsVTKWriter(mesh, fileName).apply();
}

int main() {
  omp_set_num_threads(1);

  LSType mask;
  {
    double gridDelta = 1.0;
    double extent = 10;
    double bounds[2 * D] = {-extent, extent, -extent, extent};

    typename lsDomain<NumericType, D>::BoundaryType boundaryCons[D];
    boundaryCons[0] =
        lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
    boundaryCons[1] = lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

    if constexpr (D == 3) {
      boundaryCons[1] =
          lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
      boundaryCons[2] =
          lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;
    }

    mask = lsSmartPointer<lsDomain<NumericType, D>>::New(bounds, boundaryCons,
                                                         gridDelta);
  }

  // make mask geometry ( just a simple box )
  {
    double min[D] = {-5., 1.};
    double max[D] = {5., 10.};
    if constexpr (D == 3) {
      min[1] = -5.;
      min[2] = 1.;
      max[1] = 5.;
      max[2] = 10.;
    }
    lsMakeGeometry(mask, lsSmartPointer<lsBox<NumericType, D>>::New(min, max))
        .apply();
    writeLS(mask, "mask_initial.vtp");
  }

  // now make substrate ( plane ) at the same height as the bottom of the mask
  auto substrate = LSType::New(mask->getGrid());
  {
    double origin[D] = {};
    double normal[D] = {};
    origin[D - 1] = 1.;
    normal[D - 1] = 1.;
    lsMakeGeometry(substrate,
                   lsSmartPointer<lsPlane<NumericType, D>>::New(origin, normal))
        .apply();
    writeLS(substrate, "subs_initial.vtp");
  }

  lsBooleanOperation(substrate, mask, lsBooleanOperationEnum::UNION).apply();

  writeLS(substrate, "subs_afterBool.vtp");

  return 0;
}
