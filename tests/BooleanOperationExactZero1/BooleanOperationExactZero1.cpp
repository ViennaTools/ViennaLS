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
  Check if Boolean Operation works with exact zero LS values and only slightly
  positive numbers \example BooleanOperationExactZero1.cpp
*/

namespace ls = viennals;

using NumericType = double;
constexpr int D = 2;

using LSType = ls::SmartPointer<ls::Domain<NumericType, D>>;

// void writeLS(LSType levelSet, std::string fileName) {
//   auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
//   lsToMesh(levelSet, mesh).apply();
//   ls::VTKWriter(mesh, fileName).apply();
//   lsToSurfaceMesh(levelSet, mesh).apply();
//   ls::VTKWriter(mesh, "surf_" + fileName).apply();
// }

int main() {
  omp_set_num_threads(1);

  LSType mask;
  {
    double gridDelta = 1.0;
    double extent = 10;
    double bounds[2 * D] = {-extent, extent, -extent, extent};

    typename ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
    boundaryCons[0] =
        ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
    boundaryCons[1] =
        ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

    if constexpr (D == 3) {
      boundaryCons[1] =
          ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
      boundaryCons[2] =
          ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;
    }

    mask = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);
  }

  constexpr double yVal = 0.5;

  // make 1st box
  {
    double min[D] = {-7., yVal};
    double max[D] = {0., 10.};
    // if constexpr (D == 3) {
    //   min[1] = -5.;
    //   min[2] = 1.;
    //   max[1] = 5.;
    //   max[2] = 10.;
    // }
    ls::MakeGeometry(mask,
                     ls::SmartPointer<ls::Box<NumericType, D>>::New(min, max))
        .apply();
    // writeLS(mask, "mask_initial.vtp");
  }

  // make 2nd box
  auto substrate = LSType::New(mask->getGrid());
  {
    double min[D] = {0., yVal};
    double max[D] = {7., 10.};
    // if constexpr (D == 3) {
    //   min[1] = -5.;
    //   min[2] = 1.;
    //   max[1] = 5.;
    //   max[2] = 10.;
    // }
    ls::MakeGeometry(substrate,
                     ls::SmartPointer<ls::Box<NumericType, D>>::New(min, max))
        .apply();
    // writeLS(substrate, "subs_initial.vtp");
  }

  ls::BooleanOperation<NumericType, D> boolOp(substrate, mask,
                                              ls::BooleanOperationEnum::UNION);
  boolOp.setPruneResult(false);
  boolOp.apply();

  // writeLS(substrate, "subs_afterBool.vtp");

  ls::Prune<NumericType, D> pruner(substrate);
  pruner.setRemoveStrayZeros(true);
  pruner.apply();

  // writeLS(substrate, "subs_afterPrune1.vtp");

  ls::Prune(substrate).apply();

  // writeLS(substrate, "subs_afterPrune2.vtp");

  LSTEST_ASSERT_VALID_LS(substrate, NumericType, D);

  return 0;
}
