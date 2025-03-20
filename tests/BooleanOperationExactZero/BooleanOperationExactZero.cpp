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

#include "result.hpp"

/**
  Check if Boolean Operation works with exact zero LS values
  \example BooleanOperationExactZero.cpp
*/

namespace ls = viennals;

using NumericType = double;
constexpr int D = 2;

using LSType = ls::SmartPointer<ls::Domain<NumericType, D>>;

void writeLS(LSType levelSet, std::string fileName) {
  auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
  ls::ToMesh(levelSet, mesh, false).apply();
  ls::VTKWriter(mesh, fileName).apply();
}

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
    ls::MakeGeometry(mask,
                     ls::SmartPointer<ls::Box<NumericType, D>>::New(min, max))
        .apply();
    // writeLS(mask, "mask_initial.vtp");
  }

  // now make substrate ( plane ) at the same height as the bottom of the mask
  auto substrate = LSType::New(mask->getGrid());
  {
    double origin[D] = {};
    double normal[D] = {};
    origin[D - 1] = 1.;
    normal[D - 1] = 1.;
    ls::MakeGeometry(
        substrate,
        ls::SmartPointer<ls::Plane<NumericType, D>>::New(origin, normal))
        .apply();
    // writeLS(substrate, "subs_initial.vtp");
  }

  ls::BooleanOperation(substrate, mask, ls::BooleanOperationEnum::UNION)
      .apply();

  // writeLS(substrate, "subs_afterBool.vtp");

  // iterate through all values and check if they are correct
  unsigned counter = 0;
  for (viennahrle::ConstSparseIterator<
           typename ls::Domain<double, D>::DomainType>
           it(substrate->getDomain());
       !it.isFinished(); ++it) {
    auto indices = it.getStartIndices();
    LSTEST_ASSERT(indices == resultIndices[counter]);

    auto val = it.getValue();
    if (it.isDefined()) {
      LSTEST_ASSERT((val - resultValues[counter]) < 1e-4)
    } else {
      LSTEST_ASSERT((val < 0.) == (resultValues[counter] < 0.))
    }
    ++counter;
  }

  return 0;
}
