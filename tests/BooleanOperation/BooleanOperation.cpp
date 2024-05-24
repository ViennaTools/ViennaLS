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

namespace ls = viennals;

/**
  Example of boolean operations on level sets
  using two spheres.
  \example BooleanOperation.cpp
*/

int main() {

  constexpr int D = 3;
  omp_set_num_threads(4);

  double gridDelta = 1.0;
  double bounds[2 * D] = {-20, 20, -20, 20};
  if constexpr (D == 3) {
    bounds[4] = -20;
    bounds[5] = 20;
  }
  typename ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] = ls::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto sphere1 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  auto sphere2 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin[3] = {0., 0., 0.};
  double radius = 15.7;

  ls::MakeGeometry<double, D>(
      sphere1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
      .apply();
  origin[0] = 15.0;
  radius = 9.5;
  ls::MakeGeometry<double, D>(
      sphere2, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
      .apply();

  // put some data into either LS
  {
    using ScalarDataType = ls::Domain<double, D>::PointDataType::ScalarDataType;
    ScalarDataType scalars1(sphere1->getNumberOfPoints(), 0.);
    sphere1->getPointData().insertNextScalarData(std::move(scalars1),
                                                 "originID");

    ScalarDataType scalars2(sphere2->getNumberOfPoints(), 1.);
    sphere2->getPointData().insertNextScalarData(std::move(scalars2),
                                                 "originID");
  }

  // {
  //   auto mesh1 = ls::ls::SmartPointer<ls::Mesh<>>::New();
  //   auto mesh2 = ls::ls::SmartPointer<ls::Mesh<>>::New();

  //   std::cout << "Extracting..." << std::endl;
  //   ls::ToSurfaceMesh<double, D>(sphere1, mesh1).apply();
  //   ls::ToSurfaceMesh<double, D>(sphere2, mesh2).apply();

  //   ls::VTKWriter<double>(mesh1, "sphere1.vtk").apply();
  //   ls::VTKWriter<double>(mesh2, "sphere2.vtk").apply();

  //   ls::ToMesh<double, D>(sphere1, mesh1).apply();
  //   ls::ToMesh<double, D>(sphere2, mesh2).apply();

  //   ls::VTKWriter<double>(mesh1, "LS1.vtk").apply();
  //   ls::VTKWriter<double>(mesh2, "LS2.vtk").apply();
  // }

  // Perform a boolean operation
  ls::BooleanOperation<double, D>(sphere1, sphere2,
                                  ls::BooleanOperationEnum::UNION)
      .apply();

  LSTEST_ASSERT_VALID_LS(sphere1, double, D)

  // std::cout << "Extracting..." << std::endl;
  // auto mesh = ls::ls::SmartPointer<ls::Mesh<>>::New();
  // ls::ToSurfaceMesh<double, D>(sphere1, mesh).apply();
  // ls::VTKWriter<double>(mesh, "after.vtp").apply();

  return 0;
}
