#include <iostream>

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
  Example of boolean operations on level sets
  using two spheres.
  \example BooleanOperation.cpp
*/

int main() {

  constexpr int D = 3;
  omp_set_num_threads(4);

  double gridDelta = 1.0;
  double bounds[2 * D] = {-20, 20, -20, 20};
  if (D == 3) {
    bounds[4] = -20;
    bounds[5] = 20;
  }
  typename lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto sphere1 =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);
  auto sphere2 =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

  double origin[3] = {0., 0., 0.};
  double radius = 15.7;

  lsMakeGeometry<double, D>(
      sphere1, lsSmartPointer<lsSphere<double, D>>::New(origin, radius))
      .apply();
  origin[0] = 15.0;
  radius = 9.5;
  lsMakeGeometry<double, D>(
      sphere2, lsSmartPointer<lsSphere<double, D>>::New(origin, radius))
      .apply();

  // {
  //   auto mesh1 = lsSmartPointer<lsMesh<>>::New();
  //   auto mesh2 = lsSmartPointer<lsMesh<>>::New();

  //   std::cout << "Extracting..." << std::endl;
  //   lsToSurfaceMesh<double, D>(sphere1, mesh1).apply();
  //   lsToSurfaceMesh<double, D>(sphere2, mesh2).apply();

  //   lsVTKWriter<double>(mesh1, "sphere1.vtk").apply();
  //   lsVTKWriter<double>(mesh2, "sphere2.vtk").apply();

  //   lsToMesh<double, D>(sphere1, mesh1).apply();
  //   lsToMesh<double, D>(sphere2, mesh2).apply();

  //   lsVTKWriter<double>(mesh1, "LS1.vtk").apply();
  //   lsVTKWriter<double>(mesh2, "LS2.vtk").apply();
  // }

  // Perform a boolean operation
  lsBooleanOperation<double, D>(sphere1, sphere2, lsBooleanOperationEnum::UNION)
      .apply();

  LSTEST_ASSERT_VALID_LS(sphere1, double, D)

  // std::cout << "Extracting..." << std::endl;
  // auto mesh = lsSmartPointer<lsMesh<>>::New();
  // lsToSurfaceMesh<double, D>(sphere1, mesh).apply();

  // mesh->print();

  // lsVTKWriter<double>(mesh, "after.vtk").apply();

  return 0;
}
