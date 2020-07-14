#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  This example shows how to make sure the pre-built shared library is used.
  ViennaLS is a header only project. However, certain commonly used
  template specializations can be built with the library and used during
  developement to decrease compile times. In order to make sure these
  pre-compiled types are used, typedefs for those types are included with each
  header. These typedefs will also be available when building as header-only in
  order to avoid code changes. The available specialisations are listed in
  "lsPreCompileMacros.hpp" \example SharedLib.cpp
*/

int main() {

  // do not need to define dimension
  // since we are using predefined typedefs
  // constexpr int D = 3;
  omp_set_num_threads(4);

  double gridDelta = 0.25;

  // Usually we would use lsDomain<float, D>.
  // Since we want to make sure we get an error
  // if we do not use a pre-built type, we use
  // the specialization typedef
  auto sphere1 = lsSmartPointer<lsDomain_float_3>::New(gridDelta);
  auto sphere2 = lsSmartPointer<lsDomain_float_3>::New(gridDelta);

  float origin[3] = {5., 0., 0.};
  float radius = 7.3;

  {
    // these typedefs are available for all templated classes
    auto sphere = lsSmartPointer<lsSphere_float_3>::New(origin, radius);
    lsMakeGeometry_float_3(sphere1, sphere).apply();
    origin[0] = -5.0;
    radius = 9.5;
    sphere = lsSmartPointer<lsSphere_float_3>::New(origin, radius);
    lsMakeGeometry_float_3(sphere2, sphere).apply();
  }

  {
    auto mesh1 = lsSmartPointer<lsMesh>::New();
    auto mesh2 = lsSmartPointer<lsMesh>::New();

    std::cout << "Extracting..." << std::endl;
    lsToSurfaceMesh_float_3(sphere1, mesh1).apply();
    lsToSurfaceMesh_float_3(sphere2, mesh2).apply();

    lsVTKWriter(mesh1, "sphere1.vtk").apply();
    lsVTKWriter(mesh2, "sphere2.vtk").apply();
  }

  // Perform a boolean operation
  lsBooleanOperation_float_3(sphere1, sphere2,
                             lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
      .apply();

  std::cout << "Extracting..." << std::endl;
  auto mesh = lsSmartPointer<lsMesh>::New();
  lsToSurfaceMesh_float_3(sphere1, mesh).apply();

  mesh->print();

  lsVTKWriter(mesh, "after.vtk").apply();

  return 0;
}
