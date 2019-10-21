#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToExplicitMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  This example is based on the boolean operation example.
  It shows how to make sure the pre-built shared library is used.
  ViennaLS is a header only project. However, certain commonly used
  template specializations can be built with the library and used during
  developement to decrese compile times. In order to make sure these
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

  // usually we would use lsDomain<float, D>
  // since we want to make sure we get an error
  // if we do not use a pre-built type, we use
  // the specialization typedef
  lsDomain_float_3 sphere1(gridDelta);
  lsDomain_float_3 sphere2(gridDelta);

  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  // these typedefs are available for all templated classes
  lsMakeGeometry_float_3(sphere1).makeSphere(origin, radius);
  origin[0] = -5.0;
  radius = 9.5;
  lsMakeGeometry_float_3(sphere2).makeSphere(origin, radius);

  {
    lsMesh mesh1, mesh2;

    std::cout << "Extracting..." << std::endl;
    lsToExplicitMesh_float_3(sphere1, mesh1).apply();
    lsToExplicitMesh_float_3(sphere2, mesh2).apply();

    lsVTKWriter(mesh1).writeVTKLegacy("sphere1.vtk");
    lsVTKWriter(mesh2).writeVTKLegacy("sphere2.vtk");
  }

  // Perform a boolean operation
  lsBooleanOperation_float_3(sphere1).XOR(sphere2);

  std::cout << "Extracting..." << std::endl;
  lsMesh mesh;
  lsToExplicitMesh_float_3(sphere1, mesh).apply();

  mesh.print();

  lsVTKWriter(mesh).writeVTKLegacy("after.vtk");

  return 0;
}
