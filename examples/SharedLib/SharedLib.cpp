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

namespace ls = viennals;

int main() {

  // do not need to define dimension
  // since we are using predefined typedefs
  // constexpr int D = 3;
  omp_set_num_threads(4);

  double gridDelta = 0.25;

  // Usually we would use ls::Domain<float, D>.
  // Since we want to make sure we get an error
  // if we do not use a pre-built type, we use
  // the specialization typedef
  auto sphere1 = ls::SmartPointer<ls::Domain_float_3>::New(gridDelta);
  auto sphere2 = ls::SmartPointer<ls::Domain_float_3>::New(gridDelta);

  float origin[3] = {5., 0., 0.};
  float radius = 7.3;

  {
    // these typedefs are available for all templated classes
    auto sphere = ls::SmartPointer<ls::Sphere_float_3>::New(origin, radius);
    ls::MakeGeometry_float_3(sphere1, sphere).apply();
    origin[0] = -5.0;
    radius = 9.5;
    sphere = ls::SmartPointer<ls::Sphere_float_3>::New(origin, radius);
    ls::MakeGeometry_float_3(sphere2, sphere).apply();
  }

  {
    auto mesh1 = ls::SmartPointer<ls::Mesh<float>>::New();
    auto mesh2 = ls::SmartPointer<ls::Mesh<float>>::New();

    std::cout << "Extracting..." << std::endl;
    ls::ToSurfaceMesh_float_3(sphere1, mesh1).apply();
    ls::ToSurfaceMesh_float_3(sphere2, mesh2).apply();

    ls::VTKWriter<float>(mesh1, "sphere1.vtp").apply();
    ls::VTKWriter<float>(mesh2, "sphere2.vtp").apply();
  }

  // Perform a boolean operation
  ls::BooleanOperation_float_3(sphere1, sphere2,
                               ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
      .apply();

  std::cout << "Extracting..." << std::endl;
  auto mesh = ls::SmartPointer<ls::Mesh<float>>::New();
  ls::ToSurfaceMesh_float_3(sphere1, mesh).apply();

  mesh->print();

  ls::VTKWriter<float>(mesh, "after.vtp").apply();

  return 0;
}
