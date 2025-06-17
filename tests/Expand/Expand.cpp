#include <iostream>

#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsReduce.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Simple 2D example showing how to use lsExpand in order
  to increase the number of stored grid points around
  a level set interface.
  \example Expand.cpp
*/

namespace ls = viennals;

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  double extent = 15;
  double gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  typename ls::BoundaryConditionEnum boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;

  auto sphere1 = ls::Domain<double, D>::New(bounds, boundaryCons, gridDelta);

  double origin[D] = {5., 0.};
  double radius = 7.3;

  ls::MakeGeometry<double, D>(sphere1,
                              ls::Sphere<double, D>::New(origin, radius))
      .apply();

  {
    auto mesh = ls::Mesh<>::New();
    ls::ToMesh<double, D>(sphere1, mesh).apply();
    ls::VTKWriter<double>(mesh, "sphere.vtk").apply();
  }

  {
    auto mesh = ls::Mesh<>::New();
    ls::Expand<double, D>(sphere1, 5).apply();
    ls::ToMesh<double, D>(sphere1, mesh).apply();
    ls::VTKWriter<double>(mesh, "sphereExpanded.vtk").apply();

    ls::Reduce<double, D>(sphere1, 1).apply();
    ls::ToMesh<double, D>(sphere1, mesh).apply();
    ls::VTKWriter<double>(mesh, "sphereReduced.vtk").apply();
  }

  return 0;
}
