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

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  double extent = 15;
  double gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  auto sphere1 =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

  double origin[D] = {5., 0.};
  double radius = 7.3;

  lsMakeGeometry<double, D>(
      sphere1, lsSmartPointer<lsSphere<double, D>>::New(origin, radius))
      .apply();

  {
    auto mesh = lsSmartPointer<lsMesh>::New();
    lsToMesh<double, D>(sphere1, mesh).apply();
    lsVTKWriter(mesh, "sphere.vtk").apply();
  }

  {
    auto mesh = lsSmartPointer<lsMesh>::New();
    lsExpand<double, D>(sphere1, 5).apply();
    lsToMesh<double, D>(sphere1, mesh).apply();
    lsVTKWriter(mesh, "sphereExpanded.vtk").apply();

    lsReduce<double, D>(sphere1, 1).apply();
    lsToMesh<double, D>(sphere1, mesh).apply();
    lsVTKWriter(mesh, "sphereReduced.vtk").apply();
  }

  return 0;
}
