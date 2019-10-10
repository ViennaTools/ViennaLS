#include <iostream>

#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsReduce.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  double extent = 15;
  double gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  lsDomain_double_2::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = lsDomain_double_2::BoundaryType::SYMMETRIC_BOUNDARY;
  lsDomain_double_2 sphere1(bounds, boundaryCons, gridDelta);
  // lsDomain_double_3 sphere2(bounds, gridDelta, boundaryCons);

  double origin[D] = {5., 0.};
  double radius = 7.3;

  lsMakeGeometry<double, D>(sphere1).makeSphere(origin, radius, 2);

  {
    lsMesh mesh;
    lsToMesh<double, D>(sphere1, mesh).apply();
    lsVTKWriter(mesh).writeVTKLegacy("sphere.vtk");
  }

  {
    lsMesh mesh;
    lsExpand<double, D>(sphere1).apply(5);
    lsToMesh<double, D>(sphere1, mesh).apply();
    lsVTKWriter(mesh).writeVTKLegacy("sphereExpanded.vtk");

    lsReduce<double, D>(sphere1).apply(1);
    lsToMesh<double, D>(sphere1, mesh).apply();
    lsVTKWriter(mesh).writeVTKLegacy("sphereReduced.vtk");
  }

  return 0;
}
