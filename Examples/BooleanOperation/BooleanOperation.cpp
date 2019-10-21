#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToExplicitMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Example of boolean operations on level sets
  using two spheres.
  \example BooleanOperation.cpp
*/

int main() {

  constexpr int D = 3;
  omp_set_num_threads(4);

  double gridDelta = 0.25;

  lsDomain<double, D> sphere1(gridDelta);
  lsDomain<double, D> sphere2(gridDelta);

  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  lsMakeGeometry<double, D>(sphere1).makeSphere(origin, radius);
  origin[0] = -5.0;
  radius = 9.5;
  lsMakeGeometry<double, D>(sphere2).makeSphere(origin, radius);

  {
    lsMesh mesh1, mesh2;

    std::cout << "Extracting..." << std::endl;
    lsToExplicitMesh<double, D>(sphere1, mesh1).apply();
    lsToExplicitMesh<double, D>(sphere2, mesh2).apply();

    lsVTKWriter(mesh1).writeVTKLegacy("sphere1.vtk");
    lsVTKWriter(mesh2).writeVTKLegacy("sphere2.vtk");
  }

  // Perform a boolean operation
  lsBooleanOperation<double, D>(sphere1).XOR(sphere2);

  std::cout << "Extracting..." << std::endl;
  lsMesh mesh;
  lsToExplicitMesh<double, D>(sphere1, mesh).apply();

  mesh.print();

  lsVTKWriter(mesh).writeVTKLegacy("after.vtk");

  return 0;
}
