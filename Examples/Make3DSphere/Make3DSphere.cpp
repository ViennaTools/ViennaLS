#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromExplicitMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToExplicitMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsToVoxelMesh.hpp>

int main() {

  constexpr int D = 3;

  omp_set_num_threads(1);

  double gridDelta = 0.35;

  lsDomain_double_3 sphere1(gridDelta); //, boundaryCons);

  lsDomain_double_3 sphere2(gridDelta); //, boundaryCons);
  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  lsMakeGeometry<double, D>(sphere1).makeSphere(origin, radius);
  origin[0] = -5.;
  lsMakeGeometry<double, D>(sphere2).makeSphere(origin, radius);

  std::cout << "Number of points: " << sphere1.getDomain().getNumberOfPoints()
            << std::endl;
  lsMesh mesh;

  std::cout << "Expanding..." << std::endl;
  lsExpand<double, D>(sphere1).apply(2);
  lsExpand<double, D>(sphere2).apply(2);

  std::cout << "Booling..." << std::endl;
  lsBooleanOperation<double, D>(sphere1).XOR(sphere2);

  std::cout << "Extracting..." << std::endl;
  lsToExplicitMesh<double, D>(sphere1, mesh).apply();

  mesh.print();

  lsVTKWriter(mesh).writeVTKLegacy("test-" + std::to_string(radius) + ".vtk");

  // write voxelised volume mesh
  {
    lsMesh voxelMesh;
    lsToVoxelMesh<double, D>(sphere1, voxelMesh).apply();
    std::cout << "voxelMesh: " << std::endl;
    voxelMesh.print();

    lsVTKWriter(voxelMesh).writeVTU("voxelMesh.vtu");
  }


  std::cout << "Reading mesh again: " << std::endl;

  lsDomain<double, D> levelSet(gridDelta);

  lsFromExplicitMesh<double, D>(levelSet).apply(mesh);

  std::cout << "Create second mesh: " << std::endl;
  lsMesh mesh2;
  lsToExplicitMesh<double, D>(levelSet, mesh2).apply();

  mesh.print();

  lsVTKWriter(mesh).writeVTKLegacy("test2-" + std::to_string(radius) + ".vtk");

  return 0;
}
