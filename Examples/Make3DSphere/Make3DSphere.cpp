#include <iostream>
#include <chrono>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromExplicitMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToExplicitMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example showing how to write and read different
  meshes created by the algorithms lsToVoxelMesh and lsToExplicitMesh.
  \example Make3DSphere.cpp
*/

int main() {

  constexpr int D = 3;

  omp_set_num_threads(1);

  double gridDelta = 0.4;

  lsDomain<double, D> sphere1(gridDelta); //, boundaryCons);

  lsDomain<double, D> sphere2(gridDelta); //, boundaryCons);
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
  lsBooleanOperation<double, D>(sphere1).OR(sphere2);

  std::cout << "Extracting..." << std::endl;
  lsToExplicitMesh<double, D>(sphere1, mesh).apply();

  mesh.print();

  lsVTKWriter(mesh).writeVTKLegacy("test-" + std::to_string(radius) + ".vtk");

  // write voxelised volume mesh
  {
    lsMesh voxelMesh;
    auto voxelMesher = lsToVoxelMesh<double, D>(voxelMesh);

    voxelMesher.insertNextLevelSet(sphere2);
    voxelMesher.insertNextLevelSet(sphere1);

    std::cout << "voxelising" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    voxelMesher.apply();
    auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Converting to voxel mesh took "
              << std::chrono::duration_cast<std::chrono::seconds>(stop - start)
                     .count()
              << "s" << std::endl;

    std::cout << "voxelMesh: " << std::endl;
    voxelMesh.print();

    lsVTKWriter(voxelMesh).writeVTU("voxelMesh.vtu");
  }

  return 0;
}
