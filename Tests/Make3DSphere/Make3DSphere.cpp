#include <chrono>
#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example showing how to write different
  meshes created by the algorithms lsToVoxelMesh and lsToSurfaceMesh.
  \example Make3DSphere.cpp
*/

int main() {

  constexpr int D = 3;

  omp_set_num_threads(1);

  double gridDelta = 0.4;

  auto sphere1 =
      lsSmartPointer<lsDomain<double, D>>::New(gridDelta); //, boundaryCons);

  auto sphere2 =
      lsSmartPointer<lsDomain<double, D>>::New(gridDelta); //, boundaryCons);
  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  lsMakeGeometry<double, D>(
      sphere1, lsSmartPointer<lsSphere<double, D>>::New(origin, radius))
      .apply();
  origin[0] = -5.;
  lsMakeGeometry<double, D>(
      sphere2, lsSmartPointer<lsSphere<double, D>>::New(origin, radius))
      .apply();

  std::cout << "Number of points: " << sphere1->getDomain().getNumberOfPoints()
            << std::endl;
  auto mesh = lsSmartPointer<lsMesh>::New();

  std::cout << "Expanding..." << std::endl;
  lsExpand<double, D>(sphere1, 2).apply();
  lsExpand<double, D>(sphere2, 2).apply();

  std::cout << "Booling..." << std::endl;
  lsBooleanOperation<double, D>(sphere1, sphere2, lsBooleanOperationEnum::UNION)
      .apply();

  std::cout << "Extracting..." << std::endl;
  lsToSurfaceMesh<double, D>(sphere1, mesh).apply();

  mesh->print();

  lsVTKWriter(mesh, "test-" + std::to_string(radius) + ".vtk").apply();

  // write voxelised volume mesh
  {
    auto voxelMesh = lsSmartPointer<lsMesh>::New();
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
    voxelMesh->print();

    lsVTKWriter(voxelMesh, lsFileFormatEnum::VTU, "voxelMesh.vtu").apply();
  }

  return 0;
}
