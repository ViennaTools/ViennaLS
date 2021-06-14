#include <chrono>
#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToSurfaceMesh.hpp>
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
  auto mesh = lsSmartPointer<lsMesh<>>::New();

  std::cout << "Expanding..." << std::endl;
  lsExpand<double, D>(sphere1, 2).apply();
  lsExpand<double, D>(sphere2, 2).apply();

  std::cout << "Booling..." << std::endl;
  lsBooleanOperation<double, D>(sphere1, sphere2, lsBooleanOperationEnum::UNION)
      .apply();

  //   std::cout << "Extracting..." << std::endl;
  lsToDiskMesh<double, D>(sphere1, mesh).apply();
  //   std::cout << "Disk mesh:" << std::endl;
  //   mesh->print();
  //   lsVTKWriter<double>(mesh, "disks-" + std::to_string(radius) +
  //   ".vtk").apply();

  //   lsToSurfaceMesh<double, D>(sphere1, mesh).apply();
  //   std::cout << "Surface mesh:" << std::endl;
  //   mesh->print();
  //   lsVTKWriter<double>(mesh, "surface-" + std::to_string(radius) + ".vtk")
  //       .apply();

  return 0;
}
