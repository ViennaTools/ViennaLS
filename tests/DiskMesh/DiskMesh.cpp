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

namespace ls = viennals;

int main() {

  constexpr int D = 3;

  omp_set_num_threads(1);

  double gridDelta = 0.4;

  auto sphere1 = ls::SmartPointer<ls::Domain<double, D>>::New(
      gridDelta); //, boundaryCons);

  auto sphere2 = ls::SmartPointer<ls::Domain<double, D>>::New(
      gridDelta); //, boundaryCons);
  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  ls::MakeGeometry<double, D>(
      sphere1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
      .apply();
  origin[0] = -5.;
  ls::MakeGeometry<double, D>(
      sphere2, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
      .apply();

  std::cout << "Number of points: " << sphere1->getDomain().getNumberOfPoints()
            << std::endl;
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();

  std::cout << "Expanding..." << std::endl;
  ls::Expand<double, D>(sphere1, 2).apply();
  ls::Expand<double, D>(sphere2, 2).apply();

  std::cout << "Booling..." << std::endl;
  ls::BooleanOperation<double, D>(sphere1, sphere2,
                                  ls::BooleanOperationEnum::UNION)
      .apply();

  //   std::cout << "Extracting..." << std::endl;
  ls::ToDiskMesh<double, D>(sphere1, mesh).apply();
  //   std::cout << "Disk mesh:" << std::endl;
  //   mesh->print();
  //   ls::VTKWriter<double>(mesh, "disks-" + std::to_string(radius) +
  //   ".vtk").apply();

  //   ls::ToSurfaceMesh<double, D>(sphere1, mesh).apply();
  //   std::cout << "Surface mesh:" << std::endl;
  //   mesh->print();
  //   ls::VTKWriter<double>(mesh, "surface-" + std::to_string(radius) + ".vtk")
  //       .apply();

  return 0;
}
