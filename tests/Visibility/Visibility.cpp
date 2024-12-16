#include <chrono>
#include <iostream>

#include <lsCalculateVisibilities.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

int main() {

  constexpr int D = 2;

  double gridDelta = 0.4;

  auto sphere1 = ls::SmartPointer<ls::Domain<double, D>>::New(
      gridDelta); //, boundaryCons);

  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  ls::MakeGeometry<double, D>(
      sphere1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
      .apply();

  ls::Expand<double, D>(sphere1, 5).apply();
  ls::CalculateVisibilities<double, D>(sphere1, ls::Vec3D<double>{0., -1.0, 0.})
      .apply();

  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, D>(sphere1, mesh).apply();
  ls::VTKWriter<double>(mesh, "visibility_test.vtp").apply();

  return 0;
}
