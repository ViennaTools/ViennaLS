#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsTestAsserts.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Example showing how to calculate the normal vectors for
  every grid point stored in an lsDomain and outputting
  them to an explicit mesh. This also outputs the level set values
  stored in each grid point.
  \example calculateNormalVectors.cpp
*/

namespace ls = viennals;

int main() {

  constexpr int D = 3;

  omp_set_num_threads(1);

  double extent = 15;
  double gridDelta = 0.25;

  double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
  ls::Domain<double, D>::BoundaryType boundaryCons[3];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  auto sphere1 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  ls::MakeGeometry<double, D>(
      sphere1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
      .apply();

  {
    auto sphere2 = ls::SmartPointer<ls::Domain<double, D>>::New(
        bounds, boundaryCons, gridDelta);
    origin[0] = -5.;
    ls::MakeGeometry<double, D>(
        sphere2, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
        .apply();
  }

  // std::cout << "Expanding..." << std::endl;
  ls::Expand<double, D>(sphere1, 3).apply();

  // std::cout << "Number of points: " << sphere1->getNumberOfPoints()
  // << std::endl;

  // normal vectors are only valid as long as the underlying
  // level set does not change
  ls::CalculateNormalVectors<double, 3>(sphere1).apply();

  // auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  // ls::ToMesh<double, 3>(sphere1, mesh, true, true).apply();

  // auto writer = ls::VTKWriter<double>();
  // writer.setMesh(mesh);
  // writer.setFileName("explicit.vtk");
  // writer.apply();

  LSTEST_ASSERT_VALID_LS(sphere1, double, D)

  return 0;
}
