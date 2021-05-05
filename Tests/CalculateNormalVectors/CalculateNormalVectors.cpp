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

int main() {

  constexpr int D = 3;

  omp_set_num_threads(1);

  double extent = 15;
  double gridDelta = 0.25;

  double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[3];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  auto sphere1 =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  lsMakeGeometry<double, D>(
      sphere1, lsSmartPointer<lsSphere<double, D>>::New(origin, radius))
      .apply();

  {
    auto sphere2 = lsSmartPointer<lsDomain<double, D>>::New(
        bounds, boundaryCons, gridDelta);
    origin[0] = -5.;
    lsMakeGeometry<double, D>(
        sphere2, lsSmartPointer<lsSphere<double, D>>::New(origin, radius))
        .apply();
  }

  // std::cout << "Expanding..." << std::endl;
  lsExpand<double, D>(sphere1, 3).apply();

  // std::cout << "Number of points: " << sphere1->getNumberOfPoints()
  // << std::endl;

  // normal vectors are only valid as long as the underlying
  // level set does not change
  lsCalculateNormalVectors<double, 3>(sphere1).apply();

  // auto mesh = lsSmartPointer<lsMesh<>>::New();
  // lsToMesh<double, 3>(sphere1, mesh, true, true).apply();

  // auto writer = lsVTKWriter<double>();
  // writer.setMesh(mesh);
  // writer.setFileName("explicit.vtk");
  // writer.apply();

  LSTEST_ASSERT_VALID_LS(sphere1, double, D)

  return 0;
}
