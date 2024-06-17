#include <iostream>

#include <lsAdvect.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsTestAsserts.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example showing how to set boundary conditions
  for an lsDomain.
  \example BoundaryConditions.cpp
*/

namespace ls = viennals;

int main() {

  constexpr int D = 3;
  omp_set_num_threads(4);

  double gridDelta = 0.1;
  double extent = 15;

  double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
  ls::Domain<double, D>::BoundaryType boundaryCons[3];
  for (unsigned i = 0; i < D - 1; ++i)
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  boundaryCons[2] = ls::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto levelSet = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  hrleVectorType<double, D> origin(0., 0., 0.);
  hrleVectorType<double, D> normalVector(0., 1., 1.);

  ls::MakeGeometry<double, D>(
      levelSet,
      ls::SmartPointer<ls::Plane<double, D>>::New(origin, normalVector))
      .apply();

  // {
  //   std::cout << "Extracting..." << std::endl;
  //   auto mesh = ls::ls::SmartPointer<ls::Mesh<>>::New();
  //   ls::ToSurfaceMesh<double, D>(levelSet, mesh).apply();
  //   ls::VTKWriter<double>(mesh, "plane.vtk").apply();
  // }

  LSTEST_ASSERT_VALID_LS(levelSet, double, D)

  return 0;
}
