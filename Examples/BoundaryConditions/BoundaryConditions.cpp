#include <iostream>

#include <lsAdvect.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToExplicitMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example showing how to set boundary conditions
  for an lsDomain.
  \example BoundaryConditions.cpp
*/

int main() {

  constexpr int D = 3;
  omp_set_num_threads(4);

  double extent = 15;

  double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[3];
  for (unsigned i = 0; i < D - 1; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::SYMMETRIC_BOUNDARY;

  boundaryCons[2] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  lsDomain<double, D> levelSet(bounds, boundaryCons, 0.1);

  hrleVectorType<double, D> origin(0., 0., 0.);
  hrleVectorType<double, D> normalVector(0., 1., 1.);

  lsMakeGeometry<double, D>(levelSet).makePlane(origin, normalVector);

  {
    std::cout << "Extracting..." << std::endl;
    lsMesh mesh;
    lsToExplicitMesh<double, D>(levelSet, mesh).apply();
    lsVTKWriter(mesh).writeVTKLegacy("plane.vtk");
  }

  return 0;
}
