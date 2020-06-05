#include <iostream>

#include <lsFromMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  double extent = 15;
  double gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  lsDomain<double, D> sphere1(bounds, boundaryCons, gridDelta);
  // lsDomain<double, D> sphere2(bounds, gridDelta, boundaryCons);

  double origin[D] = {5., 0.};
  double radius = 7.3;

  lsMakeGeometry<double, D>(sphere1, lsSphere<double, D>(origin, radius))
      .apply();

  std::cout << "Writing" << std::endl;
  {
    lsMesh mesh;
    lsToMesh<double, D>(sphere1, mesh).apply();
    lsVTKWriter(mesh, "sphere.vtk").apply();
  }

  std::cout << "Reading" << std::endl;
  {
    lsMesh mesh;
    lsVTKReader(mesh, "sphere.vtk").apply();
    lsDomain<double, D> newLS(bounds, boundaryCons, gridDelta);
    lsFromMesh<double, D>(newLS, mesh).apply();

    std::cout << "Writing new" << std::endl;
    lsMesh newMesh;
    lsToMesh<double, D>(newLS, newMesh).apply();
    lsVTKWriter(mesh, "newMesh.vtk").apply();
  }

  return 0;
}
