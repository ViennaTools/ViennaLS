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

  auto sphere1 =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

  double origin[D] = {5., 0.};
  double radius = 7.3;

  lsMakeGeometry<double, D>(
      sphere1, lsSmartPointer<lsSphere<double, D>>::New(origin, radius))
      .apply();

  std::cout << "Writing" << std::endl;
  {
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToMesh<double, D>(sphere1, mesh).apply();
    lsVTKWriter<double>(mesh, "sphere.vtk").apply();
  }

  std::cout << "Reading" << std::endl;
  {
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsVTKReader<double>(mesh, "sphere.vtk").apply();
    auto newLS = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons,
                                                          gridDelta);
    lsFromMesh<double, D>(newLS, mesh).apply();

    std::cout << "Writing new" << std::endl;
    auto newMesh = lsSmartPointer<lsMesh<>>::New();
    lsToMesh<double, D>(newLS, newMesh).apply();
    lsVTKWriter<double>(mesh, "newMesh.vtk").apply();
  }

  return 0;
}
