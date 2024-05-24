#include <iostream>

#include <lsFromMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  double extent = 15;
  double gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  auto sphere1 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin[D] = {5., 0.};
  double radius = 7.3;

  ls::MakeGeometry<double, D>(
      sphere1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
      .apply();

  std::cout << "Writing" << std::endl;
  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(sphere1, mesh).apply();
    ls::VTKWriter<double>(mesh, "sphere.vtk").apply();
  }

  std::cout << "Reading" << std::endl;
  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::VTKReader<double>(mesh, "sphere.vtk").apply();
    auto newLS = ls::SmartPointer<ls::Domain<double, D>>::New(
        bounds, boundaryCons, gridDelta);
    ls::FromMesh<double, D>(newLS, mesh).apply();

    std::cout << "Writing new" << std::endl;
    auto newMesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(newLS, newMesh).apply();
    ls::VTKWriter<double>(mesh, "newMesh.vtk").apply();
  }

  return 0;
}
