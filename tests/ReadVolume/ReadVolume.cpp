#include <algorithm>
#include <iostream>

#include <lsDomain.hpp>
#include <lsFromVolumeMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

int main() {
  constexpr int D = 2;
  omp_set_num_threads(1);

  double gridDelta = 5e-10;
  double bounds[2 * D] = {-3.5e-8, 3.5e-8, -5e-8, 5e-8};
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  boundaryCons[1] = ls::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto domain = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  // Read mesh
  auto initialMesh = ls::SmartPointer<ls::Mesh<>>::New();
  ls::VTKReader<double>(initialMesh, "initial.vtu").apply();
  initialMesh->print();

  ls::FromVolumeMesh<double, D> reader(domain->getGrid(), initialMesh);
  reader.apply();
  auto levelSets = reader.getLevelSets();

  for (unsigned i = 0; i < levelSets.size(); ++i) {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(levelSets[i], mesh).apply();
    ls::VTKWriter<double>(mesh, "LSsurface-" + std::to_string(i) + ".vtk")
        .apply();
  }

  return 0;
}
