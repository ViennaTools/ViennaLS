#include <algorithm>
#include <iostream>

#include <lsDomain.hpp>
#include <lsFromVolumeMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

int main() {
  constexpr int D = 2;
  omp_set_num_threads(1);

  double gridDelta = 5e-10;
  double bounds[2 * D] = {-3.5e-8, 3.5e-8, -5e-8, 5e-8};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  boundaryCons[1] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  std::vector<lsSmartPointer<lsDomain<double, D>>> levelSets(
      5, lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons,
                                                  gridDelta));

  // Read mesh
  auto initialMesh = lsSmartPointer<lsMesh<>>::New();
  lsVTKReader<double>(initialMesh, "initial.vtk").apply();
  initialMesh->print();

  lsFromVolumeMesh<double, D>(levelSets, initialMesh).apply();

  for (unsigned i = 0; i < levelSets.size(); ++i) {
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<double, D>(levelSets[i], mesh).apply();
    lsVTKWriter<double>(mesh, "LSsurface-" + std::to_string(i) + ".vtk")
        .apply();
  }

  return 0;
}
