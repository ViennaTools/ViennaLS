#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

#include <lsWriteVisualizationMesh.hpp>

namespace ls = viennals;

int main() {
  omp_set_num_threads(4);

  constexpr int D = 2;
  typedef float NumericType;

  double gridDelta = 1.0;

  double bounds[2 * D] = {-20, 20, -20, 20};
  if constexpr (D == 3) {
    bounds[4] = -20;
    bounds[5] = 20;
  }

  typename ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] =
        ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] =
      ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  // create a sphere in the level set
  NumericType origin[D] = {0., 0.};
  if constexpr (D == 3)
    origin[2] = 0;
  NumericType radius = 15.3;
  ls::MakeGeometry<NumericType, D>(
      substrate,
      ls::SmartPointer<ls::Sphere<NumericType, D>>::New(origin, radius))
      .apply();

  origin[0] = 15.0;
  radius = 8.7;
  auto secondSphere =
      ls::SmartPointer<ls::Domain<NumericType, D>>::New(substrate->getGrid());
  ls::MakeGeometry<NumericType, D>(
      secondSphere,
      ls::SmartPointer<ls::Sphere<NumericType, D>>::New(origin, radius))
      .apply();

  ls::BooleanOperation<NumericType, D>(substrate, secondSphere,
                                       ls::BooleanOperationEnum::UNION)
      .apply();

  ls::Expand<NumericType, D>(substrate, 3).apply();
  ls::Expand<NumericType, D>(secondSphere, 3).apply();

  auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
  ls::ToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
  ls::VTKWriter<NumericType>(mesh, "surface_1.vtk").apply();
  ls::ToSurfaceMesh<NumericType, D>(secondSphere, mesh).apply();
  ls::VTKWriter<NumericType>(mesh, "surface_2.vtk").apply();

  auto visualizeMesh =
      ls::SmartPointer<ls::WriteVisualizationMesh<NumericType, D>>::New();
  visualizeMesh->insertNextLevelSet(secondSphere);
  visualizeMesh->insertNextLevelSet(substrate);
  visualizeMesh->setExtractHullMesh(true);
  visualizeMesh->setFileName("myFile");

  visualizeMesh->apply();

  //   ls::BooleanOperation<NumericType, D>(substrate, secondSphere,
  //                                      ls::BooleanOperationEnum::UNION)
  //       .apply();

  //   auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  //   ls::ToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
  //   ls::VTKWriter<double>(mesh, "twoSpheres.vtk").apply();

  //   substrate->getDomain().print();

  // for(hrleConstSparseIterator<Domain<NumericType, D>::DomainType>
  // it(substrate->getDomain()); !it.isFinished(); ++it) {
  //   std::cout << it.getStartIndices() << std::endl;
  // }

  return 0;
}
