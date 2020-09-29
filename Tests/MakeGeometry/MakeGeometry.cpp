#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

int main() {
  omp_set_num_threads(1);

  constexpr int D = 2;
  typedef double NumericType;

  double gridDelta = 1.0;

  double bounds[2 * D] = {-20, 20, -20, 20};
  if (D == 3) {
    bounds[4] = -20;
    bounds[5] = 20;
  }

  typename lsDomain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] =
        lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] =
      lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = lsSmartPointer<lsDomain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  // create a sphere in the level set
  NumericType origin[D] = {0., 0.};
  if (D == 3)
    origin[2] = 0;
  NumericType radius = 15.3;
  lsMakeGeometry<NumericType, D>(
      substrate, lsSmartPointer<lsSphere<NumericType, D>>::New(origin, radius))
      .apply();

  origin[0] = 15.0;
  radius = 8.7;
  auto secondSphere =
      lsSmartPointer<lsDomain<NumericType, D>>::New(substrate->getGrid());
  lsMakeGeometry<NumericType, D>(
      secondSphere,
      lsSmartPointer<lsSphere<NumericType, D>>::New(origin, radius))
      .apply();

  lsBooleanOperation<NumericType, D>(substrate, secondSphere,
                                     lsBooleanOperationEnum::UNION)
      .apply();

  auto mesh = lsSmartPointer<lsMesh>::New();
  lsToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
  lsVTKWriter(mesh, "twoSpheres.vtk").apply();

  substrate->getDomain().print();

  // for(hrleConstSparseIterator<lsDomain<NumericType, D>::DomainType>
  // it(substrate->getDomain()); !it.isFinished(); ++it) {
  //   std::cout << it.getStartIndices() << std::endl;
  // }

  return 0;
}
