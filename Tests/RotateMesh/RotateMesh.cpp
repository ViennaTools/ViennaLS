#include <iostream>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsTransformMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example showing how to rotate an lsMesh<>
  \example RotateMesh.cpp
*/

int main() {
  using NumericType = double;
  constexpr int D = 3;

  omp_set_num_threads(4);

  auto levelSet = lsSmartPointer<lsDomain<NumericType, D>>::New();
  auto mesh = lsSmartPointer<lsMesh<>>::New();

  const NumericType radius = 7.3;
  const hrleVectorType<NumericType, D> min(-50, -25., -25.);
  const hrleVectorType<NumericType, D> max(0., 0., 0.);
  lsMakeGeometry<NumericType, D>(
      levelSet, lsSmartPointer<lsBox<double, D>>::New(min, max))
      .apply();

  // const hrleVectorType<NumericType, D> centre(5., 10., 0.);
  // lsMakeGeometry<NumericType, D>(
  //     levelSet, lsSmartPointer<lsSphere<double, D>>::New(centre, radius))
  //     .apply();

  lsToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();

  lsVTKWriter<double>(mesh, "Initial.vtk").apply();

  // auto mesh = lsSmartPointer<lsMesh<>>::New();
  // mesh->insertNextNode({1., 0., 1.});
  hrleVectorType<hrleCoordType, 3> rotAxis{0., 1., 1.};
  lsTransformMesh<double>(mesh, lsTransformEnum::ROTATION, rotAxis, M_PI_4)
      .apply();

  lsVTKWriter<double>(mesh, "Rotated.vtk").apply();

  return 0;
}
