#include <iostream>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsTransformMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example showing how to rotate an Mesh<>
  \example RotateMesh.cpp
*/

namespace ls = viennals;

int main() {
  using NumericType = double;
  constexpr int D = 3;

  omp_set_num_threads(4);

  auto levelSet = ls::SmartPointer<ls::Domain<NumericType, D>>::New();
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();

  const NumericType radius = 7.3;
  const ls::VectorType<NumericType, D> min(-50, -25., -25.);
  const ls::VectorType<NumericType, D> max(0., 0., 0.);
  ls::MakeGeometry<NumericType, D>(
      levelSet, ls::SmartPointer<ls::Box<double, D>>::New(min, max))
      .apply();

  // const hrleVectorType<NumericType, D> centre(5., 10., 0.);
  // ls::MakeGeometry<NumericType, D>(
  //     levelSet, ls::SmartPointer<ls::Sphere<double, D>>::New(centre, radius))
  //     .apply();

  ls::ToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();

  ls::VTKWriter<double>(mesh, "Initial.vtk").apply();

  // auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  // mesh->insertNextNode({1., 0., 1.});
  ls::VectorType<viennahrle::CoordType, 3> rotAxis{0., 1., 1.};
  ls::TransformMesh<double>(mesh, ls::TransformEnum::ROTATION, rotAxis, M_PI_4)
      .apply();

  ls::VTKWriter<double>(mesh, "Rotated.vtk").apply();

  return 0;
}
