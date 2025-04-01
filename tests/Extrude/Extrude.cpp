#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExtrude.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Simple example showing how to use lsExtrude in order
  to transform a 2D trench into a 3D level set.
  \example Extrude.cpp
*/

namespace ls = viennals;

int main() {

  omp_set_num_threads(4);

  double extent = 15;
  double gridDelta = 0.5;

  // 2D domain boundaries
  double bounds[2 * 2] = {-extent, extent, -extent, extent};
  ls::Domain<double, 2>::BoundaryType boundaryCons[2];
  boundaryCons[0] = ls::Domain<double, 2>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[1] = ls::Domain<double, 2>::BoundaryType::INFINITE_BOUNDARY;

  auto trench = ls::SmartPointer<ls::Domain<double, 2>>::New(
      bounds, boundaryCons, gridDelta);

  {
    double origin[2] = {0., 0.};
    double normal[2] = {0., 1.};

    ls::MakeGeometry<double, 2>(
        trench, ls::SmartPointer<ls::Plane<double, 2>>::New(origin, normal))
        .apply();

    auto cutOut = ls::SmartPointer<ls::Domain<double, 2>>::New(
        bounds, boundaryCons, gridDelta);

    double minPoint[2] = {-5., -5.};
    double maxPoint[2] = {5., gridDelta};

    ls::MakeGeometry<double, 2>(
        cutOut, ls::SmartPointer<ls::Box<double, 2>>::New(minPoint, maxPoint))
        .apply();

    ls::BooleanOperation<double, 2>(
        trench, cutOut, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, 2>(trench, mesh).apply();
    ls::VTKWriter<double>(mesh, "trench_initial.vtp").apply();
  }

  ls::Vec2D<double> extrudeExtent{-5., 5.};
  std::array<ls::BoundaryConditionEnum, 3> boundaryConds{
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY,
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY};
  auto trench_3D = ls::SmartPointer<ls::Domain<double, 3>>::New();
  ls::Extrude<double>(trench, trench_3D, extrudeExtent, 2, boundaryConds)
      .apply();

  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, 3>(trench_3D, mesh).apply();
    ls::VTKWriter<double>(mesh, "trench_extrude.vtp").apply();
  }

  return 0;
}
