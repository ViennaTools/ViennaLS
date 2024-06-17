#include <chrono>
#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsExpand.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsWriteVisualizationMesh.hpp>

namespace ls = viennals;

int main(int argc, char *argv[]) {
  constexpr int D = 3;
  typedef double NumericType;

  NumericType gridDelta = 1.e-2;
  NumericType depth = 1.2;
  NumericType topRadius = 2e-1;
  NumericType baseRadius = 1.75e-1;

  // Process parameters
  double extent = 0.5;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  if constexpr (D == 3) {
    bounds[4] = -extent;
    bounds[5] = extent;
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

  {
    NumericType origin[D] = {0.};
    NumericType planeNormal[D] = {0.};
    planeNormal[D - 1] = 1.;

    auto plane =
        ls::SmartPointer<ls::Plane<NumericType, D>>::New(origin, planeNormal);
    ls::MakeGeometry<NumericType, D>(substrate, plane).apply();
  }
  {
    std::cout << "Writing substrate" << std::endl;
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    ls::VTKWriter(mesh, ls::FileFormatEnum::VTP, "substrate.vtp").apply();
  }
  {
    // make LS from trench mesh and remove from substrate
    auto hole = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);

    NumericType origin[D] = {0., 0., -depth};
    NumericType axisDirection[D] = {0., 0., 1.};

    auto cylinder = ls::SmartPointer<ls::Cylinder<NumericType, D>>::New(
        origin, axisDirection, depth, baseRadius, topRadius);
    ls::MakeGeometry<NumericType, D>(hole, cylinder).apply();
    {
      std::cout << "Writing hole" << std::endl;
      auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
      ls::ToSurfaceMesh<NumericType, D>(hole, mesh).apply();
      ls::VTKWriter(mesh, ls::FileFormatEnum::VTP, "hole.vtp").apply();
    }

    ls::BooleanOperation<NumericType, D>(
        substrate, hole, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }
  {
    std::cout << "Writing output" << std::endl;
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    ls::VTKWriter(mesh, ls::FileFormatEnum::VTP, "surface_i.vtp").apply();
  }

  return 0;
}