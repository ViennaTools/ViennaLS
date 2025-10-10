#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsExpand.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

int main() {
  omp_set_num_threads(8);

  constexpr int D = 3;
  typedef double NumericType;
  double gridDelta = 1.0;

  double extent = 50;
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

  auto levelSet = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  // create a sphere in the level set
  NumericType origin[D] = {0., 0.};
  if constexpr (D == 3)
    origin[2] = 0;
  NumericType radius = 8.0;
  ls::MakeGeometry<NumericType, D>(
      levelSet,
      ls::SmartPointer<ls::Sphere<NumericType, D>>::New(origin, radius))
      .apply();

  auto sphere = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  origin[1] = 10.0;
  ls::MakeGeometry<NumericType, D>(
      sphere, ls::SmartPointer<ls::Sphere<NumericType, D>>::New(origin, radius))
      .apply();

  ls::BooleanOperation<NumericType, D>(levelSet, sphere,
                                       ls::BooleanOperationEnum::UNION)
      .apply();

  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();

  ls::ToMesh<NumericType, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "points.vtk").apply();
  ls::ToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "surface.vtk").apply();

  // set up spherical advection dist
  auto dist = ls::SmartPointer<ls::SphereDistribution<double, D>>::New(20.0);
  ls::GeometricAdvect<NumericType, D>(levelSet, dist).apply();

  ls::ToMesh<NumericType, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "afterDepoLS.vtk").apply();
  ls::ToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "afterDepo.vtk").apply();

  // now remove the same again using spherical distribution
  auto etch =
      ls::SmartPointer<ls::SphereDistribution<NumericType, D>>::New(-20.0);
  ls::GeometricAdvect<NumericType, D>(levelSet, etch).apply();

  ls::ToMesh<NumericType, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "afterEtchLS.vtk").apply();
  ls::ToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "afterEtch.vtk").apply();

  return 0;
}
