#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsExpand.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

int main() {
  omp_set_num_threads(8);

  constexpr int D = 3;
  typedef double NumericType;
  double gridDelta = 1.0;

  double extent = 50;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  if (D == 3) {
    bounds[4] = -extent;
    bounds[5] = extent;
  }

  typename lsDomain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] =
        lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] =
      lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  lsDomain<NumericType, D> levelSet(bounds, boundaryCons, gridDelta);
  // create a sphere in the level set
  NumericType origin[D] = {0., 0.};
  if (D == 3)
    origin[2] = 0;
  NumericType radius = 8.0;
  lsMakeGeometry<NumericType, D>(levelSet,
                                 lsSphere<NumericType, D>(origin, radius))
      .apply();

  lsDomain<NumericType, D> sphere2(bounds, boundaryCons, gridDelta);

  origin[1] = 10.0;
  lsMakeGeometry<NumericType, D>(sphere2,
                                 lsSphere<NumericType, D>(origin, radius))
      .apply();

  lsBooleanOperation<NumericType, D>(levelSet, sphere2,
                                     lsBooleanOperationEnum::UNION)
      .apply();

  lsMesh mesh;

  lsToMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "points.vtk").apply();
  lsToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "surface.vtk").apply();

  // set up spherical advection dist
  lsSphereDistribution<double, D> dist(20.0, gridDelta);
  lsGeometricAdvect<NumericType, D>(levelSet, dist).apply();

  lsToMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "afterDepoLS.vtk").apply();
  lsToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "afterDepo.vtk").apply();

  // now remove the same again using spherical distribution
  lsSphereDistribution<NumericType, D> etch(-20.0, gridDelta);
  lsGeometricAdvect<NumericType, D>(levelSet, etch).apply();

  lsToMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "afterEtchLS.vtk").apply();
  lsToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "afterEtch.vtk").apply();

  return 0;
}
