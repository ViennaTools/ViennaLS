#include <hrleDenseIterator.hpp>
#include <hrleVectorType.hpp>
#include <iostream>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsReduce.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example for creating a plane.
  \example MakePlane.cpp
*/

int main() {
  constexpr int D = 3;
  double extent = 15.7;

  double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[3];
  for (unsigned i = 0; i < D - 1; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  boundaryCons[2] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  lsDomain<double, D> levelSet(bounds, boundaryCons, 0.7);
  lsMesh mesh;

  const hrleVectorType<double, D> origin(0., 0., 0.);
  const hrleVectorType<double, D> normal(1., 1., 1.);

  lsMakeGeometry<double, D>(levelSet, lsPlane<double, D>(origin, normal))
      .apply();

  lsToSurfaceMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "Plane.vtk").apply();

  lsToMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "PlanePoints.vtk").apply();

  return 0;
}
