#include <hrleDenseIterator.hpp>
#include <hrleVectorType.hpp>
#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsReduce.hpp>
#include <lsTestAsserts.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example for creating a plane.
  \example MakePlane.cpp
*/

namespace ls = viennals;

int main() {
  constexpr int D = 3;
  double extent = 15.7;
  double gridDelta = 2.8;

  double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
  ls::Domain<double, D>::BoundaryType boundaryCons[3];
  for (unsigned i = 0; i < D - 1; ++i)
    // boundaryCons[i] = ls::Domain<double, D>::BoundaryType::PERIODIC_BOUNDARY;
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  boundaryCons[2] = ls::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto levelSet = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();

  const hrleVectorType<double, D> origin(0., 0., 0.);
  // const hrleVectorType<double, D> normal(0., 0., 1.);
  const hrleVectorType<double, D> normal(1., 1., 1.);

  ls::MakeGeometry<double, D>(
      levelSet, ls::SmartPointer<ls::Plane<double, D>>::New(origin, normal))
      .apply();

  LSTEST_ASSERT_VALID_LS(levelSet, double, D)

  // ls::ToSurfaceMesh<double, D>(levelSet, mesh).apply();
  // ls::VTKWriter<double>(mesh, "Plane.vtk").apply();

  // ls::ToMesh<double, D>(levelSet, mesh).apply();
  // ls::VTKWriter<double>(mesh, "PlanePoints.vtk").apply();

  return 0;
}
