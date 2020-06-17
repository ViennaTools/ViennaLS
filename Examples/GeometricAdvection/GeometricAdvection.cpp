#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  3D Example showing how to use the library for topography
  emulation, by creating a trench geometry. A uniform
  layer of a different material is then grown on top. It is
  the same example as Deposition but emulates the deposition
  rather than simulating a slow growth.
  \example GeometricAdvection.cpp
*/

int main() {

  constexpr int D = 3;
  omp_set_num_threads(4);

  double extent = 30;
  double gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[2] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  lsDomain<double, D> substrate(bounds, boundaryCons, gridDelta);

  double origin[3] = {0., 0., 0.};
  double planeNormal[3] = {0., 0., 1.};

  lsMakeGeometry<double, D>(substrate, lsPlane<double, D>(origin, planeNormal))
      .apply();

  lsDomain<double, D> trench(bounds, boundaryCons, gridDelta);
  // make -x and +x greater than domain for numerical stability
  double minCorner[D] = {-extent - 1, -extent / 4., -15.};
  double maxCorner[D] = {extent + 1, extent / 4., 1.};
  lsMakeGeometry<double, D>(trench, lsBox<double, D>(minCorner, maxCorner))
      .apply();

  // Create trench geometry
  lsBooleanOperation<double, D>(substrate, trench,
                                lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
      .apply();

  {
    std::cout << "Extracting..." << std::endl;
    lsMesh mesh;
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter(mesh, "trench-0.vtk").apply();
  }

  // Now grow new material isotropically

  // create new levelset for new material, which will be grown
  // since it has to wrap around the substrate, just copy it
  lsDomain<double, D> newLayer(substrate);

  std::cout << "Advecting" << std::endl;
  // Grow the layer uniformly by 4 as in deposition example
  lsSphereDistribution<double, D> dist(4.0);
  lsGeometricAdvect<double, D>(newLayer, dist).apply();

  lsMesh mesh;
  lsToSurfaceMesh<double, D>(newLayer, mesh).apply();
  lsVTKWriter(mesh, "trench-final.vtk").apply();

  return 0;
}
