#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsExpand.hpp>
#include <lsFastAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

int main() {
  constexpr int D = 3;
  typedef double NumericType;

  double extent = 30;
  double gridDelta = 0.5;
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

  lsDomain<double, D> substrate(bounds, boundaryCons, gridDelta);

  double origin[3] = {0., 0., 0.};
  double planeNormal[3] = {0., D == 2, D == 3};

  lsMakeGeometry<double, D>(substrate, lsPlane<double, D>(origin, planeNormal))
    .apply();

  {
    std::cout << "Extracting..." << std::endl;
    lsMesh mesh;
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter(mesh, "plane.vtk").apply();
  }

  {
    // create layer used for booling
    std::cout << "Creating box..." << std::endl;
    lsDomain<double, D> trench(bounds, boundaryCons, gridDelta);
    double minCorner[3] = {-extent - 1, -extent / 4., -15.};
    double maxCorner[3] = {extent + 1, extent / 4., 1.0};
    if(D==2) {
      minCorner[0] = minCorner[1];
      minCorner[1] = minCorner[2];
      maxCorner[0] = maxCorner[1];
      maxCorner[1] = maxCorner[2];
    }
    lsMakeGeometry<double, D>(trench, lsBox<double, D>(minCorner, maxCorner))
        .apply();

    {
      std::cout << "Extracting..." << std::endl;
      lsMesh mesh;
      lsToMesh<double, D>(trench, mesh).apply();
      lsVTKWriter(mesh, "box.vtk").apply();
    }

    // Create trench geometry
    std::cout << "Booling trench..." << std::endl;
    lsBooleanOperation<double, D>(substrate, trench,
                                  lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  lsMesh mesh;

  lsToMesh<NumericType, D>(substrate, mesh).apply();
  lsVTKWriter(mesh, "points.vtk").apply();
  lsToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
  lsVTKWriter(mesh, "surface.vtk").apply();

  // set up spherical advection dist
  lsSphereDistribution<NumericType, D> dist(4);

  lsDomain<double, D> newLayer(substrate);

  std::cout << "FastAdvecting" << std::endl;
  lsFastAdvect<NumericType, D> fastAdvectKernel(newLayer, dist);
  fastAdvectKernel.apply();

  lsToMesh<double, D>(newLayer, mesh).apply();
  lsVTKWriter(mesh, "FastAdvect.vtk").apply();
  lsToSurfaceMesh<double, D>(newLayer, mesh).apply();
  lsVTKWriter(mesh, "finalSurface.vtk").apply();

  return 0;
}