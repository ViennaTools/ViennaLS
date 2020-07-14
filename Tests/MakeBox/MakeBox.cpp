#include <chrono>
#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example showing how to write different
  meshes created by the algorithms lsToVoxelMesh and lsToSurfaceMesh.
  \example Make3DSphere.cpp
*/

int main() {

  constexpr int D = 3;

  double gridDelta = 0.5;

  double extent = 30;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  if (D == 3) {
    bounds[4] = -extent;
    bounds[5] = extent;
  }

  typename lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto levelSet =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

  hrleVectorType<double, D> min(-10, -10, 0);
  hrleVectorType<double, D> max(10, 10, 4);
  lsMakeGeometry<double, D>(levelSet,
                            lsSmartPointer<lsBox<double, D>>::New(min, max))
      .apply();

  std::cout << "Number of points: " << levelSet->getDomain().getNumberOfPoints()
            << std::endl;
  auto mesh = lsSmartPointer<lsMesh>::New();
  lsToMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "boxLS.vtk").apply();
  lsToSurfaceMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "box.vtk").apply();

  return 0;
}
