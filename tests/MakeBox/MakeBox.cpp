#include <chrono>
#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsTestAsserts.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example showing how to write different
  meshes created by the algorithms lsToVoxelMesh and lsToSurfaceMesh.
  \example Make3DSphere.cpp
*/

namespace ls = viennals;

int main() {

  constexpr int D = 3;

  double gridDelta = 0.5;

  double extent = 30;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  if constexpr (D == 3) {
    bounds[4] = -extent;
    bounds[5] = extent;
  }

  typename ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::PERIODIC_BOUNDARY;
  }
  boundaryCons[D - 1] = ls::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto levelSet = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  ls::VectorType<double, D> min{bounds[0] + 10, bounds[2] + extent, 0};
  ls::VectorType<double, D> max{bounds[1] - 10, bounds[3] + extent / 2., 4};
  ls::MakeGeometry<double, D>(
      levelSet, ls::SmartPointer<ls::Box<double, D>>::New(min, max))
      .apply();

  LSTEST_ASSERT_VALID_LS(levelSet, double, D)

  // std::cout << "Number of points: " <<
  // levelSet->getDomain().getNumberOfPoints()
  //           << std::endl;
  // auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  // ls::ToMesh<double, D>(levelSet, mesh).apply();
  // ls::VTKWriter<double>(mesh, "boxLS.vtk").apply();
  // ls::ToSurfaceMesh<double, D>(levelSet, mesh).apply();
  // ls::VTKWriter<double>(mesh, "box.vtk").apply();

  return 0;
}
