#include <algorithm>
#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMarkVoidPoints.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsRemoveStrayPoints.hpp>

int main() {
  constexpr int D = 2;
  omp_set_num_threads(1);

  double extent = 10;
  double gridDelta = 1;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  boundaryCons[1] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

  double origin[D] = {0., 0.};
  double normal[D] = {0., 1.};

  lsMakeGeometry<double, D>(
      substrate, lsSmartPointer<lsPlane<double, D>>::New(origin, normal))
      .apply();
  {
    auto hole = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons,
                                                         gridDelta);
    origin[1] = -5.;
    lsMakeGeometry<double, D>(
        hole, lsSmartPointer<lsSphere<double, D>>::New(origin, 3.))
        .apply();

    lsBooleanOperation<double, D>(substrate, hole,
                                  lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  {
    auto explMesh = lsSmartPointer<lsMesh<>>::New();

    std::cout << "Extracting..." << std::endl;
    lsToSurfaceMesh<double, D>(substrate, explMesh).apply();

    lsVTKWriter<double>(explMesh, "before.vtk").apply();
  }

  lsMarkVoidPoints<double, D>(substrate).apply();

  // Remove the stray points
  lsRemoveStrayPoints<double, D>(substrate).apply();

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter<double>(mesh, "after.vtk").apply();
  }

  return 0;
}
