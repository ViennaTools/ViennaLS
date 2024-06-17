#include <algorithm>
#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMarkVoidPoints.hpp>
#include <lsRemoveStrayPoints.hpp>
#include <lsTestAsserts.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

int main() {
  constexpr int D = 2;
  // omp_set_num_threads(1);

  double extent = 10;
  double gridDelta = 1;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;

  boundaryCons[1] = ls::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin[D] = {0., 0.};
  double normal[D] = {0., 1.};

  ls::MakeGeometry<double, D>(
      substrate, ls::SmartPointer<ls::Plane<double, D>>::New(origin, normal))
      .apply();
  {
    auto hole = ls::SmartPointer<ls::Domain<double, D>>::New(
        bounds, boundaryCons, gridDelta);
    origin[1] = -5.;
    ls::MakeGeometry<double, D>(
        hole, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, 3.))
        .apply();

    ls::BooleanOperation<double, D>(
        substrate, hole, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  // {
  //   lsMarkVoidPoints<double, D> marker;
  //   marker.setLevelSet(substrate);
  //   marker.setVoidTopSurface(lsVoidTopSurfaceEnum::LEX_HIGHEST);
  //   marker.setSaveComponentIds(true);
  //   marker.apply();
  //   auto explMesh = ls::SmartPointer<ls::Mesh<>>::New();
  //   ls::ToMesh<double, D>(substrate, explMesh).apply();
  //   ls::VTKWriter<double>(explMesh, ls::FileFormatEnum::VTP,
  //   "before.vtp").apply();
  // }

  // Remove the stray points
  ls::RemoveStrayPoints<double, D> cleaner;
  cleaner.setLevelSet(substrate);
  cleaner.setVoidTopSurface(ls::VoidTopSurfaceEnum::LEX_HIGHEST);
  cleaner.apply();

  // check if the correct surface was removed
  LSTEST_ASSERT(substrate->getNumberOfPoints() == 42)

  // {
  //   auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  //   ls::ToMesh<double, D>(substrate, mesh).apply();
  //   ls::VTKWriter<double>(mesh, ls::FileFormatEnum::VTP,
  //   "after.vtp").apply();
  // }

  return 0;
}
