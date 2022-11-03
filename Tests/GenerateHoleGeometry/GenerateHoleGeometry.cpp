#include <chrono>
#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsExpand.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsWriteVisualizationMesh.hpp>

int main(int argc, char *argv[]) {
  constexpr int D = 3;
  typedef double NumericType;

  NumericType gridDelta = 1.e-2;
  NumericType depth = 1.2;
  NumericType topRadius = 2e-1;
  NumericType baseRadius = 1.75e-1;

  // Process parameters
  double extent = 0.5;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  if constexpr (D == 3) {
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

  auto substrate = lsSmartPointer<lsDomain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  {
    NumericType origin[D] = {0.};
    NumericType planeNormal[D] = {0.};
    planeNormal[D - 1] = 1.;

    auto plane =
        lsSmartPointer<lsPlane<NumericType, D>>::New(origin, planeNormal);
    lsMakeGeometry<NumericType, D>(substrate, plane).apply();
  }
  {
    std::cout << "Writing substrate" << std::endl;
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    lsVTKWriter(mesh, lsFileFormatEnum::VTP, "substrate.vtp").apply();
  }
  {
    // make LS from trench mesh and remove from substrate
    auto hole = lsSmartPointer<lsDomain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);

    NumericType origin[D] = {0., 0., -depth};
    NumericType axisDirection[D] = {0., 0., 1.};

    auto cylinder = lsSmartPointer<lsCylinder<NumericType, D>>::New(
        origin, axisDirection, depth, baseRadius, topRadius);
    lsMakeGeometry<NumericType, D>(hole, cylinder).apply();
    {
    std::cout << "Writing hole" << std::endl;
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<NumericType, D>(hole, mesh).apply();
    lsVTKWriter(mesh, lsFileFormatEnum::VTP, "hole.vtp").apply();
  }

    lsBooleanOperation<NumericType, D>(
        substrate, hole, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }
  {
    std::cout << "Writing output" << std::endl;
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    lsVTKWriter(mesh, lsFileFormatEnum::VTP, "surface_i.vtp").apply();
  }

  return 0;
}