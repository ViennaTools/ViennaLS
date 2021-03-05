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
  double gridDelta = 2.0;

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

  auto mask = lsSmartPointer<lsDomain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  auto levelSet = lsSmartPointer<lsDomain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);
  ;

  // create mask
  {
    NumericType normal[3] = {0., (D == 2) ? 1. : 0., (D == 3) ? 1. : 0.};
    NumericType origin[3] = {};
    lsMakeGeometry<NumericType, D>(
        mask, lsSmartPointer<lsPlane<NumericType, D>>::New(origin, normal))
        .apply();
    normal[D - 1] = -1.0;
    origin[D - 1] = -extent / 5.0;
    auto maskBottom = lsSmartPointer<lsDomain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);
    lsMakeGeometry<NumericType, D>(
        maskBottom,
        lsSmartPointer<lsPlane<NumericType, D>>::New(origin, normal))
        .apply();
    lsBooleanOperation<NumericType, D>(mask, maskBottom,
                                       lsBooleanOperationEnum::INTERSECT)
        .apply();

    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToMesh<NumericType, D>(mask, mesh).apply();
    lsVTKWriter<double>(mesh, "Plane.vtk").apply();

    auto maskHole = lsSmartPointer<lsDomain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);

    double holeOrigin[3] = {0., 0., origin[D - 1] - 2 * gridDelta};
    double axis[3] = {0.0, 0.0, 1.0};
    lsMakeGeometry<NumericType, D>(
        maskHole,
        lsSmartPointer<lsCylinder<NumericType, D>>::New(
            holeOrigin, axis, 4 * gridDelta - origin[D - 1], extent / 1.5))
        .apply();

    // double minScalar = origin[D-1] - 2 * gridDelta;
    // double maxScalar = extent/5 * 2 * gridDelta;
    // double min[3] = {minScalar, minScalar, minScalar};
    // double max[3] = {maxScalar, maxScalar, maxScalar};
    // lsMakeGeometry<NumericType, D>(maskHole,
    // lsSmartPointer<lsBox<NumericType, D>>::New(min, max)).apply();

    lsBooleanOperation<NumericType, D>(
        mask, maskHole, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();

    lsToSurfaceMesh<NumericType, D>(mask, mesh).apply();
    lsVTKWriter<double>(mesh, "Mask.vtk").apply();

    // make substrate
    lsBooleanOperation<NumericType, D>(maskBottom,
                                       lsBooleanOperationEnum::INVERT)
        .apply();
    levelSet->deepCopy(maskBottom);
    lsBooleanOperation<NumericType, D>(levelSet, mask,
                                       lsBooleanOperationEnum::UNION)
        .apply();
  }

  auto mesh = lsSmartPointer<lsMesh<>>::New();

  lsToMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter<double>(mesh, "Surface_i_p.vtk").apply();
  lsToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter<double>(mesh, "Surface_i.vtk").apply();
  lsToMesh<NumericType, D>(mask, mesh).apply();
  lsVTKWriter<double>(mesh, "Surface_m_p.vtk").apply();
  lsToSurfaceMesh<NumericType, D>(mask, mesh).apply();
  lsVTKWriter<double>(mesh, "Surface_m.vtk").apply();

  //   set up spherical advection dist
  auto dist =
      // lsSmartPointer<lsSphereDistribution<double, D>>::New(-15.0,gridDelta);
      lsSmartPointer<lsBoxDistribution<NumericType, D>>::New(
          std::array<NumericType, 3>({-gridDelta, -gridDelta, -150.0}),
          gridDelta);

  lsGeometricAdvect<NumericType, D>(levelSet, dist, mask).apply();

  lsToMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter<double>(mesh, "afterDepoLS.vtk").apply();
  lsToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter<double>(mesh, "afterDepo.vtk").apply();

  return 0;
}