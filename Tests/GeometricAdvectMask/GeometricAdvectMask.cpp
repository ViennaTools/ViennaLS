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
  omp_set_num_threads(4);

  constexpr int D = 2;
  typedef double NumericType;
  double gridDelta = 1.0;

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

  auto levelSet =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);
  // create a sphere in the level set
  NumericType normal[3] = {0., (D==2)?1.:0., (D==3)?1.:0.};
  NumericType origin[3] = {};
  lsMakeGeometry<NumericType, D>(
      levelSet, lsSmartPointer<lsPlane<NumericType, D>>::New(origin, normal))
      .apply();

  auto mask =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

  {
      NumericType extent4 = extent / 4.0;
      NumericType min[3] = {-extent4, -extent4, -extent4};
      NumericType max[3] = {extent4, extent4, extent4};

      lsMakeGeometry<NumericType, D>(
      mask, lsSmartPointer<lsBox<NumericType, D>>::New(min, max))
      .apply();
  }

  lsBooleanOperation<NumericType, D>(levelSet, mask,
                                     lsBooleanOperationEnum::UNION)
      .apply();

  auto mesh = lsSmartPointer<lsMesh>::New();

  lsToMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "Surface_i_p.vtk").apply();
  lsToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "Surface_i.vtk").apply();
  lsToMesh<NumericType, D>(mask, mesh).apply();
  lsVTKWriter(mesh, "Surface_m_p.vtk").apply();
  lsToSurfaceMesh<NumericType, D>(mask, mesh).apply();
  lsVTKWriter(mesh, "Surface_m.vtk").apply();

//   set up spherical advection dist
  auto dist =
      // lsSmartPointer<lsSphereDistribution<double, D>>::New(-15.0, gridDelta);
      lsSmartPointer<lsBoxDistribution<NumericType, D>>::New(std::array<NumericType, 3>({-gridDelta, -15.0, -gridDelta}), gridDelta);
  lsGeometricAdvect<NumericType, D>(levelSet, dist, mask).apply();

  lsToMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "afterDepoLS.vtk").apply();
  lsToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  lsVTKWriter(mesh, "afterDepo.vtk").apply();

  return 0;
}