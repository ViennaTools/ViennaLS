#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsExpand.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

int main() {
  omp_set_num_threads(8);

  constexpr int D = 3;
  typedef double NumericType;
  double gridDelta = 2.0;

  double extent = 50;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  if constexpr (D == 3) {
    bounds[4] = -extent;
    bounds[5] = extent;
  }

  typename ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] =
        ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] =
      ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto mask = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  auto levelSet = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);
  ;

  // create mask
  {
    NumericType normal[3] = {0., (D == 2) ? 1. : 0., (D == 3) ? 1. : 0.};
    NumericType origin[3]{};
    ls::MakeGeometry<NumericType, D>(
        mask, ls::SmartPointer<ls::Plane<NumericType, D>>::New(origin, normal))
        .apply();
    normal[D - 1] = -1.0;
    origin[D - 1] = -extent / 5.0;
    auto maskBottom = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);
    ls::MakeGeometry<NumericType, D>(
        maskBottom,
        ls::SmartPointer<ls::Plane<NumericType, D>>::New(origin, normal))
        .apply();
    ls::BooleanOperation<NumericType, D>(mask, maskBottom,
                                         ls::BooleanOperationEnum::INTERSECT)
        .apply();

    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<NumericType, D>(mask, mesh).apply();
    ls::VTKWriter<double>(mesh, "Plane.vtk").apply();

    auto maskHole = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);

    double holeOrigin[3] = {0., 0., origin[D - 1] - 2 * gridDelta};
    double axis[3] = {0.0, 0.0, 1.0};
    ls::MakeGeometry<NumericType, D>(
        maskHole,
        ls::SmartPointer<ls::Cylinder<NumericType, D>>::New(
            holeOrigin, axis, 4 * gridDelta - origin[D - 1], extent / 1.5))
        .apply();

    // double minScalar = origin[D-1] - 2 * gridDelta;
    // double maxScalar = extent/5 * 2 * gridDelta;
    // double min[3] = {minScalar, minScalar, minScalar};
    // double max[3] = {maxScalar, maxScalar, maxScalar};
    // ls::MakeGeometry<NumericType, D>(maskHole,
    // ls::SmartPointer<ls::Box<NumericType, D>>::New(min, max)).apply();

    ls::BooleanOperation<NumericType, D>(
        mask, maskHole, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();

    ls::ToSurfaceMesh<NumericType, D>(mask, mesh).apply();
    ls::VTKWriter<double>(mesh, "Mask.vtk").apply();

    // make substrate
    ls::BooleanOperation<NumericType, D>(maskBottom,
                                         ls::BooleanOperationEnum::INVERT)
        .apply();
    levelSet->deepCopy(maskBottom);
    ls::BooleanOperation<NumericType, D>(levelSet, mask,
                                         ls::BooleanOperationEnum::UNION)
        .apply();
  }

  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();

  ls::ToMesh<NumericType, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "Surface_i_p.vtk").apply();
  ls::ToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "Surface_i.vtk").apply();
  ls::ToMesh<NumericType, D>(mask, mesh).apply();
  ls::VTKWriter<double>(mesh, "Surface_m_p.vtk").apply();
  ls::ToSurfaceMesh<NumericType, D>(mask, mesh).apply();
  ls::VTKWriter<double>(mesh, "Surface_m.vtk").apply();

  //   set up spherical advection dist
  auto dist =
      // ls::SmartPointer<ls::SphereDistribution<double,
      // D>>::New(-15.0,gridDelta);
      ls::SmartPointer<ls::BoxDistribution<NumericType, D>>::New(
          std::array<NumericType, 3>({-gridDelta, -gridDelta, -150.0}));

  ls::GeometricAdvect<NumericType, D>(levelSet, dist, mask).apply();

  ls::ToMesh<NumericType, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "afterDepoLS.vtk").apply();
  ls::ToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  ls::VTKWriter<double>(mesh, "afterDepo.vtk").apply();

  return 0;
}