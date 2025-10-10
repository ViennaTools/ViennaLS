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
  omp_set_num_threads(12);

  constexpr int D = 2;
  typedef double NumericType;
  double gridDelta = 1.1;

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

  auto substrate = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin[3] = {0., 0., 0.};
  double planeNormal[3] = {0., D == 2, D == 3};

  ls::MakeGeometry<double, D>(
      substrate,
      ls::SmartPointer<ls::Plane<double, D>>::New(origin, planeNormal))
      .apply();

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(substrate, mesh).apply();
    ls::VTKWriter<double>(mesh, "plane.vtk").apply();
  }

  {
    // create layer used for booling
    std::cout << "Creating box..." << std::endl;
    auto trench = ls::SmartPointer<ls::Domain<double, D>>::New(
        bounds, boundaryCons, gridDelta);
    double minCorner[3] = {-extent - 1, -extent / 4., -15.};
    double maxCorner[3] = {extent + 1, extent / 4., 1.0};
    if (D == 2) {
      minCorner[0] = minCorner[1];
      minCorner[1] = minCorner[2];
      maxCorner[0] = maxCorner[1];
      maxCorner[1] = maxCorner[2];
    }
    ls::MakeGeometry<double, D>(
        trench, ls::SmartPointer<ls::Box<double, D>>::New(minCorner, maxCorner))
        .apply();

    {
      std::cout << "Extracting..." << std::endl;
      auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
      ls::ToMesh<double, D>(trench, mesh).apply();
      ls::VTKWriter<double>(mesh, "box.vtk").apply();
    }

    // Create trench geometry
    std::cout << "Booling trench..." << std::endl;
    ls::BooleanOperation<double, D>(
        substrate, trench, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();

  ls::ToMesh<NumericType, D>(substrate, mesh).apply();
  ls::VTKWriter<double>(mesh, "points.vtk").apply();
  ls::ToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
  ls::VTKWriter<double>(mesh, "surface.vtk").apply();

  // set up spherical advection dist
  // lsSphereDistribution<NumericType, D> dist(15.0);
  std::cout << "Advecting..." << std::endl;
  std::array<NumericType, 3> box = {1.1, 15};
  if constexpr (D == 3) {
    box[1] = 1.1;
    box[2] = 15;
  }
  auto dist = ls::SmartPointer<ls::BoxDistribution<NumericType, D>>::New(box);
  ls::GeometricAdvect<NumericType, D>(substrate, dist).apply();

  std::cout << "Writing results..." << std::endl;
  ls::ToMesh<NumericType, D>(substrate, mesh).apply();
  ls::VTKWriter<double>(mesh, "finalLS.vtk").apply();

  ls::ToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
  ls::VTKWriter<double>(mesh, "finalSurface.vtk").apply();

  std::cout << "Done" << std::endl;

  return 0;
}
