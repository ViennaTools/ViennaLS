#include <chrono>
#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsExpand.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsVelocityField.hpp>

namespace ls = viennals;

class velocityField : public ls::VelocityField<double> {
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int /*material*/,
                           const std::array<double, 3> & /*normalVector*/,
                           unsigned long /*pointId*/) {
    return 1;
  }
};

int main() {
  constexpr int D = 3;
  typedef double NumericType;

  double extent = 10;
  double gridDelta = 0.25;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  if constexpr (D == 3) {
    bounds[4] = -10;
    bounds[5] = 10;
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
    double minCorner[3] = {-extent - 1, -7.5, -15.};
    double maxCorner[3] = {extent + 1, 7.5, 1.0};
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

  // Distance to advect to
  double depositionDistance = 4.0;

  // set up spherical advection dist
  auto dist = ls::SmartPointer<ls::SphereDistribution<NumericType, D>>::New(
      depositionDistance, gridDelta);

  auto newLayer = ls::SmartPointer<ls::Domain<double, D>>::New(substrate);

  std::cout << "GeometricAdvecting" << std::endl;
  ls::GeometricAdvect<NumericType, D> fastAdvectKernel(newLayer, dist);

  {
    auto start = std::chrono::high_resolution_clock::now();
    fastAdvectKernel.apply();
    auto end = std::chrono::high_resolution_clock::now();
    auto diff =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << "Fast Advect: " << diff << "ms" << std::endl;
  }

  ls::ToSurfaceMesh<double, D>(newLayer, mesh).apply();
  ls::VTKWriter<double>(mesh, "GeometricAdvect.vtk").apply();
  // ls::ToSurfaceMesh<double, D>(newLayer, mesh).apply();
  // ls::VTKWriter<double>(mesh, "finalSurface.vtk").apply();

  // now rund lsAdvect for all other advection schemes
  // last scheme is SLLFS with i == 9
  for (unsigned i = 0; i < 10; ++i) {
    if (i == 4) {
      continue;
    }
    ls::Advect<double, D> advectionKernel;
    auto nextLayer = ls::SmartPointer<ls::Domain<double, D>>::New(substrate);
    advectionKernel.insertNextLevelSet(nextLayer);

    auto velocities = ls::SmartPointer<velocityField>::New();
    advectionKernel.setVelocityField(velocities);
    advectionKernel.setAdvectionTime(depositionDistance);
    advectionKernel.setIntegrationScheme(
        static_cast<ls::IntegrationSchemeEnum>(i));
    {
      auto start = std::chrono::high_resolution_clock::now();
      advectionKernel.apply();
      auto end = std::chrono::high_resolution_clock::now();
      auto diff =
          std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
              .count();
      std::cout << "Advect " << i << ": " << diff << "ms" << std::endl;
    }

    ls::ToSurfaceMesh<double, D>(nextLayer, mesh).apply();
    ls::VTKWriter<double>(mesh, "Advect-" + std::to_string(i) + ".vtk").apply();
  }

  return 0;
}