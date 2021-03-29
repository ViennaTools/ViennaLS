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

class velocityField : public lsVelocityField<double> {
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
  if (D == 3) {
    bounds[4] = -10;
    bounds[5] = 10;
  }

  typename lsDomain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] =
        lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] =
      lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

  double origin[3] = {0., 0., 0.};
  double planeNormal[3] = {0., D == 2, D == 3};

  lsMakeGeometry<double, D>(
      substrate, lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal))
      .apply();

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter<double>(mesh, "plane.vtk").apply();
  }

  {
    // create layer used for booling
    std::cout << "Creating box..." << std::endl;
    auto trench = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons,
                                                           gridDelta);
    double minCorner[3] = {-extent - 1, -7.5, -15.};
    double maxCorner[3] = {extent + 1, 7.5, 1.0};
    if (D == 2) {
      minCorner[0] = minCorner[1];
      minCorner[1] = minCorner[2];
      maxCorner[0] = maxCorner[1];
      maxCorner[1] = maxCorner[2];
    }
    lsMakeGeometry<double, D>(
        trench, lsSmartPointer<lsBox<double, D>>::New(minCorner, maxCorner))
        .apply();

    {
      std::cout << "Extracting..." << std::endl;
      auto mesh = lsSmartPointer<lsMesh<>>::New();
      lsToMesh<double, D>(trench, mesh).apply();
      lsVTKWriter<double>(mesh, "box.vtk").apply();
    }

    // Create trench geometry
    std::cout << "Booling trench..." << std::endl;
    lsBooleanOperation<double, D>(substrate, trench,
                                  lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  auto mesh = lsSmartPointer<lsMesh<>>::New();

  lsToMesh<NumericType, D>(substrate, mesh).apply();
  lsVTKWriter<double>(mesh, "points.vtk").apply();
  lsToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
  lsVTKWriter<double>(mesh, "surface.vtk").apply();

  // Distance to advect to
  double depositionDistance = 4.0;

  // set up spherical advection dist
  auto dist = lsSmartPointer<lsSphereDistribution<NumericType, D>>::New(
      depositionDistance, gridDelta);

  auto newLayer = lsSmartPointer<lsDomain<double, D>>::New(substrate);

  std::cout << "GeometricAdvecting" << std::endl;
  lsGeometricAdvect<NumericType, D> fastAdvectKernel(newLayer, dist);

  {
    auto start = std::chrono::high_resolution_clock::now();
    fastAdvectKernel.apply();
    auto end = std::chrono::high_resolution_clock::now();
    auto diff =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << "Fast Advect: " << diff << "ms" << std::endl;
  }

  lsToSurfaceMesh<double, D>(newLayer, mesh).apply();
  lsVTKWriter<double>(mesh, "GeometricAdvect.vtk").apply();
  // lsToSurfaceMesh<double, D>(newLayer, mesh).apply();
  // lsVTKWriter<double>(mesh, "finalSurface.vtk").apply();

  // now rund lsAdvect for all other advection schemes
  // last scheme is SLLFS with i == 9
  for (unsigned i = 0; i < 10; ++i) {
    if (i == 4) {
      continue;
    }
    lsAdvect<double, D> advectionKernel;
    auto nextLayer = lsSmartPointer<lsDomain<double, D>>::New(substrate);
    advectionKernel.insertNextLevelSet(nextLayer);

    auto velocities = lsSmartPointer<velocityField>::New();
    advectionKernel.setVelocityField(velocities);
    advectionKernel.setAdvectionTime(depositionDistance);
    advectionKernel.setIntegrationScheme(
        static_cast<lsIntegrationSchemeEnum>(i));
    {
      auto start = std::chrono::high_resolution_clock::now();
      advectionKernel.apply();
      auto end = std::chrono::high_resolution_clock::now();
      auto diff =
          std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
              .count();
      std::cout << "Advect " << i << ": " << diff << "ms" << std::endl;
    }

    lsToSurfaceMesh<double, D>(nextLayer, mesh).apply();
    lsVTKWriter<double>(mesh, "Advect-" + std::to_string(i) + ".vtk").apply();
  }

  return 0;
}