#include <algorithm>
#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMarkVoidPoints.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Example showing how to use void detection.
  \example VoidDetection.cpp
*/

// implement own velocity field
class velocityField : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int /*material*/,
                           const std::array<double, 3> & /*normalVector*/) {
    // Some arbitrary velocity function of your liking
    // (try changing it and see what happens :)
    double velocity = 1.;
    return velocity;
  }

  std::array<double, 3>
  getVectorVelocity(const std::array<double, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<double, 3> & /*normalVector*/) {
    return std::array<double, 3>({});
  }
};

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
  lsDomain<double, D> substrate(bounds, boundaryCons, gridDelta);

  double origin[D] = {0., 0.};
  double normal[D] = {0., 1.};

  lsMakeGeometry<double, D>(substrate, lsPlane<double, D>(origin, normal))
      .apply();
  {
    lsDomain<double, D> hole(bounds, boundaryCons, gridDelta);
    origin[1] = -5.;
    lsMakeGeometry<double, D>(hole, lsSphere<double, D>(origin, 3.)).apply();

    lsBooleanOperation<double, D>(substrate, hole,
                                  lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  {
    lsMesh explMesh;

    std::cout << "Extracting..." << std::endl;
    lsToSurfaceMesh<double, D>(substrate, explMesh).apply();

    lsVTKWriter(explMesh, "before.vtk").apply();
  }

  lsMarkVoidPoints<double, D>(substrate).apply();

  {
    std::cout << "Extracting..." << std::endl;
    lsMesh mesh;
    lsToMesh<double, D>(substrate, mesh).apply();

    auto voidPointMarkers = substrate.getVoidPointMarkers();
    std::vector<double> isVoid(
        voidPointMarkers.size()); // 0 = not void, 1 = void
    for (unsigned i = 0; i < isVoid.size(); ++i) {
      isVoid[i] = (voidPointMarkers[i]) ? 1. : 0.;
    }

    std::cout << "Points: " << substrate.getNumberOfPoints() << std::endl;
    std::cout << "Markers: " << isVoid.size() << std::endl;

    mesh.insertNextScalarData(isVoid, "voidMarkers");

    lsVTKWriter(mesh, "after.vtk").apply();
  }

  // Advection
  velocityField velocities;
  lsAdvect<double, D> advectionKernel(substrate, velocities);
  advectionKernel.setIgnoreVoids(true);
  for (unsigned i = 0; i < 30; ++i) {
    {
      lsMesh mesh;
      lsToSurfaceMesh<double, D>(substrate, mesh).apply();
      lsVTKWriter(mesh, "out-" + std::to_string(i) + ".vtk").apply();

      lsMarkVoidPoints<double, D>(substrate).apply();
      auto voidPointMarkers = substrate.getVoidPointMarkers();
      std::vector<double> isVoid(
          voidPointMarkers.size()); // 0 = not void, 1 = void
      for (unsigned i = 0; i < isVoid.size(); ++i) {
        isVoid[i] = (voidPointMarkers[i]) ? 1. : 0.;
      }
      lsToMesh<double, D>(substrate, mesh).apply();
      mesh.insertNextScalarData(isVoid, "voidMarkers");

      lsVTKWriter(mesh, "ls-out-" + std::to_string(i) + ".vtk").apply();
    }
    advectionKernel.apply();
  }

  return 0;
}
