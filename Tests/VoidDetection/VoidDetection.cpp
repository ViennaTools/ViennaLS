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
    auto explMesh = lsSmartPointer<lsMesh>::New();

    std::cout << "Extracting..." << std::endl;
    lsToSurfaceMesh<double, D>(substrate, explMesh).apply();

    lsVTKWriter(explMesh, "before.vtk").apply();
  }

  lsMarkVoidPoints<double, D>(substrate).apply();

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = lsSmartPointer<lsMesh>::New();
    lsToMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter(mesh, "after.vtk").apply();
  }

  // Advection
  auto velocities = lsSmartPointer<velocityField>::New();
  lsAdvect<double, D> advectionKernel(substrate, velocities);
  advectionKernel.setIgnoreVoids(true);
  for (unsigned i = 0; i < 30; ++i) {
    {
      auto mesh = lsSmartPointer<lsMesh>::New();
      lsToSurfaceMesh<double, D>(substrate, mesh).apply();
      lsVTKWriter(mesh, "out-" + std::to_string(i) + ".vtk").apply();

      lsMarkVoidPoints<double, D>(substrate).apply();
      lsToMesh<double, D>(substrate, mesh).apply();

      lsVTKWriter(mesh, "ls-out-" + std::to_string(i) + ".vtk").apply();
    }
    advectionKernel.apply();
  }

  return 0;
}
