#include <iostream>

#include <lsAdvect.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  This example shows how to use lsAdvect to create an egg
  shape from a spherical level set using directional growth rates.
  \example Advection.cpp
*/

// implement own velocity field
class velocityField : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int /*material*/,
                           const std::array<double, 3> &normalVector) {
    // Some arbitrary velocity function of your liking
    // (try changing it and see what happens :)
    double velocity = 1. + ((normalVector[0] > 0) ? 2.3 : 0.5) *
                               std::abs(normalVector[0] * normalVector[0]);
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

  constexpr int D = 3;
  omp_set_num_threads(4);

  double gridDelta = 0.25;

  auto sphere1 = lsSmartPointer<lsDomain<double, D>>::New(gridDelta);

  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  lsMakeGeometry<double, D>(
      sphere1, lsSmartPointer<lsSphere<double, D>>::New(origin, radius))
      .apply();

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = lsSmartPointer<lsMesh>::New();
    lsToSurfaceMesh<double, D>(sphere1, mesh).apply();
    lsVTKWriter(mesh, "before.vtk").apply();
  }

  // instantiate velocities
  auto velocities = lsSmartPointer<velocityField>::New();

  std::cout << "Advecting" << std::endl;

  lsAdvect<double, D> advectionKernel;
  advectionKernel.insertNextLevelSet(sphere1);
  advectionKernel.setVelocityField(velocities);
  advectionKernel.setAdvectionTime(2.);
  advectionKernel.setIntegrationScheme(
      lsIntegrationSchemeEnum::ENGQUIST_OSHER_1ST_ORDER);

  advectionKernel.apply();
  double advectionSteps = advectionKernel.getNumberOfTimeSteps();
  std::cout << "Number of Advection steps taken: " << advectionSteps
            << std::endl;

  lsPrune<double, D>(sphere1).apply();
  lsExpand<double, D>(sphere1, 2).apply();

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = lsSmartPointer<lsMesh>::New();
    lsToSurfaceMesh<double, D>(sphere1, mesh).apply();
    mesh->print();
    lsVTKWriter(mesh, "after.vtk").apply();
  }

  return 0;
}
