#include <array>
#include <iostream>

#include <lsAdvect.hpp>
#include <lsCompareChamfer.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsTestAsserts.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsVelocityField.hpp>
#include <vcTimer.hpp>

namespace ls = viennals;

// Define a constant velocity field
template <class T> class ConstantVelocity : public ls::VelocityField<T> {
  std::array<T, 3> velocity;

public:
  ConstantVelocity(std::array<T, 3> vel) : velocity(vel) {}

  T getScalarVelocity(const std::array<T, 3> & /*coordinate*/, int /*material*/,
                      const std::array<T, 3> & /*normalVector*/,
                      unsigned long /*pointId*/) override {
    return 0;
  }

  std::array<T, 3> getVectorVelocity(const std::array<T, 3> & /*coordinate*/,
                                     int /*material*/,
                                     const std::array<T, 3> & /*normalVector*/,
                                     unsigned long /*pointId*/) override {
    return velocity;
  }
};

int main() {
  // Define dimension and type
  constexpr int D = 3;
  using T = double;
  omp_set_num_threads(8);

  // Define grid and domain bounds
  double gridDelta = 0.1;
  double bounds[2 * D] = {-5.0, 5.0, -5.0, 5.0, -5.0, 5.0};
  ls::Domain<T, D>::BoundaryType boundaryCons[D];
  for (int i = 0; i < D; ++i)
    boundaryCons[i] = ls::Domain<T, D>::BoundaryType::INFINITE_BOUNDARY;

  // Create initial level set (Sphere)
  auto sphere =
      ls::SmartPointer<ls::Domain<T, D>>::New(bounds, boundaryCons, gridDelta);
  T origin[3]{};
  T radius = 1.5;
  ls::MakeGeometry<T, D>(sphere, ls::Sphere<T, D>::New(origin, radius)).apply();

  // Create analytical reference for error calculation at t=2.0.
  // The sphere starts at x=0 and moves with v=1 for t=2, so it ends at x=2.
  auto sphereRef =
      ls::SmartPointer<ls::Domain<T, D>>::New(bounds, boundaryCons, gridDelta);
  T originRef[3] = {2.0, 0.0, 0.0};
  ls::MakeGeometry<T, D>(sphereRef, ls::Sphere<T, D>::New(originRef, radius))
      .apply();

  // Create copies for Forward Euler, RK2 and RK3
  auto sphereFE = ls::SmartPointer<ls::Domain<T, D>>::New(sphere);
  auto sphereRK2 = ls::SmartPointer<ls::Domain<T, D>>::New(sphere);
  auto sphereRK3 = ls::SmartPointer<ls::Domain<T, D>>::New(sphere);
  auto sphereWENO3_FE = ls::SmartPointer<ls::Domain<T, D>>::New(sphere);
  auto sphereWENO3_RK2 = ls::SmartPointer<ls::Domain<T, D>>::New(sphere);
  auto sphereWENO3_RK3 = ls::SmartPointer<ls::Domain<T, D>>::New(sphere);
  auto sphereWENO5_FE = ls::SmartPointer<ls::Domain<T, D>>::New(sphere);
  auto sphereWENO5_RK2 = ls::SmartPointer<ls::Domain<T, D>>::New(sphere);
  auto sphereWENO5_RK3 = ls::SmartPointer<ls::Domain<T, D>>::New(sphere);

  auto meshInit = ls::Mesh<T>::New();
  ls::ToSurfaceMesh<T, D>(sphere, meshInit).apply();
  ls::VTKWriter<T>(meshInit, "sphereInit.vtp").apply();

  // Define constant velocity field (moving in x-direction)
  std::array<T, 3> vel = {1.0, 0.0, 0.0};
  auto velocityField = ls::SmartPointer<ConstantVelocity<T>>::New(vel);

  // Setup Advection: Forward Euler
  ls::Advect<T, D> advectFE;
  advectFE.insertNextLevelSet(sphereFE);
  advectFE.setVelocityField(velocityField);
  advectFE.setAdvectionTime(2.0);
  advectFE.setSpatialScheme(ls::SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER);
  advectFE.setTemporalScheme(ls::TemporalSchemeEnum::FORWARD_EULER);

  // Setup Advection: Runge-Kutta 2
  ls::Advect<T, D> advectRK2;
  advectRK2.insertNextLevelSet(sphereRK2);
  advectRK2.setVelocityField(velocityField);
  advectRK2.setAdvectionTime(2.0);
  advectRK2.setSpatialScheme(ls::SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER);
  advectRK2.setTemporalScheme(ls::TemporalSchemeEnum::RUNGE_KUTTA_2ND_ORDER);

  // Setup Advection: Runge-Kutta 3
  ls::Advect<T, D> advectRK3;
  advectRK3.insertNextLevelSet(sphereRK3);
  advectRK3.setVelocityField(velocityField);
  advectRK3.setAdvectionTime(2.0);
  advectRK3.setSpatialScheme(ls::SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER);
  advectRK3.setTemporalScheme(ls::TemporalSchemeEnum::RUNGE_KUTTA_3RD_ORDER);

  // Setup Advection: WENO3 Forward Euler
  ls::Advect<T, D> advectWENO3_FE;
  advectWENO3_FE.insertNextLevelSet(sphereWENO3_FE);
  advectWENO3_FE.setVelocityField(velocityField);
  advectWENO3_FE.setAdvectionTime(2.0);
  advectWENO3_FE.setSpatialScheme(ls::SpatialSchemeEnum::WENO_3RD_ORDER);
  advectWENO3_FE.setTemporalScheme(ls::TemporalSchemeEnum::FORWARD_EULER);

  // Setup Advection: WENO3 Runge-Kutta 2
  ls::Advect<T, D> advectWENO3_RK2;
  advectWENO3_RK2.insertNextLevelSet(sphereWENO3_RK2);
  advectWENO3_RK2.setVelocityField(velocityField);
  advectWENO3_RK2.setAdvectionTime(2.0);
  advectWENO3_RK2.setSpatialScheme(ls::SpatialSchemeEnum::WENO_3RD_ORDER);
  advectWENO3_RK2.setTemporalScheme(
      ls::TemporalSchemeEnum::RUNGE_KUTTA_2ND_ORDER);

  // Setup Advection: WENO3 Runge-Kutta 3
  ls::Advect<T, D> advectWENO3_RK3;
  advectWENO3_RK3.insertNextLevelSet(sphereWENO3_RK3);
  advectWENO3_RK3.setVelocityField(velocityField);
  advectWENO3_RK3.setAdvectionTime(2.0);
  advectWENO3_RK3.setSpatialScheme(ls::SpatialSchemeEnum::WENO_3RD_ORDER);
  advectWENO3_RK3.setTemporalScheme(
      ls::TemporalSchemeEnum::RUNGE_KUTTA_3RD_ORDER);

  // Setup Advection: WENO5 Forward Euler
  ls::Advect<T, D> advectWENO5_FE;
  advectWENO5_FE.insertNextLevelSet(sphereWENO5_FE);
  advectWENO5_FE.setVelocityField(velocityField);
  advectWENO5_FE.setAdvectionTime(2.0);
  advectWENO5_FE.setSpatialScheme(ls::SpatialSchemeEnum::WENO_5TH_ORDER);
  advectWENO5_FE.setTemporalScheme(ls::TemporalSchemeEnum::FORWARD_EULER);

  // Setup Advection: WENO5 Runge-Kutta 2
  ls::Advect<T, D> advectWENO5_RK2;
  advectWENO5_RK2.insertNextLevelSet(sphereWENO5_RK2);
  advectWENO5_RK2.setVelocityField(velocityField);
  advectWENO5_RK2.setAdvectionTime(2.0);
  advectWENO5_RK2.setSpatialScheme(ls::SpatialSchemeEnum::WENO_5TH_ORDER);
  advectWENO5_RK2.setTemporalScheme(
      ls::TemporalSchemeEnum::RUNGE_KUTTA_2ND_ORDER);

  // Setup Advection: WENO5 Runge-Kutta 3
  ls::Advect<T, D> advectWENO5_RK3;
  advectWENO5_RK3.insertNextLevelSet(sphereWENO5_RK3);
  advectWENO5_RK3.setVelocityField(velocityField);
  advectWENO5_RK3.setAdvectionTime(2.0);
  advectWENO5_RK3.setSpatialScheme(ls::SpatialSchemeEnum::WENO_5TH_ORDER);
  advectWENO5_RK3.setTemporalScheme(
      ls::TemporalSchemeEnum::RUNGE_KUTTA_3RD_ORDER);

  // Run Advection
  std::cout << "Running Forward Euler Advection..." << std::endl;
  viennacore::Timer timer;
  timer.start();
  advectFE.apply();
  timer.finish();
  std::cout << "Time: " << timer.currentDuration / 1e9 << "s" << std::endl;
  LSTEST_ASSERT_VALID_LS(sphereFE, T, D);

  auto meshFE = ls::Mesh<T>::New();
  ls::ToSurfaceMesh<T, D>(sphereFE, meshFE).apply();
  ls::VTKWriter<T>(meshFE, "sphereFE.vtp").apply();

  auto chamferFE = ls::CompareChamfer<T, D>(sphereRef, sphereFE);
  chamferFE.apply();
  std::cout << "Chamfer distance: " << chamferFE.getChamferDistance()
            << std::endl;
  VC_TEST_ASSERT(chamferFE.getChamferDistance() < 0.04);

  std::cout << "Running Runge-Kutta 2 Advection..." << std::endl;
  timer.start();
  advectRK2.apply();
  timer.finish();
  std::cout << "Time: " << timer.currentDuration / 1e9 << "s" << std::endl;
  LSTEST_ASSERT_VALID_LS(sphereRK2, T, D);

  auto meshRK2 = ls::Mesh<T>::New();
  ls::ToSurfaceMesh<T, D>(sphereRK2, meshRK2).apply();
  ls::VTKWriter<T>(meshRK2, "sphereRK2.vtp").apply();

  auto chamferRK2 = ls::CompareChamfer<T, D>(sphereRef, sphereRK2);
  chamferRK2.apply();
  std::cout << "Chamfer distance: " << chamferRK2.getChamferDistance()
            << std::endl;
  VC_TEST_ASSERT(chamferRK2.getChamferDistance() < 0.07);

  std::cout << "Running Runge-Kutta 3 Advection..." << std::endl;
  timer.start();
  advectRK3.apply();
  timer.finish();
  std::cout << "Time: " << timer.currentDuration / 1e9 << "s" << std::endl;
  LSTEST_ASSERT_VALID_LS(sphereRK3, T, D);

  auto meshRK3 = ls::Mesh<T>::New();
  ls::ToSurfaceMesh<T, D>(sphereRK3, meshRK3).apply();
  ls::VTKWriter<T>(meshRK3, "sphereRK3.vtp").apply();

  auto chamferRK3 = ls::CompareChamfer<T, D>(sphereRef, sphereRK3);
  chamferRK3.apply();
  std::cout << "Chamfer distance: " << chamferRK3.getChamferDistance()
            << std::endl;
  VC_TEST_ASSERT(chamferRK3.getChamferDistance() < 0.07);

  std::cout << "Running WENO3 Forward Euler Advection..." << std::endl;
  timer.start();
  advectWENO3_FE.apply();
  timer.finish();
  std::cout << "Time: " << timer.currentDuration / 1e9 << "s" << std::endl;
  LSTEST_ASSERT_VALID_LS(sphereWENO3_FE, T, D);

  auto meshWENO3_FE = ls::Mesh<T>::New();
  ls::ToSurfaceMesh<T, D>(sphereWENO3_FE, meshWENO3_FE).apply();
  ls::VTKWriter<T>(meshWENO3_FE, "sphereWENO3_FE.vtp").apply();

  auto chamferWENO3_FE = ls::CompareChamfer<T, D>(sphereRef, sphereWENO3_FE);
  chamferWENO3_FE.apply();
  std::cout << "Chamfer distance: " << chamferWENO3_FE.getChamferDistance()
            << std::endl;
  VC_TEST_ASSERT(chamferWENO3_FE.getChamferDistance() < 0.03);

  std::cout << "Running WENO3 Runge-Kutta 2 Advection..." << std::endl;
  timer.start();
  advectWENO3_RK2.apply();
  timer.finish();
  std::cout << "Time: " << timer.currentDuration / 1e9 << "s" << std::endl;
  LSTEST_ASSERT_VALID_LS(sphereWENO3_RK2, T, D);

  auto meshWENO3_RK2 = ls::Mesh<T>::New();
  ls::ToSurfaceMesh<T, D>(sphereWENO3_RK2, meshWENO3_RK2).apply();
  ls::VTKWriter<T>(meshWENO3_RK2, "sphereWENO3_RK2.vtp").apply();

  auto chamferWENO3_RK2 = ls::CompareChamfer<T, D>(sphereRef, sphereWENO3_RK2);
  chamferWENO3_RK2.apply();
  std::cout << "Chamfer distance: " << chamferWENO3_RK2.getChamferDistance()
            << std::endl;
  VC_TEST_ASSERT(chamferWENO3_RK2.getChamferDistance() < 0.008);

  std::cout << "Running WENO3 Runge-Kutta 3 Advection..." << std::endl;
  timer.start();
  advectWENO3_RK3.apply();
  timer.finish();
  std::cout << "Time: " << timer.currentDuration / 1e9 << "s" << std::endl;
  LSTEST_ASSERT_VALID_LS(sphereWENO3_RK3, T, D);

  auto meshWENO3_RK3 = ls::Mesh<T>::New();
  ls::ToSurfaceMesh<T, D>(sphereWENO3_RK3, meshWENO3_RK3).apply();
  ls::VTKWriter<T>(meshWENO3_RK3, "sphereWENO3_RK3.vtp").apply();

  auto chamferWENO3_RK3 = ls::CompareChamfer<T, D>(sphereRef, sphereWENO3_RK3);
  chamferWENO3_RK3.apply();
  std::cout << "Chamfer distance: " << chamferWENO3_RK3.getChamferDistance()
            << std::endl;
  VC_TEST_ASSERT(chamferWENO3_RK3.getChamferDistance() < 0.008);

  std::cout << "Running WENO5 Forward Euler Advection..." << std::endl;
  timer.start();
  advectWENO5_FE.apply();
  timer.finish();
  std::cout << "Time: " << timer.currentDuration / 1e9 << "s" << std::endl;
  LSTEST_ASSERT_VALID_LS(sphereWENO5_FE, T, D);

  auto meshWENO5_FE = ls::Mesh<T>::New();
  ls::ToSurfaceMesh<T, D>(sphereWENO5_FE, meshWENO5_FE).apply();
  ls::VTKWriter<T>(meshWENO5_FE, "sphereWENO5_FE.vtp").apply();

  auto chamferWENO5_FE = ls::CompareChamfer<T, D>(sphereRef, sphereWENO5_FE);
  chamferWENO5_FE.apply();
  std::cout << "Chamfer distance: " << chamferWENO5_FE.getChamferDistance()
            << std::endl;
  VC_TEST_ASSERT(chamferWENO5_FE.getChamferDistance() < 0.018);

  std::cout << "Running WENO5 Runge-Kutta 2 Advection..." << std::endl;
  timer.start();
  advectWENO5_RK2.apply();
  timer.finish();
  std::cout << "Time: " << timer.currentDuration / 1e9 << "s" << std::endl;
  LSTEST_ASSERT_VALID_LS(sphereWENO5_RK2, T, D);

  auto meshWENO5_RK2 = ls::Mesh<T>::New();
  ls::ToSurfaceMesh<T, D>(sphereWENO5_RK2, meshWENO5_RK2).apply();
  ls::VTKWriter<T>(meshWENO5_RK2, "sphereWENO5_RK2.vtp").apply();

  auto chamferWENO5_RK2 = ls::CompareChamfer<T, D>(sphereRef, sphereWENO5_RK2);
  chamferWENO5_RK2.apply();
  std::cout << "Chamfer distance: " << chamferWENO5_RK2.getChamferDistance()
            << std::endl;
  VC_TEST_ASSERT(chamferWENO5_RK2.getChamferDistance() < 0.004);

  std::cout << "Running WENO5 Runge-Kutta 3 Advection..." << std::endl;
  timer.start();
  advectWENO5_RK3.apply();
  timer.finish();
  std::cout << "Time: " << timer.currentDuration / 1e9 << "s" << std::endl;
  LSTEST_ASSERT_VALID_LS(sphereWENO5_RK3, T, D);

  auto meshWENO5_RK3 = ls::Mesh<T>::New();
  ls::ToSurfaceMesh<T, D>(sphereWENO5_RK3, meshWENO5_RK3).apply();
  ls::VTKWriter<T>(meshWENO5_RK3, "sphereWENO5_RK3.vtp").apply();

  auto chamferWENO5_RK3 = ls::CompareChamfer<T, D>(sphereRef, sphereWENO5_RK3);
  chamferWENO5_RK3.apply();
  std::cout << "Chamfer distance: " << chamferWENO5_RK3.getChamferDistance()
            << std::endl;
  VC_TEST_ASSERT(chamferWENO5_RK3.getChamferDistance() < 0.004);

  return 0;
}
