#include <iostream>
#include <array>

#include <lsAdvectForwardEuler.hpp>
#include <lsAdvectRungeKutta3.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsTestAsserts.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsVelocityField.hpp>

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

  // Define grid and domain bounds
  double gridDelta = 0.1;
  double bounds[2 * D] = {-5.0, 5.0, -5.0, 5.0, -5.0, 5.0};
  ls::Domain<T, D>::BoundaryType boundaryCons[D];
  for (int i = 0; i < D; ++i)
    boundaryCons[i] = ls::Domain<T, D>::BoundaryType::INFINITE_BOUNDARY;

  // Create initial level set (Sphere)
  auto sphere = ls::SmartPointer<ls::Domain<T, D>>::New(bounds, boundaryCons, gridDelta);
  T origin[3] = {0.0, 0.0, 0.0};
  T radius = 1.5;
  ls::MakeGeometry<T, D>(sphere, ls::Sphere<T, D>::New(origin, radius)).apply();

  // Create copies for Forward Euler and RK3
  auto sphereFE = ls::SmartPointer<ls::Domain<T, D>>::New(sphere);
  auto sphereRK3 = ls::SmartPointer<ls::Domain<T, D>>::New(sphere);

  // Define constant velocity field (moving in x-direction)
  std::array<T, 3> vel = {1.0, 0.0, 0.0};
  auto velocityField = ls::SmartPointer<ConstantVelocity<T>>::New(vel);

  // Setup Advection: Forward Euler
  ls::AdvectForwardEuler<T, D> advectFE;
  advectFE.insertNextLevelSet(sphereFE);
  advectFE.setVelocityField(velocityField);
  advectFE.setAdvectionTime(2.0);
  advectFE.setIntegrationScheme(ls::IntegrationSchemeEnum::ENGQUIST_OSHER_1ST_ORDER);

  // Setup Advection: Runge-Kutta 3
  ls::AdvectRungeKutta3<T, D> advectRK3;
  advectRK3.insertNextLevelSet(sphereRK3);
  advectRK3.setVelocityField(velocityField);
  advectRK3.setAdvectionTime(2.0);
  advectRK3.setIntegrationScheme(ls::IntegrationSchemeEnum::ENGQUIST_OSHER_1ST_ORDER);

  // Run Advection
  std::cout << "Running Forward Euler Advection..." << std::endl;
  advectFE.apply();
  LSTEST_ASSERT_VALID_LS(sphereFE, T, D);

  auto meshFE = ls::Mesh<T>::New();
  ls::ToSurfaceMesh<T, D>(sphereFE, meshFE).apply();
  ls::VTKWriter<T>(meshFE, "sphereFE.vtp").apply();

  std::cout << "Running Runge-Kutta 3 Advection..." << std::endl;
  advectRK3.apply();
  LSTEST_ASSERT_VALID_LS(sphereRK3, T, D);

  auto meshRK3 = ls::Mesh<T>::New();
  ls::ToSurfaceMesh<T, D>(sphereRK3, meshRK3).apply();
  ls::VTKWriter<T>(meshRK3, "sphereRK3.vtp").apply();

  return 0;
}
