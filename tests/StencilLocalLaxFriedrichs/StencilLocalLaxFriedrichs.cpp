#include <array>
#include <cmath>
#include <iostream>
#include <map>

#include <lsAdvect.hpp>
#include <lsCompareChamfer.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsTestAsserts.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsVelocityField.hpp>

namespace ls = viennals;

// Define a constant scalar velocity field for expansion
template <class T> class ConstantScalarVelocity : public ls::VelocityField<T> {
public:
  T getScalarVelocity(const std::array<T, 3> & /*coordinate*/, int /*material*/,
                      const std::array<T, 3> & /*normalVector*/,
                      unsigned long /*pointId*/) override {
    return 1.0;
  }

  std::array<T, 3> getVectorVelocity(const std::array<T, 3> & /*coordinate*/,
                                     int /*material*/,
                                     const std::array<T, 3> & /*normalVector*/,
                                     unsigned long /*pointId*/) override {
    return {0.0, 0.0, 0.0};
  }
};

// speed 1 where the point id is that of the point at the coordinate, 0
// elsewhere: the scheme must ask the velocity of a point with its own id
template <class T, int D> class PointIdVelocity : public ls::VelocityField<T> {
public:
  std::map<std::array<long, D>, unsigned long> idOf;
  T gridDelta = 1.;
  T getScalarVelocity(const std::array<T, 3> &coordinate, int /*material*/,
                      const std::array<T, 3> & /*normalVector*/,
                      unsigned long pointId) override {
    std::array<long, D> index;
    for (int i = 0; i < D; ++i)
      index[i] = std::lround(coordinate[i] / gridDelta);
    const auto found = idOf.find(index);
    return (found != idOf.end() && found->second == pointId) ? 1. : 0.;
  }
};

int main() {
  // Define dimension and type
  constexpr int D = 3;
  using T = double;

  // Define grid and domain bounds
  double gridDelta = 0.1;
  double bounds[2 * D] = {-3.0, 3.0, -3.0, 3.0, -3.0, 3.0};
  ls::Domain<T, D>::BoundaryType boundaryCons[D];
  for (int i = 0; i < D; ++i)
    boundaryCons[i] = ls::Domain<T, D>::BoundaryType::INFINITE_BOUNDARY;

  // Create initial level set (Sphere)
  auto sphere =
      ls::SmartPointer<ls::Domain<T, D>>::New(bounds, boundaryCons, gridDelta);
  T origin[3]{};
  T radius = 1.0;
  ls::MakeGeometry<T, D>(sphere, ls::Sphere<T, D>::New(origin, radius)).apply();

  // Output initial geometry
  // auto mesh = ls::SmartPointer<ls::Mesh<T>>::New();
  // ls::ToSurfaceMesh<T, D>(sphere, mesh).apply();
  // ls::VTKWriter<T>(mesh, "sphere_initial.vtp").apply();

  // Define constant velocity field
  auto velocityField = ls::SmartPointer<ConstantScalarVelocity<T>>::New();

  // Setup Advection
  ls::Advect<T, D> advectionKernel;
  advectionKernel.insertNextLevelSet(sphere);
  advectionKernel.setVelocityField(velocityField);
  advectionKernel.setAdvectionTime(0.5);

  // Set the specific spatial discretization scheme
  advectionKernel.setSpatialScheme(
      ls::SpatialSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER);

  // Run Advection
  std::cout << "Running Stencil Local Lax Friedrichs Advection..." << std::endl;
  advectionKernel.apply();

  // Verify the result is a valid level set
  LSTEST_ASSERT_VALID_LS(sphere, T, D);

  // Check result against analytical solution
  auto sphereRef =
      ls::SmartPointer<ls::Domain<T, D>>::New(bounds, boundaryCons, gridDelta);
  ls::MakeGeometry<T, D>(sphereRef, ls::Sphere<T, D>::New(origin, radius + 0.5))
      .apply();

  auto chamfer = ls::CompareChamfer<T, D>(sphereRef, sphere);
  chamfer.apply();
  VC_TEST_ASSERT(chamfer.getChamferDistance() < 0.035);

  // Output final geometry
  // ls::ToSurfaceMesh<T, D>(sphere, mesh).apply();
  // ls::VTKWriter<T>(mesh, "sphere_final.vtp").apply();

  // the scheme asks the velocity of each point with that point's id: every
  // point of the band near the surface gets the speed 1 of its own id
  {
    auto circle = ls::SmartPointer<ls::Domain<T, D>>::New(bounds, boundaryCons,
                                                          gridDelta);
    ls::MakeGeometry<T, D>(circle, ls::Sphere<T, D>::New(origin, radius))
        .apply();
    lsInternal::StencilLocalLaxFriedrichsScalar<T, D, 1>::prepareLS(circle);
    auto field = ls::SmartPointer<PointIdVelocity<T, D>>::New();
    field->gridDelta = gridDelta;
    for (viennahrle::ConstSparseIterator<ls::Domain<T, D>::DomainType> it(
             circle->getDomain());
         !it.isFinished(); ++it)
      if (it.isDefined()) {
        std::array<long, D> index;
        for (int i = 0; i < D; ++i)
          index[i] = it.getStartIndices()[i];
        field->idOf[index] = it.getPointId();
      }
    lsInternal::StencilLocalLaxFriedrichsScalar<T, D, 1> scheme(circle, field);
    unsigned near = 0, moving = 0;
    for (viennahrle::ConstSparseIterator<ls::Domain<T, D>::DomainType> it(
             circle->getDomain());
         !it.isFinished(); ++it) {
      if (!it.isDefined() || std::abs(it.getValue()) > 0.5)
        continue;
      ++near;
      if (scheme(it.getStartIndices(), 0).first != 0.)
        ++moving;
    }
    std::cout << moving << " of " << near
              << " points near the surface take the speed of their own id"
              << std::endl;
    VC_TEST_ASSERT(near > 0 && moving == near);
  }

  std::cout << "Test passed!" << std::endl;

  return 0;
}