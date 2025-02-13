#include <iostream>
#include <numeric>

#include <lsAdvect.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsTestAsserts.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

/**
  This example shows how to use lsAdvect to create an egg
  shape from a spherical level set using directional growth rates.
  \example Advection.cpp
*/

// implement own velocity field
class velocityField : public ls::VelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int /*material*/,
                           const std::array<double, 3> &normalVector,
                           unsigned long /*pointId*/) final {
    // Some arbitrary velocity function of your liking
    // (try changing it and see what happens :)
    double velocity = 1. + ((normalVector[0] > 0) ? 2.3 : 0.5) *
                               std::abs(normalVector[0] * normalVector[0]);
    return velocity;
  }

  std::array<double, 3>
  getVectorVelocity(const std::array<double, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<double, 3> & /*normalVector*/,
                    unsigned long /*pointId*/) final {
    return std::array<double, 3>({});
  }
};

int main() {

  constexpr int D = 3;
  omp_set_num_threads(4);

  double gridDelta = 0.6999999;

  std::vector<ls::IntegrationSchemeEnum> integrationSchemes = {
      ls::IntegrationSchemeEnum::ENGQUIST_OSHER_1ST_ORDER,
      ls::IntegrationSchemeEnum::ENGQUIST_OSHER_2ND_ORDER,
      ls::IntegrationSchemeEnum::LAX_FRIEDRICHS_1ST_ORDER,
      ls::IntegrationSchemeEnum::LAX_FRIEDRICHS_2ND_ORDER,
      ls::IntegrationSchemeEnum::LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER,
      ls::IntegrationSchemeEnum::LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER,
      ls::IntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_1ST_ORDER,
      ls::IntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_2ND_ORDER,
      ls::IntegrationSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER};

  for (auto integrationScheme : integrationSchemes) {
    auto sphere1 = ls::SmartPointer<ls::Domain<double, D>>::New(gridDelta);

    double origin[3] = {5., 0., 0.};
    double radius = 7.3;

    ls::MakeGeometry<double, D>(
        sphere1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
        .apply();

    // {
    //   std::cout << "Extracting..." << std::endl;
    //   auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    //   ls::ToSurfaceMesh<double, D>(sphere1, mesh).apply();
    //   ls::VTKWriter<double>(mesh, "before.vtk").apply();
    // }

    // instantiate velocities
    auto velocities = ls::SmartPointer<velocityField>::New();

    // std::cout << "Advecting" << std::endl;

    // Fill point data with original point IDs to see how LS changed
    {
      typename ls::PointData<double>::ScalarDataType pointIDs(
          sphere1->getNumberOfPoints());
      std::iota(std::begin(pointIDs), std::end(pointIDs), 0);
      sphere1->getPointData().insertNextScalarData(pointIDs, "originalIDs");
    }

    ls::Advect<double, D> advectionKernel;
    advectionKernel.insertNextLevelSet(sphere1);
    advectionKernel.setVelocityField(velocities);
    advectionKernel.setIntegrationScheme(integrationScheme);
    // advectionKernel.setSaveAdvectionVelocities(true);

    double time = 0.;
    for (unsigned i = 0; time < 2.0 && i < 1e2; ++i) {
      advectionKernel.apply();
      time += advectionKernel.getAdvectedTime();

      // std::string fileName = std::to_string(i) + ".vtp";
      // auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
      // ls::ToMesh<double, D>(sphere1, mesh).apply();
      // ls::VTKWriter<double>(mesh, "points_" + fileName).apply();
      // ls::ToSurfaceMesh<double, D>(sphere1, mesh).apply();
      // ls::VTKWriter(mesh, "surface_" + fileName).apply();
    }

    LSTEST_ASSERT_VALID_LS(sphere1, double, D)
  }

  // std::cout << sphere1->getNumberOfPoints() << std::endl;

  // Prune<double, D>(sphere1).apply();
  // ls::Expand<double, D>(sphere1, 2).apply();

  // {
  //   std::cout << "Extracting..." << std::endl;
  //   auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  //   ls::ToSurfaceMesh<double, D>(sphere1, mesh).apply();
  //   mesh->print();
  //   ls::VTKWriter<double>(mesh, "after.vtk").apply();
  // }

  return 0;
}
