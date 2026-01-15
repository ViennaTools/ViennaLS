#include <iostream>
#include <numeric>

#include <lsAdvect.hpp>
#include <lsCompareCriticalDimensions.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsTestAsserts.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsVelocityField.hpp>

namespace ls = viennals;

/**
  This example shows how to use lsAdvect to create an egg
  shape from a spherical level set using directional growth rates.
  \example Advection.cpp
*/

// implement own velocity field
class velocityField : public ls::VelocityField<double> {
public:
  double getScalarVelocity(const ls::Vec3D<double> & /*coordinate*/,
                           int /*material*/,
                           const ls::Vec3D<double> &normalVector,
                           unsigned long /*pointId*/) final {
    // Some arbitrary velocity function of your liking
    // (try changing it and see what happens :)
    double velocity = 1. + ((normalVector[0] > 0) ? 2.3 : 0.5) *
                               std::abs(normalVector[0] * normalVector[0]);
    return velocity;
  }

  ls::Vec3D<double>
  getVectorVelocity(const ls::Vec3D<double> & /*coordinate*/, int /*material*/,
                    const ls::Vec3D<double> & /*normalVector*/,
                    unsigned long /*pointId*/) final {
    return {0., 0., 0.};
  }
};

int main() {

  constexpr int D = 3;
  omp_set_num_threads(8);

  double gridDelta = 0.6999999;

  std::vector<ls::SpatialSchemeEnum> spatialSchemes = {
      ls::SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER,
      ls::SpatialSchemeEnum::ENGQUIST_OSHER_2ND_ORDER,
      ls::SpatialSchemeEnum::LAX_FRIEDRICHS_1ST_ORDER,
      ls::SpatialSchemeEnum::LAX_FRIEDRICHS_2ND_ORDER,
      ls::SpatialSchemeEnum::LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER,
      ls::SpatialSchemeEnum::LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER,
      ls::SpatialSchemeEnum::LOCAL_LAX_FRIEDRICHS_1ST_ORDER,
      ls::SpatialSchemeEnum::LOCAL_LAX_FRIEDRICHS_2ND_ORDER,
      ls::SpatialSchemeEnum::WENO_3RD_ORDER,
      ls::SpatialSchemeEnum::WENO_5TH_ORDER};

  for (auto scheme : spatialSchemes) {
    auto sphere1 = ls::Domain<double, D>::New(gridDelta);

    double origin[3] = {5., 0., 0.};
    double radius = 7.3;

    ls::MakeGeometry<double, D>(sphere1,
                                ls::Sphere<double, D>::New(origin, radius))
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
    advectionKernel.setSpatialScheme(scheme);
    advectionKernel.setSaveAdvectionVelocities(true);

    double time = 0.;
    for (unsigned i = 0; time < 1.0 && i < 50; ++i) {
      advectionKernel.apply();
      time += advectionKernel.getAdvectedTime();

      // std::string fileName = std::to_string(i) + ".vtp";
      // auto mesh = ls::Mesh<>::New();
      // ls::ToMesh<double, D>(sphere1, mesh).apply();
      // ls::VTKWriter<double>(mesh, "points_" + fileName).apply();
      // ls::ToSurfaceMesh<double, D>(sphere1, mesh).apply();
      // ls::VTKWriter(mesh, "surface_" + fileName).apply();
    }

    LSTEST_ASSERT_VALID_LS(sphere1, double, D)

    // Check critical dimensions against initial state
    auto sphereRef = ls::Domain<double, D>::New(gridDelta);
    ls::MakeGeometry<double, D>(sphereRef,
                                ls::Sphere<double, D>::New(origin, radius))
        .apply();

    ls::CompareCriticalDimensions<double, D> compare(sphereRef, sphere1);
    std::array<double, D> minB, maxB;
    minB.fill(std::numeric_limits<double>::lowest());
    maxB.fill(std::numeric_limits<double>::max());

    // Measure Max X (Growth velocity ~3.3)
    compare.addRange(0, minB, maxB, true);
    // Measure Min X (Growth velocity ~1.5)
    compare.addRange(0, minB, maxB, false);
    // Measure Max Y (Growth velocity ~1.0)
    compare.addRange(1, minB, maxB, true);

    compare.apply();

    double posRef, posSample, diffx, diffy, diffz;

    compare.getCriticalDimensionResult(0, posRef, posSample, diffx);
    compare.getCriticalDimensionResult(1, posRef, posSample, diffy);
    compare.getCriticalDimensionResult(2, posRef, posSample, diffz);
    VC_TEST_ASSERT(diffx > 3.45 && diffx < 3.55 && diffy > 1.45 &&
                   diffy < 1.65 && diffz > 0.85 && diffz < 1.10);
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
