#include <iostream>
#include <numeric>

// #include <lsBooleanOperation.hpp>
#include <hrleSparseIterator.hpp>
#include <lsAdvect.hpp>
#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsTestAsserts.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

/**
  This example shows how to use lsAdvect to isotropically
  grow a 2D circle with reflective/symmetric boundary conditions.
  \example Advection2D.cpp
*/

// implement own velocity field
class velocityField : public ls::VelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int /*material*/,
                           const std::array<double, 3> & /*normalVector*/,
                           unsigned long /*pointId*/) {
    return 1.;
  }

  std::array<double, 3>
  getVectorVelocity(const std::array<double, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<double, 3> & /*normalVector*/,
                    unsigned long /*pointId*/) {
    return std::array<double, 3>({});
  }
};

int main() {

  constexpr int D = 2;

  omp_set_num_threads(4);

  double extent = 100;
  double gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i) {
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }

  auto sphere1 = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin[D] = {5., 0.};
  double radius = 7.3;

  ls::MakeGeometry<double, D>(
      sphere1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
      .apply();
  // {
  //   auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  //   ls::ToMesh<double, D>(sphere1, mesh).apply();
  //   ls::VTKWriter<double>(mesh, "sphere.vtk").apply();

  //   ls::ToSurfaceMesh<double, D>(sphere1, mesh).apply();
  //   ls::VTKWriter<double>(mesh, "before2D.vtk").apply();
  // }
  // lsExpand(sphere1, 3).apply();

  // Fill point data with original point IDs to see how LS changed
  {
    typename ls::PointData<double>::ScalarDataType pointIDs(
        sphere1->getNumberOfPoints());
    std::iota(std::begin(pointIDs), std::end(pointIDs), 0);
    sphere1->getPointData().insertNextScalarData(pointIDs, "originalIDs");
  }

  // {
  //   auto mesh = ls::SmartPointer<ls::Mesh<double>>::New();
  //   lsToMesh(sphere1, mesh, true, true).apply();
  //   ls::VTKWriter(mesh, "initial.vtp").apply();
  //   lsToSurfaceMesh(sphere1, mesh).apply();
  //   ls::VTKWriter(mesh, "surface_initial.vtp").apply();
  // }

  // Advect the sphere
  auto velocities = ls::SmartPointer<velocityField>::New();

  // std::cout << "Number of points: " <<
  // sphere1->getDomain().getNumberOfPoints()
  //           << std::endl;

  // std::cout << "Advecting" << std::endl;

  ls::Advect<double, D> advectionKernel;
  advectionKernel.insertNextLevelSet(sphere1);
  advectionKernel.setVelocityField(velocities);
  advectionKernel.setSpatialScheme(
      ls::SpatialSchemeEnum::LAX_FRIEDRICHS_1ST_ORDER);
  advectionKernel.setAdvectionTime(20.);
  advectionKernel.apply();

  LSTEST_ASSERT_VALID_LS(sphere1, double, D)

  // double advectionTime = advectionKernel.getAdvectedTime();
  // unsigned advectionSteps = advectionKernel.getNumberOfTimeSteps();
  // std::cout << "Time difference: " << advectionTime << std::endl;
  // std::cout << "Number of advection steps: " << advectionSteps << std::endl;

  // lsExpand(sphere1, 5).apply();

  // {
  //   auto mesh = ls::SmartPointer<ls::Mesh<double>>::New();
  //   lsToMesh(sphere1, mesh).apply();
  //   ls::VTKWriter(mesh, "final.vtp").apply();
  //   lsToSurfaceMesh(sphere1, mesh).apply();
  //   ls::VTKWriter(mesh, "surface_final.vtp").apply();
  // }

  // std::cout << "Pruning" << std::endl;
  // Prune<double, D>(sphere1).apply();
  // std::cout << "Expanding" << std::endl;
  // ls::Expand<double, D>(sphere1, 2).apply();

  // {
  //   auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  //   std::cout << "Extracting..." << std::endl;
  //   ls::ToSurfaceMesh<double, D>(sphere1, mesh).apply();
  //   mesh->print();
  //   ls::VTKWriter<double>(mesh, "after2D.vtk").apply();
  // }

  return 0;
}
