#include <iostream>

// #include <lsBooleanOperation.hpp>
#include <hrleSparseIterator.hpp>
#include <lsAdvect.hpp>
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

/**
  This example shows how to use lsAdvect to isotropically
  grow a 2D circle with reflective/symmetric boundary conditions.
  \example Advection2D.cpp
*/

// implement own velocity field
class velocityField : public lsVelocityField<double> {
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

  omp_set_num_threads(2);

  double extent = 100;
  double gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i) {
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }

  auto sphere1 =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

  double origin[D] = {5., 0.};
  double radius = 7.3;

  lsMakeGeometry<double, D>(
      sphere1, lsSmartPointer<lsSphere<double, D>>::New(origin, radius))
      .apply();
  // {
  //   auto mesh = lsSmartPointer<lsMesh<>>::New();
  //   lsToMesh<double, D>(sphere1, mesh).apply();
  //   lsVTKWriter<double>(mesh, "sphere.vtk").apply();

  //   lsToSurfaceMesh<double, D>(sphere1, mesh).apply();
  //   lsVTKWriter<double>(mesh, "before2D.vtk").apply();
  // }

  // Advect the sphere
  auto velocities = lsSmartPointer<velocityField>::New();

  // std::cout << "Number of points: " <<
  // sphere1->getDomain().getNumberOfPoints()
  //           << std::endl;

  // std::cout << "Advecting" << std::endl;

  lsAdvect<double, D> advectionKernel;
  advectionKernel.insertNextLevelSet(sphere1);
  advectionKernel.setVelocityField(velocities);
  advectionKernel.setIntegrationScheme(
      lsIntegrationSchemeEnum::LAX_FRIEDRICHS_1ST_ORDER);
  advectionKernel.setAdvectionTime(20.);
  advectionKernel.apply();

  LSTEST_ASSERT_VALID_LS(sphere1, double, D)

  // double advectionTime = advectionKernel.getAdvectedTime();
  // unsigned advectionSteps = advectionKernel.getNumberOfTimeSteps();
  // std::cout << "Time difference: " << advectionTime << std::endl;
  // std::cout << "Number of advection steps: " << advectionSteps << std::endl;

  // std::cout << "Pruning" << std::endl;
  // lsPrune<double, D>(sphere1).apply();
  // std::cout << "Expanding" << std::endl;
  // lsExpand<double, D>(sphere1, 2).apply();

  // {
  //   auto mesh = lsSmartPointer<lsMesh<>>::New();
  //   std::cout << "Extracting..." << std::endl;
  //   lsToSurfaceMesh<double, D>(sphere1, mesh).apply();
  //   mesh->print();
  //   lsVTKWriter<double>(mesh, "after2D.vtk").apply();
  // }

  return 0;
}
