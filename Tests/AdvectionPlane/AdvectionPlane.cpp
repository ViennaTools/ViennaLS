#include <iostream>

// #include <lsBooleanOperation.hpp>
#include <hrleSparseIterator.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsPrune.hpp>
// #include <lsFromSurfaceMesh.hpp>
#include <lsAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsTestAsserts.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Minimal example of a plane surface being moved
  using lsAdvect
  \example AdvectionPlane.cpp
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
  omp_set_num_threads(1);

  double extent = 25;
  double gridDelta = 1;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[1] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto plane =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

  double origin[D] = {0., 0.};
  double normal[D] = {2., 1.};

  lsMakeGeometry<double, D>(
      plane, lsSmartPointer<lsPlane<double, D>>::New(origin, normal))
      .apply();
  // {
  //   auto mesh = lsSmartPointer<lsMesh<>>::New();

  //   std::cout << "Extracting..." << std::endl;
  //   lsToSurfaceMesh<double, D>(plane, mesh).apply();
  //   lsVTKWriter<double>(mesh, "before.vtk").apply();

  //   lsToMesh<double, D>(plane, mesh).apply();
  //   lsVTKWriter<double>(mesh, "beforeLS.vtk").apply();

  //   mesh->print();
  // }

  auto velocities = lsSmartPointer<velocityField>::New();

  // std::cout << "number of Points: " << plane->getNumberOfPoints() <<
  // std::endl;

  // std::cout << "Advecting" << std::endl;
  lsAdvect<double, D> advectionKernel(plane, velocities);
  advectionKernel.apply();
  double advectionTime = advectionKernel.getAdvectedTime();
  // std::cout << "Time difference: " << advectionTime << std::endl;

  LSTEST_ASSERT_VALID_LS(plane, double, D)

  // lsPrune<double, D>(plane).apply();
  // lsExpand<double, D>(plane, 2).apply();

  // std::cout << "Extracting..." << std::endl;
  // auto mesh = lsSmartPointer<lsMesh<>>::New();
  // lsToSurfaceMesh<double, D>(plane, mesh).apply();

  // // mesh.print();

  // lsVTKWriter<double>(mesh, "after.vtk").apply();

  return 0;
}
