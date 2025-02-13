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

namespace ls = viennals;

/**
  Minimal example of a plane surface being moved
  using lsAdvect
  \example AdvectionPlane.cpp
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
  omp_set_num_threads(1);

  double extent = 25;
  double gridDelta = 1;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  ls::BoundaryConditionEnum boundaryCons[D];
  boundaryCons[0] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
  boundaryCons[1] = ls::BoundaryConditionEnum::INFINITE_BOUNDARY;

  auto plane = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin[D] = {0., 0.};
  double normal[D] = {2., 1.};

  ls::MakeGeometry<double, D>(
      plane, ls::SmartPointer<ls::Plane<double, D>>::New(origin, normal))
      .apply();
  // {
  //   auto mesh = ls::ls::SmartPointer<ls::Mesh<>>::New();

  //   std::cout << "Extracting..." << std::endl;
  //   ls::ToSurfaceMesh<double, D>(plane, mesh).apply();
  //   ls::VTKWriter<double>(mesh, "before.vtk").apply();

  //   ls::ToMesh<double, D>(plane, mesh).apply();
  //   ls::VTKWriter<double>(mesh, "beforeLS.vtk").apply();

  //   mesh->print();
  // }

  auto velocities = ls::SmartPointer<velocityField>::New();

  // std::cout << "number of Points: " << plane->getNumberOfPoints() <<
  // std::endl;

  // std::cout << "Advecting" << std::endl;
  ls::Advect<double, D> advectionKernel(plane, velocities);
  advectionKernel.apply();
  double advectionTime = advectionKernel.getAdvectedTime();
  // std::cout << "Time difference: " << advectionTime << std::endl;

  LSTEST_ASSERT_VALID_LS(plane, double, D)

  // Prune<double, D>(plane).apply();
  // ls::Expand<double, D>(plane, 2).apply();

  // std::cout << "Extracting..." << std::endl;
  // auto mesh = ls::ls::SmartPointer<ls::Mesh<>>::New();
  // ls::ToSurfaceMesh<double, D>(plane, mesh).apply();

  // // mesh.print();

  // ls::VTKWriter<double>(mesh, "after.vtk").apply();

  return 0;
}
