#include <chrono>
#include <iostream>

#include <lsAdvect.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

/**
  This example measures the time it takes for several advection steps to run.
  \example AdvectionBenchmark.cpp
*/

// implement own velocity field
class velocityField : public ls::VelocityField<double> {
  std::vector<double> &data_;

public:
  velocityField(std::vector<double> &data) : data_(data) {}

  double getScalarVelocity(const ls::Vec3D<double> & /*coordinate*/,
                           int /*material*/,
                           const ls::Vec3D<double> & /*normalVector*/,
                           unsigned long pointId) override {
    // Some arbitrary velocity function of your liking
    // (try changing it and see what happens :)
    // double velocity = 1. + ((normalVector[0] > 0) ? 2.3 : 0.5) *
    //                            std::abs(normalVector[0] * normalVector[0]);
    // return velocity;
    return data_[pointId];
  }
};

int main() {

  constexpr int D = 3;
  omp_set_num_threads(1);

  double gridDelta = 0.25;

  auto sphere1 = ls::SmartPointer<ls::Domain<double, D>>::New(gridDelta);

  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  ls::MakeGeometry<double, D>(
      sphere1, ls::SmartPointer<ls::Sphere<double, D>>::New(origin, radius))
      .apply();

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(sphere1, mesh).apply();
    ls::VTKWriter<double>(mesh, "before.vtk").apply();
  }

  // instantiate velocities
  std::vector<double> vels(sphere1->getNumberOfPoints() * 100, 0.31415);
  auto velocities = ls::SmartPointer<velocityField>::New(vels);

  std::cout << "Advecting" << std::endl;

  const unsigned numberOfSteps = 500;
  // run several adveciton steps with different number of threads
  for (unsigned cores = 1; cores < 17; cores *= 2) {
    omp_set_num_threads(cores);

    auto levelSet = ls::SmartPointer<ls::Domain<double, D>>::New(gridDelta);
    levelSet->deepCopy(sphere1);

    levelSet->getDomain().segment();

    ls::Advect<double, D> advectionKernel;
    advectionKernel.insertNextLevelSet(levelSet);
    advectionKernel.setVelocityField(velocities);

    const auto start = std::chrono::high_resolution_clock::now();
    for (unsigned i = 0; i < numberOfSteps; ++i) {
      advectionKernel.apply();
    }
    const auto stop = std::chrono::high_resolution_clock::now();
    std::cout << "Advection with " << cores << ": "
              << std::chrono::duration_cast<std::chrono::milliseconds>(stop -
                                                                       start)
                     .count()
              << "\n";

    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(levelSet, mesh).apply();
    ls::VTKWriter<double>(mesh, "cores" + std::to_string(cores) + ".vtk")
        .apply();
  }

  return 0;
}
