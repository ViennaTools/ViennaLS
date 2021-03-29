#include <chrono>
#include <iostream>

#include <lsAdvect.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  This example measures the time it takes for several advection steps to run.
  \example AdvectionBenchmark.cpp
*/

// implement own velocity field
class velocityField : public lsVelocityField<double> {
  std::vector<double> &data_;

public:
  velocityField(std::vector<double> &data) : data_(data) {}

  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int /*material*/,
                           const std::array<double, 3> & /*normalVector*/,
                           unsigned long pointId) {
    // Some arbitrary velocity function of your liking
    // (try changing it and see what happens :)
    // double velocity = 1. + ((normalVector[0] > 0) ? 2.3 : 0.5) *
    //                            std::abs(normalVector[0] * normalVector[0]);
    // return velocity;
    return data_[pointId];
  }

  std::array<double, 3>
  getVectorVelocity(const std::array<double, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<double, 3> & /*normalVector*/,
                    unsigned long /*pointId*/) {
    return std::array<double, 3>({});
  }

  double
  getDissipationAlpha(int /*direction*/, int /*material*/,
                      const std::array<double, 3> & /*centralDifferences*/) {
    return 0;
  }
};

int main() {

  constexpr int D = 3;
  omp_set_num_threads(1);

  double gridDelta = 0.25;

  auto sphere1 = lsSmartPointer<lsDomain<double, D>>::New(gridDelta);

  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  lsMakeGeometry<double, D>(
      sphere1, lsSmartPointer<lsSphere<double, D>>::New(origin, radius))
      .apply();

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<double, D>(sphere1, mesh).apply();
    lsVTKWriter<double>(mesh, "before.vtk").apply();
  }

  // instantiate velocities
  std::vector<double> vels(sphere1->getNumberOfPoints() * 100, 0.31415);
  auto velocities = lsSmartPointer<velocityField>::New(vels);

  std::cout << "Advecting" << std::endl;

  const unsigned numberOfSteps = 500;
  // run several adveciton steps with different number of threads
  for (unsigned cores = 1; cores < 33; cores *= 2) {
    omp_set_num_threads(cores);

    auto levelSet = lsSmartPointer<lsDomain<double, D>>::New(gridDelta);
    levelSet->deepCopy(sphere1);

    levelSet->getDomain().segment();

    lsAdvect<double, D> advectionKernel;
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

    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<double, D>(levelSet, mesh).apply();
    lsVTKWriter<double>(mesh, "cores" + std::to_string(cores) + ".vtk").apply();
  }

  return 0;
}
