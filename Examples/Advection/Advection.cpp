#include <iostream>

#include <lsAdvect.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToExplicitMesh.hpp>
#include <lsVTKWriter.hpp>

// implement own velocity field
class velocityField : public lsVelocityField<double> {
public:
  double getScalarVelocity(
      hrleVectorType<double, 3> /*coordinate*/, int /*material*/,
      hrleVectorType<double, 3> normalVector = hrleVectorType<double, 3>(0.)) {
    // Some arbitrary velocity function of your liking
    // (try changing it and see what happens :)
    double velocity = 1. + ((normalVector[0] > 0) ? 2.3 : 0.5) *
                               std::abs(normalVector[0] * normalVector[0]);
    return velocity;
  }

  hrleVectorType<double, 3> getVectorVelocity(
      hrleVectorType<double, 3> /*coordinate*/, int /*material*/,
      hrleVectorType<double,
                     3> /*normalVector = hrleVectorType<double, 3>(0.)*/) {
    return hrleVectorType<double, 3>(0.);
  }
};

int main() {

  constexpr int D = 3;
  omp_set_num_threads(4);

  double extent = 15;
  double gridDelta = 0.25;

  // including lsDomain.hpp provides typedefs for pre-built
  // template specialisations, such as lsDomain<double, 3>
  lsDomain_double_3 sphere1(gridDelta);

  double origin[3] = {5., 0., 0.};
  double radius = 7.3;

  lsMakeGeometry<double, D>(sphere1).makeSphere(origin, radius);

  {
    std::cout << "Extracting..." << std::endl;
    lsMesh mesh;
    lsToExplicitMesh<double, D>(sphere1, mesh).apply();
    lsVTKWriter(mesh).writeVTKLegacy("before.vtk");
  }

  // fill vector with lsDomain pointers
  std::vector<lsDomain_double_3 *> lsDomains;
  lsDomains.push_back(&sphere1);

  velocityField velocities;

  std::cout << "Advecting" << std::endl;
  lsAdvect<double, D> advection(lsDomains, velocities);
  // advection.setIntegrationScheme(1);
  // advection.setCalculateNormalVectors(false);
  double advectionSteps = advection.apply(2.);
  std::cout << "Number of Advection steps taken: " << advectionSteps
            << std::endl;

  lsPrune<double, D>(sphere1).apply();
  lsExpand<double, D>(sphere1).apply(2);

  {
    std::cout << "Extracting..." << std::endl;
    lsMesh mesh;
    lsToExplicitMesh<double, D>(sphere1, mesh).apply();
    mesh.print();
    lsVTKWriter(mesh).writeVTKLegacy("after.vtk");
  }

  return 0;
}
