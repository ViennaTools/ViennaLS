#include <iostream>

// #include <lsBooleanOperation.hpp>
#include <hrleSparseIterator.hpp>
#include <lsAdvect.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromExplicitMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToExplicitMesh.hpp>
#include <lsToMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

// implement own velocity field
class velocityField : public lsVelocityField<double> {
public:
  double getScalarVelocity(
      hrleVectorType<double, 3> /*coordinate*/, int /*material*/,
      hrleVectorType<double,
                     3> /*normalVector = hrleVectorType<double, 3>(0.)*/) {
    return 1.;
  }

  hrleVectorType<double, 3> getVectorVelocity(
      hrleVectorType<double, 3> /*coordinate*/, int /*material*/,
      hrleVectorType<double,
                     3> /*normalVector = hrleVectorType<double, 3>(0.)*/) {
    return hrleVectorType<double, 3>(0., 0., 0.);
  }
};

int main() {

  constexpr int D = 2;

  omp_set_num_threads(2);

  double extent = 100;
  double gridDelta = 1;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  lsDomain_double_2::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = lsDomain_double_2::BoundaryType::SYMMETRIC_BOUNDARY;
  lsDomain_double_2 sphere1(bounds, boundaryCons, gridDelta);

  double origin[D] = {5., 0.};
  double radius = 7.3;

  lsMakeGeometry<double, D>(sphere1).makeSphere(origin, radius);
  {
    lsMesh mesh;
    lsToMesh<double, D>(sphere1, mesh).apply();
    lsVTKWriter(mesh).writeVTKLegacy("sphere.vtk");

    lsToExplicitMesh<double, D>(sphere1, mesh).apply();
    lsVTKWriter(mesh).writeVTKLegacy("before2D.vtk");
  }

  // fill vector with lsDomain pointers
  std::vector<lsDomain_double_2 *> lsDomains;
  lsDomains.push_back(&sphere1);

  velocityField velocities;
  std::cout << "Number of points: " << sphere1.getDomain().getNumberOfPoints()
            << std::endl;

  std::cout << "Advecting" << std::endl;
  auto advection = lsAdvect<double, D>(lsDomains, velocities);
  advection.setIntegrationScheme(4);
  double advectionTime = advection.apply(20);
  std::cout << "Time difference: " << advectionTime << std::endl;

  std::cout << "Pruning" << std::endl;
  lsPrune<double, D>(sphere1).apply();
  std::cout << "Expanding" << std::endl;
  lsExpand<double, D>(sphere1).apply(2);

  {
    lsMesh mesh;
    std::cout << "Extracting..." << std::endl;
    lsToExplicitMesh<double, D>(sphere1, mesh).apply();
    mesh.print();
    lsVTKWriter(mesh).writeVTKLegacy("after2D.vtk");
  }

  return 0;
}
