#include <iostream>

// #include <lsBooleanOperation.hpp>
#include <hrleSparseIterator.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsPrune.hpp>
// #include <lsFromExplicitMesh.hpp>
#include <lsAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToExplicitMesh.hpp>
#include <lsToMesh.hpp>
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
  omp_set_num_threads(1);

  double extent = 25;
  double gridDelta = 1;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  lsDomain_double_2::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = lsDomain_double_2::BoundaryType::SYMMETRIC_BOUNDARY;
  lsDomain_double_2 plane(bounds, boundaryCons, gridDelta);

  double origin[D] = {0., 0.};
  double normal[D] = {1., 1.};

  lsMakeGeometry<double, D>(plane).makePlane(origin, normal);
  {
    lsMesh mesh;
    lsMesh explMesh;

    std::cout << "Extracting..." << std::endl;
    lsToExplicitMesh<double, D>(plane, explMesh).apply();
    lsToMesh<double, D>(plane, mesh).apply();

    mesh.print();
    lsVTKWriter(explMesh).writeVTKLegacy("before.vtk");
    lsVTKWriter(mesh).writeVTKLegacy("beforeLS.vtk");
  }

  // fill vector with lsDomain pointers
  std::vector<lsDomain_double_2 *> lsDomains;
  lsDomains.push_back(&plane);

  velocityField velocities;

  std::cout << "number of Points: " << plane.getDomain().getNumberOfPoints()
            << std::endl;
  plane.calculateActivePointIds();
  std::cout << "number of active Points: " << plane.getNumberOfActivePoints()
            << std::endl;

  std::cout << "Advecting" << std::endl;
  double advectionTime = lsAdvect<double, D>(lsDomains, velocities).apply();
  std::cout << "Time difference: " << advectionTime << std::endl;

  lsPrune<double, D>(plane).apply();
  lsExpand<double, D>(plane).apply(2);

  std::cout << "Extracting..." << std::endl;
  lsMesh mesh;
  lsToExplicitMesh<double, D>(plane, mesh).apply();

  // mesh.print();

  lsVTKWriter(mesh).writeVTKLegacy("after.vtk");

  return 0;
}
