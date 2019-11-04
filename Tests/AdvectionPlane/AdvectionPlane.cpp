#include <iostream>

// #include <lsBooleanOperation.hpp>
#include <hrleSparseIterator.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsPrune.hpp>
// #include <lsFromSurfaceMesh.hpp>
#include <lsAdvect.hpp>
#include <lsMakeGeometry.hpp>
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
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  lsDomain<double, D> plane(bounds, boundaryCons, gridDelta);

  double origin[D] = {0., 0.};
  double normal[D] = {1., 1.};

  lsMakeGeometry<double, D>(plane).makePlane(origin, normal);
  {
    lsMesh mesh;
    lsMesh explMesh;

    std::cout << "Extracting..." << std::endl;
    lsToSurfaceMesh<double, D>(plane, explMesh).apply();
    lsToMesh<double, D>(plane, mesh).apply();

    mesh.print();
    lsVTKWriter(explMesh).writeVTKLegacy("before.vtk");
    lsVTKWriter(mesh).writeVTKLegacy("beforeLS.vtk");
  }

  // fill vector with lsDomain pointers
  std::vector<lsDomain<double, D> *> lsDomains;
  lsDomains.push_back(&plane);

  velocityField velocities;

  std::cout << "number of Points: " << plane.getDomain().getNumberOfPoints()
            << std::endl;

  std::cout << "Advecting" << std::endl;
  lsAdvect<double, D> advectionKernel(lsDomains, velocities);
  advectionKernel.apply();
  double advectionTime = advectionKernel.getAdvectionTime();
  std::cout << "Time difference: " << advectionTime << std::endl;

  lsPrune<double, D>(plane).apply();
  lsExpand<double, D>(plane, 2).apply();

  std::cout << "Extracting..." << std::endl;
  lsMesh mesh;
  lsToSurfaceMesh<double, D>(plane, mesh).apply();

  // mesh.print();

  lsVTKWriter(mesh).writeVTKLegacy("after.vtk");

  return 0;
}
