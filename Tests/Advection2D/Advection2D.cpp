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

/**
  This example shows how to use lsAdvect to isotropically
  grow a 2D circle with reflective/symmetric boundary conditions.
  \example Advection2D.cpp
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

  omp_set_num_threads(2);

  double extent = 100;
  double gridDelta = 1;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  lsDomain<double, D> sphere1(bounds, boundaryCons, gridDelta);

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

  // Advect the sphere
  velocityField velocities;
  std::cout << "Number of points: " << sphere1.getDomain().getNumberOfPoints()
            << std::endl;

  std::cout << "Advecting" << std::endl;

  lsAdvect<double, D> advectionKernel;
  advectionKernel.insertNextLevelSet(sphere1);
  advectionKernel.setVelocityField(velocities);
  advectionKernel.setIntegrationScheme(
      lsIntegrationSchemeEnum::LAX_FRIEDRICHS_1ST_ORDER);
  advectionKernel.setAdvectionTime(20.);
  advectionKernel.apply();

  double advectionTime = advectionKernel.getAdvectionTime();
  unsigned advectionSteps = advectionKernel.getNumberOfTimeSteps();
  std::cout << "Time difference: " << advectionTime << std::endl;
  std::cout << "Number of advection steps: " << advectionSteps << std::endl;

  std::cout << "Pruning" << std::endl;
  lsPrune<double, D>(sphere1).apply();
  std::cout << "Expanding" << std::endl;
  lsExpand<double, D>(sphere1, 2).apply();

  {
    lsMesh mesh;
    std::cout << "Extracting..." << std::endl;
    lsToExplicitMesh<double, D>(sphere1, mesh).apply();
    mesh.print();
    lsVTKWriter(mesh).writeVTKLegacy("after2D.vtk");
  }

  return 0;
}
