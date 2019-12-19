#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  3D Example showing how to use the library for topography
  simulation. A uniform layer is deposited on top of a pillar
  using periodic boundary conditions.
  \example PeriodicBoundary.cpp
*/

// implement own velocity field
class velocityField : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int /*material*/,
                           const std::array<double, 3> & /*normalVector*/) {
    // isotropic etch rate
    return 1;
  }

  std::array<double, 3>
  getVectorVelocity(const std::array<double, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<double, 3> & /*normalVector*/) {
    return std::array<double, 3>({});
  }
};

int main() {

  constexpr int D = 3;
  omp_set_num_threads(6);

  double extent = 20;
  double gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = lsDomain<double, D>::BoundaryType::PERIODIC_BOUNDARY;
  boundaryCons[1] = lsDomain<double, D>::BoundaryType::PERIODIC_BOUNDARY;
  boundaryCons[2] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  lsDomain<double, D> substrate(bounds, boundaryCons, gridDelta);

  double origin[D] = {0., 0., 0.};
  double planeNormal[D] = {0., 0., 1.};

  lsMakeGeometry<double, D>(substrate, lsPlane<double, D>(origin, planeNormal))
      .apply();

  {
    // create spheres used for booling
    std::cout << "Creating pillar..." << std::endl;
    lsDomain<double, D> pillar(bounds, boundaryCons, gridDelta);
    double lowerCorner[D] = {15, 15, -1};
    double upperCorner[D] = {25, 25, 10};
    lsMakeGeometry<double, D>(pillar,
                              lsBox<double, D>(lowerCorner, upperCorner))
        .apply();
    lsMesh mesh;
    lsToSurfaceMesh<double, D>(pillar, mesh).apply();
    lsVTKWriter(mesh, "pillar.vtk").apply();
    lsBooleanOperation<double, D> boolOp(substrate, pillar,
                                         lsBooleanOperationEnum::UNION);
    boolOp.apply();
  }

  // Now etch the substrate isotropically
  velocityField velocities;

  std::cout << "Advecting" << std::endl;

  lsAdvect<double, D> advectionKernel;
  advectionKernel.insertNextLevelSet(substrate);
  advectionKernel.setVelocityField(velocities);
  advectionKernel.setIntegrationScheme(
      lsIntegrationSchemeEnum::ENGQUIST_OSHER_2ND_ORDER);

  // Now advect the level set 50 times, outputting every
  // advection step. Save the physical time that
  // passed during the advection.
  double passedTime = 0.;
  unsigned numberOfSteps = 50;
  for (unsigned i = 0; i < numberOfSteps; ++i) {
    std::cout << "\rAdvection step " + std::to_string(i) + " / "
              << numberOfSteps << std::flush;
    lsMesh mesh;
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter(mesh, "pillar-" + std::to_string(i) + ".vtk").apply();

    advectionKernel.apply();
    passedTime += advectionKernel.getAdvectedTime();
  }
  std::cout << std::endl;

  std::cout << "Time passed during advection: " << passedTime << std::endl;

  return 0;
}
