#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  3D Example showing how to use the library for topography
  simulation, by creating a trench geometry. A uniform
  layer of a different material is then grown on top.
  \example Deposition.cpp
*/

// implement own velocity field
class velocityField : public lsVelocityField<double> {
public:
  double getScalarVelocity(
      const std::array<double, 3> &/*coordinate*/, int /*material*/,
      const std::array<double,
                     3> &/*normalVector = hrleVectorType<double, 3>(0.)*/) {
    // Some arbitrary velocity function of your liking
    // (try changing it and see what happens :)
    double velocity = 1.;
    return velocity;
  }

  std::array<double, 3> getVectorVelocity(
      const std::array<double, 3> &/*coordinate*/, int /*material*/,
      const std::array<double,
                     3> &/*normalVector = hrleVectorType<double, 3>(0.)*/) {
    return std::array<double, 3>({});  // initialise to zero
  }
};

int main() {

  constexpr int D = 3;
  omp_set_num_threads(4);

  double extent = 30;
  double gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[2] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  lsDomain<double, D> substrate(bounds, boundaryCons, gridDelta);

  double origin[3] = {0., 0., 0.};
  double planeNormal[3] = {0., 0., 1.};

  lsMakeGeometry<double, D>(substrate, lsPlane<double, D>(origin, planeNormal))
      .apply();

  lsDomain<double, D> trench(bounds, boundaryCons, gridDelta);
  // make -x and +x greater than domain for numerical stability
  double minCorner[D] = {-extent - 1, -extent / 4., -15.};
  double maxCorner[D] = {extent + 1, extent / 4., 1.};
  lsMakeGeometry<double, D>(trench, lsBox<double, D>(minCorner, maxCorner))
      .apply();

  // Create trench geometry
  lsBooleanOperation<double, D>(substrate, trench,
                                lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
      .apply();

  {
    std::cout << "Extracting..." << std::endl;
    lsMesh mesh;
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter(mesh, "trench-0.vtk").apply();
  }

  // Now grow new material isotropically

  // create new levelset for new material, which will be grown
  // since it has to wrap around the substrate, just copy it
  lsDomain<double, D> newLayer(substrate);

  velocityField velocities;

  std::cout << "Advecting" << std::endl;
  lsAdvect<double, D> advectionKernel;

  // the level set to be advected has to be inserted last
  // the other could be taken as a mask layer for advection
  advectionKernel.insertNextLevelSet(substrate);
  advectionKernel.insertNextLevelSet(newLayer);

  advectionKernel.setVelocityField(velocities);
  // advectionKernel.setAdvectionTime(4.);
  unsigned counter = 1;
  for (double time = 0; time < 4.; time += advectionKernel.getAdvectionTime()) {
    advectionKernel.apply();

    lsMesh mesh;
    lsToSurfaceMesh<double, D>(newLayer, mesh).apply();
    lsVTKWriter(mesh, "trench-" + std::to_string(counter) + ".vtk").apply();

    lsToMesh<double, D>(newLayer, mesh).apply();
    lsVTKWriter(mesh, "LS-" + std::to_string(counter) + ".vtk").apply();

    ++counter;
  }

  // double advectionSteps = advectionKernel.getNumberOfTimeSteps();
  // std::cout << "Number of Advection steps taken: " << advectionSteps
  // << std::endl;

  // {
  //   std::cout << "Extracting..." << std::endl;
  //   lsMesh mesh;
  //   lsToSurfaceMesh<double, D>(newLayer, mesh).apply();
  //
  //   lsVTKWriter(mesh, "grown.vtk").apply();
  // }

  return 0;
}
