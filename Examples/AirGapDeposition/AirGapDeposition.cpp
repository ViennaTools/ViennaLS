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
  2D Example showing how to use the library for topography
  simulation, by creating a trench geometry. A layer of a different material is
  then grown directionally on top. \example AirGapDeposition.cpp
*/

// implement own velocity field
class velocityField : public lsVelocityField<double> {
public:
  double getScalarVelocity(
      hrleVectorType<double, 3> /*coordinate*/, int /*material*/,
      hrleVectorType<double, 3> normalVector = hrleVectorType<double, 3>(0.)) {
    // velocity is proportional to the normal vector
    double velocity = std::abs(normalVector[0]) + std::abs(normalVector[1]);
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

  constexpr int D = 2;
  omp_set_num_threads(2);

  double extent = 30;
  double gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[1] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  lsDomain<double, D> substrate(bounds, boundaryCons, gridDelta);

  double origin[2] = {0., 0.};
  double planeNormal[2] = {0., 1.};

  lsMakeGeometry<double, D>(substrate).makePlane(origin, planeNormal);

  {
    std::cout << "Extracting..." << std::endl;
    lsMesh mesh;
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter(mesh).writeVTP("plane.vtp");
  }

  {
    // create layer used for booling
    std::cout << "Creating box..." << std::endl;
    lsDomain<double, D> trench(bounds, boundaryCons, gridDelta);
    double minCorner[D] = {-extent / 6., -25.};
    double maxCorner[D] = {extent / 6., 1.};
    lsMakeGeometry<double, D>(trench).makeBox(minCorner, maxCorner);

    {
      std::cout << "Extracting..." << std::endl;
      lsMesh mesh;
      lsToMesh<double, D>(trench, mesh).apply();
      lsVTKWriter(mesh).writeVTP("box.vtp");
    }

    // Create trench geometry
    std::cout << "Booling trench..." << std::endl;
    lsBooleanOperation<double, D>(substrate, trench,
                                  lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  // Now grow new material

  // create new levelset for new material, which will be grown
  // since it has to wrap around the substrate, just copy it
  std::cout << "Creating new layer..." << std::endl;
  lsDomain<double, D> newLayer(substrate);

  velocityField velocities;

  std::cout << "Advecting" << std::endl;
  lsAdvect<double, D> advectionKernel;

  // the level set to be advected has to be inserted last
  // the other could be taken as a mask layer for advection
  advectionKernel.insertNextLevelSet(substrate);
  advectionKernel.insertNextLevelSet(newLayer);

  advectionKernel.setVelocityField(velocities);
  advectionKernel.setIgnoreVoids(true);

  // Now advect the level set 50 times, outputting every
  // 10th advection step. Save the physical time that
  // passed during the advection.
  double passedTime = 0.;
  unsigned numberOfSteps = 60;
  for (unsigned i = 0; i < numberOfSteps; ++i) {
    advectionKernel.apply();
    passedTime += advectionKernel.getAdvectionTime();

    std::cout << "\rAdvection step " + std::to_string(i) + " / "
              << numberOfSteps << std::flush;
    lsMesh mesh;
    lsToSurfaceMesh<double, D>(newLayer, mesh).apply();
    lsVTKWriter(mesh).writeVTKLegacy("trench" + std::to_string(i) + ".vtk");
  }
  std::cout << std::endl;
  std::cout << "Time passed during advection: " << passedTime << std::endl;

  return 0;
}
