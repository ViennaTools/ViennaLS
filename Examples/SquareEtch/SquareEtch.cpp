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
  simulation, by creating a square trench geometry. A patch
  of exposed substrate is etched directionally.
  \example SquareEtch.cpp
*/

// implement own velocity field
class velocityField : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int material,
                           const std::array<double, 3> &normalVector) {
    // if the surface of material 1 is facing upwards, etch it anisotropically
    if (material == 1 && normalVector[1] > 0.) {
      return -std::abs(normalVector[1]);
    } else
      return 0.;
  }
};

int main() {

  constexpr int D = 2;
  omp_set_num_threads(1);

  double extent = 30;
  double gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent}; //, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[D-1] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  lsDomain<double, D> substrate(bounds, boundaryCons, gridDelta);

  double origin[3] = {0., 0., 0.};
  double planeNormal[3] = {0., 1.}; //, 1.};

  lsMakeGeometry<double, D>(substrate, lsPlane<double, D>(origin, planeNormal))
      .apply();

  lsDomain<double, D> trench(bounds, boundaryCons, gridDelta);
  // make -x and +x greater than domain for numerical stability
  // trench bottom is the initial bottom of the trench
  double trenchBottom = -2.;
  double minCorner[D] = {-extent / 1.5, /*-extent / 1.5,*/ trenchBottom};
  double maxCorner[D] = {extent / 1.5, /*extent / 1.5,*/ 1.};
  lsMakeGeometry<double, D>(trench, lsBox<double, D>(minCorner, maxCorner))
      .apply();

  // Create trench geometry
  lsBooleanOperation<double, D>(substrate, trench,
                                lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
      .apply();

  // in order only to etch the bottom of the trench, we need a mask layer
  lsDomain<double, D> mask(bounds, boundaryCons, gridDelta);
  // make downward facing plane to remove bottom of trench for the mask
  // layer
  // add small offset so bottom of trench is definetly gone
  origin[D-1] = trenchBottom + 1e-9;
  planeNormal[D-1] = -1.;
  lsMakeGeometry<double, D>(mask, lsPlane<double, D>(origin, planeNormal))
      .apply();
  lsBooleanOperation<double, D>(mask, substrate,
                                lsBooleanOperationEnum::INTERSECT)
      .apply();

  {
    std::cout << "Extracting..." << std::endl;
    // output substrate layer (which wraps around mask layer)
    // wrapping is necessary for stable advection
    lsMesh mesh;
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter(mesh, "surface-0.vtk").apply();

    // output mask layer
    lsToSurfaceMesh<double, D>(mask, mesh).apply();
    lsVTKWriter(mesh, "mask.vtk").apply();
  }

  // START ADVECTION
  velocityField velocities;

  std::cout << "Advecting" << std::endl;
  lsAdvect<double, D> advectionKernel;

  // the level set to be advected has to be inserted last
  // the other is used as the mask layer for etching
  advectionKernel.insertNextLevelSet(mask);
  advectionKernel.insertNextLevelSet(substrate);
  advectionKernel.setVelocityField(velocities);
  advectionKernel.setSaveAdvectionVelocities(true);

  // Lax Friedrichs is necessary for correct integration of the given velocity function
  advectionKernel.setIntegrationScheme(lsIntegrationSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS);

  // advect the level set 30 times
  for (unsigned counter = 1; counter < 400.; ++counter) {
    advectionKernel.apply();

    std::cout << "\rAdvection step: " << counter;
    lsMesh mesh;
    lsToMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter(mesh, "slsurface-" + std::to_string(counter) + ".vtk").apply();
    std::cout << "step size: " << advectionKernel.getAdvectionTime() << std::endl;
    // break;
  }
  std::cout << std::endl;
  // std::cout << "Number of Advection steps taken: " << advectionKernel << std::endl;

  lsMesh mesh;
  lsToSurfaceMesh<double, D>(substrate, mesh).apply();
  lsVTKWriter(mesh, "final.vtk").apply();

  return 0;
}
