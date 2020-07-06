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

// Numerical velocity field.
// Advection scheme will take care of numerical
// artefacts itself.
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

// Same velocity field, but analytical
// If the dissipation alphas can be derived,
// this will produce better results than numerical
// approximations. lsLocalLaxFriedrichsAnalytical has
// to be used for advection.
class analyticalField : public lsVelocityField<double> {
  const double velocity = -1;

public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int material,
                           const std::array<double, 3> &normalVector) {
    if (material != 1)
      return 0.;

    return velocity * std::abs(normalVector[1]);
  }

  double getDissipationAlpha(int direction, int material,
                             const std::array<double, 3> &centralDifferences) {
    if (material != 1)
      return 0;

    double gradient = 0.;
    for (unsigned i = 0; i < 3; ++i) {
      gradient += centralDifferences[i] * centralDifferences[i];
    }
    gradient = std::sqrt(gradient);

    // alpha for different directions
    if (direction == 0) {
      return 0;
    } else if (direction == 1) {
      return std::abs(velocity);
    } else {
      return 0;
    }
  }
};

int main() {

  constexpr int D = 2;
  omp_set_num_threads(1);

  // Change this to use the analytical velocity field
  const bool useAnalyticalVelocity = false;

  double extent = 30;
  double gridDelta = 0.47;

  double bounds[2 * D] = {-extent, extent, -extent,
                          extent}; //, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[D - 1] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

  double origin[3] = {0., 0., 0.};
  double planeNormal[3] = {0., D == 2, D == 3};
  {
    auto plane = lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal);
    lsMakeGeometry<double, D>(substrate, plane).apply();
  }

  double trenchBottom = -2.;
  {
    auto trench = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons,
                                                           gridDelta);
    // trench bottom is the initial bottom of the trench
    double minCorner[D] = {-extent / 1.5, trenchBottom};
    double maxCorner[D] = {extent / 1.5, 1.};
    auto box = lsSmartPointer<lsBox<double, D>>::New(minCorner, maxCorner);
    lsMakeGeometry<double, D>(trench, box).apply();

    // Create trench geometry
    lsBooleanOperation<double, D>(substrate, trench,
                                  lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  // in order only to etch the bottom of the trench, we need a mask layer
  auto mask =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);
  // make downward facing plane to remove bottom of trench for the mask
  // layer
  // add small offset so bottom of trench is definetly gone
  origin[D - 1] = trenchBottom + 1e-9;
  planeNormal[D - 1] = -1.;
  lsMakeGeometry<double, D>(
      mask, lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal))
      .apply();
  lsBooleanOperation<double, D>(mask, substrate,
                                lsBooleanOperationEnum::INTERSECT)
      .apply();

  std::string fileName;
  if (useAnalyticalVelocity)
    fileName = "analytical-";
  else
    fileName = "numerical-";
  {
    std::cout << "Extracting..." << std::endl;
    // output substrate layer (which wraps around mask layer)
    // wrapping is necessary for stable advection
    auto mesh = lsSmartPointer<lsMesh>::New();
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter(mesh, fileName + "0.vtk").apply();

    // output mask layer
    lsToSurfaceMesh<double, D>(mask, mesh).apply();
    lsVTKWriter(mesh, "mask.vtk").apply();
  }

  // START ADVECTION
  auto velocities = lsSmartPointer<velocityField>::New();
  lsSmartPointer<analyticalField> analyticalVelocities;

  std::cout << "Advecting" << std::endl;
  lsAdvect<double, D> advectionKernel;

  // the level set to be advected has to be inserted last
  // the other is used as the mask layer for etching
  advectionKernel.insertNextLevelSet(mask);
  advectionKernel.insertNextLevelSet(substrate);
  advectionKernel.setSaveAdvectionVelocities(true);

  if (useAnalyticalVelocity) {
    advectionKernel.setVelocityField(analyticalVelocities);
    // Analytical velocity fields and dissipation coefficients
    // can only be used with this integration scheme
    advectionKernel.setIntegrationScheme(
        lsIntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER);
  } else {
    // for numerical velocities, just use the default
    // integration scheme, which is not accurate for certain
    // velocity functions but very fast
    advectionKernel.setVelocityField(velocities);

    // For coordinate independent velocity functions
    // this numerical scheme is superior though.
    // However, it is slower.
    // advectionKernel.setIntegrationScheme(
    //     lsIntegrationSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER);
  }

  // advect the level set until 50s have passed
  double finalTime = 50;
  unsigned counter = 1;
  for (double time = 0.; time < finalTime;
       time += advectionKernel.getAdvectedTime()) {
    advectionKernel.apply();
    std::cout << "Advection step: " << counter
              << ", time: " << advectionKernel.getAdvectedTime() << std::endl;

    auto mesh = lsSmartPointer<lsMesh>::New();
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter(mesh, fileName + std::to_string(counter) + ".vtk").apply();
    ++counter;
  }
  std::cout << std::endl;
  std::cout << "Number of Advection steps taken: " << counter << std::endl;

  auto mesh = lsSmartPointer<lsMesh>::New();
  lsToSurfaceMesh<double, D>(substrate, mesh).apply();
  lsVTKWriter(mesh, "final.vtk").apply();

  return 0;
}
