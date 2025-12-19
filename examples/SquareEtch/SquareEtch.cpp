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

namespace ls = viennals;

// Numerical velocity field.
// Advection scheme will take care of numerical
// artefacts itself.
class velocityField : public ls::VelocityField<double> {
  int D;

public:
  velocityField(int D) : D(D) {}

  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int material,
                           const std::array<double, 3> &normalVector,
                           unsigned long /*pointId*/) {
    // if the surface of material 1 is facing upwards, etch it anisotropically
    if (material == 1 && normalVector[D - 1] > 0.) {
      return -std::abs(normalVector[D - 1]);
    } else
      return 0.;
  }
};

// Same velocity field, but analytical
// If the dissipation alphas can be derived,
// this will produce better results than numerical
// approximations. lsLocalLaxFriedrichsAnalytical has
// to be used for advection.
class analyticalField : public ls::VelocityField<double> {
  const double velocity = -1;
  int D;

public:
  analyticalField(int D) : D(D) {}

  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int material,
                           const std::array<double, 3> &normalVector,
                           unsigned long /*pointId*/) {
    if (material != 1)
      return 0.;

    return velocity * std::abs(normalVector[D - 1]);
  }

  double getDissipationAlpha(
      int direction, int material,
      const std::array<double, 3> &centralDifferences) {
    if (material != 1)
      return 0;

    double gradient = 0.;
    for (unsigned i = 0; i < 3; ++i) {
      gradient += centralDifferences[i] * centralDifferences[i];
    }
    gradient = std::sqrt(gradient);

    // alpha for different directions
    if (direction == D - 1) {
      return std::abs(velocity);
    } else {
      return 0;
    }
  }
};

int main() {

  constexpr int D = 3;
  omp_set_num_threads(8);

  // Change this to use the analytical velocity field
  const bool useAnalyticalVelocity = false;

  double extent = 30;
  double gridDelta = 0.47;

  double bounds[2 * D] = {-extent, extent, -extent,
                          extent}; //, -extent, extent};
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i)
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[D - 1] = ls::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin[3] = {0., 0., 0.};
  double planeNormal[3] = {0., D == 2, D == 3};
  {
    auto plane =
        ls::SmartPointer<ls::Plane<double, D>>::New(origin, planeNormal);
    ls::MakeGeometry<double, D>(substrate, plane).apply();
  }

  double trenchBottom = -2.;
  {
    auto trench = ls::SmartPointer<ls::Domain<double, D>>::New(
        bounds, boundaryCons, gridDelta);
    // trench bottom is the initial bottom of the trench
    double minCorner[D];
    double maxCorner[D];
    for (int i = 0; i < D - 1; ++i) {
      minCorner[i] = -extent / 1.5;
      maxCorner[i] = extent / 1.5;
    }
    minCorner[D - 1] = trenchBottom;
    maxCorner[D - 1] = 1.;
    auto box = ls::SmartPointer<ls::Box<double, D>>::New(minCorner, maxCorner);
    ls::MakeGeometry<double, D>(trench, box).apply();

    // Create trench geometry
    ls::BooleanOperation<double, D>(
        substrate, trench, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  // in order only to etch the bottom of the trench, we need a mask layer
  auto mask = ls::SmartPointer<ls::Domain<double, D>>::New(bounds, boundaryCons,
                                                           gridDelta);
  // make downward facing plane to remove bottom of trench for the mask
  // layer
  // add small offset so bottom of trench is definetly gone
  origin[D - 1] = trenchBottom + 1e-9;
  planeNormal[D - 1] = -1.;
  ls::MakeGeometry<double, D>(
      mask, ls::SmartPointer<ls::Plane<double, D>>::New(origin, planeNormal))
      .apply();
  ls::BooleanOperation<double, D>(mask, substrate,
                                  ls::BooleanOperationEnum::INTERSECT)
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
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(substrate, mesh).apply();
    ls::VTKWriter<double>(mesh, fileName + "0.vtp").apply();

    // output mask layer
    ls::ToSurfaceMesh<double, D>(mask, mesh).apply();
    ls::VTKWriter<double>(mesh, "mask.vtp").apply();
  }

  // START ADVECTION
  velocityField velocities(D);
  analyticalField analyticalVelocities(D);

  std::cout << "Advecting" << std::endl;
  ls::Advect<double, D> advectionKernel;

  // the level set to be advected has to be inserted last
  // the other is used as the mask layer for etching
  advectionKernel.insertNextLevelSet(mask);
  advectionKernel.insertNextLevelSet(substrate);
  advectionKernel.setSaveAdvectionVelocities(true);

  if (useAnalyticalVelocity) {
    advectionKernel.setVelocityField(
        ls::SmartPointer<analyticalField>(&analyticalVelocities, [](analyticalField *) {}));
    // Analytical velocity fields and dissipation coefficients
    // can only be used with this integration scheme
    advectionKernel.setIntegrationScheme(
        // ls::IntegrationSchemeEnum::LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER);
        ls::IntegrationSchemeEnum::WENO_5TH_ORDER);
  } else {
    // for numerical velocities, just use the default
    // integration scheme, which is not accurate for certain
    // velocity functions but very fast
    advectionKernel.setVelocityField(
        ls::SmartPointer<velocityField>(&velocities, [](velocityField *) {}));
    // advectionKernel.setIntegrationScheme(
    //     ls::IntegrationSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER);
  }


  // advect the level set until 50s have passed
  double finalTime = 50;
  unsigned counter = 1;
  for (double time = 0.; time < finalTime;
       time += advectionKernel.getAdvectedTime()) {
    advectionKernel.apply();
    std::cout << "Advection step: " << counter
              << ", time: " << advectionKernel.getAdvectedTime() << std::endl;

    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(substrate, mesh).apply();
    ls::VTKWriter<double>(mesh, fileName + std::to_string(counter) + ".vtp")
        .apply();
    ++counter;
  }
  std::cout << std::endl;
  std::cout << "Number of Advection steps taken: " << counter << std::endl;

  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToSurfaceMesh<double, D>(substrate, mesh).apply();
  ls::VTKWriter<double>(mesh, "final.vtp").apply();

  return 0;
}
