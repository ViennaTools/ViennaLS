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
// This class defines the etch rate based on the surface normal.
// It models a directional etch where only upward-facing surfaces of the
// material are removed. The advection scheme will handle numerical dissipation.
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

// Analytical velocity field.
// This class provides the scalar velocity and the dissipation coefficients
// (alphas) analytically. Using analytical dissipation with the Local
// Lax-Friedrichs scheme (lsLocalLaxFriedrichsAnalytical) typically yields
// better accuracy than numerical approximations.
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

  // Toggle between numerical and analytical velocity fields.
  const bool useAnalyticalVelocity = false;

  double extent = 30;
  double gridDelta = 0.47;

  double bounds[2 * D];
  for (int i = 0; i < D; ++i) {
    bounds[2 * i] = -extent;
    bounds[2 * i + 1] = extent;
  }

  // Define boundary conditions:
  // Reflective boundaries for the sides (simulating an infinite array).
  // Infinite boundary at the top (open direction).
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i)
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[D - 1] = ls::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  // Create substrate level set (initially an empty domain with bounds).
  auto substrate = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  // Initialize the substrate as a flat surface at z=0 (or y=0 in 2D).
  double origin[3] = {0., 0., 0.};
  double planeNormal[3] = {0., D == 2, D == 3};
  {
    auto plane =
        ls::SmartPointer<ls::Plane<double, D>>::New(origin, planeNormal);
    ls::MakeGeometry<double, D>(substrate, plane).apply();
  }

  double trenchBottom = -2.;
  {
    // Create the trench geometry.
    // Define a box for the trench volume and subtract it from the substrate.
    auto trench = ls::SmartPointer<ls::Domain<double, D>>::New(
        bounds, boundaryCons, gridDelta);
    // Define the box dimensions for the trench.
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

    // Subtract the trench box from the substrate (Substrate \ Trench).
    ls::BooleanOperation<double, D>(
        substrate, trench, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  // Create a mask layer to restrict etching to the trench bottom.
  // The mask will cover the top surface and sidewalls, but expose the bottom.
  auto mask = ls::SmartPointer<ls::Domain<double, D>>::New(bounds, boundaryCons,
                                                           gridDelta);

  // Create a half-space for the mask that exists above the trench bottom.
  // Add a small offset to ensure the mask ends above the trench floor.
  origin[D - 1] = trenchBottom + 1e-9;
  planeNormal[D - 1] = -1.;
  ls::MakeGeometry<double, D>(
      mask, ls::SmartPointer<ls::Plane<double, D>>::New(origin, planeNormal))
      .apply();
  // Intersect the mask half-space with the substrate geometry.
  // Result: Mask exists where Substrate exists AND z > TrenchBottom.
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
  auto velocities = ls::SmartPointer<velocityField>::New(D);
  auto analyticalVelocities = ls::SmartPointer<analyticalField>::New(D);

  std::cout << "Advecting" << std::endl;
  ls::Advect<double, D> advectionKernel;

  // Insert level sets.
  // The order determines the material ID: Mask is Material 0, Substrate is
  // Material 1. The velocity field is configured to only move Material 1.
  advectionKernel.insertNextLevelSet(mask);
  advectionKernel.insertNextLevelSet(substrate);
  advectionKernel.setSaveAdvectionVelocities(true);

  if (useAnalyticalVelocity) {
    advectionKernel.setVelocityField(analyticalVelocities);
    // Analytical velocity fields and dissipation coefficients
    // can only be used with this spatial discretization scheme
    advectionKernel.setSpatialScheme(
        // ls::SpatialSchemeEnum::LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER);
        ls::SpatialSchemeEnum::WENO_5TH_ORDER);
  } else {
    // for numerical velocities, just use the default
    // spatial discretization scheme, which is not accurate for certain
    // velocity functions but very fast
    advectionKernel.setVelocityField(velocities);

    // For coordinate independent velocity functions
    // this numerical scheme is superior though.
    // However, it is slower.
    // advectionKernel.setSpatialScheme(
    //     ls::SpatialSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER);
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
