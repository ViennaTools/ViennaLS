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

namespace ls = viennals;

using NumericType = float;

// implement own velocity field
class velocityField : public ls::VelocityField<NumericType> {
public:
  NumericType
  getScalarVelocity(const std::array<NumericType, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<NumericType, 3> & /*normalVector*/,
                    unsigned long /*pointId*/) {
    // Some arbitrary velocity function of your liking
    // (try changing it and see what happens :)
    NumericType velocity = 1.;
    return velocity;
  }

  std::array<NumericType, 3>
  getVectorVelocity(const std::array<NumericType, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<NumericType, 3> & /*normalVector*/,
                    unsigned long /*pointId*/) {
    return std::array<NumericType, 3>({}); // initialise to zero
  }
};

int main() {

  constexpr int D = 3;
  omp_set_num_threads(4);

  NumericType extent = 30;
  NumericType gridDelta = 0.5;

  double bounds[2 * D];
  for (int i = 0; i < D; ++i) {
    bounds[2 * i] = -extent;
    bounds[2 * i + 1] = extent;
  }

  ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i)
    boundaryCons[i] =
        ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[D - 1] =
      ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  NumericType origin[D] = {0.};
  NumericType planeNormal[D] = {0.};
  planeNormal[D - 1] = 1.;

  {
    auto plane =
        ls::SmartPointer<ls::Plane<NumericType, D>>::New(origin, planeNormal);
    ls::MakeGeometry<NumericType, D>(substrate, plane).apply();
  }

  {
    auto trench = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);
    // make -x and +x greater than domain for numerical stability
    NumericType ylimit = extent / 4.;
    NumericType minCorner[D];
    NumericType maxCorner[D];
    if constexpr (D == 2) {
      minCorner[0] = -ylimit;
      maxCorner[0] = ylimit;
    } else {
      minCorner[0] = -extent - 1;
      maxCorner[0] = extent + 1;
      minCorner[1] = -ylimit;
      maxCorner[1] = ylimit;
    }
    minCorner[D - 1] = -15.;
    maxCorner[D - 1] = 1.;
    auto box =
        ls::SmartPointer<ls::Box<NumericType, D>>::New(minCorner, maxCorner);
    ls::MakeGeometry<NumericType, D>(trench, box).apply();

    // Create trench geometry
    ls::BooleanOperation<NumericType, D>(
        substrate, trench, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
    ls::ToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    ls::VTKWriter<NumericType>(mesh, "trench-0.vtp").apply();
  }

  // Now grow new material isotropically

  // create new levelset for new material, which will be grown
  // since it has to wrap around the substrate, just copy it
  auto newLayer = ls::SmartPointer<ls::Domain<NumericType, D>>::New(substrate);

  auto velocities = ls::SmartPointer<velocityField>::New();

  std::cout << "Advecting" << std::endl;
  ls::Advect<NumericType, D> advectionKernel;

  // the level set to be advected has to be inserted last
  // the other could be taken as a mask layer for advection
  advectionKernel.insertNextLevelSet(substrate);
  advectionKernel.insertNextLevelSet(newLayer);
  advectionKernel.setVelocityField(velocities);
  // advectionKernel.setAdvectionTime(4.);
  unsigned counter = 1;
  advectionKernel.setIntegrationScheme(
      ls::IntegrationSchemeEnum::WENO_5TH_ORDER);
  for (NumericType time = 0; time < 4.;
       time += advectionKernel.getAdvectedTime()) {
    advectionKernel.apply();

    auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
    ls::ToSurfaceMesh<NumericType, D>(newLayer, mesh).apply();
    ls::VTKWriter<NumericType>(mesh,
                               "trench-" + std::to_string(counter) + ".vtp")
        .apply();

    ls::ToMesh<NumericType, D>(newLayer, mesh).apply();
    ls::VTKWriter<NumericType>(mesh, "LS-" + std::to_string(counter) + ".vtp")
        .apply();

    ++counter;
  }

  // NumericType advectionSteps = advectionKernel.getNumberOfTimeSteps();
  // std::cout << "Number of Advection steps taken: " << advectionSteps
  // << std::endl;

  return 0;
}
