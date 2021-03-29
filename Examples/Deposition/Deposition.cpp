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

using NumericType = float;

// implement own velocity field
class velocityField : public lsVelocityField<NumericType> {
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

  double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
  lsDomain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i)
    boundaryCons[i] =
        lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[2] = lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = lsSmartPointer<lsDomain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  NumericType origin[3] = {0., 0., 0.};
  NumericType planeNormal[3] = {0., 0., 1.};

  {
    auto plane =
        lsSmartPointer<lsPlane<NumericType, D>>::New(origin, planeNormal);
    lsMakeGeometry<NumericType, D>(substrate, plane).apply();
  }

  {
    auto trench = lsSmartPointer<lsDomain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);
    // make -x and +x greater than domain for numerical stability
    NumericType ylimit = extent / 4.;
    NumericType minCorner[D] = {-extent - 1, -ylimit, -15.};
    NumericType maxCorner[D] = {extent + 1, ylimit, 1.};
    auto box = lsSmartPointer<lsBox<NumericType, D>>::New(minCorner, maxCorner);
    lsMakeGeometry<NumericType, D>(trench, box).apply();

    // Create trench geometry
    lsBooleanOperation<NumericType, D>(
        substrate, trench, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();
    lsToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    lsVTKWriter<NumericType>(mesh, "trench-0.vtk").apply();
  }

  // Now grow new material isotropically

  // create new levelset for new material, which will be grown
  // since it has to wrap around the substrate, just copy it
  auto newLayer = lsSmartPointer<lsDomain<NumericType, D>>::New(substrate);

  auto velocities = lsSmartPointer<velocityField>::New();

  std::cout << "Advecting" << std::endl;
  lsAdvect<NumericType, D> advectionKernel;

  // the level set to be advected has to be inserted last
  // the other could be taken as a mask layer for advection
  advectionKernel.insertNextLevelSet(substrate);
  advectionKernel.insertNextLevelSet(newLayer);

  advectionKernel.setVelocityField(velocities);
  // advectionKernel.setAdvectionTime(4.);
  unsigned counter = 1;
  for (NumericType time = 0; time < 4.;
       time += advectionKernel.getAdvectedTime()) {
    advectionKernel.apply();

    auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();
    lsToSurfaceMesh<NumericType, D>(newLayer, mesh).apply();
    lsVTKWriter<NumericType>(mesh, "trench-" + std::to_string(counter) + ".vtk")
        .apply();

    lsToMesh<NumericType, D>(newLayer, mesh).apply();
    lsVTKWriter<NumericType>(mesh, "LS-" + std::to_string(counter) + ".vtk")
        .apply();

    ++counter;
  }

  // NumericType advectionSteps = advectionKernel.getNumberOfTimeSteps();
  // std::cout << "Number of Advection steps taken: " << advectionSteps
  // << std::endl;

  return 0;
}
