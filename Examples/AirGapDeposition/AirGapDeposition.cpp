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

using NumericType = float;

// implement own velocity field
class velocityField : public lsVelocityField<NumericType> {
public:
  NumericType
  getScalarVelocity(const std::array<NumericType, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<NumericType, 3> &normalVector,
                    unsigned long /*pointId*/) {
    // velocity is proportional to the normal vector
    NumericType velocity =
        std::abs(normalVector[0]) + std::abs(normalVector[1]);
    return velocity;
  }

  std::array<NumericType, 3>
  getVectorVelocity(const std::array<NumericType, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<NumericType, 3> & /*normalVector*/,
                    unsigned long /*pointId*/) {
    return std::array<NumericType, 3>({});
  }
};

int main() {

  constexpr int D = 2;
  omp_set_num_threads(2);

  NumericType extent = 30;
  NumericType gridDelta = 0.5;

  hrleCoordType bounds[2 * D] = {-extent, extent, -extent, extent};
  lsDomain<NumericType, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[1] = lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = lsSmartPointer<lsDomain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  NumericType origin[2] = {0., 0.};
  NumericType planeNormal[2] = {0., 1.};

  {
    auto plane =
        lsSmartPointer<lsPlane<NumericType, D>>::New(origin, planeNormal);
    lsMakeGeometry<NumericType, D>(substrate, plane).apply();
  }

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();
    lsToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    lsVTKWriter<NumericType>(mesh, "plane.vtk").apply();
  }

  {
    // create layer used for booling
    std::cout << "Creating box..." << std::endl;
    auto trench = lsSmartPointer<lsDomain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);
    NumericType xlimit = extent / 6.;
    NumericType minCorner[D] = {-xlimit, -25.};
    NumericType maxCorner[D] = {xlimit, 1.};
    auto box = lsSmartPointer<lsBox<NumericType, D>>::New(minCorner, maxCorner);
    lsMakeGeometry<NumericType, D>(trench, box).apply();

    {
      std::cout << "Extracting..." << std::endl;
      auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();
      lsToMesh<NumericType, D>(trench, mesh).apply();
      lsVTKWriter<NumericType>(mesh, "box.vtk").apply();
    }

    // Create trench geometry
    std::cout << "Booling trench..." << std::endl;
    lsBooleanOperation<NumericType, D>(
        substrate, trench, lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  // Now grow new material

  // create new levelset for new material, which will be grown
  // since it has to wrap around the substrate, just copy it
  std::cout << "Creating new layer..." << std::endl;
  auto newLayer = lsSmartPointer<lsDomain<NumericType, D>>::New(substrate);

  auto velocities = lsSmartPointer<velocityField>::New();

  std::cout << "Advecting" << std::endl;
  lsAdvect<NumericType, D> advectionKernel;

  // the level set to be advected has to be inserted last
  // the other could be taken as a mask layer for advection
  advectionKernel.insertNextLevelSet(substrate);
  advectionKernel.insertNextLevelSet(newLayer);

  advectionKernel.setVelocityField(velocities);
  advectionKernel.setIgnoreVoids(true);

  // Now advect the level set 50 times, outputting every
  // advection step. Save the physical time that
  // passed during the advection.
  NumericType passedTime = 0.;
  unsigned numberOfSteps = 60;
  for (unsigned i = 0; i < numberOfSteps; ++i) {
    advectionKernel.apply();
    passedTime += advectionKernel.getAdvectedTime();

    std::cout << "\rAdvection step " + std::to_string(i) + " / "
              << numberOfSteps << std::flush;
    auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();
    lsToSurfaceMesh<NumericType, D>(newLayer, mesh).apply();
    lsVTKWriter<NumericType>(mesh, "trench" + std::to_string(i) + ".vtk")
        .apply();
  }
  std::cout << std::endl;
  std::cout << "Time passed during advection: " << passedTime << std::endl;

  return 0;
}
