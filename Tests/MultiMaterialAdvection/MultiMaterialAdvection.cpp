#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Example showing how to grow/shrink different neighbouring materials
  at different speeds.
  \example MultiMaterialAdvection.cpp
*/

// implement own velocity field for advection
// in this case just grow one of the materials
class velocityField : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int material,
                           const std::array<double, 3> & /*normalVector*/) {
    // Note that only the top material grows, so having two different,
    // positive velocities will only apply in the first advection step.
    // In the next step, the levelSets of the materials will not overlap
    // anymore, so the velocity of the top layer will be used.
    // For some applications, this problem can be solved by advecting
    // the levelsets individually.
    // Grow the wrapped top material and etch the lower material.
    return ((material == 1) ? 0.5 : -0.2);
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
  omp_set_num_threads(4);

  // set up simulation domains and geometry
  double gridDelta = 0.25;

  auto sphere1 = lsSmartPointer<lsDomain<double, D>>::New(gridDelta);
  auto sphere2 = lsSmartPointer<lsDomain<double, D>>::New(gridDelta);

  double origin[3] = {5., 0., 0.};
  double radius = 9.5;

  lsMakeGeometry<double, D>(
      sphere1, lsSmartPointer<lsSphere<double, D>>::New(origin, radius))
      .apply();
  origin[0] = -5.0;
  radius = 7.3;
  lsMakeGeometry<double, D>(
      sphere2, lsSmartPointer<lsSphere<double, D>>::New(origin, radius))
      .apply();

  // Perform a boolean operation
  // Sphere2 is now the union of both original spheres.
  // This is required for the advection kernel to correctly
  // consider both materials.
  // Higher materials must always "wrap" ALL lower materials
  lsBooleanOperation<double, D>(sphere2, sphere1, lsBooleanOperationEnum::UNION)
      .apply();

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = lsSmartPointer<lsMesh>::New();
    lsToSurfaceMesh<double, D>(sphere1, mesh).apply();
    lsVTKWriter(mesh, "lower_0.vtk").apply();

    lsToSurfaceMesh<double, D>(sphere2, mesh).apply();
    lsVTKWriter(mesh, "union_0.vtk").apply();
  }

  // ADVECTION
  // fill vector with lsDomain pointers
  std::vector<lsSmartPointer<lsDomain<double, D>>> lsDomains;
  lsDomains.push_back(sphere1);
  lsDomains.push_back(sphere2);

  auto velocities = lsSmartPointer<velocityField>::New();

  std::cout << "Advecting" << std::endl;
  lsAdvect<double, D> advection(lsDomains, velocities);
  // We do not need normal vectors since our velocityField does not use them
  // This could be left on, but would decrease efficiency
  advection.setCalculateNormalVectors(false);
  advection.setAdvectionTime(5.);
  advection.apply();
  double advectionSteps = advection.getNumberOfTimeSteps();
  std::cout << "Number of Advection steps taken: " << advectionSteps
            << std::endl;

  // Output result
  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = lsSmartPointer<lsMesh>::New();
    lsToSurfaceMesh<double, D>(sphere1, mesh).apply();
    lsVTKWriter(mesh, "lower_1.vtk").apply();

    lsToSurfaceMesh<double, D>(sphere2, mesh).apply();
    mesh->print();
    lsVTKWriter(mesh, "union_1.vtk").apply();
  }

  return 0;
}
