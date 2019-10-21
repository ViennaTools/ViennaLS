#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToExplicitMesh.hpp>
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
      hrleVectorType<double, 3> /*coordinate*/, int /*material*/,
      hrleVectorType<double,
                     3> /*normalVector = hrleVectorType<double, 3>(0.)*/) {
    // Some arbitrary velocity function of your liking
    // (try changing it and see what happens :)
    double velocity = 1.;
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

  constexpr int D = 3;
  omp_set_num_threads(4);

  double extent = 30;
  double gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::SYMMETRIC_BOUNDARY;
  boundaryCons[2] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  lsDomain<double, D> substrate(bounds, boundaryCons, gridDelta);

  double origin[3] = {0., 0., 0.};
  double planeNormal[3] = {0., 0., 1.};

  lsMakeGeometry<double, D>(substrate).makePlane(origin, planeNormal);

  {
    std::cout << "Extracting..." << std::endl;
    lsMesh mesh;
    lsToExplicitMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter(mesh).writeVTKLegacy("plane.vtk");
  }

  lsDomain<double, D> trench(bounds, boundaryCons, gridDelta);
  double minCorner[D] = {-extent, -extent / 4., -15.};
  double maxCorner[D] = {extent, extent / 4., 1.};
  lsMakeGeometry<double, D>(trench).makeBox(minCorner, maxCorner);

  // Create trench geometry
  lsBooleanOperation<double, D>(substrate).XOR(trench);

  {
    std::cout << "Extracting..." << std::endl;
    lsMesh mesh;
    lsToExplicitMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter(mesh).writeVTKLegacy("trench.vtk");
  }

  // Now grow new material isotropically

  // fill vector with lsDomain pointers
  std::vector<lsDomain<double, D> *> lsDomains;
  lsDomains.push_back(&substrate);

  // create new levelset for new material, which will be grown
  // since it has to wrap around the substrate, just copy it
  lsDomain<double, D> newLayer(substrate);
  lsDomains.push_back(&newLayer);

  velocityField velocities;

  std::cout << "Advecting" << std::endl;
  lsAdvect<double, D> advection(lsDomains, velocities);
  // advection.setIntegrationScheme(1);
  // advection.setCalculateNormalVectors(false);
  double advectionSteps = advection.apply(4.);
  std::cout << "Number of Advection steps taken: " << advectionSteps
            << std::endl;

  {
    std::cout << "Extracting..." << std::endl;
    for (unsigned i = 0; i < lsDomains.size(); ++i) {
      lsMesh mesh;
      lsToExplicitMesh<double, D>(*(lsDomains[i]), mesh).apply();
      lsVTKWriter(mesh).writeVTKLegacy("grown-" + std::to_string(i) + ".vtk");
      lsVTKWriter(mesh).writeVTP("grown-" + std::to_string(i) + ".vtp");
    }
  }

  return 0;
}
