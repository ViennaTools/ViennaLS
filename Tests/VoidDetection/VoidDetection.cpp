#include <algorithm>
#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMarkVoidPoints.hpp>
#include <lsToExplicitMesh.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  Example showing how to use void detection.
  \example VoidDetection.cpp
*/

// implement own velocity field
class velocityField : public lsVelocityField<double> {
public:
  double getScalarVelocity(
      hrleVectorType<double, 3> /*coordinate*/, int /*material*/,
      hrleVectorType<double, 3> /*normalVector = hrleVectorType<double, 3>(0.)*/) {
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
  constexpr int D = 2;
  omp_set_num_threads(1);

  double extent = 10;
  double gridDelta = 1;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i)
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::SYMMETRIC_BOUNDARY;

  boundaryCons[1] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;
  lsDomain<double, D> substrate(bounds, boundaryCons, gridDelta);

  double origin[D] = {0., 0.};
  double normal[D] = {0., 1.};

  lsMakeGeometry<double, D>(substrate).makePlane(origin, normal);
  {
    lsDomain<double, D> hole(bounds, boundaryCons, gridDelta);
    origin[1] = -5.;
    lsMakeGeometry<double, D>(hole).makeSphere(origin, 3.);

    lsBooleanOperation<double, D>(substrate, hole,
                                  lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  {
    lsMesh explMesh;

    std::cout << "Extracting..." << std::endl;
    lsToExplicitMesh<double, D>(substrate, explMesh).apply();

    lsVTKWriter(explMesh).writeVTKLegacy("before.vtk");
  }

  lsMarkVoidPoints<double, D>(substrate).apply();

  {
    std::cout << "Extracting..." << std::endl;
    lsMesh mesh;
    lsToMesh<double, D>(substrate, mesh).apply();

    auto voidPointMarkers = substrate.getVoidPointMarkers();
    std::vector<double> isVoid(voidPointMarkers.size()); // 0 = not void, 1 = void
    for (unsigned i = 0; i < isVoid.size(); ++i) {
      isVoid[i] = (voidPointMarkers[i]) ? 1. : 0.;
    }

    std::cout << "Points: " << substrate.getNumberOfPoints() << std::endl;
    std::cout << "Markers: " << isVoid.size() << std::endl;

    mesh.insertNextScalarData(isVoid, "voidMarkers");

    lsVTKWriter(mesh).writeVTKLegacy("after.vtk");
  }


  // Advection
  velocityField velocities;
  lsAdvect<double, D> advectionKernel(substrate, velocities);
  advectionKernel.setIgnoreVoids(true);
  for(unsigned i=0; i< 30; ++i) {
    {
      lsMesh mesh;
      lsToExplicitMesh<double, D>(substrate, mesh).apply();
      lsVTKWriter(mesh).writeVTKLegacy("out-" + std::to_string(i) + ".vtk");
    }
    advectionKernel.apply();
  }

  return 0;
}
