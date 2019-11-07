#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsConvexHull.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  3D Example showing how to use the library for topography
  simulation. A uniform layer is deposited on top of a pillar
  using periodic boundary conditions.
  \example PeriodicBoundary.cpp
*/

// implement own velocity field
class velocityField : public lsVelocityField<double> {
public:
  double getScalarVelocity(
      hrleVectorType<double, 3> /*coordinate*/, int /*material*/,
      hrleVectorType<double,
                     3> /*normalVector = hrleVectorType<double, 3>(0.)*/) {
    // isotropic etch rate
    return 1;
  }

  hrleVectorType<double, 3> getVectorVelocity(
      hrleVectorType<double, 3> /*coordinate*/, int /*material*/,
      hrleVectorType<double,
                     3> /*normalVector = hrleVectorType<double, 3>(0.)*/) {
    return hrleVectorType<double, 3>(0.);
  }
};

// create a rounded cone as described in the paper and save it in the levelSet
// define a pointcloud and create hull mesh using lsConvexHull
void makeRoundCone(lsMesh &mesh, hrleVectorType<double, 3> center,
                   double radius, double height) {
  // cone is just a circle with a point above the center
  lsPointCloud<double, 3> cloud;
  // frist inside top point
  {
    hrleVectorType<double, 3> topPoint = center;
    topPoint[2] += height;
    cloud.insertNextPoint(topPoint);
  }

  // now create all points of the base
  unsigned numberOfBasePoints = 20;
  for (unsigned i = 0; i < numberOfBasePoints; ++i) {
    double angle = double(i) / double(numberOfBasePoints) * 2. * 3.141592;
    double x = center[0] + radius * cos(angle);
    double y = center[1] + radius * sin(angle);
    cloud.insertNextPoint(hrleVectorType<double, 3>(x, y, center[2]));
  }

  lsMesh pointMesh;
  for (unsigned i = 0; i < cloud.points.size(); ++i) {
    pointMesh.nodes.push_back(hrleVectorType<double, 3>(
        cloud.points[i][0], cloud.points[i][1], cloud.points[i][2]));
    pointMesh.vertices.push_back(hrleVectorType<unsigned, 1>(i));
  }
  lsVTKWriter(pointMesh).writeVTP("points.vtp");

  lsConvexHull<double, 3>(mesh, cloud).apply();
}

int main() {

  constexpr int D = 3;
  omp_set_num_threads(6);

  // scale in micrometers
  double coneRadius = 3.5;
  double yExtent = 5 * std::sqrt(3) * coneRadius / 4;

  double gridDelta = 0.2;

  double bounds[2 * D] = {-10.5, 10.5, -yExtent, yExtent, -5, 5};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = lsDomain<double, D>::BoundaryType::PERIODIC_BOUNDARY;
  boundaryCons[1] = lsDomain<double, D>::BoundaryType::PERIODIC_BOUNDARY;
  boundaryCons[2] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  lsDomain<double, D> substrate(bounds, boundaryCons, gridDelta);

  double origin[D] = {0., 0., 0.};
  double planeNormal[D] = {0., 0., 1.};

  lsMakeGeometry<double, D>(substrate, lsPlane<double, D>(origin, planeNormal))
      .apply();

  {
    // create spheres used for booling
    std::cout << "Creating pillar..." << std::endl;
    lsDomain<double, D> cone(bounds, boundaryCons, gridDelta);

    hrleVectorType<double, D> coneCenter(0., 0., 0.);

    lsMesh coneMesh;
    makeRoundCone(coneMesh, coneCenter, 1.5, 1.5);

    // lsMesh mesh;
    // lsToSurfaceMesh<double, D>(cone, mesh).apply();
    lsVTKWriter(coneMesh).writeVTP("pillar.vtp");
    // lsToMesh<double, D>(pillar, mesh).apply();
    // lsVTKWriter(mesh).writeVTP("LS.vtp");
    // lsBooleanOperation<double, D> boolOp(substrate, pillar,
    //                                      lsBooleanOperationEnum::UNION);
    // boolOp.apply();
    // lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    // lsVTKWriter(mesh).writeVTP("substrate.vtp");
  }

  // // Now etch the substrate isotropically
  // velocityField velocities;
  //
  // std::cout << "Advecting" << std::endl;
  //
  // lsAdvect<double, D> advectionKernel;
  // advectionKernel.insertNextLevelSet(substrate);
  // advectionKernel.setVelocityField(velocities);
  // // advectionKernel.setIntegrationScheme(
  // //     lsIntegrationSchemeEnum::ENGQUIST_OSHER_2ND_ORDER);
  //
  // // Now advect the level set 50 times, outputting every
  // // advection step. Save the physical time that
  // // passed during the advection.
  // double passedTime = 0.;
  // unsigned numberOfSteps = 50;
  // for (unsigned i = 0; i < numberOfSteps; ++i) {
  //   std::cout << "\rAdvection step " + std::to_string(i) + " / "
  //             << numberOfSteps << std::flush;
  //   lsMesh mesh;
  //   lsToSurfaceMesh<double, D>(substrate, mesh).apply();
  //   lsVTKWriter(mesh).writeVTP("pillar-" + std::to_string(i) + ".vtp");
  //
  //   advectionKernel.apply();
  //   passedTime += advectionKernel.getAdvectionTime();
  // }
  // std::cout << std::endl;
  //
  // std::cout << "Time passed during advection: " << passedTime << std::endl;

  return 0;
}
