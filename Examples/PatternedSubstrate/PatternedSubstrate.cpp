#include <iostream>
#include <random>

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
  simulation. A hexagonal pattern of rounded cones is formed.
  These cones are then used as masks for etching. A uniform
  layer is then deposited on top creating voids in the structure.
  \example PatternedSubstrate.cpp
*/

// implement velocity field describing a directional etch
class directionalEtch : public lsVelocityField<double> {
public:
  double getScalarVelocity(
      hrleVectorType<double, 3> /*coordinate*/, int material,
      hrleVectorType<double, 3> normalVector = hrleVectorType<double, 3>(0.)) {
    // etch directionally
    if (material > 0) {
      return (normalVector[2] > 0.) ? -normalVector[2] : 0;
    } else {
      return 0;
    }
  }

  hrleVectorType<double, 3> getVectorVelocity(
      hrleVectorType<double, 3> /*coordinate*/, int /*material*/,
      hrleVectorType<double,
                     3> /*normalVector = hrleVectorType<double, 3>(0.)*/) {
    return hrleVectorType<double, 3>(0.);
  }
};

// implement velocity field describing an isotropic deposition
class isotropicDepo : public lsVelocityField<double> {
public:
  double getScalarVelocity(
      hrleVectorType<double, 3> /*coordinate*/, int /*material*/,
      hrleVectorType<double,
                     3> /*normalVector = hrleVectorType<double, 3>(0.)*/) {
    // deposit isotropically everywhere
    return 1;
  }

  hrleVectorType<double, 3> getVectorVelocity(
      hrleVectorType<double, 3> /*coordinate*/, int /*material*/,
      hrleVectorType<double,
                     3> /*normalVector = hrleVectorType<double, 3>(0.)*/) {
    return hrleVectorType<double, 3>(0.);
  }
};

// create a rounded cone as the primitive pattern.
// Define a pointcloud and create a hull mesh using lsConvexHull.
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
  unsigned numberOfBasePoints = 40;
  unsigned numberOfEdgePoints = 7;
  for (unsigned i = 0; i < numberOfBasePoints; ++i) {
    double angle = double(i) / double(numberOfBasePoints) * 2. * 3.141592;
    for (unsigned j = 1; j <= numberOfEdgePoints; ++j) {
      double distance = double(j) / double(numberOfEdgePoints) * radius;
      double pointHeight = std::sqrt(double(numberOfEdgePoints - j) /
                                     double(numberOfEdgePoints)) *
                           height;
      double x = center[0] + distance * cos(angle);
      double y = center[1] + distance * sin(angle);
      cloud.insertNextPoint(
          hrleVectorType<double, 3>(x, y, center[2] + pointHeight));
    }
  }

  lsConvexHull<double, 3>(mesh, cloud).apply();
}

int main() {

  constexpr int D = 3;
  omp_set_num_threads(6);

  // scale in micrometers
  double coneDistance = 3.5;
  double xExtent = 21;
  double yConeDelta = std::sqrt(3) * coneDistance / 2;
  double yExtent = 6 * yConeDelta;

  double gridDelta = 0.15;

  double bounds[2 * D] = {-xExtent / 2., xExtent / 2., -yExtent / 2.,
                          yExtent / 2.,  -5,           5};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = lsDomain<double, D>::BoundaryType::PERIODIC_BOUNDARY;
  boundaryCons[1] = lsDomain<double, D>::BoundaryType::PERIODIC_BOUNDARY;
  boundaryCons[2] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  lsDomain<double, D> substrate(bounds, boundaryCons, gridDelta);

  double origin[D] = {0., 0., 0.001};
  double planeNormal[D] = {0., 0., 1.};

  lsMakeGeometry<double, D>(substrate, lsPlane<double, D>(origin, planeNormal))
      .apply();

  // copy the structure to add the pattern on top
  lsDomain<double, D> pattern(bounds, boundaryCons, gridDelta);
  pattern.setLevelSetWidth(2);

  // Create varying cones and put them in hexagonal pattern ---------
  {
    std::cout << "Creating pattern..." << std::endl;

    // need to place cone one grid delta below surface to avoid rounding
    hrleVectorType<double, D> coneCenter(-xExtent / 2.0 + coneDistance / 2.0,
                                         -3 * yConeDelta, -gridDelta);
    double coneRadius = 1.4;
    double coneHeight = 1.5;
    // adjust since cone is slightly below the surface
    {
      double gradient = coneHeight / coneRadius;
      coneRadius += gridDelta / gradient;
      coneHeight += gridDelta * gradient;
    }

    // random radius cones
    double variation = 0.1;
    std::mt19937 gen(532132432);
    std::uniform_real_distribution<> dis(1 - variation, 1 + variation);

    // for each row
    for (unsigned j = 0; j < 6; ++j) {
      // for each cone in a row
      for (unsigned i = 0; i < 6; ++i) {
        // make ls from cone mesh and add to substrate
        lsDomain<double, D> cone(bounds, boundaryCons, gridDelta);
        // create cone
        lsMesh coneMesh;
        makeRoundCone(coneMesh, coneCenter, coneRadius * dis(gen),
                      coneHeight * dis(gen));

        lsFromSurfaceMesh<double, D>(cone, coneMesh, false).apply();
        lsBooleanOperation<double, D> boolOp(pattern, cone,
                                             lsBooleanOperationEnum::UNION);
        boolOp.apply();

        // now shift mesh for next bool
        coneCenter[0] += coneDistance;
      }
      coneCenter[0] = -xExtent / 2. + ((j % 2) ? coneDistance / 2.0 : 0);
      coneCenter[1] += yConeDelta;
    }
  }

  lsBooleanOperation<double, D>(substrate, pattern,
                                lsBooleanOperationEnum::UNION)
      .apply();

  // Etch the substrate under the pattern ---------------------------
  unsigned numberOfEtchSteps = 30;
  std::cout << "Advecting" << std::endl;

  lsAdvect<double, D> advectionKernel;
  advectionKernel.insertNextLevelSet(pattern);
  advectionKernel.insertNextLevelSet(substrate);
  {
    directionalEtch velocities;
    advectionKernel.setVelocityField(velocities);

    // Now advect the level set, outputting every
    // advection step. Save the physical time that
    // passed during the advection.
    double passedTime = 0.;
    for (unsigned i = 0; i < numberOfEtchSteps; ++i) {
      std::cout << "\rEtch step " + std::to_string(i) + " / "
                << numberOfEtchSteps << std::flush;
      lsMesh mesh;
      lsToSurfaceMesh<double, D>(substrate, mesh).apply();
      lsVTKWriter(mesh, "substrate-" + std::to_string(i) +
                                       ".vtp").apply();

      advectionKernel.apply();
      passedTime += advectionKernel.getAdvectionTime();
    }
    std::cout << std::endl;

    lsMesh mesh;
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter(mesh, "substrate-" + std::to_string(numberOfEtchSteps) + ".vtp").apply(
        );

    std::cout << "Time passed during directional etch: " << passedTime
              << std::endl;
  }

  // Deposit new layer ----------------------------------------------
  // new level set for new layer
  lsDomain<double, D> fillLayer(substrate);
  {
    isotropicDepo velocities;
    advectionKernel.setVelocityField(velocities);

    advectionKernel.insertNextLevelSet(fillLayer);

    // stop advection in voids, which will form
    advectionKernel.setIgnoreVoids(true);

    double passedTime = 0.;
    unsigned numberOfDepoSteps = 30;
    for (unsigned i = 0; i < numberOfDepoSteps; ++i) {
      std::cout << "\rDepo step " + std::to_string(i) + " / "
                << numberOfDepoSteps << std::flush;
      lsMesh mesh;
      lsToSurfaceMesh<double, D>(fillLayer, mesh).apply();
      lsVTKWriter(mesh, "fillLayer-" + std::to_string(numberOfEtchSteps + 1 + i) + ".vtp").apply(
          );

      advectionKernel.apply();
      passedTime += advectionKernel.getAdvectionTime();
    }
    std::cout << std::endl;

    lsMesh mesh;
    lsToSurfaceMesh<double, D>(fillLayer, mesh).apply();
    lsVTKWriter(mesh, "fillLayer-" + std::to_string(numberOfEtchSteps + numberOfDepoSteps) +
    ".vtp").apply(
        );

    std::cout << "Time passed during isotropic deposition: " << passedTime
              << std::endl;
  }

  // now output the final level sets
  {
    lsMesh mesh;
    lsToSurfaceMesh<double, D>(substrate, mesh).apply();
    lsVTKWriter(mesh, "final-substrate.vtp").apply();

    lsToSurfaceMesh<double, D>(fillLayer, mesh).apply();
    lsVTKWriter(mesh, "final-fillLayer.vtp").apply();
  }

  return 0;
}
