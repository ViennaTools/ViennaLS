#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
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
  simulation. A substrate with big air inclusions is etched isotropically.
  Voids are detected automatically and only etched once they are exposed.
  \example VoidEtching.cpp
*/

// implement own velocity field
class velocityField : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int /*material*/,
                           const std::array<double, 3> & /*normalVector*/,
                           unsigned long /*pointId*/) {
    // isotropic etch rate
    return 0;
  }

  std::array<double, 3>
  getVectorVelocity(const std::array<double, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<double, 3> & /*normalVector*/,
                    unsigned long /*pointId*/) {
    return std::array<double, 3>({1., 0.});
  }
};

int main() {

  constexpr int D = 2;
  omp_set_num_threads(1);

  double extent = 20;
  double gridDelta = 1;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = lsDomain<double, D>::BoundaryType::PERIODIC_BOUNDARY;
  boundaryCons[1] = lsDomain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

  double origin[D] = {0., 0.};
  double planeNormal[D] = {0., 1.};

  lsMakeGeometry<double, D>(
      substrate, lsSmartPointer<lsPlane<double, D>>::New(origin, planeNormal))
      .apply();

  std::cout << substrate->getGrid().getMinGridPoint() << std::endl;
  std::cout << substrate->getGrid().getMaxGridPoint() << std::endl;

  // for(hrleConstSparseStarIterator<lsDomain<double, D>::DomainType>
  // it(substrate.getDomain()); !it.isFinished(); it.next()) {
  //   std::cout << it.getIndices() << ": " << it.getCenter().getEndIndices() <<
  //   ", "; for(unsigned i = 0; i< 2*D; ++i) std::cout <<
  //   it.getNeighbor(i).getEndIndices() << ", "; std::cout << std::endl;
  //
  //   // << it.getNeighbor(2).getValue() << " -- " << it.getCenter().getValue()
  //   << " -- " << it.getNeighbor(0).getValue() << std::endl;
  // }

  // for(hrleConstSparseIterator<lsDomain<double, D>::DomainType>
  // it(substrate.getDomain()); !it.isFinished(); it.next()) {
  //   std::cout << it.getStartIndices() << " -- " << it.getEndIndices() <<
  //   std::endl;
  // }

  auto mesh = lsSmartPointer<lsMesh<>>::New();
  lsToMesh<double, D>(substrate, mesh).apply();
  lsVTKWriter<double>(mesh, lsFileFormatEnum::VTP, "normal.vtp").apply();

  // lsExpand<double, D>(substrate, 4).apply();
  // lsToMesh<double, D>(substrate, mesh).apply();
  // lsVTKWriter<double>(mesh, lsFileFormatEnum::VTP, "expanded.vtp").apply();

  // lsPrune<double, D>(substrate).apply();
  // lsToMesh<double, D>(substrate, mesh).apply();
  // lsVTKWriter<double>(mesh, lsFileFormatEnum::VTP, "pruned.vtp").apply();
  // -----------------------------------------------------

  {
    // create spheres used for booling
    std::cout << "Creating pillar..." << std::endl;
    auto pillar = lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons,
                                                           gridDelta);
    double lowerCorner[D] = {5, -1};
    double upperCorner[D] = {15, 10};
    lsMakeGeometry<double, D>(
        pillar, lsSmartPointer<lsBox<double, D>>::New(lowerCorner, upperCorner))
        .apply();
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<double, D>(pillar, mesh).apply();
    lsVTKWriter<double>(mesh, lsFileFormatEnum::VTP, "pillar.vtp").apply();
    lsBooleanOperation<double, D> boolOp(substrate, pillar,
                                         lsBooleanOperationEnum::UNION);
    boolOp.apply();
  }

  // Now etch the substrate isotropically
  auto velocities = lsSmartPointer<velocityField>::New();

  std::cout << "Advecting" << std::endl;

  lsAdvect<double, D> advectionKernel;
  advectionKernel.insertNextLevelSet(substrate);
  advectionKernel.setVelocityField(velocities);

  // Now advect the level set 50 times, outputting every
  // advection step. Save the physical time that
  // passed during the advection.
  double passedTime = 0.;
  unsigned numberOfSteps = 500;
  for (unsigned i = 0; i < numberOfSteps; ++i) {
    if (true) {
      std::cout << "\rAdvection step " + std::to_string(i) + " / "
                << numberOfSteps << std::flush;
      auto mesh = lsSmartPointer<lsMesh<>>::New();
      lsToSurfaceMesh<double, D>(substrate, mesh).apply();
      lsVTKWriter<double>(mesh, lsFileFormatEnum::VTP,
                          "pillar-" + std::to_string(i) + ".vtp")
          .apply();

      lsToMesh<double, D>(substrate, mesh).apply();
      lsVTKWriter<double>(mesh, lsFileFormatEnum::VTP,
                          "LS-" + std::to_string(i) + ".vtp")
          .apply();
    }

    advectionKernel.apply();
    passedTime += advectionKernel.getAdvectedTime();
  }
  std::cout << std::endl;

  std::cout << "Time passed during advection: " << passedTime << std::endl;

  return 0;
}
