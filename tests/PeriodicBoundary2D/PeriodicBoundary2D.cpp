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

namespace ls = viennals;

// implement own velocity field
class velocityField : public ls::VelocityField<double> {
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
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = ls::Domain<double, D>::BoundaryType::PERIODIC_BOUNDARY;
  boundaryCons[1] = ls::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  double origin[D] = {0., 0.};
  double planeNormal[D] = {0., 1.};

  ls::MakeGeometry<double, D>(
      substrate,
      ls::SmartPointer<ls::Plane<double, D>>::New(origin, planeNormal))
      .apply();

  std::cout << substrate->getGrid().getMinGridPoint() << std::endl;
  std::cout << substrate->getGrid().getMaxGridPoint() << std::endl;

  // for(hrleConstSparseStarIterator<Domain<double, D>::DomainType, 1>
  // it(substrate.getDomain()); !it.isFinished(); it.next()) {
  //   std::cout << it.getIndices() << ": " << it.getCenter().getEndIndices() <<
  //   ", "; for(unsigned i = 0; i< 2*D; ++i) std::cout <<
  //   it.getNeighbor(i).getEndIndices() << ", "; std::cout << std::endl;
  //
  //   // << it.getNeighbor(2).getValue() << " -- " << it.getCenter().getValue()
  //   << " -- " << it.getNeighbor(0).getValue() << std::endl;
  // }

  // for(hrleConstSparseIterator<Domain<double, D>::DomainType>
  // it(substrate.getDomain()); !it.isFinished(); it.next()) {
  //   std::cout << it.getStartIndices() << " -- " << it.getEndIndices() <<
  //   std::endl;
  // }

  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, D>(substrate, mesh).apply();
  ls::VTKWriter<double>(mesh, ls::FileFormatEnum::VTP, "normal.vtp").apply();

  // ls::Expand<double, D>(substrate, 4).apply();
  // ls::ToMesh<double, D>(substrate, mesh).apply();
  // ls::VTKWriter<double>(mesh, ls::FileFormatEnum::VTP,
  // "expanded.vtp").apply();

  // Prune<double, D>(substrate).apply();
  // ls::ToMesh<double, D>(substrate, mesh).apply();
  // ls::VTKWriter<double>(mesh, ls::FileFormatEnum::VTP, "pruned.vtp").apply();
  // -----------------------------------------------------

  {
    // create spheres used for booling
    std::cout << "Creating pillar..." << std::endl;
    auto pillar = ls::SmartPointer<ls::Domain<double, D>>::New(
        bounds, boundaryCons, gridDelta);
    double lowerCorner[D] = {5, -1};
    double upperCorner[D] = {15, 10};
    ls::MakeGeometry<double, D>(
        pillar,
        ls::SmartPointer<ls::Box<double, D>>::New(lowerCorner, upperCorner))
        .apply();
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(pillar, mesh).apply();
    ls::VTKWriter<double>(mesh, ls::FileFormatEnum::VTP, "pillar.vtp").apply();
    ls::BooleanOperation<double, D> boolOp(substrate, pillar,
                                           ls::BooleanOperationEnum::UNION);
    boolOp.apply();
  }

  // Now etch the substrate isotropically
  auto velocities = ls::SmartPointer<velocityField>::New();

  std::cout << "Advecting" << std::endl;

  ls::Advect<double, D> advectionKernel;
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
      auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
      ls::ToSurfaceMesh<double, D>(substrate, mesh).apply();
      ls::VTKWriter<double>(mesh, ls::FileFormatEnum::VTP,
                            "pillar-" + std::to_string(i) + ".vtp")
          .apply();

      ls::ToMesh<double, D>(substrate, mesh).apply();
      ls::VTKWriter<double>(mesh, ls::FileFormatEnum::VTP,
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
