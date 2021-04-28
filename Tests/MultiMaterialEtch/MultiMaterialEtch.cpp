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
  Example showing how to grow/shrink different neighbouring materials
  at different speeds.
  \example MultiMaterialAdvection.cpp
*/

// implement own velocity field for advection
// in this case just grow one of the materials
class depositionVel : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int /*material*/,
                           const std::array<double, 3> & /*normalVector*/,
                           unsigned long /*pointId*/) {
    return 0.1;
  }
};

class etchingVel : public lsVelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int material,
                           const std::array<double, 3> & /*normalVector*/,
                           unsigned long /*pointId*/) {
    return (material == 1) ? -0.3 : 0;
  }

  // std::array<double, 3>
  // getVectorVelocity(const std::array<double, 3> & /*coordinate*/,
  //                   int /*material*/,
  //                   const std::array<double, 3> & /*normalVector*/,  unsigned
  //                   long /*pointId*/) {
  //   return std::array<double, 3>({});
  // }
};

int main() {
  omp_set_num_threads(1);

  constexpr int D = 2;
  typedef double NumericType;
  double gridDelta = 1.1;

  double extent = 10;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  if (D == 3) {
    bounds[4] = -extent;
    bounds[5] = extent;
  }

  typename lsDomain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] =
        lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] =
      lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);
  auto mask =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

  double planeOrigin[3] = {0., 0., 0.};
  double planeNormal[3] = {0., D == 2, D == 3};

  lsMakeGeometry<double, D>(substrate, lsSmartPointer<lsPlane<double, D>>::New(
                                           planeOrigin, planeNormal))
      .apply();

  double maskOrigin[3] = {0., -10 * (D == 2), -10 * (D == 3)};
  double maskNormal[3] = {0, -(D == 2), -(D == 3)};

  lsMakeGeometry<double, D>(
      mask, lsSmartPointer<lsPlane<double, D>>::New(maskOrigin, maskNormal))
      .apply();

  lsBooleanOperation<double, D>(mask, substrate,
                                lsBooleanOperationEnum::INTERSECT)
      .apply();

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<double, D>(mask, mesh).apply();
    lsVTKWriter<double>(mesh, "maskPlane.vtk").apply();
  }

  // {
  //   // create box used for mask
  //   std::cout << "Creating box..." << std::endl;
  //   lsDomain<double, D> trench(bounds, boundaryCons, gridDelta);
  //   double minCorner[3] = {-extent / 4., (D==3)?(-extent / 4.):-15, -15.};
  //   double maxCorner[3] = {extent / 4., (D==3)?(extent / 4.):5.0, 5.0};
  //   lsMakeGeometry<double, D>(trench, lsBox<double, D>(minCorner, maxCorner))
  //       .apply();

  //   {
  //     std::cout << "Extracting..." << std::endl;
  //     lsMesh<> mesh;
  //     lsToMesh<double, D>(trench, mesh).apply();
  //     lsVTKWriter<double>(mesh, "box.vtk").apply();
  //   }

  //   // Create trench geometry
  //   std::cout << "Booling trench..." << std::endl;
  //   lsBooleanOperation<double, D>(mask, trench,
  //                                 lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
  //       .apply();

  //   {
  //     std::cout << "Extracting..." << std::endl;
  //     lsMesh<> mesh;
  //     lsToMesh<double, D>(mask, mesh).apply();
  //     lsVTKWriter<double>(mesh, "mask.vtk").apply();
  //   }
  // }

  {
    auto mesh = lsSmartPointer<lsMesh<>>::New();

    lsToMesh<NumericType, D>(substrate, mesh).apply();
    lsVTKWriter<double>(mesh, "points.vtk").apply();
    lsToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    lsVTKWriter<double>(mesh, "surface.vtk").apply();
  }

  auto depoVel = lsSmartPointer<depositionVel>::New();
  auto etchVel = lsSmartPointer<etchingVel>::New();

  std::cout << "Advecting" << std::endl;
  lsAdvect<double, D> deposition(depoVel);
  deposition.insertNextLevelSet(mask);
  deposition.insertNextLevelSet(substrate);
  deposition.setAdvectionTime(1);

  lsAdvect<double, D> etching(etchVel);
  etching.insertNextLevelSet(mask);
  etching.insertNextLevelSet(substrate);
  etching.setAdvectionTime(1);

  {
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<NumericType, D>(mask, mesh).apply();
    lsVTKWriter<double>(mesh, "mask0.vtk").apply();
    lsToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    lsVTKWriter<double>(mesh, "surface0.vtk").apply();
  }

  for (unsigned i = 1; i < 10; ++i) {
    lsAdvect<double, D> deposition(depoVel);
    deposition.insertNextLevelSet(mask);
    deposition.insertNextLevelSet(substrate);
    deposition.setAdvectionTime(1);

    lsAdvect<double, D> etching(etchVel);
    etching.insertNextLevelSet(mask);
    etching.insertNextLevelSet(substrate);
    etching.setAdvectionTime(1);

    deposition.apply();

    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<NumericType, D>(mask, mesh).apply();
    lsVTKWriter<double>(mesh, "mask" + std::to_string(2 * i) + ".vtk").apply();
    lsToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    lsVTKWriter<double>(mesh, "surface" + std::to_string(2 * i) + ".vtk")
        .apply();
    std::cout << "DepoSteps: " << deposition.getNumberOfTimeSteps()
              << std::endl;

    etching.apply();

    lsToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    lsVTKWriter<double>(mesh, "surface" + std::to_string(2 * i + 1) + ".vtk")
        .apply();
    lsToSurfaceMesh<NumericType, D>(mask, mesh).apply();
    lsVTKWriter<double>(mesh, "mask" + std::to_string(2 * i + 1) + ".vtk")
        .apply();
    std::cout << "EtchSteps: " << etching.getNumberOfTimeSteps() << std::endl;
  }

  return 0;
}
