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

namespace ls = viennals;

// implement own velocity field for advection
// in this case just grow one of the materials
class depositionVel : public ls::VelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int /*material*/,
                           const std::array<double, 3> & /*normalVector*/,
                           unsigned long /*pointId*/) {
    return 0.1;
  }
};

class etchingVel : public ls::VelocityField<double> {
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
  if constexpr (D == 3) {
    bounds[4] = -extent;
    bounds[5] = extent;
  }

  typename ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] =
        ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] =
      ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  auto mask = ls::SmartPointer<ls::Domain<double, D>>::New(bounds, boundaryCons,
                                                           gridDelta);

  double planeOrigin[3] = {0., 0., 0.};
  double planeNormal[3] = {0., D == 2, D == 3};

  ls::MakeGeometry<double, D>(
      substrate,
      ls::SmartPointer<ls::Plane<double, D>>::New(planeOrigin, planeNormal))
      .apply();

  double maskOrigin[3] = {0., -10 * (D == 2), -10 * (D == 3)};
  double maskNormal[3] = {0, -1.0 * (D == 2), -1.0 * (D == 3)};

  ls::MakeGeometry<double, D>(
      mask, ls::SmartPointer<ls::Plane<double, D>>::New(maskOrigin, maskNormal))
      .apply();

  ls::BooleanOperation<double, D>(mask, substrate,
                                  ls::BooleanOperationEnum::INTERSECT)
      .apply();

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(mask, mesh).apply();
    ls::VTKWriter<double>(mesh, "maskPlane.vtk").apply();
  }

  // {
  //   // create box used for mask
  //   std::cout << "Creating box..." << std::endl;
  //   ls::Domain<double, D> trench(bounds, boundaryCons, gridDelta);
  //   double minCorner[3] = {-extent / 4., (D==3)?(-extent / 4.):-15, -15.};
  //   double maxCorner[3] = {extent / 4., (D==3)?(extent / 4.):5.0, 5.0};
  //   ls::MakeGeometry<double, D>(trench, ls::Box<double, D>(minCorner,
  //   maxCorner))
  //       .apply();

  //   {
  //     std::cout << "Extracting..." << std::endl;
  //     Mesh<> mesh;
  //     ls::ToMesh<double, D>(trench, mesh).apply();
  //     ls::VTKWriter<double>(mesh, "box.vtk").apply();
  //   }

  //   // Create trench geometry
  //   std::cout << "Booling trench..." << std::endl;
  //   ls::BooleanOperation<double, D>(mask, trench,
  //                                 ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
  //       .apply();

  //   {
  //     std::cout << "Extracting..." << std::endl;
  //     Mesh<> mesh;
  //     ls::ToMesh<double, D>(mask, mesh).apply();
  //     ls::VTKWriter<double>(mesh, "mask.vtk").apply();
  //   }
  // }

  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();

    ls::ToMesh<NumericType, D>(substrate, mesh).apply();
    ls::VTKWriter<double>(mesh, "points.vtk").apply();
    ls::ToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    ls::VTKWriter<double>(mesh, "surface.vtk").apply();
  }

  auto depoVel = ls::SmartPointer<depositionVel>::New();
  auto etchVel = ls::SmartPointer<etchingVel>::New();

  std::cout << "Advecting" << std::endl;
  ls::Advect<double, D> deposition;
  deposition.setVelocityField(depoVel);
  deposition.insertNextLevelSet(mask);
  deposition.insertNextLevelSet(substrate);
  deposition.setAdvectionTime(1);

  ls::Advect<double, D> etching;
  etching.setVelocityField(etchVel);
  etching.insertNextLevelSet(mask);
  etching.insertNextLevelSet(substrate);
  etching.setAdvectionTime(1);

  {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<NumericType, D>(mask, mesh).apply();
    ls::VTKWriter<double>(mesh, "mask0.vtk").apply();
    ls::ToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    ls::VTKWriter<double>(mesh, "surface0.vtk").apply();
  }

  for (unsigned i = 1; i < 10; ++i) {
    ls::Advect<double, D> deposition;
    deposition.setVelocityField(depoVel);
    deposition.insertNextLevelSet(mask);
    deposition.insertNextLevelSet(substrate);
    deposition.setAdvectionTime(1);

    ls::Advect<double, D> etching;
    etching.setVelocityField(etchVel);
    etching.insertNextLevelSet(mask);
    etching.insertNextLevelSet(substrate);
    etching.setAdvectionTime(1);

    deposition.apply();

    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<NumericType, D>(mask, mesh).apply();
    ls::VTKWriter<double>(mesh, "mask" + std::to_string(2 * i) + ".vtk")
        .apply();
    ls::ToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    ls::VTKWriter<double>(mesh, "surface" + std::to_string(2 * i) + ".vtk")
        .apply();
    std::cout << "DepoSteps: " << deposition.getNumberOfTimeSteps()
              << std::endl;

    etching.apply();

    ls::ToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    ls::VTKWriter<double>(mesh, "surface" + std::to_string(2 * i + 1) + ".vtk")
        .apply();
    ls::ToSurfaceMesh<NumericType, D>(mask, mesh).apply();
    ls::VTKWriter<double>(mesh, "mask" + std::to_string(2 * i + 1) + ".vtk")
        .apply();
    std::cout << "EtchSteps: " << etching.getNumberOfTimeSteps() << std::endl;
  }

  return 0;
}
