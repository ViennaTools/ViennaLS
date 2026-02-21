#include <iostream>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsToMesh.hpp>
#include <lsToMultiSurfaceMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsWriteVisualizationMesh.hpp>

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
    return 1.;
  }
};

class etchingVel : public ls::VelocityField<double> {
public:
  double getScalarVelocity(const std::array<double, 3> & /*coordinate*/,
                           int material,
                           const std::array<double, 3> & /*normalVector*/,
                           unsigned long /*pointId*/) {
    if (material == 0)
      return 0.;
    else if (material == 1)
      return -10.;
    else
      return -0.1;
  }
};

int main() {
  omp_set_num_threads(1);

  constexpr int D = 2;
  double gridDelta = 0.59;

  double extent = 10;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  if constexpr (D == 3) {
    bounds[4] = -extent;
    bounds[5] = extent;
  }

  typename ls::Domain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] = ls::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

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
    // create box used for mask
    std::cout << "Creating box..." << std::endl;
    auto trench = ls::Domain<double, D>::New(bounds, boundaryCons, gridDelta);
    double minCorner[3] = {-extent / 2., (D == 3) ? (-extent / 4.) : -15, -15.};
    double maxCorner[3] = {extent / 2., (D == 3) ? (extent / 4.) : 5.0, 5.0};
    ls::MakeGeometry<double, D>(trench,
                                ls::Box<double, D>::New(minCorner, maxCorner))
        .apply();

    // Create trench geometry
    std::cout << "Booling trench..." << std::endl;
    ls::BooleanOperation<double, D>(
        mask, trench, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();

    {
      std::cout << "Extracting..." << std::endl;
      auto mesh = ls::Mesh<>::New();
      ls::ToSurfaceMesh<double, D>(mask, mesh).apply();
      ls::VTKWriter<double>(mesh, "mask.vtp").apply();
    }
  }

  auto depoVel = ls::SmartPointer<depositionVel>::New();
  auto etchVel = ls::SmartPointer<etchingVel>::New();

  maskOrigin[D] += gridDelta;
  ls::MakeGeometry<double, D>(
      substrate,
      ls::SmartPointer<ls::Plane<double, D>>::New(maskOrigin, planeNormal))
      .apply();

  ls::BooleanOperation<double, D>(substrate, mask,
                                  ls::BooleanOperationEnum::UNION)
      .apply();

  auto polymer = ls::Domain<double, D>::New(substrate);

  std::cout << "Depositing..." << std::endl;
  ls::Advect<double, D> deposition;
  deposition.setVelocityField(depoVel);
  deposition.insertNextLevelSet(mask);
  deposition.insertNextLevelSet(substrate);
  deposition.insertNextLevelSet(polymer);
  deposition.setAdvectionTime(0.3);
  deposition.apply();

  auto visualizeMesh =
      ls::SmartPointer<ls::WriteVisualizationMesh<double, D>>::New();
  visualizeMesh->insertNextLevelSet(mask);
  visualizeMesh->insertNextLevelSet(substrate);
  visualizeMesh->insertNextLevelSet(polymer);
  visualizeMesh->setSharpCorners(true);
  visualizeMesh->setExtractHullMesh(true);
  visualizeMesh->setFileName("visualizationMesh");

  visualizeMesh->apply();

  {
    auto mesh = ls::Mesh<>::New();
    ls::ToMultiSurfaceMesh<double, D> mesher;
    mesher.insertNextLevelSet(mask);
    mesher.insertNextLevelSet(substrate);
    mesher.insertNextLevelSet(polymer);
    mesher.setMesh(mesh);
    mesher.setSharpCorners(true);
    mesher.apply();
    ls::VTKWriter<double>(mesh, "multiSurfaceMesh").apply();

    ls::ToMesh<double, D>(mask, mesh).apply();
    ls::VTKWriter<double>(mesh, "mask.vtp").apply();
    ls::ToMesh<double, D>(substrate, mesh).apply();
    ls::VTKWriter<double>(mesh, "substrate.vtp").apply();
    ls::ToMesh<double, D>(polymer, mesh).apply();
    ls::VTKWriter<double>(mesh, "polymer.vtp").apply();
  }

  double etchTime = 10.5;

  for (int i = 0; i < 3; ++i) {
    auto m = ls::Domain<double, D>::New(mask);
    auto s = ls::Domain<double, D>::New(substrate);
    auto p = ls::Domain<double, D>::New(polymer);

    ls::Advect<double, D> etching;
    etching.setVelocityField(etchVel);
    etching.insertNextLevelSet(m);
    etching.insertNextLevelSet(s);
    etching.insertNextLevelSet(p);
    etching.setAdvectionTime(etchTime);
    etching.setTimeStepRatio(0.4999 / double(i + 1));
    etching.setAdaptiveTimeStepping(true, 20);

    std::cout << "Etching..." << std::endl;
    etching.apply();

    std::cout << "Time: " << etching.getAdvectedTime()
              << ", Number of steps: " << etching.getNumberOfTimeSteps()
              << std::endl;

    auto mesh = ls::Mesh<>::New();
    ls::ToMultiSurfaceMesh<double, D> mesher;
    mesher.insertNextLevelSet(m);
    mesher.insertNextLevelSet(s);
    mesher.insertNextLevelSet(p);
    mesher.setSharpCorners(true);
    mesher.setMesh(mesh);
    mesher.apply();
    ls::VTKWriter<double>(mesh, "etching_ada_" + std::to_string(i) + ".vtp")
        .apply();
  }

  for (int i = 0; i < 3; ++i) {
    auto m = ls::Domain<double, D>::New(mask);
    auto s = ls::Domain<double, D>::New(substrate);
    auto p = ls::Domain<double, D>::New(polymer);

    ls::Advect<double, D> etching;
    etching.setVelocityField(etchVel);
    etching.insertNextLevelSet(m);
    etching.insertNextLevelSet(s);
    etching.insertNextLevelSet(p);
    etching.setAdvectionTime(etchTime);
    etching.setTimeStepRatio(0.4999 / double(i + 1));

    std::cout << "Etching..." << std::endl;
    etching.apply();

    std::cout << "Time: " << etching.getAdvectedTime()
              << ", Number of steps: " << etching.getNumberOfTimeSteps()
              << std::endl;

    {
      auto mesh = ls::Mesh<>::New();
      ls::ToMultiSurfaceMesh<double, D> mesher;
      mesher.insertNextLevelSet(m);
      mesher.insertNextLevelSet(s);
      mesher.insertNextLevelSet(p);
      mesher.setMesh(mesh);
      mesher.apply();
      ls::VTKWriter<double>(mesh, "etching_noada_" + std::to_string(i) + ".vtp")
          .apply();
    }
  }

  if constexpr (false) {
    auto m = ls::Domain<double, D>::New(mask);
    auto s = ls::Domain<double, D>::New(substrate);
    auto p = ls::Domain<double, D>::New(polymer);

    ls::Advect<double, D> etching;
    etching.setVelocityField(etchVel);
    etching.insertNextLevelSet(m);
    etching.insertNextLevelSet(s);
    etching.insertNextLevelSet(p);
    etching.setAdvectionTime(etchTime);
    etching.setTimeStepRatio(0.25);
    // etching.setAdaptiveTimeStepping(true, 20);
    etching.setSingleStep(true);

    std::cout << "Etching..." << std::endl;
    double time = 0;
    int i = 0;
    while (time < etchTime) {
      etching.apply();
      time += etching.getAdvectedTime();
      auto mesh = ls::Mesh<>::New();
      ls::ToMultiSurfaceMesh<double, D> mesher;
      mesher.insertNextLevelSet(m);
      mesher.insertNextLevelSet(s);
      mesher.insertNextLevelSet(p);
      mesher.setMesh(mesh);
      mesher.apply();
      ls::VTKWriter<double>(mesh, "t_etching_" + std::to_string(i) + ".vtp")
          .apply();
      i++;
    }
  }

  return 0;
}