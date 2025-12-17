#include <iostream>

#include <lsAdvect.hpp>
#include <lsAdvectForwardEuler.hpp>
#include <lsAdvectRungeKutta3.hpp>
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
  2D Example showing how to use the library for topography
  simulation, by creating a trench geometry. A layer of a different material is
  then grown directionally on top. \example AirGapDeposition.cpp
*/

namespace ls = viennals;

using NumericType = float;

// implement own velocity field
class velocityField : public ls::VelocityField<NumericType> {
public:
  NumericType
  getScalarVelocity(const std::array<NumericType, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<NumericType, 3> &normalVector,
                    unsigned long /*pointId*/) {
    // velocity is proportional to the normal vector
    NumericType velocity =
        std::abs(normalVector[0]) + std::abs(normalVector[1]);
    return velocity;
  }

  std::array<NumericType, 3>
  getVectorVelocity(const std::array<NumericType, 3> & /*coordinate*/,
                    int /*material*/,
                    const std::array<NumericType, 3> & /*normalVector*/,
                    unsigned long /*pointId*/) {
    return std::array<NumericType, 3>({});
  }
};

template <class AdvectKernelType, class NumericType, int D>
double runSimulation(AdvectKernelType &kernel,
                     ls::SmartPointer<ls::Domain<NumericType, D>> newLayer,
                     double totalTime, double outputInterval,
                     std::string name) {
  double passedTime = 0.0;
  unsigned stepCounter = 0;

  while (passedTime < totalTime) {
    double dt = outputInterval;
    if (passedTime + dt > totalTime) {
      dt = totalTime - passedTime;
    }

    kernel.setAdvectionTime(dt);
    kernel.apply();
    passedTime += kernel.getAdvectedTime();

    std::cout << "\r" << name << " Advection time: " << std::fixed
              << std::setprecision(1) << std::setw(4) << std::setfill('0')
              << passedTime << " / " << std::setw(4) << std::setfill('0')
              << totalTime << "s" << std::flush;

    auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
    ls::ToSurfaceMesh<NumericType, D>(newLayer, mesh).apply();
    ls::VTKWriter<NumericType> writer(
        mesh, "trench_" + name + "_" + std::to_string(stepCounter) + ".vtp");
    writer.addMetaData("time", passedTime);
    writer.apply();

    ++stepCounter;
  }
  std::cout << std::endl;
  return passedTime;
}

int main() {

  constexpr int D = 2;
  omp_set_num_threads(8);

  NumericType extent = 30;
  NumericType gridDelta = 0.5;

  double bounds[2 * D] = {-extent, extent, -extent, extent};
  ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] =
      ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[1] = ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  NumericType origin[2] = {0., 0.};
  NumericType planeNormal[2] = {0., 1.};

  {
    auto plane =
        ls::SmartPointer<ls::Plane<NumericType, D>>::New(origin, planeNormal);
    ls::MakeGeometry<NumericType, D>(substrate, plane).apply();
  }

  {
    std::cout << "Extracting..." << std::endl;
    auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
    ls::ToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
    ls::VTKWriter<NumericType>(mesh, "plane.vtp").apply();
  }

  {
    // create layer used for booling
    std::cout << "Creating box..." << std::endl;
    auto trench = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);
    NumericType xlimit = extent / 6.;
    NumericType minCorner[D] = {-xlimit, -25.};
    NumericType maxCorner[D] = {xlimit, 1.};
    auto box =
        ls::SmartPointer<ls::Box<NumericType, D>>::New(minCorner, maxCorner);
    ls::MakeGeometry<NumericType, D>(trench, box).apply();

    {
      std::cout << "Extracting..." << std::endl;
      auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
      ls::ToMesh<NumericType, D>(trench, mesh).apply();
      ls::VTKWriter<NumericType>(mesh, "box.vtp").apply();
    }

    // Create trench geometry
    std::cout << "Booling trench..." << std::endl;
    ls::BooleanOperation<NumericType, D>(
        substrate, trench, ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }

  // Now grow new material

  // Create copies for FE and RK simulations
  std::cout << "Creating new layers..." << std::endl;
  auto substrateFE =
      ls::SmartPointer<ls::Domain<NumericType, D>>::New(substrate);
  auto newLayerFE =
      ls::SmartPointer<ls::Domain<NumericType, D>>::New(substrateFE);

  auto substrateRK =
      ls::SmartPointer<ls::Domain<NumericType, D>>::New(substrate);
  auto newLayerRK =
      ls::SmartPointer<ls::Domain<NumericType, D>>::New(substrateRK);

  auto velocities = ls::SmartPointer<velocityField>::New();

  double totalSimulationTime = 16.5;
  double outputInterval = 0.5;

  std::cout << "Advecting" << std::endl;

  // FE Kernel
  ls::AdvectForwardEuler<NumericType, D> advectionKernelFE;
  advectionKernelFE.insertNextLevelSet(substrateFE);
  advectionKernelFE.insertNextLevelSet(newLayerFE);
  advectionKernelFE.setVelocityField(velocities);
  advectionKernelFE.setIgnoreVoids(true);

  double passedTimeFE = runSimulation(
      advectionKernelFE, newLayerFE, totalSimulationTime, outputInterval, "FE");

  // RK Kernel
  ls::AdvectRungeKutta3<NumericType, D> advectionKernelRK;
  advectionKernelRK.insertNextLevelSet(substrateRK);
  advectionKernelRK.insertNextLevelSet(newLayerRK);
  advectionKernelRK.setVelocityField(velocities);
  advectionKernelRK.setIgnoreVoids(true);

  double passedTimeRK = runSimulation(
      advectionKernelRK, newLayerRK, totalSimulationTime, outputInterval, "RK");

  std::cout << "Time passed FE: " << passedTimeFE << std::endl;
  std::cout << "Time passed RK: " << passedTimeRK << std::endl;

  // FE Output
  {
    ls::WriteVisualizationMesh<NumericType, D> writer;
    writer.insertNextLevelSet(substrateFE);
    writer.insertNextLevelSet(newLayerFE);
    writer.addMetaData("time", passedTimeFE);
    writer.setFileName("airgap_FE");
    writer.setExtractHullMesh(true);
    writer.apply();

    ls::ToMultiSurfaceMesh<NumericType, D> multiMeshKernel;
    multiMeshKernel.insertNextLevelSet(substrateFE);
    multiMeshKernel.insertNextLevelSet(newLayerFE);
    auto multiMesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
    multiMeshKernel.setMesh(multiMesh);
    multiMeshKernel.apply();
    ls::VTKWriter<NumericType>(multiMesh, "multimesh_FE.vtp").apply();
  }

  // RK Output
  {
    ls::WriteVisualizationMesh<NumericType, D> writer;
    writer.insertNextLevelSet(substrateRK);
    writer.insertNextLevelSet(newLayerRK);
    writer.addMetaData("time", passedTimeRK);
    writer.setFileName("airgap_RK");
    writer.setExtractHullMesh(true);
    writer.apply();

    ls::ToMultiSurfaceMesh<NumericType, D> multiMeshKernel;
    multiMeshKernel.insertNextLevelSet(substrateRK);
    multiMeshKernel.insertNextLevelSet(newLayerRK);
    auto multiMesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
    multiMeshKernel.setMesh(multiMesh);
    multiMeshKernel.apply();
    ls::VTKWriter<NumericType>(multiMesh, "multimesh_RK.vtp").apply();
  }

  return 0;
}
