#include <array>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <string>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsOxidationMaterials.hpp>
#include <lsOxidationModel.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

using NumericType = double;

template <int D>
using LevelSet = ls::SmartPointer<ls::Domain<NumericType, D>>;

template <int D>
LevelSet<D> makeBoxLevelSet(const double *bounds,
                         typename ls::Domain<NumericType, D>::BoundaryType *boundaryCons,
                         NumericType gridDelta,
                         const ls::VectorType<NumericType, D> &minCorner,
                         const ls::VectorType<NumericType, D> &maxCorner) {
  auto levelSet = ls::Domain<NumericType, D>::New(bounds, boundaryCons, gridDelta);
  ls::MakeGeometry<NumericType, D> makeGeometry(
      levelSet, ls::Box<NumericType, D>::New(minCorner, maxCorner));
  std::array<bool, D> ignoreBCs;
  ignoreBCs.fill(false);
  ignoreBCs[1] = true;
  makeGeometry.setIgnoreBoundaryConditions(ignoreBCs);
  makeGeometry.apply();
  return levelSet;
}

template <int D>
LevelSet<D> makeStepLevelSet(const double *bounds,
                          typename ls::Domain<NumericType, D>::BoundaryType *boundaryCons,
                          NumericType gridDelta, NumericType xMin,
                          NumericType xMax, NumericType yMin,
                          NumericType leftTop, NumericType rightTop,
                          NumericType stepX, NumericType zMin = 0.,
                          NumericType zMax = 0., NumericType stepZ = 0.) {
  (void)xMin;
  (void)yMin;
  (void)zMin;

  auto step = ls::Domain<NumericType, D>::New(bounds, boundaryCons, gridDelta);
  ls::VectorType<NumericType, D> planeOrigin(0.);
  planeOrigin[1] = leftTop;
  ls::VectorType<NumericType, D> planeNormal(0.);
  planeNormal[1] = 1.;
  ls::MakeGeometry<NumericType, D>(
      step, ls::Plane<NumericType, D>::New(planeOrigin, planeNormal))
      .apply();

  ls::VectorType<NumericType, D> rightMin(0.);
  rightMin[0] = stepX; rightMin[1] = leftTop;
  if constexpr (D == 3) rightMin[2] = stepZ;

  ls::VectorType<NumericType, D> rightMax(0.);
  rightMax[0] = xMax; rightMax[1] = rightTop;
  if constexpr (D == 3) rightMax[2] = zMax;
  auto rightBlock =
      makeBoxLevelSet(bounds, boundaryCons, gridDelta, rightMin, rightMax);

  ls::BooleanOperation<NumericType, D>(
      step, rightBlock, ls::BooleanOperationEnum::UNION)
      .apply();
  return step;
}

template <int D>
void writeSurface(LevelSet<D> levelSet, const std::string &fileName) {
  auto mesh = ls::Mesh<NumericType>::New();
  auto surfaceMesh = ls::ToSurfaceMesh<NumericType, D>(levelSet, mesh);
  surfaceMesh.setSharpCorners(false);
  surfaceMesh.apply();
  ls::VTKWriter<NumericType>(mesh, fileName).apply();
}

template <int D>
void writeDeformationCSV(
    const ls::SmartPointer<
        ls::OxidationDeformation<NumericType, D>> &deformation,
    const viennahrle::Index<D> &minIndex, const viennahrle::Index<D> &maxIndex,
    NumericType gridDelta, const std::string &fileName) {
  std::ofstream file(fileName);
  if constexpr (D == 2) {
    file << "x,y,velocity_x,velocity_y,pressure,strain_trace,"
            "stress_xx,stress_xy,stress_yy,von_mises_stress\n";
    for (auto j = minIndex[1]; j <= maxIndex[1]; ++j) {
      for (auto i = minIndex[0]; i <= maxIndex[0]; ++i) {
        ls::Vec3D<NumericType> coordinate{i * gridDelta, j * gridDelta, 0.};
        const auto velocity = deformation->getVelocity(coordinate);
        const auto stress = deformation->getStressTensor(coordinate);
        file << coordinate[0] << ',' << coordinate[1] << ',' << velocity[0]
             << ',' << velocity[1] << ',' << deformation->getPressure(coordinate)
             << ',' << deformation->getStrainTrace(coordinate) << ','
             << stress[0] << ',' << stress[1] << ',' << stress[4] << ','
             << deformation->getVonMisesStress(coordinate) << '\n';
      }
    }
  } else {
    file << "x,y,z,velocity_x,velocity_y,velocity_z,pressure,strain_trace,"
            "stress_xx,stress_xy,stress_xz,stress_yy,stress_yz,stress_zz,von_mises_stress\n";
    for (auto k = minIndex[2]; k <= maxIndex[2]; ++k) {
      for (auto j = minIndex[1]; j <= maxIndex[1]; ++j) {
        for (auto i = minIndex[0]; i <= maxIndex[0]; ++i) {
          ls::Vec3D<NumericType> coordinate{i * gridDelta, j * gridDelta, k * gridDelta};
          const auto velocity = deformation->getVelocity(coordinate);
          const auto stress = deformation->getStressTensor(coordinate);
          file << coordinate[0] << ',' << coordinate[1] << ',' << coordinate[2] << ',' 
               << velocity[0] << ',' << velocity[1] << ',' << velocity[2] << ','
               << deformation->getPressure(coordinate) << ',' 
               << deformation->getStrainTrace(coordinate) << ','
               << stress[0] << ',' << stress[1] << ',' << stress[2] << ',' 
               << stress[4] << ',' << stress[5] << ',' << stress[8] << ','
               << deformation->getVonMisesStress(coordinate) << '\n';
        }
      }
    }
  }
}

template <int D>
int runSimulation() {
  constexpr NumericType gridDelta = 0.05; // micrometer
  constexpr NumericType xExtent = 5.;     // micrometer
  constexpr NumericType yMin = -2.;       // micrometer
  constexpr NumericType yMax = 4.;        // micrometer
  constexpr NumericType zExtent = 5.;     // micrometer
  constexpr NumericType stepX = 0.;
  constexpr NumericType stepZ = 0.;
  constexpr NumericType leftSiTop = 0.;
  constexpr NumericType rightSiTop = 1.;
  constexpr NumericType oxideThickness = 0.2;
  constexpr NumericType advectionTime = 0.1; // hr

  double bounds[2 * D] = {-xExtent, xExtent, yMin, yMax};
  if constexpr (D == 3) {
    bounds[4] = -zExtent; bounds[5] = zExtent;
  }
  typename ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[1] = ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;
  if constexpr (D == 3) {
    boundaryCons[2] = ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }

  auto siInterface =
      makeStepLevelSet(bounds, boundaryCons, gridDelta, -xExtent, xExtent, yMin,
                       leftSiTop, rightSiTop, stepX, -zExtent, zExtent, stepZ);

  auto ambientInterface = ls::Domain<NumericType, D>::New(siInterface);
  auto initialOxide =
      ls::SmartPointer<ls::SphereDistribution<viennahrle::CoordType, D>>::New(
          oxideThickness);
  ls::GeometricAdvect<NumericType, D>(ambientInterface, initialOxide).apply();

  writeSurface(siInterface, "step_oxidation_si_initial.vtp");
  writeSurface(ambientInterface, "step_oxidation_ambient_initial.vtp");

  auto params =
      ls::OxidationMaterials<NumericType>::wet1000CDealGrove100();
  // The silicon step is represented as the negative/inside side of this level
  // set. A negative velocity consumes silicon and moves the Si/SiO2 interface
  // into the step; the oxide/ambient interface is moved by the deformation
  // velocity field below.
  params.velocitySign = -1.0;
  params.maxIterations = 10000;
  params.tolerance = 1e-7;

  auto oxidationVelocity =
      ls::OxidationDiffusion<NumericType, D>::New(
          siInterface, ambientInterface, params);
  viennahrle::Index<D> minIndex;
  viennahrle::Index<D> maxIndex;
  minIndex[0] = -100; maxIndex[0] = 100;
  minIndex[1] = -40;  maxIndex[1] = 80;
  if constexpr (D == 3) {
    minIndex[2] = -100;
    maxIndex[2] = 100;
  }
  oxidationVelocity->setSolveBounds(minIndex, maxIndex);
  auto deformationParams =
      ls::OxidationMaterials<NumericType>::oxideMechanics1000C(
          advectionTime);

  auto deformationVelocity =
      ls::OxidationDeformation<NumericType, D>::New(
          siInterface, ambientInterface, oxidationVelocity, params,
          deformationParams);
  deformationVelocity->setSolveBounds(minIndex, maxIndex);

  ls::OxidationCouplingParameters<NumericType> couplingParams;
  couplingParams.maxIterations = 8;
  couplingParams.tolerance = 1e-6;
  couplingParams.relaxation = 1.0;
  auto coupledOxidation =
      ls::OxidationModel<NumericType, D>::New(
          oxidationVelocity, deformationVelocity, couplingParams);
  coupledOxidation->setSolveBounds(minIndex, maxIndex);
  coupledOxidation->apply();

  std::cout << "Coupled oxidation iterations: "
            << coupledOxidation->getIterations()
            << ", relative pressure residual: "
            << coupledOxidation->getResidual()
            << std::endl;
  std::cout << "Diffusion solve nodes: "
            << oxidationVelocity->getNumberOfSolutionNodes()
            << ", iterations: " << oxidationVelocity->getIterations()
            << ", residual: " << oxidationVelocity->getResidual()
            << std::endl;

  std::cout << "Deformation solve nodes: "
            << deformationVelocity->getNumberOfSolutionNodes()
            << ", iterations: " << deformationVelocity->getIterations()
            << ", residual: " << deformationVelocity->getResidual()
            << std::endl;
  const auto ambientFlatSpeed = deformationVelocity->avgExpansionSpeed();
  const auto siliconFlatSpeed =
      ambientFlatSpeed / (params.expansionCoefficient - 1.);
  std::cout << "Flat split benchmark: silicon fraction "
            << siliconFlatSpeed / (siliconFlatSpeed + ambientFlatSpeed)
            << ", oxide/ambient fraction "
            << ambientFlatSpeed / (siliconFlatSpeed + ambientFlatSpeed)
            << std::endl;
  ls::Vec3D<NumericType> flatReactionPoint{2., rightSiTop + gridDelta, 0.};
  ls::Vec3D<NumericType> flatAmbientPoint{2., rightSiTop + oxideThickness, 0.};
  if constexpr (D == 3) {
    flatReactionPoint[2] = 2.;
    flatAmbientPoint[2] = 2.;
  }
  const NumericType localSiliconSpeed =
      std::abs(oxidationVelocity->getScalarVelocity(flatReactionPoint, 0,
                                                    {0., 1., 0.}, 0));
  const NumericType localAmbientSpeed =
      std::abs(deformationVelocity->getScalarVelocity(flatAmbientPoint, 0,
                                                      {0., 1., 0.}, 0));
  std::cout << "Flat top velocity split: silicon "
            << localSiliconSpeed / (localSiliconSpeed + localAmbientSpeed)
            << ", oxide/ambient "
            << localAmbientSpeed / (localSiliconSpeed + localAmbientSpeed)
            << "; expected flat displacements after " << advectionTime
            << " hr: silicon " << localSiliconSpeed * advectionTime
            << " um, oxide/ambient " << localAmbientSpeed * advectionTime
            << " um"
            << std::endl;

  writeDeformationCSV(deformationVelocity, minIndex, maxIndex, gridDelta,
                      "step_oxidation_deformation.csv");

  ls::Advect<NumericType, D> ambientAdvection;
  ambientAdvection.insertNextLevelSet(ambientInterface);
  ambientAdvection.setVelocityField(deformationVelocity);
  ambientAdvection.setSpatialScheme(
      ls::SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER);
  ambientAdvection.setTemporalScheme(
      ls::TemporalSchemeEnum::FORWARD_EULER);
  ambientAdvection.setAdvectionTime(advectionTime);
  ambientAdvection.apply();

  ls::Advect<NumericType, D> reactionAdvection;
  reactionAdvection.insertNextLevelSet(siInterface);
  reactionAdvection.setVelocityField(oxidationVelocity);
  reactionAdvection.setSpatialScheme(
      ls::SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER);
  reactionAdvection.setTemporalScheme(
      ls::TemporalSchemeEnum::FORWARD_EULER);
  reactionAdvection.setAdvectionTime(advectionTime);
  reactionAdvection.apply();

  writeSurface(siInterface, "step_oxidation_si_after.vtp");
  writeSurface(ambientInterface, "step_oxidation_ambient_after.vtp");

  std::cout << "Wrote step_oxidation_si_initial.vtp, "
               "step_oxidation_ambient_initial.vtp, "
               "step_oxidation_si_after.vtp, "
               "step_oxidation_ambient_after.vtp, and "
               "step_oxidation_deformation.csv"
            << std::endl;
  std::cout << "The ambient interface is moved by the oxide deformation field; "
               "pressure, strain_trace, stress components, and von Mises "
               "stress are finite-difference mechanics diagnostics."
            << std::endl;

  return 0;
}

int main(int argc, char **argv) {
  omp_set_num_threads(4);
  int dim = 2; // Default to 2D
  if (argc > 1) {
    dim = std::stoi(argv[1]);
  }

  if (dim == 2) return runSimulation<2>();
  if (dim == 3) return runSimulation<3>();
  
  std::cerr << "Unsupported dimension: " << dim << std::endl;
  return 1;
}
