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
#include <lsOxidationModel.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

using NumericType = double;
constexpr int D = 2;

using LevelSet = ls::SmartPointer<ls::Domain<NumericType, D>>;

LevelSet makeBoxLevelSet(const double *bounds,
                         ls::Domain<NumericType, D>::BoundaryType *boundaryCons,
                         NumericType gridDelta,
                         const ls::VectorType<NumericType, D> &minCorner,
                         const ls::VectorType<NumericType, D> &maxCorner) {
  auto levelSet = ls::Domain<NumericType, D>::New(bounds, boundaryCons, gridDelta);
  ls::MakeGeometry<NumericType, D> makeGeometry(
      levelSet, ls::Box<NumericType, D>::New(minCorner, maxCorner));
  makeGeometry.setIgnoreBoundaryConditions(std::array<bool, D>{false, true});
  makeGeometry.apply();
  return levelSet;
}

LevelSet makeStepLevelSet(const double *bounds,
                          ls::Domain<NumericType, D>::BoundaryType *boundaryCons,
                          NumericType gridDelta, NumericType xMin,
                          NumericType xMax, NumericType yMin,
                          NumericType leftTop, NumericType rightTop,
                          NumericType stepX) {
  (void)xMin;
  (void)yMin;

  auto step = ls::Domain<NumericType, D>::New(bounds, boundaryCons, gridDelta);
  const ls::VectorType<NumericType, D> planeOrigin{0., leftTop};
  const ls::VectorType<NumericType, D> planeNormal{0., 1.};
  ls::MakeGeometry<NumericType, D>(
      step, ls::Plane<NumericType, D>::New(planeOrigin, planeNormal))
      .apply();

  ls::VectorType<NumericType, D> rightMin{stepX, leftTop};
  ls::VectorType<NumericType, D> rightMax{xMax, rightTop};
  auto rightBlock =
      makeBoxLevelSet(bounds, boundaryCons, gridDelta, rightMin, rightMax);

  ls::BooleanOperation<NumericType, D>(
      step, rightBlock, ls::BooleanOperationEnum::UNION)
      .apply();
  return step;
}

void writeSurface(LevelSet levelSet, const std::string &fileName) {
  auto mesh = ls::Mesh<NumericType>::New();
  auto surfaceMesh = ls::ToSurfaceMesh<NumericType, D>(levelSet, mesh);
  surfaceMesh.setSharpCorners(true);
  surfaceMesh.apply();
//   ls::ToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  ls::VTKWriter<NumericType>(mesh, fileName).apply();
}

void writeDeformationCSV(
    const ls::SmartPointer<
        ls::OxidationDeformationVelocityField<NumericType, D>> &deformation,
    const viennahrle::Index<D> &minIndex, const viennahrle::Index<D> &maxIndex,
    NumericType gridDelta, const std::string &fileName) {
  std::ofstream file(fileName);
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
}

int main() {
  omp_set_num_threads(4);

  constexpr NumericType gridDelta = 0.05; // micrometer
  constexpr NumericType xExtent = 5.;     // micrometer
  constexpr NumericType yMin = -2.;       // micrometer
  constexpr NumericType yMax = 4.;        // micrometer
  constexpr NumericType stepX = 0.;
  constexpr NumericType leftSiTop = 0.;
  constexpr NumericType rightSiTop = 1.;
  constexpr NumericType oxideThickness = 0.2;
  constexpr NumericType advectionTime = 0.1; // hr

  double bounds[2 * D] = {-xExtent, xExtent, yMin, yMax};
  ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[1] = ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto siInterface =
      makeStepLevelSet(bounds, boundaryCons, gridDelta, -xExtent, xExtent, yMin,
                       leftSiTop, rightSiTop, stepX);

  auto ambientInterface = ls::Domain<NumericType, D>::New(siInterface);
  auto initialOxide =
      ls::SmartPointer<ls::SphereDistribution<viennahrle::CoordType, D>>::New(
          oxideThickness);
  ls::GeometricAdvect<NumericType, D>(ambientInterface, initialOxide).apply();

  writeSurface(siInterface, "step_oxidation_si_initial.vtk");
  writeSurface(ambientInterface, "step_oxidation_ambient_initial.vtk");

  ls::OxidationParameters<NumericType> params;
  // Wet oxidation of <100> Si at 1000 C, represented in Deal-Grove form.
  // With C*/N normalized to 1: B = 2D and B/A ~= k for large gas transfer h.
  // Literature values give B ~= 0.31 um^2/hr and B/A ~= 0.74 um/hr.
  params.diffusionCoefficient = 0.157;
  params.reactionRate = 0.74;
  params.transferCoefficient = 100.0;
  params.equilibriumConcentration = 1.0;
  params.oxidantMoleculeDensity = 1.0;
  params.expansionCoefficient = 2.27;
  params.stressCouplingCoefficient = 1.0e-15; // 1/Pa
  params.minStressRateFactor = 0.25;
  params.maxStressRateFactor = 4.0;
  // The silicon step is represented as the negative/inside side of this level
  // set. A negative velocity consumes silicon and moves the Si/SiO2 interface
  // into the step; the oxide/ambient interface is moved by the deformation
  // velocity field below.
  params.velocitySign = -1.0;
  params.maxIterations = 10000;
  params.tolerance = 1e-7;

  auto oxidationVelocity =
      ls::OxidationDiffusionVelocityField<NumericType, D>::New(
          siInterface, ambientInterface, params);
  viennahrle::Index<D> minIndex{-100, -40};
  viennahrle::Index<D> maxIndex{100, 80};
  oxidationVelocity->setSolveBounds(minIndex, maxIndex);
  ls::OxidationDeformationParameters<NumericType> deformationParams;
  deformationParams.viscosity = 1.0e7;    // Pa hr, diagnostic scale
  deformationParams.bulkModulus = 7.5e8;  // Pa, pressure reference from paper
  deformationParams.shearModulus = 3.0e10;
  deformationParams.stressTimeStep = advectionTime;
  deformationParams.freeSurfaceTractionScale = 1.0;
  deformationParams.substrateNormalStiffness = 1.0e9; // Pa/um elastic support
  deformationParams.pressureGradientScale = 1e-3;
  deformationParams.mechanicsIterations = 2;
  deformationParams.mechanicsTolerance = 1e-7;
  deformationParams.pressureIterations = 500;
  deformationParams.stokesIterations = 100;
  deformationParams.pressureTolerance = 1e-6;
  deformationParams.stokesTolerance = 1e-7;
  // Advect the free oxide surface by the full Stokes vector field. The
  // freeSurfaceVelocityScale kinematic-split path is disabled here; both
  // approaches are physically consistent but the Stokes path captures lateral
  // flow and shear that the scalar kinematic split misses.
  deformationParams.freeSurfaceVelocityScale = 0.0;
  deformationParams.vectorVelocityScale = 1.0;
  deformationParams.maxIterations = 10000;
  deformationParams.tolerance = 1e-7;

  auto deformationVelocity =
      ls::OxidationDeformationVelocityField<NumericType, D>::New(
          siInterface, ambientInterface, oxidationVelocity, params,
          deformationParams);
  deformationVelocity->setSolveBounds(minIndex, maxIndex);

  ls::OxidationCouplingParameters<NumericType> couplingParams;
  couplingParams.maxIterations = 8;
  couplingParams.tolerance = 1e-6;
  couplingParams.relaxation = 1.0;
  auto coupledOxidation =
      ls::OxidationCoupledModel<NumericType, D>::New(
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
  const auto ambientFlatSpeed =
      deformationVelocity->getAverageBoundaryExpansionVelocity();
  const auto siliconFlatSpeed =
      ambientFlatSpeed / (params.expansionCoefficient - 1.);
  std::cout << "Flat split benchmark: silicon fraction "
            << siliconFlatSpeed / (siliconFlatSpeed + ambientFlatSpeed)
            << ", oxide/ambient fraction "
            << ambientFlatSpeed / (siliconFlatSpeed + ambientFlatSpeed)
            << std::endl;
  const ls::Vec3D<NumericType> flatReactionPoint{2., rightSiTop + gridDelta,
                                                 0.};
  const ls::Vec3D<NumericType> flatAmbientPoint{2., rightSiTop + oxideThickness,
                                                0.};
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
      ls::TemporalSchemeEnum::RUNGE_KUTTA_2ND_ORDER);
  ambientAdvection.setAdvectionTime(advectionTime);
  ambientAdvection.apply();

  ls::Advect<NumericType, D> reactionAdvection;
  reactionAdvection.insertNextLevelSet(siInterface);
  reactionAdvection.setVelocityField(oxidationVelocity);
  reactionAdvection.setSpatialScheme(
      ls::SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER);
  reactionAdvection.setTemporalScheme(
      ls::TemporalSchemeEnum::RUNGE_KUTTA_2ND_ORDER);
  reactionAdvection.setAdvectionTime(advectionTime);
  reactionAdvection.apply();

  writeSurface(siInterface, "step_oxidation_si_after.vtk");
  writeSurface(ambientInterface, "step_oxidation_ambient_after.vtk");

  std::cout << "Wrote step_oxidation_si_initial.vtk, "
               "step_oxidation_ambient_initial.vtk, "
               "step_oxidation_si_after.vtk, "
               "step_oxidation_ambient_after.vtk, and "
               "step_oxidation_deformation.csv"
            << std::endl;
  std::cout << "The ambient interface is moved by the oxide deformation field; "
               "pressure, strain_trace, stress components, and von Mises "
               "stress are finite-difference mechanics diagnostics."
            << std::endl;

  return 0;
}
