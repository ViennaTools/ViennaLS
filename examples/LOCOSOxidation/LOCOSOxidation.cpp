#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include <lsAdvect.hpp>
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

void writeSurface(LevelSet levelSet, const std::string &fileName) {
  auto mesh = ls::Mesh<NumericType>::New();
  ls::ToSurfaceMesh<NumericType, D>(levelSet, mesh).apply();
  ls::VTKWriter<NumericType>(mesh, fileName).apply();
}

LevelSet makePlane(const double *bounds,
                   ls::Domain<NumericType, D>::BoundaryType *boundaryCons,
                   NumericType gridDelta, NumericType y) {
  auto levelSet = ls::Domain<NumericType, D>::New(bounds, boundaryCons, gridDelta);
  const ls::VectorType<NumericType, D> origin{0., y};
  const ls::VectorType<NumericType, D> normal{0., 1.};
  ls::MakeGeometry<NumericType, D>(
      levelSet, ls::Plane<NumericType, D>::New(origin, normal))
      .apply();
  return levelSet;
}

LevelSet makeMask(const double *bounds,
                  ls::Domain<NumericType, D>::BoundaryType *boundaryCons,
                  NumericType gridDelta, NumericType xMin, NumericType xMax,
                  NumericType yMin, NumericType yMax) {
  auto mask = ls::Domain<NumericType, D>::New(bounds, boundaryCons, gridDelta);
  const ls::VectorType<NumericType, D> minCorner{xMin, yMin};
  const ls::VectorType<NumericType, D> maxCorner{xMax, yMax};
  ls::MakeGeometry<NumericType, D> makeGeometry(
      mask, ls::Box<NumericType, D>::New(minCorner, maxCorner));
  makeGeometry.setIgnoreBoundaryConditions(std::array<bool, D>{false, true});
  makeGeometry.apply();
  return mask;
}

void writeDiagnostics(
    const ls::SmartPointer<ls::OxidationDiffusionVelocityField<NumericType, D>>
        &diffusion,
    const ls::SmartPointer<ls::OxidationDeformationVelocityField<NumericType, D>>
        &deformation,
    const viennahrle::Index<D> &minIndex, const viennahrle::Index<D> &maxIndex,
    NumericType gridDelta, const std::string &fileName) {
  std::ofstream file(fileName);
  file << "x,y,concentration,velocity_x,velocity_y,pressure,strain_trace,"
          "von_mises_stress\n";
  for (auto j = minIndex[1]; j <= maxIndex[1]; ++j) {
    for (auto i = minIndex[0]; i <= maxIndex[0]; ++i) {
      ls::Vec3D<NumericType> coordinate{i * gridDelta, j * gridDelta, 0.};
      const auto velocity = deformation->getVelocity(coordinate);
      file << coordinate[0] << ',' << coordinate[1] << ','
           << diffusion->getConcentration(coordinate) << ',' << velocity[0]
           << ',' << velocity[1] << ',' << deformation->getPressure(coordinate)
           << ',' << deformation->getStrainTrace(coordinate) << ','
           << deformation->getVonMisesStress(coordinate) << '\n';
    }
  }
}

int main() {
  omp_set_num_threads(4);

  constexpr NumericType gridDelta = 0.05;       // um
  constexpr NumericType xExtent = 4.;           // um
  constexpr NumericType yMin = -1.;             // um
  constexpr NumericType yMax = 2.;              // um
  constexpr NumericType padOxideThickness = 0.15;
  constexpr NumericType maskThickness = 0.2;     // um
  constexpr NumericType maskEdge = 0.;          // open window is x > 0
  constexpr NumericType advectionTime = 0.35;   // hr
  constexpr NumericType maskContactEpsilon = 1.e-6; // um

  double bounds[2 * D] = {-xExtent, xExtent, yMin, yMax};
  ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[1] = ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto siInterface = makePlane(bounds, boundaryCons, gridDelta, 0.);
  auto ambientInterface = ls::Domain<NumericType, D>::New(siInterface);
  auto initialOxide =
      ls::SmartPointer<ls::SphereDistribution<viennahrle::CoordType, D>>::New(
          padOxideThickness);
  ls::GeometricAdvect<NumericType, D>(ambientInterface, initialOxide).apply();

  // The mask bottom is numerically offset by a tiny contact epsilon from the
  // pad-oxide/ambient interface. This is far below the grid resolution, so the
  // generated mask lies flat on the oxide surface while Cartesian stencils
  // unambiguously hit the mask boundary in the covered region.
  auto maskInterface =
      makeMask(bounds, boundaryCons, gridDelta, -xExtent, maskEdge,
               padOxideThickness - maskContactEpsilon,
               padOxideThickness + maskThickness);

  writeSurface(siInterface, "locos_si_initial.vtk");
  writeSurface(ambientInterface, "locos_ambient_initial.vtk");
  writeSurface(maskInterface, "locos_mask.vtk");

  ls::OxidationParameters<NumericType> params;
  params.diffusionCoefficient = 0.157;
  params.reactionRate = 0.74;
  params.transferCoefficient = 100.;
  params.equilibriumConcentration = 1.;
  params.oxidantMoleculeDensity = 1.;
  params.expansionCoefficient = 2.27;
  params.velocitySign = -1.;
  params.stressCouplingCoefficient = 1.e-15;
  params.minStressRateFactor = 0.25;
  params.maxStressRateFactor = 4.;
  params.maskTransferCoefficient = 0.; // nitride-like oxidant blocking
  params.maskConcentration = 0.;
  params.maxIterations = 10000;
  params.tolerance = 1.e-7;

  auto oxidationVelocity =
      ls::OxidationDiffusionVelocityField<NumericType, D>::New(
          siInterface, ambientInterface, params);
  oxidationVelocity->setMaskInterface(maskInterface, -1);
  viennahrle::Index<D> minIndex{-80, -20};
  viennahrle::Index<D> maxIndex{80, 40};
  oxidationVelocity->setSolveBounds(minIndex, maxIndex);

  ls::OxidationDeformationParameters<NumericType> deformationParams;
  deformationParams.viscosity = 1.e7;
  deformationParams.bulkModulus = 7.5e8;
  deformationParams.shearModulus = 3.e10;
  deformationParams.stressTimeStep = advectionTime;
  deformationParams.freeSurfaceTractionScale = 1.;
  deformationParams.substrateNormalStiffness = 1.e9;
  deformationParams.maskVelocityScale = 0.15;
  deformationParams.maskNormalStiffness = 2.e9;
  deformationParams.pressureGradientScale = 1.e-3;
  deformationParams.mechanicsIterations = 2;
  deformationParams.mechanicsTolerance = 1.e-7;
  deformationParams.pressureIterations = 500;
  deformationParams.stokesIterations = 100;
  deformationParams.pressureTolerance = 1.e-6;
  deformationParams.stokesTolerance = 1.e-7;
  deformationParams.freeSurfaceVelocityScale = 1.;
  deformationParams.vectorVelocityScale = 0.;
  deformationParams.maxIterations = 10000;
  deformationParams.tolerance = 1.e-7;

  auto deformationVelocity =
      ls::OxidationDeformationVelocityField<NumericType, D>::New(
          siInterface, ambientInterface, oxidationVelocity, params,
          deformationParams);
  deformationVelocity->setMaskInterface(maskInterface, -1);
  deformationVelocity->setSolveBounds(minIndex, maxIndex);

  ls::OxidationCouplingParameters<NumericType> couplingParams;
  couplingParams.maxIterations = 8;
  couplingParams.tolerance = 1.e-6;
  couplingParams.relaxation = 1.;
  auto coupledOxidation =
      ls::OxidationCoupledModel<NumericType, D>::New(
          oxidationVelocity, deformationVelocity, couplingParams);
  coupledOxidation->setSolveBounds(minIndex, maxIndex);
  coupledOxidation->apply();

  const ls::Vec3D<NumericType> openReactionPoint{1.5, gridDelta, 0.};
  const ls::Vec3D<NumericType> maskedReactionPoint{-1.5, gridDelta, 0.};
  const NumericType openConcentration =
      oxidationVelocity->getConcentration(openReactionPoint);
  const NumericType maskedConcentration =
      oxidationVelocity->getConcentration(maskedReactionPoint);
  const NumericType openSiliconSpeed =
      std::abs(oxidationVelocity->getScalarVelocity(openReactionPoint, 0,
                                                    {0., 1., 0.}, 0));
  const NumericType maskedSiliconSpeed =
      std::abs(oxidationVelocity->getScalarVelocity(maskedReactionPoint, 0,
                                                    {0., 1., 0.}, 0));

  std::cout << "Coupled oxidation iterations: "
            << coupledOxidation->getIterations()
            << ", residual: " << coupledOxidation->getResidual() << '\n';
  std::cout << "Diffusion nodes: "
            << oxidationVelocity->getNumberOfSolutionNodes()
            << ", iterations: " << oxidationVelocity->getIterations()
            << ", residual: " << oxidationVelocity->getResidual() << '\n';
  std::cout << "Deformation nodes: "
            << deformationVelocity->getNumberOfSolutionNodes()
            << ", iterations: " << deformationVelocity->getIterations()
            << ", residual: " << deformationVelocity->getResidual() << '\n';
  std::cout << "Average oxide expansion speed: "
            << deformationVelocity->getAverageBoundaryExpansionVelocity()
            << " um/hr\n";
  std::cout << "Open-window concentration: " << openConcentration
            << ", masked concentration: " << maskedConcentration << '\n';
  std::cout << "Open-window Si speed: " << openSiliconSpeed
            << " um/hr, masked Si speed: " << maskedSiliconSpeed
            << " um/hr, suppression ratio: "
            << (openSiliconSpeed > 0. ? maskedSiliconSpeed / openSiliconSpeed
                                      : 0.)
            << '\n';

  if (!std::isfinite(openSiliconSpeed) || !std::isfinite(maskedSiliconSpeed) ||
      openSiliconSpeed <= 0. ||
      maskedSiliconSpeed / openSiliconSpeed > 0.05) {
    std::cerr << "LOCOS mask sanity check failed: masked oxidation is not "
                 "sufficiently suppressed.\n";
    return 1;
  }

  writeDiagnostics(oxidationVelocity, deformationVelocity, minIndex, maxIndex,
                   gridDelta, "locos_oxidation_diagnostics.csv");

  ls::OxidationMaskParameters<NumericType> maskBendingParams;
  maskBendingParams.youngModulus = 2.5e11;       // Pa, LPCVD Si3N4-like
  maskBendingParams.poissonRatio = 0.27;
  maskBendingParams.thickness = maskThickness;
  maskBendingParams.referenceThickness = 0.25;   // um
  maskBendingParams.velocityScale = 0.35;
  maskBendingParams.pressureVelocityScale = 2.e-7; // um / (hr Pa)
  maskBendingParams.maxVelocity = 0.25;          // um/hr
  maskBendingParams.relaxation = 0.9;
  maskBendingParams.tolerance = 5.e-6;
  auto maskBendingVelocity =
      ls::OxidationMaskBendingVelocityField<NumericType, D>::New(
          deformationVelocity, maskInterface, maskBendingParams, 1);
  viennahrle::Index<D> maskMinIndex{-80, 3};
  viennahrle::Index<D> maskMaxIndex{0, 8};
  maskBendingVelocity->setSolveBounds(maskMinIndex, maskMaxIndex);
  maskBendingVelocity->apply();
  const ls::Vec3D<NumericType> maskBottomSample{-1.5, padOxideThickness, 0.};
  const ls::Vec3D<NumericType> maskContactNodeSample{-1.5, 0.2, 0.};
  const ls::Vec3D<NumericType> maskTopSample{
      -1.5, padOxideThickness + maskThickness, 0.};
  const auto maskBottomVelocity = maskBendingVelocity->getVectorVelocity(
      maskBottomSample, 0, {0., -1., 0.}, 0);
  const auto maskTopVelocity = maskBendingVelocity->getVectorVelocity(
      maskTopSample, 0, {0., 1., 0.}, 0);
  const auto maskContactNodeVelocity = maskBendingVelocity->getVectorVelocity(
      maskContactNodeSample, 0, {0., -1., 0.}, 0);
  const NumericType maskBottomPressure =
      deformationVelocity->getPressure(maskBottomSample);
  std::cout << "Mask elasticity nodes: "
            << maskBendingVelocity->getNumberOfSolutionNodes()
            << ", contact nodes: "
            << maskBendingVelocity->getNumberOfContactNodes()
            << ", iterations: " << maskBendingVelocity->getIterations()
            << ", residual: " << maskBendingVelocity->getResidual() << '\n';
  std::cout << "Mask bottom velocity: (" << maskBottomVelocity[0] << ", "
            << maskBottomVelocity[1] << ") um/hr, top velocity: ("
            << maskTopVelocity[0] << ", " << maskTopVelocity[1]
            << ") um/hr, oxide pressure: " << maskBottomPressure << " Pa\n";
  std::cout << "Mask contact-node velocity: (" << maskContactNodeVelocity[0]
            << ", " << maskContactNodeVelocity[1] << ") um/hr\n";

  auto constrainedAmbientVelocity =
      ls::OxidationConstrainedAmbientVelocityField<NumericType, D>::New(
          deformationVelocity, maskBendingVelocity, maskInterface, -1, 1.);
  const auto maskedAmbientVelocity =
      constrainedAmbientVelocity->getVectorVelocity(maskBottomSample, 0,
                                                    {0., 1., 0.}, 0);
  const NumericType maskedAmbientScalar =
      constrainedAmbientVelocity->getScalarVelocity(maskBottomSample, 0,
                                                    {0., 1., 0.}, 0);
  std::cout << "Masked oxide/mask velocity: (" << maskedAmbientVelocity[0]
            << ", " << maskedAmbientVelocity[1]
            << ") um/hr, scalar growth: " << maskedAmbientScalar << " um/hr\n";

  ls::Advect<NumericType, D> ambientAdvection;
  ambientAdvection.insertNextLevelSet(ambientInterface);
  ambientAdvection.setVelocityField(constrainedAmbientVelocity);
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

  ls::Advect<NumericType, D> maskAdvection;
  maskAdvection.insertNextLevelSet(maskInterface);
  maskAdvection.setVelocityField(maskBendingVelocity);
  maskAdvection.setSpatialScheme(
      ls::SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER);
  maskAdvection.setTemporalScheme(
      ls::TemporalSchemeEnum::RUNGE_KUTTA_2ND_ORDER);
  maskAdvection.setAdvectionTime(advectionTime);
  maskAdvection.apply();

  writeSurface(siInterface, "locos_si_after.vtk");
  writeSurface(ambientInterface, "locos_ambient_after.vtk");
  writeSurface(maskInterface, "locos_mask_after.vtk");

  std::cout << "Wrote locos_si_initial.vtk, locos_ambient_initial.vtk, "
               "locos_mask.vtk, locos_si_after.vtk, locos_ambient_after.vtk, "
               "locos_mask_after.vtk, and locos_oxidation_diagnostics.csv\n";

  return 0;
}
