#include <lsDomain.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsOxidationModel.hpp>
#include <lsTestAsserts.hpp>

#include <cmath>

namespace ls = viennals;

int main() {
  constexpr int D = 2;
  constexpr double gridDelta = 1.;
  constexpr double extent = 8.;
  constexpr double oxideThickness = 6.;

  double bounds[2 * D] = {-extent, extent, -1., oxideThickness + 1.};
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = ls::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[1] = ls::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto reactionInterface =
      ls::Domain<double, D>::New(bounds, boundaryCons, gridDelta);
  auto ambientInterface =
      ls::Domain<double, D>::New(bounds, boundaryCons, gridDelta);

  const ls::VectorType<double, D> reactionOrigin{0., 0.};
  const ls::VectorType<double, D> ambientOrigin{0., oxideThickness};
  const ls::VectorType<double, D> normal{0., 1.};

  ls::MakeGeometry<double, D>(
      reactionInterface,
      ls::Plane<double, D>::New(reactionOrigin, normal))
      .apply();
  ls::MakeGeometry<double, D>(
      ambientInterface,
      ls::Plane<double, D>::New(ambientOrigin, normal))
      .apply();

  ls::OxidationParameters<double> parameters;
  parameters.diffusionCoefficient = 1.;
  parameters.reactionRate = 1.;
  parameters.transferCoefficient = 1.;
  parameters.equilibriumConcentration = 1.;
  parameters.expansionCoefficient = 2.27;
  parameters.maxIterations = 20000;
  parameters.tolerance = 1e-10;

  auto oxidation = ls::OxidationDiffusionVelocityField<double, D>::New(
      reactionInterface, ambientInterface, parameters);
  viennahrle::Index<D> minIndex{-8, -1};
  viennahrle::Index<D> maxIndex{8, 7};
  oxidation->setSolveBounds(minIndex, maxIndex);
  oxidation->apply();

  VC_TEST_ASSERT(oxidation->getNumberOfSolutionNodes() > 0)
  VC_TEST_ASSERT(oxidation->getResidual() < 1e-6)

  ls::Vec3D<double> nearReaction{0., 1., 0.};
  ls::Vec3D<double> nearAmbient{0., oxideThickness - 1., 0.};

  const double reactionConcentration =
      oxidation->getConcentration(nearReaction);
  const double ambientConcentration = oxidation->getConcentration(nearAmbient);

  VC_TEST_ASSERT(reactionConcentration > 0.)
  VC_TEST_ASSERT(ambientConcentration > reactionConcentration)
  VC_TEST_ASSERT(ambientConcentration < parameters.equilibriumConcentration)

  const double velocity =
      oxidation->getScalarVelocity(nearReaction, 0, {0., 1., 0.}, 0);
  VC_TEST_ASSERT(velocity > 0.)

  ls::OxidationDeformationParameters<double> deformationParameters;
  deformationParameters.viscosity = 1e3;
  deformationParameters.mechanicsIterations = 10;
  deformationParameters.pressureIterations = 500;
  deformationParameters.stokesIterations = 100;
  deformationParameters.tolerance = 1e-8;

  auto deformation =
      ls::OxidationDeformationVelocityField<double, D>::New(
          reactionInterface, ambientInterface, oxidation, parameters,
          deformationParameters);
  deformation->setSolveBounds(minIndex, maxIndex);
  deformation->apply();

  VC_TEST_ASSERT(deformation->getNumberOfSolutionNodes() > 0)
  VC_TEST_ASSERT(std::isfinite(deformation->getResidual()))
  VC_TEST_ASSERT(deformation->getResidual() < 2.)

  // The (gamma-1)/gamma expansion ratio must be checked at the same y location
  // as the deformation boundary velocity: the reaction interface node at y=0.
  // Querying at y=1 vs y=0 uses different concentrations and breaks the ratio.
  ls::Vec3D<double> reactionBoundary{0., 0., 0.};
  const double siliconSpeed =
      std::abs(oxidation->getScalarVelocity(reactionBoundary, 0, {0., 1., 0.}, 0));
  const double ambientSpeed =
      std::abs(deformation->getVectorVelocity(reactionBoundary, 0, {0., 1., 0.}, 0)[1]);
  const double oxidePressure = deformation->getPressure(nearReaction);
  const auto stressTensor = deformation->getStressTensor(nearReaction);
  const double vonMisesStress = deformation->getVonMisesStress(nearReaction);
  VC_TEST_ASSERT(std::isfinite(oxidePressure))
  VC_TEST_ASSERT(std::abs(oxidePressure) > 0.)
  VC_TEST_ASSERT(std::isfinite(stressTensor[0]))
  VC_TEST_ASSERT(std::isfinite(stressTensor[4]))
  VC_TEST_ASSERT(std::isfinite(vonMisesStress))
  VC_TEST_ASSERT(vonMisesStress >= 0.)

  const double ambientFraction = ambientSpeed / (ambientSpeed + siliconSpeed);
  const double siliconFraction = siliconSpeed / (ambientSpeed + siliconSpeed);
  const double expectedAmbientFraction =
      (parameters.expansionCoefficient - 1.) / parameters.expansionCoefficient;
  const double expectedSiliconFraction = 1. / parameters.expansionCoefficient;

  VC_TEST_ASSERT(std::abs(ambientFraction - expectedAmbientFraction) < 0.03)
  VC_TEST_ASSERT(std::abs(siliconFraction - expectedSiliconFraction) < 0.03)

  auto maskInterface = ls::Domain<double, D>::New(bounds, boundaryCons, gridDelta);
  double maskMin[D] = {2., -1.};
  double maskMax[D] = {extent, oxideThickness + 1.};
  ls::MakeGeometry<double, D>(
      maskInterface, ls::Box<double, D>::New(maskMin, maskMax))
      .apply();

  auto maskedOxidation = ls::OxidationDiffusionVelocityField<double, D>::New(
      reactionInterface, ambientInterface, parameters);
  maskedOxidation->setMaskInterface(maskInterface, -1);
  maskedOxidation->setSolveBounds(minIndex, maxIndex);
  maskedOxidation->apply();

  VC_TEST_ASSERT(maskedOxidation->getNumberOfSolutionNodes() <
                 oxidation->getNumberOfSolutionNodes())
  const ls::Vec3D<double> nearMaskedAmbient{-1., oxideThickness - 1., 0.};
  VC_TEST_ASSERT(maskedOxidation->getConcentration(nearMaskedAmbient) > 0.)

  ls::OxidationDeformationParameters<double> maskedDeformationParameters =
      deformationParameters;
  maskedDeformationParameters.maskNormalStiffness = 1e6;
  auto maskedDeformation =
      ls::OxidationDeformationVelocityField<double, D>::New(
          reactionInterface, ambientInterface, maskedOxidation, parameters,
          maskedDeformationParameters);
  maskedDeformation->setMaskInterface(maskInterface, -1);
  maskedDeformation->setSolveBounds(minIndex, maxIndex);
  maskedDeformation->apply();
  VC_TEST_ASSERT(maskedDeformation->getNumberOfSolutionNodes() <
                 deformation->getNumberOfSolutionNodes())
  VC_TEST_ASSERT(std::isfinite(maskedDeformation->getResidual()))

  ls::OxidationParameters<double> stressParameters = parameters;
  stressParameters.stressCouplingCoefficient = 1e-9;
  auto stressOxidation = ls::OxidationDiffusionVelocityField<double, D>::New(
      reactionInterface, ambientInterface, stressParameters);
  stressOxidation->setSolveBounds(minIndex, maxIndex);
  stressOxidation->apply();
  const double unstressedRate =
      stressOxidation->getEffectiveReactionRate(nearReaction);
  stressOxidation->setPressure(nearReaction, 1e9);
  const double compressedRate =
      stressOxidation->getEffectiveReactionRate(nearReaction);
  VC_TEST_ASSERT(compressedRate < unstressedRate)

  ls::OxidationCouplingParameters<double> couplingParameters;
  couplingParameters.maxIterations = 2;
  couplingParameters.tolerance = 1e12;
  auto coupled = ls::OxidationCoupledModel<double, D>::New(
      oxidation, deformation, couplingParameters);
  coupled->setSolveBounds(minIndex, maxIndex);
  coupled->apply();
  VC_TEST_ASSERT(coupled->getIterations() < couplingParameters.maxIterations)

  LSTEST_ASSERT_VALID_LS(reactionInterface, double, D)
  LSTEST_ASSERT_VALID_LS(ambientInterface, double, D)

  return 0;
}
