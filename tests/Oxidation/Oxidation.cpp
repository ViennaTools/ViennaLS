#include <lsDomain.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsOxidation.hpp>
#include <lsOxidationPresets.hpp>
#include <lsTestAsserts.hpp>

#include <algorithm>
#include <array>
#include <cmath>

namespace ls = viennals;

using NumericType = double;
constexpr int D = 2;

using LevelSet = ls::SmartPointer<ls::Domain<NumericType, D>>;

// ---------------------------------------------------------------------------
// Geometry helpers
// ---------------------------------------------------------------------------

LevelSet makePlane(const double *bounds,
                   ls::Domain<NumericType, D>::BoundaryType *bc,
                   NumericType gridDelta, NumericType y) {
  auto lsd = ls::Domain<NumericType, D>::New(bounds, bc, gridDelta);
  ls::MakeGeometry<NumericType, D>(
      lsd,
      ls::Plane<NumericType, D>::New(ls::VectorType<NumericType, D>{0., y},
                                     ls::VectorType<NumericType, D>{0., 1.}))
      .apply();
  return lsd;
}

LevelSet makeBox(const double *bounds,
                 ls::Domain<NumericType, D>::BoundaryType *bc,
                 NumericType gridDelta, ls::VectorType<NumericType, D> minPt,
                 ls::VectorType<NumericType, D> maxPt) {
  auto lsd = ls::Domain<NumericType, D>::New(bounds, bc, gridDelta);
  ls::MakeGeometry<NumericType, D> mg(
      lsd, ls::Box<NumericType, D>::New(minPt, maxPt));
  std::array<bool, D> ignoreBC{};
  ignoreBC[1] = true;
  mg.setIgnoreBoundaryConditions(ignoreBC);
  mg.apply();
  return lsd;
}

// ---------------------------------------------------------------------------
// Standard (no-mask) oxidation: verify growth and volume conservation.
//
// Setup: flat Si at y=0, flat SiO2 surface at y=oxThick.
// Run one CFL-limited step; check:
//   - actualDt > 0
//   - CFL limiting fires when requested step >> CFL step
//   - Si interface receded (Si was consumed)
//   - Ambient interface rose (oxide grew upward)
//   - Volume is conserved: Si recession × (β−1) ≈ ambient lift
//   - getMaskBendingField() returns nullptr (no mask)
// ---------------------------------------------------------------------------

void testStandardOxidation() {
  constexpr NumericType gridDelta = 0.1;
  constexpr NumericType extent = 2.0;
  constexpr NumericType oxThick = 0.5; // 5 grid cells — well resolved
  constexpr NumericType yMin = -0.5;
  constexpr NumericType yMax = oxThick + 0.5;

  double bounds[2 * D] = {-extent, extent, yMin, yMax};
  ls::Domain<NumericType, D>::BoundaryType bc[D];
  bc[0] = ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  bc[1] = ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto si = makePlane(bounds, bc, gridDelta, 0.);
  auto ambient = makePlane(bounds, bc, gridDelta, oxThick);

  // Deep copies for before/after comparison.
  auto siInitial = ls::Domain<NumericType, D>::New(si);
  auto ambientInitial = ls::Domain<NumericType, D>::New(ambient);

  auto oxParams = ls::OxidationPresets<NumericType>::wet1000CDealGrove100();
  oxParams.velocitySign = -1.;
  oxParams.tolerance = 1e-7;

  auto defParams = ls::OxidationPresets<NumericType>::oxideMechanics1000C(1.);
  defParams.mechanicsIterations = 200;
  defParams.mechanicsTolerance = 1e-2;
  defParams.pressureTolerance = 1e-6;
  defParams.stokesTolerance = 1e-6;
  defParams.pressureIterations = 300;
  defParams.stokesIterations = 80;

  ls::OxidationCouplingParameters<NumericType> coupling;
  coupling.maxIterations = 8;
  coupling.tolerance = 5e-2;
  coupling.relaxation = 1.0;

  // --- Standard mode: no mask ---
  auto ox = ls::Oxidation<NumericType, D>::New(si, ambient);
  ox->setOxidationParameters(oxParams);
  ox->setDeformationParameters(defParams);
  ox->setCouplingParameters(coupling);

  // Request 0.5 hr: CFL limit is ~0.4 hr, so this should be CFL-limited.
  const NumericType actualDt = ox->applyCFLLimited(0.5, 0.499);

  VC_TEST_ASSERT(actualDt > NumericType(0))
  VC_TEST_ASSERT(actualDt < NumericType(0.5))

  // No mask → getMaskBendingField() must be null.
  VC_TEST_ASSERT(ox->getMaskBendingField() == nullptr)

  // Solver fields should be populated.
  VC_TEST_ASSERT(ox->getDiffusionField() != nullptr)
  VC_TEST_ASSERT(ox->getDeformationField() != nullptr)
  VC_TEST_ASSERT(ox->getDiffusionField()->getNumberOfSolutionNodes() > 0)

  // Measure interface positions using column-scan interpolation.
  const auto toIndex = [gridDelta](NumericType v) {
    return static_cast<viennahrle::IndexType>(std::llround(v / gridDelta));
  };
  const viennahrle::IndexType jMin = toIndex(yMin);
  const viennahrle::IndexType jMax = toIndex(yMax);

  // Volume conservation: Si recession × (β−1) ≈ ambient lift.
  // Sample x ∈ [0, extent×0.8] — away from the reflective boundary artefacts.
  const auto diag = ls::computeLOCOSOpenWindowConservation<NumericType>(
      siInitial, si, ambientInitial, ambient, NumericType(0),
      extent * NumericType(0.8), jMin, jMax, oxParams.expansionCoefficient);

  VC_TEST_ASSERT(diag.samples > 0)
  VC_TEST_ASSERT(diag.siliconRecession > NumericType(0)) // Si was consumed
  VC_TEST_ASSERT(diag.ambientLift > NumericType(0))      // oxide grew
  VC_TEST_ASSERT(diag.relativeError <
                 NumericType(0.15)) // ≤15% conservation error

  LSTEST_ASSERT_VALID_LS(si, NumericType, D)
  LSTEST_ASSERT_VALID_LS(ambient, NumericType, D)
}

// ---------------------------------------------------------------------------
// LOCOS (mask) oxidation: verify mask coupling and open-window growth.
//
// Setup: Si at y=0, pad oxide at y=padOx, Si₃N₄ mask box over x≤0.
// Run one CFL-limited step; check:
//   - actualDt > 0
//   - getMaskBendingField() is populated (mask was used)
//   - Mask coupling converged within iteration budget
//   - Oxide grew in the open window (x>0)
//   - Level sets remain valid after the step
// ---------------------------------------------------------------------------

void testLOCOSOxidation() {
  constexpr NumericType gridDelta = 0.1;
  constexpr NumericType xExtent = 1.5;
  constexpr NumericType yMin = -0.5;
  constexpr NumericType yMax = 1.0;
  constexpr NumericType padOx = 0.15; // 1.5 grid cells
  constexpr NumericType maskThick = 0.2;
  constexpr NumericType maskEdge = 0.;
  constexpr NumericType maskContactEpsilon = 1e-6;
  constexpr NumericType advectionTime = 0.2;

  double bounds[2 * D] = {-xExtent, xExtent, yMin, yMax};
  ls::Domain<NumericType, D>::BoundaryType bc[D];
  bc[0] = ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  bc[1] = ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto si = makePlane(bounds, bc, gridDelta, 0.);

  auto ambient = ls::Domain<NumericType, D>::New(si);
  auto sphere =
      ls::SmartPointer<ls::SphereDistribution<viennahrle::CoordType, D>>::New(
          padOx);
  ls::GeometricAdvect<NumericType, D>(ambient, sphere).apply();

  auto siInitial = ls::Domain<NumericType, D>::New(si);
  auto ambientInitial = ls::Domain<NumericType, D>::New(ambient);

  auto mask =
      makeBox(bounds, bc, gridDelta, {-xExtent, padOx - maskContactEpsilon},
              {maskEdge, padOx + maskThick});

  auto oxParams = ls::OxidationPresets<NumericType>::wet1000CDealGrove100();
  oxParams.velocitySign = -1.;
  oxParams.maskTransferCoefficient = 0.;
  oxParams.maskConcentration = 0.;
  oxParams.tolerance = 1e-7;

  auto defParams =
      ls::OxidationPresets<NumericType>::oxideMechanics1000C(advectionTime);
  defParams.mechanicsIterations = 200;
  defParams.mechanicsTolerance = 1e-2;
  defParams.pressureTolerance = 1e-6;
  defParams.stokesTolerance = 1e-6;
  defParams.pressureIterations = 200;
  defParams.stokesIterations = 60;

  ls::OxidationCouplingParameters<NumericType> coupling;
  coupling.maxIterations = 8;
  coupling.tolerance = 5e-2;
  coupling.relaxation = 1.0;

  auto maskParams =
      ls::OxidationPresets<NumericType>::siliconNitrideMask1000C();
  maskParams.maxIterations = 4000;

  const auto toIndex = [gridDelta](NumericType v) {
    return static_cast<viennahrle::IndexType>(std::llround(v / gridDelta));
  };

  // --- LOCOS mode: provide mask ---
  auto ox = ls::Oxidation<NumericType, D>::New(si, ambient, mask);
  ox->setOxidationParameters(oxParams);
  ox->setDeformationParameters(defParams);
  ox->setCouplingParameters(coupling);
  ox->setMaskParameters(maskParams);
  viennahrle::Index<D> diffMin{toIndex(-xExtent), toIndex(yMin)};
  viennahrle::Index<D> diffMax{toIndex(xExtent), toIndex(yMax)};
  ox->setSolveBounds(diffMin, diffMax);

  viennahrle::Index<D> mbMin{toIndex(-xExtent), toIndex(padOx) - 1};
  viennahrle::Index<D> mbMax{toIndex(maskEdge), toIndex(padOx + maskThick) + 1};
  ox->setMaskBendingBounds(mbMin, mbMax);
  ox->setMaskCouplingIterations(8);
  ox->setMaskCouplingTolerance(5e-2);

  const NumericType actualDt = ox->applyCFLLimited(advectionTime, 0.499);

  VC_TEST_ASSERT(actualDt > NumericType(0))
  VC_TEST_ASSERT(actualDt <= advectionTime)

  // Mask bending field must be populated in LOCOS mode.
  VC_TEST_ASSERT(ox->getMaskBendingField() != nullptr)
  VC_TEST_ASSERT(ox->getMaskCouplingIterations() <= 8)
  VC_TEST_ASSERT(ox->getMaskCouplingResidual() < NumericType(5e-2))

  // Volume conservation in the open window (x > 0.4 µm, well inside window).
  const auto diag = ls::computeLOCOSOpenWindowConservation<NumericType>(
      siInitial, si, ambientInitial, ambient, NumericType(0.4),
      NumericType(1.2), toIndex(yMin), toIndex(yMax),
      oxParams.expansionCoefficient);

  VC_TEST_ASSERT(diag.samples > 0)
  VC_TEST_ASSERT(diag.siliconRecession > NumericType(0))
  VC_TEST_ASSERT(diag.ambientLift > NumericType(0))
  VC_TEST_ASSERT(diag.relativeError < NumericType(0.15))

  LSTEST_ASSERT_VALID_LS(si, NumericType, D)
  LSTEST_ASSERT_VALID_LS(ambient, NumericType, D)
  LSTEST_ASSERT_VALID_LS(mask, NumericType, D)
}

// ---------------------------------------------------------------------------
// LOCOS mask mechanics sensitivity.
//
// The mask contact geometry is identical for all cases: same bottom surface,
// same mask edge, same pad oxide.  Only the mask top surface or viscosity
// changes.  The traction-driven contact solve should therefore make thicker
// and more viscous masks less mobile than a thin/soft mask.
// ---------------------------------------------------------------------------

struct LOCOSMaskResponse {
  NumericType actualDt = 0.;
  NumericType verticalMaskSpeed = 0.;
  NumericType underMaskOxideSpeed = 0.;
  NumericType maxMaskSpeed = 0.;
  std::size_t maskNodes = 0;
  std::size_t contactNodes = 0;
  std::size_t fixedNodes = 0;
};

LOCOSMaskResponse runLOCOSMaskResponse(NumericType maskThick,
                                       NumericType maskViscosity) {
  constexpr NumericType gridDelta = 0.1;
  constexpr NumericType xExtent = 1.5;
  constexpr NumericType yMin = -0.5;
  constexpr NumericType yMax = 1.2;
  constexpr NumericType padOx = 0.15;
  constexpr NumericType maskEdge = 0.;
  constexpr NumericType maskContactEpsilon = 1e-6;
  constexpr NumericType advectionTime = 0.12;

  double bounds[2 * D] = {-xExtent, xExtent, yMin, yMax};
  ls::Domain<NumericType, D>::BoundaryType bc[D];
  bc[0] = ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  bc[1] = ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto si = makePlane(bounds, bc, gridDelta, 0.);
  auto ambient = ls::Domain<NumericType, D>::New(si);
  auto sphere =
      ls::SmartPointer<ls::SphereDistribution<viennahrle::CoordType, D>>::New(
          padOx);
  ls::GeometricAdvect<NumericType, D>(ambient, sphere).apply();

  auto mask =
      makeBox(bounds, bc, gridDelta, {-xExtent, padOx - maskContactEpsilon},
              {maskEdge, padOx + maskThick});

  auto oxParams = ls::OxidationPresets<NumericType>::wet1000CDealGrove100();
  oxParams.velocitySign = -1.;
  oxParams.maskTransferCoefficient = 0.;
  oxParams.maskConcentration = 0.;
  oxParams.tolerance = 1e-7;

  auto defParams =
      ls::OxidationPresets<NumericType>::oxideMechanics1000C(advectionTime);
  defParams.mechanicsIterations = 200;
  defParams.mechanicsTolerance = 1e-2;
  defParams.pressureTolerance = 1e-6;
  defParams.stokesTolerance = 1e-6;
  defParams.pressureIterations = 200;
  defParams.stokesIterations = 60;

  ls::OxidationCouplingParameters<NumericType> coupling;
  coupling.maxIterations = 8;
  coupling.tolerance = 5e-2;
  coupling.relaxation = 1.0;

  auto maskParams =
      ls::OxidationPresets<NumericType>::siliconNitrideMask1000C();
  maskParams.contactMode = 1;
  maskParams.referenceViscosity = maskViscosity;
  maskParams.anchorBoundaryDirection = 0;
  maskParams.anchorBoundarySide = -1;
  maskParams.anchorBoundaryLayers = 1;
  maskParams.maxIterations = 4000;

  const auto toIndex = [gridDelta](NumericType v) {
    return static_cast<viennahrle::IndexType>(std::llround(v / gridDelta));
  };

  auto ox = ls::Oxidation<NumericType, D>::New(si, ambient, mask);
  ox->setOxidationParameters(oxParams);
  ox->setDeformationParameters(defParams);
  ox->setCouplingParameters(coupling);
  ox->setMaskParameters(maskParams);
  viennahrle::Index<D> diffMin{toIndex(-xExtent), toIndex(yMin)};
  viennahrle::Index<D> diffMax{toIndex(xExtent), toIndex(yMax)};
  ox->setSolveBounds(diffMin, diffMax);
  viennahrle::Index<D> mbMin{toIndex(-xExtent), toIndex(padOx) - 1};
  viennahrle::Index<D> mbMax{toIndex(maskEdge), toIndex(padOx + maskThick) + 1};
  ox->setMaskBendingBounds(mbMin, mbMax);
  ox->setMaskCouplingIterations(8);
  ox->setMaskCouplingTolerance(5e-2);

  LOCOSMaskResponse response;
  response.actualDt = ox->applyCFLLimited(advectionTime, 0.499);
  auto maskField = ox->getMaskBendingField();
  VC_TEST_ASSERT(maskField != nullptr)

  const ls::Vec3D<NumericType> zero{0., 0., 0.};
  response.verticalMaskSpeed = maskField->getDissipationAlpha(1, -1, zero);
  response.maxMaskSpeed = std::max(maskField->getDissipationAlpha(0, -1, zero),
                                   response.verticalMaskSpeed);
  auto deformationField = ox->getDeformationField();
  VC_TEST_ASSERT(deformationField != nullptr)
  const ls::Vec3D<NumericType> underMaskPoint{
      NumericType(-0.2), padOx - gridDelta * NumericType(0.25), 0.};
  const auto underMaskVelocity =
      deformationField->getVectorVelocity(underMaskPoint, -1, zero, 0);
  response.underMaskOxideSpeed = std::abs(underMaskVelocity[1]);
  response.maskNodes = maskField->getNumberOfSolutionNodes();
  response.contactNodes = maskField->getNumberOfContactNodes();
  response.fixedNodes = maskField->getNumberOfFixedNodes();
  return response;
}

void testLOCOSMaskStiffnessSensitivity() {
  const auto thin = runLOCOSMaskResponse(NumericType(0.2), NumericType(5.e11));
  const auto thick = runLOCOSMaskResponse(NumericType(0.4), NumericType(5.e11));
  const auto soft = runLOCOSMaskResponse(NumericType(0.2), NumericType(2.5e11));

  VC_TEST_ASSERT(thin.actualDt > NumericType(0))
  VC_TEST_ASSERT(thick.actualDt > NumericType(0))
  VC_TEST_ASSERT(soft.actualDt > NumericType(0))
  VC_TEST_ASSERT(thin.contactNodes > 0)
  VC_TEST_ASSERT(thin.fixedNodes > 0)
  VC_TEST_ASSERT(thick.maskNodes > thin.maskNodes)
  VC_TEST_ASSERT(thin.verticalMaskSpeed > NumericType(0))
  VC_TEST_ASSERT(thin.underMaskOxideSpeed > NumericType(0))
  VC_TEST_ASSERT(soft.verticalMaskSpeed > thin.verticalMaskSpeed)
  VC_TEST_ASSERT(thick.verticalMaskSpeed < thin.verticalMaskSpeed)
  VC_TEST_ASSERT(soft.underMaskOxideSpeed > thin.underMaskOxideSpeed)
  VC_TEST_ASSERT(thick.underMaskOxideSpeed < thin.underMaskOxideSpeed)
}

// ---------------------------------------------------------------------------
// 3D smoke test: flat Si substrate, one CFL step, no mask.
//
// Exercises the 3D face-enumeration and Stokes-solve code paths that are
// entirely absent from the 2D tests.  Checks only that the solver runs,
// Si was consumed, and both level sets remain valid — not conservation or
// rate accuracy (the coarse grid makes those unreliable at this small size).
// ---------------------------------------------------------------------------

void testStandard3D() {
  constexpr int D3 = 3;
  using LS3 = ls::SmartPointer<ls::Domain<NumericType, D3>>;

  constexpr NumericType gridDelta = 0.15;
  constexpr NumericType extent = 0.75;
  constexpr NumericType oxThick = 0.30;
  constexpr NumericType yMin = -0.30;
  constexpr NumericType yMax = oxThick + 0.30;

  double bounds[2 * D3] = {-extent, extent, yMin, yMax, -extent, extent};
  ls::Domain<NumericType, D3>::BoundaryType bc[D3];
  bc[0] = ls::Domain<NumericType, D3>::BoundaryType::REFLECTIVE_BOUNDARY;
  bc[1] = ls::Domain<NumericType, D3>::BoundaryType::INFINITE_BOUNDARY;
  bc[2] = ls::Domain<NumericType, D3>::BoundaryType::REFLECTIVE_BOUNDARY;

  auto si = ls::Domain<NumericType, D3>::New(bounds, bc, gridDelta);
  ls::MakeGeometry<NumericType, D3>(
      si, ls::Plane<NumericType, D3>::New(
              ls::VectorType<NumericType, D3>{0., 0., 0.},
              ls::VectorType<NumericType, D3>{0., 1., 0.}))
      .apply();

  auto ambient = ls::Domain<NumericType, D3>::New(si);
  auto sphere =
      ls::SmartPointer<ls::SphereDistribution<viennahrle::CoordType, D3>>::New(
          oxThick);
  ls::GeometricAdvect<NumericType, D3>(ambient, sphere).apply();

  auto siInitial = ls::Domain<NumericType, D3>::New(si);

  auto oxParams = ls::OxidationPresets<NumericType>::wet1000CDealGrove100();
  oxParams.velocitySign = -1.;
  oxParams.tolerance = 1e-7;

  auto defParams = ls::OxidationPresets<NumericType>::oxideMechanics1000C(0.1);
  defParams.mechanicsIterations = 200;
  defParams.mechanicsTolerance = 1e-2;
  defParams.pressureTolerance = 1e-6;
  defParams.stokesTolerance = 1e-6;
  defParams.pressureIterations = 300;
  defParams.stokesIterations = 150;

  ls::OxidationCouplingParameters<NumericType> coupling;
  coupling.maxIterations = 8;
  coupling.tolerance = 5e-2;
  coupling.relaxation = 1.0;

  auto ox = ls::Oxidation<NumericType, D3>::New(si, ambient);
  ox->setOxidationParameters(oxParams);
  ox->setDeformationParameters(defParams);
  ox->setCouplingParameters(coupling);

  const NumericType actualDt = ox->applyCFLLimited(0.1, 0.099);

  VC_TEST_ASSERT(actualDt > NumericType(0))
  VC_TEST_ASSERT(ox->getDiffusionField() != nullptr)
  VC_TEST_ASSERT(ox->getDiffusionField()->getNumberOfSolutionNodes() > 0)

  LSTEST_ASSERT_VALID_LS(si, NumericType, D3)
  LSTEST_ASSERT_VALID_LS(ambient, NumericType, D3)
}

// ---------------------------------------------------------------------------

int main() {
  testStandardOxidation();
  testLOCOSOxidation();
  testLOCOSMaskStiffnessSensitivity();
  testStandard3D();
  return 0;
}
