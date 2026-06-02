#include <lsDomain.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsOxidation.hpp>
#include <lsOxidationPresets.hpp>
#include <lsTestAsserts.hpp>

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
      lsd, ls::Plane<NumericType, D>::New(
               ls::VectorType<NumericType, D>{0., y},
               ls::VectorType<NumericType, D>{0., 1.}))
      .apply();
  return lsd;
}

LevelSet makeBox(const double *bounds,
                 ls::Domain<NumericType, D>::BoundaryType *bc,
                 NumericType gridDelta,
                 ls::VectorType<NumericType, D> minPt,
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
  constexpr NumericType gridDelta   = 0.1;
  constexpr NumericType extent      = 2.0;
  constexpr NumericType oxThick     = 0.5;  // 5 grid cells — well resolved
  constexpr NumericType yMin        = -0.5;
  constexpr NumericType yMax        = oxThick + 0.5;

  double bounds[2 * D] = {-extent, extent, yMin, yMax};
  ls::Domain<NumericType, D>::BoundaryType bc[D];
  bc[0] = ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  bc[1] = ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto si      = makePlane(bounds, bc, gridDelta, 0.);
  auto ambient = makePlane(bounds, bc, gridDelta, oxThick);

  // Deep copies for before/after comparison.
  auto siInitial      = ls::Domain<NumericType, D>::New(si);
  auto ambientInitial = ls::Domain<NumericType, D>::New(ambient);

  auto oxParams =
      ls::OxidationPresets<NumericType>::wet1000CDealGrove100();
  oxParams.velocitySign = -1.;
  oxParams.tolerance    = 1e-7;

  auto defParams =
      ls::OxidationPresets<NumericType>::oxideMechanics1000C(1.);
  defParams.mechanicsIterations = 10;
  defParams.pressureIterations  = 300;
  defParams.stokesIterations    = 80;

  ls::OxidationCouplingParameters<NumericType> coupling;
  coupling.maxIterations = 6;
  coupling.tolerance     = 1e-4;

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
      siInitial, si, ambientInitial, ambient,
      NumericType(0), extent * NumericType(0.8),
      jMin, jMax, oxParams.expansionCoefficient);

  VC_TEST_ASSERT(diag.samples > 0)
  VC_TEST_ASSERT(diag.siliconRecession > NumericType(0))  // Si was consumed
  VC_TEST_ASSERT(diag.ambientLift      > NumericType(0))  // oxide grew
  VC_TEST_ASSERT(diag.relativeError < NumericType(0.15))  // ≤15% conservation error

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
  constexpr NumericType gridDelta          = 0.1;
  constexpr NumericType xExtent            = 1.5;
  constexpr NumericType yMin               = -0.5;
  constexpr NumericType yMax               = 1.0;
  constexpr NumericType padOx              = 0.15;  // 1.5 grid cells
  constexpr NumericType maskThick          = 0.2;
  constexpr NumericType maskEdge           = 0.;
  constexpr NumericType maskContactEpsilon = 1e-6;
  constexpr NumericType advectionTime      = 0.2;

  double bounds[2 * D] = {-xExtent, xExtent, yMin, yMax};
  ls::Domain<NumericType, D>::BoundaryType bc[D];
  bc[0] = ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  bc[1] = ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto si = makePlane(bounds, bc, gridDelta, 0.);

  auto ambient = ls::Domain<NumericType, D>::New(si);
  auto sphere  = ls::SmartPointer<ls::SphereDistribution<
      viennahrle::CoordType, D>>::New(padOx);
  ls::GeometricAdvect<NumericType, D>(ambient, sphere).apply();

  auto siInitial      = ls::Domain<NumericType, D>::New(si);
  auto ambientInitial = ls::Domain<NumericType, D>::New(ambient);

  auto mask = makeBox(bounds, bc, gridDelta,
                      {-xExtent, padOx - maskContactEpsilon},
                      {maskEdge,  padOx + maskThick});

  auto oxParams =
      ls::OxidationPresets<NumericType>::wet1000CDealGrove100();
  oxParams.velocitySign              = -1.;
  oxParams.maskTransferCoefficient   = 0.;
  oxParams.maskConcentration         = 0.;
  oxParams.tolerance                 = 1e-7;

  auto defParams =
      ls::OxidationPresets<NumericType>::oxideMechanics1000C(advectionTime);
  defParams.pressureIterations = 200;
  defParams.stokesIterations   = 60;

  ls::OxidationCouplingParameters<NumericType> coupling;
  coupling.maxIterations = 6;
  coupling.tolerance     = 1e-6;

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
  viennahrle::Index<D> diffMax{toIndex( xExtent), toIndex(yMax)};
  ox->setSolveBounds(diffMin, diffMax);

  viennahrle::Index<D> mbMin{toIndex(-xExtent), toIndex(padOx) - 1};
  viennahrle::Index<D> mbMax{toIndex( maskEdge), toIndex(padOx + maskThick) + 1};
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
      siInitial, si, ambientInitial, ambient,
      NumericType(0.4), NumericType(1.2),
      toIndex(yMin), toIndex(yMax),
      oxParams.expansionCoefficient);

  VC_TEST_ASSERT(diag.samples > 0)
  VC_TEST_ASSERT(diag.siliconRecession > NumericType(0))
  VC_TEST_ASSERT(diag.ambientLift      > NumericType(0))
  VC_TEST_ASSERT(diag.relativeError < NumericType(0.15))

  LSTEST_ASSERT_VALID_LS(si, NumericType, D)
  LSTEST_ASSERT_VALID_LS(ambient, NumericType, D)
  LSTEST_ASSERT_VALID_LS(mask, NumericType, D)
}

// ---------------------------------------------------------------------------

int main() {
  testStandardOxidation();
  testLOCOSOxidation();
  return 0;
}
