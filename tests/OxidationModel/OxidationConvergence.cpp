#include <lsDomain.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsGeometries.hpp>
#include <lsLOCOSOxidation.hpp>
#include <lsMakeGeometry.hpp>
#include <lsOxidationMaterials.hpp>
#include <lsOxidationModel.hpp>
#include <lsTestAsserts.hpp>

#include <array>
#include <cmath>

namespace ls = viennals;

using NumericType = double;
constexpr int D = 2;
using LevelSet = ls::SmartPointer<ls::Domain<NumericType, D>>;

namespace {

LevelSet makePlane(const double *bounds,
                   ls::Domain<NumericType, D>::BoundaryType *boundaryCons,
                   NumericType gridDelta, NumericType y) {
  auto levelSet =
      ls::Domain<NumericType, D>::New(bounds, boundaryCons, gridDelta);
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

NumericType flatAmbientFraction(NumericType minMechanicsBoundaryDistance,
                                NumericType tolerance) {
  constexpr NumericType gridDelta = 0.5;
  constexpr NumericType extent = 4.;
  constexpr NumericType oxideThickness = 3.;
  double bounds[2 * D] = {-extent, extent, -1., oxideThickness + 1.};
  ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] =
      ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[1] =
      ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto reactionInterface = makePlane(bounds, boundaryCons, gridDelta, 0.);
  auto ambientInterface =
      makePlane(bounds, boundaryCons, gridDelta, oxideThickness);

  ls::OxidationParameters<NumericType> oxParams;
  oxParams.diffusionCoefficient = 1.;
  oxParams.reactionRate = 1.;
  oxParams.transferCoefficient = 1.;
  oxParams.equilibriumConcentration = 1.;
  oxParams.expansionCoefficient = 2.27;
  oxParams.maxIterations = 20000;
  oxParams.tolerance = tolerance;

  viennahrle::Index<D> minIndex{-8, -1};
  viennahrle::Index<D> maxIndex{8, 4};

  auto diffusion = ls::OxidationDiffusionVelocityField<NumericType, D>::New(
      reactionInterface, ambientInterface, oxParams);
  diffusion->setSolveBounds(minIndex, maxIndex);
  diffusion->apply();

  ls::OxidationDeformationParameters<NumericType> defParams;
  defParams.viscosity = 1.e3;
  defParams.bulkModulus = 1.;
  defParams.minMechanicsBoundaryDistance = minMechanicsBoundaryDistance;
  defParams.mechanicsIterations = 10;
  defParams.pressureIterations = 300;
  defParams.stokesIterations = 80;
  defParams.tolerance = tolerance;
  defParams.pressureTolerance = tolerance;
  defParams.stokesTolerance = tolerance;

  auto deformation =
      ls::OxidationDeformationVelocityField<NumericType, D>::New(
          reactionInterface, ambientInterface, diffusion, oxParams, defParams);
  deformation->setSolveBounds(minIndex, maxIndex);
  deformation->apply();

  const ls::Vec3D<NumericType> reactionBoundary{0., 0., 0.};
  const NumericType siliconSpeed = std::abs(diffusion->getScalarVelocity(
      reactionBoundary, 0, {0., 1., 0.}, 0));
  const NumericType ambientSpeed = std::abs(deformation->getVectorVelocity(
      reactionBoundary, 0, {0., 1., 0.}, 0)[1]);
  return ambientSpeed / (ambientSpeed + siliconSpeed);
}

void testMechanicsConvergence() {
  const NumericType reference = flatAmbientFraction(1.e-6, 1.e-9);
  const NumericType regularized = flatAmbientFraction(1.e-3, 1.e-9);
  const NumericType loose = flatAmbientFraction(1.e-6, 1.e-7);
  const NumericType expected = (2.27 - 1.) / 2.27;

  VC_TEST_ASSERT(std::abs(reference - expected) < 0.03)
  VC_TEST_ASSERT(std::abs(regularized - reference) < 0.01)
  VC_TEST_ASSERT(std::abs(loose - reference) < 0.01)
}

void testLOCOSInterfaceConvergence() {
  constexpr NumericType gridDelta = 0.1;
  constexpr NumericType xExtent = 1.5;
  constexpr NumericType yMin = -0.5;
  constexpr NumericType yMax = 1.0;
  constexpr NumericType padOxideThickness = 0.15;
  constexpr NumericType maskThickness = 0.2;
  constexpr NumericType maskEdge = 0.;
  constexpr NumericType advectionTime = 0.2;
  constexpr NumericType maskContactEpsilon = 1.e-6;

  double bounds[2 * D] = {-xExtent, xExtent, yMin, yMax};
  ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] =
      ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[1] =
      ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  auto siInterface = makePlane(bounds, boundaryCons, gridDelta, 0.);
  auto ambientInterface = ls::Domain<NumericType, D>::New(siInterface);
  auto initialOxide =
      ls::SmartPointer<ls::SphereDistribution<viennahrle::CoordType, D>>::New(
          padOxideThickness);
  ls::GeometricAdvect<NumericType, D>(ambientInterface, initialOxide).apply();
  auto siInitial = ls::Domain<NumericType, D>::New(siInterface);
  auto ambientInitial = ls::Domain<NumericType, D>::New(ambientInterface);

  auto maskInterface =
      makeMask(bounds, boundaryCons, gridDelta, -xExtent, maskEdge,
               padOxideThickness - maskContactEpsilon,
               padOxideThickness + maskThickness);

  auto oxParams =
      ls::OxidationProcessPresets<NumericType>::wet1000CDealGrove100();
  oxParams.velocitySign = -1.;
  oxParams.maskTransferCoefficient = 0.;
  oxParams.maskConcentration = 0.;
  oxParams.tolerance = 1.e-7;

  auto defParams =
      ls::OxidationProcessPresets<NumericType>::oxideMechanics1000C(
          advectionTime);
  defParams.pressureIterations = 200;
  defParams.stokesIterations = 60;

  ls::OxidationCouplingParameters<NumericType> couplingParams;
  couplingParams.maxIterations = 6;
  couplingParams.tolerance = 1.e-6;

  auto maskParams =
      ls::OxidationProcessPresets<NumericType>::siliconNitrideMask1000C();
  maskParams.maxIterations = 4000;

  const auto toIndex = [](NumericType value) {
    return static_cast<viennahrle::IndexType>(std::llround(value / gridDelta));
  };
  viennahrle::Index<D> diffMinIndex{toIndex(-xExtent), toIndex(yMin)};
  viennahrle::Index<D> diffMaxIndex{toIndex(xExtent), toIndex(yMax)};
  viennahrle::Index<D> maskMinIndex{toIndex(-xExtent),
                                    toIndex(padOxideThickness) - 1};
  viennahrle::Index<D> maskMaxIndex{
      toIndex(maskEdge), toIndex(padOxideThickness + maskThickness) + 1};

  auto locos =
      ls::LOCOSOxidation<NumericType, D>::New(siInterface, ambientInterface,
                                              maskInterface);
  locos->setOxidationParameters(oxParams);
  locos->setDeformationParameters(defParams);
  locos->setCouplingParameters(couplingParams);
  locos->setMaskParameters(maskParams);
  locos->setSolveBounds(diffMinIndex, diffMaxIndex);
  locos->setMaskBendingBounds(maskMinIndex, maskMaxIndex);
  locos->setMaskCouplingIterations(8);
  locos->setMaskCouplingTolerance(5.e-2);
  locos->apply(advectionTime);

  VC_TEST_ASSERT(locos->getMaskCouplingIterations() <= 8)
  VC_TEST_ASSERT(locos->getMaskCouplingResidual() < 5.e-2)

  const auto conservation = ls::computeLOCOSOpenWindowConservation<NumericType>(
      siInitial, siInterface, ambientInitial, ambientInterface, 0.4, 1.2,
      toIndex(yMin), toIndex(yMax), oxParams.expansionCoefficient);
  VC_TEST_ASSERT(conservation.samples > 0)
  VC_TEST_ASSERT(conservation.relativeError < 0.1)
}

} // namespace

int main() {
  testMechanicsConvergence();
  testLOCOSInterfaceConvergence();
  return 0;
}
