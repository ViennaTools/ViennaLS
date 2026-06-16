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
  boundaryCons[1] = ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

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

  auto oxParams = ls::OxidationPresets::wet1000CDealGrove100();
  oxParams.velocitySign = -1.;
  oxParams.maskTransferCoefficient = 0.;
  oxParams.maskConcentration = 0.;
  oxParams.tolerance = 1.e-7;

  auto defParams = ls::OxidationPresets::oxideMechanics1000C(advectionTime);
  defParams.mechanicsIterations = 200;
  defParams.mechanicsTolerance = 1e-2;
  defParams.pressureTolerance = 1e-6;
  defParams.stokesTolerance = 1e-6;
  defParams.pressureIterations = 200;
  defParams.stokesIterations = 60;

  ls::OxidationCouplingParameters couplingParams;
  couplingParams.maxIterations = 8;
  couplingParams.tolerance = 5e-2;
  couplingParams.relaxation = 1.0;

  auto maskParams = ls::OxidationPresets::siliconNitrideMask1000C();
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

  auto locos = ls::Oxidation<NumericType, D>::New(siInterface, ambientInterface,
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
  testLOCOSInterfaceConvergence();
  return 0;
}
