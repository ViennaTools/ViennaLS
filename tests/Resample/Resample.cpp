// Resample onto a coarser grid when the reflective extent is not a multiple
// of the new spacing (0.4 um on 1.2 nm): the new grid's edge must lie inside
// the old extent and its edge column must hold the surface's distances.  The
// new grid used to end beyond the old one; its edge column held the
// undefined sentinel (+-1.8e308) as defined values, and an oxidation run
// regridded that way aborted with a non-finite diffusion matrix.
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsResample.hpp>
#include <lsTestAsserts.hpp>

#include <cmath>
#include <cstdio>

namespace ls = viennals;
using T = double;
constexpr int D = 2;
using LevelSet = ls::SmartPointer<ls::Domain<T, D>>;

constexpr T oldDelta = 0.0006, newDelta = 0.0012;

LevelSet plane(T extent, T height, T slope) {
  double bounds[2 * D] = {-extent, extent, -0.02, 0.02};
  ls::Domain<T, D>::BoundaryType boundaries[D] = {
      ls::Domain<T, D>::BoundaryType::REFLECTIVE_BOUNDARY,
      ls::Domain<T, D>::BoundaryType::INFINITE_BOUNDARY};
  auto levelSet = LevelSet::New(bounds, boundaries, oldDelta);
  ls::MakeGeometry<T, D>(levelSet,
                         ls::Plane<T, D>::New(ls::VectorType<T, D>{0., height},
                                              ls::VectorType<T, D>{-slope, 1.}))
      .apply();
  return levelSet;
}

// largest |phi| of a defined point, and the largest deviation from the
// plane's exact signed distance (in new grid units) at the two edge columns
void check(T extent, T height, T slope) {
  auto source = plane(extent, height, slope);
  const T oldEdge = source->getGrid().getMaxGridPoint(0) * oldDelta;
  ls::Resample<T, D> resample(source, newDelta);
  resample.apply();
  auto result = resample.getResult();
  const auto &grid = result->getGrid();
  const int low = grid.getMinGridPoint(0), high = grid.getMaxGridPoint(0);
  T largest = 0., deviation = 0.;
  const T norm = std::sqrt(1. + slope * slope);
  viennahrle::ConstSparseIterator<typename ls::Domain<T, D>::DomainType> it(
      result->getDomain());
  for (; !it.isFinished(); ++it) {
    if (!it.isDefined())
      continue;
    const T value = it.getValue();
    largest = std::max(largest, std::abs(value));
    const auto index = it.getStartIndices();
    if ((index[0] == low || index[0] == high) && std::abs(value) <= 1.) {
      const T x = index[0] * newDelta, y = index[1] * newDelta;
      const T exact = (y - height - slope * x) / norm / newDelta;
      deviation = std::max(deviation, std::abs(value - exact));
    }
  }
  std::printf("extent %.3f, slope %.2f: old edge %.4f, new edge %.4f um, "
              "max |phi| %.3g, edge columns off the plane by %.3g cells\n",
              extent, slope, oldEdge, high * newDelta, largest, deviation);
  VC_TEST_ASSERT(high * newDelta <= oldEdge + 1e-12)
  VC_TEST_ASSERT(oldEdge - high * newDelta < newDelta)
  VC_TEST_ASSERT(largest < 10.)
  VC_TEST_ASSERT(deviation < 0.05)
}

int main() {
  viennacore::Logger::setLogLevel(viennacore::LogLevel::ERROR);
  for (T slope : {0., 0.3}) {
    check(0.4, 0.0024, slope);   // 0.4002 um / 1.2 nm = 333.5 cells
    check(0.675, 0.0024, slope); // 0.6756 um / 1.2 nm = 563 cells
  }
  return 0;
}
