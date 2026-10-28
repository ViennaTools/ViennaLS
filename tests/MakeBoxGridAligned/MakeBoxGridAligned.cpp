#include <iostream>

#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMakeGeometry.hpp>
#include <lsTestAsserts.hpp>

/**
  A box whose face lies on a grid line, where the face coordinate divided by
  the grid spacing comes out one ulp off the integer (-0.2975 / 0.0025 =
  -118.99999999999999).  The box overhangs the reflective boundary, so the
  rows inside it have no defined point on their left.  FromSurfaceMesh used
  to drop the band point inside the face; the face point (value 0) then set
  the sign of the next rows' undefined runs and the box came out hollow
  (undefined neighbours of opposite sign), and Expand read past the point
  data.  A LOCOS nitride mask on a 2.5 nm grid crashed this way.
*/

namespace ls = viennals;

int main() {
  constexpr int D = 2;
  using T = double;
  const T gridDelta = 0.0025;
  double bounds[2 * D] = {-1.2, 1.2, -1., 1.5};
  ls::BoundaryConditionEnum bcs[D] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};

  for (const T bottom : {T(-1e-6), T(0.005 - 1e-6)}) {
    auto box = ls::SmartPointer<ls::Domain<T, D>>::New(bounds, bcs, gridDelta);
    ls::MakeGeometry<T, D> geometry(
        box, ls::SmartPointer<ls::Box<T, D>>::New(
                 ls::VectorType<T, D>{-2.4, bottom},
                 ls::VectorType<T, D>{-0.2975, bottom + 0.14}));
    geometry.setIgnoreBoundaryConditions(std::array<bool, D>{false, true});
    geometry.apply();
    LSTEST_ASSERT_VALID_LS(box, T, D)

    ls::Expand<T, D>(box, 9).apply();
    LSTEST_ASSERT_VALID_LS(box, T, D)
  }
  return 0;
}
