#include <cmath>
#include <iostream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMaterialStack.hpp>
#include <vcTestAsserts.hpp>

/**
  The material of points around the corners of a Euclidean trench whose
  corners lie on grid nodes, as ViennaPS's MakeTrench builds it (a plane less
  a box, grid spacing 2.5, 70 wide and 70 deep): every point 0.01 and 0.1
  grid spacings beside a face near a corner, on both sides of it, takes the
  material of its side. At a concave corner on a grid node, the node holds
  the corner with the normal of one face and the node beyond the corner
  along the other face holds it with the other's normal; the two faces of
  the cell then pass through one point.
  \example MaterialStack.cpp
*/

namespace ls = viennals;

int main() {
  constexpr int D = 2;
  using T = double;
  omp_set_num_threads(1);

  const T gridDelta = 2.5;
  T bounds[2 * D] = {-55., 55., -100., 100.};
  ls::BoundaryConditionEnum boundaryCons[D] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
  auto newLevelSet = [&]() {
    auto domain = ls::SmartPointer<ls::Domain<T, D>>::New(bounds, boundaryCons,
                                                          gridDelta);
    domain->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
    domain->setLevelSetWidth(12);
    return domain;
  };
  auto substrate = newLevelSet();
  {
    T origin[D] = {0., 0.};
    T normal[D] = {0., 1.};
    ls::MakeGeometry<T, D>(
        substrate, ls::SmartPointer<ls::Plane<T, D>>::New(origin, normal))
        .apply();
  }
  auto trench = newLevelSet();
  {
    T minCorner[D] = {-35., -70.};
    T maxCorner[D] = {35., 10.};
    ls::MakeGeometry<T, D>(
        trench, ls::SmartPointer<ls::Box<T, D>>::New(minCorner, maxCorner))
        .apply();
  }
  ls::BooleanOperation<T, D>(substrate, trench,
                             ls::BooleanOperationEnum::RELATIVE_COMPLEMENT)
      .apply();

  ls::MaterialStack<T, D> stack({substrate});
  // inside the substrate: below the top beside the trench, or below its floor
  auto solid = [](T x, T y) {
    return (y < 0. && std::abs(x) > 35.) || y < -70.;
  };
  int wrong = 0, probes = 0;
  // the corners, and the directions along their faces
  const T corners[4][2] = {{-35., -70.}, {35., -70.}, {-35., 0.}, {35., 0.}};
  for (const auto &c : corners)
    for (T along : {0.01, 0.1, 0.5})
      for (T off : {0.01, 0.1})
        for (int axis = 0; axis < D; ++axis)
          for (int s1 : {-1, 1})
            for (int s2 : {-1, 1}) {
              // a point beside the face through the corner normal to the
              // other axis, off it by off and along it by along (um)
              T p[D] = {c[0], c[1]};
              p[axis] += s1 * off * gridDelta;
              p[1 - axis] += s2 * along * gridDelta;
              const bool inside = solid(p[0], p[1]);
              const int expected = inside ? 0 : 1;
              const int found = stack.material(
                  ls::Vec3D<T>{p[0] / gridDelta, p[1] / gridDelta, T(0)});
              ++probes;
              if (found != expected) {
                ++wrong;
                std::cout << "point (" << p[0] << ", " << p[1] << "): material "
                          << found << ", value "
                          << stack.value(0,
                                         ls::Vec3D<T>{p[0] / gridDelta,
                                                      p[1] / gridDelta, T(0)})
                          << " grid spacings, expected " << expected
                          << std::endl;
              }
            }
  std::cout << "points around the trench corners: " << wrong << " of " << probes
            << " in the wrong material" << std::endl;
  VC_TEST_ASSERT(wrong == 0);
  return 0;
}
