#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <hrleSparseIterator.hpp>

#include <lsAdvect.hpp>
#include <lsCheck.hpp>
#include <lsCompareSurfaceDistance.hpp>
#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMesh.hpp>
#include <lsToSurfaceMesh.hpp>

/**
  Advect on Euclidean level sets (setEuclideanRenormalization), values only:
  after every step the points next to the surface keep their values and the
  rest of the band takes exact distances to the surface. Every case runs once
  on a Manhattan and once on a Euclidean level set, from the same surface
  with the same scheme, and measures the surface each delivers (marching
  squares or cubes for both) against the exact one with
  CompareSurfaceDistance.
  The Euclidean run lies at least as close to the exact surface as the
  Manhattan run, in the Chamfer and in the maximum distance:
  - a substrate etched directionally beside a mask that does not advect
    (exact: a trench whose walls go straight down from the mask's edges),
  - a rotated rectangle etched from above by a vector velocity (exact: the
    rectangle intersected with its translate by the etched depth),
  - isotropic growth on a substrate beside a mask of zero rate (exact: the
    mask unchanged, the substrate surface lifted by the grown thickness).
  A sphere grown and a sphere etched isotropically (exact: the spheres of the
  radius plus and minus the distance), where both normalizations share the
  first-order bias of the scheme, lie within half a grid spacing of the exact
  sphere, the resolution of the grid, with both.
  \example AdvectEuclidean.cpp
*/

namespace ls = viennals;
using Polygon = std::vector<std::array<double, 2>>;

// the substrate (material 1) etched at the upward component of its normal
class Etch : public ls::VelocityField<double> {
public:
  double getScalarVelocity(const ls::Vec3D<double> &, int material,
                           const ls::Vec3D<double> &normal,
                           unsigned long) override {
    return material == 1 ? -std::max(normal[1], 0.) : 0.;
  }
};

// etching from above: the surface moves down at unit speed where it faces up
class FromAbove : public ls::VelocityField<double> {
public:
  ls::Vec3D<double> getVectorVelocity(const ls::Vec3D<double> &, int,
                                      const ls::Vec3D<double> &normal,
                                      unsigned long) override {
    return normal[1] > 0 ? ls::Vec3D<double>{0., -1., 0.}
                         : ls::Vec3D<double>{0., 0., 0.};
  }
};

// growth at unit speed on the substrate (material 1), none on the mask
class OnSubstrate : public ls::VelocityField<double> {
public:
  double getScalarVelocity(const ls::Vec3D<double> &, int material,
                           const ls::Vec3D<double> &, unsigned long) override {
    return material == 1 ? 1. : 0.;
  }
};

// isotropic motion at a constant speed
class Constant : public ls::VelocityField<double> {
  double speed;

public:
  explicit Constant(double s) : speed(s) {}
  double getScalarVelocity(const ls::Vec3D<double> &, int,
                           const ls::Vec3D<double> &, unsigned long) override {
    return speed;
  }
};

// a polyline as a 2D surface mesh (closed: the last point joins the first)
ls::SmartPointer<ls::Mesh<>> polyline(const Polygon &P, bool closed) {
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  for (const auto &p : P)
    mesh->insertNextNode({p[0], p[1], 0.});
  const unsigned n = unsigned(P.size());
  for (unsigned k = 0; k + 1 < n; ++k)
    mesh->insertNextLine({k, k + 1});
  if (closed)
    mesh->insertNextLine({n - 1, 0});
  return mesh;
}

// the corners of the rectangle of half extents hx, hy about (cx, cy), turned
// by angle, counter-clockwise
Polygon rectangle(double cx, double cy, double hx, double hy, double angle) {
  const double c = std::cos(angle), s = std::sin(angle);
  Polygon P;
  for (const auto &q :
       {std::array<double, 2>{-hx, -hy}, {hx, -hy}, {hx, hy}, {-hx, hy}})
    P.push_back({cx + c * q[0] - s * q[1], cy + s * q[0] + c * q[1]});
  return P;
}

// the part of convex polygon P inside convex polygon Q, both
// counter-clockwise (Sutherland-Hodgman)
Polygon clip(Polygon P, const Polygon &Q) {
  for (std::size_t k = 0; k < Q.size() && !P.empty(); ++k) {
    const auto &a = Q[k], &b = Q[(k + 1) % Q.size()];
    auto side = [&](const std::array<double, 2> &p) {
      return (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]);
    };
    Polygon out;
    for (std::size_t i = 0; i < P.size(); ++i) {
      const auto &p = P[i], &q = P[(i + 1) % P.size()];
      const double sp = side(p), sq = side(q);
      if (sp >= 0)
        out.push_back(p);
      if ((sp >= 0) != (sq >= 0)) {
        const double t = sp / (sp - sq);
        out.push_back({p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])});
      }
    }
    P = out;
  }
  return P;
}

// a sphere of radius r about the origin, the icosahedron divided four times,
// its points on the sphere
ls::SmartPointer<ls::Mesh<>> sphere(double r) {
  const double t = (1. + std::sqrt(5.)) / 2.;
  std::vector<std::array<double, 3>> V = {{-1, t, 0},  {1, t, 0},   {-1, -t, 0},
                                          {1, -t, 0},  {0, -1, t},  {0, 1, t},
                                          {0, -1, -t}, {0, 1, -t},  {t, 0, -1},
                                          {t, 0, 1},   {-t, 0, -1}, {-t, 0, 1}};
  std::vector<std::array<unsigned, 3>> F = {
      {0, 11, 5}, {0, 5, 1},  {0, 1, 7},   {0, 7, 10}, {0, 10, 11},
      {1, 5, 9},  {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
      {3, 9, 4},  {3, 4, 2},  {3, 2, 6},   {3, 6, 8},  {3, 8, 9},
      {4, 9, 5},  {2, 4, 11}, {6, 2, 10},  {8, 6, 7},  {9, 8, 1}};
  auto onSphere = [](std::array<double, 3> p) {
    const double n = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    return std::array<double, 3>{p[0] / n, p[1] / n, p[2] / n};
  };
  for (auto &p : V)
    p = onSphere(p);
  for (int level = 0; level < 4; ++level) {
    std::vector<std::array<unsigned, 3>> next;
    std::map<std::pair<unsigned, unsigned>, unsigned> middle;
    auto mid = [&](unsigned a, unsigned b) {
      const auto key = std::minmax(a, b);
      const auto found = middle.find(key);
      if (found != middle.end())
        return found->second;
      V.push_back(
          onSphere({V[a][0] + V[b][0], V[a][1] + V[b][1], V[a][2] + V[b][2]}));
      return middle[key] = unsigned(V.size() - 1);
    };
    for (const auto &f : F) {
      const unsigned ab = mid(f[0], f[1]), bc = mid(f[1], f[2]),
                     ca = mid(f[2], f[0]);
      next.push_back({f[0], ab, ca});
      next.push_back({f[1], bc, ab});
      next.push_back({f[2], ca, bc});
      next.push_back({ab, bc, ca});
    }
    F = next;
  }
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  for (const auto &p : V)
    mesh->insertNextNode({r * p[0], r * p[1], r * p[2]});
  for (const auto &f : F)
    mesh->insertNextTriangle(f);
  return mesh;
}

struct Result {
  double chamfer = 0., worst = 0.;
  bool valid = false;
};

// the surface a level set delivers against the exact mesh, in a box
template <int D>
Result measure(ls::SmartPointer<ls::Domain<double, D>> levelSet,
               ls::SmartPointer<ls::Mesh<>> exact, double gridDelta,
               const std::array<double, D> *lo,
               const std::array<double, D> *hi) {
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToSurfaceMesh<double, D>(levelSet, mesh).apply();
  ls::CompareSurfaceDistance<double, D> compare;
  compare.setMeshTarget(exact);
  compare.setMeshSample(mesh);
  compare.setSamplingSpacing(gridDelta / 8);
  if (lo && hi)
    compare.setBox(*lo, *hi);
  compare.apply();
  // both bands follow the layer rules of Check: a Euclidean band's values,
  // exact distances capped by the kept values next to the surface, change
  // by at most one grid spacing between axis neighbours
  bool valid =
      !mesh->getNodes().empty() && std::isfinite(compare.getMaxDistance());
  ls::Check<double, D> check(levelSet);
  check.apply();
  valid = valid && check.isValid();
  return {compare.getChamferDistance(), compare.getMaxDistance(), valid};
}

// one case on a Manhattan and a Euclidean level set: make(euclidean) builds
// the stack, the top last; true where the Euclidean run lies at least as
// close to the exact surface, or with a bound (grid spacings) where both
// runs lie within it of the exact surface
template <int D>
bool compareCase(
    const std::string &name,
    const std::function<
        std::vector<ls::SmartPointer<ls::Domain<double, D>>>(bool)> &make,
    ls::SmartPointer<ls::VelocityField<double>> field, double time,
    ls::SmartPointer<ls::Mesh<>> exact, double gridDelta,
    const std::array<double, D> *lo, const std::array<double, D> *hi,
    double bound = 0.) {
  Result result[2];
  unsigned steps[2];
  for (int euclidean = 0; euclidean < 2; ++euclidean) {
    auto stack = make(euclidean);
    ls::Advect<double, D> advect;
    for (auto &levelSet : stack)
      advect.insertNextLevelSet(levelSet);
    advect.setVelocityField(field);
    advect.setAdvectionTime(time);
    advect.setEuclideanRenormalization(euclidean);
    advect.apply();
    steps[euclidean] = advect.getNumberOfTimeSteps();
    result[euclidean] = measure<D>(stack.back(), exact, gridDelta, lo, hi);
  }
  const bool good = result[0].valid && result[1].valid &&
                    (bound > 0. ? result[0].worst <= bound * gridDelta &&
                                      result[1].worst <= bound * gridDelta
                                : result[1].chamfer <= result[0].chamfer &&
                                      result[1].worst <= result[0].worst);
  std::printf("%s: Manhattan %u steps, Chamfer %.4g, max %.4g; Euclidean %u "
              "steps, Chamfer %.4g, max %.4g (grid spacings)%s%s\n",
              name.c_str(), steps[0], result[0].chamfer / gridDelta,
              result[0].worst / gridDelta, steps[1],
              result[1].chamfer / gridDelta, result[1].worst / gridDelta,
              result[0].valid && result[1].valid ? "" : ", INVALID",
              good ? "" : ", FAILED");
  return good;
}

int main() {
  bool failed = false;
  const double gridDelta = 0.25, width = 6.;
  const double a = 2.37, ys = 0.29, h = 1.13;
  // 2D level sets of the reflective strip |x| <= width
  auto newLevelSet = [&](bool euclidean, double gs, bool reflective) {
    double bounds[4] = {-width, width, -6., 6.};
    ls::BoundaryConditionEnum bc[2] = {
        reflective ? ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY
                   : ls::BoundaryConditionEnum::INFINITE_BOUNDARY,
        ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
    auto levelSet =
        ls::SmartPointer<ls::Domain<double, 2>>::New(bounds, bc, gs);
    if (euclidean)
      levelSet->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
    return levelSet;
  };
  // the mask |x| <= a, reaching a little below the substrate surface, and the
  // surface of substrate and mask, as ViennaPS builds them
  auto maskAndTop = [&](bool euclidean) {
    auto mask = newLevelSet(euclidean, gridDelta, true);
    ls::FromSurfaceMesh<double, 2>(
        mask,
        polyline({{a, ys - 1e-4}, {-a, ys - 1e-4}, {-a, ys + h}, {a, ys + h}},
                 true))
        .apply();
    auto top = newLevelSet(euclidean, gridDelta, true);
    ls::FromSurfaceMesh<double, 2>(top, polyline({{-width, ys},
                                                  {-a, ys},
                                                  {-a, ys + h},
                                                  {a, ys + h},
                                                  {a, ys},
                                                  {width, ys}},
                                                 false))
        .apply();
    return std::vector<ls::SmartPointer<ls::Domain<double, 2>>>{mask, top};
  };
  const std::array<double, 2> lo{-width + gridDelta, -1e3},
      hi{width - gridDelta, 1e3};

  // directional etching beside the mask by r
  {
    const double r = 1.;
    auto exact = polyline({{-width, ys - r},
                           {-a, ys - r},
                           {-a, ys + h},
                           {a, ys + h},
                           {a, ys - r},
                           {width, ys - r}},
                          false);
    failed |= !compareCase<2>("2D directional etching beside a mask",
                              maskAndTop, ls::SmartPointer<Etch>::New(), r,
                              exact, gridDelta, &lo, &hi);
  }
  // growth by d beside the mask of zero rate
  {
    const double d = 0.75;
    auto exact = polyline({{-width, ys + d},
                           {-a, ys + d},
                           {-a, ys + h},
                           {a, ys + h},
                           {a, ys + d},
                           {width, ys + d}},
                          false);
    failed |= !compareCase<2>("2D growth beside a mask of zero rate",
                              maskAndTop, ls::SmartPointer<OnSubstrate>::New(),
                              d, exact, gridDelta, &lo, &hi);
  }
  // growth by d over the mask of zero rate, which the film covers from its
  // edges (exact: the film at ys + d beside the mask, arcs of radius d - h
  // around the mask's corners, the mask top between them); where the film
  // meets the mask, the values measured to the one and to the other used
  // to leave the Euclidean band invalid
  {
    const double d = 1.5, rc = d - h;
    Polygon P{{-width, ys + d}, {-a, ys + d}};
    for (int k = 1; k <= 32; ++k) {
      const double t = M_PI / 2 * (1. - k / 32.);
      P.push_back({-a + rc * std::cos(t), ys + h + rc * std::sin(t)});
    }
    for (int k = 0; k < 32; ++k) {
      const double t = M_PI / 2 * (2. - k / 32.);
      P.push_back({a + rc * std::cos(t), ys + h + rc * std::sin(t)});
    }
    P.push_back({a, ys + d});
    P.push_back({width, ys + d});
    failed |= !compareCase<2>("2D growth over a mask of zero rate", maskAndTop,
                              ls::SmartPointer<OnSubstrate>::New(), d,
                              polyline(P, false), gridDelta, &lo, &hi);
  }
  // a rotated rectangle etched from above by depth
  {
    const double depth = 2.5;
    const Polygon box = rectangle(0.3, -0.2, 2.1, 1.4, M_PI / 6);
    Polygon lowered = box;
    for (auto &p : lowered)
      p[1] -= depth;
    auto exact = polyline(clip(box, lowered), true);
    auto make = [&](bool euclidean) {
      auto levelSet = newLevelSet(euclidean, gridDelta, false);
      Polygon clockwise(box.rbegin(), box.rend());
      ls::FromSurfaceMesh<double, 2>(levelSet, polyline(clockwise, true))
          .apply();
      return std::vector<ls::SmartPointer<ls::Domain<double, 2>>>{levelSet};
    };
    failed |= !compareCase<2>("2D rectangle etched from above", make,
                              ls::SmartPointer<FromAbove>::New(), depth, exact,
                              gridDelta, nullptr, nullptr);
  }
  // a sphere of radius 1.5 grown and etched by 0.5, within half a grid
  // spacing of the exact sphere with both normalizations
  for (const double speed : {1., -1.}) {
    const double gs = 0.1, radius = 1.5, distance = 0.5;
    auto make = [&](bool euclidean) {
      double bounds[6] = {-3., 3., -3., 3., -3., 3.};
      ls::BoundaryConditionEnum bc[3] = {
          ls::BoundaryConditionEnum::INFINITE_BOUNDARY,
          ls::BoundaryConditionEnum::INFINITE_BOUNDARY,
          ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
      auto levelSet =
          ls::SmartPointer<ls::Domain<double, 3>>::New(bounds, bc, gs);
      if (euclidean)
        levelSet->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
      ls::FromSurfaceMesh<double, 3>(levelSet, sphere(radius)).apply();
      return std::vector<ls::SmartPointer<ls::Domain<double, 3>>>{levelSet};
    };
    failed |= !compareCase<3>(
        speed > 0 ? "3D sphere grown" : "3D sphere etched", make,
        ls::SmartPointer<Constant>::New(speed), distance,
        sphere(radius + speed * distance), gs, nullptr, nullptr, 0.5);
  }
  return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
