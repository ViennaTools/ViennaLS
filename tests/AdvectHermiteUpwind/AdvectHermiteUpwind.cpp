#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <lsAdvectHermite.hpp>
#include <lsCheck.hpp>
#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMesh.hpp>
#include <lsToDualContouringMesh.hpp>

/**
  AdvectHermite with the upwind update (setUpwind), the update ViennaPS's
  plasma etching models take on Euclidean domains: a substrate, the topmost
  material, etched directionally (the rate the upward component of the
  normal, as ions from above give) beside a mask that does not advect, the
  level set below it in the stack as ViennaPS builds a mask. The exact
  boundary is a trench whose walls go straight down from the mask's edges to
  the floor, without undercut. The mask keeps its faces and corners, the
  upwind values capped by its level set, the trench follows the exact
  boundary to first order in the grid spacing at its bottom corners (the
  bounds per scheme below), and the band stays a valid level set; with
  first- and second-order
  Engquist-Osher, second-order Lax-Friedrichs and WENO 5, and with
  third-order Runge-Kutta.

  Two more cases check the update against exact solutions where the rate
  depends on the normal or on the material below the surface: a rotated
  rectangle etched from above by a vector velocity (the exact body is the
  rectangle intersected with its translate by the etched depth, with no
  solid left beside its lower corners), and isotropic growth on the floor of
  an opening of a mask of zero rate (the exact surface keeps the mask's top
  and walls and lifts the floor by the grown thickness). \example
  AdvectHermiteUpwind.cpp
*/

namespace ls = viennals;

// the substrate (material 1) etched at the upward component of its normal,
// the mask (material 0) not
class Etch : public ls::VelocityField<double> {
public:
  double getScalarVelocity(const ls::Vec3D<double> &, int material,
                           const ls::Vec3D<double> &normal,
                           unsigned long) override {
    return material == 1 ? -std::max(normal[1], 0.) : 0.;
  }
};

double segmentDistance(double x, double y, double ax, double ay, double bx,
                       double by) {
  const double ex = bx - ax, ey = by - ay;
  const double t =
      std::clamp(((x - ax) * ex + (y - ay) * ey) / (ex * ex + ey * ey), 0., 1.);
  return std::hypot(x - ax - t * ex, y - ay - t * ey);
}

// distance to the mask |x| <= a on the substrate surface y = ys, its top and
// walls
double maskDistance(double x, double y, double a, double ys, double h) {
  x = std::abs(x);
  return std::min(segmentDistance(x, y, 0., ys + h, a, ys + h),
                  segmentDistance(x, y, a, ys + h, a, ys));
}

// distance to the exact boundary after directional etching by r beside the
// mask: its top, its walls continued straight down to the floor at ys - r,
// and the floor
double exactDistance(double x, double y, double a, double ys, double h,
                     double r, double width) {
  x = std::abs(x);
  double d = segmentDistance(x, y, 0., ys + h, a, ys + h);
  d = std::min(d, segmentDistance(x, y, a, ys + h, a, ys - r));
  return std::min(d, segmentDistance(x, y, a, ys - r, width, ys - r));
}

using Polygon = std::vector<std::array<double, 2>>;

// distance of (x, y) to the boundary of polygon P
double polygonDistance(double x, double y, const Polygon &P) {
  double d = std::numeric_limits<double>::max();
  for (std::size_t k = 0; k < P.size(); ++k) {
    const auto &a = P[k], &b = P[(k + 1) % P.size()];
    d = std::min(d, segmentDistance(x, y, a[0], a[1], b[0], b[1]));
  }
  return d;
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

// a closed polygon as lines with the outside on the left (clockwise), as
// FromSurfaceMesh takes a solid
ls::SmartPointer<ls::Mesh<>> outline(const Polygon &P) {
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  for (auto it = P.rbegin(); it != P.rend(); ++it)
    mesh->insertNextNode({(*it)[0], (*it)[1], 0.});
  for (unsigned k = 0; k < P.size(); ++k)
    mesh->insertNextLine({k, unsigned((k + 1) % P.size())});
  return mesh;
}

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

// the worst distance of the vertices of the dual contouring mesh of
// levelSet from the exact boundary (grid spacings), whether the band is a
// valid level set
template <class F>
bool measure(ls::SmartPointer<ls::Domain<double, 2>> levelSet, F exact,
             double gridDelta, double &worst, double &worstX, double &worstY) {
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToDualContouringMesh<double, 2>(levelSet, mesh).apply();
  worst = 0.;
  for (const auto &n : mesh->getNodes()) {
    const double e = exact(n[0], n[1]) / gridDelta;
    if (e > worst) {
      worst = e;
      worstX = n[0];
      worstY = n[1];
    }
  }
  ls::Check<double, 2> check(levelSet);
  check.apply();
  return !mesh->getNodes().empty() && check.isValid();
}

int main() {
  const double gridDelta = 0.25, r = 1., width = 6.;
  const double a = 2.37, ys = 0.29, h = 1.13;
  struct Case {
    std::string name;
    ls::SpatialSchemeEnum scheme;
    ls::TemporalSchemeEnum time;
    // the worst vertex of the trench and of the mask, grid spacings: the
    // first-order rounding of the trench's bottom corners, up to three
    // quarters of a grid spacing (one with Runge-Kutta 3); the mask does not
    // move in the exact solution and keeps its faces with every scheme
    double bound, maskBound;
  };
  const Case cases[] = {
      {"Engquist-Osher", ls::SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER,
       ls::TemporalSchemeEnum::FORWARD_EULER, 0.75, 1e-3},
      {"Engquist-Osher, Runge-Kutta 3",
       ls::SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER,
       ls::TemporalSchemeEnum::RUNGE_KUTTA_3RD_ORDER, 1.0, 1e-3},
      {"Engquist-Osher 2nd order",
       ls::SpatialSchemeEnum::ENGQUIST_OSHER_2ND_ORDER,
       ls::TemporalSchemeEnum::FORWARD_EULER, 0.75, 1e-3},
      {"Lax-Friedrichs 2nd order",
       ls::SpatialSchemeEnum::LAX_FRIEDRICHS_2ND_ORDER,
       ls::TemporalSchemeEnum::FORWARD_EULER, 0.75, 1e-3},
      {"WENO 5", ls::SpatialSchemeEnum::WENO_5TH_ORDER,
       ls::TemporalSchemeEnum::RUNGE_KUTTA_3RD_ORDER, 0.75, 1e-3}};
  bool failed = false;
  for (const auto &c : cases) {
    double bounds[4] = {-width, width, -6., 6.};
    ls::BoundaryConditionEnum boundaryCons[2] = {
        ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
        ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
    auto newLevelSet = [&]() {
      auto levelSet = ls::SmartPointer<ls::Domain<double, 2>>::New(
          bounds, boundaryCons, gridDelta);
      levelSet->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
      levelSet->setLevelSetWidth(12);
      return levelSet;
    };
    // the mask, reaching a little below the substrate surface as in
    // ViennaPS, below in the stack
    auto mask = newLevelSet();
    {
      auto outline = ls::SmartPointer<ls::Mesh<>>::New();
      const double pts[4][2] = {
          {a, ys - 1e-4}, {-a, ys - 1e-4}, {-a, ys + h}, {a, ys + h}};
      for (const auto &q : pts)
        outline->insertNextNode({q[0], q[1], 0.});
      for (unsigned k = 0; k < 4; ++k)
        outline->insertNextLine({k, (k + 1) % 4});
      ls::FromSurfaceMesh<double, 2>(mask, outline).apply();
    }
    // the surface of substrate and mask
    auto top = newLevelSet();
    {
      auto surface = ls::SmartPointer<ls::Mesh<>>::New();
      const double pts[6][2] = {{-width, ys}, {-a, ys}, {-a, ys + h},
                                {a, ys + h},  {a, ys},  {width, ys}};
      for (const auto &q : pts)
        surface->insertNextNode({q[0], q[1], 0.});
      for (unsigned k = 0; k < 5; ++k)
        surface->insertNextLine({k, k + 1});
      ls::FromSurfaceMesh<double, 2>(top, surface).apply();
    }

    ls::AdvectHermite<double, 2> advect(
        std::vector<ls::SmartPointer<ls::Domain<double, 2>>>{mask, top},
        ls::SmartPointer<Etch>::New());
    advect.setUpwind(true);
    advect.setSpatialScheme(c.scheme);
    advect.setTemporalScheme(c.time);
    advect.setAdvectionTime(r);
    advect.apply();

    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToDualContouringMesh<double, 2>(top, mesh).apply();
    double deviation = 0., maskDeviation = 0., worstX = 0., worstY = 0.;
    for (const auto &n : mesh->getNodes()) {
      if (std::abs(n[0]) > width - gridDelta)
        continue; // the reflective boundary
      const double e = exactDistance(n[0], n[1], a, ys, h, r, width);
      if (e > deviation) {
        deviation = e;
        worstX = n[0];
        worstY = n[1];
      }
      // the vertices on the mask's top and walls
      if (n[1] > ys + 0.5 * gridDelta && std::abs(n[0]) < a + 0.5 * gridDelta)
        maskDeviation =
            std::max(maskDeviation, maskDistance(n[0], n[1], a, ys, h));
    }
    ls::Check<double, 2> check(top);
    check.apply();
    std::cout << "2D directional etching beside a mask, upwind update ("
              << c.name << "): " << advect.getNumberOfTimeSteps()
              << " steps, worst vertex " << deviation / gridDelta << " at ("
              << worstX << ", " << worstY << "), on the mask "
              << maskDeviation / gridDelta << " grid spacings, "
              << (check.isValid() ? "valid" : "INVALID") << std::endl;
    if (mesh->getNodes().empty() || !check.isValid() ||
        deviation / gridDelta > c.bound ||
        maskDeviation / gridDelta > c.maskBound) {
      std::cout << "FAILED: " << c.name << std::endl;
      failed = true;
    }
  }

  // a rotated rectangle etched from above by a vector velocity: the exact
  // body is the rectangle intersected with its translate down by the etched
  // depth; its lower faces stand still while the upper ones sweep past
  // their nodes, which leaves no solid beside the lower corners
  {
    const double depth = 2.5;
    const Polygon box = rectangle(0.3, -0.2, 2.1, 1.4, M_PI / 6);
    Polygon lowered = box;
    for (auto &p : lowered)
      p[1] -= depth;
    const Polygon exact = clip(box, lowered);
    double bounds[4] = {-width, width, -width, width};
    ls::BoundaryConditionEnum boundaryCons[2] = {
        ls::BoundaryConditionEnum::INFINITE_BOUNDARY,
        ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
    auto levelSet = ls::SmartPointer<ls::Domain<double, 2>>::New(
        bounds, boundaryCons, gridDelta);
    levelSet->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
    levelSet->setLevelSetWidth(12);
    ls::FromSurfaceMesh<double, 2>(levelSet, outline(box)).apply();
    ls::AdvectHermite<double, 2> advect(
        std::vector<ls::SmartPointer<ls::Domain<double, 2>>>{levelSet},
        ls::SmartPointer<FromAbove>::New());
    advect.setUpwind(true);
    advect.setAdvectionTime(depth);
    advect.apply();
    double worst = 0., worstX = 0., worstY = 0.;
    const bool valid = measure(
        levelSet,
        [&](double x, double y) { return polygonDistance(x, y, exact); },
        gridDelta, worst, worstX, worstY);
    std::cout << "2D rectangle etched from above, upwind update: "
              << advect.getNumberOfTimeSteps() << " steps, worst vertex "
              << worst << " grid spacings at (" << worstX << ", " << worstY
              << "), " << (valid ? "valid" : "INVALID") << std::endl;
    // first-order rounding at the corners, up to three quarters of a grid
    // spacing
    if (!valid || worst > 0.75) {
      std::cout << "FAILED: rectangle etched from above" << std::endl;
      failed = true;
    }
  }

  // isotropic growth by d on the substrate beside the mask |x| <= a, whose
  // rate is zero, the mask lowest in the stack: the exact surface keeps the
  // mask's top and walls and lifts the substrate surface to ys + d
  {
    const double d = 0.75;
    double bounds[4] = {-width, width, -6., 6.};
    ls::BoundaryConditionEnum boundaryCons[2] = {
        ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
        ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
    auto newLevelSet = [&]() {
      auto levelSet = ls::SmartPointer<ls::Domain<double, 2>>::New(
          bounds, boundaryCons, gridDelta);
      levelSet->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
      levelSet->setLevelSetWidth(12);
      return levelSet;
    };
    auto mask = newLevelSet();
    ls::FromSurfaceMesh<double, 2>(
        mask,
        outline({{-a, ys - 1e-4}, {a, ys - 1e-4}, {a, ys + h}, {-a, ys + h}}))
        .apply();
    auto top = newLevelSet();
    {
      auto surface = ls::SmartPointer<ls::Mesh<>>::New();
      const double pts[6][2] = {{-width, ys}, {-a, ys}, {-a, ys + h},
                                {a, ys + h},  {a, ys},  {width, ys}};
      for (const auto &q : pts)
        surface->insertNextNode({q[0], q[1], 0.});
      for (unsigned k = 0; k < 5; ++k)
        surface->insertNextLine({k, k + 1});
      ls::FromSurfaceMesh<double, 2>(top, surface).apply();
    }
    ls::AdvectHermite<double, 2> advect(
        std::vector<ls::SmartPointer<ls::Domain<double, 2>>>{mask, top},
        ls::SmartPointer<OnSubstrate>::New());
    advect.setUpwind(true);
    advect.setAdvectionTime(d);
    advect.apply();
    auto exact = [&](double x, double y) {
      x = std::abs(x);
      double e = segmentDistance(x, y, 0., ys + h, a, ys + h);
      e = std::min(e, segmentDistance(x, y, a, ys + h, a, ys + d));
      return std::min(e, segmentDistance(x, y, a, ys + d, width, ys + d));
    };
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToDualContouringMesh<double, 2>(top, mesh).apply();
    double worst = 0., worstX = 0., worstY = 0., maskWorst = 0.;
    for (const auto &n : mesh->getNodes()) {
      if (std::abs(n[0]) > width - gridDelta)
        continue; // the reflective boundary
      const double e = exact(n[0], n[1]) / gridDelta;
      if (e > worst) {
        worst = e;
        worstX = n[0];
        worstY = n[1];
      }
      // the vertices on the mask's top and walls above the grown surface
      if (n[1] > ys + d + gridDelta && std::abs(n[0]) < a + 0.5 * gridDelta)
        maskWorst = std::max(maskWorst, e);
    }
    ls::Check<double, 2> check(top);
    check.apply();
    std::cout << "2D growth beside a mask of zero rate, upwind update: "
              << advect.getNumberOfTimeSteps() << " steps, worst vertex "
              << worst << " grid spacings at (" << worstX << ", " << worstY
              << "), on the mask " << maskWorst << ", "
              << (check.isValid() ? "valid" : "INVALID") << std::endl;
    // the corners where the grown surface meets the mask's walls to first
    // order, three quarters of a grid spacing; the mask does not move
    if (mesh->getNodes().empty() || !check.isValid() || worst > 0.75 ||
        maskWorst > 1e-3) {
      std::cout << "FAILED: growth beside a mask" << std::endl;
      failed = true;
    }
  }
  return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
