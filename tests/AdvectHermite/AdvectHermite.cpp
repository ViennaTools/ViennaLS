#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <lsAdvectHermite.hpp>
#include <lsBooleanOperation.hpp>
#include <lsCheck.hpp>
#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMesh.hpp>
#include <lsToDualContouringMesh.hpp>

/**
  Advection of rotated boxes with AdvectHermite in 2D and 3D by 2.25 grid
  spacings, compared with the exact results: uniform etching keeps the
  shrunk box with its sharp corners exactly, uniform deposition rounds the
  convex corners into arcs (the dual contouring polygon around them stays
  within a quarter grid spacing), directional etching moves the faces facing
  up while the others stay, and directional deposition sweeps the box up,
  with new faces at the edges between moving and static faces. Every process
  runs once in one call and once in five calls with the rate changed between
  them (same total motion), which continues from the band state the previous
  call left. A ridge across
  a periodic boundary (2D), and a substrate etched beside a mask (2D) and
  through a square opening of a mask (3D), check periodic boundaries and
  material stacks, also with the mask below in the stack and in single steps
  as ViennaPS runs them, the mask reaching a little below the substrate
  surface as MakeTrench builds it; their vertices include those at the
  reflective boundaries. A mask on a sloped substrate (2D) checks the cones that
  the interface splits at the feet of the mask, and the slope alone a Manhattan
  substrate below. A 2D trench etched with a noisy rate must not reach below
  the mask from a folded mesh at the corner of its bottom. The rotated boxes run
  once more in single precision. A trench through a mask, etched with a noisy
  view factor in single steps, must end as a valid level set in double and
  single precision, also with the mask below in the stack, and keep its mesh out
  of the mask. \example AdvectHermite.cpp
*/

namespace ls = viennals;

using Point = std::array<double, 3>;

int fail(int D, const std::string &message) {
  std::cout << D << "D: " << message << std::endl;
  return EXIT_FAILURE;
}

template <int D> struct RotatedBox {
  Point half, centre;
  std::array<Point, 3> R;

  RotatedBox(double scale) {
    if constexpr (D == 2) {
      half = {5.15 * scale, 3.85 * scale, 0.};
      centre = {0.31 * scale, -0.27 * scale, 0.};
    } else {
      half = {5.15 * scale, 4.35 * scale, 3.55 * scale};
      centre = {0.31 * scale, -0.27 * scale, 0.43 * scale};
    }
    // Rz(30°), followed by Rx(20°) in 3D, as in EuclideanLS/advect_prototype.py
    const double a = 30. * M_PI / 180., b = (D == 3) ? 20. * M_PI / 180. : 0.;
    const double ca = std::cos(a), sa = std::sin(a), cb = std::cos(b),
                 sb = std::sin(b);
    R = {{{ca, -sa, 0.}, {cb * sa, cb * ca, -sb}, {sb * sa, sb * ca, cb}}};
  }

  double sdf(const Point &p, double inset = 0.) const {
    double outside = 0., inside = -1e300;
    for (int j = 0; j < D; ++j) {
      double l = 0.;
      for (int i = 0; i < 3; ++i)
        l += R[i][j] * (p[i] - centre[i]);
      const double q = std::abs(l) - (half[j] - inset);
      outside += (q > 0.) ? q * q : 0.;
      inside = std::max(inside, q);
    }
    return std::sqrt(outside) + std::min(inside, 0.);
  }

  std::vector<Point> corners(double inset = 0.) const {
    std::vector<Point> result;
    for (int n = 0; n < (1 << D); ++n) {
      Point p = centre;
      for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j)
          p[i] +=
              R[i][j] * (((n >> j) & 1) ? half[j] - inset : inset - half[j]);
      result.push_back(p);
    }
    return result;
  }

  // closed mesh with the corner order and orientation of MakeGeometry::makeBox
  // signed distance to the box swept up by r (directional deposition), a
  // convex body: outside, the distance to the union of the shifted boxes,
  // convex in the shift; inside, the largest signed distance to its
  // supporting planes, among them its facets with the normals of the box and
  // those across an edge of the box and y
  double sweptSdf(const Point &p, double r) const {
    auto shifted = [&](double s) { return sdf({p[0], p[1] - s, p[2]}); };
    const double g = (std::sqrt(5.) - 1.) / 2.;
    double lo = 0., hi = r, x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
    double f1 = shifted(x1), f2 = shifted(x2);
    while (hi - lo > 1e-13 * std::max(1., r)) {
      if (f1 < f2) {
        hi = x2;
        x2 = x1;
        f2 = f1;
        x1 = hi - g * (hi - lo);
        f1 = shifted(x1);
      } else {
        lo = x1;
        x1 = x2;
        f1 = f2;
        x2 = lo + g * (hi - lo);
        f2 = shifted(x2);
      }
    }
    const double outside = std::min({shifted(0.), shifted(r), f1, f2});
    if (outside > 0.)
      return outside;
    std::vector<Point> vertices;
    for (const auto &c : corners()) {
      vertices.push_back(c);
      vertices.push_back({c[0], c[1] + r, c[2]});
    }
    std::vector<Point> normals;
    for (int j = 0; j < D; ++j) {
      const Point a{R[0][j], R[1][j], R[2][j]};
      Point across{1., 0., 0.};
      if constexpr (D == 3) {
        const double l = std::hypot(a[0], a[2]);
        across = {-a[2] / l, 0., a[0] / l};
      }
      for (const Point &n : {a, across}) {
        normals.push_back(n);
        normals.push_back({-n[0], -n[1], -n[2]});
      }
    }
    double inside = -1e300;
    for (const auto &n : normals) {
      double support = -1e300;
      for (const auto &v : vertices)
        support = std::max(support, n[0] * v[0] + n[1] * v[1] + n[2] * v[2]);
      inside =
          std::max(inside, n[0] * p[0] + n[1] * p[1] + n[2] * p[2] - support);
    }
    return inside;
  }

  ls::SmartPointer<ls::Mesh<>> mesh() const {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    for (const auto &c : corners())
      mesh->insertNextNode({c[0], c[1], c[2]});
    if constexpr (D == 2) {
      const std::array<unsigned, 2> lines[4] = {{0, 2}, {2, 3}, {3, 1}, {1, 0}};
      for (const auto &l : lines)
        mesh->insertNextLine(l);
    } else {
      const std::array<unsigned, 3> triangles[12] = {
          {0, 3, 1}, {0, 2, 3}, {0, 1, 5}, {0, 5, 4}, {0, 4, 2}, {4, 6, 2},
          {7, 6, 4}, {7, 4, 5}, {7, 2, 6}, {7, 3, 2}, {1, 3, 5}, {3, 7, 5}};
      for (const auto &t : triangles)
        mesh->insertNextTriangle(t);
    }
    return mesh;
  }
};

class Rated : public ls::VelocityField<double> {
public:
  double rate = 1.;
};

class Uniform : public Rated {
  double v;

public:
  explicit Uniform(double velocity) : v(velocity) {}
  double getScalarVelocity(const ls::Vec3D<double> &, int,
                           const ls::Vec3D<double> &, unsigned long) override {
    return rate * v;
  }
};

// etching of the faces looking up, speed proportional to that component
class Directional : public Rated {
public:
  double getScalarVelocity(const ls::Vec3D<double> &, int,
                           const ls::Vec3D<double> &normal,
                           unsigned long) override {
    return -rate * std::max(normal[1], 0.);
  }
};

// deposition on the faces looking up, speed proportional to that component
class DirectionalDeposition : public Rated {
public:
  double getScalarVelocity(const ls::Vec3D<double> &, int,
                           const ls::Vec3D<double> &normal,
                           unsigned long) override {
    return rate * std::max(normal[1], 0.);
  }
};

template <int D>
ls::SmartPointer<ls::Domain<double, D>> levelSet(const RotatedBox<D> &box,
                                                 double gridDelta) {
  double bounds[2 * D];
  ls::BoundaryConditionEnum boundaryCons[D];
  for (int i = 0; i < D; ++i) {
    bounds[2 * i] = -18. * gridDelta;
    bounds[2 * i + 1] = 18. * gridDelta;
    boundaryCons[i] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
  }
  auto domain = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  domain->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
  domain->setLevelSetWidth(12);
  ls::FromSurfaceMesh<double, D> mesher(domain, box.mesh());
  mesher.setUseWindingNumber(true);
  mesher.apply();
  return domain;
}

template <int D> int check(double gridDelta) {
  const RotatedBox<D> box(gridDelta);
  const double r = 2.25 * gridDelta;
  struct Case {
    std::string name;
    ls::SmartPointer<Rated> velocity;
    double tolerance; // worst vertex, grid spacings
  };
  const std::vector<Case> processes = {
      {"etching", ls::SmartPointer<Uniform>::New(-1.), 1e-6},
      {"deposition", ls::SmartPointer<Uniform>::New(1.), 0.3},
      {"directional etching", ls::SmartPointer<Directional>::New(), 0.05},
      {"directional deposition", ls::SmartPointer<DirectionalDeposition>::New(),
       D == 2 ? 1e-6 : 0.1}};
  // one call with rate 1, or five calls with changing rates and the same
  // total motion
  const std::vector<std::vector<double>> schedules = {
      {1.}, {1.0, 0.6, 1.4, 0.8, 1.2}};
  for (const auto &schedule : schedules)
    for (auto c : processes) {
      const std::string kind = c.name;
      if (schedule.size() > 1)
        c.name += " in " + std::to_string(schedule.size()) + " calls";
      auto domain = levelSet<D>(box, gridDelta);
      const auto start = std::chrono::steady_clock::now();
      unsigned steps = 0, rebuilds = 0;
      for (double rate : schedule) {
        c.velocity->rate = rate;
        ls::AdvectHermite<double, D> advect(domain, c.velocity);
        advect.setAdvectionTime(2.25 * gridDelta / schedule.size());
        // the calls continue from the band state the previous one left
        advect.setRebuildAtEnd(schedule.size() == 1);
        advect.apply();
        steps += advect.getNumberOfTimeSteps();
        rebuilds += advect.getNumberOfRebuilds();
      }
      const double seconds = std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() - start)
                                 .count();
      auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
      ls::ToDualContouringMesh<double, D>(domain, mesh).apply();
      if (mesh->getNodes().empty())
        return fail(D, c.name + ": no surface");

      double deviation = 0.;
      for (const auto &n : mesh->getNodes()) {
        const Point p{n[0], n[1], n[2]};
        double exact;
        if (kind == "etching")
          exact = box.sdf(p, r);
        else if (kind == "deposition")
          exact = box.sdf(p) - r;
        else
          exact = (kind == "directional deposition")
                      ? box.sweptSdf(p, r)
                      : std::max(box.sdf(p), box.sdf({p[0], p[1] + r, p[2]}));
        deviation = std::max(deviation, std::abs(exact));
      }
      double cornerError = 0.;
      if (kind == "etching")
        for (const auto &corner : box.corners(r)) {
          double best = 1e300;
          for (const auto &n : mesh->getNodes())
            best = std::min(best, std::hypot(n[0] - corner[0], n[1] - corner[1],
                                             n[2] - corner[2]));
          cornerError = std::max(cornerError, best);
        }
      std::cout << D << "D " << c.name << ": " << steps << " steps, "
                << rebuilds << " rebuild(s), " << seconds << " s, worst vertex "
                << deviation / gridDelta << ", corners "
                << cornerError / gridDelta << " grid spacings" << std::endl;
      if (deviation / gridDelta > c.tolerance ||
          cornerError / gridDelta > std::max(c.tolerance, 1e-6))
        return fail(D, c.name + " off the exact result");
    }
  return EXIT_SUCCESS;
}

// 2D surface periodic in x with a trapezoidal ridge across the periodic
// boundary, material below
struct PeriodicRidge {
  double L, y0, y1, top, slope;
  std::vector<std::array<double, 2>> points() const {
    return {{-L, y1},
            {-L + top, y1},
            {-L + top + slope, y0},
            {L - top - slope, y0},
            {L - top, y1},
            {L, y1}};
  }
  double height(double x) const {
    x = x - 2 * L * std::floor((x + L) / (2 * L)); // into [-L, L)
    const auto pts = points();
    for (std::size_t k = 0; k + 1 < pts.size(); ++k)
      if (x >= pts[k][0] && x <= pts[k + 1][0])
        return pts[k][1] + (pts[k + 1][1] - pts[k][1]) * (x - pts[k][0]) /
                               (pts[k + 1][0] - pts[k][0]);
    return y1;
  }
  // exact signed distance, periodic images included
  double sdf(double x, double y) const {
    const auto pts = points();
    double d = 1e300;
    for (int image = -1; image <= 1; ++image)
      for (std::size_t k = 0; k + 1 < pts.size(); ++k) {
        const double ax = pts[k][0] + 2 * L * image, ay = pts[k][1];
        const double bx = pts[k + 1][0] + 2 * L * image, by = pts[k + 1][1];
        const double ex = bx - ax, ey = by - ay;
        const double t = std::clamp(
            ((x - ax) * ex + (y - ay) * ey) / (ex * ex + ey * ey), 0., 1.);
        d = std::min(d, std::hypot(x - ax - t * ex, y - ay - t * ey));
      }
    return (y < height(x)) ? -d : d;
  }
};

int checkPeriodic() {
  const double gridDelta = 0.5, r = 2.25 * gridDelta;
  const PeriodicRidge ridge{6., 0.37, 0.37 + 1.33, 0.83, 0.71};
  struct Case {
    std::string name;
    ls::SmartPointer<Rated> velocity;
    double tolerance;
  };
  const std::vector<Case> cases = {
      {"etching", ls::SmartPointer<Uniform>::New(-1.), 0.3},
      {"deposition", ls::SmartPointer<Uniform>::New(1.), 0.3},
      {"directional etching", ls::SmartPointer<Directional>::New(), 1e-6}};
  for (const auto &c : cases) {
    double bounds[4] = {-ridge.L, ridge.L, -5., 8.};
    ls::BoundaryConditionEnum boundaryCons[2] = {
        ls::BoundaryConditionEnum::PERIODIC_BOUNDARY,
        ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
    auto domain = ls::SmartPointer<ls::Domain<double, 2>>::New(
        bounds, boundaryCons, gridDelta);
    domain->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
    domain->setLevelSetWidth(12);
    auto surface = ls::SmartPointer<ls::Mesh<>>::New();
    for (const auto &q : ridge.points())
      surface->insertNextNode({q[0], q[1], 0.});
    for (unsigned k = 0; k + 1 < ridge.points().size(); ++k)
      surface->insertNextLine({k, k + 1});
    ls::FromSurfaceMesh<double, 2>(domain, surface).apply();

    ls::AdvectHermite<double, 2> advect(domain, c.velocity);
    advect.setAdvectionTime(r);
    advect.apply();
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToDualContouringMesh<double, 2>(domain, mesh).apply();
    double deviation = 0., xMin = 1e300, xMax = -1e300;
    for (const auto &n : mesh->getNodes()) {
      double exact;
      if (c.name == "etching")
        exact = ridge.sdf(n[0], n[1]) + r;
      else if (c.name == "deposition")
        exact = ridge.sdf(n[0], n[1]) - r;
      else
        exact = ridge.sdf(n[0], n[1] + r);
      deviation = std::max(deviation, std::abs(exact));
      xMin = std::min(xMin, n[0]);
      xMax = std::max(xMax, n[0]);
    }
    std::cout << "2D periodic " << c.name << ": "
              << advect.getNumberOfTimeSteps() << " steps, "
              << advect.getNumberOfRebuilds() << " rebuild(s), worst vertex "
              << deviation / gridDelta << " grid spacings, x from " << xMin
              << " to " << xMax << std::endl;
    // the surface reaches across the periodic boundary on both sides
    if (deviation / gridDelta > c.tolerance || xMin > -ridge.L + gridDelta ||
        xMax < ridge.L - gridDelta)
      return fail(2, "periodic " + c.name + " off the exact result");
  }
  return EXIT_SUCCESS;
}

// substrate (material 0) etched beside a mask (material 1, not etched)
class MaskedEtch : public ls::VelocityField<double> {
  bool directional;

public:
  explicit MaskedEtch(bool passedDirectional)
      : directional(passedDirectional) {}
  double getScalarVelocity(const ls::Vec3D<double> &, int material,
                           const ls::Vec3D<double> &normal,
                           unsigned long) override {
    if (material != 0)
      return 0.;
    return directional ? -std::max(normal[1], 0.) : -1.;
  }
};

double segmentDistance(double x, double y, double ax, double ay, double bx,
                       double by) {
  const double ex = bx - ax, ey = by - ay;
  const double t =
      std::clamp(((x - ax) * ex + (y - ay) * ey) / (ex * ex + ey * ey), 0., 1.);
  return std::hypot(x - ax - t * ex, y - ay - t * ey);
}

// distance to the exact boundary after etching by r beside a mask |x| <= a
// on the substrate surface y = ys: a trench with vertical walls
// (directional), or a trench undercutting the mask with quarter circles of
// radius r around its bottom edges (isotropic)
double maskedDistance(double x, double y, double a, double ys, double h,
                      double r, double width, bool directional) {
  x = std::abs(x);                                         // symmetric
  double d = segmentDistance(x, y, 0., ys + h, a, ys + h); // mask top
  if (directional) {
    d = std::min(d, segmentDistance(x, y, a, ys + h, a, ys - r));
  } else {
    d = std::min(d, segmentDistance(x, y, a, ys + h, a, ys)); // mask wall
    d = std::min(d, segmentDistance(x, y, a, ys, a - r, ys)); // underside
    // quarter circle around (a, ys) from (a - r, ys) to (a, ys - r)
    const double dx = x - a, dy = y - ys;
    if (dx <= 0 && dy <= 0)
      d = std::min(d, std::abs(std::hypot(dx, dy) - r));
    else
      d = std::min(d, std::min(std::hypot(x - (a - r), y - ys),
                               std::hypot(x - a, y - (ys - r))));
  }
  return std::min(d, segmentDistance(x, y, a, ys - r, width, ys - r));
}

int checkMaterials() {
  const double gridDelta = 0.5, r = 2.25 * gridDelta, width = 6.;
  const double a = 2.37, ys = 0.29, h = 1.13;
  for (bool directional : {true, false}) {
    double bounds[4] = {-width, width, -6., 6.};
    ls::BoundaryConditionEnum boundaryCons[2] = {
        ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
        ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
    // substrate below, Manhattan as usual
    auto substrate = ls::SmartPointer<ls::Domain<double, 2>>::New(
        bounds, boundaryCons, gridDelta);
    auto plane = ls::SmartPointer<ls::Mesh<>>::New();
    plane->insertNextNode({-width, ys, 0.});
    plane->insertNextNode({width, ys, 0.});
    plane->insertNextLine({0, 1});
    ls::FromSurfaceMesh<double, 2>(substrate, plane).apply();
    // surface of substrate and mask, Euclidean with Hermite data
    auto top = ls::SmartPointer<ls::Domain<double, 2>>::New(
        bounds, boundaryCons, gridDelta);
    top->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
    top->setLevelSetWidth(12);
    auto surface = ls::SmartPointer<ls::Mesh<>>::New();
    const double pts[6][2] = {{-width, ys}, {-a, ys}, {-a, ys + h},
                              {a, ys + h},  {a, ys},  {width, ys}};
    for (const auto &q : pts)
      surface->insertNextNode({q[0], q[1], 0.});
    for (unsigned k = 0; k < 5; ++k)
      surface->insertNextLine({k, k + 1});
    ls::FromSurfaceMesh<double, 2>(top, surface).apply();

    ls::AdvectHermite<double, 2> advect(
        std::vector<ls::SmartPointer<ls::Domain<double, 2>>>{substrate, top},
        ls::SmartPointer<MaskedEtch>::New(directional));
    advect.setAdvectionTime(r);
    advect.apply();
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToDualContouringMesh<double, 2>(top, mesh).apply();
    double deviation = 0.;
    for (const auto &n : mesh->getNodes()) {
      if (std::abs(n[0]) > width - gridDelta)
        continue; // the reflective boundary
      deviation = std::max(deviation, maskedDistance(n[0], n[1], a, ys, h, r,
                                                     width, directional));
    }
    const std::string name =
        directional ? "directional etching" : "isotropic etching";
    std::cout << "2D " << name
              << " beside a mask: " << advect.getNumberOfTimeSteps()
              << " steps, " << advect.getNumberOfRebuilds()
              << " rebuild(s), worst vertex " << deviation / gridDelta
              << " grid spacings" << std::endl;
    if (mesh->getNodes().empty() ||
        deviation / gridDelta > (directional ? 0.05 : 0.3))
      return fail(2, name + " beside a mask off the exact result");
  }
  return EXIT_SUCCESS;
}

// the stack in the order of ViennaPS's MakeTrench: the mask (material 0, not
// etched) below, substrate and mask (material 1, etched) on top
class ReversedEtch : public ls::VelocityField<double> {
  bool directional;

public:
  explicit ReversedEtch(bool passedDirectional)
      : directional(passedDirectional) {}
  double getScalarVelocity(const ls::Vec3D<double> &, int material,
                           const ls::Vec3D<double> &normal,
                           unsigned long) override {
    if (material != 1)
      return 0.;
    return directional ? -std::max(normal[1], 0.) : -1.;
  }
};

// distance to the exact boundary after etching by r through an opening
// |x| < a of a mask on the substrate surface y = ys: a pit of the points
// within r of the opening, undercutting the mask
double openingProfileDistance(double x, double y, double a, double ys, double h,
                              double r) {
  x = std::abs(x);
  const double sd = x - a;
  double d = std::hypot(y - (ys + h), std::max(-sd, 0.));
  d = std::min(d, std::hypot(sd, std::max({ys - y, y - (ys + h), 0.})));
  d = std::min(d, std::hypot(y - ys, std::max({-sd, sd - r, 0.})));
  double pit;
  if (y <= ys) {
    pit = std::abs(std::hypot(y - ys, std::max(sd, 0.)) - r);
  } else {
    pit = std::hypot(y - ys, sd - r);
    if (sd <= 0.)
      pit = std::min(pit, y - ys + r);
  }
  return std::min(d, pit);
}

// distance to the exact boundary after directional etching by r through an
// opening |x| < a of a mask with its top at yt on the substrate surface
// y = ys: a trench with vertical walls
double openingTrenchDistance(double x, double y, double a, double ys, double yt,
                             double r) {
  x = std::abs(x);
  double d = std::hypot(y - yt, std::max(a - x, 0.)); // mask top
  d = std::min(d, std::hypot(x - a, std::max({ys - r - y, y - yt, 0.})));
  return std::min(d, std::hypot(y - (ys - r), std::max(x - a, 0.)));
}

// a mask reaching eps below the substrate surface leaves a sliver of mask
// just behind the foot of its wall; the undercut below it must not take the
// mask's material from there (1.28 grid spacings off at grid spacing 0.125
// before, and 1.12 at 0.5 with eps 0.01), and with directional etching the
// substrate below the sliver keeps the wall's direction, a sliver splitting
// nothing in the cone of the fan there
int checkReversedStack(double gridDelta, double eps, bool directional = false) {
  const double r = 1.125, width = 6.;
  const double a = 2.37, ys = 0.29, h = 1.13;
  double bounds[4] = {-width, width, -6., 6.};
  ls::BoundaryConditionEnum boundaryCons[2] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
  auto levelSetOf = [&](const std::vector<std::array<double, 2>> &points) {
    auto domain = ls::SmartPointer<ls::Domain<double, 2>>::New(
        bounds, boundaryCons, gridDelta);
    domain->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
    domain->setLevelSetWidth(12);
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    for (const auto &q : points)
      mesh->insertNextNode({q[0], q[1], 0.});
    for (unsigned k = 0; k + 1 < points.size(); ++k)
      mesh->insertNextLine({k, k + 1});
    ls::FromSurfaceMesh<double, 2>(domain, mesh).apply();
    return domain;
  };
  // mask slabs |x| >= a reaching eps below the substrate surface, as
  // MakeTrench builds them
  auto mask = levelSetOf({{-width, ys - eps + h},
                          {-a, ys - eps + h},
                          {-a, ys - eps},
                          {-width, ys - eps}});
  ls::BooleanOperation<double, 2>(mask,
                                  levelSetOf({{width, ys - eps},
                                              {a, ys - eps},
                                              {a, ys - eps + h},
                                              {width, ys - eps + h}}),
                                  ls::BooleanOperationEnum::UNION)
      .apply();
  auto top = levelSetOf({{-width, ys - eps + h},
                         {-a, ys - eps + h},
                         {-a, ys},
                         {a, ys},
                         {a, ys - eps + h},
                         {width, ys - eps + h}});
  // single steps with velocities given between them, as in ViennaPS
  ls::AdvectHermite<double, 2> advect(
      std::vector<ls::SmartPointer<ls::Domain<double, 2>>>{mask, top},
      ls::SmartPointer<ReversedEtch>::New(directional));
  advect.setSingleStep(true);
  advect.setTimeStepRatio(0.4999);
  double time = 0.;
  unsigned calls = 0;
  while (time < r - 1e-12 && calls < 100) {
    advect.setAdvectionTime(r - time);
    advect.apply();
    time += advect.getAdvectedTime();
    ++calls;
  }
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToDualContouringMesh<double, 2>(top, mesh).apply();
  double deviation = 0.;
  for (const auto &n : mesh->getNodes())
    deviation = std::max(
        deviation,
        directional ? openingTrenchDistance(n[0], n[1], a, ys, ys - eps + h, r)
                    : openingProfileDistance(n[0], n[1], a, ys, h, r));
  std::cout << "2D " << (directional ? "directional" : "isotropic")
            << " etching through a mask opening, mask below in "
               "the stack, single steps, grid spacing "
            << gridDelta << ", mask " << eps
            << " below the substrate: " << calls << " calls, worst vertex "
            << deviation / gridDelta << " grid spacings" << std::endl;
  if (mesh->getNodes().empty() ||
      deviation / gridDelta > (directional ? 0.05 : 0.3))
    return fail(2, "etching with the mask below in the stack off the exact "
                   "result");
  return EXIT_SUCCESS;
}

// distance to the arc of radius r around (cx, cy) from angle t0
// counterclockwise to t1
double arcDistance(double x, double y, double cx, double cy, double r,
                   double t0, double t1) {
  double t = std::atan2(y - cy, x - cx) - t0;
  t -= 2 * M_PI * std::floor(t / (2 * M_PI));
  if (t <= t1 - t0)
    return std::abs(std::hypot(x - cx, y - cy) - r);
  return std::min(
      std::hypot(x - cx - r * std::cos(t0), y - cy - r * std::sin(t0)),
      std::hypot(x - cx - r * std::cos(t1), y - cy - r * std::sin(t1)));
}

// distance to the exact boundary after isotropic etching by r of a substrate
// whose surface rises with y = ys + x tan(alpha), beside a mask |x| <= a with
// its top at yt: the mask, its underside laid open along the interface for r
// next to each foot, arcs of radius r around the feet and the slope moved
// down by r; without the mask the moved slope
double slopeDistance(double x, double y, double a, double ys, double yt,
                     double alpha, double r, double width, bool masked) {
  const double t = std::tan(alpha), c = std::cos(alpha), s = std::sin(alpha);
  if (!masked)
    return std::abs((y - ys - t * x) * c + r);
  const double ly = ys - t * a, ry = ys + t * a, far = 3. * width;
  double d = segmentDistance(x, y, -a, yt, a, yt);
  d = std::min(d, segmentDistance(x, y, -a, ly, -a, yt));
  d = std::min(d, segmentDistance(x, y, a, ry, a, yt));
  // the lower foot, where the wall leans over the slope
  d = std::min(d, segmentDistance(x, y, -a, ly, -a + r * c, ly + r * s));
  d = std::min(d, arcDistance(x, y, -a, ly, r, alpha - M_PI / 2, alpha));
  d = std::min(d, segmentDistance(x, y, -a + r * s, ly - r * c,
                                  -a + r * s - far * c, ly - r * c - far * s));
  // the higher foot
  d = std::min(d, segmentDistance(x, y, a, ry, a - r * c, ry - r * s));
  d = std::min(d,
               arcDistance(x, y, a, ry, r, M_PI + alpha, 1.5 * M_PI + alpha));
  return std::min(d,
                  segmentDistance(x, y, a + r * s, ry - r * c,
                                  a + r * s + far * c, ry - r * c + far * s));
}

// a mask on a slope, etched isotropically beside it: at the lower foot the
// wall leans over the slope and the substrate below the mask etches past the
// wall's normal up to the interface, at the higher foot the interface splits
// the cone into the mask's part, which stays, and the substrate's (0.71 and
// 0.87 grid spacings off at 20 and 30 degrees before). The end of the cavity
// below the mask lies in a cell crossed by the interface and the arc, and
// dual contouring places its vertex where the interface meets the arc's
// tangent at the other crossing, up to a third of a grid spacing beyond the
// end.
// Without the mask the moved slope is exact, also with the substrate below
// Manhattan, whose band holds values up to a grid spacing (3.3 grid spacings
// off before, where the value taken for a corner outside that band moved the
// interface below the slope and left parts of the slope static)
int checkSlope(double degrees, bool masked) {
  const double gridDelta = 0.5, r = 2.25 * gridDelta, width = 6.;
  const double a = 2.37, ys = 0.29, h = 1.13;
  const double alpha = degrees * M_PI / 180., t = std::tan(alpha);
  // the wall at the higher foot keeps the height h
  const double yt = ys + t * a + h;
  double bounds[4] = {-width, width, -8., 8.};
  ls::BoundaryConditionEnum boundaryCons[2] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
  // substrate below, Manhattan as usual
  auto substrate = ls::SmartPointer<ls::Domain<double, 2>>::New(
      bounds, boundaryCons, gridDelta);
  auto plane = ls::SmartPointer<ls::Mesh<>>::New();
  plane->insertNextNode({-width, ys - t * width, 0.});
  plane->insertNextNode({width, ys + t * width, 0.});
  plane->insertNextLine({0, 1});
  ls::FromSurfaceMesh<double, 2>(substrate, plane).apply();
  auto top = ls::SmartPointer<ls::Domain<double, 2>>::New(bounds, boundaryCons,
                                                          gridDelta);
  top->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
  top->setLevelSetWidth(12);
  auto surface = ls::SmartPointer<ls::Mesh<>>::New();
  std::vector<std::array<double, 2>> pts = {{-width, ys - t * width},
                                            {width, ys + t * width}};
  if (masked)
    pts = {{-width, ys - t * width}, {-a, ys - t * a},       {-a, yt}, {a, yt},
           {a, ys + t * a},          {width, ys + t * width}};
  for (const auto &q : pts)
    surface->insertNextNode({q[0], q[1], 0.});
  for (unsigned k = 0; k + 1 < pts.size(); ++k)
    surface->insertNextLine({k, k + 1});
  ls::FromSurfaceMesh<double, 2>(top, surface).apply();
  ls::AdvectHermite<double, 2> advect(
      std::vector<ls::SmartPointer<ls::Domain<double, 2>>>{substrate, top},
      ls::SmartPointer<MaskedEtch>::New(false));
  advect.setAdvectionTime(r);
  advect.apply();
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToDualContouringMesh<double, 2>(top, mesh).apply();
  double deviation = 0.;
  for (const auto &n : mesh->getNodes()) {
    if (std::abs(n[0]) > width - 2 * gridDelta)
      continue; // the slope meets the reflective boundaries at an angle
    deviation = std::max(deviation, slopeDistance(n[0], n[1], a, ys, yt, alpha,
                                                  r, width, masked));
  }
  const std::string name = "isotropic etching of a slope of " +
                           std::to_string(std::lround(degrees)) + " degrees " +
                           (masked ? "beside a mask" : "without a mask");
  std::cout << "2D " << name << ": " << advect.getNumberOfTimeSteps()
            << " steps, " << advect.getNumberOfRebuilds()
            << " rebuild(s), worst vertex " << deviation / gridDelta
            << " grid spacings" << std::endl;
  if (mesh->getNodes().empty() ||
      deviation / gridDelta > (masked ? 0.35 : 1e-6))
    return fail(2, name + " off the exact result");
  return EXIT_SUCCESS;
}

// relative noise in [-1, 1], fixed per cell of 0.05, seed and step
double stepNoise(const ls::Vec3D<double> &p, unsigned seed, unsigned step) {
  unsigned long long k = 1469598103934665603ULL ^ (seed * 2654435761ULL + step);
  for (int i = 0; i < 2; ++i) {
    k ^= static_cast<unsigned long long>(std::llround(p[i] / 0.05));
    k *= 1099511628211ULL;
  }
  k ^= k >> 33;
  k *= 0xff51afd7ed558ccdULL;
  k ^= k >> 33;
  return 2. * double(k % 1000003ULL) / 1000003. - 1.;
}

// etching of the substrate (material 1 with the mask below in the stack):
// downwards with the view factor of the opening, a little everywhere, and a
// noise that changes every step, as the flux of a ray-traced process
class NoisyRate : public ls::VelocityField<double> {
public:
  double a = 0., top = 0.;
  unsigned seed = 0, step = 0;
  double getScalarVelocity(const ls::Vec3D<double> &p, int material,
                           const ls::Vec3D<double> &normal,
                           unsigned long) override {
    if (material != 1)
      return 0.;
    const double d = top - p[1];
    const double view = (d <= 1e-9)
                            ? 1.
                            : 0.5 * (std::sin(std::atan((a - p[0]) / d)) +
                                     std::sin(std::atan((a + p[0]) / d)));
    return -(std::max(normal[1], 0.) * view + 0.03) *
           (1. + 0.3 * stepNoise(p, seed, step));
  }
};

// a 2D trench through a mask with the mask below in the stack, etched with a
// noisy rate in single steps at grid spacing 0.125. Next to the corner of the
// trench bottom two planes of the floor some 15 degrees apart put the vertex
// of a cell far beyond the wall of the mask, the mesh folded there, and the
// rebuild took the fold for a corner of the surface and etched from it into
// the substrate below the mask (0.11 and 0.22 beyond the mask edge for seeds
// 487 and 254 before, where the slow rate everywhere reaches 0.04)
int checkNoisyTrench(unsigned seed) {
  const double gridDelta = 0.125, width = 6., eps = 1e-4;
  const double a = 2.37, ys = 0.29, h = 1.13;
  double bounds[4] = {-width, width, -6., 6.};
  ls::BoundaryConditionEnum boundaryCons[2] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
  auto levelSetOf = [&](const std::vector<std::array<double, 2>> &points) {
    auto domain = ls::SmartPointer<ls::Domain<double, 2>>::New(
        bounds, boundaryCons, gridDelta);
    domain->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
    domain->setLevelSetWidth(12);
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    for (const auto &q : points)
      mesh->insertNextNode({q[0], q[1], 0.});
    for (unsigned k = 0; k + 1 < points.size(); ++k)
      mesh->insertNextLine({k, k + 1});
    ls::FromSurfaceMesh<double, 2>(domain, mesh).apply();
    return domain;
  };
  auto mask = levelSetOf({{-width, ys - eps + h},
                          {-a, ys - eps + h},
                          {-a, ys - eps},
                          {-width, ys - eps}});
  ls::BooleanOperation<double, 2>(mask,
                                  levelSetOf({{width, ys - eps},
                                              {a, ys - eps},
                                              {a, ys - eps + h},
                                              {width, ys - eps + h}}),
                                  ls::BooleanOperationEnum::UNION)
      .apply();
  auto top = levelSetOf({{-width, ys - eps + h},
                         {-a, ys - eps + h},
                         {-a, ys},
                         {a, ys},
                         {a, ys - eps + h},
                         {width, ys - eps + h}});
  auto velocity = ls::SmartPointer<NoisyRate>::New();
  velocity->a = a;
  velocity->top = ys - eps + h;
  velocity->seed = seed;
  ls::AdvectHermite<double, 2> advect(
      std::vector<ls::SmartPointer<ls::Domain<double, 2>>>{mask, top},
      velocity);
  advect.setSingleStep(true);
  advect.setTimeStepRatio(0.4999);
  double time = 0.;
  unsigned calls = 0;
  while (time < 1. - 1e-12 && calls < 100) {
    velocity->step = calls;
    advect.setAdvectionTime(1. - time);
    advect.apply();
    time += advect.getAdvectedTime();
    ++calls;
  }
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToDualContouringMesh<double, 2>(top, mesh).apply();
  double reach = 0.;
  for (const auto &n : mesh->getNodes())
    if (n[1] < ys - 0.02 && std::abs(n[0]) > a)
      reach = std::max(reach, std::abs(n[0]) - a);
  std::cout << "2D trench etched with a noisy rate, mask below in the stack, "
               "single steps, grid spacing 0.125, seed "
            << seed << ": " << calls << " calls, farthest vertex " << reach
            << " beyond the mask edge" << std::endl;
  if (mesh->getNodes().empty() || reach > 0.06)
    return fail(2, "trench etched with a noisy rate, seed " +
                       std::to_string(seed) + ", reaches below the mask");
  return EXIT_SUCCESS;
}

// a planar quad with corners in cyclic order as two triangles facing outward
void insertQuad(ls::Mesh<> &mesh, const std::array<unsigned, 4> &q,
                const Point &outward) {
  const auto &n = mesh.getNodes();
  const auto &a = n[q[0]], &b = n[q[1]], &c = n[q[2]];
  const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
  const double v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
  const double w[3] = {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                       u[0] * v[1] - u[1] * v[0]};
  if (w[0] * outward[0] + w[1] * outward[1] + w[2] * outward[2] < 0.) {
    mesh.insertNextTriangle({q[0], q[2], q[1]});
    mesh.insertNextTriangle({q[0], q[3], q[2]});
  } else {
    mesh.insertNextTriangle({q[0], q[1], q[2]});
    mesh.insertNextTriangle({q[0], q[2], q[3]});
  }
}

// distance to the exact boundary after etching by r through a square opening
// |x|, |z| < a of a mask (top yt = ys + h) on the substrate surface y = ys:
// a pit with vertical walls (directional), or a pit of the points within r
// of the opening, with quarter cylinders along its edges and sphere octants
// at its corners, undercutting the mask (isotropic)
double openingDistance(const Point &p, double a, double ys, double h, double r,
                       bool directional) {
  const double x = std::abs(p[0]), y = p[1], z = std::abs(p[2]);
  const double yt = ys + h, ox = x - a, oz = z - a;
  // signed distance in the plane to the opening, negative inside
  const double sd = (ox > 0. || oz > 0.)
                        ? std::hypot(std::max(ox, 0.), std::max(oz, 0.))
                        : std::max(ox, oz);
  double d = std::hypot(y - yt, std::max(-sd, 0.)); // mask top
  // walls of the opening, down to the bottom of the pit when directional
  const double bottom = directional ? ys - r : ys;
  const double dy = std::max({bottom - y, y - yt, 0.});
  d = std::min(d, std::hypot(ox, std::max(oz, 0.), dy));
  d = std::min(d, std::hypot(oz, std::max(ox, 0.), dy));
  if (directional)
    return std::min(d, std::hypot(y - (ys - r), std::max(sd, 0.)));
  // underside of the mask over the undercut, 0 <= sd <= r
  d = std::min(d, std::hypot(y - ys, std::max({-sd, sd - r, 0.})));
  // pit surface; from above, its nearest points lie on the rim or straight
  // below on the flat bottom
  double pit;
  if (y <= ys) {
    pit = std::abs(std::hypot(y - ys, std::max(sd, 0.)) - r);
  } else {
    pit = std::hypot(y - ys, sd - r);
    if (sd <= 0.)
      pit = std::min(pit, y - ys + r);
  }
  return std::min(d, pit);
}

int checkOpening() {
  const double gridDelta = 1., r = 2.25 * gridDelta, width = 8.;
  const double a = 3.37, ys = 0.29, h = 2.13;
  for (bool directional : {true, false}) {
    double bounds[6] = {-width, width, -6., 6., -width, width};
    ls::BoundaryConditionEnum boundaryCons[3] = {
        ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
        ls::BoundaryConditionEnum::INFINITE_BOUNDARY,
        ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY};
    // substrate below, Manhattan as usual
    auto substrate = ls::SmartPointer<ls::Domain<double, 3>>::New(
        bounds, boundaryCons, gridDelta);
    auto plane = ls::SmartPointer<ls::Mesh<>>::New();
    for (const auto &q : {Point{-width, ys, -width}, Point{width, ys, -width},
                          Point{width, ys, width}, Point{-width, ys, width}})
      plane->insertNextNode(q);
    insertQuad(*plane, {0, 1, 2, 3}, {0., 1., 0.});
    ls::FromSurfaceMesh<double, 3>(substrate, plane).apply();
    // surface of substrate and mask, Euclidean with Hermite data: mask top
    // (outer corners 0 to 3, inner 4 to 7), walls of the opening and its
    // floor (corners 8 to 11)
    auto top = ls::SmartPointer<ls::Domain<double, 3>>::New(
        bounds, boundaryCons, gridDelta);
    top->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
    top->setLevelSetWidth(12);
    auto surface = ls::SmartPointer<ls::Mesh<>>::New();
    const double sx[4] = {-1., 1., 1., -1.}, sz[4] = {-1., -1., 1., 1.};
    for (const auto &[size, y] :
         {std::pair{width, ys + h}, std::pair{a, ys + h}, std::pair{a, ys}})
      for (int k = 0; k < 4; ++k)
        surface->insertNextNode({sx[k] * size, y, sz[k] * size});
    for (unsigned k = 0; k < 4; ++k) {
      const unsigned l = (k + 1) % 4;
      insertQuad(*surface, {k, l, 4 + l, 4 + k}, {0., 1., 0.});
      insertQuad(*surface, {4 + k, 4 + l, 8 + l, 8 + k},
                 {-(sx[k] + sx[l]), 0., -(sz[k] + sz[l])});
    }
    insertQuad(*surface, {8, 9, 10, 11}, {0., 1., 0.});
    ls::FromSurfaceMesh<double, 3>(top, surface).apply();

    ls::AdvectHermite<double, 3> advect(
        std::vector<ls::SmartPointer<ls::Domain<double, 3>>>{substrate, top},
        ls::SmartPointer<MaskedEtch>::New(directional));
    advect.setAdvectionTime(r);
    advect.apply();
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToDualContouringMesh<double, 3>(top, mesh).apply();
    double deviation = 0.;
    for (const auto &n : mesh->getNodes()) {
      if (std::abs(n[0]) > width - gridDelta ||
          std::abs(n[2]) > width - gridDelta)
        continue; // the reflective boundaries
      deviation = std::max(deviation, openingDistance({n[0], n[1], n[2]}, a, ys,
                                                      h, r, directional));
    }
    const std::string name =
        directional ? "directional etching" : "isotropic etching";
    std::cout << "3D " << name
              << " through a square opening: " << advect.getNumberOfTimeSteps()
              << " steps, " << advect.getNumberOfRebuilds()
              << " rebuild(s), worst vertex " << deviation / gridDelta
              << " grid spacings" << std::endl;
    if (mesh->getNodes().empty() ||
        deviation / gridDelta > (directional ? 0.05 : 0.3))
      return fail(3, name + " through a square opening off the exact result");
  }
  return EXIT_SUCCESS;
}

// a trench through a mask in 3D, z up, both level sets Euclidean, etched in
// single steps with the view factor of the opening and a reproducible noise
// of 10 % or 30 %, like the flux of a ViennaPS process; the band must end as
// a valid level set. A moving concave fan at a corner of the trench bottom
// had offered its plane to nodes below the floor beside it, libigl in single
// precision had taken the closest point of a node next to an edge on the
// wrong element, with its sign, and with 30 % the rebuild had left folded
// quads out of its mesh together with the floor beneath them, so that the
// nodes above and below the hole took their distances from its rim. With the
// mask below in the stack and the time step of ViennaPS the last rebuild
// leaves the band, and dual contouring had put vertices into the mask above
// the thin tip of the cavity under it: an edge crossed by the ceiling of the
// cavity took the plane of the floor at its interpolated crossing
double trenchViewFactor(double x, double z, double a, double top) {
  const double h = top - z;
  if (h <= 1e-9)
    return 1.;
  return 0.5 *
         (std::sin(std::atan((a - x) / h)) + std::sin(std::atan((a + x) / h)));
}

// relative noise in [-1, 1], fixed per cell of 0.05 and per seed
double trenchNoise(const ls::Vec3D<double> &p, unsigned seed) {
  unsigned long long h = 1469598103934665603ULL ^ seed;
  for (int i = 0; i < 3; ++i) {
    h ^= static_cast<unsigned long long>(std::llround(p[i] / 0.05));
    h *= 1099511628211ULL;
  }
  h ^= h >> 33;
  h *= 0xff51afd7ed558ccdULL;
  h ^= h >> 33;
  return 2. * double(h % 1000003ULL) / 1000003. - 1.;
}

template <class T> class TrenchEtch : public ls::VelocityField<T> {
public:
  double a = 0., top = 0., noise = 0.1;
  unsigned seed = 0;
  int etched = 0; // the material that etches
  T getScalarVelocity(const ls::Vec3D<T> &p, int material, const ls::Vec3D<T> &,
                      unsigned long) override {
    if (material != etched)
      return T(0);
    const ls::Vec3D<double> q{double(p[0]), double(p[1]), double(p[2])};
    return T(-trenchViewFactor(q[0], q[2], a, top) *
             (1. + noise * trenchNoise(q, seed)));
  }
};

template <class T>
int checkTrench(unsigned seed, double noise = 0.1, bool maskBelow = false) {
  // with the mask below in the stack it reaches eps below the substrate
  // surface, as MakeTrench builds it
  const double gridDelta = 0.5, X = 6., Y = 1.5, a = 2.37, base = 0.29,
               eps = maskBelow ? 1e-4 : 0., top = 1.42 - eps;
  auto newLevelSet = [&]() {
    double bounds[6] = {-X, X, -Y, Y, -5., 5.};
    ls::BoundaryConditionEnum boundaryCons[3] = {
        ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
        ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
        ls::BoundaryConditionEnum::INFINITE_BOUNDARY};
    auto levelSet = ls::SmartPointer<ls::Domain<T, 3>>::New(
        bounds, boundaryCons, T(gridDelta));
    levelSet->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
    levelSet->setLevelSetWidth(12);
    return levelSet;
  };
  auto convert = [](const ls::Mesh<> &mesh) {
    auto result = ls::SmartPointer<ls::Mesh<T>>::New();
    for (const auto &n : mesh.getNodes())
      result->insertNextNode({T(n[0]), T(n[1]), T(n[2])});
    for (const auto &e : mesh.triangles)
      result->insertNextTriangle(e);
    return result;
  };
  // below in the stack the substrate below the plane z = base, or the mask:
  // slabs |x| >= a from eps below the substrate surface to the top, open at
  // the domain boundaries
  ls::SmartPointer<ls::Domain<T, 3>> lower;
  if (!maskBelow) {
    ls::Mesh<> plane;
    for (const auto &q : {Point{-X, -Y, base}, Point{X, -Y, base},
                          Point{X, Y, base}, Point{-X, Y, base}})
      plane.insertNextNode(q);
    insertQuad(plane, {0, 1, 2, 3}, {0., 0., 1.});
    lower = newLevelSet();
    ls::FromSurfaceMesh<T, 3>(lower, convert(plane)).apply();
  } else {
    for (int side : {-1, 1}) {
      ls::Mesh<> slab;
      const double xo = side * X, xi = side * a, z0 = base - eps;
      const Point corners[8] = {{xo, -Y, top}, {xi, -Y, top}, {xi, Y, top},
                                {xo, Y, top},  {xi, -Y, z0},  {xo, -Y, z0},
                                {xo, Y, z0},   {xi, Y, z0}};
      for (const auto &q : corners)
        slab.insertNextNode(q);
      insertQuad(slab, {0, 1, 2, 3}, {0., 0., 1.});
      insertQuad(slab, {1, 4, 7, 2}, {-double(side), 0., 0.});
      insertQuad(slab, {4, 5, 6, 7}, {0., 0., -1.});
      auto part = newLevelSet();
      ls::FromSurfaceMesh<T, 3>(part, convert(slab)).apply();
      if (!lower)
        lower = part;
      else
        ls::BooleanOperation<T, 3>(lower, part, ls::BooleanOperationEnum::UNION)
            .apply();
    }
  }
  // mask tops, walls of the opening |x| < a and its floor
  ls::Mesh<> profile;
  const double xs[6] = {-X, -a, -a, a, a, X};
  const double zs[6] = {top, top, base, base, top, top};
  for (int k = 0; k < 6; ++k)
    for (double y : {-Y, Y})
      profile.insertNextNode({xs[k], y, zs[k]});
  const Point outward[5] = {
      {0., 0., 1.}, {1., 0., 0.}, {0., 0., 1.}, {-1., 0., 0.}, {0., 0., 1.}};
  for (unsigned k = 0; k < 5; ++k)
    insertQuad(profile, {2 * k, 2 * k + 1, 2 * k + 3, 2 * k + 2}, outward[k]);
  auto surface = newLevelSet();
  ls::FromSurfaceMesh<T, 3>(surface, convert(profile)).apply();

  auto velocity = ls::SmartPointer<TrenchEtch<T>>::New();
  velocity->a = a;
  velocity->top = top;
  velocity->seed = seed;
  velocity->noise = noise;
  velocity->etched = maskBelow ? 1 : 0;
  ls::AdvectHermite<T, 3> advect(
      std::vector<ls::SmartPointer<ls::Domain<T, 3>>>{lower, surface},
      velocity);
  advect.setSingleStep(true);
  // with the mask below, the time step of ViennaPS, after which the last
  // rebuild leaves the band
  if (maskBelow)
    advect.setTimeStepRatio(T(0.4999));
  const double time = 2.25 * gridDelta;
  double done = 0.;
  unsigned calls = 0;
  while (done < time - 1e-9 && calls < 100) {
    advect.setAdvectionTime(T(time - done));
    advect.apply();
    done += advect.getAdvectedTime();
    ++calls;
  }
  ls::Check<T, 3> check(surface);
  check.apply();
  auto mesh = ls::SmartPointer<ls::Mesh<T>>::New();
  if (check.isValid())
    ls::ToDualContouringMesh<T, 3>(surface, mesh).apply();
  double lowest = 1e300;
  // vertices inside the body of the mask, more than 0.08 beside its wall,
  // 0.06 above its bottom and 0.07 below its top
  int inMask = 0;
  for (const auto &n : mesh->getNodes()) {
    lowest = std::min(lowest, double(n[2]));
    const double x = std::abs(double(n[0])), z = double(n[2]);
    if (x > a + 0.08 && z > base + 0.06 && z < top - 0.07)
      ++inMask;
  }
  std::string name = std::string(sizeof(T) == 4 ? "single" : "double") +
                     " precision, seed " + std::to_string(seed);
  if (noise != 0.1)
    name += ", noise " + std::to_string(std::lround(100. * noise)) + " %";
  if (maskBelow)
    name += ", mask below";
  std::cout << "3D trench etched with a noisy view factor, " << name << ": "
            << calls << " calls, valid " << check.isValid()
            << ", deepest vertex " << lowest << std::endl;
  // the floor sinks by at most the largest speed times the time
  if (!check.isValid() || mesh->getNodes().empty() ||
      lowest < base - (1. + noise) * time)
    return fail(3, "trench etched with a noisy view factor, " + name);
  if (inMask > 0)
    return fail(3, "trench etched with a noisy view factor, " + name + ": " +
                       std::to_string(inMask) + " vertices in the mask");
  return EXIT_SUCCESS;
}

// the rotated boxes of check() in single precision, one call each: the exact
// cases exact to the rounding of float, the others within the bounds of the
// double runs
template <class T> class Kinds : public ls::VelocityField<T> {
public:
  int kind = 0; // etching, deposition, directional etching and deposition
  T getScalarVelocity(const ls::Vec3D<T> &, int, const ls::Vec3D<T> &normal,
                      unsigned long) override {
    if (kind < 2)
      return (kind == 0) ? T(-1) : T(1);
    return (kind == 2 ? T(-1) : T(1)) * std::max(normal[1], T(0));
  }
};

template <int D> int checkSinglePrecision(double gridDelta) {
  const RotatedBox<D> box(gridDelta);
  const double r = 2.25 * gridDelta;
  const std::string names[4] = {"etching", "deposition", "directional etching",
                                "directional deposition"};
  const double tolerances[4] = {1e-4, 0.3, 1e-4, D == 2 ? 0.01 : 0.1};
  for (int kind = 0; kind < 4; ++kind) {
    double bounds[2 * D];
    ls::BoundaryConditionEnum boundaryCons[D];
    for (int i = 0; i < D; ++i) {
      bounds[2 * i] = -18. * gridDelta;
      bounds[2 * i + 1] = 18. * gridDelta;
      boundaryCons[i] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
    }
    auto domain = ls::SmartPointer<ls::Domain<float, D>>::New(
        bounds, boundaryCons, float(gridDelta));
    domain->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
    domain->setLevelSetWidth(12);
    auto boxMesh = box.mesh();
    auto mesh = ls::SmartPointer<ls::Mesh<float>>::New();
    for (const auto &n : boxMesh->getNodes())
      mesh->insertNextNode({float(n[0]), float(n[1]), float(n[2])});
    for (const auto &e : boxMesh->template getElements<D>())
      mesh->insertNextElement(e);
    ls::FromSurfaceMesh<float, D> mesher(domain, mesh);
    mesher.setUseWindingNumber(true);
    mesher.apply();
    auto velocity = ls::SmartPointer<Kinds<float>>::New();
    velocity->kind = kind;
    ls::AdvectHermite<float, D> advect(domain, velocity);
    advect.setAdvectionTime(float(r));
    advect.apply();
    auto surface = ls::SmartPointer<ls::Mesh<float>>::New();
    ls::ToDualContouringMesh<float, D>(domain, surface).apply();
    if (surface->getNodes().empty())
      return fail(D, names[kind] + " in single precision: no surface");
    double deviation = 0.;
    for (const auto &n : surface->getNodes()) {
      const Point p{n[0], n[1], n[2]};
      double exact;
      if (kind == 0)
        exact = box.sdf(p, r);
      else if (kind == 1)
        exact = box.sdf(p) - r;
      else if (kind == 2)
        exact = std::max(box.sdf(p), box.sdf({p[0], p[1] + r, p[2]}));
      else
        exact = box.sweptSdf(p, r);
      deviation = std::max(deviation, std::abs(exact));
    }
    std::cout << D << "D " << names[kind]
              << " in single precision: worst vertex " << deviation / gridDelta
              << " grid spacings" << std::endl;
    if (deviation / gridDelta > tolerances[kind])
      return fail(D, names[kind] + " in single precision off");
  }
  return EXIT_SUCCESS;
}

int main() {
  omp_set_num_threads(4);
  if (check<2>(0.5) != EXIT_SUCCESS || check<3>(1.) != EXIT_SUCCESS ||
      checkPeriodic() != EXIT_SUCCESS || checkMaterials() != EXIT_SUCCESS ||
      checkReversedStack(0.5, 1e-4) != EXIT_SUCCESS ||
      checkReversedStack(0.125, 1e-4) != EXIT_SUCCESS ||
      checkReversedStack(0.5, 1e-2) != EXIT_SUCCESS ||
      checkReversedStack(0.5, 1e-4, true) != EXIT_SUCCESS ||
      checkReversedStack(0.5, 1e-2, true) != EXIT_SUCCESS ||
      checkSlope(10., true) != EXIT_SUCCESS ||
      checkSlope(20., true) != EXIT_SUCCESS ||
      checkSlope(30., true) != EXIT_SUCCESS ||
      checkSlope(20., false) != EXIT_SUCCESS ||
      checkSlope(30., false) != EXIT_SUCCESS ||
      checkNoisyTrench(487) != EXIT_SUCCESS ||
      checkNoisyTrench(254) != EXIT_SUCCESS || checkOpening() != EXIT_SUCCESS ||
      checkSinglePrecision<2>(0.5) != EXIT_SUCCESS ||
      checkSinglePrecision<3>(1.) != EXIT_SUCCESS ||
      checkTrench<double>(0) != EXIT_SUCCESS ||
      checkTrench<double>(1) != EXIT_SUCCESS ||
      checkTrench<float>(0) != EXIT_SUCCESS ||
      checkTrench<float>(1) != EXIT_SUCCESS ||
      checkTrench<double>(533, 0.3) != EXIT_SUCCESS ||
      checkTrench<float>(533, 0.3) != EXIT_SUCCESS ||
      checkTrench<double>(0, 0.3, true) != EXIT_SUCCESS ||
      checkTrench<float>(0, 0.3, true) != EXIT_SUCCESS)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
