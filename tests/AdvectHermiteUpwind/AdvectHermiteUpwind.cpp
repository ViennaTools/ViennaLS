#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
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
  third-order Runge-Kutta. \example AdvectHermiteUpwind.cpp
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

int main() {
  const double gridDelta = 0.25, r = 1., width = 6.;
  const double a = 2.37, ys = 0.29, h = 1.13;
  struct Case {
    std::string name;
    ls::SpatialSchemeEnum scheme;
    ls::TemporalSchemeEnum time;
    // the worst vertex of the trench and of the mask, grid spacings: the
    // first-order upwind update rounds the trench's bottom corners (0.47,
    // 0.67 with Runge-Kutta 3, 0.41 with Lax-Friedrichs), the second-order
    // and WENO schemes bulge the lower wall below the mask (0.47, 0.53) and
    // move the mask's vertices by up to 0.01 through their wider stencils
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
       ls::TemporalSchemeEnum::FORWARD_EULER, 0.75, 0.02},
      {"Lax-Friedrichs 2nd order",
       ls::SpatialSchemeEnum::LAX_FRIEDRICHS_2ND_ORDER,
       ls::TemporalSchemeEnum::FORWARD_EULER, 0.75, 1e-3},
      {"WENO 5", ls::SpatialSchemeEnum::WENO_5TH_ORDER,
       ls::TemporalSchemeEnum::RUNGE_KUTTA_3RD_ORDER, 0.75, 0.02}};
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
  return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
