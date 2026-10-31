// Oxide growth around a Si disk must not depend on the grid orientation.
//
// With isotropic kinetics and a linear viscous oxide the exact flow is the
// radial source flow v = c / r, the Si front and the free surface stay
// round and centred.  Two grid effects broke this:
// - the Si-front normal that orients the expansion velocity came from a
//   central difference of phi clamped to one grid unit, which turned the
//   normals of tilted fronts towards the diagonal by up to 7 degrees at any
//   grid spacing: the oxide flowed towards the diagonals (tangential
//   velocity 9 % of the radial one at the Si front, radial velocity +-6 % at
//   the free surface) and bulged there;
// - the gas-side velocity of the free surface came from the FIRST of
//   equidistant nearest oxide nodes in scan order (-y first), which shifted
//   the free surface in y.
#include <lsGeometricAdvect.hpp>
#include <lsGeometricAdvectDistributions.hpp>
#include <lsMakeGeometry.hpp>
#include <lsOxidation.hpp>
#include <lsTestAsserts.hpp>
#include <lsToSurfaceMesh.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include <omp.h>

namespace ls = viennals;
using T = double;
constexpr int D = 2;

struct Shape {
  T cx, cy, mean, spread;
};

Shape shapeOf(ls::SmartPointer<ls::Domain<T, D>> levelSet) {
  auto mesh = ls::SmartPointer<ls::Mesh<T>>::New();
  ls::ToSurfaceMesh<T, D>(levelSet, mesh).apply();
  const auto &nodes = mesh->getNodes();
  Shape s{0., 0., 0., 0.};
  for (const auto &p : nodes) {
    s.cx += p[0];
    s.cy += p[1];
  }
  s.cx /= nodes.size();
  s.cy /= nodes.size();
  T lo = 1e9, hi = -1e9;
  for (const auto &p : nodes) {
    const T r = std::hypot(p[0], p[1]);
    s.mean += r;
    lo = std::min(lo, r);
    hi = std::max(hi, r);
  }
  s.mean /= nodes.size();
  s.spread = hi - lo;
  return s;
}

int main() {
  omp_set_num_threads(2);
  ls::Logger::setLogLevel(ls::LogLevel::WARNING);

  const T a0 = 0.3, h = 0.01, L = a0 + 0.6;
  double bounds[2 * D] = {-L, L, -L, L};
  ls::Domain<T, D>::BoundaryType bc[D] = {
      ls::Domain<T, D>::BoundaryType::INFINITE_BOUNDARY,
      ls::Domain<T, D>::BoundaryType::INFINITE_BOUNDARY};
  auto si = ls::SmartPointer<ls::Domain<T, D>>::New(bounds, bc, h);
  ls::MakeGeometry<T, D>(
      si, ls::Sphere<T, D>::New(ls::VectorType<T, D>{0., 0.}, a0))
      .apply();
  auto ambient = ls::SmartPointer<ls::Domain<T, D>>::New(si);
  ls::GeometricAdvect<T, D>(
      ambient,
      ls::SmartPointer<ls::SphereDistribution<viennahrle::CoordType, D>>::New(
          3 * h))
      .apply();

  // wet 1000 C Deal-Grove rates, isotropic, no stress kinetics
  ls::OxidationParameters ox;
  ox.diffusionCoefficient = 0.315463 / 2;
  ox.reactionRate = 0.744430;
  ox.transferCoefficient = 100.;
  ox.equilibriumConcentration = 1.;
  ox.oxidantMoleculeDensity = 1.;
  ox.expansionCoefficient = 2.27;
  ox.velocitySign = -1.;
  ox.temperature = 1273.15;
  ox.reactionRateRatio111 = 1.;
  auto model = ls::Oxidation<T, D>::New(si, ambient);
  model->setOxidationParameters(ox);

  T t = 0., dt = 0.004;
  const T tEnd = 0.3;
  while (tEnd - t > 1e-9) {
    const T requested = std::min(dt, tEnd - t);
    ls::OxidationDeformationParameters p;
    p.viscosity = 3.0e9; // linear viscous oxide
    p.stressTimeStep = requested;
    model->setDeformationParameters(p);
    const T actual = model->applyCFLLimited(requested, 0.499);
    t += actual;
    const T vmax = model->getLastMaxVelocity();
    dt = std::min(vmax > 0 ? 0.499 * h / vmax : requested, 2 * actual);
  }

  const auto siShape = shapeOf(si);
  const auto oxShape = shapeOf(ambient);
  std::printf("Si: centroid (%+.3f, %+.3f) nm, spread %.2f nm; "
              "oxide: centroid (%+.3f, %+.3f) nm, spread %.2f nm\n",
              siShape.cx * 1e3, siShape.cy * 1e3, siShape.spread * 1e3,
              oxShape.cx * 1e3, oxShape.cy * 1e3, oxShape.spread * 1e3);
  std::fflush(stdout);

  // velocity of the last solve on rings near the Si front and near the
  // free surface, every 5 degrees
  auto field = model->getDeformationField();
  auto ring = [&](T fraction, T &tangential, T &radialSpread) {
    const T r = siShape.mean + fraction * (oxShape.mean - siShape.mean);
    std::vector<T> vr, vt;
    for (int k = 0; k < 72; ++k) {
      const T a = k * 5. * M_PI / 180.;
      const auto v = field->getVelocity(
          ls::Vec3D<T>{r * std::cos(a), r * std::sin(a), 0.});
      vr.push_back(v[0] * std::cos(a) + v[1] * std::sin(a));
      vt.push_back(-v[0] * std::sin(a) + v[1] * std::cos(a));
    }
    T mean = 0.;
    for (T x : vr)
      mean += x;
    mean /= vr.size();
    tangential = 0.;
    for (T x : vt)
      tangential = std::max(tangential, std::abs(x) / mean);
    radialSpread = (*std::max_element(vr.begin(), vr.end()) -
                    *std::min_element(vr.begin(), vr.end())) /
                   mean;
  };
  T tNear, sNear, tFar, sFar;
  ring(0.1, tNear, sNear);
  ring(0.8, tFar, sFar);
  std::printf("near the Si front: max |v_theta| / v_r %.4f, v_r spread %.4f; "
              "near the free surface: %.4f, %.4f\n",
              tNear, sNear, tFar, sFar);
  std::fflush(stdout);

  // Measured here: clamped normals 0.066 near the Si front and a radial
  // spread of 0.13 near the free surface (0.019 and 0.031 after that fix,
  // 0.0070 and 0.0099 with the front evaluated at its interface points);
  // with the first-found nearest node the oxide centroid sits at -4.05 nm
  // in y.
  VC_TEST_ASSERT(tNear < 0.03)
  VC_TEST_ASSERT(sFar < 0.04)
  VC_TEST_ASSERT(std::abs(oxShape.cx) < 3e-4 && std::abs(oxShape.cy) < 3e-4)
  VC_TEST_ASSERT(std::abs(siShape.cx) < 1e-4 && std::abs(siShape.cy) < 1e-4)
  return 0;
}
