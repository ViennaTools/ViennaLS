// Stress-dependent oxidation of a Si cylinder (Kao) against the exact 1-D
// radial model: the radial source flow of an Eyring oxide in the annulus, the
// reaction law on the interface normal stress (V_k) and the diffusivity law on
// the pressure (V_D), each law alone.  Reference radii from
// patches/analysis/kao/kao1d.py (normal-stress law, eta0 3e9 Pa hr, Vc 300 A^3,
// 30 nm seed, run_cyl(..., 0.3, "convex", 0.15)).  The oxide is purely
// viscous here, like the reference: a Maxwell oxide (G 30 GPa, relaxation
// 0.1 hr in this 0.15 hr run) builds up its stress over the relaxation time,
// since new oxide is born stress-free (V_D alone: 0.9 nm more Si consumed).
// The elastic run had also hidden the third defect below.
//
// Three defects made the 2-D cylinder miss these, at any grid spacing:
// - the stress laws read the MAC stress of the first node layer at the Si
//   front, which is O(1) wrong there (V_k alone: 1.2 nm too much Si consumed);
// - each diffusion node used its own D on its faces, not conservative with a
//   pressure-dependent D (V_D alone: 2.9 nm too little Si consumed);
// - the reaction rate and front speed were per-node values copied to the
//   front, not evaluated at its interface points (V_k alone: 0.8 nm too
//   little Si consumed; 0.09 nm off with the interface points).
#include <lsGeometricAdvect.hpp>
#include <lsGeometricAdvectDistributions.hpp>
#include <lsMakeGeometry.hpp>
#include <lsOxidation.hpp>
#include <lsTestAsserts.hpp>
#include <lsToSurfaceMesh.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>

#include <omp.h>

namespace ls = viennals;
using T = double;
constexpr int D = 2;

T meanRadius(ls::SmartPointer<ls::Domain<T, D>> levelSet) {
  auto mesh = ls::SmartPointer<ls::Mesh<T>>::New();
  ls::ToSurfaceMesh<T, D>(levelSet, mesh).apply();
  T sum = 0.;
  for (const auto &p : mesh->getNodes())
    sum += std::hypot(p[0], p[1]);
  return sum / mesh->getNodes().size();
}

// Si radius after oxidising a cylinder of radius a0 for tEnd
T siRadius(T Vk, T Vd) {
  const T a0 = 0.3, h = 0.01, tEnd = 0.15, L = a0 + 1.0;
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

  // wet 1000 C, isotropic kinetics
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
  ox.reactionActivationVolume = Vk * 1e-30;
  ox.diffusionActivationVolume = Vd * 1e-30;
  ls::OxidationCouplingParameters coupling;
  // a viscous oxide's pressure does not shrink with the step, so halving it
  // does not help the pressure-concentration coupling: some steps need more
  // than 30 iterations
  coupling.maxIterations = 100;
  coupling.tolerance = 5e-2;
  auto model = ls::Oxidation<T, D>::New(si, ambient);
  model->setOxidationParameters(ox);
  model->setCouplingParameters(coupling);

  T t = 0., dt = 0.004;
  while (tEnd - t > 1e-9) {
    const T requested = std::min(dt, tEnd - t);
    ls::OxidationDeformationParameters p;
    p.viscosity = 3.0e9; // G = 0: purely viscous
    p.shearActivationVolume = 300e-30;
    p.stressTimeStep = requested;
    model->setDeformationParameters(p);
    const T actual = model->applyCFLLimited(requested, 0.499);
    t += actual;
    const T vmax = model->getLastMaxVelocity();
    dt = std::min(vmax > 0 ? 0.499 * h / vmax : requested, 2 * actual);
  }
  return meanRadius(si);
}

int main() {
  omp_set_num_threads(2);
  ls::Logger::setLogLevel(ls::LogLevel::WARNING);
  const T reactionOnly = siRadius(25., 0.);
  const T diffusivityOnly = siRadius(0., 40.);
  std::printf("Si radius after 0.15 hr: V_k only %.5f um (1-D 0.26916), "
              "V_D only %.5f um (1-D 0.25640)\n",
              reactionOnly, diffusivityOnly);
  std::fflush(stdout);
  VC_TEST_ASSERT(std::abs(reactionOnly - 0.26916) < 8e-4)
  VC_TEST_ASSERT(std::abs(diffusivityOnly - 0.25640) < 8e-4)
  return 0;
}
