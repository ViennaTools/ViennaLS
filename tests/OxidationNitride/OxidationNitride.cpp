// Maxwell nitride of the monolithic oxide + nitride mechanics against exact
// solutions (Senez et al., IEEE TED 43, 720 (1996): the deviatoric stress
// relaxes over tau = eta_N / mu, the mean stress stays elastic).
//
// 1. Relaxation: a nitride film on a pad oxide carrying a uniform deviatoric
//    in-plane pre-stress (sigma_xx = s0, sigma_zz = -s0, sigma_yy = 0) is in
//    equilibrium and keeps its volume, so nothing moves and the stress decays
//    as exp(-t / tau).  The Eyring nitride (shear-thinned viscosity) relaxes
//    faster.
// 2. Elastic limit: with a huge viscosity the nitride is the elastic one, in
//    a case that moves (pre-stress in half of the film).
// 3. Preset: the Si3N4 preset reproduces Senez's Table II.
#include <lsDomain.hpp>
#include <lsInterior.hpp>
#include <lsMakeGeometry.hpp>
#include <lsOxidationDeformation.hpp>
#include <lsOxidationDiffusion.hpp>
#include <lsOxidationPresets.hpp>
#include <lsTestAsserts.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>

#include <omp.h>

namespace ls = viennals;

using T = double;
constexpr int D = 2;
using LevelSet = ls::SmartPointer<ls::Domain<T, D>>;

constexpr T gridDelta = 0.005, padTop = 0.04, nitrideTop = 0.08, s0 = 1e9;

ls::Domain<T, D>::BoundaryType boundaries[D] = {
    ls::Domain<T, D>::BoundaryType::REFLECTIVE_BOUNDARY,
    ls::Domain<T, D>::BoundaryType::INFINITE_BOUNDARY};
double bounds[2 * D] = {-0.05, 0.05, -0.03, 0.13};

LevelSet makePlane(T height) {
  auto levelSet = LevelSet::New(bounds, boundaries, gridDelta);
  ls::MakeGeometry<T, D>(levelSet,
                         ls::Plane<T, D>::New(ls::VectorType<T, D>{0., height},
                                              ls::VectorType<T, D>{0., 1.}))
      .apply();
  return levelSet;
}

// Fill the mask and store the nitride stress of every point ("NitStress",
// the step-start stress the mechanics reads).
void storeNitrideStress(
    LevelSet mask,
    const std::function<std::array<T, 9>(const viennahrle::Index<D> &)>
        &stress) {
  ls::Interior<T, D>(mask).apply();
  using VectorData = typename ls::PointData<T>::VectorDataType;
  VectorData rows[3];
  viennahrle::ConstSparseIterator<typename ls::Domain<T, D>::DomainType> it(
      mask->getDomain());
  for (; !it.isFinished(); ++it)
    if (it.isDefined()) {
      const auto s = stress(it.getStartIndices());
      for (unsigned r = 0; r < 3; ++r)
        rows[r].push_back({s[3 * r], s[3 * r + 1], s[3 * r + 2]});
    }
  for (unsigned r = 0; r < 3; ++r)
    mask->getPointData().insertReplaceVectorData(std::move(rows[r]),
                                                 "NitStressR" +
                                                     std::to_string(r));
}

// One mechanics solve over dt of the Si / pad oxide / nitride stack with the
// given nitride pre-stress and creep (viscosity 0 = elastic).
ls::SmartPointer<ls::OxidationDeformation<T, D>>
solveStack(const std::function<std::array<T, 9>(const viennahrle::Index<D> &)>
               &stress,
           T viscosity, T activationVolume, T dt) {
  auto si = makePlane(0.);
  auto ambient = makePlane(nitrideTop);
  auto mask = LevelSet::New(bounds, boundaries, gridDelta);
  double lo[D] = {-1., padTop}, hi[D] = {1., 1.};
  ls::MakeGeometry<T, D>(mask, ls::Box<T, D>::New(lo, hi)).apply();
  storeNitrideStress(mask, stress);

  ls::OxidationParameters oxidation; // no reaction: nothing grows
  oxidation.reactionRate = 0.;
  oxidation.velocitySign = -1.;
  oxidation.temperature = 1273.15;
  ls::OxidationDeformationParameters mechanics;
  mechanics.viscosity = 1e12;
  mechanics.shearModulus = 3e10;
  mechanics.stressTimeStep = dt;
  const auto nitride = ls::OxidationPresets::siliconNitrideMask1000C();
  mechanics.nitrideYoungModulus = nitride.youngModulus;
  mechanics.nitridePoissonRatio = nitride.poissonRatio;
  mechanics.nitrideViscosity = viscosity;
  mechanics.nitrideShearActivationVolume = activationVolume;
  auto diffusion = ls::OxidationDiffusion<T, D>::New(si, ambient, oxidation);
  diffusion->setMaskInterface(mask, -1);
  diffusion->apply();
  auto deformation = ls::OxidationDeformation<T, D>::New(
      si, ambient, diffusion, oxidation, mechanics);
  deformation->setMaskInterface(mask, -1);
  deformation->apply();
  VC_TEST_ASSERT(deformation->lastSolveConverged())
  return deformation;
}

std::array<T, 9> inPlane(T s) {
  std::array<T, 9> stress{};
  stress[0] = s;  // xx
  stress[8] = -s; // zz
  return stress;
}

// sigma_xx range of the nitride interior and the largest speed anywhere
void measure(ls::SmartPointer<ls::OxidationDeformation<T, D>> deformation,
             T &minXX, T &maxXX, T &maxSpeed) {
  minXX = std::numeric_limits<T>::max();
  maxXX = -minXX;
  maxSpeed = 0.;
  deformation->forEachSolutionNode([&](const viennahrle::Index<D> &index, T,
                                       const std::array<T, 9> &s) {
    const auto v = deformation->getVelocity(index);
    maxSpeed = std::max(maxSpeed, std::hypot(v[0], v[1]));
    const T y = index[1] * gridDelta;
    if (y > padTop + 1.5 * gridDelta && y < nitrideTop - 1.5 * gridDelta) {
      minXX = std::min(minXX, s[0]);
      maxXX = std::max(maxXX, s[0]);
    }
  });
}

void testRelaxation() {
  const auto nitride = ls::OxidationPresets::siliconNitrideMask1000C();
  const T mu = nitride.youngModulus / (2. * (1. + nitride.poissonRatio));
  const T tau = nitride.viscosity / mu, dt = 0.5;
  const T exact = s0 * std::exp(-dt / tau);
  const auto uniform = [](const viennahrle::Index<D> &) { return inPlane(s0); };

  T minXX, maxXX, maxSpeed;
  measure(solveStack(uniform, nitride.viscosity, 0., dt), minXX, maxXX,
          maxSpeed);
  std::printf("Maxwell nitride, tau %.3f hr: sigma_xx %.6f..%.6f of exact "
              "s0 exp(-t/tau) (%.4f s0), max speed %.2e um/hr\n",
              tau, minXX / exact, maxXX / exact, exact / s0, maxSpeed);
  VC_TEST_ASSERT(std::abs(minXX / exact - 1.) < 1e-6 &&
                 std::abs(maxXX / exact - 1.) < 1e-6)
  VC_TEST_ASSERT(maxSpeed < 1e-9)

  T minE, maxE, speedE;
  measure(solveStack(uniform, nitride.viscosity, nitride.shearActivationVolume,
                     dt),
          minE, maxE, speedE);
  std::printf("Eyring nitride (V_N %.0f A^3): sigma_xx %.4f..%.4f s0\n",
              nitride.shearActivationVolume * 1e30, minE / s0, maxE / s0);
  VC_TEST_ASSERT(maxE < 0.9 * exact && minE > 0.)
}

void testElasticLimit() {
  // pre-stress in the left half only: the film moves
  const auto half = [](const viennahrle::Index<D> &index) {
    return inPlane(index[0] < 0 ? s0 : 0.);
  };
  const T dt = 0.01;
  auto elastic = solveStack(half, 0., 0., dt);
  auto stiff = solveStack(half, 1e30, 0., dt);
  T maxSpeed = 0., maxDifference = 0.;
  elastic->forEachSolutionNode(
      [&](const viennahrle::Index<D> &index, T, const std::array<T, 9> &) {
        const auto a = elastic->getVelocity(index);
        const auto b = stiff->getVelocity(index);
        maxSpeed = std::max(maxSpeed, std::hypot(a[0], a[1]));
        maxDifference =
            std::max(maxDifference, std::hypot(a[0] - b[0], a[1] - b[1]));
      });
  std::printf("elastic limit: max speed %.3e um/hr, eta_N 1e30 differs by "
              "%.2e of it\n",
              maxSpeed, maxDifference / maxSpeed);
  VC_TEST_ASSERT(maxSpeed > 0.)
  VC_TEST_ASSERT(maxDifference < 1e-9 * maxSpeed)
}

void testPreset() {
  // Senez Table II: poise -> Pa hr is / 10 / 3600
  const T points[4][3] = {{800., 6.0e16, 120.},
                          {900., 3.1e16, 130.},
                          {1000., 5.0e15, 170.},
                          {1100., 1.1e15, 350.}};
  for (const auto &p : points) {
    const auto m = ls::OxidationPresets::siliconNitrideMask(p[0]);
    VC_TEST_ASSERT(std::abs(m.viscosity / (p[1] / 36000.) - 1.) < 1e-12)
    VC_TEST_ASSERT(std::abs(m.shearActivationVolume / (p[2] * 1e-30) - 1.) <
                   1e-12)
  }
  const auto mid = ls::OxidationPresets::siliconNitrideMask(950.);
  VC_TEST_ASSERT(mid.viscosity < 3.1e16 / 36000. &&
                 mid.viscosity > 5.0e15 / 36000.)
}

int main() {
  omp_set_num_threads(2);
  viennacore::Logger::setLogLevel(viennacore::LogLevel::WARNING);
  testPreset();
  testRelaxation();
  testElasticLimit();
  return 0;
}
