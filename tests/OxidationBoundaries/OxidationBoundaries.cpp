// Boundary treatments of the oxidation solvers against exact solutions.
//
// 1. Free surface (OxidationDeformation): a film on a tilted Si surface
//    carrying a uniform in-plane pre-stress is in equilibrium, so the exact
//    velocity is zero and the exact pressure uniform.  A free surface must
//    keep the tangential stress (only sigma.n vanishes there).
// 2. Free surface through the grid points: a one-cell seed oxide around a
//    convex Si corner (nodes on the surface, the level set undefined beyond)
//    must give a regular system.
// 2b. Growth on a tilted Si front (Oxidation): a flat film on a planar front
//    at 5.7 degrees to the grid translates rigidly along the front normal,
//    so p = 0 and there is no shear.  The front's rate and speed must come
//    from its own interface points (the face crossings): per-node values
//    copied to the stair-stepped front, and Si-side ghosts with no oxide
//    neighbour, gave GPa stress and a film sliding down the slope.
// 2c. Stress-free birth (Oxidation): a flat dry film carrying a uniform
//    deviatoric pre-stress grows on a flat front.  The flow is a rigid
//    translation, so the old oxide keeps its pre-stress (decayed) and the
//    oxide formed during the run, born stress-free, has none.  Carrying the
//    history with weights renormalised over the oxide nodes handed the first
//    oxide layer its own old stress every step.
// 2d. Interior fill to a guide (Interior): a 150-cell film filled from its
//    surface with the Si as guide must be defined down to the Si.  The fill
//    stopped after 100 layers, so a thick film kept no stress history next
//    to the Si front.  The fill must keep the level-set width: it reset it
//    to 1, and meshing the filled oxide (ToSurfaceMesh expands a level set
//    narrower than 2) then dropped the fill and the oxide's stress history.
// 3. Reflective boundary (OxidationDiffusion): a strip whose centre lies on
//    the mirror plane must give the same oxidant concentration as the full
//    strip.
// 3b. Surface through grid points (OxidationDiffusion): dry oxidant through
//    the oxide around a Si disk whose oxide surface passes through grid
//    points must give the exact front speed of the cylinder.  The nodes on
//    the surface have rows of 1 / distance; the convergence test was
//    relative to their right-hand side (2e10), so the solve stopped at the
//    initial guess, C = 0 in the oxide, and the Si did not oxidise.
// 4. Nitride advection (OxidationMaskTracker): a nitride edge lifted rigidly
//    by 19 cells must keep its corners, and a nitride box turned rigidly
//    through 0.3 rad its shape; plain level-set advection rounds the
//    corners by numerical diffusion.  A lifted and bent strip on the mirror
//    plane must stay a valid level set over several process calls.
#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsCheck.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsGeometricAdvectDistributions.hpp>
#include <lsGeometries.hpp>
#include <lsInterior.hpp>
#include <lsMakeGeometry.hpp>
#include <lsOxidation.hpp>
#include <lsOxidationDeformation.hpp>
#include <lsOxidationDiffusion.hpp>
#include <lsOxidationPresets.hpp>
#include <lsTestAsserts.hpp>
#include <lsToSurfaceMesh.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>

#include <omp.h>

namespace ls = viennals;

using T = double;
constexpr int D = 2;
using LevelSet = ls::SmartPointer<ls::Domain<T, D>>;

ls::Domain<T, D>::BoundaryType boundaries[D] = {
    ls::Domain<T, D>::BoundaryType::REFLECTIVE_BOUNDARY,
    ls::Domain<T, D>::BoundaryType::INFINITE_BOUNDARY};

LevelSet makePlane(const double *bounds, T gridDelta,
                   ls::VectorType<T, D> origin, ls::VectorType<T, D> normal) {
  auto levelSet = LevelSet::New(bounds, boundaries, gridDelta);
  ls::MakeGeometry<T, D>(levelSet, ls::Plane<T, D>::New(origin, normal))
      .apply();
  return levelSet;
}

// Fill the level set's interior and store one stress tensor on every point.
void storeStress(LevelSet levelSet, const std::array<T, 9> &stress,
                 const std::string &name) {
  ls::Interior<T, D>(levelSet).apply();
  using VectorData = typename ls::PointData<T>::VectorDataType;
  VectorData rows[3];
  viennahrle::ConstSparseIterator<typename ls::Domain<T, D>::DomainType> it(
      levelSet->getDomain());
  for (; !it.isFinished(); ++it)
    if (it.isDefined())
      for (unsigned r = 0; r < 3; ++r)
        rows[r].push_back(
            {stress[3 * r], stress[3 * r + 1], stress[3 * r + 2]});
  for (unsigned r = 0; r < 3; ++r)
    levelSet->getPointData().insertReplaceVectorData(
        std::move(rows[r]), name + "R" + std::to_string(r));
}

void testTiltedFilm(T slope) {
  constexpr T gridDelta = 0.01, halfWidth = 1.5, thickness = 0.1, s0 = 1e9;
  const T angle = std::atan(slope);
  const ls::VectorType<T, D> normal{-std::sin(angle), std::cos(angle)};
  const ls::VectorType<T, D> tangent{std::cos(angle), std::sin(angle)};
  double bounds[2 * D] = {-halfWidth, halfWidth, -halfWidth * slope - 0.2,
                          halfWidth * slope + 0.3};
  auto si = makePlane(bounds, gridDelta, {0., 0.}, normal);
  auto ambient =
      makePlane(bounds, gridDelta,
                {normal[0] * thickness, normal[1] * thickness}, normal);
  // deviatoric pre-stress s0 (t t - n n)
  std::array<T, 9> stress{};
  for (unsigned i = 0; i < D; ++i)
    for (unsigned j = 0; j < D; ++j)
      stress[3 * i + j] =
          s0 * (tangent[i] * tangent[j] - normal[i] * normal[j]);
  storeStress(ambient, stress, "OxStress");

  ls::OxidationParameters oxidation; // no reaction: the Si front is fixed
  oxidation.reactionRate = 0.;
  oxidation.velocitySign = -1.;
  ls::OxidationDeformationParameters mechanics;
  mechanics.viscosity = 1e12;
  mechanics.shearModulus = 3e10;
  mechanics.stressTimeStep = 0.01;
  auto diffusion = ls::OxidationDiffusion<T, D>::New(si, ambient, oxidation);
  diffusion->apply();
  auto deformation = ls::OxidationDeformation<T, D>::New(si, ambient, diffusion,
                                                         oxidation, mechanics);
  deformation->apply();
  VC_TEST_ASSERT(deformation->lastSolveConverged())

  const T decay = std::exp(-mechanics.stressTimeStep * mechanics.shearModulus /
                           mechanics.viscosity);
  const T implicitViscosity = (1. - decay) * mechanics.viscosity;
  const T velocityScale = s0 * gridDelta / (2. * implicitViscosity);
  T maxVelocity = 0., maxPressureError = 0.;
  deformation->forEachSolutionNode([&](const viennahrle::Index<D> &index,
                                       T pressure, const std::array<T, 9> &) {
    if (std::abs(index[0] * gridDelta) > 0.3) // away from the mirrors
      return;
    const auto v = deformation->getVelocity(index);
    maxVelocity = std::max(maxVelocity, std::hypot(v[0], v[1]));
    maxPressureError =
        std::max(maxPressureError, std::abs(pressure + decay * s0));
  });
  std::printf("tilted film, slope %.2f: max|v| %.2e of the release scale, "
              "max pressure error %.2e of s0\n",
              slope, maxVelocity / velocityScale, maxPressureError / s0);
  std::fflush(stdout);
  VC_TEST_ASSERT(maxVelocity < 1e-2 * velocityScale)
  VC_TEST_ASSERT(maxPressureError < 1e-3 * s0)
}

void testCornerSeed() {
  constexpr T gridDelta = 0.01;
  double bounds[2 * D] = {0., 1., -1., 2.};
  auto si = makePlane(bounds, gridDelta, {0., 0.}, {0., 1.});
  auto fin = LevelSet::New(bounds, boundaries, gridDelta);
  ls::MakeGeometry<T, D>(fin,
                         ls::Box<T, D>::New(ls::VectorType<T, D>{-0.25, -0.1},
                                            ls::VectorType<T, D>{0.25, 1.}))
      .apply();
  ls::BooleanOperation<T, D>(si, fin, ls::BooleanOperationEnum::UNION).apply();
  auto ambient = LevelSet::New(si);
  ls::GeometricAdvect<T, D>(
      ambient,
      ls::SmartPointer<ls::SphereDistribution<viennahrle::CoordType, D>>::New(
          gridDelta))
      .apply();
  auto parameters = ls::OxidationPresets::wet1000CDealGrove100();
  parameters.velocitySign = -1.;
  auto mechanics = ls::OxidationPresets::oxideMechanics1000C(0.004);
  auto diffusion = ls::OxidationDiffusion<T, D>::New(si, ambient, parameters);
  diffusion->apply();
  auto deformation = ls::OxidationDeformation<T, D>::New(si, ambient, diffusion,
                                                         parameters, mechanics);
  deformation->apply();
  VC_TEST_ASSERT(deformation->lastSolveConverged())
}

// Dry 1000 C oxidation (Drob's Deal-Grove fit), isotropic, no stress
// kinetics, linear Maxwell oxide, run to tEnd with CFL-limited steps.
ls::SmartPointer<ls::Oxidation<T, D>> growDryFilm(LevelSet si, LevelSet ambient,
                                                  T gridDelta, T tEnd) {
  ls::OxidationParameters ox;
  ox.diffusionCoefficient = 0.064 / 2;
  ox.reactionRate = 0.0564;
  ox.transferCoefficient = 100.;
  ox.equilibriumConcentration = 1.;
  ox.oxidantMoleculeDensity = 1.;
  ox.expansionCoefficient = 2.27;
  ox.velocitySign = -1.;
  ox.temperature = 1273.15;
  ox.reactionRateRatio111 = 1.;
  auto model = ls::Oxidation<T, D>::New(si, ambient);
  model->setOxidationParameters(ox);
  ls::OxidationCouplingParameters coupling;
  coupling.maxIterations = 100;
  coupling.tolerance = 2e-2;
  model->setCouplingParameters(coupling);
  const auto preset = ls::OxidationPresets::oxideMechanics1000C(1.0);
  T t = 0., dt = 0.004;
  while (tEnd - t > 1e-9) {
    const T requested = std::min(dt, tEnd - t);
    ls::OxidationDeformationParameters p;
    p.viscosity = 6.5e11;
    p.shearModulus = preset.shearModulus;
    p.stressTimeStep = requested;
    model->setDeformationParameters(p);
    const T actual = model->applyCFLLimited(requested, 0.499);
    t += actual;
    const T vmax = model->getLastMaxVelocity();
    dt = std::min(vmax > 0 ? 0.499 * gridDelta / vmax : requested, 2 * actual);
  }
  return model;
}

// sqrt(J2) of the in-plane deviatoric part of a 2-D Cauchy stress
T deviatoricMagnitude(const std::array<T, 9> &s) {
  const T mean = 0.5 * (s[0] + s[4]);
  const T dxx = s[0] - mean, dyy = s[4] - mean;
  return std::sqrt(0.5 * (dxx * dxx + dyy * dyy) + s[1] * s[1]);
}

void testTiltedFrontGrowth() {
  constexpr T gridDelta = 0.0024, halfWidth = 0.25, thickness = 0.02,
              slope = 0.1, window = 0.1;
  const T angle = std::atan(slope);
  const ls::VectorType<T, D> normal{-std::sin(angle), std::cos(angle)};
  double bounds[2 * D] = {-halfWidth, halfWidth, -halfWidth * slope - 0.1,
                          halfWidth * slope + 0.15};
  auto si = makePlane(bounds, gridDelta, {0., 0.}, normal);
  auto ambient =
      makePlane(bounds, gridDelta,
                {normal[0] * thickness, normal[1] * thickness}, normal);
  auto model = growDryFilm(si, ambient, gridDelta, 0.05);
  auto field = model->getDeformationField();
  T maxPressure = 0., maxShear = 0., meanNormal = 0., meanTangential = 0.;
  int count = 0;
  field->forEachSolutionNode([&](const viennahrle::Index<D> &index, T pressure,
                                 const std::array<T, 9> &s) {
    if (std::abs(index[0] * gridDelta) > window)
      return;
    maxPressure = std::max(maxPressure, std::abs(pressure));
    maxShear = std::max(maxShear, deviatoricMagnitude(s));
    const auto v = field->getVelocity(index);
    meanNormal += v[0] * normal[0] + v[1] * normal[1];
    meanTangential += v[0] * normal[1] - v[1] * normal[0];
    ++count;
  });
  meanNormal /= count;
  meanTangential /= count;
  std::printf("tilted front growth: max|p| %.3f GPa, max shear %.3f GPa, "
              "mean tangential / normal velocity %.4f\n",
              maxPressure / 1e9, maxShear / 1e9, meanTangential / meanNormal);
  std::fflush(stdout);
  // measured: 0.98 / 2.47 GPa and -0.047 with per-node rates and front
  // velocities, 0.002 / 0.002 GPa and 0.0002 with the interface points
  VC_TEST_ASSERT(maxPressure < 0.05e9)
  VC_TEST_ASSERT(maxShear < 0.05e9)
  VC_TEST_ASSERT(std::abs(meanTangential) < 5e-3 * meanNormal)
}

void testStressFreeBirth() {
  constexpr T gridDelta = 0.001, halfWidth = 0.05, thickness = 0.02, s0 = 1e9,
              tEnd = 0.2;
  double bounds[2 * D] = {-halfWidth, halfWidth, -0.02, thickness + 0.02};
  auto si = makePlane(bounds, gridDelta, {0., 0.}, {0., 1.});
  auto ambient = makePlane(bounds, gridDelta, {0., thickness}, {0., 1.});
  std::array<T, 9> stress{};
  stress[0] = s0;
  stress[4] = -s0;
  storeStress(ambient, stress, "OxStress");
  auto model = growDryFilm(si, ambient, gridDelta, tEnd);
  const auto preset = ls::OxidationPresets::oxideMechanics1000C(1.0);
  const T decay = std::exp(-tEnd * preset.shearModulus / 6.5e11);
  // The front receded to y_f < 0; the oxide made from the consumed Si
  // reaches up to (alpha - 1) |y_f|, the old oxide lies above it.  The
  // semi-Lagrangian transport smears the jump between them over ~3 cells.
  auto mesh = ls::SmartPointer<ls::Mesh<T>>::New();
  ls::ToSurfaceMesh<T, D>(si, mesh).apply();
  T front = 0.;
  for (const auto &node : mesh->getNodes())
    front += node[1];
  front /= mesh->getNodes().size();
  const T interface = (2.27 - 1.) * -front;
  T maxNew = 0., minOld = std::numeric_limits<T>::max(), maxOld = 0.;
  int countNew = 0;
  model->getDeformationField()->forEachSolutionNode(
      [&](const viennahrle::Index<D> &index, T, const std::array<T, 9> &s) {
        const T y = index[1] * gridDelta, magnitude = deviatoricMagnitude(s);
        if (y < interface - 4 * gridDelta) {
          maxNew = std::max(maxNew, magnitude);
          ++countNew;
        } else if (y > interface + 4 * gridDelta &&
                   y < interface + thickness - 4 * gridDelta) {
          minOld = std::min(minOld, magnitude);
          maxOld = std::max(maxOld, magnitude);
        }
      });
  std::printf("stress-free birth: front at %.2f nm, %d new oxide nodes, max "
              "deviatoric stress %.3f of s0 there; old oxide %.4f..%.4f of "
              "decayed s0\n",
              front * 1e3, countNew, maxNew / s0, minOld / (decay * s0),
              maxOld / (decay * s0));
  std::fflush(stdout);
  // measured: 0.991 of s0 in the new oxide with the renormalised transport,
  // 0.007 born stress-free; old oxide 0.9988..1.0000
  VC_TEST_ASSERT(countNew > 0)
  VC_TEST_ASSERT(maxNew < 0.02 * s0)
  VC_TEST_ASSERT(minOld > 0.98 * decay * s0 && maxOld < 1.02 * decay * s0)
}

void testGuidedFill() {
  constexpr T gridDelta = 0.001, thickness = 0.15;
  double bounds[2 * D] = {-0.01, 0.01, -0.02, thickness + 0.02};
  auto si = makePlane(bounds, gridDelta, {0., 0.}, {0., 1.});
  auto ambient = makePlane(bounds, gridDelta, {0., thickness}, {0., 1.});
  ls::Expand<T, D>(ambient, 2).apply(); // the oxide after Advect
  ls::Interior<T, D> fill(ambient);
  fill.setGuide(si);
  fill.apply();
  const int top = static_cast<int>(std::lround(thickness / gridDelta));
  auto undefined = [&]() {
    viennahrle::ConstSparseIterator<typename ls::Domain<T, D>::DomainType> it(
        ambient->getDomain());
    int missing = 0;
    for (int j = 1; j < top; ++j) {
      it.goToIndices(viennahrle::Index<D>(0, j));
      if (!it.isDefined())
        ++missing;
    }
    return missing;
  };
  const int missing = undefined();
  std::printf("guided fill: %d of %d film points undefined, width %d\n",
              missing, top - 1, ambient->getLevelSetWidth());
  // measured: 48 of 149 undefined with the 100-layer stop, 0 filled to the
  // guide
  VC_TEST_ASSERT(missing == 0)
  VC_TEST_ASSERT(ambient->getLevelSetWidth() == 2)
  // meshing the filled level set must leave it filled
  auto mesh = ls::SmartPointer<ls::Mesh<T>>::New();
  ls::ToSurfaceMesh<T, D>(ambient, mesh).apply();
  const int missingAfterMesh = undefined();
  std::printf("guided fill: %d undefined after ToSurfaceMesh\n",
              missingAfterMesh);
  std::fflush(stdout);
  // measured: 148 of 149 before the fix
  VC_TEST_ASSERT(missingAfterMesh == 0)
}

ls::SmartPointer<ls::OxidationDiffusion<T, D>> solveStrip(bool mirrored) {
  // pad oxide 30 nm, nitride 100 nm, strip half width 0.305 um (edges
  // between grid points); the mirrored run has the strip centre on its
  // lower reflective boundary
  constexpr T gridDelta = 0.01, pad = 0.03, halfStrip = 0.305;
  const T extent = mirrored ? 0.6 : 1.2;
  const T maskLow = mirrored ? -2. : -halfStrip;
  const T maskHigh = mirrored ? -0.6 + halfStrip : halfStrip;
  double bounds[2 * D] = {-extent, extent, -0.5, 0.5};
  auto si = makePlane(bounds, gridDelta, {0., 0.}, {0., 1.});
  auto ambient = makePlane(bounds, gridDelta, {0., pad}, {0., 1.});
  auto mask = LevelSet::New(bounds, boundaries, gridDelta);
  ls::MakeGeometry<T, D> box(
      mask, ls::Box<T, D>::New(ls::VectorType<T, D>{maskLow, pad - 1e-6},
                               ls::VectorType<T, D>{maskHigh, pad + 0.1}));
  box.setIgnoreBoundaryConditions(std::array<bool, 3>{false, true, false});
  box.apply();
  auto parameters = ls::OxidationPresets::wet1000CDealGrove100();
  parameters.velocitySign = -1.;
  parameters.maskTransferCoefficient = 0.;
  parameters.maskConcentration = 0.;
  parameters.tolerance = 1e-12;
  parameters.maxIterations = 100000;
  auto diffusion = ls::OxidationDiffusion<T, D>::New(si, ambient, parameters);
  diffusion->setMaskInterface(mask, -1);
  diffusion->apply();
  return diffusion;
}

void testDiffusionMirror() {
  auto mirrored = solveStrip(true);
  auto full = solveStrip(false);
  for (int offset : {0, 1, 2, 5, 10, 20, 30, 31, 40}) {
    const T y = 0.015;
    const T cm =
        mirrored->getConcentration(ls::Vec3D<T>{-0.6 + offset * 0.01, y, 0.});
    const T cf = full->getConcentration(ls::Vec3D<T>{offset * 0.01, y, 0.});
    VC_TEST_ASSERT(cf > 0.)
    VC_TEST_ASSERT(std::abs(cm - cf) < 1e-4 * cf)
  }
}

void testSurfaceThroughGridPoints() {
  // dry 1000 C kinetics (Deal-Grove, Si(110)), no stress laws
  constexpr T gridDelta = 0.01, a = 0.36, b = 0.39;
  double bounds[2 * D] = {-0.6, 0.6, -0.6, 0.6};
  auto si = LevelSet::New(bounds, boundaries, gridDelta);
  ls::MakeGeometry<T, D>(si,
                         ls::Sphere<T, D>::New(ls::VectorType<T, D>{0., 0.}, a))
      .apply();
  auto ambient = LevelSet::New(bounds, boundaries, gridDelta);
  ls::MakeGeometry<T, D>(ambient,
                         ls::Sphere<T, D>::New(ls::VectorType<T, D>{0., 0.}, b))
      .apply();
  ls::OxidationParameters parameters;
  parameters.diffusionCoefficient = 0.00521938;
  parameters.reactionRate = 0.0651269;
  parameters.transferCoefficient = 100.;
  parameters.equilibriumConcentration = 1.;
  parameters.oxidantMoleculeDensity = 1.;
  parameters.expansionCoefficient = 2.27;
  parameters.velocitySign = -1.;
  parameters.reactionRateRatio111 = 1.;
  parameters.tolerance = 1e-7; // as ViennaPS
  auto diffusion = ls::OxidationDiffusion<T, D>::New(si, ambient, parameters);
  diffusion->apply();
  // flux per radian through the annulus, the front concentration F / (a k)
  const T flux = 1. / (1. / (b * parameters.transferCoefficient) +
                       std::log(b / a) / parameters.diffusionCoefficient +
                       1. / (a * parameters.reactionRate));
  const T exact = flux / a / parameters.expansionCoefficient;
  for (T angle : {0., 0.5, 1., 1.5707963267948966}) {
    const auto front = diffusion->frontVelocityAt(
        ls::Vec3D<T>{a * std::cos(angle), a * std::sin(angle), 0.});
    std::printf("surface through grid points: front speed at %.2f rad %.5f "
                "um/hr, exact %.5f\n",
                angle, front.speed, exact);
    VC_TEST_ASSERT(front.found)
    // measured: 0.2 % (0 before the fix)
    VC_TEST_ASSERT(std::abs(front.speed / exact - 1.) < 0.01)
  }
}

class Lift final : public ls::VelocityField<T> {
public:
  static constexpr T speed = 0.193; // um/hr
  ls::Vec3D<T> getVectorVelocity(const ls::Vec3D<T> &, int,
                                 const ls::Vec3D<T> &, unsigned long) final {
    return {0., speed, 0.};
  }
};

// Outermost x of the surface at the given height above its lowest point.
T wallAbove(LevelSet levelSet, T height) {
  auto mesh = ls::SmartPointer<ls::Mesh<T>>::New();
  ls::ToSurfaceMesh<T, D>(levelSet, mesh).apply();
  const auto &nodes = mesh->getNodes();
  T yMin = std::numeric_limits<T>::max();
  for (const auto &node : nodes)
    yMin = std::min(yMin, node[1]);
  const T y = yMin + height;
  T x = std::numeric_limits<T>::lowest();
  for (const auto &line : mesh->lines) {
    const auto &a = nodes[line[0]];
    const auto &b = nodes[line[1]];
    if (a[1] != b[1] && y >= std::min(a[1], b[1]) && y <= std::max(a[1], b[1]))
      x = std::max(x, a[0] + (b[0] - a[0]) * (y - a[1]) / (b[1] - a[1]));
  }
  return x;
}

void testMaskTracker() {
  // the edge of a 140 nm nitride on a 5 nm grid, lifted 0.1 um in 84 steps
  constexpr T gridDelta = 0.005, edge = -0.7375, dt = 0.006;
  double bounds[2 * D] = {-0.7625, 0.7625, -1., 1.5};
  for (bool tracked : {false, true}) {
    auto mask = LevelSet::New(bounds, boundaries, gridDelta);
    ls::MakeGeometry<T, D>(
        mask, ls::Box<T, D>::New(ls::VectorType<T, D>{-0.8625, 0.015},
                                 ls::VectorType<T, D>{edge, 0.155}))
        .apply();
    ls::OxidationMaskTracker<T, D> tracker(mask);
    auto lift = ls::SmartPointer<Lift>::New();
    for (int step = 0; step < 84; ++step) {
      ls::Advect<T, D> advect;
      advect.insertNextLevelSet(mask);
      advect.setVelocityField(lift);
      advect.setAdvectionTime(dt);
      advect.apply();
      if (tracked)
        tracker.correct(
            mask,
            [](const ls::Vec3D<T> &) {
              return ls::Vec3D<T>{0., Lift::speed, 0.};
            },
            dt);
    }
    const T inward = edge - wallAbove(mask, 0.0025);
    std::printf("nitride edge %s: side wall 2.5 nm above the bottom moved in "
                "by %.1f nm\n",
                tracked ? "tracked" : "advected", inward * 1e3);
    if (tracked)
      VC_TEST_ASSERT(inward < 0.002)
  }
}

class Spin final : public ls::VelocityField<T> {
public:
  static constexpr T omega = 0.6; // rad/hr about the origin
  ls::Vec3D<T> getVectorVelocity(const ls::Vec3D<T> &x, int,
                                 const ls::Vec3D<T> &, unsigned long) final {
    return {-omega * x[1], omega * x[0], 0.};
  }
};

// Largest distance of the surface from the box |x| <= a, |y| <= b turned by
// angle about the origin.
T turnedBoxDeviation(LevelSet levelSet, T a, T b, T angle) {
  auto mesh = ls::SmartPointer<ls::Mesh<T>>::New();
  ls::ToSurfaceMesh<T, D>(levelSet, mesh).apply();
  const T c = std::cos(angle), s = std::sin(angle);
  T worst = 0.;
  for (const auto &node : mesh->getNodes()) {
    const T x = c * node[0] + s * node[1], y = -s * node[0] + c * node[1];
    const T qx = std::abs(x) - a, qy = std::abs(y) - b;
    const T outside = std::hypot(std::max(qx, T(0)), std::max(qy, T(0)));
    const T inside = std::min(std::max(qx, qy), T(0));
    worst = std::max(worst, std::abs(outside + inside));
  }
  return worst;
}

void testMaskRotation() {
  constexpr T gridDelta = 0.005, a = 0.1, b = 0.07, dt = 0.006;
  constexpr int steps = 84; // 0.3 rad
  double bounds[2 * D] = {-0.5, 0.5, -0.5, 0.5};
  for (bool tracked : {false, true}) {
    auto mask = LevelSet::New(bounds, boundaries, gridDelta);
    ls::MakeGeometry<T, D>(mask,
                           ls::Box<T, D>::New(ls::VectorType<T, D>{-a, -b},
                                              ls::VectorType<T, D>{a, b}))
        .apply();
    ls::OxidationMaskTracker<T, D> tracker(mask);
    auto spin = ls::SmartPointer<Spin>::New();
    for (int step = 0; step < steps; ++step) {
      ls::Advect<T, D> advect;
      advect.insertNextLevelSet(mask);
      advect.setVelocityField(spin);
      advect.setAdvectionTime(dt);
      advect.apply();
      if (tracked)
        tracker.correct(
            mask,
            [](const ls::Vec3D<T> &x) {
              return ls::Vec3D<T>{-Spin::omega * x[1], Spin::omega * x[0], 0.};
            },
            dt);
    }
    const T deviation =
        turnedBoxDeviation(mask, a, b, Spin::omega * steps * dt);
    std::printf("nitride box %s through 0.3 rad: largest distance from the "
                "turned box %.1f nm\n",
                tracked ? "tracked" : "advected", deviation * 1e3);
    if (tracked)
      VC_TEST_ASSERT(deviation < 0.002)
  }
}

// Lift growing toward the strip edge plus an outward drift: the nitride
// bends and stretches.  Before the tracker normalised its band this left
// neighbours like -0.504 / +0.527 across the surface (invalid for the
// level-set algorithms; LOCOS runs crashed).
ls::Vec3D<T> bendVelocity(const ls::Vec3D<T> &x) {
  const T s = (x[0] + 0.9125) / 0.325;
  return {0.02 * s, 0.05 + 0.25 * s * s, 0.};
}

class Bend final : public ls::VelocityField<T> {
public:
  ls::Vec3D<T> getVectorVelocity(const ls::Vec3D<T> &x, int,
                                 const ls::Vec3D<T> &, unsigned long) final {
    return bendVelocity(x);
  }
};

void testMaskTrackerBending() {
  constexpr T gridDelta = 0.005, dt = 0.006;
  double bounds[2 * D] = {-0.9125, 0.9125, -1., 1.5};
  auto mask = LevelSet::New(bounds, boundaries, gridDelta);
  ls::MakeGeometry<T, D>(
      mask, ls::Box<T, D>::New(ls::VectorType<T, D>{-1.0125, 0.015},
                               ls::VectorType<T, D>{-0.5875, 0.155}))
      .apply();
  auto bend = ls::SmartPointer<Bend>::New();
  for (int call = 0; call < 3; ++call) {
    ls::OxidationMaskTracker<T, D> tracker(mask);
    for (int step = 0; step < 42; ++step) {
      ls::Advect<T, D> advect;
      advect.insertNextLevelSet(mask);
      advect.setVelocityField(bend);
      advect.setAdvectionTime(dt);
      advect.apply();
      tracker.correct(mask, bendVelocity, dt);
      ls::Check<T, D> check(mask);
      check.apply();
      VC_TEST_ASSERT(check.isValid())
      ls::Interior<T, D>(mask).apply();
    }
  }
}

int main() {
  omp_set_num_threads(2);
  viennacore::Logger::setLogLevel(viennacore::LogLevel::WARNING);
  for (T slope : {0., 0.5, 1.})
    testTiltedFilm(slope);
  testCornerSeed();
  testTiltedFrontGrowth();
  testStressFreeBirth();
  testGuidedFill();
  testDiffusionMirror();
  testSurfaceThroughGridPoints();
  testMaskTracker();
  testMaskRotation();
  testMaskTrackerBending();
  return 0;
}
