#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <string>

#include <lsDomain.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsOxidation.hpp>
#include <lsOxidationPresets.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

using NumericType = double;
constexpr int D = 2;

using LevelSet = ls::SmartPointer<ls::Domain<NumericType, D>>;

void writeSurface(LevelSet levelSet, const std::string &fileName) {
  if (levelSet->getLevelSetWidth() < 2)
    ls::Expand<NumericType, D>(levelSet, 2).apply();
  auto mesh = ls::Mesh<NumericType>::New();
  auto surfaceMesh = ls::ToSurfaceMesh<NumericType, D>(levelSet, mesh);
  surfaceMesh.setSharpCorners(false);
  surfaceMesh.apply();
  ls::VTKWriter<NumericType>(mesh, fileName).apply();
}

LevelSet makePlane(const double *bounds,
                   ls::Domain<NumericType, D>::BoundaryType *boundaryCons,
                   NumericType gridDelta, NumericType y) {
  auto levelSet =
      ls::Domain<NumericType, D>::New(bounds, boundaryCons, gridDelta);
  const ls::VectorType<NumericType, D> origin{0., y};
  const ls::VectorType<NumericType, D> normal{0., 1.};
  ls::MakeGeometry<NumericType, D>(
      levelSet, ls::Plane<NumericType, D>::New(origin, normal))
      .apply();
  return levelSet;
}

LevelSet makeMask(const double *bounds,
                  ls::Domain<NumericType, D>::BoundaryType *boundaryCons,
                  NumericType gridDelta, NumericType xMin, NumericType xMax,
                  NumericType yMin, NumericType yMax) {
  auto mask = ls::Domain<NumericType, D>::New(bounds, boundaryCons, gridDelta);
  const ls::VectorType<NumericType, D> minCorner{xMin, yMin};
  const ls::VectorType<NumericType, D> maxCorner{xMax, yMax};
  ls::MakeGeometry<NumericType, D> makeGeometry(
      mask, ls::Box<NumericType, D>::New(minCorner, maxCorner));
  makeGeometry.setIgnoreBoundaryConditions(std::array<bool, D>{false, true});
  makeGeometry.apply();
  return mask;
}

void writeDiagnostics(
    const ls::SmartPointer<ls::OxidationDiffusion<NumericType, D>> &diffusion,
    const ls::SmartPointer<ls::OxidationDeformation<NumericType, D>>
        &deformation,
    const viennahrle::Index<D> &minIndex, const viennahrle::Index<D> &maxIndex,
    NumericType gridDelta, const std::string &fileName) {
  std::ofstream file(fileName);
  file << "x,y,concentration,velocity_x,velocity_y,pressure,strain_trace,"
          "von_mises_stress\n";
  for (auto j = minIndex[1]; j <= maxIndex[1]; ++j) {
    for (auto i = minIndex[0]; i <= maxIndex[0]; ++i) {
      ls::Vec3D<NumericType> coordinate{i * gridDelta, j * gridDelta, 0.};
      const auto velocity = deformation->getVelocity(coordinate);
      file << coordinate[0] << ',' << coordinate[1] << ','
           << diffusion->getConcentration(coordinate) << ',' << velocity[0]
           << ',' << velocity[1] << ',' << deformation->getPressure(coordinate)
           << ',' << deformation->getStrainTrace(coordinate) << ','
           << deformation->getVonMisesStress(coordinate) << '\n';
    }
  }
}

int main(int argc, char *argv[]) {
  omp_set_num_threads(4);

  // Pass "gpu" as first argument to use the GPU BiCGSTAB solver, e.g.:
  //   ./LOCOSOxidation gpu
  // Requires a GPU build (-DVIENNALS_USE_GPU=ON). CPU is the default.
  const bool useGpu = argc > 1 && std::string(argv[1]) == "gpu";

  // Geometry — matches ViennaPS locosOxidation defaults.
  constexpr NumericType gridDelta = 0.01;           // um
  constexpr NumericType xExtent = 1.;               // um  (REFLECTIVE at x=0)
  constexpr NumericType yMin = -1.;                 // um
  constexpr NumericType yMax = 2.;                  // um
  constexpr NumericType padOxideThickness = 0.03;   // um
  constexpr NumericType maskThickness = 0.05;       // um
  constexpr NumericType maskEdge = 0.;              // open window: x > 0
  constexpr NumericType advectionTime = 0.1;        // hr
  constexpr NumericType timeStep = 0.01;            // hr per outer step
  constexpr NumericType maskContactEpsilon = 1.e-6; // um

  double bounds[2 * D] = {-xExtent, xExtent, yMin, yMax};
  ls::Domain<NumericType, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] =
      ls::Domain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[1] = ls::Domain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  // Si substrate at y = 0.
  auto siInterface = makePlane(bounds, boundaryCons, gridDelta, 0.);

  // Pad SiO2: geometric offset of Si by padOxideThickness.
  auto ambientInterface = ls::Domain<NumericType, D>::New(siInterface);
  auto initialOxide =
      ls::SmartPointer<ls::SphereDistribution<viennahrle::CoordType, D>>::New(
          padOxideThickness);
  ls::GeometricAdvect<NumericType, D>(ambientInterface, initialOxide).apply();

  // Capture pre-oxidation state for volume conservation check.
  auto siInitialForConservation = ls::Domain<NumericType, D>::New(siInterface);
  auto ambientInitialForConservation =
      ls::Domain<NumericType, D>::New(ambientInterface);

  // Si3N4 mask box: covers x ∈ [−xExtent, maskEdge], sitting on the pad oxide.
  // The tiny contact epsilon ensures Cartesian stencils unambiguously cross the
  // mask boundary at the contact face.
  auto maskInterface =
      makeMask(bounds, boundaryCons, gridDelta, -xExtent, maskEdge,
               padOxideThickness - maskContactEpsilon,
               padOxideThickness + maskThickness);

  writeSurface(siInterface, "locos_si_initial.vtp");
  writeSurface(ambientInterface, "locos_ambient_initial.vtp");
  writeSurface(maskInterface, "locos_mask.vtp");

  // ── Oxidation parameters ──────────────────────────────────────────────────

  auto oxParams = ls::OxidationPresets::wet1000CDealGrove100();
  oxParams.velocitySign = -1.;
  oxParams.maskTransferCoefficient = 0.; // nitride is a perfect oxidant block
  oxParams.maskConcentration = 0.;
  oxParams.maxIterations = 10000;
  oxParams.tolerance = 1.e-7;

  auto defParams = ls::OxidationPresets::oxideMechanics1000C(timeStep);

  // Coupling parameters — match ViennaPS locosOxidation config.
  ls::OxidationCouplingParameters couplingParams;
  couplingParams.maxIterations = 100;
  couplingParams.tolerance = 2.e-2;
  couplingParams.relaxation = 1.;

  auto maskParams = ls::OxidationPresets::siliconNitrideMask1000C();

  // ── Solve bounds ──────────────────────────────────────────────────────────

  const auto toIndex = [&](NumericType x) {
    return static_cast<viennahrle::IndexType>(std::llround(x / gridDelta));
  };

  // Diffusion/deformation: covers the entire oxide-band region.
  viennahrle::Index<D> diffMinIndex{toIndex(-xExtent), toIndex(yMin)};
  viennahrle::Index<D> diffMaxIndex{toIndex(xExtent), toIndex(yMax)};

  // Mask bending: brackets the mask geometry with a one-cell margin.
  viennahrle::Index<D> maskMinIndex{toIndex(-xExtent),
                                    toIndex(padOxideThickness) - 1};
  viennahrle::Index<D> maskMaxIndex{
      toIndex(maskEdge), toIndex(padOxideThickness + maskThickness) + 1};

  // ── LOCOS solver ──────────────────────────────────────────────────────────

  auto locos = ls::Oxidation<NumericType, D>::New(siInterface, ambientInterface,
                                                  maskInterface);
  locos->setOxidationParameters(oxParams);
  locos->setDeformationParameters(defParams);
  locos->setCouplingParameters(couplingParams);
  locos->setMaskParameters(maskParams);
  locos->setSolveBounds(diffMinIndex, diffMaxIndex);
  locos->setMaskBendingBounds(maskMinIndex, maskMaxIndex);
  locos->setMaskCouplingIterations(30);
  locos->setMaskCouplingTolerance(NumericType(1.e-2));
  if (useGpu)
    locos->setGpuMode(ls::GpuMode::Gpu);

  // ── Time-stepping loop ────────────────────────────────────────────────────

  NumericType elapsed = 0.;
  const NumericType timeEps = advectionTime * NumericType(1e-8);
  while (advectionTime - elapsed > timeEps) {
    const NumericType dt = std::min(timeStep, advectionTime - elapsed);
    elapsed += locos->applyCFLLimited(dt, NumericType(0.499));
    std::cout << "t = " << elapsed << " hr\n";
  }

  // ── Diagnostics ───────────────────────────────────────────────────────────

  const auto diffusion = locos->getDiffusionField();
  const auto deformation = locos->getDeformationField();
  const auto maskBending = locos->getMaskBendingField();

  // Sample points in open window and under mask.
  const ls::Vec3D<NumericType> openReactionPoint{0.5, gridDelta, 0.};
  const ls::Vec3D<NumericType> maskedReactionPoint{-0.5, gridDelta, 0.};
  const NumericType openConcentration =
      diffusion->getConcentration(openReactionPoint);
  const NumericType maskedConcentration =
      diffusion->getConcentration(maskedReactionPoint);
  const NumericType openSiliconSpeed = std::abs(
      diffusion->getScalarVelocity(openReactionPoint, 0, {0., 1., 0.}, 0));
  const NumericType maskedSiliconSpeed = std::abs(
      diffusion->getScalarVelocity(maskedReactionPoint, 0, {0., 1., 0.}, 0));

  std::cout << "Diffusion nodes: " << diffusion->getNumberOfSolutionNodes()
            << ", iterations: " << diffusion->getIterations()
            << ", residual: " << diffusion->getResidual() << '\n';
  std::cout << "Deformation nodes: " << deformation->getNumberOfSolutionNodes()
            << ", iterations: " << deformation->getIterations()
            << ", residual: " << deformation->getResidual() << '\n';
  std::cout << "Average oxide expansion speed: "
            << deformation->avgExpansionSpeed() << " um/hr\n";
  std::cout << "Open-window concentration: " << openConcentration
            << ", masked concentration: " << maskedConcentration << '\n';
  std::cout << "Open-window Si speed: " << openSiliconSpeed
            << " um/hr, masked Si speed: " << maskedSiliconSpeed
            << " um/hr, suppression ratio: "
            << (openSiliconSpeed > 0. ? maskedSiliconSpeed / openSiliconSpeed
                                      : 0.)
            << '\n';

  if (!std::isfinite(openSiliconSpeed) || !std::isfinite(maskedSiliconSpeed) ||
      openSiliconSpeed <= 0. || maskedSiliconSpeed / openSiliconSpeed > 0.05) {
    std::cerr << "LOCOS mask sanity check failed: masked oxidation is not "
                 "sufficiently suppressed.\n";
    return 1;
  }

  writeDiagnostics(diffusion, deformation, diffMinIndex, diffMaxIndex,
                   gridDelta, "locos_oxidation_diagnostics.csv");

  const ls::Vec3D<NumericType> maskBottomSample{-0.5, padOxideThickness, 0.};
  const ls::Vec3D<NumericType> maskAnchorSample{
      -xExtent + gridDelta, padOxideThickness + maskThickness * 0.5, 0.};
  const ls::Vec3D<NumericType> maskMidSample{
      -0.5, padOxideThickness + maskThickness * 0.5, 0.};
  const ls::Vec3D<NumericType> maskEdgeSample{
      -gridDelta, padOxideThickness + maskThickness * 0.5, 0.};
  const ls::Vec3D<NumericType> maskContactNodeSample{-gridDelta,
                                                     padOxideThickness, 0.};
  const ls::Vec3D<NumericType> maskTopSample{
      -0.5, padOxideThickness + maskThickness, 0.};
  const auto maskBottomVelocity =
      maskBending->getVectorVelocity(maskBottomSample, 0, {0., -1., 0.}, 0);
  const auto maskTopVelocity =
      maskBending->getVectorVelocity(maskTopSample, 0, {0., 1., 0.}, 0);
  const auto maskContactNodeVelocity = maskBending->getVectorVelocity(
      maskContactNodeSample, 0, {0., -1., 0.}, 0);
  const auto maskAnchorVelocity =
      maskBending->getVectorVelocity(maskAnchorSample, 0, {0., 1., 0.}, 0);
  const auto maskMidVelocity =
      maskBending->getVectorVelocity(maskMidSample, 0, {0., 1., 0.}, 0);
  const auto maskEdgeVelocity =
      maskBending->getVectorVelocity(maskEdgeSample, 0, {0., 1., 0.}, 0);
  const NumericType maskBottomPressure =
      deformation->getPressure(maskBottomSample);
  std::cout << "Mask elasticity nodes: "
            << maskBending->getNumberOfSolutionNodes()
            << ", contact nodes: " << maskBending->getNumberOfContactNodes()
            << ", iterations: " << maskBending->getIterations()
            << ", residual: " << maskBending->getResidual() << '\n';
  std::cout << "Oxide/mask interface solve iterations: "
            << locos->getMaskCouplingIterations()
            << ", residual: " << locos->getMaskCouplingResidual() << '\n';
  std::cout << "Mask bottom velocity: (" << maskBottomVelocity[0] << ", "
            << maskBottomVelocity[1] << ") um/hr, top velocity: ("
            << maskTopVelocity[0] << ", " << maskTopVelocity[1]
            << ") um/hr, oxide pressure: " << maskBottomPressure << " Pa\n";
  std::cout << "Mask contact-node velocity: (" << maskContactNodeVelocity[0]
            << ", " << maskContactNodeVelocity[1] << ") um/hr\n";
  std::cout << "Mask lateral bending samples: anchor (" << maskAnchorVelocity[0]
            << ", " << maskAnchorVelocity[1] << "), mid (" << maskMidVelocity[0]
            << ", " << maskMidVelocity[1] << "), edge (" << maskEdgeVelocity[0]
            << ", " << maskEdgeVelocity[1] << ") um/hr\n";

  // Volume conservation in the open window.
  const auto conservation = ls::computeLOCOSOpenWindowConservation<NumericType>(
      siInitialForConservation, siInterface, ambientInitialForConservation,
      ambientInterface, 0.1, 0.9, toIndex(yMin), toIndex(yMax),
      oxParams.expansionCoefficient);
  std::cout << "Open-window conservation samples: " << conservation.samples
            << ", Si consumed area: " << conservation.siliconRecession
            << " um^2, ambient lift area: " << conservation.ambientLift
            << " um^2, expected ambient lift: "
            << conservation.expectedAmbientLift
            << " um^2, lift/Si ratio: " << conservation.ambientLiftRatio
            << " (expected " << (oxParams.expansionCoefficient - 1.)
            << "), relative error: " << conservation.relativeError << '\n';

  writeSurface(siInterface, "locos_si_after.vtp");
  writeSurface(ambientInterface, "locos_ambient_after.vtp");
  writeSurface(maskInterface, "locos_mask_after.vtp");

  std::cout << "Wrote locos_si_initial.vtp, locos_ambient_initial.vtp, "
               "locos_mask.vtp, locos_si_after.vtp, locos_ambient_after.vtp, "
               "locos_mask_after.vtp, and locos_oxidation_diagnostics.csv\n";

  return 0;
}
