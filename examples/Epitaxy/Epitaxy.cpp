#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

#include <vcTimer.hpp>

using namespace viennals;
constexpr int D = 3;
using T = double;

// Implement Epitaxy velocity field
class epitaxy final : public VelocityField<T> {
  std::vector<double> velocities;
  static constexpr double R111 = 0.5;
  static constexpr double R100 = 1.;
  static constexpr double low =
      (D > 2) ? 0.5773502691896257 : 0.7071067811865476;
  static constexpr double high = 1.0;

public:
  epitaxy(std::vector<double> vel) : velocities(vel){};

  double getScalarVelocity(const std::array<T, 3> & /*coordinate*/,
                           int material, const std::array<T, 3> &normal,
                           unsigned long /* pointID */) final {
    double vel = std::max(std::abs(normal[0]), std::abs(normal[2]));
    constexpr double factor = (R100 - R111) / (high - low);
    vel = (vel - low) * factor + R111;
    if (std::abs(normal[0]) < std::abs(normal[2])) {
      vel *= 2.;
    }

    return vel *
           ((material < int(velocities.size())) ? velocities[material] : 0);
  }
};

void writeSurface(SmartPointer<Domain<T, D>> domain,
                  const std::string &filename) {
  auto mesh = Mesh<T>::New();
  ToSurfaceMesh<T, D>(domain, mesh).apply();
  VTKWriter<T>(mesh, filename).apply();
}

int main(int argc, char *argv[]) {

  // Create hole geometry
  double bounds[2 * D] = {-1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
  BoundaryConditionEnum boundaryConditions[2 * D] = {
      BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      BoundaryConditionEnum::INFINITE_BOUNDARY};
  const T gridDelta = 0.03;

  const T finWidth = 0.5;
  const T finHeight = 0.2;

  auto mask = Domain<T, D>::New(bounds, boundaryConditions, gridDelta);
  MakeGeometry<T, D>(mask, Plane<T, D>::New(VectorType<T, D>{0.0, 0.0, 0.0},
                                            VectorType<T, D>{0.0, 0.0, 1.0}))
      .apply();

  auto substrate = Domain<T, D>::New(bounds, boundaryConditions, gridDelta);
  T minPoint[D] = {-finWidth / 2, -finWidth / 2, 0.};
  T maxPoint[D] = {finWidth / 2, finWidth / 2, finHeight};
  MakeGeometry<T, D>(substrate, Box<T, D>::New(minPoint, maxPoint)).apply();
  BooleanOperation(substrate, mask, BooleanOperationEnum::UNION).apply();

  writeSurface(mask, "mask.vtp");
  writeSurface(substrate, "substrate.vtp");

  std::vector<SmartPointer<Domain<T, D>>> levelSets;
  levelSets.push_back(mask);
  levelSets.push_back(substrate);

  PrepareStencilLocalLaxFriedrichs<T, D>(levelSets,
                                         std::vector<bool>{false, true});

  auto velocityField = SmartPointer<epitaxy>::New(std::vector<double>{0., -.5});

  Advect<T, D> advectionKernel(levelSets, velocityField);
  advectionKernel.setSpatialScheme(
      SpatialSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER);
  advectionKernel.setAdvectionTime(.5);

  Timer timer;
  timer.start();
  advectionKernel.apply();
  timer.finish();

  std::cout << "Epitaxy took " << timer.currentDuration / 1e9 << "s\n";

  FinalizeStencilLocalLaxFriedrichs<T, D>(levelSets);

  writeSurface(substrate, "epitaxy.vtp");
}