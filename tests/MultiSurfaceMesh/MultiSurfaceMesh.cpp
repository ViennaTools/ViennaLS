#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToMesh.hpp>
#include <lsToMultiSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

using namespace viennals;

class velocityField : public VelocityField<double> {
public:
  Vec3D<double> getVectorVelocity(const Vec3D<double> &, int material,
                                  const Vec3D<double> &, unsigned long) final {
    if (material == 3)
      return {0., 0., 0.};
    return {0., -1., 0.};
  }
};

int main() {

  double gridDelta = 0.17;
  double extent[4] = {-10.0, 10.0, -10.0, 10.0};
  BoundaryConditionEnum bc[2] = {BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
                                 BoundaryConditionEnum::INFINITE_BOUNDARY};

  std::vector<SmartPointer<Domain<double, 2>>> layers;
  double origin[2] = {0.0, 0.0};
  double normal[2] = {0.0, 1.0};

  {
    auto substrate = Domain<double, 2>::New(extent, bc, gridDelta);
    MakeGeometry<double, 2>(substrate, Plane<double, 2>::New(origin, normal))
        .apply();
    layers.push_back(substrate);
  }

  for (unsigned i = 0; i < 2; ++i) {
    auto layer = Domain<double, 2>::New(extent, bc, gridDelta);
    origin[1] += 1.0;
    MakeGeometry<double, 2>(layer, Plane<double, 2>::New(origin, normal))
        .apply();
    layers.push_back(layer);
  }

  {
    auto box = Domain<double, 2>::New(extent, bc, gridDelta);
    double minPoint[2] = {-2.5, origin[1] - 0.5 * gridDelta};
    double maxPoint[2] = {2.5, origin[1] + 5.0};
    MakeGeometry<double, 2>(box, Box<double, 2>::New(minPoint, maxPoint))
        .apply();
    BooleanOperation<double, 2>(box, layers.back(), BooleanOperationEnum::UNION)
        .apply();
    layers.push_back(box);
  }

  {
    Advect<double, 2> advectionKernel;
    for (const auto &layer : layers) {
      advectionKernel.insertNextLevelSet(layer);
    }
    advectionKernel.setVelocityField(SmartPointer<velocityField>::New());
    advectionKernel.setAdvectionTime(3.0);
    advectionKernel.apply();
  }

  auto mesh = Mesh<double>::New();
  ToMultiSurfaceMesh<double, 2> mesher(layers, mesh);
  mesher.apply();

  std::cout << "Writing regular marching cubes mesh..." << std::endl;
  VTKWriter<double>(mesh, "multi_surface_mesh_no_sharp.vtp").apply();

  mesher.setSharpCorners(true);
  mesher.apply();

  std::cout << "Writing sharp corner marching cubes mesh..." << std::endl;
  VTKWriter<double>(mesh, "multi_surface_mesh_sharp.vtp").apply();

  std::cout << "Writing level set meshes..." << std::endl;
  int i = 0;
  for (const auto &layer : layers) {
    auto lsmesh = Mesh<double>::New();
    ToMesh<double, 2>(layer, lsmesh).apply();
    VTKWriter<double>(lsmesh, "layer_" + std::to_string(i) + ".vtp").apply();
    ++i;
  }
}