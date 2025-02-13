#include <iostream>

#include <lsDomain.hpp>
#include <lsFromVolumeMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

int main(int argc, char *argv[]) {

  using NumericType = double;
  constexpr int D = 3;
  double gridDelta = 0.00023;

  std::string fileName;
  if (argc > 1) {
    fileName = std::string(argv[1]);
  } else {
    fileName = "volumeInitial.vtu";
  }

  auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
  ls::VTKReader(mesh, ls::FileFormatEnum::VTU, fileName).apply();

  // reorder numbering
  {
    typename ls::PointData<NumericType>::ScalarDataType *materialData =
        mesh->getCellData().getScalarData("Material");

    std::vector<int> translator = {3, 2, 4, 7, 7, 6, 5, 7, 1, 0};
    if (materialData != nullptr) {
      for (auto &cell : *materialData) {
        cell = translator[std::round(cell)];
      }
    }
  }

  mesh->print();

  ls::VTKWriter(mesh, ls::FileFormatEnum::VTU, "ReadVolumeMesh.vtu").apply();

  double bounds[2 * D] = {-6, 6, 1e-10, 0.078, -0.034, 0.034};
  ls::BoundaryConditionEnum boundaryCons[D];
  for (unsigned i = 0; i < D; ++i) {
    boundaryCons[i] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[0] = ls::BoundaryConditionEnum::INFINITE_BOUNDARY;

  auto domain = ls::SmartPointer<ls::Domain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  // read in as LS
  ls::FromVolumeMesh<double, D> reader(domain->getGrid(), mesh);
  reader.apply();
  auto levelSets = reader.getLevelSets();

  for (unsigned i = 0; i < levelSets.size(); ++i) {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(levelSets[i], mesh).apply();
    ls::VTKWriter<double>(mesh, "LSsurface-" + std::to_string(i) + ".vtp")
        .apply();
  }

  return 0;
}