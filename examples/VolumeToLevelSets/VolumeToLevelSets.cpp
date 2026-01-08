#include <iostream>

#include <lsDomain.hpp>
#include <lsFromVolumeMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

namespace ls = viennals;

int main(int argc, char *argv[]) {

  using NumericType = double;
  constexpr int D = 3;
  double gridDelta = 0.1;

  std::string fileName;
  if (argc > 1) {
    fileName = std::string(argv[1]);
  } else {
    fileName = "volumeInitial.vtu";
    // Generate a simple volume mesh if no file is provided
    auto mesh = ls::SmartPointer<ls::Mesh<NumericType>>::New();
    // Create a simple cube [-1, 1]^3 consisting of 5 tetrahedra
    mesh->nodes = {{-1., -1., -1.}, {1., -1., -1.}, {1., 1., -1.},
                   {-1., 1., -1.},  {-1., -1., 1.}, {1., -1., 1.},
                   {1., 1., 1.},    {-1., 1., 1.}};
    // 5-tetra decomposition of a cube
    mesh->tetras = {
        {0, 1, 3, 4}, {1, 2, 3, 6}, {1, 4, 5, 6}, {3, 4, 6, 7}, {1, 3, 4, 6}};
    // Assign materials
    std::vector<double> materials = {0, 0, 1, 1, 1};
    mesh->cellData.insertNextScalarData(materials, "Material");

    ls::VTKWriter<NumericType>(mesh, ls::FileFormatEnum::VTU, fileName).apply();
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
        int id = std::round(cell);
        if (id >= 0 && static_cast<std::size_t>(id) < translator.size()) {
          cell = translator[id];
        }
      }
    }
  }

  mesh->print();

  ls::VTKWriter(mesh, ls::FileFormatEnum::VTU, "ReadVolumeMesh.vtu").apply();

  double bounds[2 * D] = {-2, 2, -2, 2, -2, 2};
  ls::BoundaryConditionEnum boundaryCons[D];
  for (unsigned i = 0; i < D; ++i) {
    boundaryCons[i] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
  }
  // boundaryCons[0] = ls::BoundaryConditionEnum::INFINITE_BOUNDARY;

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

    std::cout << "---------------------------------" << std::endl;
    std::cout << "LevelSet " << i << std::endl;
    levelSets[i]->print();

    auto meshVol = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToMesh<double, D>(levelSets[i], meshVol).apply();
    ls::VTKWriter<double>(meshVol, "LSvolume-" + std::to_string(i) + ".vtu")
        .apply();
  }

  return 0;
}