#include <iostream>

#include <lsDomain.hpp>
#include <lsFromVolumeMesh.hpp>
#include <lsSmartPointer.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKReader.hpp>
#include <lsVTKWriter.hpp>

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

  auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();
  lsVTKReader(mesh, lsFileFormatEnum::VTU, fileName).apply();

  // reorder numbering
  {
    typename lsPointData<NumericType>::ScalarDataType *materialData =
        mesh->getCellData().getScalarData("Material");

    std::vector<int> translator = {3, 2, 4, 7, 7, 6, 5, 7, 1, 0};
    if (materialData == nullptr) {
      std::cout << "Could not get material data" << std::endl;
    } else {
      for (auto &cell : *materialData) {
        cell = translator[std::round(cell)];
      }
    }
  }

  mesh->print();

  lsVTKWriter(mesh, lsFileFormatEnum::VTU, "ReadVolumeMesh.vtu").apply();

  double bounds[2 * D] = {-6, 6, 1e-10, 0.078, -0.034, 0.034};
  lsBoundaryConditionEnum<D> boundaryCons[D];
  for (unsigned i = 0; i < D; ++i) {
    boundaryCons[i] = lsBoundaryConditionEnum<D>::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[0] = lsBoundaryConditionEnum<D>::INFINITE_BOUNDARY;

  auto domain = lsSmartPointer<lsDomain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);

  // read in as LS
  lsFromVolumeMesh<double, D> reader(domain->getGrid(), mesh);
  reader.apply();
  auto levelSets = reader.getLevelSets();

  for (unsigned i = 0; i < levelSets.size(); ++i) {
    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToSurfaceMesh<double, D>(levelSets[i], mesh).apply();
    lsVTKWriter<double>(mesh, "LSsurface-" + std::to_string(i) + ".vtk")
        .apply();
  }

  return 0;
}