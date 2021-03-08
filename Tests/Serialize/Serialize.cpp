#include <fstream>
#include <iostream>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPointData.hpp>
#include <lsReader.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsWriter.hpp>

/**
  Minimal example showing how to serialize an lsDomain and deserialize.
  \example Serialize.cpp
*/

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();
  auto mesh = lsSmartPointer<lsMesh<>>::New();

  const double radius = 7.3;
  const hrleVectorType<double, D> centre(5., 0.);

  lsMakeGeometry<double, 2>(
      levelSet, lsSmartPointer<lsSphere<double, D>>::New(centre, radius))
      .apply();

  lsPointData<double> &data = levelSet->getPointData();
  typename lsPointData<double>::ScalarDataType scalars;
  typename lsPointData<double>::VectorDataType vectors;
  for (unsigned i = 0; i < levelSet->getNumberOfPoints(); ++i) {
    scalars.push_back(i);
    vectors.push_back(
        typename lsPointData<double>::VectorDataType::value_type({double(i)}));
  }

  data.insertNextScalarData(scalars, "myScalars");
  data.insertNextVectorData(vectors, "myVectors");

  // {
  //   std::ofstream fout("test.lvst", std::ofstream::binary);
  //   levelSet->serialize(fout);
  //   fout.close();
  // }
  lsWriter<double, D>(levelSet, "test.lvst").apply();

  {
    auto newLevelSet = lsSmartPointer<lsDomain<double, D>>::New();
    // std::ifstream fin("test.lvst", std::ofstream::binary);
    lsReader<double, D>(newLevelSet, "test.lvst").apply();
    // newLevelSet->deserialize(fin);
    lsPointData<double> &newData = newLevelSet->getPointData();
    std::cout << newData.getScalarDataSize() << std::endl;
    auto newScalars = newData.getScalarData(0);
    std::cout << newData.getScalarDataLabel(0) << std::endl;
    for (auto i : *newScalars) {
      std::cout << i << std::endl;
    }
    auto newVectors = newData.getVectorData(0);
    std::cout << newData.getVectorDataLabel(0) << std::endl;
    for (auto i : *newVectors) {
      std::cout << i[0] << ", " << i[1] << ", " << i[2] << std::endl;
    }
    // fin.close();

    auto mesh = lsSmartPointer<lsMesh<>>::New();
    lsToMesh<double, D>(newLevelSet, mesh).apply();
    lsVTKWriter<double>(mesh, "test.vtk").apply();
  }

  return 0;
}
