#include <iostream>
#include <sstream>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPointData.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

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

  std::stringstream stream;
  levelSet->serialize(stream);

  {
    auto newLevelSet = lsSmartPointer<lsDomain<double, D>>::New();
    newLevelSet->deserialize(stream);

    if (newLevelSet->getNumberOfPoints() != levelSet->getNumberOfPoints()) {
      std::cout << "Levelset points were not properly deserialized"
                << std::endl;
      return EXIT_FAILURE;
    }

    lsPointData<double> &newData = newLevelSet->getPointData();
    auto scalarDataSize = newData.getScalarDataSize();
    if (scalarDataSize != levelSet->getPointData().getScalarDataSize()) {
      std::cout << "Scalar data was not properly deserialized" << std::endl;
      return EXIT_FAILURE;
    }

    auto newScalars = newData.getScalarData(0);
    std::cout << newData.getScalarDataLabel(0) << std::endl;
    for (auto i : *newScalars) {
      std::cout << i << std::endl;
    }

    auto vectorDataSize = newData.getVectorDataSize();
    if (vectorDataSize != levelSet->getPointData().getVectorDataSize()) {
      std::cout << "Vector data was not properly deserialized" << std::endl;
      return EXIT_FAILURE;
    }

    auto newVectors = newData.getVectorData(0);
    std::cout << newData.getVectorDataLabel(0) << std::endl;
    for (auto i : *newVectors) {
      std::cout << i[0] << ", " << i[1] << ", " << i[2] << std::endl;
    }
  }

  return EXIT_SUCCESS;
}
