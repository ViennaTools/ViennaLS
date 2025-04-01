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

namespace ls = viennals;

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  auto levelSet = ls::SmartPointer<ls::Domain<double, D>>::New();
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();

  const double radius = 7.3;
  const ls::VectorType<double, D> centre{5., 0.};

  ls::MakeGeometry<double, 2>(
      levelSet, ls::SmartPointer<ls::Sphere<double, D>>::New(centre, radius))
      .apply();

  ls::PointData<double> &data = levelSet->getPointData();
  typename ls::PointData<double>::ScalarDataType scalars;
  typename ls::PointData<double>::VectorDataType vectors;
  for (unsigned i = 0; i < levelSet->getNumberOfPoints(); ++i) {
    scalars.push_back(i);
    vectors.push_back(
        typename ls::PointData<double>::VectorDataType::value_type(
            {double(i)}));
  }

  data.insertNextScalarData(scalars, "myScalars");
  data.insertNextVectorData(vectors, "myVectors");

  std::stringstream stream;
  levelSet->serialize(stream);

  {
    auto newLevelSet = ls::SmartPointer<ls::Domain<double, D>>::New();
    newLevelSet->deserialize(stream);

    if (newLevelSet->getNumberOfPoints() != levelSet->getNumberOfPoints()) {
      std::cout << "Levelset points were not properly deserialized"
                << std::endl;
      return EXIT_FAILURE;
    }

    ls::PointData<double> &newData = newLevelSet->getPointData();
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
