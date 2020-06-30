#include <fstream>
#include <iostream>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPointData.hpp>

/**
  Minimal example showing how to serialize an lsDomain and deserialize.
  \example Serialize.cpp
*/

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();
  auto mesh = lsSmartPointer<lsMesh>::New();

  const double radius = 7.3;
  const hrleVectorType<double, D> centre(5., 0.);

  lsMakeGeometry<double, 2>(
      levelSet, lsSmartPointer<lsSphere<double, D>>::New(centre, radius))
      .apply();

  lsPointData &data = levelSet->getPointData();
  typename lsPointData::ScalarDataType scalars;
  typename lsPointData::VectorDataType vectors;
  for (unsigned i = 0; i < levelSet->getNumberOfPoints(); ++i) {
    scalars.push_back(i);
    vectors.push_back(
        typename lsPointData::VectorDataType::value_type({double(i)}));
  }

  data.insertNextScalarData(scalars, "myScalars");
  data.insertNextVectorData(vectors, "myVectors");

  {
    std::ofstream fout("test.lvst", std::ofstream::binary);
    levelSet->serialize(fout);
    fout.close();
  }

  {
    auto newLevelSet = lsSmartPointer<lsDomain<double, D>>::New();
    std::ifstream fin("test.lvst", std::ofstream::binary);
    newLevelSet->deserialize(fin);
    lsPointData &newData = newLevelSet->getPointData();
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
    fin.close();
  }

  return 0;
}
