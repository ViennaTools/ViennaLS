#include <iostream>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPointData.hpp>
#include <lsReader.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsWriter.hpp>

/**
  Minimal example showing how to serialize an lsDomain and deserialize.
  \example Serialize.cpp
*/

int main() {
  constexpr int D = 3;

  omp_set_num_threads(4);

  auto levelSet = lsSmartPointer<lsDomain<double, D>>::New();

  const double radius = 7.3;
  const hrleVectorType<double, D> centre(5., 0.);

  lsMakeGeometry<double, D>(
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

  lsWriter<double, D>(levelSet, "test.lvst").apply();

  // read it in again
  auto newLevelSet = lsSmartPointer<lsDomain<double, D>>::New();
  lsReader<double, D>(newLevelSet, "test.lvst").apply();

  auto mesh = lsSmartPointer<lsMesh<>>::New();
  lsToSurfaceMesh<double, D>(levelSet, mesh).apply();
  lsVTKWriter<double>(mesh, "test.vtk").apply();

  return 0;
}
