#include <iostream>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPointData.hpp>
#include <lsReader.hpp>
#include <lsTestAsserts.hpp>
#include <lsWriter.hpp>

/**
  Minimal example showing how to serialize an lsDomain and deserialize.
  \example Serialize.cpp
*/

int main() {
  constexpr int D = 2;

  omp_set_num_threads(4);

  double extent = 10;
  double gridDelta = 1;
  double bounds[2 * D] = {-extent, extent, -extent, extent};
  lsDomain<double, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i) {
    boundaryCons[i] = lsDomain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] = lsBoundaryConditionEnum<D>::INFINITE_BOUNDARY;

  auto levelSet =
      lsSmartPointer<lsDomain<double, D>>::New(bounds, boundaryCons, gridDelta);

  const double radius = 7.3;
  const hrleVectorType<double, D> centre(5., 0.);
  // const double centre[D] = {0,0.1};
  // const double normal[D] = {0, 1};

  // lsMakeGeometry<double, D>(levelSet, lsSmartPointer<lsPlane<double,
  // D>>::New(centre, normal)).apply();
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

  // TODO: using the output of hrleDomain::print here to compare data in LS
  // there should definitely be a better way
  std::string domainString;
  {
    std::ostringstream oss;
    levelSet->getDomain().print(oss);
    domainString = oss.str();
  }

  lsWriter<double, D>(levelSet, "test.lvst").apply();

  // read it in again
  auto newLevelSet = lsSmartPointer<lsDomain<double, D>>::New();
  lsReader<double, D>(newLevelSet, "test.lvst").apply();

  {
    std::ostringstream oss;
    newLevelSet->getDomain().print(oss);
    LSTEST_ASSERT(oss.str() == domainString)
  }

  return 0;
}
