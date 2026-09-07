#include <lsMesh.hpp>
#include <vcTestAsserts.hpp>

#include <cmath>
#include <limits>

template <class T> void testEmptyAndUnique() {
  viennals::Mesh<T> mesh;
  mesh.removeDuplicateNodes();
  VC_TEST_ASSERT(mesh.nodes.empty());

  mesh.nodes = {{1, 2, 3}};
  mesh.pointData.insertNextScalarData({T(7)}, "Values");
  mesh.removeDuplicateNodes();
  VC_TEST_ASSERT(mesh.nodes.size() == 1);
  VC_TEST_ASSERT(mesh.pointData.getScalarData("Values")->at(0) == T(7));

  mesh.nodes.push_back({4, 5, 6});
  mesh.pointData.getScalarData("Values")->push_back(T(8));
  mesh.lines = {{1, 0}};
  const auto originalNodes = mesh.nodes;
  const auto originalLines = mesh.lines;
  const auto originalValues = *mesh.pointData.getScalarData("Values");
  mesh.removeDuplicateNodes();
  VC_TEST_ASSERT(mesh.nodes == originalNodes);
  VC_TEST_ASSERT(mesh.lines == originalLines);
  VC_TEST_ASSERT(*mesh.pointData.getScalarData("Values") == originalValues);
}

template <class T> void testConnectivityAndData() {
  viennals::Mesh<T> mesh;
  mesh.nodes = {{0, 0, 0}, {1, 0, 0}, {0, 0, 0}, {0, 1, 0}, {1, 0, 0},
                {0, 0, 1}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
  mesh.vertices = {{2}, {9}};
  mesh.lines = {{4, 3}, {0, 2}};
  mesh.triangles = {{2, 3, 9}};
  mesh.tetras = {{0, 4, 5, 9}};
  mesh.hexas = {{0, 1, 3, 5, 6, 7, 8, 9}};
  mesh.minimumExtent = {0, 0, 0};
  mesh.maximumExtent = {1, 1, 1};
  mesh.pointData.insertNextScalarData({10, 20, 99, 30, 98, 40, 50, 60, 70, 80},
                                      "Values");
  std::vector<viennals::Vec3D<T>> vectors;
  for (const T value : *mesh.pointData.getScalarData("Values"))
    vectors.push_back({value, 0, -value});
  mesh.pointData.insertNextVectorData(vectors, "Vectors");
  const std::vector<T> cellValues = {1, 2, 3, 4, 5, 6, 7};
  mesh.cellData.insertNextScalarData(cellValues, "Materials");

  viennals::Mesh<T> expected;
  expected.nodes = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
                    {1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
  expected.vertices = {{0}, {7}};
  expected.lines = {{1, 2}, {0, 0}};
  expected.triangles = {{0, 2, 7}};
  expected.tetras = {{0, 1, 3, 7}};
  expected.hexas = {{0, 1, 2, 3, 4, 5, 6, 7}};
  const std::vector<T> expectedValues = {10, 20, 30, 40, 50, 60, 70, 80};
  std::vector<viennals::Vec3D<T>> expectedVectors;
  for (const T value : expectedValues)
    expectedVectors.push_back({value, 0, -value});

  // Verify both the first application and idempotence.
  for (int pass = 0; pass < 2; ++pass) {
    mesh.removeDuplicateNodes();
    VC_TEST_ASSERT(mesh.nodes == expected.nodes);
    VC_TEST_ASSERT(mesh.vertices == expected.vertices);
    VC_TEST_ASSERT(mesh.lines == expected.lines);
    VC_TEST_ASSERT(mesh.triangles == expected.triangles);
    VC_TEST_ASSERT(mesh.tetras == expected.tetras);
    VC_TEST_ASSERT(mesh.hexas == expected.hexas);
    VC_TEST_ASSERT(*mesh.pointData.getScalarData("Values") == expectedValues);
    VC_TEST_ASSERT(*mesh.pointData.getVectorData("Vectors") == expectedVectors);
    VC_TEST_ASSERT(*mesh.cellData.getScalarData("Materials") == cellValues);
    VC_TEST_ASSERT(mesh.minimumExtent == expected.nodes.front());
    VC_TEST_ASSERT(mesh.maximumExtent == expected.nodes.back());
  }
}

template <class T> void testFloatingPointEquality() {
  const T inf = std::numeric_limits<T>::infinity();
  const T nan = std::numeric_limits<T>::quiet_NaN();
  viennals::Mesh<T> mesh;
  mesh.nodes = {{-T(0), 0, 0},
                {T(0), 0, 0},
                {inf, 0, 0},
                {inf, 0, 0},
                {-inf, 0, 0},
                {nan, 0, 0},
                {nan, 0, 0},
                {0, nan, 0},
                {0, nan, 0},
                {0, 0, nan},
                {0, 0, nan},
                {1, 0, 0},
                {std::nextafter(T(1), T(2)), 0, 0}};
  for (unsigned i = 0; i < mesh.nodes.size(); ++i)
    mesh.vertices.push_back({i});
  const std::vector<std::array<unsigned, 1>> expectedVertices = {
      {0}, {0}, {1}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}};
  mesh.removeDuplicateNodes();
  VC_TEST_ASSERT(mesh.nodes.size() == 11);
  VC_TEST_ASSERT(mesh.vertices == expectedVertices);
  VC_TEST_ASSERT(std::signbit(mesh.nodes[0][0]));
  VC_TEST_ASSERT(mesh.nodes[1][0] == inf);
  VC_TEST_ASSERT(mesh.nodes[2][0] == -inf);
  VC_TEST_ASSERT(std::isnan(mesh.nodes[3][0]));
  VC_TEST_ASSERT(std::isnan(mesh.nodes[5][1]));
  VC_TEST_ASSERT(std::isnan(mesh.nodes[7][2]));
  VC_TEST_ASSERT(mesh.nodes[9][0] != mesh.nodes[10][0]);
}

template <class T> void testAllDuplicatesAndInvalidData() {
  viennals::Mesh<T> mesh;
  mesh.nodes.assign(100, {1, 2, 3});
  mesh.lines = {{99, 50}};
  // Reject incomplete data without partially modifying the mesh.
  mesh.pointData.insertNextScalarData({T(7)}, "Values");
  const auto originalNodes = mesh.nodes;
  const auto originalLines = mesh.lines;
  bool rejected = false;
  try {
    mesh.removeDuplicateNodes();
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  VC_TEST_ASSERT(rejected);
  VC_TEST_ASSERT(mesh.nodes == originalNodes);
  VC_TEST_ASSERT(mesh.lines == originalLines);
  VC_TEST_ASSERT(mesh.pointData.getScalarData("Values")->size() == 1);

  mesh.pointData.getScalarData("Values")->resize(100, T(99));
  mesh.removeDuplicateNodes();
  VC_TEST_ASSERT(mesh.nodes.size() == 1);
  VC_TEST_ASSERT(mesh.lines.front()[0] == 0 && mesh.lines.front()[1] == 0);
  VC_TEST_ASSERT(mesh.pointData.getScalarData("Values")->size() == 1);
  VC_TEST_ASSERT(mesh.pointData.getScalarData("Values")->front() == T(7));
}

template <class T> void runTests() {
  testEmptyAndUnique<T>();
  testConnectivityAndData<T>();
  testFloatingPointEquality<T>();
  testAllDuplicatesAndInvalidData<T>();
}

int main() {
  runTests<float>();
  runTests<double>();
}
