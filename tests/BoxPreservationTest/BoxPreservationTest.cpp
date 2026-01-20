#include <iostream>
#include <string>

#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMarchingCubes.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

template <int D> void runTest() {
  using T = double;

  // 1. Setup Domain
  double gridDelta = (D == 2) ? 0.1 : 0.1;
  // Define bounds large enough for the box
  double bounds[2 * D];
  for (int i = 0; i < 2 * D; ++i)
    bounds[i] = (i % 2 == 0) ? -3.0 : 3.0;

  typename viennals::Domain<T, D>::BoundaryType boundaryCons[D];
  for (int i = 0; i < D; ++i)
    boundaryCons[i] = viennals::Domain<T, D>::BoundaryType::INFINITE_BOUNDARY;

  auto domain = viennals::SmartPointer<viennals::Domain<T, D>>::New(
      bounds, boundaryCons, gridDelta);

  // 2. Create Box Geometry
  // Box from (-1, -1, -1) to (1, 1, 1)
  T minCorner[D];
  T maxCorner[D];
  for (int i = 0; i < D; ++i) {
    minCorner[i] = -1.0;
    maxCorner[i] = 1.0;
  }

  auto box =
      viennals::SmartPointer<viennals::Box<T, D>>::New(minCorner, maxCorner);

  std::cout << "--- Running " << D << "D Test ---" << std::endl;
  std::cout << "Initial Box: Min(-1...), Max(1...)" << std::endl;

  std::cout << "Creating Box Level Set..." << std::endl;
  viennals::MakeGeometry<T, D>(domain, box).apply();

  std::cout << "Saving Level Set..." << std::endl;
  auto lsMesh = viennals::SmartPointer<viennals::Mesh<T>>::New();
  viennals::ToMesh<T, D>(domain, lsMesh).apply();
  std::string vtuName = "BoxLevelSet_" + std::to_string(D) + "D.vtu";
  viennals::VTKWriter<T>(lsMesh, vtuName).apply();

  // 3. Convert Level Set back to Mesh
  std::cout << "Converting Level Set to Mesh..." << std::endl;
  auto mesh = viennals::SmartPointer<viennals::Mesh<T>>::New();
  viennals::ToSurfaceMesh<T, D>(domain, mesh).apply();

  // 4. Output Mesh statistics
  std::cout << "Mesh Nodes: " << mesh->getNodes().size() << std::endl;
  if constexpr (D == 2) {
    std::cout << "Mesh Lines: " << mesh->template getElements<2>().size()
              << std::endl;
  } else {
    std::cout << "Mesh Triangles: " << mesh->template getElements<3>().size()
              << std::endl;
  }

  // 5. Write to file
  std::string filename = "BoxFinal_" + std::to_string(D) + "D.vtp";
  viennals::VTKWriter<T>(mesh, filename).apply();
  std::cout << "Written mesh to " << filename << std::endl;

  // 6. Test the exposed Marching Cubes edge tables
  if constexpr (D == 2) {
    // Case 1: Only corner 0 is inside.
    // Edges connected to corner 0 are 0 and 3.
    // Expected bitmask: (1<<0) | (1<<3) = 1 | 8 = 9.
    unsigned int signs = 1; // Corner 0 inside
    unsigned int edges =
        lsInternal::MarchingCubes::getIntersectedEdges2d(signs);

    if (edges == 9) {
      std::cout << "Marching Squares Edge Table Test: PASSED" << std::endl;
    } else {
      std::cout << "Marching Squares Edge Table Test: FAILED (Expected 9, got "
                << edges << ")" << std::endl;
    }
  } else {
    // Case 1: Only corner 0 is inside.
    // Edges connected to corner 0 are 0, 3, 8.
    // Expected bitmask: (1<<0) | (1<<3) | (1<<8) = 1 | 8 | 256 = 265.
    unsigned int signs = 1; // Corner 0 inside
    unsigned int edges =
        lsInternal::MarchingCubes::getIntersectedEdges3d(signs);

    if (edges == 265) {
      std::cout << "Marching Cubes Edge Table Test: PASSED" << std::endl;
    } else {
      std::cout << "Marching Cubes Edge Table Test: FAILED (Expected 265, got "
                << edges << ")" << std::endl;
    }
  }
  std::cout << std::endl;
}

int main() {
  runTest<2>();
  runTest<3>();
  return 0;
}
