#include <iostream>
#include <string>

#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMarchingCubes.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

template <int D> void runTest() {
  using T = double;

  // 1. Setup Domain
  double gridDelta = (D == 2) ? 0.0485 : 0.0485;
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
  std::string vtuName = "BoxFinal_" + std::to_string(D) + "D.vtu";
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

  // 6. Create Sphere Geometry
  auto domainSphere = viennals::SmartPointer<viennals::Domain<T, D>>::New(
      bounds, boundaryCons, gridDelta);
  T origin[D];
  for (int i = 0; i < D; ++i)
    origin[i] = 0;
  T radius = 1.0;
  auto sphere =
      viennals::SmartPointer<viennals::Sphere<T, D>>::New(origin, radius);

  std::cout << "Creating Sphere Level Set..." << std::endl;
  viennals::MakeGeometry<T, D>(domainSphere, sphere).apply();

  std::cout << "Converting Sphere Level Set to Mesh..." << std::endl;
  auto meshSphere = viennals::SmartPointer<viennals::Mesh<T>>::New();
  viennals::ToSurfaceMesh<T, D>(domainSphere, meshSphere).apply();

  std::string filenameSphere = "SphereFinal_" + std::to_string(D) + "D.vtp";
  viennals::VTKWriter<T>(meshSphere, filenameSphere).apply();
  viennals::ToMesh<T, D>(domainSphere, meshSphere).apply();

  filenameSphere = "SphereFinal_" + std::to_string(D) + "D.vtu";
  viennals::VTKWriter<T>(meshSphere, filenameSphere).apply();
  std::cout << "Written mesh to " << filenameSphere << std::endl;

  // 7. Test the exposed Marching Cubes edge tables
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

  // 8. Create Plane with Sphere Cavity
  for (int i = 0; i < D - 1; ++i)
    boundaryCons[i] =
        viennals::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  boundaryCons[D - 1] =
      viennals::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

  auto substrate = viennals::SmartPointer<viennals::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  // double origin[3] = {0., 0., 0.};
  {
    T planeNormal[D];
    for (int i = 0; i < D; ++i)
      planeNormal[i] = 0.0;
    planeNormal[D - 1] = 1.0;
    auto plane =
        viennals::SmartPointer<viennals::Plane<double, D>>::New(origin,
                                                                planeNormal);
    viennals::MakeGeometry<double, D>(substrate, plane).apply();
  }

  auto meshPlane = viennals::SmartPointer<viennals::Mesh<T>>::New();
  viennals::ToSurfaceMesh<T, D>(substrate, meshPlane).apply();
  std::string filenamePlane = "CavityPlane_" + std::to_string(D) + "D.vtp";
  viennals::VTKWriter<T>(meshPlane, filenamePlane).apply();
  viennals::ToMesh<T, D>(substrate, meshPlane).apply();
  filenamePlane = "CavityPlane_" + std::to_string(D) + "D.vtu";
  viennals::VTKWriter<T>(meshPlane, filenamePlane).apply();

  {
    // create spheres used for booling
    std::cout << "Creating sphere..." << std::endl;
    auto sphere = viennals::SmartPointer<viennals::Domain<double, D>>::New(
        bounds, boundaryCons, gridDelta);

    origin[D - 1] = -0.6;
    double radius = 1.0;
    viennals::MakeGeometry<double, D>(
        sphere, viennals::SmartPointer<viennals::Sphere<double, D>>::New(origin, radius))
        .apply();
  
    auto meshSphereCavity = viennals::SmartPointer<viennals::Mesh<T>>::New();
    viennals::ToSurfaceMesh<T, D>(sphere, meshSphereCavity).apply();
    std::string filenameSphereCavity = "CavitySphere_" + std::to_string(D) + "D.vtp";
    viennals::VTKWriter<T>(meshSphereCavity, filenameSphereCavity).apply();
    viennals::ToMesh<T, D>(sphere, meshSphereCavity).apply();
    filenameSphereCavity = "CavitySphere_" + std::to_string(D) + "D.vtu";
    viennals::VTKWriter<T>(meshSphereCavity, filenameSphereCavity).apply();

        
    viennals::BooleanOperation<double, D> boolOp(
        substrate, sphere, viennals::BooleanOperationEnum::RELATIVE_COMPLEMENT);
    boolOp.apply();

    std::cout << "Converting Cavity Level Set to Mesh..." << std::endl;
    auto meshCavity = viennals::SmartPointer<viennals::Mesh<T>>::New();
    viennals::ToSurfaceMesh<T, D>(substrate, meshCavity).apply();

    std::string filenameCavity = "CavityFinal_" + std::to_string(D) + "D.vtp";
    viennals::VTKWriter<T>(meshCavity, filenameCavity).apply();
    viennals::ToMesh<T, D>(substrate, meshCavity).apply();

    filenameCavity = "CavityFinal_" + std::to_string(D) + "D.vtu";
    viennals::VTKWriter<T>(meshCavity, filenameCavity).apply();
    std::cout << "Written mesh to " << filenameCavity << std::endl;

    std::cout << std::endl;
  }

  // 9. Create Plane with Box Cavity
  {
    std::cout << "--- Box Cavity Test ---" << std::endl;
    // Reset substrate
    substrate = viennals::SmartPointer<viennals::Domain<double, D>>::New(
        bounds, boundaryCons, gridDelta);

    T planeNormal[D];
    for (int i = 0; i < D; ++i)
      planeNormal[i] = 0.0;
    planeNormal[D - 1] = 1.0;

    // Reset origin
    for (int i = 0; i < D; ++i) origin[i] = 0.0;

    auto plane =
        viennals::SmartPointer<viennals::Plane<double, D>>::New(origin,
                                                                planeNormal);
    viennals::MakeGeometry<double, D>(substrate, plane).apply();

    // Create box used for booling
    std::cout << "Creating box..." << std::endl;
    auto boxDomain = viennals::SmartPointer<viennals::Domain<double, D>>::New(
        bounds, boundaryCons, gridDelta);

    T minBox[D];
    T maxBox[D];
    for (int i = 0; i < D; ++i) {
      minBox[i] = -1.0;
      maxBox[i] = 1.0;
    }

    viennals::MakeGeometry<double, D>(
        boxDomain, viennals::SmartPointer<viennals::Box<double, D>>::New(minBox, maxBox))
        .apply();

    viennals::BooleanOperation<double, D> boolOp(
        substrate, boxDomain, viennals::BooleanOperationEnum::RELATIVE_COMPLEMENT);
    boolOp.apply();

    std::cout << "Converting Box Cavity Level Set to Mesh..." << std::endl;
    auto meshBoxCavity = viennals::SmartPointer<viennals::Mesh<T>>::New();
    viennals::ToSurfaceMesh<T, D>(substrate, meshBoxCavity).apply();

    std::string filenameBoxCavity = "CavityBoxFinal_" + std::to_string(D) + "D.vtp";
    viennals::VTKWriter<T>(meshBoxCavity, filenameBoxCavity).apply();

    viennals::ToMesh<T, D>(substrate, meshBoxCavity).apply();
    filenameBoxCavity = "CavityBoxFinal_" + std::to_string(D) + "D.vtu";
    viennals::VTKWriter<T>(meshBoxCavity, filenameBoxCavity).apply();
    std::cout << "Written mesh to " << filenameBoxCavity << std::endl;
    std::cout << std::endl;
  }
}

int main() {
  // omp_set_num_threads(1);
  runTest<2>();
  // runTest<3>();
  return 0;
}
