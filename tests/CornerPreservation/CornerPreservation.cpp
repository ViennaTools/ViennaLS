#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <lsBooleanOperation.hpp>
#include <lsCompareChamfer.hpp>
#include <lsDomain.hpp>
#include <lsMakeGeometry.hpp>
#include <lsMarchingCubes.hpp>
#include <lsTestAsserts.hpp>
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

  std::cout << "--- Running " << D << "D Test ---" << std::endl;
  // 2. Create Box Geometry
  {
    std::cout << std::endl << "----- Box Test -----" << std::endl;
    typename viennals::Domain<T, D>::BoundaryType boundaryCons[D];
    for (int i = 0; i < D; ++i)
      boundaryCons[i] = viennals::Domain<T, D>::BoundaryType::INFINITE_BOUNDARY;

    auto domain = viennals::SmartPointer<viennals::Domain<T, D>>::New(
        bounds, boundaryCons, gridDelta);

    // Box from (-1, -1, -1) to (1, 1, 1)
    T minCorner[D];
    T maxCorner[D];
    for (int i = 0; i < D; ++i) {
      minCorner[i] = -1.0;
      maxCorner[i] = 1.0;
    }

    auto box =
        viennals::SmartPointer<viennals::Box<T, D>>::New(minCorner, maxCorner);

    std::cout << "Creating Box Level Set..." << std::endl;
    viennals::MakeGeometry<T, D>(domain, box).apply();

    // 3. Convert Level Set back to Mesh
    std::cout << "Converting Level Set to Mesh..." << std::endl;
    auto mesh = viennals::SmartPointer<viennals::Mesh<T>>::New();
    {
      viennals::ToSurfaceMesh<T, D> toSurfaceMesh(domain, mesh);
      toSurfaceMesh.setSharpCorners(true);
      toSurfaceMesh.apply();
    }

    // 5. Write to file
    std::string filename = "BoxFinal_" + std::to_string(D) + "D.vtp";
    viennals::VTKWriter<T>(mesh, filename).apply();
    std::cout << "Written mesh to " << filename << std::endl;

    std::vector<std::vector<T>> expectedCorners;
    for (int i = 0; i < (1 << D); ++i) {
      std::vector<T> corner(D);
      for (int j = 0; j < D; ++j) {
        corner[j] = ((i >> j) & 1) ? 1.0 : -1.0;
      }
      expectedCorners.push_back(corner);
    }
    LSTEST_ASSERT_MESH_CORNERS(mesh, expectedCorners, D);
  }

  // 6. Create Sphere Geometry
  {
    std::cout << std::endl << "----- Sphere Test -----" << std::endl;
    typename viennals::Domain<T, D>::BoundaryType boundaryCons[D];
    for (int i = 0; i < D; ++i)
      boundaryCons[i] = viennals::Domain<T, D>::BoundaryType::INFINITE_BOUNDARY;

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
    {
      viennals::ToSurfaceMesh<T, D> toSurfaceMesh(domainSphere, meshSphere);
      toSurfaceMesh.setSharpCorners(true);
      toSurfaceMesh.apply();
    }
    std::string filenameSphere = "SphereFinal_" + std::to_string(D) + "D.vtp";
    viennals::VTKWriter<T>(meshSphere, filenameSphere).apply();
    std::cout << "Written mesh to " << filenameSphere << std::endl;
  }

  // 7. Create Plane with Sphere Cavity
  {
    std::cout << std::endl << "----- Cavity Test -----" << std::endl;
    typename viennals::Domain<T, D>::BoundaryType boundaryCons[D];
    for (int i = 0; i < D - 1; ++i)
      boundaryCons[i] =
          viennals::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
    boundaryCons[D - 1] =
        viennals::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

    auto substrate = viennals::SmartPointer<viennals::Domain<double, D>>::New(
        bounds, boundaryCons, gridDelta);

    T origin[D];
    for (int i = 0; i < D; ++i)
      origin[i] = 0.0;

    {
      T planeNormal[D];
      for (int i = 0; i < D; ++i)
        planeNormal[i] = 0.0;
      planeNormal[D - 1] = 1.0;
      auto plane = viennals::SmartPointer<viennals::Plane<double, D>>::New(
          origin, planeNormal);
      std::cout << "Creating Cavity Level Set..." << std::endl;
      viennals::MakeGeometry<double, D>(substrate, plane).apply();
    }

    {
      // create spheres used for booling
      auto sphere = viennals::SmartPointer<viennals::Domain<double, D>>::New(
          bounds, boundaryCons, gridDelta);

      origin[D - 1] = -0.6;
      double radius = 1.0;
      viennals::MakeGeometry<double, D>(
          sphere, viennals::SmartPointer<viennals::Sphere<double, D>>::New(
                      origin, radius))
          .apply();

      viennals::BooleanOperation<double, D> boolOp(
          substrate, sphere,
          viennals::BooleanOperationEnum::RELATIVE_COMPLEMENT);
      boolOp.apply();

      std::cout << "Converting Cavity Level Set to Mesh..." << std::endl;
      auto meshCavity = viennals::SmartPointer<viennals::Mesh<T>>::New();
      {
        viennals::ToSurfaceMesh<T, D> toSurfaceMesh(substrate, meshCavity);
        toSurfaceMesh.setSharpCorners(true);
        toSurfaceMesh.apply();
      }
      std::string filenameCavity = "CavityFinal_" + std::to_string(D) + "D.vtp";
      viennals::VTKWriter<T>(meshCavity, filenameCavity).apply();
      std::cout << "Written mesh to " << filenameCavity << std::endl;

      std::cout << std::endl;
    }
  }

  // 8. Create Plane with Box Cavity
  {
    std::cout << std::endl << "----- Box Cavity Test -----" << std::endl;

    typename viennals::Domain<T, D>::BoundaryType boundaryCons[D];
    for (int i = 0; i < D - 1; ++i)
      boundaryCons[i] =
          viennals::Domain<double, D>::BoundaryType::REFLECTIVE_BOUNDARY;
    boundaryCons[D - 1] =
        viennals::Domain<double, D>::BoundaryType::INFINITE_BOUNDARY;

    // Reset substrate
    auto substrate = viennals::SmartPointer<viennals::Domain<double, D>>::New(
        bounds, boundaryCons, gridDelta);

    T origin[D];
    for (int i = 0; i < D; ++i)
      origin[i] = 0.0;
    origin[D - 1] = 0.025;

    T planeNormal[D];
    for (int i = 0; i < D; ++i)
      planeNormal[i] = 0.0;
    planeNormal[D - 1] = 1.0;

    auto plane = viennals::SmartPointer<viennals::Plane<double, D>>::New(
        origin, planeNormal);
    std::cout << "Creating Box Cavity Level Set..." << std::endl;
    viennals::MakeGeometry<double, D>(substrate, plane).apply();

    // Create box used for booling
    auto boxDomain = viennals::SmartPointer<viennals::Domain<double, D>>::New(
        bounds, boundaryCons, gridDelta);

    T minBox[D];
    T maxBox[D];
    for (int i = 0; i < D; ++i) {
      minBox[i] = -1.0;
      maxBox[i] = 1.0;
    }

    viennals::MakeGeometry<double, D>(
        boxDomain,
        viennals::SmartPointer<viennals::Box<double, D>>::New(minBox, maxBox))
        .apply();

    viennals::BooleanOperation<double, D> boolOp(
        substrate, boxDomain,
        viennals::BooleanOperationEnum::RELATIVE_COMPLEMENT);
    boolOp.apply();

    std::cout << "Converting Box Cavity Level Set to Mesh..." << std::endl;
    auto meshBoxCavity = viennals::SmartPointer<viennals::Mesh<T>>::New();
    {
      viennals::ToSurfaceMesh<T, D> toSurfaceMesh(substrate, meshBoxCavity);
      toSurfaceMesh.setSharpCorners(true);
      toSurfaceMesh.apply();
    }
    std::string filenameBoxCavity =
        "CavityBoxFinal_" + std::to_string(D) + "D.vtp";
    viennals::VTKWriter<T>(meshBoxCavity, filenameBoxCavity).apply();

    std::vector<std::vector<T>> expectedCornersCavity;
    for (int i = 0; i < (1 << D); ++i) {
      std::vector<T> corner(D);
      for (int j = 0; j < D - 1; ++j) {
        corner[j] = ((i >> j) & 1) ? 1.0 : -1.0;
      }
      corner[D - 1] = ((i >> (D - 1)) & 1) ? 0.025 : -1.0;
      expectedCornersCavity.push_back(corner);
    }
    LSTEST_ASSERT_MESH_CORNERS(meshBoxCavity, expectedCornersCavity, D);

    viennals::ToMesh<T, D>(substrate, meshBoxCavity).apply();
    filenameBoxCavity = "CavityBoxFinal_" + std::to_string(D) + "D.vtu";
    viennals::VTKWriter<T>(meshBoxCavity, filenameBoxCavity).apply();
    std::cout << "Written mesh to " << filenameBoxCavity << std::endl;
    std::cout << std::endl;
  }
}

int main() {
  omp_set_num_threads(8);
  runTest<2>();
  runTest<3>();
  return 0;
}
