#include <iostream>

#include <lsDomain.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsToMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsBooleanOperation.hpp>

#include <lsExtrude.hpp>

#include <hrleSparseIterator.hpp>
#include <hrleDenseIterator.hpp>

/**
  Example showing how to create a 2D level set domain from
  a polygon. The polygon is read from a CSV file and converted
  to a level set domain.
  The intended use is to identify points which serve as exposure
  locations for subsequent Gaussian convolution.
  \example Lithography.cpp
*/

namespace ls = viennals;

void readPolygonCSV(const std::string &filename,
                    ls::SmartPointer<ls::Mesh<>> &mesh,
                    double scaleFactor = 1.0,
                    double shiftX = 0.0,
                    double shiftY = 0.0,
                    double shiftZ = 0.0,
                    bool is2D = true) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  std::vector<std::array<double, 3>> points;
  std::string line;

  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string value;
    std::vector<double> coords;

    while (std::getline(ss, value, ',')) {
      coords.push_back(std::stod(value));
    }

    double x = 0., y = 0., z = 0.;
    if (coords.size() >= 2) {
      x = coords[0];
      y = coords[1];
      if (!is2D && coords.size() >= 3)
        z = coords[2];
    } else {
      std::cerr << "Invalid line in CSV: " << line << std::endl;
      continue;
    }

    points.push_back({
      x * scaleFactor + shiftX,
      y * scaleFactor + shiftY,
      z * scaleFactor + shiftZ
    });
  }
  file.close();

  if (points.size() < 2) return;

  // --- Determine winding order using signed area ---
  double signedArea = 0.0;
  for (size_t i = 0; i < points.size(); ++i) {
    const auto &p1 = points[i];
    const auto &p2 = points[(i + 1) % points.size()];
    signedArea += (p1[0] * p2[1] - p2[0] * p1[1]);
  }
  bool isCCW = (signedArea > 0);  // counter-clockwise if true

  // --- Insert nodes and lines with correct order ---
  std::vector<unsigned> indices;
  for (const auto &p : points) {
    indices.push_back(mesh->insertNextNode(p));
  }

  for (size_t i = 1; i < indices.size(); ++i) {
    if (isCCW)
      mesh->insertNextLine({indices[i], indices[i - 1]});
    else
      mesh->insertNextLine({indices[i - 1], indices[i]});
  }

  // Close the polygon if not already closed
  if (points.front() != points.back()) {
    if (isCCW)
      mesh->insertNextLine({indices.front(), indices.back()});
    else
      mesh->insertNextLine({indices.back(), indices.front()});
  }
}

int main(int argc, char* argv[]) {
  ls::Logger::getInstance().setLogLevel(ls::LogLevel::DEBUG);

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << "<polygonFile.csv>" << std::endl;
    return 1;
  }

  std::string csvFilename = argv[1];
  constexpr int D = 2;

  // scale in micrometers
  double xExtent = 50;
  double yExtent = 50;
  double gridDelta = 0.199;

  double bounds[2 * D] = {-xExtent / 2., xExtent / 2., -yExtent / 2., yExtent / 2.};
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = ls::Domain<double, D>::BoundaryType::PERIODIC_BOUNDARY;
  boundaryCons[1] = ls::Domain<double, D>::BoundaryType::PERIODIC_BOUNDARY;

  // copy the structure to add the pattern on top
  auto pattern = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  // Create pattern from CSV polygon
  {
    std::cout << "== Reading polygon from CSV... ";

    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    readPolygonCSV(csvFilename, mesh);

    // Convert to domain
    ls::FromSurfaceMesh<double, D>(pattern, mesh).apply();
    auto writemesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(pattern, writemesh).apply();
    ls::VTKWriter<double>(writemesh, "writtenMesh.vtp").apply();

    std::cout << "done!" << std::endl;
  }

  // === Iterate over level set grid ===
  using LevelSetType = ls::Domain<double, D>;
  auto &levelSet = pattern->getDomain();

  // Sparse iterator
  hrleSparseIterator<typename LevelSetType::DomainType> its(levelSet);
  for (; !its.isFinished(); ++its) {
    const auto &coord = its.getStartIndices();
    const double value = its.getValue();
    if (value <= 0.0) {
      // Do whatever you want with the points
      // Take care of runs (see HRLE)
      // std::cout << "Point (" << coord[0] << ", " << coord[1] << ") is exposed; SDF = " << value << std::endl;
    }
  }

  // Dense iterator
  using LevelSetType = ls::Domain<double, D>;
  hrleDenseIterator<typename LevelSetType::DomainType> itd(levelSet);
  for (; !itd.isFinished(); ++itd) {
    auto idx = itd.getIndices();
    double value = itd.getValue();
    if (value <= 0.) {
      // Do whatever you want with the points
      // std::cout << "Point (" << idx[0] << ", " << idx[1] << ") is exposed; SDF = " << value << std::endl;
    }
  }

  // === Extrude 2D pattern into 3D ===
  std::cout << "== Extruding 2D into 3D... ";

  constexpr int D3 = 3;
  using Domain3D = ls::Domain<double, D3>;
  auto pattern3D = ls::SmartPointer<Domain3D>::New();

  // Define extrusion extent in Z (0 to 5)
  std::array<double, 2> extrudeExtent = {0. - gridDelta, 5. + gridDelta}; // add buffer or +/- gridDelta

  // Boundary conditions in 3D (x, y periodic, z reflective)
  std::array<ls::BoundaryConditionEnum, 3> boundaryConds = {
      ls::BoundaryConditionEnum::PERIODIC_BOUNDARY,
      ls::BoundaryConditionEnum::PERIODIC_BOUNDARY,
      ls::BoundaryConditionEnum::INFINITE_BOUNDARY};

  // Run extrusion using the available mesh-based implementation
  // ls::Extrude<double>(pattern, pattern3D, extrudeExtent, 2, boundaryConds).apply();
  ls::Extrude<double>(pattern, pattern3D, extrudeExtent, 2, boundaryConds).apply();

  // === Export extruded mesh to SDF VTK grid for visualization ===
  auto sdf3D_preCap = ls::SmartPointer<ls::Mesh<double>>::New();
  ls::ToMesh<double, D3>(pattern3D, sdf3D_preCap).apply();
  ls::VTKWriter<double>(sdf3D_preCap, "extrudedSDF_preCap.vtp").apply();

  // === Add capping layers (bottom and top planes) ===
  std::cout << "adding capping planes... ";

  // Create bottom substrate (z >= 0)
  auto bottomLS = ls::SmartPointer<Domain3D>::New(pattern3D->getGrid());
  double originLow[3] = {0., 0., 0.};
  double normalLow[3] = {0., 0., -1.}; // z >= 0
  auto bottomPlane = ls::SmartPointer<ls::Plane<double, D3>>::New(originLow, normalLow);
  ls::MakeGeometry<double, 3>(bottomLS, bottomPlane).apply();

  // Create top cap (z <= 5 µm)
  auto topLS = ls::SmartPointer<Domain3D>::New(pattern3D->getGrid());
  double originHigh[3] = {0., 0., 5.0}; // Adjust to match extrusion
  double normalHigh[3] = {0., 0., 1.}; // z <= 5
  auto topPlane = ls::SmartPointer<ls::Plane<double, D3>>::New(originHigh, normalHigh);
  ls::MakeGeometry<double, D3>(topLS, topPlane).apply();

  // Intersect with bottom
  ls::BooleanOperation<double, D3>(pattern3D, bottomLS, ls::BooleanOperationEnum::INTERSECT).apply();
  // Intersect with top
  ls::BooleanOperation<double, D3>(pattern3D, topLS, ls::BooleanOperationEnum::INTERSECT).apply();
  std::cout << "done!" << std::endl;

  // === Export extruded mesh to SDF VTK grid for visualization ===
  auto sdf3D = ls::SmartPointer<ls::Mesh<double>>::New();
  ls::ToMesh<double, D3>(pattern3D, sdf3D).apply();
  ls::VTKWriter<double>(sdf3D, "extrudedSDF.vtp").apply();

  // === Export Surface Mesh (safe if a surface exists)
  try {
    auto mesh3D = ls::SmartPointer<ls::Mesh<double>>::New();
    ls::ToSurfaceMesh<double, D3>(pattern3D, mesh3D).apply();
    if (!mesh3D->getNodes().empty()) {
      ls::VTKWriter<double>(mesh3D, "extrudedMesh.vtp").apply();
    } else {
      std::cerr << "Warning: No surface found to export!" << std::endl;
    }
  } catch (...) {
    std::cerr << "Exception: Failed to extract surface mesh from extruded level set." << std::endl;
  }

  return 0;
}
