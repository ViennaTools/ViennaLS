#include <iostream>

#include <lsDomain.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

#include <hrleSparseIterator.hpp>

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

  std::string line;
  unsigned lastIdx = -1;
  unsigned firstIdx = -1;
  std::array<double, 3> firstPoint = {0., 0., 0.};
  std::array<double, 3> lastPoint = {0., 0., 0.};

  // Track bounding box
  double minx = 1e100, maxx = -1e100;
  double miny = 1e100, maxy = -1e100;
  double minz = 1e100, maxz = -1e100;

  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string value;
    std::vector<double> coords;

    while (std::getline(ss, value, ',')) {
      coords.push_back(std::stod(value));
    }

    double x = 0.0, y = 0.0, z = 0.0;
    if (coords.size() >= 2) {
      x = coords[0];
      y = coords[1];
      if (!is2D && coords.size() >= 3)
        z = coords[2];
    } else {
      std::cerr << "Invalid line in CSV: " << line << std::endl;
      continue;
    }

    x = x * scaleFactor + shiftX;
    y = y * scaleFactor + shiftY;
    z = z * scaleFactor + shiftZ;

    // Update bounds
    minx = std::min(minx, x); maxx = std::max(maxx, x);
    miny = std::min(miny, y); maxy = std::max(maxy, y);
    minz = std::min(minz, z); maxz = std::max(maxz, z);

    // Insert node
    unsigned idx = mesh->insertNextNode({x, y, z});
    if (firstIdx == -1) {
      firstIdx = idx;
      firstPoint = {x, y, z};
    } else {
      mesh->insertNextLine({idx, lastIdx});
    }

    lastIdx = idx;
    lastPoint = {x, y, z};
  }
  file.close();

  // Close the polygon if last point != first point
  if ((lastIdx != -1) && (firstPoint != lastPoint)) {
    mesh->insertNextLine({firstIdx, lastIdx});
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
  double gridDelta = 1.0;

  double bounds[2 * D] = {-xExtent / 2., xExtent / 2., -yExtent / 2., yExtent / 2.};
  ls::Domain<double, D>::BoundaryType boundaryCons[D];
  boundaryCons[0] = ls::Domain<double, D>::BoundaryType::PERIODIC_BOUNDARY;
  boundaryCons[1] = ls::Domain<double, D>::BoundaryType::PERIODIC_BOUNDARY;

  auto substrate = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  // copy the structure to add the pattern on top
  auto pattern = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);

  // Create pattern from CSV polygon
  {
    std::cout << "Reading polygon from CSV..." << std::endl;

    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    readPolygonCSV(csvFilename, mesh);

    // Convert to domain
    ls::FromSurfaceMesh<double, D>(pattern, mesh).apply();
    auto writemesh = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToSurfaceMesh<double, D>(pattern, writemesh).apply();
    ls::VTKWriter<double>(writemesh, "writtenMesh.vtp").apply();

    std::cout << "Done reading polygon from CSV." << std::endl;
  }

  std::cout << "Iterating over level set grid points with negative signed distance..." << std::endl;

  using LevelSetType = ls::Domain<double, D>;
  auto &levelSet = pattern->getDomain();

  hrleSparseIterator<typename LevelSetType::DomainType> it(levelSet);
  for (; !it.isFinished(); ++it) {
    const auto &coord = it.getStartIndices();
    const double value = it.getValue();
    if (value <= 0.0) {
      // Gaussians applied here
      std::cout << "Point (" << coord[0] << ", " << coord[1];
      if constexpr (D == 3) std::cout << ", " << coord[2];
      std::cout << ") is exposed; SDF = " << value << std::endl;
    }
  }

  return 0;
}
