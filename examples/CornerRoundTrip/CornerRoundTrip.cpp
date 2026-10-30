#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMesh.hpp>
#include <lsToDualContouringMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>

/**
  A rotated square (2D) and a rotated cube (3D) are converted from a surface
  mesh to a level set and back several times: marching cubes on Manhattan
  and on Euclidean level set values, and dual contouring on Euclidean values
  with the closest points and normals stored by FromSurfaceMesh. Every round
  trip writes the extracted surface to a VTP file and prints how far the true
  corners and edges lie from it, together with the change of area or
  volume.
  All lengths are in units of the grid spacing.
  \example CornerRoundTrip.cpp
*/

namespace ls = viennals;
using namespace viennacore; // vector arithmetic on Vec3D
using T = double;
using Vec = ls::Vec3D<T>;

constexpr T gridDelta = 1.0;
constexpr int numberOfTrips = 5;

T dot(const Vec &a, const Vec &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vec closestOnSegment(const Vec &p, const Vec &a, const Vec &b) {
  const Vec ab = b - a;
  const T len2 = dot(ab, ab);
  if (len2 <= 0.)
    return a;
  const T t = std::clamp(dot(p - a, ab) / len2, T(0), T(1));
  return a + ab * t;
}

// closest point on triangle abc (Ericson, Real-Time Collision Detection 5.1.5)
Vec closestOnTriangle(const Vec &p, const Vec &a, const Vec &b, const Vec &c) {
  const Vec ab = b - a, ac = c - a, ap = p - a;
  const T d1 = dot(ab, ap), d2 = dot(ac, ap);
  if (d1 <= 0 && d2 <= 0)
    return a;
  const Vec bp = p - b;
  const T d3 = dot(ab, bp), d4 = dot(ac, bp);
  if (d3 >= 0 && d4 <= d3)
    return b;
  const T vc = d1 * d4 - d3 * d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0)
    return a + ab * (d1 / (d1 - d3));
  const Vec cp = p - c;
  const T d5 = dot(ab, cp), d6 = dot(ac, cp);
  if (d6 >= 0 && d5 <= d6)
    return c;
  const T vb = d5 * d2 - d1 * d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0)
    return a + ac * (d2 / (d2 - d6));
  const T va = d3 * d6 - d5 * d4;
  if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0)
    return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
  const T denom = va + vb + vc;
  if (std::abs(denom) < 1e-300)
    return a;
  return a + ab * (vb / denom) + ac * (vc / denom);
}

template <int D>
T distanceToMesh(const Vec &p, ls::SmartPointer<ls::Mesh<T>> mesh) {
  const auto &nodes = mesh->getNodes();
  T best = std::numeric_limits<T>::max();
  for (const auto &e : mesh->template getElements<D>()) {
    Vec q;
    if constexpr (D == 2)
      q = closestOnSegment(p, nodes[e[0]], nodes[e[1]]);
    else
      q = closestOnTriangle(p, nodes[e[0]], nodes[e[1]], nodes[e[2]]);
    const Vec d = p - q;
    best = std::min(best, dot(d, d));
  }
  return std::sqrt(best);
}

// signed area (2D) or volume (3D) enclosed by an oriented surface mesh
template <int D> T enclosedMeasure(ls::SmartPointer<ls::Mesh<T>> mesh) {
  const auto &nodes = mesh->getNodes();
  T sum = 0.;
  for (const auto &e : mesh->template getElements<D>()) {
    if constexpr (D == 2) {
      const Vec &a = nodes[e[0]], &b = nodes[e[1]];
      sum += 0.5 * (a[0] * b[1] - b[0] * a[1]);
    } else {
      const Vec &a = nodes[e[0]], &b = nodes[e[1]], &c = nodes[e[2]];
      sum += dot(a, CrossProduct(b, c)) / 6.;
    }
  }
  return sum;
}

// box with half sizes `half`, rotated by R (rows) about its centre
template <int D> struct RotatedBox {
  Vec half, centre;
  std::array<Vec, 3> R;

  Vec toWorld(const Vec &local) const {
    return Vec{dot(R[0], local), dot(R[1], local), dot(R[2], local)} + centre;
  }

  T sdf(const Vec &p) const {
    const Vec d = p - centre;
    T outside = 0., inside = std::numeric_limits<T>::lowest();
    for (int j = 0; j < D; ++j) {
      const T local = R[0][j] * d[0] + R[1][j] * d[1] + R[2][j] * d[2];
      const T q = std::abs(local) - half[j];
      outside += std::max(q, T(0)) * std::max(q, T(0));
      inside = std::max(inside, q);
    }
    return std::sqrt(outside) + std::min(inside, T(0));
  }

  // corners in the order of MakeGeometry::makeBox: bit 0 is x, 1 is y, 2 is z
  std::vector<Vec> corners() const {
    std::vector<Vec> result;
    for (int n = 0; n < (1 << D); ++n) {
      Vec local{0., 0., 0.};
      for (int j = 0; j < D; ++j)
        local[j] = ((n >> j) & 1) ? half[j] : -half[j];
      result.push_back(toWorld(local));
    }
    return result;
  }

  std::vector<std::pair<Vec, Vec>> edges() const {
    std::vector<std::pair<Vec, Vec>> result;
    const auto c = corners();
    for (int n = 0; n < (1 << D); ++n)
      for (int j = 0; j < D; ++j)
        if (!((n >> j) & 1))
          result.emplace_back(c[n], c[n | (1 << j)]);
    return result;
  }

  // closed surface mesh with the orientation used by MakeGeometry::makeBox
  ls::SmartPointer<ls::Mesh<T>> mesh() const {
    auto mesh = ls::Mesh<T>::New();
    for (const auto &c : corners())
      mesh->insertNextNode(c);
    if constexpr (D == 2) {
      const std::array<unsigned, 2> lines[4] = {{0, 2}, {2, 3}, {3, 1}, {1, 0}};
      for (const auto &l : lines)
        mesh->insertNextLine(l);
    } else {
      const std::array<unsigned, 3> triangles[12] = {
          {0, 3, 1}, {0, 2, 3}, {0, 1, 5}, {0, 5, 4}, {0, 4, 2}, {4, 6, 2},
          {7, 6, 4}, {7, 4, 5}, {7, 2, 6}, {7, 3, 2}, {1, 3, 5}, {3, 7, 5}};
      for (const auto &t : triangles)
        mesh->insertNextTriangle(t);
    }
    return mesh;
  }
};

template <int D>
void roundTrips(const RotatedBox<D> &box, ls::NormalizationEnum normalization,
                bool dualContouring) {
  const bool euclidean = normalization == ls::NormalizationEnum::EUCLIDEAN;
  const std::string prefix = std::to_string(D) + "D_" +
                             (euclidean ? "euclidean_" : "manhattan_") +
                             (dualContouring ? "dc" : "mc");

  // domain around the box with a margin of four cells
  T bounds[2 * D];
  ls::BoundaryConditionEnum boundaryCons[D];
  for (int j = 0; j < D; ++j) {
    bounds[2 * j] = std::numeric_limits<T>::max();
    bounds[2 * j + 1] = std::numeric_limits<T>::lowest();
    for (const auto &c : box.corners()) {
      bounds[2 * j] = std::min(bounds[2 * j], c[j] - 4 * gridDelta);
      bounds[2 * j + 1] = std::max(bounds[2 * j + 1], c[j] + 4 * gridDelta);
    }
    boundaryCons[j] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
  }

  // sample points along the true edges, endpoints excluded
  std::vector<Vec> edgePoints;
  if constexpr (D == 3)
    for (const auto &[a, b] : box.edges())
      for (int i = 1; i < 50; ++i)
        edgePoints.push_back(a + (b - a) * (i / 50.));

  auto mesh = box.mesh();
  const T exactMeasure = enclosedMeasure<D>(mesh);
  ls::VTKWriter<T>(mesh, std::to_string(D) + "D_exact.vtp").apply();

  std::cout << "\n"
            << D << "D, " << (euclidean ? "Euclidean" : "Manhattan")
            << " values, "
            << (dualContouring ? "dual contouring" : "marching cubes")
            << "\n trip  corner max  corner mean"
            << (D == 3 ? "  edge max  edge mean" : "") << "  node deviation  "
            << (D == 2 ? "area" : "volume") << " change\n";

  for (int trip = 1; trip <= numberOfTrips; ++trip) {
    auto levelSet = ls::SmartPointer<ls::Domain<T, D>>::New(
        bounds, boundaryCons, gridDelta);
    levelSet->setNormalization(normalization);
    ls::FromSurfaceMesh<T, D> mesher(levelSet, mesh);
    mesher.setUseWindingNumber(dualContouring); // closed meshes, robust sign
    mesher.apply();

    if (trip == 1) {
      auto points = ls::Mesh<T>::New();
      ls::ToMesh<T, D>(levelSet, points).apply();
      ls::VTKWriter<T>(points, prefix + "_levelset.vtp").apply();
    }

    auto extracted = ls::Mesh<T>::New();
    if (dualContouring)
      ls::ToDualContouringMesh<T, D>(levelSet, extracted).apply();
    else
      ls::ToSurfaceMesh<T, D>(levelSet, extracted).apply();
    ls::VTKWriter<T>(extracted,
                     prefix + "_trip" + std::to_string(trip) + ".vtp")
        .apply();

    T cornerMax = 0., cornerSum = 0.;
    for (const auto &c : box.corners()) {
      const T d = distanceToMesh<D>(c, extracted);
      cornerMax = std::max(cornerMax, d);
      cornerSum += d;
    }
    T edgeMax = 0., edgeSum = 0.;
    for (const auto &p : edgePoints) {
      const T d = distanceToMesh<D>(p, extracted);
      edgeMax = std::max(edgeMax, d);
      edgeSum += d;
    }
    T nodeDeviation = 0.;
    for (const auto &n : extracted->getNodes())
      nodeDeviation = std::max(nodeDeviation, std::abs(box.sdf(n)));
    const T measure = enclosedMeasure<D>(extracted);

    std::cout << std::fixed << std::setprecision(3) << std::setw(5) << trip
              << std::setw(12) << cornerMax / gridDelta << std::setw(13)
              << cornerSum / (1 << D) / gridDelta;
    if constexpr (D == 3)
      std::cout << std::setw(10) << edgeMax / gridDelta << std::setw(11)
                << edgeSum / edgePoints.size() / gridDelta;
    std::cout << std::setw(16) << nodeDeviation / gridDelta << std::setw(10)
              << std::setprecision(2) << 100. * (measure / exactMeasure - 1.)
              << " %" << std::endl;

    mesh = extracted;
  }
}

int main() {
  omp_set_num_threads(4);

  const T a = 30. * M_PI / 180.; // rotation about z
  const T b = 20. * M_PI / 180.; // rotation about x (3D only)

  RotatedBox<2> square;
  square.half = Vec{5.15, 3.85, 0.};
  square.centre = Vec{0.31, -0.27, 0.};
  square.R = {Vec{std::cos(a), -std::sin(a), 0.},
              Vec{std::sin(a), std::cos(a), 0.}, Vec{0., 0., 1.}};

  RotatedBox<3> cube;
  cube.half = Vec{5.15, 4.35, 3.55};
  cube.centre = Vec{0.31, -0.27, 0.43};
  // R = Rx(b) * Rz(a)
  cube.R = {
      Vec{std::cos(a), -std::sin(a), 0.},
      Vec{std::cos(b) * std::sin(a), std::cos(b) * std::cos(a), -std::sin(b)},
      Vec{std::sin(b) * std::sin(a), std::sin(b) * std::cos(a), std::cos(b)}};

  for (auto [normalization, dualContouring] :
       {std::pair{ls::NormalizationEnum::MANHATTAN, false},
        std::pair{ls::NormalizationEnum::EUCLIDEAN, false},
        std::pair{ls::NormalizationEnum::EUCLIDEAN, true}}) {
    roundTrips<2>(square, normalization, dualContouring);
    roundTrips<3>(cube, normalization, dualContouring);
  }
  return EXIT_SUCCESS;
}
