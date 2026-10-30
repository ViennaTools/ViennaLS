#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <string>

#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMesh.hpp>
#include <lsReinitialize.hpp>
#include <lsToDualContouringMesh.hpp>

/**
  Dual contouring of rotated boxes in 2D and 3D. With the closest points and
  normals stored by the Euclidean FromSurfaceMesh every round trip reproduces
  the box exactly, with a closed mesh whose triangles all have area and whose
  lines all have length. The box swept along y, as directional deposition
  grows it, has faces meeting at smaller angles and in 3D a face narrower
  than a cell, and it is reproduced exactly as well. A Manhattan level set
  without stored normals is contoured with estimated normals.
  \example DualContouring.cpp
*/

namespace ls = viennals;

using Point = std::array<double, 3>;

int fail(int D, const std::string &message) {
  std::cout << D << "D: " << message << std::endl;
  return EXIT_FAILURE;
}

template <int D> struct RotatedBox {
  Point half{3.15, 2.35, 1.85}, centre{0.31, -0.27, D == 3 ? 0.43 : 0.};
  std::array<Point, 3> R;

  RotatedBox() {
    const double a = 0.52, b = 0.35;
    if constexpr (D == 2)
      R = {{{std::cos(a), -std::sin(a), 0.},
            {std::sin(a), std::cos(a), 0.},
            {0., 0., 1.}}};
    else
      R = {
          {{std::cos(a), -std::sin(a), 0.},
           {std::cos(b) * std::sin(a), std::cos(b) * std::cos(a), -std::sin(b)},
           {std::sin(b) * std::sin(a), std::sin(b) * std::cos(a),
            std::cos(b)}}};
  }

  double sdf(const Point &p) const {
    double outside = 0., inside = -1e300;
    for (int j = 0; j < D; ++j) {
      double l = 0.;
      for (int i = 0; i < 3; ++i)
        l += R[i][j] * (p[i] - centre[i]);
      const double q = std::abs(l) - half[j];
      outside += (q > 0.) ? q * q : 0.;
      inside = std::max(inside, q);
    }
    return std::sqrt(outside) + std::min(inside, 0.);
  }

  std::vector<Point> corners() const {
    std::vector<Point> result;
    for (int n = 0; n < (1 << D); ++n) {
      Point p = centre;
      for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j)
          p[i] += R[i][j] * (((n >> j) & 1) ? half[j] : -half[j]);
      result.push_back(p);
    }
    return result;
  }

  // closed mesh with the corner order and orientation of MakeGeometry::makeBox
  ls::SmartPointer<ls::Mesh<>> mesh() const {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    for (const auto &c : corners())
      mesh->insertNextNode({c[0], c[1], c[2]});
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

  // the box swept up by r, as directional deposition grows it: the faces
  // facing up moved by r, the others in place, and between them the faces
  // swept by the edges (3D) or corners (2D) where the two kinds meet
  ls::SmartPointer<ls::Mesh<>> sweptMesh(double r) const {
    auto box = mesh();
    const auto &n = box->getNodes();
    const unsigned shift = n.size();
    auto swept = ls::SmartPointer<ls::Mesh<>>::New();
    for (const auto &c : n)
      swept->insertNextNode(c);
    for (const auto &c : n)
      swept->insertNextNode({c[0], c[1] + r, c[2]});
    // whether the outward normal of an element points up, by its edges
    std::map<std::pair<unsigned, unsigned>, bool> up;
    for (const auto &e : box->template getElements<D>()) {
      bool facesUp;
      if constexpr (D == 2)
        facesUp = n[e[1]][0] > n[e[0]][0];
      else
        facesUp = (n[e[1]][2] - n[e[0]][2]) * (n[e[2]][0] - n[e[0]][0]) >
                  (n[e[1]][0] - n[e[0]][0]) * (n[e[2]][2] - n[e[0]][2]);
      const unsigned s = facesUp ? shift : 0;
      if constexpr (D == 2) {
        swept->insertNextLine({e[0] + s, e[1] + s});
        up[{e[0], e[1]}] = facesUp;
      } else {
        swept->insertNextTriangle({e[0] + s, e[1] + s, e[2] + s});
        for (int k = 0; k < 3; ++k)
          up[{e[k], e[(k + 1) % 3]}] = facesUp;
      }
    }
    if constexpr (D == 2) {
      // a corner between an incoming and an outgoing line
      for (const auto &[in, inUp] : up)
        for (const auto &[out, outUp] : up)
          if (in.second == out.first && inUp != outUp) {
            const unsigned v = in.second;
            if (inUp)
              swept->insertNextLine({v + shift, v});
            else
              swept->insertNextLine({v, v + shift});
          }
    } else {
      // an edge a-b of a face facing down whose other face faces up
      for (const auto &[edge, edgeUp] : up) {
        const auto other = up.find({edge.second, edge.first});
        if (edgeUp || other == up.end() || !other->second)
          continue;
        const unsigned a = edge.first, b = edge.second;
        swept->insertNextTriangle({b, a, a + shift});
        swept->insertNextTriangle({b, a + shift, b + shift});
      }
    }
    return swept;
  }
};

// how far a point lies off a closed convex mesh: the largest signed distance
// to the planes of its elements, zero exactly on the surface
template <int D>
double offConvex(ls::SmartPointer<ls::Mesh<>> mesh, const Point &p) {
  const auto &n = mesh->getNodes();
  double largest = -1e300;
  for (const auto &e : mesh->template getElements<D>()) {
    ls::Vec3D<double> normal;
    if constexpr (D == 2)
      normal = {n[e[0]][1] - n[e[1]][1], n[e[1]][0] - n[e[0]][0], 0.};
    else
      normal = ls::CrossProduct(
          ls::Vec3D<double>{n[e[1]][0] - n[e[0]][0], n[e[1]][1] - n[e[0]][1],
                            n[e[1]][2] - n[e[0]][2]},
          ls::Vec3D<double>{n[e[2]][0] - n[e[0]][0], n[e[2]][1] - n[e[0]][1],
                            n[e[2]][2] - n[e[0]][2]});
    const double length = ls::Norm(normal);
    double distance = 0.;
    for (int i = 0; i < D; ++i)
      distance += normal[i] / length * (p[i] - n[e[0]][i]);
    largest = std::max(largest, distance);
  }
  return std::abs(largest);
}

// signed area (2D) or volume (3D) enclosed by an oriented mesh
template <int D> double enclosed(ls::SmartPointer<ls::Mesh<>> mesh) {
  const auto &n = mesh->getNodes();
  double sum = 0.;
  for (const auto &e : mesh->template getElements<D>()) {
    if constexpr (D == 2)
      sum += 0.5 * (n[e[0]][0] * n[e[1]][1] - n[e[1]][0] * n[e[0]][1]);
    else
      sum +=
          (n[e[0]][0] * (n[e[1]][1] * n[e[2]][2] - n[e[1]][2] * n[e[2]][1]) -
           n[e[0]][1] * (n[e[1]][0] * n[e[2]][2] - n[e[1]][2] * n[e[2]][0]) +
           n[e[0]][2] * (n[e[1]][0] * n[e[2]][1] - n[e[1]][1] * n[e[2]][0])) /
          6.;
  }
  return sum;
}

template <int D>
ls::SmartPointer<ls::Domain<double, D>> domain(double gridDelta,
                                               ls::NormalizationEnum norm) {
  double bounds[2 * D];
  ls::BoundaryConditionEnum boundaryCons[D];
  for (int i = 0; i < D; ++i) {
    bounds[2 * i] = -8.;
    bounds[2 * i + 1] = 8.;
    boundaryCons[i] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
  }
  auto levelSet = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  levelSet->setNormalization(norm);
  levelSet->setLevelSetWidth(4);
  return levelSet;
}

template <int D> int checkRoundTrips() {
  const RotatedBox<D> box;
  const double gridDelta = 0.5;
  auto mesh = box.mesh();
  // signed: ViennaLS closes 2D boxes clockwise, the extraction must agree
  const double exact = enclosed<D>(mesh);
  for (int trip = 1; trip <= 3; ++trip) {
    auto levelSet = domain<D>(gridDelta, ls::NormalizationEnum::EUCLIDEAN);
    ls::FromSurfaceMesh<double, D> mesher(levelSet, mesh);
    mesher.setUseWindingNumber(true);
    mesher.apply();
    auto extracted = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToDualContouringMesh<double, D>(levelSet, extracted).apply();

    // nodes on the box, corners reproduced, area or volume kept
    double nodeError = 0.;
    for (const auto &n : extracted->getNodes())
      nodeError = std::max(nodeError, std::abs(box.sdf({n[0], n[1], n[2]})));
    double cornerError = 0.;
    for (const auto &c : box.corners()) {
      double best = 1e300;
      for (const auto &n : extracted->getNodes())
        best =
            std::min(best, std::hypot(n[0] - c[0], n[1] - c[1], n[2] - c[2]));
      cornerError = std::max(cornerError, best);
    }
    const double measureError = std::abs(enclosed<D>(extracted) / exact - 1.);
    // elements without length or area, where vertices of neighbouring cells
    // coincide or three vertices lie on one sharp edge, are avoided by
    // welding and by the choice of the diagonal
    unsigned flat = 0;
    const auto &n = extracted->getNodes();
    for (const auto &e : extracted->template getElements<D>()) {
      if constexpr (D == 2) {
        if (std::hypot(n[e[1]][0] - n[e[0]][0], n[e[1]][1] - n[e[0]][1]) <=
            1e-10 * gridDelta)
          ++flat;
      } else {
        const ls::Vec3D<double> ab{n[e[1]][0] - n[e[0]][0],
                                   n[e[1]][1] - n[e[0]][1],
                                   n[e[1]][2] - n[e[0]][2]},
            ac{n[e[2]][0] - n[e[0]][0], n[e[2]][1] - n[e[0]][1],
               n[e[2]][2] - n[e[0]][2]};
        const double longest = std::max(
            {ls::Norm(ab), ls::Norm(ac),
             std::hypot(n[e[2]][0] - n[e[1]][0], n[e[2]][1] - n[e[1]][1],
                        n[e[2]][2] - n[e[1]][2])});
        if (ls::Norm(ls::CrossProduct(ab, ac)) <= 1e-10 * longest * longest)
          ++flat;
      }
    }
    const bool closed = ls::Reinitialize<double, D>::isClosed(*extracted);
    std::cout << D << "D trip " << trip << ": " << extracted->getNodes().size()
              << " nodes, node error " << nodeError / gridDelta
              << ", corner error " << cornerError / gridDelta << ", relative "
              << (D == 2 ? "area" : "volume") << " error " << measureError
              << ", "
              << (D == 2 ? "lines without length " : "triangles without area ")
              << flat << (closed ? ", closed" : ", open") << std::endl;
    if (flat > 0)
      return fail(D,
                  "elements without measure in trip " + std::to_string(trip));
    if (!closed)
      return fail(D, "open mesh in trip " + std::to_string(trip));
    if (cornerError / gridDelta > 1e-9)
      return fail(D, "corner not reproduced in trip " + std::to_string(trip));
    if (nodeError / gridDelta > 1e-9 || measureError > 1e-9)
      return fail(D, "round trip " + std::to_string(trip) + " leaves the box");
    mesh = extracted;
  }

  // Manhattan values without stored normals: estimated Hermite data
  auto levelSet = domain<D>(gridDelta, ls::NormalizationEnum::MANHATTAN);
  ls::FromSurfaceMesh<double, D>(levelSet, box.mesh()).apply();
  auto extracted = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToDualContouringMesh<double, D>(levelSet, extracted).apply();
  double nodeError = 0.;
  for (const auto &n : extracted->getNodes())
    nodeError = std::max(nodeError, std::abs(box.sdf({n[0], n[1], n[2]})));
  const double measureError = std::abs(enclosed<D>(extracted) / exact - 1.);
  std::cout << D << "D estimated normals: node error " << nodeError / gridDelta
            << ", relative " << (D == 2 ? "area" : "volume") << " error "
            << measureError << std::endl;
  if (extracted->getNodes().empty() || nodeError / gridDelta > 1. ||
      measureError > 0.05)
    return fail(D, "dual contouring with estimated normals is off");
  return EXIT_SUCCESS;
}

// the swept box, whose swept faces meet the box faces at normals 22.8 to
// 67.2 degrees apart (2D 29.8 and 60.2); the threshold 0.1 of earlier
// versions merged those below 35 degrees (2D 0.06 grid spacings off, 3D
// volume 1.6e-3)
template <int D> int checkSwept() {
  const RotatedBox<D> box;
  const double gridDelta = 0.5;
  auto mesh = box.sweptMesh(2.25 * gridDelta);
  auto levelSet = domain<D>(gridDelta, ls::NormalizationEnum::EUCLIDEAN);
  ls::FromSurfaceMesh<double, D> mesher(levelSet, mesh);
  mesher.setUseWindingNumber(true);
  mesher.apply();
  auto extracted = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToDualContouringMesh<double, D>(levelSet, extracted).apply();
  double nodeError = 0.;
  for (const auto &n : extracted->getNodes())
    nodeError = std::max(nodeError, offConvex<D>(mesh, {n[0], n[1], n[2]}));
  const double measureError =
      std::abs(enclosed<D>(extracted) / enclosed<D>(mesh) - 1.);
  const bool closed = ls::Reinitialize<double, D>::isClosed(*extracted);
  std::cout << D << "D swept box: " << extracted->getNodes().size()
            << " nodes, node error " << nodeError / gridDelta << ", relative "
            << (D == 2 ? "area" : "volume") << " error " << measureError
            << (closed ? ", closed" : ", open") << std::endl;
  // in 3D one cell holds both edges of a swept face 1.3 grid spacings wide
  // and gets a vertex on each (a single vertex was 0.31 grid spacings off)
  if (!closed || nodeError / gridDelta > 1e-9 || measureError > 1e-9)
    return fail(D, "swept box off");
  return EXIT_SUCCESS;
}

int main() {
  omp_set_num_threads(4);
  if (checkRoundTrips<2>() != EXIT_SUCCESS ||
      checkRoundTrips<3>() != EXIT_SUCCESS || checkSwept<2>() != EXIT_SUCCESS ||
      checkSwept<3>() != EXIT_SUCCESS)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
