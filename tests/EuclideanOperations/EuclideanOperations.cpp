#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <lsAdvect.hpp>
#include <lsAdvectHermite.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsExtrude.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsMesh.hpp>
#include <lsPrune.hpp>
#include <lsReduce.hpp>
#include <lsReinitialize.hpp>
#include <lsRemoveStrayPoints.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToDualContouringMesh.hpp>

/**
  Operations on Euclidean level sets with Hermite data, on rotated boxes in 2D
  and 3D with exact references. Reduce and Prune keep the exact values and the
  normalization. Expand and the Boolean operations take their values from the
  dual contouring mesh of the stored planes (Reinitialize): union,
  intersection and relative complement of two boxes give surfaces on the
  exact result, with the sharp edges where the boxes cut each other, and the
  exact signed distance wherever it equals the minimum or maximum of the two.
  INVERT turns the stored normals. The stray-point removal of ViennaPS keeps
  the exact values of the remaining box. ToDiskMesh places the disks on the
  exact surface. AdvectHermite widens a reduced band, runs single steps and
  carries other point data through its rebuilds. A union into an empty level
  set keeps the closest points and normals of the other, and Extrude those of
  a Euclidean level set. Advect returns a Euclidean level set whose closest
  points match its values. GeometricAdvect with a sphere or box distribution
  grows and shrinks the box exactly where the result has sharp edges.
  \example EuclideanOperations.cpp
*/

namespace ls = viennals;

using Point = std::array<double, 3>;

int fail(int D, const std::string &message) {
  std::cout << D << "D: " << message << std::endl;
  return EXIT_FAILURE;
}

template <int D> struct RotatedBox {
  Point half, centre;
  std::array<Point, 3> R;

  RotatedBox(Point passedHalf, Point passedCentre, double degreesA,
             double degreesB)
      : half(passedHalf), centre(passedCentre) {
    // Rz(a), followed by Rx(b) in 3D
    const double a = degreesA * M_PI / 180.,
                 b = (D == 3) ? degreesB * M_PI / 180. : 0.;
    const double ca = std::cos(a), sa = std::sin(a), cb = std::cos(b),
                 sb = std::sin(b);
    R = {{{ca, -sa, 0.}, {cb * sa, cb * ca, -sb}, {sb * sa, sb * ca, cb}}};
  }

  RotatedBox(Point passedHalf, Point passedCentre, std::array<Point, 3> passedR)
      : half(passedHalf), centre(passedCentre), R(passedR) {}

  // a point given in the frame of the box
  Point toWorld(const Point &local) const {
    Point p = centre;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        p[i] += R[i][j] * local[j];
    return p;
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

  // closed mesh with the corner order and orientation of MakeGeometry::makeBox
  ls::SmartPointer<ls::Mesh<>> mesh() const {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    for (int n = 0; n < (1 << D); ++n) {
      Point p = centre;
      for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j)
          p[i] += R[i][j] * (((n >> j) & 1) ? half[j] : -half[j]);
      mesh->insertNextNode(p);
    }
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
ls::SmartPointer<ls::Domain<double, D>>
levelSet(const std::vector<RotatedBox<D>> &boxes, double gridDelta,
         double extent = 18.) {
  double bounds[2 * D];
  ls::BoundaryConditionEnum boundaryCons[D];
  for (int i = 0; i < D; ++i) {
    bounds[2 * i] = -extent * gridDelta;
    bounds[2 * i + 1] = extent * gridDelta;
    boundaryCons[i] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
  }
  auto domain = ls::SmartPointer<ls::Domain<double, D>>::New(
      bounds, boundaryCons, gridDelta);
  domain->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
  domain->setLevelSetWidth(12);
  // all boxes in one mesh
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  for (const auto &box : boxes) {
    const unsigned offset = mesh->getNodes().size();
    auto part = box.mesh();
    for (const auto &n : part->getNodes())
      mesh->insertNextNode(n);
    for (auto e : part->template getElements<D>()) {
      for (auto &v : e)
        v += offset;
      mesh->insertNextElement(e);
    }
  }
  ls::FromSurfaceMesh<double, D> mesher(domain, mesh);
  mesher.setUseWindingNumber(true);
  mesher.apply();
  return domain;
}

// largest difference of the level set values from f on the nodes where
// where(f) holds, grid spacings
template <int D, class F, class W>
double valueError(ls::SmartPointer<ls::Domain<double, D>> domain, F f,
                  W where) {
  const double gridDelta = domain->getGrid().getGridDelta();
  double error = 0.;
  for (viennahrle::ConstSparseIterator<
           typename ls::Domain<double, D>::DomainType>
           it(domain->getDomain());
       !it.isFinished(); it.next()) {
    if (!it.isDefined())
      continue;
    Point p{0., 0., 0.};
    for (int i = 0; i < D; ++i)
      p[i] = it.getStartIndices()[i] * gridDelta;
    const double exact = f(p);
    if (where(exact))
      error = std::max(error, std::abs(it.getValue() * gridDelta - exact));
  }
  return error / gridDelta;
}

// largest |f| at the vertices of the dual contouring mesh, grid spacings
template <int D, class F>
double surfaceError(ls::SmartPointer<ls::Domain<double, D>> domain, F f) {
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToDualContouringMesh<double, D>(domain, mesh).apply();
  if (mesh->getNodes().empty())
    return 1e300;
  double error = 0.;
  for (const auto &n : mesh->getNodes())
    error = std::max(error, std::abs(f(Point{n[0], n[1], n[2]})));
  return error / domain->getGrid().getGridDelta();
}

template <int D> int check(double gridDelta) {
  const double s = gridDelta;
  const RotatedBox<D> A(Point{5.15 * s, 3.85 * s, 4.35 * s},
                        Point{0.31 * s, -0.27 * s, 0.43 * s}, 30., 20.);
  // Boolean operations on boxes P and Q: a corner of Q pokes two grid
  // spacings out through a face of P, with its diagonal along the face
  // normal (P's local y in 2D, x in 3D). Corners and the points where edges
  // cross the other surface stay farther apart than a cell diagonal, so that
  // in 2D every cell holds at most one feature of the results, which are
  // exact. In 3D the edges of Q leave P's face at 35 degrees: within about
  // 1.7 grid spacings of the points where they pierce it, a cell holds P's
  // face and an edge of Q. Its dual contouring vertex lies on the result, and
  // the mesh between such vertices cuts across the edge.
  const RotatedBox<D> P =
      (D == 2) ? A
               : RotatedBox<D>(Point{6.5 * s, 8. * s, 8. * s},
                               Point{0.31 * s, -0.27 * s, 0.43 * s}, 30., 20.);
  const RotatedBox<D> Q = [&]() {
    if constexpr (D == 2) {
      const double h = 2.5 * s;
      return RotatedBox<D>(
          Point{h, h, 0.},
          P.toWorld({0., P.half[1] - (std::sqrt(2.) * h - 2. * s), 0.}), 75.,
          0.);
    } else {
      // the cube diagonal (1, 1, 1) turned onto P's local x
      const double h = 3.3 * s, a = 1. / std::sqrt(3.), b = 1. / std::sqrt(2.),
                   c = 1. / std::sqrt(6.);
      const std::array<Point, 3> turn = {
          {{a, a, a}, {b, -b, 0.}, {c, c, -2. * c}}};
      std::array<Point, 3> R{};
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          for (int k = 0; k < 3; ++k)
            R[i][j] += P.R[i][k] * turn[k][j];
      return RotatedBox<D>(
          Point{h, h, h},
          P.toWorld({P.half[0] - (std::sqrt(3.) * h - 2. * s), 0., 0.}), R);
    }
  }();
  const double extent = (D == 2) ? 18. : 22.;
  // the surfaces of the Boolean operations are exact, and their values in
  // 2D. In 3D the values of the relative complement next to the edges of Q,
  // where the mesh cuts across, stay within 0.05 grid spacings
  const double booleanSurfaceBound = 1e-9,
               booleanValueBound = (D == 2) ? 1e-9 : 0.1;
  const RotatedBox<D> C(Point{1.3 * s, 1.1 * s, 1.2 * s},
                        Point{-11.2 * s, 10.7 * s, 9.6 * s}, 8., -5.);
  const auto sdfA = [&](const Point &p) { return A.sdf(p); };
  const auto everywhere = [](double) { return true; };
  const double tolerance = 1e-9;
  auto report = [&](const std::string &name, double error,
                    double bound = 1e-9) {
    std::cout << D << "D " << name << ": " << error << " grid spacings"
              << std::endl;
    return error <= bound;
  };
  const auto euclidean = [](ls::SmartPointer<ls::Domain<double, D>> domain) {
    return domain->getNormalization() == ls::NormalizationEnum::EUCLIDEAN;
  };

  // Reduce and Prune keep the values and the normalization
  {
    auto a = levelSet<D>({A}, gridDelta);
    ls::Reduce<double, D>(a, 2).apply();
    if (!euclidean(a) || !report("Reduce", valueError<D>(a, sdfA, everywhere)))
      return fail(D, "Reduce");
    ls::Prune<double, D>(a).apply();
    if (!euclidean(a) || !report("Prune", valueError<D>(a, sdfA, everywhere)))
      return fail(D, "Prune");
    // Expand fills the band again from the planes
    ls::Expand<double, D>(a, 10).apply();
    if (!euclidean(a) || a->getLevelSetWidth() != 10 ||
        !report("Expand", valueError<D>(a, sdfA, everywhere)))
      return fail(D, "Expand");
  }

  // a union into an empty level set, as a layout is assembled from its
  // shapes, keeps the closest points and normals of the shape
  {
    auto a = levelSet<D>({A}, gridDelta);
    auto empty = ls::SmartPointer<ls::Domain<double, D>>::New(a->getGrid());
    empty->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
    empty->setLevelSetWidth(12);
    ls::BooleanOperation<double, D>(empty, a, ls::BooleanOperationEnum::UNION)
        .apply();
    if (!euclidean(empty) ||
        !report("union into an empty level set",
                valueError<D>(empty, sdfA, everywhere)) ||
        !report("union into an empty level set, surface",
                surfaceError<D>(empty, sdfA)))
      return fail(D, "union into an empty level set");
  }

  // Extrude keeps a Euclidean level set: the values of the prism, and the
  // closest point and normal of the node with the same x and y, at the
  // node's own z
  if constexpr (D == 2) {
    auto a = levelSet<D>({A}, gridDelta);
    auto prism = ls::SmartPointer<ls::Domain<double, 3>>::New();
    ls::Extrude<double>(a, prism, {-3. * gridDelta, 3. * gridDelta}, 2).apply();
    const auto sdfPrism = [&](const Point &p) {
      return A.sdf({p[0], p[1], 0.});
    };
    std::map<viennahrle::Index<2>, unsigned> ids;
    for (viennahrle::ConstSparseIterator<
             typename ls::Domain<double, 2>::DomainType>
             it(a->getDomain());
         !it.isFinished(); it.next())
      if (it.isDefined())
        ids.emplace(it.getStartIndices(), it.getPointId());
    const auto &data = a->getPointData();
    const auto &prismData = prism->getPointData();
    const auto *closest =
        data.getVectorData(ls::FromSurfaceMesh<double, 2>::closestPointsLabel);
    const auto *normals =
        data.getVectorData(ls::FromSurfaceMesh<double, 2>::surfaceNormalsLabel);
    const auto *prismClosest = prismData.getVectorData(
        ls::FromSurfaceMesh<double, 3>::closestPointsLabel, true);
    const auto *prismNormals = prismData.getVectorData(
        ls::FromSurfaceMesh<double, 3>::surfaceNormalsLabel, true);
    double mismatch = (prismClosest && prismNormals) ? 0. : 1e300;
    for (viennahrle::ConstSparseIterator<
             typename ls::Domain<double, 3>::DomainType>
             it(prism->getDomain());
         mismatch < 1e300 && !it.isFinished(); it.next()) {
      if (!it.isDefined())
        continue;
      const auto &k = it.getStartIndices();
      const auto found = ids.find(viennahrle::Index<2>(k[0], k[1]));
      if (found == ids.end()) {
        mismatch = 1e300;
        break;
      }
      const auto &c = (*closest)[found->second], &n = (*normals)[found->second];
      const auto &pc = (*prismClosest)[it.getPointId()],
                 &pn = (*prismNormals)[it.getPointId()];
      mismatch =
          std::max({mismatch, std::abs(pc[0] - c[0]), std::abs(pc[1] - c[1]),
                    std::abs(pc[2] - k[2] * gridDelta), std::abs(pn[0] - n[0]),
                    std::abs(pn[1] - n[1]), std::abs(pn[2])});
    }
    if (prism->getNormalization() != ls::NormalizationEnum::EUCLIDEAN ||
        prism->getLevelSetWidth() != a->getLevelSetWidth() ||
        !report("Extrude", valueError<3>(prism, sdfPrism, everywhere)) ||
        !report("Extrude closest points and normals", mismatch / gridDelta))
      return fail(D, "Extrude");
  }

  // Boolean operations: the surface vanishes on the exact result, and the
  // values are exact wherever the minimum or maximum is the signed distance
  struct Case {
    std::string name;
    ls::BooleanOperationEnum operation;
  };
  for (const auto &c :
       std::vector<Case>{{"union", ls::BooleanOperationEnum::UNION},
                         {"intersection", ls::BooleanOperationEnum::INTERSECT},
                         {"relative complement",
                          ls::BooleanOperationEnum::RELATIVE_COMPLEMENT}}) {
    auto a = levelSet<D>({P}, gridDelta, extent);
    auto b = levelSet<D>({Q}, gridDelta, extent);
    ls::BooleanOperation<double, D>(a, b, c.operation).apply();
    const auto exact = [&](const Point &p) {
      if (c.operation == ls::BooleanOperationEnum::UNION)
        return std::min(P.sdf(p), Q.sdf(p));
      if (c.operation == ls::BooleanOperationEnum::INTERSECT)
        return std::max(P.sdf(p), Q.sdf(p));
      return std::max(P.sdf(p), -Q.sdf(p));
    };
    // outside a union and inside an intersection the minimum and maximum
    // are the signed distance
    const auto valid = [&](double e) {
      return (c.operation == ls::BooleanOperationEnum::UNION) ? e > 0. : e < 0.;
    };
    if (!euclidean(a) ||
        !report(c.name + " surface", surfaceError<D>(a, exact),
                booleanSurfaceBound) ||
        !report(c.name + " values", valueError<D>(a, exact, valid),
                booleanValueBound))
      return fail(D, c.name);
  }

  // INVERT turns values and normals, closest points stay
  {
    auto a = levelSet<D>({A}, gridDelta);
    const auto normals = *a->getPointData().getVectorData(
        ls::FromSurfaceMesh<double, D>::surfaceNormalsLabel);
    ls::BooleanOperation<double, D>(a, ls::BooleanOperationEnum::INVERT)
        .apply();
    const auto &turned = *a->getPointData().getVectorData(
        ls::FromSurfaceMesh<double, D>::surfaceNormalsLabel);
    double normalError = 0.;
    for (std::size_t i = 0; i < normals.size(); ++i)
      for (int k = 0; k < 3; ++k)
        normalError =
            std::max(normalError, std::abs(normals[i][k] + turned[i][k]));
    const auto inverted = [&](const Point &p) { return -A.sdf(p); };
    if (!euclidean(a) ||
        !report("invert values", valueError<D>(a, inverted, everywhere)) ||
        !report("invert normals", normalError))
      return fail(D, "INVERT");
  }

  // stray-point removal as in ViennaPS: the small box goes, the large one
  // keeps its exact values
  {
    auto a = levelSet<D>({A, C}, gridDelta);
    ls::BooleanOperation<double, D>(a, ls::BooleanOperationEnum::INVERT)
        .apply();
    ls::Expand<double, D>(a, 3).apply();
    ls::RemoveStrayPoints<double, D>(a).apply();
    ls::Prune<double, D>(a).apply();
    ls::BooleanOperation<double, D>(a, ls::BooleanOperationEnum::INVERT)
        .apply();
    if (!euclidean(a) ||
        !report("stray points removed", valueError<D>(a, sdfA, everywhere)) ||
        !report("stray points surface", surfaceError<D>(a, sdfA)))
      return fail(D, "stray-point removal");
  }

  // disks on the exact surface with unit normals
  {
    auto a = levelSet<D>({A}, gridDelta);
    auto disks = ls::SmartPointer<ls::Mesh<>>::New();
    ls::ToDiskMesh<double, D>(a, disks).apply();
    double positionError = 0., normalError = 0.;
    const auto *normals = disks->getCellData().getVectorData("Normals");
    for (std::size_t i = 0; i < disks->getNodes().size(); ++i) {
      const auto &n = disks->getNodes()[i];
      positionError =
          std::max(positionError, std::abs(A.sdf({n[0], n[1], n[2]})));
      const auto &v = (*normals)[i];
      normalError = std::max(
          normalError,
          std::abs(std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) - 1.));
    }
    if (disks->getNodes().empty() ||
        !report("disk positions", positionError / gridDelta) ||
        !report("disk normals", normalError))
      return fail(D, "ToDiskMesh");
  }

  // AdvectHermite: a reduced band is widened, single steps, other point data
  // carried through the rebuilds
  {
    class Etching : public ls::VelocityField<double> {
    public:
      double getScalarVelocity(const ls::Vec3D<double> &, int,
                               const ls::Vec3D<double> &,
                               unsigned long) override {
        return -1.;
      }
    };
    auto a = levelSet<D>({A}, gridDelta);
    ls::Reduce<double, D>(a, 2).apply();
    // a label holding the x coordinate of every node
    std::vector<double> x(a->getNumberOfPoints());
    for (viennahrle::ConstSparseIterator<
             typename ls::Domain<double, D>::DomainType>
             it(a->getDomain());
         !it.isFinished(); it.next())
      if (it.isDefined())
        x[it.getPointId()] = it.getStartIndices()[0] * gridDelta;
    a->getPointData().insertNextScalarData(std::move(x), "Coordinate");
    ls::AdvectHermite<double, D> advect(a, ls::SmartPointer<Etching>::New());
    advect.setAdvectionTime(2.25 * gridDelta);
    advect.setSingleStep(true);
    unsigned calls = 0;
    double time = 0.;
    while (time < 2.25 * gridDelta - 1e-12 && calls < 100) {
      advect.setAdvectionTime(2.25 * gridDelta - time);
      advect.apply();
      if (advect.getNumberOfTimeSteps() != 1)
        return fail(D, "AdvectHermite single step");
      time += advect.getAdvectedTime();
      ++calls;
    }
    // shrunk box, exact
    const double r = 2.25 * gridDelta;
    const RotatedBox<D> shrunk(
        Point{A.half[0] - r, A.half[1] - r, A.half[2] - r}, A.centre,
        (D == 3) ? 30. : 30., 20.);
    const auto exact = [&](const Point &p) { return shrunk.sdf(p); };
    // the coordinate label: exact at the nodes of the reduced band it was
    // written on, present at every node
    const auto *coordinate =
        a->getPointData().getScalarData("Coordinate", true);
    double coordinateError = 0.;
    std::size_t kept = 0;
    if (coordinate != nullptr && coordinate->size() == a->getNumberOfPoints()) {
      for (viennahrle::ConstSparseIterator<
               typename ls::Domain<double, D>::DomainType>
               it(a->getDomain());
           !it.isFinished(); it.next()) {
        if (!it.isDefined())
          continue;
        Point p{0., 0., 0.};
        for (int i = 0; i < D; ++i)
          p[i] = it.getStartIndices()[i] * gridDelta;
        if (std::abs(A.sdf(p)) > gridDelta)
          continue;
        coordinateError = std::max(
            coordinateError, std::abs((*coordinate)[it.getPointId()] - p[0]));
        ++kept;
      }
    } else {
      coordinateError = 1e300;
    }
    std::cout << D << "D AdvectHermite from a reduced band: " << calls
              << " single steps, coordinate label on " << kept
              << " nodes of the reduced band" << std::endl;
    if (!report("AdvectHermite surface", surfaceError<D>(a, exact)) ||
        !report("AdvectHermite coordinate label",
                coordinateError / gridDelta) ||
        kept == 0)
      return fail(D, "AdvectHermite");
  }

  // Advect: advected without the closest points and normals and made
  // Euclidean again with its band width, the closest points at the distance
  // of the values, on the result with rounded corners
  {
    class Etching : public ls::VelocityField<double> {
    public:
      double getScalarVelocity(const ls::Vec3D<double> &, int,
                               const ls::Vec3D<double> &,
                               unsigned long) override {
        return -1.;
      }
    };
    auto a = levelSet<D>({A}, gridDelta);
    ls::Advect<double, D> advect(a, ls::SmartPointer<Etching>::New());
    advect.setAdvectionTime(2.25 * gridDelta);
    advect.apply();
    const double r = 2.25 * gridDelta;
    const RotatedBox<D> shrunk(
        Point{A.half[0] - r, A.half[1] - r, A.half[2] - r}, A.centre, 30., 20.);
    const auto exact = [&](const Point &p) { return shrunk.sdf(p); };
    const auto *closest = a->getPointData().getVectorData(
        ls::FromSurfaceMesh<double, D>::closestPointsLabel, true);
    double consistency =
        (closest != nullptr && closest->size() == a->getNumberOfPoints())
            ? 0.
            : 1e300;
    for (viennahrle::ConstSparseIterator<
             typename ls::Domain<double, D>::DomainType>
             it(a->getDomain());
         consistency < 1e300 && !it.isFinished(); it.next()) {
      if (!it.isDefined())
        continue;
      const auto &c = (*closest)[it.getPointId()];
      double distance = 0.;
      for (int i = 0; i < D; ++i)
        distance += std::pow(it.getStartIndices()[i] * gridDelta - c[i], 2);
      consistency =
          std::max(consistency, std::abs(std::sqrt(distance) -
                                         std::abs(it.getValue()) * gridDelta));
    }
    if (!euclidean(a) || a->getLevelSetWidth() != 12 ||
        !report("Advect, closest points at the values",
                consistency / gridDelta) ||
        !report("Advect surface", surfaceError<D>(a, exact), 0.5))
      return fail(D, "Advect");
  }

  // GeometricAdvect with a sphere or box distribution: every face moves by
  // the support of the distribution and the edges and corners sweep it,
  // exactly where the result has sharp edges and within dual contouring of
  // the tangent planes where it is rounded. A ball grows the box into the
  // box with rounded edges, whose signed distance is the box's less the
  // radius, and shrinks it into the box smaller by the radius; a box shrinks
  // every face by its support and in 2D grows the box into the convex hull of
  // the sums of their corners. The band keeps its width, and its closest
  // points lie at the distance of its values.
  {
    const double r = 2.3 * gridDelta, e = 1.7 * gridDelta;
    const std::array<double, 3> b{1.9 * gridDelta, 1.3 * gridDelta,
                                  1.6 * gridDelta};
    Point inward{0., 0., 0.};
    for (int j = 0; j < D; ++j)
      for (int i = 0; i < D; ++i)
        inward[j] += b[i] * std::abs(A.R[i][j]);
    const RotatedBox<D> shrunk(
        Point{A.half[0] - e, A.half[1] - e, A.half[2] - e}, A.centre, A.R);
    const RotatedBox<D> eroded(Point{A.half[0] - inward[0],
                                     A.half[1] - inward[1],
                                     A.half[2] - inward[2]},
                               A.centre, A.R);
    // 2D: convex hull of the sums of the corners, counterclockwise
    std::vector<std::array<double, 2>> hull;
    if constexpr (D == 2) {
      std::vector<std::array<double, 2>> sums;
      for (int m = 0; m < 4; ++m) {
        const Point corner = A.toWorld({((m & 1) ? 1. : -1.) * A.half[0],
                                        ((m & 2) ? 1. : -1.) * A.half[1], 0.});
        for (int k = 0; k < 4; ++k)
          sums.push_back({corner[0] + ((k & 1) ? b[0] : -b[0]),
                          corner[1] + ((k & 2) ? b[1] : -b[1])});
      }
      std::sort(sums.begin(), sums.end());
      auto turn = [](const auto &o, const auto &p, const auto &q) {
        return (p[0] - o[0]) * (q[1] - o[1]) - (p[1] - o[1]) * (q[0] - o[0]);
      };
      for (int pass = 0; pass < 2; ++pass) {
        const std::size_t start = hull.size();
        for (const auto &p : sums) {
          while (hull.size() >= start + 2 &&
                 turn(hull[hull.size() - 2], hull.back(), p) <= 0.)
            hull.pop_back();
          hull.push_back(p);
        }
        hull.pop_back();
        std::reverse(sums.begin(), sums.end());
      }
    }
    const auto hullDistance = [&](const Point &p) {
      double inside = -1e300, outside = 1e300;
      for (std::size_t k = 0; k < hull.size(); ++k) {
        const auto &u = hull[k], &v = hull[(k + 1) % hull.size()];
        const double ex = v[0] - u[0], ey = v[1] - u[1],
                     length = std::hypot(ex, ey);
        inside = std::max(inside,
                          (ey * (p[0] - u[0]) - ex * (p[1] - u[1])) / length);
        const double t = std::clamp(((p[0] - u[0]) * ex + (p[1] - u[1]) * ey) /
                                        (length * length),
                                    0., 1.);
        outside = std::min(
            outside, std::hypot(p[0] - u[0] - t * ex, p[1] - u[1] - t * ey));
      }
      return (inside > 0.) ? outside : inside;
    };

    struct Case {
      std::string name;
      ls::SmartPointer<ls::GeometricAdvectDistribution<double, D>> distribution;
      std::function<double(const Point &)> exact;
      double bound;
      RotatedBox<D> box;
    };
    std::vector<Case> cases = {
        {"sphere deposition",
         ls::SmartPointer<ls::SphereDistribution<double, D>>::New(r),
         // the rounded edges and corners are within dual contouring of
         // their tangent planes: 0.10 (2D) and 0.23 (3D) for this radius
         [&](const Point &p) { return A.sdf(p) - r; }, (D == 2) ? 0.15 : 0.3,
         A},
        {"sphere etching",
         ls::SmartPointer<ls::SphereDistribution<double, D>>::New(-e),
         [&](const Point &p) { return shrunk.sdf(p); }, 1e-9, A},
        {"box etching",
         ls::SmartPointer<ls::BoxDistribution<double, D>>::New(
             std::array<double, 3>{-b[0], -b[1], -b[2]}),
         [&](const Point &p) { return eroded.sdf(p); }, 1e-9, A}};
    if constexpr (D == 2)
      cases.push_back({"box deposition",
                       ls::SmartPointer<ls::BoxDistribution<double, D>>::New(b),
                       hullDistance, 1e-9, A});
    // 3D: the sum of the two boxes is the zonotope of their six half axes,
    // bounded by the planes across every two of them; the largest signed
    // distance to these planes is the signed distance on and inside it. With
    // A two faces of the sum meet at 10 degrees, which dual contouring takes
    // as one plane (singular value threshold 0.01, 11.4 degrees), so the box
    // turns by 35 and 40 degrees, where the faces meet at 24 degrees or more.
    // Two cells still hold two corners of the sum, whose one vertex lies
    // within 0.09 grid spacings, as for the exact sum converted directly
    const RotatedBox<D> turned(A.half, A.centre, 35., 40.);
    std::vector<Point> generators;
    for (int j = 0; j < 3; ++j)
      generators.push_back(Point{turned.R[0][j] * turned.half[j],
                                 turned.R[1][j] * turned.half[j],
                                 turned.R[2][j] * turned.half[j]});
    for (int i = 0; i < 3; ++i) {
      Point g{0., 0., 0.};
      g[i] = b[i];
      generators.push_back(g);
    }
    const auto zonotopeDistance = [&](const Point &p) {
      double distance = -1e300;
      for (std::size_t i = 0; i < generators.size(); ++i)
        for (std::size_t j = i + 1; j < generators.size(); ++j) {
          const Point &u = generators[i], &v = generators[j];
          Point n{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                  u[0] * v[1] - u[1] * v[0]};
          const double length =
              std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
          if (length < 1e-12)
            continue;
          for (auto &c : n)
            c /= length;
          double support = 0., along = 0.;
          for (const auto &g : generators)
            support += std::abs(g[0] * n[0] + g[1] * n[1] + g[2] * n[2]);
          for (int k = 0; k < 3; ++k)
            along += n[k] * (p[k] - turned.centre[k]);
          distance = std::max(distance, std::abs(along) - support);
        }
      return distance;
    };
    if constexpr (D == 3)
      cases.push_back({"box deposition",
                       ls::SmartPointer<ls::BoxDistribution<double, D>>::New(b),
                       zonotopeDistance, 0.1, turned});
    for (const auto &c : cases) {
      auto a = levelSet<D>({c.box}, gridDelta);
      ls::GeometricAdvect<double, D>(a, c.distribution).apply();
      const auto *closest = a->getPointData().getVectorData(
          ls::FromSurfaceMesh<double, D>::closestPointsLabel, true);
      double consistency =
          (closest != nullptr && closest->size() == a->getNumberOfPoints())
              ? 0.
              : 1e300;
      for (viennahrle::ConstSparseIterator<
               typename ls::Domain<double, D>::DomainType>
               it(a->getDomain());
           consistency < 1e300 && !it.isFinished(); it.next()) {
        if (!it.isDefined())
          continue;
        const auto &q = (*closest)[it.getPointId()];
        double distance = 0.;
        for (int i = 0; i < D; ++i)
          distance += std::pow(it.getStartIndices()[i] * gridDelta - q[i], 2);
        consistency = std::max(consistency,
                               std::abs(std::sqrt(distance) -
                                        std::abs(it.getValue()) * gridDelta));
      }
      if (!euclidean(a) || a->getLevelSetWidth() != 12 ||
          !report("GeometricAdvect " + c.name +
                      ", closest points at the values",
                  consistency / gridDelta) ||
          !report("GeometricAdvect " + c.name, surfaceError<D>(a, c.exact),
                  c.bound))
        return fail(D, "GeometricAdvect " + c.name);
    }
  }
  return EXIT_SUCCESS;
}

int main() {
  omp_set_num_threads(4);
  if (check<2>(0.5) != EXIT_SUCCESS || check<3>(1.) != EXIT_SUCCESS)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
