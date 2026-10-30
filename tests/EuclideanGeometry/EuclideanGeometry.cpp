#include <array>
#include <cmath>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

#include <lsCheck.hpp>
#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsFromVolumeMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsReader.hpp>
#include <lsToDualContouringMesh.hpp>
#include <lsToMesh.hpp>
#include <lsWriter.hpp>

/**
  Checks the normalization tag of a level set (its default, copies, the
  serialization round trip, streams written before format version 1 and the
  Writer/Reader file round trip) and the Euclidean construction against
  exact signed distances: a rotated box from a surface mesh with both sign
  methods, spheres and tilted planes from MakeGeometry, a two-material
  volume mesh, and a box with its faces through grid nodes at the default
  width, whose band also holds the nodes one grid spacing away.
  \example EuclideanGeometry.cpp
*/

namespace ls = viennals;

template <int D> using LevelSetType = ls::Domain<double, D>;

int fail(int D, const std::string &message) {
  std::cout << D << "D: " << message << std::endl;
  return EXIT_FAILURE;
}

template <int D>
std::string serialized(ls::SmartPointer<LevelSetType<D>> levelSet) {
  std::stringstream stream;
  levelSet->serialize(stream);
  return stream.str();
}

template <int D>
ls::SmartPointer<LevelSetType<D>> deserialized(const std::string &bytes) {
  auto levelSet = ls::SmartPointer<LevelSetType<D>>::New();
  std::stringstream stream(bytes);
  levelSet->deserialize(stream);
  return levelSet;
}

template <int D> int checkNormalization() {
  double bounds[2 * D];
  ls::BoundaryConditionEnum boundaryCons[D];
  for (int i = 0; i < D; ++i) {
    bounds[2 * i] = -10.;
    bounds[2 * i + 1] = 10.;
    boundaryCons[i] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
  }
  auto levelSet =
      ls::SmartPointer<LevelSetType<D>>::New(bounds, boundaryCons, 0.5);
  ls::VectorType<double, D> centre{};
  for (int i = 0; i < D; ++i)
    centre[i] = 0.3 * (i + 1);
  ls::MakeGeometry<double, D>(
      levelSet, ls::SmartPointer<ls::Sphere<double, D>>::New(centre, 4.1))
      .apply();

  if (levelSet->getNormalization() != ls::NormalizationEnum::MANHATTAN)
    return fail(D, "default normalization is not MANHATTAN");

  levelSet->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
  if (levelSet->getNormalization() != ls::NormalizationEnum::EUCLIDEAN)
    return fail(D, "setNormalization had no effect");

  auto copy = ls::SmartPointer<LevelSetType<D>>::New(levelSet);
  if (copy->getNormalization() != ls::NormalizationEnum::EUCLIDEAN)
    return fail(D, "copy constructor lost the normalization");
  auto target = ls::SmartPointer<LevelSetType<D>>::New();
  target->deepCopy(levelSet);
  if (target->getNormalization() != ls::NormalizationEnum::EUCLIDEAN)
    return fail(D, "deepCopy lost the normalization");

  // the tag is written after the point data, so check both stream layouts
  for (bool withPointData : {false, true}) {
    auto source = ls::SmartPointer<LevelSetType<D>>::New(levelSet);
    if (withPointData) {
      typename ls::PointData<double>::ScalarDataType scalars(
          source->getNumberOfPoints(), 1.);
      source->getPointData().insertNextScalarData(scalars, "scalars");
    }
    const std::string label =
        withPointData ? " (with point data)" : " (without point data)";
    const std::string bytes = serialized<D>(source);

    auto roundTrip = deserialized<D>(bytes);
    if (roundTrip->getNormalization() != ls::NormalizationEnum::EUCLIDEAN)
      return fail(D, "round trip lost the normalization" + label);
    if (serialized<D>(roundTrip) != bytes)
      return fail(D, "round trip changed the level set" + label);

    // a stream of format version 0 has no tag and must read as MANHATTAN
    std::string legacy = bytes;
    legacy[8] = 0;
    legacy.pop_back();
    auto fromLegacy = deserialized<D>(legacy);
    if (fromLegacy->getNormalization() != ls::NormalizationEnum::MANHATTAN)
      return fail(D, "version 0 stream did not read as MANHATTAN" + label);
    std::string expected = bytes;
    expected.back() = static_cast<char>(ls::NormalizationEnum::MANHATTAN);
    if (serialized<D>(fromLegacy) != expected)
      return fail(D, "version 0 stream changed the level set" + label);
  }

  const std::string fileName =
      "EuclideanGeometry" + std::to_string(D) + "D.lvst";
  ls::Writer<double, D>(levelSet, fileName).apply();
  auto fromFile = ls::SmartPointer<LevelSetType<D>>::New();
  ls::Reader<double, D>(fromFile, fileName).apply();
  if (fromFile->getNormalization() != ls::NormalizationEnum::EUCLIDEAN)
    return fail(D, "file round trip lost the normalization");
  if (serialized<D>(fromFile) != serialized<D>(levelSet))
    return fail(D, "file round trip changed the level set");

  return EXIT_SUCCESS;
}

// Rotated box with half sizes `half` about `centre`, rows of R are the rotation
template <int D> struct RotatedBox {
  std::array<double, 3> half, centre;
  std::array<std::array<double, 3>, 3> R;

  std::array<double, 3> local(const std::array<double, 3> &p) const {
    std::array<double, 3> l{0., 0., 0.};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        l[j] += R[i][j] * (p[i] - centre[i]);
    return l;
  }

  double sdf(const std::array<double, 3> &p) const {
    const auto l = local(p);
    double outside = 0., inside = -1e300;
    for (int j = 0; j < D; ++j) {
      const double q = std::abs(l[j]) - half[j];
      outside += (q > 0.) ? q * q : 0.;
      inside = std::max(inside, q);
    }
    return std::sqrt(outside) + std::min(inside, 0.);
  }

  // closed surface mesh with the corner order and orientation of makeBox
  ls::SmartPointer<ls::Mesh<>> mesh() const {
    auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
    for (int n = 0; n < (1 << D); ++n) {
      std::array<double, 3> l{0., 0., 0.};
      for (int j = 0; j < D; ++j)
        l[j] = ((n >> j) & 1) ? half[j] : -half[j];
      ls::Vec3D<double> p{centre[0], centre[1], D == 3 ? centre[2] : 0.};
      for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j)
          p[i] += R[i][j] * l[j];
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

using Point = std::array<double, 3>;

// Euclidean domain on [-8, 8] per direction with reflective boundaries; with
// an infinite last direction the bounds there are one grid spacing, as in
// ViennaPS domains
template <int D>
ls::SmartPointer<LevelSetType<D>> euclideanDomain(double gridDelta, int width,
                                                  bool infiniteLast = false) {
  double bounds[2 * D];
  ls::BoundaryConditionEnum boundaryCons[D];
  for (int i = 0; i < D; ++i) {
    bounds[2 * i] = -8.;
    bounds[2 * i + 1] = 8.;
    boundaryCons[i] = ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY;
  }
  if (infiniteLast) {
    bounds[2 * D - 2] = -gridDelta;
    bounds[2 * D - 1] = gridDelta;
    boundaryCons[D - 1] = ls::BoundaryConditionEnum::INFINITE_BOUNDARY;
  }
  auto levelSet =
      ls::SmartPointer<LevelSetType<D>>::New(bounds, boundaryCons, gridDelta);
  levelSet->setNormalization(ls::NormalizationEnum::EUCLIDEAN);
  levelSet->setLevelSetWidth(width);
  return levelSet;
}

// Compares a Euclidean level set with the exact signed distance `sdf`: every
// stored value inside `inRegion` equals the exact distance in grid spacings,
// no value exceeds the band, and every grid node between the indices lo and
// hi that lies inside the band and inside `inRegion` is stored.
template <int D, class SdfType, class RegionType>
int compareWithExact(const std::string &name,
                     ls::SmartPointer<LevelSetType<D>> levelSet, SdfType sdf,
                     RegionType inRegion, std::array<long, 3> lo,
                     std::array<long, 3> hi) {
  const double gridDelta = levelSet->getGrid().getGridDelta();
  const double halfWidth = 0.5 * levelSet->getLevelSetWidth();
  if (levelSet->getNormalization() != ls::NormalizationEnum::EUCLIDEAN)
    return fail(D, name + ": normalization is not EUCLIDEAN");
  ls::Check<double, D> check(levelSet);
  check.apply();
  if (!check.isValid())
    return fail(D, name + ": level set is not valid: " + check.what());

  auto points = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, D>(levelSet, points).apply();
  const auto &values = *points->getCellData().getScalarData("LSValues");
  std::set<std::array<long, 3>> stored;
  double maxError = 0.;
  for (std::size_t n = 0; n < values.size(); ++n) {
    const auto &x = points->getNodes()[n];
    const Point p{x[0], x[1], x[2]};
    if (std::abs(values[n]) > halfWidth + 1e-12)
      return fail(D, name + ": value outside the band " +
                         std::to_string(values[n]));
    stored.insert({std::lround(x[0] / gridDelta), std::lround(x[1] / gridDelta),
                   D == 3 ? std::lround(x[2] / gridDelta) : 0L});
    if (inRegion(p))
      maxError = std::max(maxError, std::abs(values[n] - sdf(p) / gridDelta));
  }
  if (maxError > 1e-9)
    return fail(D, name + ": largest deviation from the exact distance " +
                       std::to_string(maxError));

  for (long i = lo[0]; i <= hi[0]; ++i)
    for (long j = lo[1]; j <= hi[1]; ++j)
      for (long k = lo[2]; k <= hi[2]; ++k) {
        const Point p{i * gridDelta, j * gridDelta, k * gridDelta};
        if (inRegion(p) && std::abs(sdf(p)) / gridDelta < halfWidth - 1e-6 &&
            !stored.count({i, j, k}))
          return fail(D, name + ": grid node inside the band is missing at (" +
                             std::to_string(i) + ", " + std::to_string(j) +
                             ", " + std::to_string(k) + ")");
      }

  std::cout << D << "D " << name << ": " << values.size()
            << " values, largest deviation " << maxError << std::endl;
  return EXIT_SUCCESS;
}

template <int D> std::array<long, 3> nodeRange(long n) {
  return {n, n, D == 3 ? n : 0L};
}

auto everywhere = [](const Point &) { return true; };

template <int D> int checkEuclideanBox(bool useWindingNumber) {
  const double a = 0.52, b = 0.35; // rotation about z, then about x
  RotatedBox<D> box;
  box.half = {3.15, 2.35, 1.85};
  box.centre = {0.31, -0.27, D == 3 ? 0.43 : 0.};
  if constexpr (D == 2)
    box.R = {{{std::cos(a), -std::sin(a), 0.},
              {std::sin(a), std::cos(a), 0.},
              {0., 0., 1.}}};
  else
    box.R = {
        {{std::cos(a), -std::sin(a), 0.},
         {std::cos(b) * std::sin(a), std::cos(b) * std::cos(a), -std::sin(b)},
         {std::sin(b) * std::sin(a), std::sin(b) * std::cos(a), std::cos(b)}}};

  auto levelSet = euclideanDomain<D>(0.5, 8);
  ls::FromSurfaceMesh<double, D> mesher(levelSet, box.mesh());
  mesher.setUseWindingNumber(useWindingNumber);
  mesher.apply();
  auto sdf = [&box](const Point &p) { return box.sdf(p); };
  const auto hi = nodeRange<D>(16);
  const std::string name =
      useWindingNumber ? "box, winding number sign" : "box, pseudo-normal sign";
  if (compareWithExact<D>(name, levelSet, sdf, everywhere,
                          {-hi[0], -hi[1], -hi[2]}, hi) != EXIT_SUCCESS)
    return EXIT_FAILURE;

  // Hermite data: closest point on the surface at the node's distance, and
  // the outward unit normal of a face containing it
  auto points = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, D>(levelSet, points).apply();
  const auto &values = *points->getCellData().getScalarData("LSValues");
  const auto *closest = points->getCellData().getVectorData(
      ls::FromSurfaceMesh<double, D>::closestPointsLabel);
  const auto *normals = points->getCellData().getVectorData(
      ls::FromSurfaceMesh<double, D>::surfaceNormalsLabel);
  if (closest == nullptr || normals == nullptr ||
      closest->size() != values.size() || normals->size() != values.size())
    return fail(D, name + ": closest points or normals missing");
  const double gridDelta = 0.5;
  for (std::size_t n = 0; n < values.size(); ++n) {
    const auto &x = points->getNodes()[n];
    const Point c{(*closest)[n][0], (*closest)[n][1], (*closest)[n][2]};
    const auto &nrm = (*normals)[n];
    double dist = 0., along = 0., unit = 0.;
    for (int i = 0; i < 3; ++i) {
      dist += (x[i] - c[i]) * (x[i] - c[i]);
      along += (x[i] - c[i]) * nrm[i];
      unit += nrm[i] * nrm[i];
    }
    if (std::abs(std::sqrt(dist) - std::abs(values[n]) * gridDelta) > 1e-9 ||
        std::abs(box.sdf(c)) > 1e-9)
      return fail(D, name + ": closest point inconsistent with the value");
    if (std::abs(unit - 1.) > 1e-9 || along * values[n] < -1e-9)
      return fail(D, name + ": normal is not an outward unit normal");
    // the normal is a face normal of the box and its face plane contains c
    const auto l = box.local(c);
    bool onFace = false;
    for (int j = 0; j < D; ++j) {
      double nj = 0.;
      for (int i = 0; i < 3; ++i)
        nj += box.R[i][j] * nrm[i];
      if (std::abs(std::abs(nj) - 1.) < 1e-9 &&
          std::abs(nj * l[j] - box.half[j]) < 1e-9)
        onFace = true;
    }
    if (!onFace)
      return fail(D, name + ": normal does not belong to the face at the "
                            "closest point");
  }
  std::cout << D << "D " << name << ": closest points and normals exact"
            << std::endl;
  return EXIT_SUCCESS;
}

template <int D> int checkEuclideanSphere() {
  auto levelSet = euclideanDomain<D>(0.5, 8);
  ls::VectorType<double, D> centre{};
  centre[0] = 0.31;
  centre[1] = -0.27;
  if constexpr (D == 3)
    centre[2] = 0.43;
  const double radius = 4.1;
  ls::MakeGeometry<double, D>(
      levelSet, ls::SmartPointer<ls::Sphere<double, D>>::New(centre, radius))
      .apply();
  auto sdf = [&](const Point &p) {
    double r2 = 0.;
    for (int i = 0; i < D; ++i)
      r2 += (p[i] - centre[i]) * (p[i] - centre[i]);
    return std::sqrt(r2) - radius;
  };
  const auto hi = nodeRange<D>(16);
  if (compareWithExact<D>("sphere", levelSet, sdf, everywhere,
                          {-hi[0], -hi[1], -hi[2]}, hi) != EXIT_SUCCESS)
    return EXIT_FAILURE;

  // stored closest points on the sphere and radial normals
  auto points = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToMesh<double, D>(levelSet, points).apply();
  const auto &values = *points->getCellData().getScalarData("LSValues");
  const auto *closest = points->getCellData().getVectorData(
      ls::FromSurfaceMesh<double, D>::closestPointsLabel);
  const auto *normals = points->getCellData().getVectorData(
      ls::FromSurfaceMesh<double, D>::surfaceNormalsLabel);
  if (closest == nullptr || normals == nullptr ||
      closest->size() != values.size())
    return fail(D, "sphere: closest points or normals missing");
  for (std::size_t n = 0; n < values.size(); ++n) {
    const auto &x = points->getNodes()[n];
    double toCentre = 0., toNode = 0., normalError = 0.;
    for (int i = 0; i < D; ++i) {
      toCentre +=
          ((*closest)[n][i] - centre[i]) * ((*closest)[n][i] - centre[i]);
      toNode += (x[i] - (*closest)[n][i]) * (x[i] - (*closest)[n][i]);
      const double radial = ((*closest)[n][i] - centre[i]) / radius;
      normalError = std::max(normalError, std::abs((*normals)[n][i] - radial));
    }
    if (std::abs(std::sqrt(toCentre) - radius) > 1e-9 ||
        std::abs(std::sqrt(toNode) - std::abs(values[n]) * 0.5) > 1e-9 ||
        normalError > 1e-9)
      return fail(D, "sphere: closest point or normal inconsistent");
  }
  std::cout << D << "D sphere: closest points and normals exact" << std::endl;
  return EXIT_SUCCESS;
}

// tilted plane through MakeGeometry: an open mesh reaching past the domain,
// with an infinite last direction
template <int D> int checkEuclideanPlane() {
  const double gridDelta = 0.5;
  auto levelSet = euclideanDomain<D>(gridDelta, 8, true);
  ls::VectorType<double, D> origin{}, normal{};
  origin[0] = 0.3;
  normal[0] = 0.3;
  if constexpr (D == 3) {
    origin[1] = -0.2;
    normal[1] = -0.2;
  }
  origin[D - 1] = 0.37;
  normal[D - 1] = 1.;
  double norm = 0.;
  for (int i = 0; i < D; ++i)
    norm += normal[i] * normal[i];
  for (int i = 0; i < D; ++i)
    normal[i] /= std::sqrt(norm);
  ls::MakeGeometry<double, D>(
      levelSet, ls::SmartPointer<ls::Plane<double, D>>::New(origin, normal))
      .apply();
  auto sdf = [&](const Point &p) {
    double d = 0.;
    for (int i = 0; i < D; ++i)
      d += normal[i] * (p[i] - origin[i]);
    return d;
  };
  // away from the reflective walls, where the closest points leave the mesh
  auto inner = [](const Point &p) {
    for (int i = 0; i < D - 1; ++i)
      if (std::abs(p[i]) > 7.)
        return false;
    return true;
  };
  std::array<long, 3> lo{-16, -16, 0}, hi{16, 16, 0};
  lo[D - 1] = -20;
  hi[D - 1] = 20;
  return compareWithExact<D>("plane, infinite boundary", levelSet, sdf, inner,
                             lo, hi);
}

// two materials in a 2D triangle mesh, split at y = 0.3
int checkEuclideanVolume() {
  constexpr int D = 2;
  const double gridDelta = 0.5;
  double bounds[2 * D] = {-4., 4., -4., 4.};
  ls::BoundaryConditionEnum boundaryCons[D] = {
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY,
      ls::BoundaryConditionEnum::REFLECTIVE_BOUNDARY};
  auto domain =
      ls::SmartPointer<LevelSetType<D>>::New(bounds, boundaryCons, gridDelta);

  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  const double nodes[6][2] = {{-4., -4.}, {4., -4.}, {4., 0.3},
                              {-4., 0.3}, {4., 4.},  {-4., 4.}};
  for (const auto &n : nodes)
    mesh->insertNextNode({n[0], n[1], 0.});
  const std::array<unsigned, 3> triangles[4] = {
      {0, 1, 2}, {0, 2, 3}, {3, 2, 4}, {3, 4, 5}};
  for (const auto &t : triangles)
    mesh->insertNextTriangle(t);
  mesh->getCellData().insertNextScalarData({0., 0., 1., 1.}, "Material");

  ls::FromVolumeMesh<double, D> reader(domain->getGrid(), mesh);
  reader.setNormalization(ls::NormalizationEnum::EUCLIDEAN);
  reader.setLevelSetWidth(8);
  reader.apply();
  auto levelSets = reader.getLevelSets();
  if (levelSets.empty())
    return fail(D, "volume mesh: no level sets created");
  for (const auto &l : levelSets)
    if (l->getNormalization() != ls::NormalizationEnum::EUCLIDEAN)
      return fail(D, "volume mesh: normalization was not passed on");

  // the first level set encloses material 0, below the interface
  auto sdf = [](const Point &p) { return p[1] - 0.3; };
  return compareWithExact<D>("volume mesh, material 0", levelSets[0], sdf,
                             everywhere, {-8, -8, 0}, {8, 8, 0});
}

// A box from MakeGeometry with its faces through grid nodes, at the default
// width 2: the band holds the nodes on the faces and those exactly one grid
// spacing away, so that dual contouring finds the surface
template <int D> int checkAlignedBox() {
  const double gridDelta = 0.5;
  auto levelSet = euclideanDomain<D>(gridDelta, 2);
  ls::VectorType<double, D> minCorner, maxCorner;
  for (int i = 0; i < D; ++i) {
    minCorner[i] = -2. + 0.5 * i;
    maxCorner[i] = 2. - 0.5 * i;
  }
  ls::MakeGeometry<double, D>(
      levelSet, ls::SmartPointer<ls::Box<double, D>>::New(minCorner, maxCorner))
      .apply();
  auto sdf = [&](const Point &p) {
    double outside = 0., inside = -1e300;
    for (int i = 0; i < D; ++i) {
      const double q = std::max(minCorner[i] - p[i], p[i] - maxCorner[i]);
      outside += (q > 0.) ? q * q : 0.;
      inside = std::max(inside, q);
    }
    return std::sqrt(outside) + std::min(inside, 0.);
  };
  if (compareWithExact<D>(
          "aligned box at width 2", levelSet, sdf,
          [](const Point &) { return true; }, nodeRange<D>(-6),
          nodeRange<D>(6)) != EXIT_SUCCESS)
    return EXIT_FAILURE;
  std::size_t away = 0;
  for (viennahrle::ConstSparseIterator<typename LevelSetType<D>::DomainType> it(
           levelSet->getDomain());
       !it.isFinished(); it.next())
    if (it.isDefined() && std::abs(std::abs(it.getValue()) - 1.) < 1e-12)
      ++away;
  auto mesh = ls::SmartPointer<ls::Mesh<>>::New();
  ls::ToDualContouringMesh<double, D>(levelSet, mesh).apply();
  double surfaceError = 0.;
  for (const auto &n : mesh->getNodes())
    surfaceError = std::max(surfaceError, std::abs(sdf({n[0], n[1], n[2]})));
  std::cout << D << "D aligned box at width 2: " << away
            << " nodes one grid spacing away, " << mesh->getNodes().size()
            << " dual contouring vertices, largest deviation "
            << surfaceError / gridDelta << std::endl;
  if (away == 0 || mesh->getNodes().empty() || surfaceError / gridDelta > 1e-9)
    return fail(D, "aligned box at width 2");
  return EXIT_SUCCESS;
}

int main() {
  omp_set_num_threads(4);

  if (checkNormalization<2>() != EXIT_SUCCESS ||
      checkNormalization<3>() != EXIT_SUCCESS)
    return EXIT_FAILURE;
  std::cout << "Normalization checks passed" << std::endl;

  for (bool useWindingNumber : {false, true})
    if (checkEuclideanBox<2>(useWindingNumber) != EXIT_SUCCESS ||
        checkEuclideanBox<3>(useWindingNumber) != EXIT_SUCCESS)
      return EXIT_FAILURE;
  if (checkEuclideanSphere<2>() != EXIT_SUCCESS ||
      checkEuclideanSphere<3>() != EXIT_SUCCESS)
    return EXIT_FAILURE;
  if (checkEuclideanPlane<2>() != EXIT_SUCCESS ||
      checkEuclideanPlane<3>() != EXIT_SUCCESS)
    return EXIT_FAILURE;
  if (checkEuclideanVolume() != EXIT_SUCCESS)
    return EXIT_FAILURE;
  if (checkAlignedBox<2>() != EXIT_SUCCESS ||
      checkAlignedBox<3>() != EXIT_SUCCESS)
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
