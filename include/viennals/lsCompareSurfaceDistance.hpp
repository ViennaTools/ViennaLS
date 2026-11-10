#pragma once

#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsMesh.hpp>
#include <lsPreCompileMacros.hpp>
#include <lsToDualContouringMesh.hpp>
#include <lsToSurfaceMesh.hpp>

#include <igl/AABB.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace viennals {

using namespace viennacore;

/// Distances between two surfaces, each given as a surface mesh (lines in 2D,
/// triangles in 3D, in the units of the domain) or as a level set of any grid
/// spacing and normalization. A Euclidean level set is extracted by dual
/// contouring and a Manhattan one by ToSurfaceMesh, so each level set is
/// measured by the surface it delivers.
///
/// Both surfaces are sampled uniformly, every segment in equal pieces and every
/// triangle by a grid along its longest edge clipped to it, no piece longer or
/// wider than the sampling spacing, each sample at the centroid of its piece
/// and weighted by its length or area. Every
/// sample and every vertex that an element uses is measured to the nearest
/// point of the other surface, a segment in 2D and a triangle in 3D, through an
/// AABB tree. The results are
/// - the forward distance, the mean distance of the target's samples to the
///   sample surface,
/// - the backward distance, the mean distance of the sample's samples to the
///   target surface,
/// - the Chamfer distance, the mean of the forward and the backward distance,
/// - the RMS distance, the root of the mean of both mean squares, and
/// - the maximum distance over the samples and vertices of both surfaces (the
///   Hausdorff distance up to the sampling).
/// These do not depend on the grids or the vertices of the two surfaces, so a
/// coarse run compares with a fine reference or with an exact surface. An
/// optional box restricts the samples and vertices of both surfaces, which are
/// still measured to the whole other surface.
template <class T, int D = 2> class CompareSurfaceDistance {
  using Point = Eigen::Matrix<double, 1, D>;

  SmartPointer<Domain<T, D>> levelSetTarget = nullptr;
  SmartPointer<Domain<T, D>> levelSetSample = nullptr;
  SmartPointer<Mesh<T>> meshTarget = nullptr;
  SmartPointer<Mesh<T>> meshSample = nullptr;
  T samplingSpacing = 0.;
  bool useBox = false;
  std::array<T, D> boxMin{};
  std::array<T, D> boxMax{};

  // results
  T forwardDistance = 0.;  // target to sample
  T backwardDistance = 0.; // sample to target
  T chamferDistance = 0.;
  T rmsChamferDistance = 0.;
  T maxDistance = 0.;
  std::size_t numTargetPoints = 0; // samples of the target
  std::size_t numSamplePoints = 0; // samples of the sample

  // optional output: the surfaces with the distance of every vertex
  SmartPointer<Mesh<T>> outputMeshTarget = nullptr;
  SmartPointer<Mesh<T>> outputMeshSample = nullptr;

  // a surface as the matrices of libigl, its vertices that an element uses
  // and the tree of its elements
  struct Surface {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    std::vector<int> usedVertices;
    igl::AABB<Eigen::MatrixXd, D> tree;
  };

  // the mean, the mean square and the largest distance of the samples of one
  // surface to the other, and their number
  struct Measure {
    double mean = 0., meanSquare = 0., maximum = 0.;
    std::size_t samples = 0;
  };

  static SmartPointer<Mesh<T>> surfaceOf(SmartPointer<Domain<T, D>> levelSet) {
    auto mesh = SmartPointer<Mesh<T>>::New();
    if (levelSet->getNormalization() == NormalizationEnum::EUCLIDEAN) {
      ToDualContouringMesh<T, D>(levelSet, mesh).apply();
    } else {
      auto working = levelSet;
      if (levelSet->getLevelSetWidth() < 2) {
        working = SmartPointer<Domain<T, D>>::New(levelSet);
        Expand<T, D>(working, 2).apply();
      }
      ToSurfaceMesh<T, D>(working, mesh).apply();
    }
    return mesh;
  }

  static void fill(const Mesh<T> &mesh, Surface &surface) {
    const auto &nodes = mesh.getNodes();
    surface.V.resize(nodes.size(), D);
    for (std::size_t p = 0; p < nodes.size(); ++p)
      for (int i = 0; i < D; ++i)
        surface.V(p, i) = nodes[p][i];
    std::vector<char> used(nodes.size(), 0);
    if constexpr (D == 2) {
      surface.F.resize(mesh.lines.size(), 2);
      for (std::size_t e = 0; e < mesh.lines.size(); ++e)
        for (int k = 0; k < 2; ++k) {
          surface.F(e, k) = mesh.lines[e][k];
          used[mesh.lines[e][k]] = 1;
        }
    } else {
      surface.F.resize(mesh.triangles.size(), 3);
      for (std::size_t e = 0; e < mesh.triangles.size(); ++e)
        for (int k = 0; k < 3; ++k) {
          surface.F(e, k) = mesh.triangles[e][k];
          used[mesh.triangles[e][k]] = 1;
        }
    }
    surface.usedVertices.clear();
    for (std::size_t p = 0; p < nodes.size(); ++p)
      if (used[p])
        surface.usedVertices.push_back(int(p));
    if (surface.F.rows() > 0)
      surface.tree.init(surface.V, surface.F);
  }

  bool inBox(const Point &p) const {
    if (!useBox)
      return true;
    for (int i = 0; i < D; ++i)
      if (p(i) < boxMin[i] || p(i) > boxMax[i])
        return false;
    return true;
  }

  // the samples of a triangle: in the frame of its longest edge (x along it, y
  // across it in its plane) a grid of at most the spacing in both directions,
  // every cell clipped to the triangle and sampled at the centroid of the
  // piece, weighted by its area. The pieces cover the triangle, and a long thin
  // triangle gets a single row of cells.
  void triangleSamples(const Point &p0, const Point &p1, const Point &p2,
                       double spacing, std::vector<Point> &points,
                       std::vector<double> &weights) const {
    if constexpr (D == 3) {
      Point a = p0, b = p1, c = p2;
      const double l0 = (p1 - p0).squaredNorm(), l1 = (p2 - p1).squaredNorm(),
                   l2 = (p0 - p2).squaredNorm();
      if (l1 >= l0 && l1 >= l2) {
        a = p1;
        b = p2;
        c = p0;
      } else if (l2 >= l0 && l2 >= l1) {
        a = p2;
        b = p0;
        c = p1;
      }
      const double length = (b - a).norm();
      if (!(length > 0.))
        return;
      const Point ex = (b - a) / length;
      const Point across = (c - a) - ex * (c - a).dot(ex);
      const double height = across.norm();
      if (!(height > 0.))
        return;
      const Point ey = across / height;
      // the triangle (0, 0), (length, 0), (cx, height), counterclockwise
      const double cx = (c - a).dot(ex);
      const double x0 = std::min(0., cx), x1 = std::max(length, cx);
      const long nx = std::max(1L, long(std::ceil((x1 - x0) / spacing)));
      const long ny = std::max(1L, long(std::ceil(height / spacing)));
      const double dx = (x1 - x0) / nx, dy = height / ny;
      using Corner = std::array<double, 2>;
      // the part of a polygon left of the directed line from u to v
      auto clip = [](const std::vector<Corner> &polygon, const Corner &u,
                     const Corner &v) {
        std::vector<Corner> result;
        auto side = [&](const Corner &q) {
          return (v[0] - u[0]) * (q[1] - u[1]) - (v[1] - u[1]) * (q[0] - u[0]);
        };
        for (std::size_t k = 0; k < polygon.size(); ++k) {
          const Corner &p = polygon[k], &q = polygon[(k + 1) % polygon.size()];
          const double sp = side(p), sq = side(q);
          if (sp >= 0.)
            result.push_back(p);
          if ((sp >= 0.) != (sq >= 0.)) {
            const double t = sp / (sp - sq);
            result.push_back(
                {p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])});
          }
        }
        return result;
      };
      const Corner A{0., 0.}, B{length, 0.}, C{cx, height};
      for (long j = 0; j < ny; ++j) {
        const double ya = j * dy, yb = (j + 1) * dy;
        // the cells of the row that the triangle reaches
        const double xl = std::min(cx * ya / height, cx * yb / height);
        const double xr = std::max(length + (cx - length) * ya / height,
                                   length + (cx - length) * yb / height);
        const long i0 = std::max(0L, long(std::floor((xl - x0) / dx)));
        const long i1 = std::min(nx - 1, long(std::floor((xr - x0) / dx)));
        for (long i = i0; i <= i1; ++i) {
          const double xa = x0 + i * dx, xb = x0 + (i + 1) * dx;
          std::vector<Corner> cell{{xa, ya}, {xb, ya}, {xb, yb}, {xa, yb}};
          cell = clip(cell, B, C);
          if (cell.size() >= 3)
            cell = clip(cell, C, A);
          if (cell.size() < 3)
            continue;
          double area = 0., gx = 0., gy = 0.;
          for (std::size_t k = 0; k < cell.size(); ++k) {
            const Corner &p = cell[k], &q = cell[(k + 1) % cell.size()];
            const double cross = p[0] * q[1] - q[0] * p[1];
            area += cross;
            gx += (p[0] + q[0]) * cross;
            gy += (p[1] + q[1]) * cross;
          }
          area *= 0.5;
          if (!(area > 0.))
            continue;
          const Point point =
              a + ex * (gx / (6. * area)) + ey * (gy / (6. * area));
          if (inBox(point)) {
            points.push_back(point);
            weights.push_back(area);
          }
        }
      }
    }
  }

  // the samples of the elements of a surface and their weights
  void samplesOf(const Surface &s, double spacing, std::vector<Point> &points,
                 std::vector<double> &weights) const {
    points.clear();
    weights.clear();
    for (Eigen::Index e = 0; e < s.F.rows(); ++e) {
      if constexpr (D == 2) {
        const Point a = s.V.row(s.F(e, 0)), b = s.V.row(s.F(e, 1));
        const double length = (b - a).norm();
        const long n = std::max(1L, long(std::ceil(length / spacing)));
        for (long k = 0; k < n; ++k) {
          const Point p = a + (b - a) * ((k + 0.5) / n);
          if (inBox(p)) {
            points.push_back(p);
            weights.push_back(length / n);
          }
        }
      } else {
        triangleSamples(s.V.row(s.F(e, 0)), s.V.row(s.F(e, 1)),
                        s.V.row(s.F(e, 2)), spacing, points, weights);
      }
    }
  }

  // the samples of one surface measured to the other, and the distance of
  // every vertex of the first surface to the other (NaN where unused)
  Measure measure(const Surface &from, const Surface &to, double spacing,
                  std::vector<double> &vertexDistances) const {
    std::vector<Point> points;
    std::vector<double> weights;
    samplesOf(from, spacing, points, weights);
    Measure m;
    m.samples = points.size();
    double sum = 0., sumSquares = 0., sumWeights = 0., largest = 0.;
    const long n = long(points.size());
#pragma omp parallel for reduction(+ : sum, sumSquares, sumWeights)            \
    reduction(max : largest) schedule(dynamic, 1024)
    for (long k = 0; k < n; ++k) {
      int nearest = -1;
      Point closest;
      const double d = std::sqrt(
          to.tree.squared_distance(to.V, to.F, points[k], nearest, closest));
      sum += weights[k] * d;
      sumSquares += weights[k] * d * d;
      sumWeights += weights[k];
      largest = std::max(largest, d);
    }
    vertexDistances.assign(from.V.rows(),
                           std::numeric_limits<double>::quiet_NaN());
    const long nv = long(from.usedVertices.size());
#pragma omp parallel for reduction(max : largest)
    for (long k = 0; k < nv; ++k) {
      const int p = from.usedVertices[k];
      const Point x = from.V.row(p);
      int nearest = -1;
      Point closest;
      const double d =
          std::sqrt(to.tree.squared_distance(to.V, to.F, x, nearest, closest));
      vertexDistances[p] = d;
      if (inBox(x))
        largest = std::max(largest, d);
    }
    if (sumWeights > 0.) {
      m.mean = sum / sumWeights;
      m.meanSquare = sumSquares / sumWeights;
    }
    m.maximum = largest;
    return m;
  }

  static void writeOutput(SmartPointer<Mesh<T>> output, const Mesh<T> &mesh,
                          const std::vector<double> &distances,
                          const std::string &name) {
    output->clear();
    output->nodes = mesh.nodes;
    if constexpr (D == 2)
      output->lines = mesh.lines;
    else
      output->triangles = mesh.triangles;
    std::vector<T> data(distances.begin(), distances.end());
    output->pointData.insertNextScalarData(std::move(data), name);
    output->minimumExtent = mesh.minimumExtent;
    output->maximumExtent = mesh.maximumExtent;
  }

  void setFailed(T value) {
    forwardDistance = backwardDistance = chamferDistance = rmsChamferDistance =
        maxDistance = value;
    numTargetPoints = numSamplePoints = 0;
  }

public:
  CompareSurfaceDistance() = default;

  CompareSurfaceDistance(SmartPointer<Domain<T, D>> passedLevelSetTarget,
                         SmartPointer<Domain<T, D>> passedLevelSetSample)
      : levelSetTarget(passedLevelSetTarget),
        levelSetSample(passedLevelSetSample) {}

  CompareSurfaceDistance(SmartPointer<Mesh<T>> passedMeshTarget,
                         SmartPointer<Mesh<T>> passedMeshSample)
      : meshTarget(passedMeshTarget), meshSample(passedMeshSample) {}

  /// Set the target as a level set, extracted when applied
  void setLevelSetTarget(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetTarget = passedLevelSet;
    meshTarget = nullptr;
  }

  /// Set the sample as a level set, extracted when applied
  void setLevelSetSample(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSetSample = passedLevelSet;
    meshSample = nullptr;
  }

  /// Set the target as a surface mesh (lines in 2D, triangles in 3D)
  void setMeshTarget(SmartPointer<Mesh<T>> passedMesh) {
    meshTarget = passedMesh;
    levelSetTarget = nullptr;
  }

  /// Set the sample as a surface mesh (lines in 2D, triangles in 3D)
  void setMeshSample(SmartPointer<Mesh<T>> passedMesh) {
    meshSample = passedMesh;
    levelSetSample = nullptr;
  }

  /// Set the largest length (2D) or edge (3D) of a sample element. With 0,
  /// the default, it is a fifth of the smallest grid spacing of the level
  /// sets given, or a fifth of the mean edge length of the meshes.
  void setSamplingSpacing(T passedSpacing) { samplingSpacing = passedSpacing; }

  /// Restrict the samples and vertices of both surfaces to a box
  void setBox(const std::array<T, D> &minCorner,
              const std::array<T, D> &maxCorner) {
    boxMin = minCorner;
    boxMax = maxCorner;
    useBox = true;
  }

  /// Measure the whole surfaces again
  void clearBox() { useBox = false; }

  /// Set output mesh for the target surface with the distance of every
  /// vertex to the sample surface
  void setOutputMeshTarget(SmartPointer<Mesh<T>> passedMesh) {
    outputMeshTarget = passedMesh;
  }

  /// Set output mesh for the sample surface with the distance of every
  /// vertex to the target surface
  void setOutputMeshSample(SmartPointer<Mesh<T>> passedMesh) {
    outputMeshSample = passedMesh;
  }

  /// Compute the distances
  void apply() {
    auto target = meshTarget;
    if (target == nullptr && levelSetTarget != nullptr)
      target = surfaceOf(levelSetTarget);
    auto sample = meshSample;
    if (sample == nullptr && levelSetSample != nullptr)
      sample = surfaceOf(levelSetSample);
    if (target == nullptr || sample == nullptr) {
      VIENNACORE_LOG_WARNING("CompareSurfaceDistance: a surface is missing.");
      setFailed(std::numeric_limits<T>::quiet_NaN());
      return;
    }

    Surface t, s;
    fill(*target, t);
    fill(*sample, s);
    if (t.F.rows() == 0 || s.F.rows() == 0) {
      VIENNACORE_LOG_WARNING(
          "CompareSurfaceDistance: a surface has no elements.");
      setFailed(std::numeric_limits<T>::infinity());
      return;
    }

    double spacing = samplingSpacing;
    if (spacing <= 0.) {
      double grid = std::numeric_limits<double>::infinity();
      for (const auto &ls : {levelSetTarget, levelSetSample})
        if (ls != nullptr)
          grid = std::min(grid, double(ls->getGrid().getGridDelta()));
      if (std::isfinite(grid)) {
        spacing = grid / 5.;
      } else {
        double sum = 0.;
        long count = 0;
        for (const Surface *surface : {&t, &s})
          for (Eigen::Index e = 0; e < surface->F.rows(); ++e)
            for (int k = 0; k < D; ++k) {
              sum += (surface->V.row(surface->F(e, (k + 1) % D)) -
                      surface->V.row(surface->F(e, k)))
                         .norm();
              ++count;
            }
        spacing = sum / std::max(1L, count) / 5.;
      }
    }
    if (!(spacing > 0.)) {
      VIENNACORE_LOG_WARNING(
          "CompareSurfaceDistance: no positive sampling spacing.");
      setFailed(std::numeric_limits<T>::quiet_NaN());
      return;
    }

    std::vector<double> targetDistances, sampleDistances;
    const Measure forward = measure(t, s, spacing, targetDistances);
    const Measure backward = measure(s, t, spacing, sampleDistances);
    numTargetPoints = forward.samples;
    numSamplePoints = backward.samples;
    if (forward.samples == 0 || backward.samples == 0) {
      VIENNACORE_LOG_WARNING(
          "CompareSurfaceDistance: no samples, the box misses a surface.");
      setFailed(std::numeric_limits<T>::quiet_NaN());
      numTargetPoints = forward.samples;
      numSamplePoints = backward.samples;
      return;
    }
    forwardDistance = T(forward.mean);
    backwardDistance = T(backward.mean);
    chamferDistance = T(0.5 * (forward.mean + backward.mean));
    rmsChamferDistance =
        T(std::sqrt(0.5 * (forward.meanSquare + backward.meanSquare)));
    maxDistance = T(std::max(forward.maximum, backward.maximum));

    if (outputMeshTarget != nullptr)
      writeOutput(outputMeshTarget, *target, targetDistances,
                  "DistanceToSample");
    if (outputMeshSample != nullptr)
      writeOutput(outputMeshSample, *sample, sampleDistances,
                  "DistanceToTarget");
  }

  /// Get the forward distance (mean distance of the target to the sample)
  T getForwardDistance() const { return forwardDistance; }

  /// Get the backward distance (mean distance of the sample to the target)
  T getBackwardDistance() const { return backwardDistance; }

  /// Get the Chamfer distance (mean of forward and backward)
  T getChamferDistance() const { return chamferDistance; }

  /// Get the RMS distance (root of the mean of both mean squares)
  T getRMSChamferDistance() const { return rmsChamferDistance; }

  /// Get the maximum distance of a sample or vertex in either direction
  T getMaxDistance() const { return maxDistance; }

  /// Get the number of samples of the target surface
  std::size_t getNumTargetPoints() const { return numTargetPoints; }

  /// Get the number of samples of the sample surface
  std::size_t getNumSamplePoints() const { return numSamplePoints; }
};

// Add template specializations for this class
PRECOMPILE_PRECISION_DIMENSION(CompareSurfaceDistance)

} // namespace viennals
