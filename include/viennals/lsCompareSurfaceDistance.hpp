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
/// Both surfaces are sampled uniformly, every segment in pieces and every
/// triangle in congruent sub-triangles no longer than the sampling spacing,
/// each sample weighted by its length or area. Every sample and every vertex
/// that an element uses is measured to the nearest point of the other surface,
/// a segment in 2D and a triangle in 3D, through an AABB tree. The results are
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
        const Point a = s.V.row(s.F(e, 0)), b = s.V.row(s.F(e, 1)),
                    c = s.V.row(s.F(e, 2));
        const Point u = b - a, v = c - a;
        const double area = 0.5 * u.cross(v).norm();
        const double longest = std::max({u.norm(), v.norm(), (c - b).norm()});
        const long n = std::max(1L, long(std::ceil(longest / spacing)));
        const double w = area / double(n * n);
        // the centroids of the n * n congruent sub-triangles, those with the
        // orientation of the triangle and the inverted ones between them
        for (long i = 0; i < n; ++i)
          for (long j = 0; i + j < n; ++j) {
            const Point p =
                a + u * ((i + 1. / 3.) / n) + v * ((j + 1. / 3.) / n);
            if (inBox(p)) {
              points.push_back(p);
              weights.push_back(w);
            }
            if (i + j + 1 < n) {
              const Point q =
                  a + u * ((i + 2. / 3.) / n) + v * ((j + 2. / 3.) / n);
              if (inBox(q)) {
                points.push_back(q);
                weights.push_back(w);
              }
            }
          }
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
