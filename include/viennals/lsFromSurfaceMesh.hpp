#pragma once

#include <lsPreCompileMacros.hpp>

#include <limits>
#include <vector>

#include <hrleTypes.hpp>

#include <lsDomain.hpp>
#include <lsMesh.hpp>

#include <Eigen/Core>
#include <igl/signed_distance.h>

namespace viennals {

using namespace viennacore;

/// Construct a level set from an explicit mesh.
template <class T, int D> class FromSurfaceMesh {
  using hrleIndexType = viennahrle::IndexType;

  /// Class defining a box used in ray tracing optimisation
  class box {
    VectorType<hrleIndexType, D - 1> xMin, xMax;

  public:
    /// Sets xMin, xMax to the lowest and highest value of the two vectors for
    /// each dimension respectively.
    box(const VectorType<hrleIndexType, D - 1> &idx0,
        const VectorType<hrleIndexType, D - 1> &idx1) {
      xMin = Min(idx0, idx1);
      xMax = Max(idx0, idx1);
    }

    /// Creates a copy of the passed box.
    box(const box &box) {
      xMin = box.xMin;
      xMax = box.xMax;
    }

    /// Sets both xMin and Xmax to the passed vector.
    explicit box(const VectorType<hrleIndexType, D - 1> &idx) {
      xMin = idx;
      xMax = idx;
    }

    box &operator=(const box &b) = default;

    /// Checks whether there are points in the box and returns result.
    bool is_empty() const {
      for (int i = 0; i < D - 1; i++)
        if (xMax[i] < xMin[i])
          return true;
      return false;
    }

    /// Returns vector of lowest point of box in each dimension.
    const VectorType<hrleIndexType, D - 1> &min() const { return xMin; }
    /// Returns vector of highest point of box in each dimension.
    const VectorType<hrleIndexType, D - 1> &max() const { return xMax; }

    /// Iterator over all grid points, contained by a box.
    class iterator {
      VectorType<hrleIndexType, D - 1> pos;
      const box &b;

    public:
      explicit iterator(const box &bx) : pos(bx.min()), b(bx) {}

      iterator &operator++() {
        int i;
        for (i = 0; i < D - 2; i++) {
          ++pos[i];
          if (pos[i] <= b.xMax[i]) {
            break;
          } else {
            pos[i] = b.xMin[i];
          }
        }
        if (i == D - 2)
          ++pos[i];
        return *this;
      }

      iterator operator++(int) {
        iterator tmp = *this;
        ++(*this);
        return tmp;
      }

      bool is_finished() const { return (pos[D - 2] > b.xMax[D - 2]); }

      const VectorType<hrleIndexType, D - 1> &operator*() const { return pos; }
    };
  };

  SmartPointer<Domain<T, D>> levelSet = SmartPointer<Domain<T, D>>::New();
  SmartPointer<Mesh<T>> mesh = SmartPointer<Mesh<T>>::New();
  // bool removeBoundaryTriangles = true;
  std::array<bool, 3> removeBoundaryTriangles{true, true, true};
  T boundaryEps = 1e-5;
  T distanceEps = 1e-4;
  T signEps = 1e-6;
  // rounding of the precision, below the fixed tolerances with double
  static constexpr T roundOff = T(1024) * std::numeric_limits<T>::epsilon();
  bool useWindingNumber = false;

  /// this function checks if a line containing the point with coordinates
  /// "Point" and parallel to axis "dir" (x=0, y=1) intersects the surface
  /// element given by the nodes "c"
  /// If there is
  /// an intersection this function returns true, otherwise false the
  /// intersection coordinate is returned by "intersection"
  int calculateGridlineTriangleIntersection(VectorType<T, 2> Point,
                                            const VectorType<T, 2> *c, int dir,
                                            T &intersection) {

    bool inside_pos = true;
    bool inside_neg = true;

    T A[2];

    const int dirA = (dir + 1) % 2;

    for (int k = 0; k < 2; k++) {
      A[k] = c[(k + 1) % 2][dirA] - Point[dirA];
      if (k == dir)
        A[k] = -A[k];
      if (A[k] < T(0))
        inside_pos = false;
      if (A[k] > T(0))
        inside_neg = false;
    }

    if ((!inside_pos && inside_neg) || (inside_pos && !inside_neg)) {
      T sum = A[0] + A[1];
      int k = 0;
      for (int i = 1; i < 2; ++i) {
        if (inside_pos) {
          if (A[i] > A[k])
            k = i;
        } else {
          if (A[i] < A[k])
            k = i;
        }
      }
      intersection = c[k][dir] +
                     (c[(k + 1) % 2][dir] - c[k][dir]) * (A[(k + 1) % 2] / sum);
      return (inside_pos) ? 1 : -1;
    } else {
      return 0;
    }
  }

  /// this function checks if a line containing the point with coordinates
  /// "Point" and parallel to axis "dir" (x=0, y=1, z=2) intersects the surface
  /// element given by the nodes "c".
  /// If there is
  /// an intersection this function returns true, otherwise false the
  /// intersection coordinate is returned by "intersection"
  int calculateGridlineTriangleIntersection(VectorType<T, 3> Point,
                                            const VectorType<T, 3> *c, int dir,
                                            T &intersection) {

    bool inside_pos = true;
    bool inside_neg = true;

    T A[3];

    const int dirA = (dir + 1) % 3;
    const int dirB = (dir + 2) % 3;

    for (int k = 0; k < 3; k++) {
      bool swapped =
          (c[(k + 1) % 3] <
           c[(k + 2) % 3]); // necessary to guarantee anti commutativity
      const VectorType<T, 3> &v1 = (swapped) ? c[(k + 2) % 3] : c[(k + 1) % 3];
      const VectorType<T, 3> &v2 = (swapped) ? c[(k + 1) % 3] : c[(k + 2) % 3];

      A[k] = (v1[dirA] - Point[dirA]) * (v2[dirB] - Point[dirB]) -
             (v2[dirA] - Point[dirA]) * (v1[dirB] - Point[dirB]);

      if (swapped)
        A[k] = -A[k];

      if (A[k] < T(0))
        inside_pos = false;
      if (A[k] > T(0))
        inside_neg = false;
    }

    if ((!inside_pos && inside_neg) || (inside_pos && !inside_neg)) {
      T sum = A[0] + A[1] + A[2];
      int k = 0;
      for (int i = 1; i < 3; ++i) {
        if (inside_pos) {
          if (A[i] > A[k])
            k = i;
        } else {
          if (A[i] < A[k])
            k = i;
        }
      }
      intersection =
          c[k][dir] +
          (c[(k + 1) % 3][dir] - c[k][dir]) * (A[(k + 1) % 3] / sum) +
          (c[(k + 2) % 3][dir] - c[k][dir]) * (A[(k + 2) % 3] / sum);
      return (inside_pos) ? 1 : -1;
    } else {
      return 0;
    }
  }

  /// Euclidean normalization: every grid node within half the level set width
  /// of the mesh receives its exact signed distance to the mesh (libigl, sign
  /// from angle-weighted pseudo-normals or from the generalized winding
  /// number). Elements are filtered at the domain boundary with the same
  /// rules as the grid-line construction.
  void applyEuclidean(const bool removeBoundaries[D]) {
    using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    const auto &grid = levelSet->getGrid();
    const T gridDelta = grid.getGridDelta();
    const int width = std::max(levelSet->getLevelSetWidth(), 2);
    const T halfWidth = T(0.5) * width; // in grid spacings

    VectorType<T, D> gridMin, gridMax;
    for (unsigned i = 0; i < D; ++i) {
      gridMin[i] = grid.getMinIndex(i) * gridDelta;
      gridMax[i] = grid.getMaxIndex(i) * gridDelta;
    }

    // node coordinates in grid spacings, snapped onto the domain boundary
    const auto &nodes = mesh->getNodes();
    MatrixType V(nodes.size(), D);
    for (std::size_t n = 0; n < nodes.size(); ++n) {
      for (int i = 0; i < D; ++i) {
        T x = nodes[n][i];
        if (std::abs(x - gridMin[i]) < boundaryEps * gridDelta)
          x = gridMin[i];
        if (std::abs(x - gridMax[i]) < boundaryEps * gridDelta)
          x = gridMax[i];
        V(n, i) = x / gridDelta;
      }
    }

    // drop elements lying in a removed domain boundary
    const auto &elements = mesh->template getElements<D>();
    std::vector<std::size_t> kept;
    for (std::size_t e = 0; e < elements.size(); ++e) {
      bool removeElement = false;
      for (int i = 0; i < D; ++i) {
        bool atMin = true, atMax = true;
        for (int q = 0; q < D; ++q) {
          const T x = V(elements[e][q], i) * gridDelta;
          if (x > gridMin[i])
            atMin = false;
          if (x < gridMax[i])
            atMax = false;
        }
        if (removeBoundaries[i] && (atMin || atMax))
          removeElement = true;
      }
      if (removeElement)
        continue;
      // elements with a vertex at no finite position are left out
      bool finite = true;
      for (int q = 0; q < D; ++q)
        for (int i = 0; i < D; ++i)
          finite = finite && std::isfinite(V(elements[e][q], i));
      if (!finite)
        continue;
      // elements without area add no surface and their normals are round-off
      if constexpr (D == 2) {
        if (std::hypot(V(elements[e][1], 0) - V(elements[e][0], 0),
                       V(elements[e][1], 1) - V(elements[e][0], 1)) <
            std::max(T(1e-10), roundOff))
          continue;
      } else {
        Vec3D<T> edge[3];
        for (int i = 0; i < 3; ++i) {
          edge[0][i] = V(elements[e][1], i) - V(elements[e][0], i);
          edge[1][i] = V(elements[e][2], i) - V(elements[e][0], i);
          edge[2][i] = V(elements[e][2], i) - V(elements[e][1], i);
        }
        const T longest =
            std::max({Norm(edge[0]), Norm(edge[1]), Norm(edge[2])});
        if (Norm(CrossProduct(edge[0], edge[1])) <=
            std::max(T(1e-10), roundOff) * longest * longest)
          continue;
      }
      kept.push_back(e);
    }
    if (kept.empty()) {
      // nothing inside the domain, empty level set as in the grid-line path
      std::vector<std::pair<viennahrle::Index<D>, T>> noPoints;
      levelSet->insertPoints(noPoints);
      levelSet->getDomain().segment();
      levelSet->finalize(width);
      return;
    }

    // libigl rotates 2D segments clockwise to get their normal, ViennaLS lines
    // have the outside on their left, so 2D segments are reversed
    Eigen::MatrixXi F(kept.size(), D);
    for (std::size_t k = 0; k < kept.size(); ++k)
      for (int q = 0; q < D; ++q)
        F(k, q) = elements[kept[k]][(D == 2) ? 1 - q : q];

    // grid nodes within halfWidth of any element
    viennahrle::Index<D> lo, hi;
    std::vector<std::pair<viennahrle::Index<D>, viennahrle::Index<D>>> boxes;
    for (std::size_t k = 0; k < kept.size(); ++k) {
      viennahrle::Index<D> bMin, bMax;
      for (int i = 0; i < D; ++i) {
        T eMin = V(F(k, 0), i), eMax = V(F(k, 0), i);
        for (int q = 1; q < D; ++q) {
          eMin = std::min(eMin, V(F(k, q), i));
          eMax = std::max(eMax, V(F(k, q), i));
        }
        bMin[i] = static_cast<hrleIndexType>(std::ceil(eMin - halfWidth));
        bMax[i] = static_cast<hrleIndexType>(std::floor(eMax + halfWidth));
        if (removeBoundaries[i]) {
          bMin[i] = std::max(bMin[i], grid.getMinIndex(i));
          bMax[i] = std::min(bMax[i], grid.getMaxIndex(i));
        }
        lo[i] = (k == 0) ? bMin[i] : std::min(lo[i], bMin[i]);
        hi[i] = (k == 0) ? bMax[i] : std::max(hi[i], bMax[i]);
      }
      boxes.emplace_back(bMin, bMax);
    }
    std::size_t total = 1;
    std::array<std::size_t, D> extent{};
    for (int i = 0; i < D; ++i) {
      extent[i] = (hi[i] >= lo[i]) ? std::size_t(hi[i] - lo[i] + 1) : 0;
      total *= extent[i];
    }
    std::vector<char> marked(total, 0);
    for (const auto &[bMin, bMax] : boxes) {
      viennahrle::Index<D> idx = bMin;
      bool empty = false;
      for (int i = 0; i < D; ++i)
        empty |= bMax[i] < bMin[i];
      while (!empty) {
        std::size_t flat = 0;
        for (int i = D - 1; i >= 0; --i)
          flat = flat * extent[i] + std::size_t(idx[i] - lo[i]);
        marked[flat] = 1;
        int i = 0;
        for (; i < D; ++i) {
          if (++idx[i] <= bMax[i])
            break;
          idx[i] = bMin[i];
        }
        if (i == D)
          break;
      }
    }
    std::vector<viennahrle::Index<D>> candidates;
    for (std::size_t flat = 0; flat < total; ++flat) {
      if (!marked[flat])
        continue;
      viennahrle::Index<D> idx;
      std::size_t rest = flat;
      for (int i = 0; i < D; ++i) {
        idx[i] = lo[i] + hrleIndexType(rest % extent[i]);
        rest /= extent[i];
      }
      candidates.push_back(idx);
    }

    MatrixType P(candidates.size(), D);
    for (std::size_t c = 0; c < candidates.size(); ++c)
      for (int i = 0; i < D; ++i)
        P(c, i) = candidates[c][i];

    // libigl leaves out points at its bound itself; nodes at exactly half
    // the width, as next to a face through grid nodes, belong to the band.
    // It runs in double precision also for float level sets: in single
    // precision a closest point next to an edge of the mesh can fall on the
    // wrong element, and the sign with it
    const double searchBound =
        halfWidth * (1. + std::max(1e-9, double(roundOff)));
    Eigen::Matrix<double, Eigen::Dynamic, 1> iglDistances;
    Eigen::VectorXi I;
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Pd =
        P.template cast<double>();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> iglClosest,
        iglNormals;
    igl::signed_distance(
        Pd, V.template cast<double>().eval(), F,
        useWindingNumber ? igl::SIGNED_DISTANCE_TYPE_WINDING_NUMBER
                         : igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL,
        -searchBound, searchBound, iglDistances, I, iglClosest, iglNormals);
    // with winding-number signs libigl scales the distance by 1 - 2w, whose
    // rounding differs from call to call: the distance to the closest point,
    // with the sign alone
    if (useWindingNumber)
      for (Eigen::Index c = 0; c < iglDistances.size(); ++c)
        if (std::isfinite(iglDistances(c)))
          iglDistances(c) = std::copysign(
              (Pd.row(c) - iglClosest.row(c)).norm(), iglDistances(c));
    const Eigen::Matrix<T, Eigen::Dynamic, 1> S =
        iglDistances.template cast<T>();
    const MatrixType C = iglClosest.template cast<T>();

    // outward unit normal of an element
    auto normalOf = [&](int f) {
      Vec3D<T> normal;
      if constexpr (D == 2) {
        // segments are reversed above, (dy, -dx) points outward
        normal = Vec3D<T>{V(F(f, 1), 1) - V(F(f, 0), 1),
                          V(F(f, 0), 0) - V(F(f, 1), 0), T(0)};
      } else {
        Vec3D<T> e1, e2;
        for (int i = 0; i < 3; ++i) {
          e1[i] = V(F(f, 1), i) - V(F(f, 0), i);
          e2[i] = V(F(f, 2), i) - V(F(f, 0), i);
        }
        normal = CrossProduct(e1, e2);
      }
      const T length = Norm(normal);
      if (length > T(0))
        for (int i = 0; i < 3; ++i)
          normal[i] /= length;
      return normal;
    };
    std::vector<std::vector<int>> elementsAt(V.rows());
    for (int f = 0; f < int(F.rows()); ++f)
      for (int q = 0; q < D; ++q)
        elementsAt[F(f, q)].push_back(f);

    // a closest point on an edge or vertex of the mesh is shared by several
    // elements; the one whose normal points most to the side of the node is
    // chosen, so that the stored normal does not depend on the element order
    auto sharedElement = [&](std::size_t c) {
      int f = I(c);
      Vec3D<T> dir{};
      for (int i = 0; i < D; ++i)
        dir[i] = P(c, i) - C(c, i);
      const T distance = Norm(dir);
      if (distance < std::max(T(1e-12), roundOff))
        return f;
      const T tolerance = std::max(T(1e-9), roundOff);
      std::vector<int> at; // vertices of f the closest point lies on or between
      if constexpr (D == 2) {
        T t = 0., len2 = 0.;
        for (int i = 0; i < 2; ++i) {
          const T e = V(F(f, 1), i) - V(F(f, 0), i);
          t += (C(c, i) - V(F(f, 0), i)) * e;
          len2 += e * e;
        }
        t /= std::max(len2, T(1e-300));
        if (t <= tolerance)
          at = {F(f, 0)};
        else if (t >= 1 - tolerance)
          at = {F(f, 1)};
      } else {
        Vec3D<T> v0, v1, v2;
        for (int i = 0; i < 3; ++i) {
          v0[i] = V(F(f, 1), i) - V(F(f, 0), i);
          v1[i] = V(F(f, 2), i) - V(F(f, 0), i);
          v2[i] = C(c, i) - V(F(f, 0), i);
        }
        const T d00 = DotProduct(v0, v0), d01 = DotProduct(v0, v1),
                d11 = DotProduct(v1, v1), d20 = DotProduct(v2, v0),
                d21 = DotProduct(v2, v1);
        const T denominator = d00 * d11 - d01 * d01;
        if (std::abs(denominator) < T(1e-300))
          return f;
        const T l1 = (d11 * d20 - d01 * d21) / denominator;
        const T l2 = (d00 * d21 - d01 * d20) / denominator;
        const T lambda[3] = {1 - l1 - l2, l1, l2};
        for (int q = 0; q < 3; ++q)
          if (lambda[q] > tolerance)
            at.push_back(F(f, q));
        if (at.size() == 3)
          at.clear();
      }
      if (at.empty())
        return f;
      const T side = (S(c) < 0) ? T(-1) : T(1);
      T best = -std::numeric_limits<T>::max();
      for (int g : elementsAt[at[0]]) {
        bool sharesAll = true;
        for (int v : at) {
          bool found = false;
          for (int q = 0; q < D; ++q)
            found = found || F(g, q) == v;
          sharesAll = sharesAll && found;
        }
        if (!sharesAll)
          continue;
        const T score = side * DotProduct(normalOf(g), dir) / distance;
        if (score > best) {
          best = score;
          f = g;
        }
      }
      return f;
    };

    // every node keeps its closest point on the mesh and the outward normal
    // of the closest element: the Hermite data of the surface
    struct Entry {
      viennahrle::Index<D> index;
      T value;
      Vec3D<T> closestPoint, normal;
    };
    std::vector<Entry> entries;
    entries.reserve(candidates.size());
    for (std::size_t c = 0; c < candidates.size(); ++c) {
      if (!std::isfinite(S(c)) || std::abs(S(c)) > halfWidth)
        continue;
      Entry entry{grid.globalIndices2LocalIndices(candidates[c]),
                  (S(c) == T(0)) ? T(0) : S(c), Vec3D<T>{}, Vec3D<T>{}};
      // a node wrapped across a periodic boundary keeps its closest point in
      // its own frame
      for (int i = 0; i < D; ++i)
        entry.closestPoint[i] =
            (C(c, i) + (grid.isBoundaryPeriodic(i)
                            ? T(entry.index[i] - candidates[c][i])
                            : T(0))) *
            gridDelta;
      entry.normal = normalOf(sharedElement(c));
      entries.push_back(entry);
    }

    // periodic images can map two nodes onto one, keep the smaller distance
    std::sort(entries.begin(), entries.end(), [](const auto &a, const auto &b) {
      if (a.index != b.index)
        return a.index < b.index;
      return std::abs(a.value) < std::abs(b.value);
    });
    entries.erase(std::unique(entries.begin(), entries.end(),
                              [](const auto &a, const auto &b) {
                                return a.index == b.index;
                              }),
                  entries.end());

    std::vector<std::pair<viennahrle::Index<D>, T>> points;
    typename PointData<T>::VectorDataType closestPoints, normals;
    points.reserve(entries.size());
    closestPoints.reserve(entries.size());
    normals.reserve(entries.size());
    for (const auto &entry : entries) {
      points.emplace_back(entry.index, entry.value);
      closestPoints.push_back(entry.closestPoint);
      normals.push_back(entry.normal);
    }

    // points are sorted, so their order is the point id order of the domain
    levelSet->insertPoints(points);
    levelSet->getDomain().segment();
    levelSet->finalize(width);
    auto &pointData = levelSet->getPointData();
    pointData.clear();
    pointData.insertNextVectorData(std::move(closestPoints),
                                   closestPointsLabel);
    pointData.insertNextVectorData(std::move(normals), surfaceNormalsLabel);
  }

public:
  /// Point data written by the Euclidean construction: the closest point on
  /// the mesh and the outward normal of the closest mesh element of every
  /// level set point.
  static constexpr char closestPointsLabel[] = "ClosestPoints";
  static constexpr char surfaceNormalsLabel[] = "SurfaceNormals";

  FromSurfaceMesh() = default;

  FromSurfaceMesh(SmartPointer<Domain<T, D>> passedLevelSet,
                  SmartPointer<Mesh<T>> passedMesh,
                  bool passedRemoveBoundaryTriangles = true)
      : levelSet(passedLevelSet), mesh(passedMesh),
        removeBoundaryTriangles{passedRemoveBoundaryTriangles,
                                passedRemoveBoundaryTriangles,
                                passedRemoveBoundaryTriangles} {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  void setMesh(SmartPointer<Mesh<T>> passedMesh) { mesh = passedMesh; }

  /// Sign of Euclidean distances: angle-weighted pseudo-normals (default,
  /// also for open meshes) or the generalized winding number (closed meshes,
  /// robust against non-manifold and self-intersecting parts).
  void setUseWindingNumber(bool passedUseWindingNumber) {
    useWindingNumber = passedUseWindingNumber;
  }

  /// Set whether all triangles outside the domain should be ignored (=true)
  /// or whether boundary conditions should be applied correctly to such
  /// triangles(=false). Defaults to true.
  void setRemoveBoundaryTriangles(bool passedRemoveBoundaryTriangles) {
    removeBoundaryTriangles.fill(passedRemoveBoundaryTriangles);
  }

  /// Set whether all triangles outside the domain should be ignored (=true)
  /// or whether boundary conditions should be applied correctly to such
  /// triangles(=false), for each direction. Defaults to true for all
  /// directions.
  template <std::size_t N>
  void setRemoveBoundaryTriangles(
      std::array<bool, N> passedRemoveBoundaryTriangles) {
    for (unsigned i = 0; i < D && i < N; ++i) {
      removeBoundaryTriangles[i] = passedRemoveBoundaryTriangles[i];
    }
  }

  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addError("No level set was passed to FromSurfaceMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      Logger::getInstance()
          .addError("No mesh was passed to FromSurfaceMesh.")
          .print();
      return;
    }

    const auto &grid = levelSet->getGrid();

    // caluclate which directions should apply removeBoundaryTriangles
    // it only makes sense to keep boundary triangles for periodic BNC
    bool removeBoundaries[D];
    for (unsigned i = 0; i < D; ++i) {
      if (!removeBoundaryTriangles[i] && grid.isBoundaryPeriodic(i)) {
        removeBoundaries[i] = false;
      } else {
        removeBoundaries[i] = true;
      }
    }

    if (levelSet->getNormalization() == NormalizationEnum::EUCLIDEAN) {
      applyEuclidean(removeBoundaries);
      return;
    }

    std::vector<std::pair<viennahrle::Index<D>, T>> points2;

    // setup list of grid points with distances to surface elements
    {
      typedef std::vector<std::pair<viennahrle::Index<D>, std::pair<T, T>>>
          point_vector;
      point_vector points;
      T gridDelta = grid.getGridDelta();

      VectorType<T, D> gridMin, gridMax;
      for (unsigned i = 0; i < D; ++i) {
        gridMin[i] = grid.getMinIndex(i) * gridDelta;
        gridMax[i] = grid.getMaxIndex(i) * gridDelta;
      }

      // for each surface element do
      std::vector<std::array<unsigned, D>> &elements =
          mesh->template getElements<D>();

      for (unsigned currentElement = 0; currentElement < elements.size();
           currentElement++) {
        VectorType<T, D> nodes[D]; // nodes of element
        VectorType<T, D> center;   // center point of triangle
        std::fill(center.begin(), center.end(), 0);

        std::bitset<2 * D> flags;
        flags.set();
        bool removeElement = false;

        for (int dim = 0; dim < D; dim++) {
          for (int q = 0; q < D; q++) {
            nodes[q][dim] = mesh->nodes[elements[currentElement][q]][dim];
            if (std::abs(nodes[q][dim] - gridMin[dim]) <
                boundaryEps * gridDelta)
              nodes[q][dim] = gridMin[dim];
            if (std::abs(nodes[q][dim] - gridMax[dim]) <
                boundaryEps * gridDelta)
              nodes[q][dim] = gridMax[dim];

            if (nodes[q][dim] > gridMin[dim])
              flags.reset(dim);
            if (nodes[q][dim] < gridMax[dim])
              flags.reset(dim + D);

            center[dim] += nodes[q][dim]; // center point calculation
          }
          if (removeBoundaries[dim] && (flags[dim] || flags[dim + D])) {
            removeElement = true;
          }
        }

        if (removeElement)
          continue;

        // center point calculation
        for (int q = 0; q < D; q++)
          center[q] /= static_cast<T>(D);

        // VectorType<T,D> normal=NormalVector(c);  //normal vector
        // calculation

        // determine the minimum and maximum nodes of the element, based on
        // their coordinates
        VectorType<T, D> minNode = nodes[0];
        VectorType<T, D> maxNode = nodes[0];
        for (int i = 1; i < D; ++i) {
          minNode = Min(minNode, nodes[i]);
          maxNode = Max(maxNode, nodes[i]);
        }

        // find indices which describe the triangle
        viennahrle::Index<D> minIndex, maxIndex;
        for (int q = 0; q < D; q++) {
          minIndex[q] =
              static_cast<hrleIndexType>(std::ceil(minNode[q] / gridDelta));
          maxIndex[q] =
              static_cast<hrleIndexType>(std::floor(maxNode[q] / gridDelta));
        }

        // Only iterate over indices which are actually inside the domain
        // if boundary conditions should be ignored
        for (unsigned i = 0; i < D; ++i) {
          if (removeBoundaries[i]) {
            minIndex[i] = std::max(minIndex[i], grid.getMinIndex(i));
            // for Periodic BNC, MaxGridPoint is MaxIndex-1
            maxIndex[i] = std::min(maxIndex[i], grid.getMaxGridPoint(i));
          }
        }

        // each cartesian direction
        for (int z = 0; z < D; z++) {
          VectorType<hrleIndexType, D - 1> min_bb, max_bb;

          for (int h = 0; h < D - 1; ++h) {
            min_bb[h] = minIndex[(z + h + 1) % D];
            max_bb[h] = maxIndex[(z + h + 1) % D];
          }

          box bb(min_bb, max_bb);

          if (bb.is_empty())
            continue;

          for (typename box::iterator it_bb(bb); !it_bb.is_finished();
               ++it_bb) {

            viennahrle::Index<D> it_b;
            for (int h = 0; h < D - 1; ++h)
              it_b[(z + h + 1) % D] = (*it_bb)[h];

            VectorType<T, D> p;
            for (int k = 1; k < D; k++)
              p[(k + z) % D] = grid.gridPositionOfGlobalIndex(
                  (k + z) % D, it_b[(k + z) % D]);

            T intersection;
            int intersection_status = calculateGridlineTriangleIntersection(
                p, nodes, z, intersection);

            if (intersection_status != 0) {
              // if there is an intersection
              assert(intersection >= minNode[z] &&
                     "Intersection should be inside Domain!");
              assert(intersection <= maxNode[z] &&
                     "Intersection should be inside Domain!");
              // intersection = std::max(intersection, minNode[z]);
              // intersection = std::min(intersection, maxNode[z]);

              if (removeBoundaries[z] &&
                  intersection > gridDelta * (grid.getMaxIndex(z)))
                continue;
              if (removeBoundaries[z] &&
                  intersection < grid.getMinLocalCoordinate(z))
                continue;

              T intersection2 =
                  grid.globalCoordinate2LocalIndex(z, intersection);

              auto floor = static_cast<hrleIndexType>(
                  std::floor(intersection2 - distanceEps));
              auto ceil = static_cast<hrleIndexType>(
                  std::ceil(intersection2 + distanceEps));

              floor = std::max(floor, minIndex[z] - 1);
              ceil = std::min(ceil, maxIndex[z] + 1);
              floor = std::max(floor, grid.getMinIndex(z));
              ceil = std::min(ceil, grid.getMaxIndex(z));

              if (!removeBoundaries[z]) {
                floor = grid.globalIndex2LocalIndex(z, floor);
                ceil = grid.globalIndex2LocalIndex(z, ceil);
              }

              VectorType<T, D> t = center;
              t[z] -= intersection;
              for (int k = 1; k < D; k++)
                t[(z + k) % D] -= p[(z + k) % D];
              Normalize(t);

              for (it_b[z] = floor; it_b[z] <= ceil; ++(it_b[z])) {

                T RealDistance = it_b[z] - intersection2;

                T SignDistance = RealDistance;

                SignDistance -= signEps * t[z];

                if (intersection_status < 0) {
                  RealDistance = -RealDistance;
                  SignDistance = -SignDistance;
                }

                if (RealDistance < 0.)
                  SignDistance = -SignDistance;

                RealDistance *= (1. - distanceEps * 1e-3);
                if (RealDistance > 1.)
                  RealDistance = 1.;
                if (RealDistance < -1.)
                  RealDistance = -1.;

                if (RealDistance == 0.)
                  RealDistance = 0.; // to avoid zeros with negative sign

                points.push_back(
                    std::make_pair(grid.globalIndices2LocalIndices(it_b),
                                   std::make_pair(SignDistance, RealDistance)));
              }
            }
          }
        }
      }

      std::sort(points.begin(), points.end()); // sort points lexicographically

      // setup list of index/distance pairs for level set initialization
      auto it_points = points.begin();
      while (it_points != points.end()) {
        viennahrle::Index<D> tmp = it_points->first;
        points2.push_back(
            std::make_pair(it_points->first, it_points->second.second));

        do {
          ++it_points;
        } while ((it_points != points.end()) && (it_points->first == tmp));
      }
    }

    levelSet->insertPoints(points2); // initialize level set function
    levelSet->getDomain().segment(); // distribute points evenly across threads
    levelSet->finalize(2);
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(FromSurfaceMesh)

} // namespace viennals
