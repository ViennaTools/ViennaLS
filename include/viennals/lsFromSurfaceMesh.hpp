#pragma once

#include <lsPreCompileMacros.hpp>

#include <hrleTypes.hpp>

#include <lsDomain.hpp>
#include <lsMesh.hpp>

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

public:
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
          .addWarning("No level set was passed to FromSurfaceMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      Logger::getInstance()
          .addWarning("No mesh was passed to FromSurfaceMesh.")
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
