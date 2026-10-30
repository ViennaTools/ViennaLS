#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <unordered_map>
#include <vector>

#include <hrleSparseIterator.hpp>

#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>

#include <vcVectorType.hpp>

namespace viennals {

using namespace viennacore;

/// Material of points from the level sets below a surface in a material
/// stack, lowest first, as in Advect: a point belongs to the lowest level set
/// that contains it strictly, and to the surface above all of them. A point
/// on an interface belongs to the material above. Points are in grid
/// spacings. Values of a Euclidean level set come from its Hermite data: the
/// node of the point's cell whose closest point lies nearest to the point,
/// and the distance along that closest point's outward direction, which keeps
/// the sign right next to corners, where the multilinear interpolation used
/// for other level sets misplaces the surface by up to a tenth of a cell.
template <class T, int D> class MaterialStack {
  using hrleIndex = viennahrle::Index<D>;
  struct IndexHash {
    std::size_t operator()(const hrleIndex &index) const {
      std::size_t h = 0;
      for (int k = 0; k < D; ++k)
        h = h * 1000003u ^ std::hash<long long>()(index[k]);
      return h;
    }
  };

  std::vector<SmartPointer<Domain<T, D>>> levelSets;
  // defined values of every level set
  std::vector<std::unordered_map<hrleIndex, T, IndexHash>> values;
  // Hermite data of every node of a Euclidean level set, grid spacings: its
  // closest point and the outward direction there; for the planes
  // AdvectHermite stores, a point on the node's plane and its normal
  struct Plane {
    Vec3D<T> point, direction;
  };
  std::vector<std::unordered_map<hrleIndex, Plane, IndexHash>> planes;
  // rounding of the precision, below the fixed tolerances with double
  static constexpr T roundOff = T(1024) * std::numeric_limits<T>::epsilon();
  std::vector<char> nodePlanes;

  // value of level set l at a grid node; just outside its band one more than
  // the nearest defined neighbour along an axis on the same side, which
  // continues a planar surface exactly where a cell it crosses has a corner
  // outside the band (a Manhattan band at a slope), further away only the
  // side matters
  T nodeValue(std::size_t l, const hrleIndex &index) const {
    auto found = values[l].find(index);
    if (found != values[l].end())
      return found->second;
    viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType> it(
        levelSets[l]->getDomain(), index);
    const T side = (it.getValue() < 0) ? T(-1) : T(1);
    const T far = T(levelSets[l]->getLevelSetWidth()) / 2 + 1;
    T nearest = far;
    const auto &grid = levelSets[l]->getGrid();
    for (int d = 0; d < D; ++d)
      for (int step : {-1, 1}) {
        hrleIndex q = index;
        q[d] += step;
        auto neighbour = values[l].find(grid.globalIndices2LocalIndices(q));
        if (neighbour != values[l].end() && neighbour->second * side > 0)
          nearest = std::min(nearest, std::abs(neighbour->second) + 1);
      }
    return side * nearest;
  }

public:
  MaterialStack() = default;

  explicit MaterialStack(
      std::vector<SmartPointer<Domain<T, D>>> passedLevelSets)
      : levelSets(std::move(passedLevelSets)) {
    values.assign(levelSets.size(), {});
    planes.assign(levelSets.size(), {});
    nodePlanes.assign(levelSets.size(), 0);
    for (std::size_t l = 0; l < levelSets.size(); ++l) {
      const auto &pointData = levelSets[l]->getPointData();
      const auto *closest = pointData.getVectorData(
          FromSurfaceMesh<T, D>::closestPointsLabel, true);
      const auto *normals = pointData.getVectorData(
          FromSurfaceMesh<T, D>::surfaceNormalsLabel, true);
      const auto points = levelSets[l]->getNumberOfPoints();
      const bool hermite =
          levelSets[l]->getNormalization() == NormalizationEnum::EUCLIDEAN &&
          closest != nullptr && normals != nullptr &&
          closest->size() == points && normals->size() == points;
      // the band state of AdvectHermite marks planes instead of closest points
      nodePlanes[l] = pointData.getScalarData("HermiteValue", true) != nullptr;
      const T gridDelta = levelSets[l]->getGrid().getGridDelta();
      for (viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>
               it(levelSets[l]->getDomain());
           !it.isFinished(); it.next()) {
        if (!it.isDefined())
          continue;
        values[l][it.getStartIndices()] = it.getValue();
        if (!hermite)
          continue;
        Plane plane;
        Vec3D<T> r{};
        T length = 0.;
        for (int i = 0; i < 3; ++i) {
          plane.point[i] = (*closest)[it.getPointId()][i] / gridDelta;
          plane.direction[i] = (*normals)[it.getPointId()][i];
          if (i < D)
            r[i] = it.getStartIndices()[i] - plane.point[i];
          length += r[i] * r[i];
        }
        length = std::sqrt(length);
        if (!nodePlanes[l] && length > std::max(T(1e-6), roundOff)) {
          const T side = (it.getValue() < 0) ? T(-1) : T(1);
          for (int i = 0; i < 3; ++i)
            plane.direction[i] = side * r[i] / length;
        }
        planes[l][it.getStartIndices()] = plane;
      }
    }
  }

  /// Number of level sets below the surface, also the material of the
  /// surface above them.
  std::size_t size() const { return levelSets.size(); }

  /// Hermite data of level set l at the grid node p: its closest point, or
  /// for stored planes the plane point, and the outward direction there;
  /// false without Hermite data at p.
  bool nodePlane(std::size_t l, const Vec3D<T> &p, Vec3D<T> &point,
                 Vec3D<T> &direction) const {
    if (planes[l].empty())
      return false;
    const auto &grid = levelSets[l]->getGrid();
    hrleIndex q;
    for (int i = 0; i < D; ++i)
      q[i] = static_cast<viennahrle::IndexType>(std::lround(p[i]));
    const auto local = grid.globalIndices2LocalIndices(q);
    auto found = planes[l].find(local);
    if (found == planes[l].end())
      return false;
    point = found->second.point;
    for (int i = 0; i < D; ++i)
      if (grid.isBoundaryPeriodic(i))
        point[i] += T(q[i] - local[i]);
    direction = found->second.direction;
    return true;
  }

  /// Value of level set l at p. With Hermite data, from the distinct face
  /// planes of the nodes of p's cell: the largest signed distance where they
  /// meet convexly (each face inside the other planes), the smallest where
  /// concavely, the plane itself where there is one; otherwise, as at saddle
  /// corners, the distance along the outward direction of the closest point
  /// nearest to p. Without Hermite data, multilinear.
  T value(std::size_t l, const Vec3D<T> &p) const {
    if (!planes[l].empty()) {
      const auto &grid = levelSets[l]->getGrid();
      hrleIndex base;
      for (int i = 0; i < D; ++i)
        base[i] = static_cast<viennahrle::IndexType>(std::floor(p[i]));
      auto along = [&](const Vec3D<T> &point, const Vec3D<T> &direction) {
        T result = 0.;
        for (int i = 0; i < D; ++i)
          result += (p[i] - point[i]) * direction[i];
        return result;
      };
      std::vector<Plane> faces;
      T nearest = std::numeric_limits<T>::max(), nearestValue = 0.;
      for (int c = 0; c < (1 << D); ++c) {
        hrleIndex q = base;
        for (int i = 0; i < D; ++i)
          q[i] += (c >> i) & 1;
        const auto local = grid.globalIndices2LocalIndices(q);
        auto found = planes[l].find(local);
        if (found == planes[l].end())
          continue;
        // a node across a periodic boundary holds its point in its own frame
        Plane plane = found->second;
        for (int i = 0; i < D; ++i)
          if (grid.isBoundaryPeriodic(i))
            plane.point[i] += T(q[i] - local[i]);
        T distance = 0., offset = 0., cross = 0.;
        for (int i = 0; i < D; ++i) {
          distance += (p[i] - plane.point[i]) * (p[i] - plane.point[i]);
          // the node's closest point along its normal: a face
          const T r = T(q[i]) - plane.point[i];
          offset += r * plane.direction[i];
          cross += r * r;
        }
        if (distance < nearest) {
          nearest = distance;
          nearestValue = along(plane.point, plane.direction);
        }
        const bool face = nodePlanes[l] || cross - offset * offset <=
                                               std::max(T(1e-12), roundOff) *
                                                   std::max(cross, T(1));
        if (!face)
          continue;
        bool seen = false;
        for (const auto &f : faces) {
          T dn = 0., dp = 0.;
          for (int i = 0; i < D; ++i) {
            dn += (f.direction[i] - plane.direction[i]) *
                  (f.direction[i] - plane.direction[i]);
            dp += (plane.point[i] - f.point[i]) * f.direction[i];
          }
          seen = seen || (dn < std::max(T(1e-12), roundOff * roundOff) &&
                          std::abs(dp) < std::max(T(1e-9), roundOff));
        }
        if (!seen)
          faces.push_back(plane);
      }
      if (faces.size() == 1)
        return along(faces[0].point, faces[0].direction);
      if (!faces.empty()) {
        bool convex = true, concave = true;
        for (std::size_t a = 0; a < faces.size(); ++a)
          for (std::size_t b = 0; b < faces.size(); ++b) {
            if (a == b)
              continue;
            T side = 0.;
            for (int i = 0; i < D; ++i)
              side += (faces[a].point[i] - faces[b].point[i]) *
                      faces[b].direction[i];
            convex = convex && side <= std::max(T(1e-9), roundOff);
            concave = concave && side >= -std::max(T(1e-9), roundOff);
          }
        if (convex != concave) {
          T result = convex ? -std::numeric_limits<T>::max()
                            : std::numeric_limits<T>::max();
          for (const auto &f : faces)
            result = convex ? std::max(result, along(f.point, f.direction))
                            : std::min(result, along(f.point, f.direction));
          return result;
        }
      }
      if (nearest < std::numeric_limits<T>::max())
        return nearestValue;
    }
    const auto &grid = levelSets[l]->getGrid();
    hrleIndex base;
    T fraction[3] = {0, 0, 0};
    for (int i = 0; i < D; ++i) {
      base[i] = static_cast<viennahrle::IndexType>(std::floor(p[i]));
      fraction[i] = p[i] - T(base[i]);
    }
    T result = 0.;
    for (int c = 0; c < (1 << D); ++c) {
      hrleIndex q = base;
      T weight = 1.;
      for (int i = 0; i < D; ++i) {
        const int bit = (c >> i) & 1;
        q[i] += bit;
        weight *= bit ? fraction[i] : 1 - fraction[i];
      }
      if (weight != 0)
        result += weight * nodeValue(l, grid.globalIndices2LocalIndices(q));
    }
    return result;
  }

  /// Whether level set l holds a value at the grid node p.
  bool defined(std::size_t l, const Vec3D<T> &p) const {
    hrleIndex q;
    for (int i = 0; i < D; ++i)
      q[i] = static_cast<viennahrle::IndexType>(std::lround(p[i]));
    return values[l].count(
               levelSets[l]->getGrid().globalIndices2LocalIndices(q)) > 0;
  }

  /// Material at p: the lowest level set that contains it strictly, or size()
  /// above all of them.
  int material(const Vec3D<T> &p) const {
    for (std::size_t l = 0; l < levelSets.size(); ++l)
      if (value(l, p) < T(0))
        return int(l);
    return int(levelSets.size());
  }
};

} // namespace viennals
