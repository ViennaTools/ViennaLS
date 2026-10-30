#pragma once

#include <hrleSparseIterator.hpp>

#include <lsBooleanOperation.hpp>
#include <lsConcepts.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromMesh.hpp>
#include <lsGeometricAdvectDistributions.hpp>
#include <lsPreCompileMacros.hpp>
#include <lsReinitialize.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToDualContouringMesh.hpp>

#include <vcLogger.hpp>
#include <vcSmartPointer.hpp>
#include <vcVectorType.hpp>

#include <igl/AABB.h>

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <set>

#ifndef NDEBUG // if in debug build
#include <lsCheck.hpp>
#include <lsToMesh.hpp>
#include <lsVTKWriter.hpp>
#endif

namespace viennals {

using namespace viennacore;

/// This class advects the level set according to a given distribution.
/// This distribution is overlayed at every grid point of the old surface. All
/// cells within this distribution are then filled, with cells at the edge
/// marked with the correct level set values. Therefore, the surface can be
/// shifted long distances in one step. This algorithm is therefore preferable
/// to normal advection if there is growth/reduction by a purely geometric
/// directional distribution.
template <class T, int D> class GeometricAdvect {
  using hrleIndexType = viennahrle::IndexType;
  using hrleCoordType = viennahrle::CoordType;

  SmartPointer<Domain<T, D>> levelSet = nullptr;
  SmartPointer<Domain<T, D>> maskLevelSet = nullptr;
  SmartPointer<GeometricAdvectDistribution<T, D>> dist = nullptr;
  static constexpr T cutoffValue =
      T(1.) + std::numeric_limits<T>::epsilon() * T(100);

  static void incrementIndices(viennahrle::Index<D> &indices,
                               const viennahrle::Index<D> &min,
                               const viennahrle::Index<D> &max) {
    int dim = 0;
    for (; dim < D - 1; ++dim) {
      if (indices[dim] < max[dim])
        break;
      indices[dim] = min[dim];
    }
    ++indices[dim];
  }

  template <class K, class V, template <class...> class MapType, class... Ts>
  MapType<V, K> inverseTranslator(MapType<K, V, Ts...> &map) {
    MapType<V, K> inv;
    std::for_each(map.begin(), map.end(), [&inv](const std::pair<K, V> &p) {
      inv.insert(std::make_pair(p.second, p.first));
    });
    return inv;
  }

  using Point = Eigen::Matrix<double, 1, D>;
  using Element = std::array<Point, D>;

  // whether the box of half extents h around x meets the segment (2D) or
  // triangle (3D) e, by the separating axes of the two
  static bool boxMeets(const Point &x, const std::array<double, D> &h,
                       const Element &e) {
    std::array<Point, D> v;
    for (int k = 0; k < D; ++k)
      v[k] = e[k] - x;
    for (int i = 0; i < D; ++i) {
      double lo = v[0](i), hi = v[0](i);
      for (int k = 1; k < D; ++k) {
        lo = std::min(lo, v[k](i));
        hi = std::max(hi, v[k](i));
      }
      if (lo > h[i] || hi < -h[i])
        return false;
    }
    if constexpr (D == 2) {
      const double nx = v[0](1) - v[1](1), ny = v[1](0) - v[0](0);
      return std::abs(nx * v[0](0) + ny * v[0](1)) <=
             h[0] * std::abs(nx) + h[1] * std::abs(ny);
    } else {
      const Eigen::Vector3d a = v[0].transpose(), b = v[1].transpose(),
                            c = v[2].transpose();
      const Eigen::Vector3d extent(h[0], h[1], h[2]);
      const Eigen::Vector3d edges[3] = {b - a, c - b, a - c};
      const Eigen::Vector3d n = edges[0].cross(edges[1]);
      if (std::abs(n.dot(a)) > extent.dot(n.cwiseAbs()))
        return false;
      for (int i = 0; i < 3; ++i)
        for (const auto &f : edges) {
          const Eigen::Vector3d axis = Eigen::Vector3d::Unit(i).cross(f);
          const double pa = axis.dot(a), pb = axis.dot(b), pc = axis.dot(c);
          const double r = extent.dot(axis.cwiseAbs());
          if (std::min({pa, pb, pc}) > r || std::max({pa, pb, pc}) < -r)
            return false;
        }
      return true;
    }
  }

  template <class Elements>
  static void toMatrices(const std::vector<Point> &points,
                         const Elements &elements, Eigen::MatrixXd &V,
                         Eigen::MatrixXi &F) {
    V.resize(points.size(), D);
    for (std::size_t p = 0; p < points.size(); ++p)
      V.row(p) = points[p];
    F.resize(elements.size(), D);
    for (std::size_t e = 0; e < elements.size(); ++e)
      for (int k = 0; k < D; ++k)
        F(e, k) = elements[e][k];
  }

  // A Euclidean level set with a sphere or box distribution. The new solid
  // is the old one joined with (deposition) or reduced by (etching) the
  // Minkowski sum of the distribution with the contributing surface, the
  // elements of the dual contouring mesh of the old surface that do not lie
  // on the mask. Every node of the new band takes the exact distance to that
  // sum or to the old surface, whichever bounds the new solid nearer, and
  // the closest point and normal of that surface as its plane wherever the
  // point lies on the new surface; Reinitialize builds the exact band from
  // the dual contouring mesh of these planes, so that the sharp edges of the
  // result stay sharp. Returns false for other distributions.
  bool applyEuclidean() {
    // the shape of the distribution, a ball or an axis-aligned box; signs as
    // in getBounds: deposition only if all its extents are positive
    bool ball = false, deposition = true;
    double radius = 0.;
    std::array<double, D> half{};
    if (auto *sphere = dynamic_cast<SphereDistribution<T, D> *>(dist.get())) {
      ball = true;
      radius = std::abs(double(sphere->radius));
      deposition = sphere->radius > 0;
    } else if (auto *box = dynamic_cast<BoxDistribution<T, D> *>(dist.get())) {
      for (int i = 0; i < D; ++i) {
        half[i] = std::abs(double(box->posExtent[i]));
        deposition = deposition && box->posExtent[i] > 0;
      }
    } else {
      return false;
    }
    std::array<double, D> reach;
    double largest = 0.;
    for (int i = 0; i < D; ++i) {
      reach[i] = ball ? radius : half[i];
      largest = std::max(largest, reach[i]);
    }
    if (largest == 0.)
      return true;

    auto &pointData = levelSet->getPointData();
    const auto *oldClosest = pointData.getVectorData(
        FromSurfaceMesh<T, D>::closestPointsLabel, true);
    const auto *oldNormals = pointData.getVectorData(
        FromSurfaceMesh<T, D>::surfaceNormalsLabel, true);
    if (oldClosest == nullptr || oldNormals == nullptr ||
        oldClosest->size() != levelSet->getNumberOfPoints() ||
        oldNormals->size() != levelSet->getNumberOfPoints())
      return false;

    const auto &grid = levelSet->getGrid();
    const double gridDelta = grid.getGridDelta();
    const int width = std::max(levelSet->getLevelSetWidth(), 2);
    // nodes up to this distance from the new surface enter the new band
    const double band = (0.5 * width + 1.) * gridDelta;
    const double tolerance = 1e-6 * gridDelta;
    constexpr double infinity = std::numeric_limits<double>::infinity();

    // the contributing surface
    std::vector<Element> elements;
    {
      auto surface = SmartPointer<Mesh<T>>::New();
      ToDualContouringMesh<T, D>(levelSet, surface).apply();
      std::vector<Element> all;
      auto add = [&](const auto &ids) {
        Element e;
        for (int k = 0; k < D; ++k)
          for (int i = 0; i < D; ++i)
            e[k](i) = surface->nodes[ids[k]][i];
        all.push_back(e);
      };
      if constexpr (D == 3)
        for (const auto &t : surface->triangles)
          add(t);
      else
        for (const auto &l : surface->lines)
          add(l);
      // an element whose centroid lies on the surface of the mask does not
      // contribute
      std::vector<char> onMask(all.size(), 0);
      if (maskLevelSet != nullptr) {
        auto maskSurface = SmartPointer<Mesh<T>>::New();
        ToDualContouringMesh<T, D>(maskLevelSet, maskSurface).apply();
        std::vector<Point> maskNodes;
        for (const auto &n : maskSurface->nodes) {
          Point p;
          for (int i = 0; i < D; ++i)
            p(i) = n[i];
          maskNodes.push_back(p);
        }
        Eigen::MatrixXd MV;
        Eigen::MatrixXi MF;
        if constexpr (D == 3)
          toMatrices(maskNodes, maskSurface->triangles, MV, MF);
        else
          toMatrices(maskNodes, maskSurface->lines, MV, MF);
        if (MF.rows() > 0) {
          igl::AABB<Eigen::MatrixXd, D> maskTree;
          maskTree.init(MV, MF);
          const double onSurface = 1e-3 * gridDelta;
#pragma omp parallel for
          for (long e = 0; e < long(all.size()); ++e) {
            Point centroid = Point::Zero();
            for (int k = 0; k < D; ++k)
              centroid += all[e][k] / double(D);
            int nearest = -1;
            Point c;
            maskTree.squared_distance(MV, MF, centroid, onSurface * onSurface,
                                      nearest, c);
            onMask[e] = nearest >= 0;
          }
        }
      }
      for (std::size_t e = 0; e < all.size(); ++e)
        if (!onMask[e])
          elements.push_back(all[e]);
    }
    // images of the elements near reflective and periodic boundaries, so
    // that the sum reaches across them
    for (int i = 0; i < D; ++i) {
      if (!grid.isBoundaryReflective(i) && !grid.isBoundaryPeriodic(i))
        continue;
      const double lo = grid.getMinGridPoint(i) * gridDelta;
      const double hi =
          (grid.getMaxGridPoint(i) + (grid.isBoundaryPeriodic(i) ? 1 : 0)) *
          gridDelta;
      const double near = reach[i] + band + gridDelta;
      const std::size_t count = elements.size();
      for (std::size_t e = 0; e < count; ++e) {
        double emin = elements[e][0](i), emax = emin;
        for (int k = 1; k < D; ++k) {
          emin = std::min(emin, elements[e][k](i));
          emax = std::max(emax, elements[e][k](i));
        }
        for (const bool low : {true, false}) {
          if (low ? emin > lo + near : emax < hi - near)
            continue;
          Element image = elements[e];
          for (auto &p : image)
            p(i) = grid.isBoundaryPeriodic(i) ? p(i) + (low ? hi - lo : lo - hi)
                                              : 2. * (low ? lo : hi) - p(i);
          elements.push_back(image);
        }
      }
    }
    if (elements.empty())
      return true;

    // outward unit normals of the elements, for the planes of their offsets
    std::vector<Point> elementNormals(elements.size(), Point::Zero());
    for (std::size_t e = 0; e < elements.size(); ++e) {
      Point n;
      if constexpr (D == 2) {
        n(0) = elements[e][0](1) - elements[e][1](1);
        n(1) = elements[e][1](0) - elements[e][0](0);
      } else {
        n = (elements[e][1] - elements[e][0])
                .cross(elements[e][2] - elements[e][0]);
      }
      if (n.norm() > 0.)
        elementNormals[e] = n / n.norm();
    }

    // the sum: for a ball the distance to the contributing surface; for a
    // box the facets of the sum of each element with the box, which contain
    // the boundary of the whole sum, and the elements by the cells of a
    // bucket grid their bounding boxes, widened by the box, overlap
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::AABB<Eigen::MatrixXd, D> tree;
    std::array<double, D> bucketOrigin{}, bucketSize{};
    std::array<long, D> bucketCount{};
    std::vector<std::vector<unsigned>> buckets;
    // the elements whose bounding boxes, widened by the box, contain x
    const std::vector<unsigned> noElements;
    auto bucketOf = [&](const Point &x) -> const std::vector<unsigned> & {
      long id = 0;
      for (int i = D - 1; i >= 0; --i) {
        const long c =
            long(std::floor((x(i) - bucketOrigin[i]) / bucketSize[i]));
        if (c < 0 || c >= bucketCount[i])
          return noElements;
        id = id * bucketCount[i] + c;
      }
      return buckets[id];
    };
    // whether x lies in the sum, for a box of half extents h
    auto inBox = [&](const Point &x, const std::array<double, D> &h) {
      for (const unsigned e : bucketOf(x))
        if (boxMeets(x, h, elements[e]))
          return true;
      return false;
    };
    // the projection of x onto the intersection of the half-spaces beyond
    // the offsets, by the support of the box, of the elements whose sums
    // contain x: inward for an etching, outward for a deposition. The normal
    // is the outward normal of the offset most active at the projection, a
    // face of the result also where the projection lies on an edge or corner
    auto projectBeyondOffsets = [&](const Point &x, bool outward,
                                    const std::array<double, D> &h,
                                    const std::vector<unsigned> &candidates,
                                    Point &projection, Point &faceNormal) {
      std::vector<std::pair<Point, double>> planes; // a.p <= b
      for (const unsigned e : candidates) {
        if (!boxMeets(x, h, elements[e]))
          continue;
        const Point &n = elementNormals[e];
        double support = 0.;
        for (int i = 0; i < D; ++i)
          support += h[i] * std::abs(n(i));
        const double offset = n.dot(elements[e][0]);
        const Point a = outward ? Point(-n) : n;
        const double b = outward ? -(offset + support) : offset - support;
        bool seen = false;
        for (const auto &[a2, b2] : planes)
          seen = seen || ((a - a2).norm() < 1e-9 &&
                          std::abs(b - b2) < 1e-9 * (gridDelta + std::abs(b)));
        if (!seen)
          planes.emplace_back(a, b);
      }
      const int count = int(planes.size());
      if (count == 0 || count > 12)
        return false;
      double best = infinity;
      for (int mask = 1; mask < (1 << count); ++mask) {
        int active[D], size = 0;
        bool small = true;
        for (int k = 0; k < count; ++k)
          if ((mask >> k) & 1) {
            if (size == D) {
              small = false;
              break;
            }
            active[size++] = k;
          }
        if (!small)
          continue;
        // nearest point of the intersection of the active planes
        Eigen::MatrixXd A(size, D);
        Eigen::VectorXd r(size);
        for (int k = 0; k < size; ++k) {
          A.row(k) = planes[active[k]].first;
          r(k) = planes[active[k]].first.dot(x) - planes[active[k]].second;
        }
        const Eigen::MatrixXd G = A * A.transpose();
        if (std::abs(G.determinant()) < 1e-12)
          continue;
        const Eigen::VectorXd lambda = G.ldlt().solve(r);
        if ((lambda.array() < -1e-12).any())
          continue;
        const Point p = x - (A.transpose() * lambda).transpose();
        bool feasible = true;
        for (const auto &[a, b] : planes)
          feasible = feasible && a.dot(p) <= b + tolerance;
        if (feasible && (p - x).norm() < best) {
          best = (p - x).norm();
          projection = p;
          int strongest = 0;
          for (int k = 1; k < size; ++k)
            if (lambda(k) > lambda(strongest))
              strongest = k;
          const Point &a = planes[active[strongest]].first;
          faceNormal = outward ? Point(-a) : a;
        }
      }
      return best < infinity;
    };
    if (ball) {
      std::vector<Point> points;
      std::vector<std::array<unsigned, D>> simplices;
      for (const auto &e : elements) {
        std::array<unsigned, D> s;
        for (int k = 0; k < D; ++k) {
          s[k] = points.size();
          points.push_back(e[k]);
        }
        simplices.push_back(s);
      }
      toMatrices(points, simplices, V, F);
    } else {
      std::array<double, D> lo, hi;
      lo.fill(infinity);
      hi.fill(-infinity);
      for (const auto &e : elements)
        for (const auto &p : e)
          for (int i = 0; i < D; ++i) {
            lo[i] = std::min(lo[i], p(i) - half[i]);
            hi[i] = std::max(hi[i], p(i) + half[i]);
          }
      long total = 1;
      for (int i = 0; i < D; ++i) {
        bucketOrigin[i] = lo[i] - gridDelta;
        bucketSize[i] = std::max(2. * half[i], 2. * gridDelta);
        bucketCount[i] =
            long((hi[i] + gridDelta - bucketOrigin[i]) / bucketSize[i]) + 1;
        total *= bucketCount[i];
      }
      buckets.resize(total);
      for (unsigned e = 0; e < elements.size(); ++e) {
        std::array<long, D> from, to;
        for (int i = 0; i < D; ++i) {
          double emin = elements[e][0](i), emax = emin;
          for (int k = 1; k < D; ++k) {
            emin = std::min(emin, elements[e][k](i));
            emax = std::max(emax, elements[e][k](i));
          }
          from[i] = long((emin - half[i] - tolerance - bucketOrigin[i]) /
                         bucketSize[i]);
          to[i] = long((emax + half[i] + tolerance - bucketOrigin[i]) /
                       bucketSize[i]);
        }
        std::array<long, D> c = from;
        while (true) {
          long id = 0;
          for (int i = D - 1; i >= 0; --i)
            id = id * bucketCount[i] + c[i];
          buckets[id].push_back(e);
          int i = 0;
          for (; i < D; ++i) {
            if (c[i] < to[i]) {
              ++c[i];
              break;
            }
            c[i] = from[i];
          }
          if (i == D)
            break;
        }
      }
      constexpr int corners = 1 << D;
      std::array<Point, corners> corner;
      for (int m = 0; m < corners; ++m)
        for (int i = 0; i < D; ++i)
          corner[m](i) = ((m >> i) & 1) ? half[i] : -half[i];
      std::map<std::array<double, D>, unsigned> ids;
      std::vector<Point> unique;
      std::vector<std::array<unsigned, D>> elementIds;
      for (const auto &e : elements) {
        std::array<unsigned, D> s;
        for (int k = 0; k < D; ++k) {
          std::array<double, D> key;
          for (int i = 0; i < D; ++i)
            key[i] = e[k](i);
          const auto [it, inserted] = ids.try_emplace(key, unique.size());
          if (inserted)
            unique.push_back(e[k]);
          s[k] = it->second;
        }
        elementIds.push_back(s);
      }
      std::vector<Point> points;
      for (const auto &p : unique)
        for (int m = 0; m < corners; ++m)
          points.push_back(p + corner[m]);
      std::vector<std::array<unsigned, D>> simplices;
      auto area = [&](const std::array<unsigned, D> &s) {
        if constexpr (D == 2)
          return (points[s[1]] - points[s[0]]).norm();
        else
          return (points[s[1]] - points[s[0]])
              .cross(points[s[2]] - points[s[0]])
              .norm();
      };
      const double degenerate =
          (D == 2) ? 1e-9 * gridDelta : 1e-18 * gridDelta * gridDelta;
      auto addSimplex = [&](const std::array<unsigned, D> &s) {
        if (area(s) > degenerate)
          simplices.push_back(s);
      };
      // each element at every corner of the box
      for (const auto &s : elementIds)
        for (int m = 0; m < corners; ++m) {
          std::array<unsigned, D> c;
          for (int k = 0; k < D; ++k)
            c[k] = s[k] * corners + m;
          addSimplex(c);
        }
      // edges of the box at every vertex (2D), parallelograms of every edge
      // with every edge of the box and faces of the box at every vertex (3D)
      std::vector<std::pair<int, int>> boxEdges;
      for (int m = 0; m < corners; ++m)
        for (int i = 0; i < D; ++i)
          if (!((m >> i) & 1))
            boxEdges.emplace_back(m, m | (1 << i));
      if constexpr (D == 2) {
        for (unsigned p = 0; p < unique.size(); ++p)
          for (const auto &[a, b] : boxEdges)
            addSimplex({p * corners + a, p * corners + b});
      } else {
        std::set<std::pair<unsigned, unsigned>> edges;
        for (const auto &s : elementIds)
          for (int k = 0; k < 3; ++k)
            edges.emplace(std::min(s[k], s[(k + 1) % 3]),
                          std::max(s[k], s[(k + 1) % 3]));
        for (const auto &[p, q] : edges)
          for (const auto &[a, b] : boxEdges) {
            addSimplex({p * corners + a, q * corners + a, q * corners + b});
            addSimplex({p * corners + a, q * corners + b, p * corners + b});
          }
        for (unsigned p = 0; p < unique.size(); ++p)
          for (int i = 0; i < 3; ++i)
            for (int side = 0; side < 2; ++side) {
              const int j = (i + 1) % 3, l = (i + 2) % 3, base = side << i;
              const unsigned c0 = p * corners + base,
                             c1 = p * corners + (base | (1 << j)),
                             c2 = p * corners + (base | (1 << j) | (1 << l)),
                             c3 = p * corners + (base | (1 << l));
              addSimplex({c0, c1, c2});
              addSimplex({c0, c2, c3});
            }
      }
      // pieces inside the sum are no part of its boundary: they would give
      // nodes inside the sum a nearest point that is not on the new surface
      std::array<double, D> inner;
      for (int i = 0; i < D; ++i)
        inner[i] = std::max(half[i] - tolerance, 0.);
      std::vector<char> interior(simplices.size(), 0);
#pragma omp parallel for schedule(dynamic, 256)
      for (long q = 0; q < long(simplices.size()); ++q) {
        const auto &sx = simplices[q];
        Point centroid = Point::Zero();
        for (int k = 0; k < D; ++k)
          centroid += points[sx[k]] / double(D);
        bool inside = inBox(centroid, inner);
        for (int k = 0; inside && k < D; ++k) {
          const Point &p = points[sx[k]], &n = points[sx[(k + 1) % D]];
          inside =
              inBox(p + 1e-3 * (centroid - p), inner) &&
              inBox(0.5 * (p + n) + 1e-3 * (centroid - 0.5 * (p + n)), inner);
        }
        interior[q] = inside;
      }
      std::vector<std::array<unsigned, D>> boundary;
      for (std::size_t q = 0; q < simplices.size(); ++q)
        if (!interior[q])
          boundary.push_back(simplices[q]);
      toMatrices(points, boundary, V, F);
    }
    if (F.rows() == 0)
      return true;
    tree.init(V, F);

    // whether a point of the old surface lies within the sum, not on its
    // boundary
    auto covered = [&](const Point &c) {
      if (ball) {
        if (radius <= tolerance)
          return false;
        int nearest = -1;
        Point y;
        tree.squared_distance(
            V, F, c, (radius - tolerance) * (radius - tolerance), nearest, y);
        return nearest >= 0;
      }
      std::array<double, D> h;
      for (int i = 0; i < D; ++i)
        h[i] = std::max(half[i] - tolerance, 0.);
      return inBox(c, h);
    };

    // nodes near the sum: the bounding box of the elements widened by the
    // reach of the distribution and the band, within reflective and periodic
    // boundaries
    viennahrle::Index<D> regionMin, regionMax;
    {
      std::array<double, D> lo, hi;
      lo.fill(infinity);
      hi.fill(-infinity);
      for (const auto &e : elements)
        for (const auto &p : e)
          for (int i = 0; i < D; ++i) {
            lo[i] = std::min(lo[i], p(i));
            hi[i] = std::max(hi[i], p(i));
          }
      for (int i = 0; i < D; ++i) {
        regionMin[i] =
            hrleIndexType(std::floor((lo[i] - reach[i] - band) / gridDelta));
        regionMax[i] =
            hrleIndexType(std::ceil((hi[i] + reach[i] + band) / gridDelta));
        if (!grid.isNegBoundaryInfinite(i))
          regionMin[i] = std::max(regionMin[i], grid.getMinGridPoint(i));
        if (!grid.isPosBoundaryInfinite(i))
          regionMax[i] = std::min(regionMax[i], grid.getMaxGridPoint(i));
      }
    }
    auto inRegion = [&](const viennahrle::Index<D> &index) {
      for (int i = 0; i < D; ++i)
        if (index[i] < regionMin[i] || index[i] > regionMax[i])
          return false;
      return true;
    };

    struct Node {
      viennahrle::Index<D> index;
      double value; // old signed distance, infinite outside the old band
      Point closest, normal;
      bool defined;
      // the mask, which an etching leaves
      double maskValue;
      Point maskClosest, maskNormal;
      bool maskPlane;
    };
    struct Record {
      viennahrle::Index<D> index;
      T value;
      Vec3D<T> closest, normal;
      T priority;
    };
    auto toVec = [](const Point &p) {
      Vec3D<T> v{T(0), T(0), T(0)};
      for (int i = 0; i < D; ++i)
        v[i] = T(p(i));
      return v;
    };
    auto fromVec = [](const Vec3D<T> &v) {
      Point p;
      for (int i = 0; i < D; ++i)
        p(i) = v[i];
      return p;
    };

    // old band nodes away from the sum keep value and plane; the nodes of
    // the region that can lie in the new band are evaluated
    std::vector<Record> records;
    std::vector<Node> nodes;
    using DomainType = typename Domain<T, D>::DomainType;
    for (viennahrle::ConstSparseIterator<DomainType> it(levelSet->getDomain());
         !it.isFinished(); ++it) {
      if (!it.isDefined() || inRegion(it.getStartIndices()))
        continue;
      const auto id = it.getPointId();
      records.push_back({it.getStartIndices(), it.getValue(), (*oldClosest)[id],
                         (*oldNormals)[id], T(1)});
    }
    // an etching leaves the mask: each node joins it before the band is
    // rebuilt, so that the surface between the mask and the etched region
    // takes the planes of both
    const bool keepMask = !deposition && maskLevelSet != nullptr;
    const typename PointData<T>::VectorDataType *maskClosest = nullptr,
                                                *maskNormals = nullptr;
    if (keepMask) {
      auto &maskData = maskLevelSet->getPointData();
      maskClosest = maskData.getVectorData(
          FromSurfaceMesh<T, D>::closestPointsLabel, true);
      maskNormals = maskData.getVectorData(
          FromSurfaceMesh<T, D>::surfaceNormalsLabel, true);
      if (maskClosest == nullptr || maskNormals == nullptr ||
          maskClosest->size() != maskLevelSet->getNumberOfPoints() ||
          maskNormals->size() != maskLevelSet->getNumberOfPoints()) {
        maskClosest = nullptr;
        maskNormals = nullptr;
      }
    }
    {
      viennahrle::ConstSparseIterator<DomainType> it(levelSet->getDomain());
      std::unique_ptr<viennahrle::ConstSparseIterator<DomainType>> maskIt;
      if (keepMask)
        maskIt = std::make_unique<viennahrle::ConstSparseIterator<DomainType>>(
            maskLevelSet->getDomain());
      for (viennahrle::Index<D> index = regionMin; index <= regionMax;
           incrementIndices(index, regionMin, regionMax)) {
        it.goToIndicesSequential(index);
        const bool defined = it.isDefined();
        double value = it.getValue() * gridDelta;
        if (!defined)
          value = (value < 0.) ? -infinity : infinity;
        // far inside the old solid in a deposition, or far outside it in an
        // etching, a node stays far from the new surface
        if (deposition ? value < -band : value > band)
          continue;
        Node node{index,         value,         Point::Zero(),
                  Point::Zero(), defined,       infinity,
                  Point::Zero(), Point::Zero(), false};
        if (defined) {
          node.closest = fromVec((*oldClosest)[it.getPointId()]);
          node.normal = fromVec((*oldNormals)[it.getPointId()]);
        }
        if (keepMask) {
          maskIt->goToIndicesSequential(index);
          node.maskValue = maskIt->getValue() * gridDelta;
          if (!maskIt->isDefined()) {
            node.maskValue = (node.maskValue < 0.) ? -infinity : infinity;
          } else if (maskClosest != nullptr) {
            node.maskClosest = fromVec((*maskClosest)[maskIt->getPointId()]);
            node.maskNormal = fromVec((*maskNormals)[maskIt->getPointId()]);
            node.maskPlane = true;
          }
        }
        nodes.push_back(node);
      }
    }

    std::vector<Record> evaluated(nodes.size());
    std::vector<char> kept(nodes.size(), 0);
#pragma omp parallel for schedule(dynamic, 64)
    for (long k = 0; k < long(nodes.size()); ++k) {
      const Node &node = nodes[k];
      Point x;
      for (int i = 0; i < D; ++i)
        x(i) = node.index[i] * gridDelta;

      // signed distance to the sum, exact outside it and a lower bound of
      // the depth inside, with the closest point and the normal pointing out
      // of the sum; inside, the point counts only where it lies on the
      // boundary of the sum
      bool inside, onBoundary = false;
      double sum = infinity;
      Point sumClosest = x, sumNormal = Point::Zero();
      int nearest = -1;
      Point y;
      if (ball) {
        const double sqrd = tree.squared_distance(
            V, F, x, (radius + band) * (radius + band), nearest, y);
        if (nearest >= 0) {
          const double d = std::sqrt(sqrd);
          inside = d <= radius;
          sum = d - radius;
          if (d > 0.) {
            sumNormal = (x - y) / d;
            sumClosest = y + radius * sumNormal;
            onBoundary = !inside;
            if (inside) {
              int other = -1;
              Point z;
              tree.squared_distance(V, F, sumClosest,
                                    (radius - tolerance) * (radius - tolerance),
                                    other, z);
              onBoundary = other < 0;
              // next to a sharp edge of the result the point lies just
              // beyond it, on the offset of a face that still bounds the
              // result there; the plane of that offset holds
              if (!onBoundary && std::abs(sumNormal.dot(
                                     elementNormals[nearest])) > 1. - 1e-9) {
                const double within = std::max(radius - gridDelta, 0.);
                other = -1;
                if (within > 0.)
                  tree.squared_distance(V, F, sumClosest, within * within,
                                        other, z);
                onBoundary = other < 0;
              }
            }
          }
        } else {
          inside = false;
        }
      } else {
        inside = inBox(x, half);
        const double sqrd =
            tree.squared_distance(V, F, x, band * band, nearest, y);
        const double d = (nearest >= 0) ? std::sqrt(sqrd) : infinity;
        sum = inside ? -d : d;
        if (nearest >= 0 && d > 0.) {
          sumNormal = (inside ? y - x : x - y) / d;
          sumClosest = y;
          onBoundary = !inside || !inBox(y + tolerance * sumNormal, half);
          // the plane of the nearest piece, out of the sum: at an edge or
          // corner of the sum a face of it, as FromSurfaceMesh gives, where
          // the direction to the node could lie within the singular value
          // threshold of a face and blur the edge
          Point m;
          if constexpr (D == 2) {
            m(0) = V(F(nearest, 0), 1) - V(F(nearest, 1), 1);
            m(1) = V(F(nearest, 1), 0) - V(F(nearest, 0), 0);
          } else {
            const Point a = V.row(F(nearest, 0)), b = V.row(F(nearest, 1)),
                        c = V.row(F(nearest, 2));
            m = (b - a).cross(c - a);
          }
          if (m.norm() > 0.) {
            m /= m.norm();
            sumNormal = (m.dot(sumNormal) < 0.) ? Point(-m) : m;
          }
          // next to a sharp edge of the result the nearest piece runs on
          // beyond it into the sum. There the side of the boundary away
          // from the node is locally the intersection of the half-spaces
          // beyond the offsets of the elements whose sums contain the node:
          // the projection onto it is the nearest point of the boundary,
          // where it lies on the boundary
          if (!onBoundary) {
            Point p, n;
            if (projectBeyondOffsets(x, deposition, half, bucketOf(x), p, n)) {
              const double length = (p - x).norm();
              if (length > tolerance) {
                const Point u = (p - x) / length;
                if (!inBox(p + tolerance * u, half) &&
                    inBox(p - tolerance * u, half)) {
                  // out of the sum, so that the node takes the face normal
                  sumNormal = deposition ? n : Point(-n);
                  sumClosest = p;
                  onBoundary = true;
                }
              }
            }
          }
        }
      }

      // the new solid and the plane of the node, where its point lies on
      // the new surface
      const double old = node.value;
      double value;
      Point closest = x, normal = Point::Zero();
      bool plane = false;
      auto oldPlane = [&]() {
        if (node.defined && !covered(node.closest)) {
          closest = node.closest;
          normal = node.normal;
          plane = true;
        }
      };
      if (deposition) {
        if (old > 0.) {
          if (inside) {
            value = std::min(sum, -tolerance);
            closest = sumClosest;
            normal = sumNormal;
            plane = onBoundary;
          } else if (sum < old) {
            value = sum;
            closest = sumClosest;
            normal = sumNormal;
            plane = sum > tolerance;
          } else {
            value = old;
            oldPlane();
          }
        } else {
          value = std::min(old, inside ? std::min(sum, -tolerance) : old);
          oldPlane();
        }
      } else {
        if (old < 0.) {
          if (inside) {
            value = std::max(-sum, tolerance);
            closest = sumClosest;
            normal = -sumNormal;
            plane = onBoundary;
          } else if (sum < -old) {
            value = -sum;
            closest = sumClosest;
            normal = -sumNormal;
            plane = sum > tolerance;
          } else {
            value = old;
            oldPlane();
          }
        } else {
          value = old;
          oldPlane();
        }
      }
      if (keepMask && node.maskValue < value) {
        value = node.maskValue;
        closest = node.maskClosest;
        normal = node.maskNormal;
        plane = node.maskPlane;
      }
      if (!(std::abs(value) <= band))
        continue;
      evaluated[k] = {node.index, T(value / gridDelta), toVec(closest),
                      toVec(normal), T(plane ? 1 : 0)};
      kept[k] = 1;
    }
    for (std::size_t k = 0; k < nodes.size(); ++k)
      if (kept[k])
        records.push_back(evaluated[k]);
    std::sort(
        records.begin(), records.end(),
        [](const Record &a, const Record &b) { return a.index < b.index; });

    typename Domain<T, D>::PointValueVectorType points;
    typename PointData<T>::VectorDataType closestPoints, normals;
    typename PointData<T>::ScalarDataType priorities;
    for (const auto &r : records) {
      points.emplace_back(r.index, r.value);
      closestPoints.push_back(r.closest);
      normals.push_back(r.normal);
      priorities.push_back(r.priority);
    }
    // points are sorted, so their order is the point id order of the domain
    levelSet->insertPoints(points);
    levelSet->getDomain().segment();
    levelSet->finalize(width);
    pointData.clear();
    pointData.insertNextVectorData(std::move(closestPoints),
                                   FromSurfaceMesh<T, D>::closestPointsLabel);
    pointData.insertNextVectorData(std::move(normals),
                                   FromSurfaceMesh<T, D>::surfaceNormalsLabel);
    pointData.insertNextScalarData(std::move(priorities),
                                   ToDualContouringMesh<T, D>::priorityLabel);
    Reinitialize<T, D>(levelSet, width).apply();
    return true;
  }

public:
  GeometricAdvect() = default;

  GeometricAdvect(SmartPointer<Domain<T, D>> passedLevelSet,
                  SmartPointer<GeometricAdvectDistribution<T, D>> passedDist,
                  SmartPointer<Domain<T, D>> passedMaskLevelSet = nullptr)
      : levelSet(passedLevelSet), maskLevelSet(passedMaskLevelSet),
        dist(passedDist) {}

  /// Set the levelset which should be advected.
  void setLevelSet(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  /// Set which advection distribution to use. Must be derived from
  /// GeometricAdvectDistribution.
  void setAdvectionDistribution(
      SmartPointer<GeometricAdvectDistribution<T, D>> passedDist) {
    dist = passedDist;
  }

  /// Set the levelset, which should be used as a mask. This level set
  /// has to be wrapped by the levelset set by setLevelSet, so the mask
  /// is entirely inside the advected level set.
  void setMaskLevelSet(SmartPointer<Domain<T, D>> passedMaskLevelSet) {
    maskLevelSet = passedMaskLevelSet;
  }

  /// Perform geometrical advection.
  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addError("No level set passed to GeometricAdvect. Not Advecting.")
          .print();
      return;
    }
    if (dist == nullptr) {
      Logger::getInstance()
          .addError("No GeometricAdvectDistribution passed to "
                    "GeometricAdvect. Not "
                    "Advecting.")
          .print();
      return;
    }

    // a Euclidean level set keeps its sharp edges with sphere and box
    // distributions; with others it is rebuilt from the values below
    const bool euclidean =
        levelSet->getNormalization() == NormalizationEnum::EUCLIDEAN;
    const int euclideanWidth = levelSet->getLevelSetWidth();
    if (euclidean) {
      if (applyEuclidean())
        return;
      VIENNACORE_LOG_WARNING(
          "GeometricAdvect: this distribution is applied to the values of "
          "the Euclidean level set, which rounds its sharp edges; sphere and "
          "box distributions keep them.");
    }

    // levelSet must have at least a width of 3
    Expand<T, D>(levelSet, 3).apply();

    if (maskLevelSet != nullptr) {
      Expand<T, D>(maskLevelSet, 3).apply();
    }

    dist->prepare(levelSet);

    typedef typename Domain<T, D>::DomainType DomainType;

    auto &domain = levelSet->getDomain();

    auto &grid = levelSet->getGrid();
    const auto gridDelta = grid.getGridDelta();
    const bool useSurfacePointId = dist->useSurfacePointId();

    // Extract the original surface as a point cloud of grid
    // points shifted to the surface (disk mesh)
    auto surfaceMesh = SmartPointer<Mesh<hrleCoordType>>::New();
    auto pointIdTranslator =
        SmartPointer<typename ToDiskMesh<T, D>::TranslatorType>::New();
    ToDiskMesh<T, D, hrleCoordType>(levelSet, surfaceMesh, pointIdTranslator)
        .apply();
    if (!useSurfacePointId)
      *pointIdTranslator = inverseTranslator(*pointIdTranslator);

    // find bounds of distribution
    auto distBounds = dist->getBounds();

    // TODO: need to add support for periodic boundary conditions!
    viennahrle::Index<D> distMin, distMax;

    bool minPointNegative = domain.getDomainSegment(0).definedValues[0] < 0.;
    bool maxPointNegative =
        domain.getDomainSegment(domain.getNumberOfSegments() - 1)
            .definedValues.back() < 0.;
    bool distIsPositive = true;

    // find bounding box of old domain
    hrleIndexType bounds[6];
    domain.getDomainBounds(bounds);
    viennahrle::Index<D> min, max;
    for (unsigned i = 0; i < D; ++i) {
      // translate from coords to indices
      distMin[i] =
          distBounds[2 * i] / gridDelta + ((distBounds[2 * i] < 0) ? -2 : 2);
      distMax[i] = distBounds[2 * i + 1] / gridDelta +
                   ((distBounds[2 * i + 1] < 0) ? -2 : 2);
      if (distBounds[2 * i] >= 0) {
        distIsPositive = false;
      }

      // use the extent of the diskMesh to identify bounding box of new
      // level set
      // TODO: respect periodic boundary condition
      min[i] = surfaceMesh->minimumExtent[i] / gridDelta;
      // TODO also do the same thing for positive point and etching
      if (grid.isNegBoundaryInfinite(i) && minPointNegative && distMin[i] < 0) {
        min[i] -= 2;
      } else {
        if (distIsPositive) {
          min[i] += distMin[i];
        } else {
          min[i] -= distMin[i];
        }
      }
      // if calculated index is out of bounds, set the extent
      // TODO: need to add periodic BNC handling here
      if (min[i] < grid.getMinGridPoint(i)) {
        min[i] = grid.getMinGridPoint(i);
      }

      max[i] = surfaceMesh->maximumExtent[i] / gridDelta;
      if (grid.isPosBoundaryInfinite(i) && maxPointNegative && distMax[i] > 0) {
        max[i] += 2;
      } else {
        if (distIsPositive) {
          max[i] += distMax[i];
        } else {
          max[i] -= distMax[i];
        }
      }
      if (max[i] > grid.getMaxGridPoint(i)) {
        max[i] = grid.getMaxGridPoint(i);
      }
    }

    // Remove contribute points if they are part of the mask
    // If a mask is supplied, remove all contribute points which
    // lie on (or inside) the mask
    if (maskLevelSet != nullptr) {
      // Go over all contribute points and see if they are on the mask surface
      auto &maskDomain = maskLevelSet->getDomain();
      auto values = surfaceMesh->cellData.getScalarData("LSValues");
      auto valueIt = values->begin();

      auto newSurfaceMesh = SmartPointer<Mesh<hrleCoordType>>::New();
      PointData<hrleCoordType>::ScalarDataType newValues;
      viennahrle::ConstSparseIterator<DomainType> maskIt(maskDomain);
      for (auto &node : surfaceMesh->getNodes()) {
        viennahrle::Index<D> index;
        for (unsigned i = 0; i < D; ++i) {
          index[i] = std::round(node[i] / gridDelta);
        }
        // can do sequential, because surfaceNodes are lexicographically sorted
        // from lsToDiskMesh
        maskIt.goToIndicesSequential(index);
        // if it is a mask point, mark it to maybe use it in new level set
        if (!maskIt.isDefined() || !(maskIt.getValue() < *valueIt + 1e-5)) {
          newSurfaceMesh->insertNextNode(node);
          newValues.push_back(*valueIt);
          // insert vertex
          std::array<unsigned, 1> vertex{};
          vertex[0] = newSurfaceMesh->nodes.size();
          newSurfaceMesh->insertNextVertex(vertex);
        }
        ++valueIt;
      }
      newSurfaceMesh->cellData.insertNextScalarData(newValues, "LSValues");
      // use new mesh as surfaceMesh
      newSurfaceMesh->minimumExtent = surfaceMesh->minimumExtent;
      newSurfaceMesh->maximumExtent = surfaceMesh->maximumExtent;
      surfaceMesh = newSurfaceMesh;
    }

#ifndef NDEBUG // if in debug build
    {
      VIENNACORE_LOG_DEBUG("GeomAdvect: Writing debug meshes");
      VTKWriter<hrleCoordType>(surfaceMesh, FileFormatEnum::VTP,
                               "DEBUG_lsGeomAdvectMesh_contributewoMask.vtp")
          .apply();
      auto mesh = SmartPointer<Mesh<T>>::New();
      if (maskLevelSet != nullptr) {
        ToMesh<T, D>(maskLevelSet, mesh).apply();
        VTKWriter<T>(mesh, FileFormatEnum::VTP,
                     "DEBUG_lsGeomAdvectMesh_mask.vtp")
            .apply();
      }
      ToMesh<T, D>(levelSet, mesh).apply();
      VTKWriter<T>(mesh, FileFormatEnum::VTP,
                   "DEBUG_lsGeomAdvectMesh_initial.vtp")
          .apply();
    }

#endif

    const auto &surfaceNodes = surfaceMesh->getNodes();

    // initialize with segmentation for whole range
    typename viennahrle::Domain<T, D>::IndexPoints segmentation;

    {
      unsigned long long numPoints = 1;
      unsigned long long pointsPerDimension[D];
      for (unsigned i = 0; i < D; ++i) {
        pointsPerDimension[i] = numPoints;
        numPoints *= max[i] - min[i];
      }
      unsigned long numberOfSegments = domain.getNumberOfSegments();
      unsigned long long pointsPerSegment = numPoints / numberOfSegments;
      unsigned long long pointId = 0;
      for (unsigned i = 0; i < numberOfSegments - 1; ++i) {
        pointId = pointsPerSegment * (i + 1);
        viennahrle::Index<D> segmentPoint;
        for (int j = D - 1; j >= 0; --j) {
          segmentPoint[j] = pointId / (pointsPerDimension[j]) + min[j];
          pointId %= pointsPerDimension[j];
        }
        segmentation.push_back(segmentPoint);
      }
    }

    typedef std::vector<std::pair<viennahrle::Index<D>, T>> PointValueVector;
    std::vector<PointValueVector> newPoints;
    newPoints.resize(domain.getNumberOfSegments());

    const T initialDistance = (distIsPositive)
                                  ? std::numeric_limits<double>::max()
                                  : std::numeric_limits<double>::lowest();

#ifndef NDEBUG
    {
      std::ostringstream oss;
      oss << "GeomAdvect: Min: " << min << ", Max: " << max << std::endl;
      VIENNACORE_LOG_DEBUG(oss.str());
    }
#endif
// set up multithreading
#pragma omp parallel for
    for (unsigned p = 0; p < domain.getNumberOfSegments(); ++p) {

      viennahrle::Index<D> startVector;
      if (p == 0) {
        startVector = min;
      } else {
        startVector = segmentation[p - 1];
        incrementIndices(startVector, min, max);
      }

      viennahrle::Index<D> endVector =
          (p != static_cast<int>(domain.getNumberOfSegments() - 1))
              ? segmentation[p]
              : grid.incrementIndices(max);

      viennahrle::ConstSparseIterator<DomainType> checkIt(levelSet->getDomain(),
                                                          startVector);

      // Mask iterator for checking whether inside mask or not
      std::unique_ptr<viennahrle::ConstSparseIterator<DomainType>> maskIt =
          nullptr;
      if (maskLevelSet != nullptr) {
        maskIt = std::make_unique<viennahrle::ConstSparseIterator<DomainType>>(
            maskLevelSet->getDomain(), startVector);
      }

      // Iterate through the bounds of new lsDomain lexicographically
      for (viennahrle::Index<D> currentIndex = startVector;
           currentIndex <= endVector;
           incrementIndices(currentIndex, min, max)) {
        // if point is already full in old level set, skip it
        checkIt.goToIndicesSequential(currentIndex);
        T oldValue = checkIt.getValue();
        // if run is already negative undefined, just ignore the point
        if (distIsPositive) {
          if (oldValue < -cutoffValue) {
            continue;
          }
        } else if (oldValue > cutoffValue) {
          continue;
        }

        VectorType<hrleCoordType, 3> currentCoords{};
        VectorType<hrleCoordType, 3> currentDistMin{};
        VectorType<hrleCoordType, 3> currentDistMax{};

        for (unsigned i = 0; i < D; ++i) {
          currentCoords[i] = currentIndex[i] * gridDelta;

          currentDistMin[i] = currentIndex[i] - std::abs(distMin[i]);
          if (currentDistMin[i] < grid.getMinGridPoint(i)) {
            currentDistMin[i] = grid.getMinGridPoint(i);
          }
          currentDistMin[i] *= gridDelta;

          currentDistMax[i] = currentIndex[i] + std::abs(distMax[i]);
          if (currentDistMin[i] > grid.getMaxGridPoint(i)) {
            currentDistMin[i] = grid.getMaxGridPoint(i);
          }
          currentDistMax[i] *= gridDelta;
        }

        T distance = initialDistance;

        unsigned long currentPointId = 0;
        // now check which surface points contribute to currentIndex
        for (auto surfIt = surfaceNodes.begin(); surfIt != surfaceNodes.end();
             ++surfIt, ++currentPointId) {

          auto &currentNode = *surfIt;

          // if we are outside min/max go to next index inside
          {
            bool outside = false;
            for (unsigned i = 0; i < D; ++i) {
              if ((currentNode[i] < currentDistMin[i]) ||
                  (currentNode[i] > currentDistMax[i])) {
                outside = true;
                break;
              }
            }
            if (outside) {
              continue;
            }
          }

          // TODO: does this really save time? Try without it.
          if (!dist->isInside(currentNode, currentCoords, 2 * gridDelta)) {
            continue;
          }

          // get filling fraction from distance to dist surface
          auto pointId = currentPointId;
          if (!useSurfacePointId) {
            pointId = pointIdTranslator->find(currentPointId)->second;
          }
          T tmpDistance =
              dist->getSignedDistance(currentNode, currentCoords, pointId) /
              gridDelta;

          // if cell is far within a distribution, set it filled
          if (distIsPositive) {
            if (tmpDistance <= -cutoffValue) {
              distance = std::numeric_limits<T>::lowest();
              break;
            }

            if (tmpDistance < distance) {
              distance = tmpDistance;
            }
          } else {
            if (tmpDistance >= cutoffValue) {
              distance = std::numeric_limits<T>::max();
              break;
            }

            if (tmpDistance > distance) {
              distance = tmpDistance;
            }
          }
        }

        // TODO: There are still issues with positive box distributions
        // if there is a mask used!
        // if point is part of the mask, keep smaller value
        if (maskLevelSet != nullptr) {
          maskIt->goToIndicesSequential(currentIndex);

          // if dist is positive, flip logic of comparison
          if (distIsPositive ^
              (std::abs(oldValue - maskIt->getValue()) < 1e-6)) {
            if (!distIsPositive && std::abs(oldValue) <= cutoffValue) {
              newPoints[p].push_back(std::make_pair(currentIndex, oldValue));
              continue;
            }
          } else {
            if (distance != initialDistance) {
              distance = std::min(maskIt->getValue(), distance);
            } else if (distIsPositive || oldValue >= 0.) {
              newPoints[p].push_back(std::make_pair(currentIndex, oldValue));
              continue;
            }
          }
        }

        if (std::abs(distance) <= cutoffValue) {
          // avoid using distribution in wrong direction
          if (distIsPositive && oldValue >= 0.) {
            newPoints[p].push_back(std::make_pair(currentIndex, distance));
          } else if (!distIsPositive && oldValue <= 0.) {
            // if we are etching, need to make sure, we are not inside mask
            if (maskIt == nullptr || maskIt->getValue() > -cutoffValue) {
              newPoints[p].push_back(std::make_pair(currentIndex, distance));
            }
          } else {
            // this only happens if distribution is very small, < 2 * gridDelta
            newPoints[p].push_back(std::make_pair(currentIndex, oldValue));
          }
        }
      }
    }

    // copy all points into the first vector
    {
      unsigned long long numberOfPoints = newPoints[0].size();
      for (unsigned i = 1; i < domain.getNumberOfSegments(); ++i) {
        numberOfPoints += newPoints[i].size();
      }
      newPoints[0].reserve(numberOfPoints);
      for (unsigned i = 1; i < domain.getNumberOfSegments(); ++i) {
        std::move(std::begin(newPoints[i]), std::end(newPoints[i]),
                  std::back_inserter(newPoints[0]));
      }
    }

    auto mesh = SmartPointer<Mesh<T>>::New();
    // output all points directly to mesh
    {
      std::vector<T> scalarData;
      for (auto it = newPoints[0].begin(); it != newPoints[0].end(); ++it) {
        Vec3D<T> node{};
        for (unsigned i = 0; i < D; ++i) {
          node[i] = T((it->first)[i]) * gridDelta;
        }

        mesh->insertNextNode(node);
        std::array<unsigned, 1> vertex{};
        vertex[0] = mesh->vertices.size();
        mesh->insertNextVertex(vertex);
        scalarData.push_back(it->second);
      }
      mesh->cellData.insertNextScalarData(scalarData, "LSValues");
    }

#ifndef NDEBUG // if in debug build
    VIENNACORE_LOG_DEBUG("GeomAdvect: Writing final mesh...");
    VTKWriter<T>(mesh, FileFormatEnum::VTP, "DEBUG_lsGeomAdvectMesh_final.vtp")
        .apply();
#endif

    FromMesh<T, D>(levelSet, mesh).apply();
    // the closest points and normals of the old surface no longer apply
    if (euclidean)
      levelSet->getPointData().clear();

#ifndef NDEBUG // if in debug build
    VIENNACORE_LOG_DEBUG("GeomAdvect: Writing final LS...");
    ToMesh<T, D>(levelSet, mesh).apply();
    VTKWriter<T>(mesh, FileFormatEnum::VTP, "DEBUG_lsGeomAdvectLS_final.vtp")
        .apply();
#endif

    Prune<T, D>(levelSet).apply();

    levelSet->getDomain().segment();
    levelSet->finalize(1);

    // a Euclidean level set returns to its width, with the closest points of
    // the dual contouring mesh of its new values
    Expand<T, D>(levelSet, euclidean ? std::max(euclideanWidth, 2) : 2).apply();

    dist->finalize();
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(GeometricAdvect)

} // namespace viennals
