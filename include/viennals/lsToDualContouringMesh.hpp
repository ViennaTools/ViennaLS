#pragma once

#include <lsPreCompileMacros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <numeric>
#include <vector>

#include <hrleSparseCellIterator.hpp>
#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMesh.hpp>

#include <Eigen/Dense>

namespace viennals {

using namespace viennacore;

/// Extracts the surface of a level set with dual contouring (Ju et al. 2002).
/// Every cell crossed by the surface receives one vertex at the least-squares
/// intersection of the tangent planes on its crossing edges, and every
/// crossing edge becomes a quad (3D) or a line (2D) between the vertices of
/// the cells around it. With the closest points and normals stored by the
/// Euclidean FromSurfaceMesh the tangent planes are exact, so sharp edges and
/// corners of the converted mesh are reproduced. Without them the crossings
/// are interpolated linearly and the normals estimated from the values of the
/// cell. Planes with a priority (AdvectHermite) are chosen by priority first.
/// In 3D a cell crossed by a face narrower than the cell, with both of its
/// edges, gets a vertex on each edge, and polygons on that face close the
/// mesh between them. Vertices of neighbouring cells that coincide are
/// welded, and each quad is split along the diagonal that leaves no triangle
/// without area or facing against the tangent plane of its crossing edge,
/// where the other does not. In 3D a vertex that its planes place beyond its
/// cell, onto the vertex of the cell it falls into, welds to it, so that the
/// mesh folds nowhere along sharp edges (setAvoidFolds).
template <class T, int D> class ToDualContouringMesh {
  using hrleDomainType = typename Domain<T, D>::DomainType;
  using hrleIndex = viennahrle::Index<D>;
  using ConstSparseCellIterator =
      viennahrle::ConstSparseCellIterator<hrleDomainType>;
  using VectorD = Eigen::Matrix<T, D, 1>;
  using MatrixD = Eigen::Matrix<T, D, D>;

  SmartPointer<Domain<T, D>> levelSet = nullptr;
  SmartPointer<Mesh<T>> mesh = nullptr;
  T vertexMargin = 1.;
  bool avoidFolds = true;
  // rounding of the precision, below the fixed tolerances with double
  static constexpr T roundOff = T(1024) * std::numeric_limits<T>::epsilon();
  T singularValueThreshold = 0.01;

  // position of corner c of the cell with lowest corner idx, grid spacings
  static VectorD cornerPosition(const hrleIndex &idx, int c) {
    VectorD p;
    for (int i = 0; i < D; ++i)
      p(i) = idx[i] + ((c >> i) & 1);
    return p;
  }

  // gradient of the (bi/tri)linear interpolant of the cell at local position
  // u in [0,1]^D
  static VectorD interpolantGradient(const ConstSparseCellIterator &cellIt,
                                     const VectorD &u) {
    VectorD g = VectorD::Zero();
    for (int c = 0; c < (1 << D); ++c) {
      const T value = cellIt.getCorner(c).getValue();
      for (int i = 0; i < D; ++i) {
        T w = ((c >> i) & 1) ? T(1) : T(-1);
        for (int j = 0; j < D; ++j)
          if (j != i)
            w *= ((c >> j) & 1) ? u(j) : T(1) - u(j);
        g(i) += w * value;
      }
    }
    return g;
  }

  // Hermite point x and unit normal n of the crossing on the edge from
  // corner a to corner b (grid spacings). The plane is the stored plane of
  // highest positive priority that cuts the edge, of an end point or else of
  // another node of the cell. An edge without one gets, with priorities, the
  // interpolated crossing and no plane (n = 0), and otherwise the plane of
  // its nearer end point. source is the corner of the cell whose plane it
  // took, -1 for none.
  void hermiteData(const ConstSparseCellIterator &cellIt, int a, int b,
                   const typename PointData<T>::VectorDataType *closest,
                   const typename PointData<T>::VectorDataType *normals,
                   const typename PointData<T>::ScalarDataType *priority,
                   T gridDelta, VectorD &x, VectorD &n, int &source) const {
    source = -1;
    const hrleIndex &idx = cellIt.getIndices();
    const VectorD pa = cornerPosition(idx, a), pb = cornerPosition(idx, b);
    const T va = cellIt.getCorner(a).getValue();
    const T vb = cellIt.getCorner(b).getValue();
    const T tLinear = va / (va - vb);

    if (closest != nullptr && normals != nullptr) {
      // the stored plane of node c of the cell, in grid spacings, and its
      // priority; a corner across a periodic boundary holds its closest
      // point in the frame of its wrapped node
      const auto &grid = levelSet->getGrid();
      auto storedPlane = [&](int c, VectorD &cp, VectorD &nrm, T &rank) {
        const auto &corner = cellIt.getCorner(c);
        if (!corner.isDefined())
          return false;
        rank = (priority != nullptr) ? (*priority)[corner.getPointId()] : T(1);
        for (int i = 0; i < D; ++i) {
          cp(i) = (*closest)[corner.getPointId()][i] / gridDelta;
          nrm(i) = (*normals)[corner.getPointId()][i];
        }
        for (int i = 0; i < D; ++i)
          if (grid.isBoundaryPeriodic(i)) {
            const auto global = idx[i] + ((c >> i) & 1);
            cp(i) += T(global - grid.globalIndex2LocalIndex(i, global));
          }
        return rank > T(0);
      };
      const T tolerance = std::max(T(1e-9), roundOff);
      const T rise = vb - va;
      // the stored plane of highest priority that cuts the edge, nearest to
      // the interpolated crossing, of an end point or else of another node of
      // the cell. At the thin tip of a cavity the end points store the floor
      // below and a wall beside, and a node above stores the ceiling that
      // crosses the edge. The plane of another node points towards the
      // positive end, as the outward normal of the surface crossing the edge
      // does, and cuts the edge within a tenth of it from the interpolated
      // crossing, within a quarter with priorities, where the alternative is
      // no plane
      auto cuttingPlane = [&](bool ends) {
        T bestT = 0., bestGap = std::numeric_limits<T>::max(), bestRank = 0.;
        VectorD bestN = VectorD::Zero();
        int bestC = -1;
        for (int c = 0; c < (1 << D); ++c) {
          VectorD cp, nrm;
          T rank;
          if ((c == a || c == b) != ends || !storedPlane(c, cp, nrm, rank))
            continue;
          const T den = nrm.dot(pb - pa);
          if (std::abs(den) < T(1e-12) || (!ends && den * rise <= T(0)))
            continue;
          const T t = nrm.dot(cp - pa) / den;
          if (t < -tolerance || t > T(1) + tolerance ||
              (!ends && std::abs(t - tLinear) >
                            (priority != nullptr ? T(0.25) : T(0.1))))
            continue;
          if (rank > bestRank ||
              (rank == bestRank && std::abs(t - tLinear) < bestGap)) {
            bestRank = rank;
            bestGap = std::abs(t - tLinear);
            bestT = std::clamp(t, T(0), T(1));
            bestN = nrm;
            bestC = c;
          }
        }
        if (bestRank <= T(0))
          return false;
        x = pa + bestT * (pb - pa);
        n = bestN;
        source = bestC;
        return true;
      };
      // a plane of priority above 2 that a node of the cell stores, as the
      // plane of the end of a cavity below a mask past which the node lies,
      // goes first wherever it cuts the edge
      if (priority != nullptr) {
        T bestGap = std::numeric_limits<T>::max();
        bool found = false;
        for (int c = 0; c < (1 << D); ++c) {
          VectorD cp, nrm;
          T rank;
          if (!storedPlane(c, cp, nrm, rank) || rank <= T(2))
            continue;
          const T den = nrm.dot(pb - pa);
          if (std::abs(den) < T(1e-12) || den * rise <= T(0))
            continue;
          const T t = nrm.dot(cp - pa) / den;
          if (t < -tolerance || t > T(1) + tolerance)
            continue;
          if (std::abs(t - tLinear) < bestGap) {
            bestGap = std::abs(t - tLinear);
            x = pa + std::clamp(t, T(0), T(1)) * (pb - pa);
            n = nrm;
            source = c;
            found = true;
          }
        }
        if (found)
          return;
      }
      if (cuttingPlane(true) || cuttingPlane(false))
        return;
      // with priorities an edge without such a plane gets the interpolated
      // crossing and no plane
      if (priority != nullptr) {
        x = pa + tLinear * (pb - pa);
        n = VectorD::Zero();
        return;
      }
      // else the plane an end point stores, the nearer one first, through
      // its closest point, with the foot of the interpolated crossing on it
      // as its point
      const int nearer = (std::abs(va) <= std::abs(vb)) ? a : b;
      for (int end : {nearer, a + b - nearer}) {
        VectorD cp, nrm;
        T rank;
        if (!storedPlane(end, cp, nrm, rank))
          continue;
        const VectorD crossing = pa + tLinear * (pb - pa);
        x = crossing - nrm.dot(crossing - cp) * nrm;
        n = nrm;
        source = end;
        return;
      }
    }

    // estimated Hermite data: linear crossing, normal of the interpolant
    x = pa + tLinear * (pb - pa);
    VectorD u;
    for (int i = 0; i < D; ++i)
      u(i) = x(i) - idx[i];
    n = interpolantGradient(cellIt, u);
    const T length = n.norm();
    // corners outside the band carry the background value, which can make
    // the gradient overflow; such an edge gets no plane
    if (length > T(0) && std::isfinite(length))
      n /= length;
    else
      n = VectorD::Zero();
  }

public:
  /// Optional scalar point data: plane priority of every node, 0 for none.
  static constexpr char priorityLabel[] = "HermitePriority";

  ToDualContouringMesh() = default;

  ToDualContouringMesh(SmartPointer<Domain<T, D>> passedLevelSet,
                       SmartPointer<Mesh<T>> passedMesh)
      : levelSet(passedLevelSet), mesh(passedMesh) {}

  void setLevelSet(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSet = passedLevelSet;
  }

  void setMesh(SmartPointer<Mesh<T>> passedMesh) { mesh = passedMesh; }

  /// Vertices stay within their cell plus this margin, in grid spacings.
  /// Defaults to 1, which lets a vertex reach a corner in a neighbouring cell.
  void setVertexMargin(T passedMargin) { vertexMargin = passedMargin; }

  /// Singular values of a cell's quadric (squared singular values of its
  /// planes) below this fraction of the largest one are dropped. Defaults to
  /// 0.01, the threshold 0.1 of Ju et al. on the planes, which keeps two
  /// planes whose normals differ by more than 11.4 degrees as a sharp
  /// feature.
  void setSingularValueThreshold(T passedThreshold) {
    singularValueThreshold = passedThreshold;
  }

  /// Whether a 3D mesh folds nowhere (default): a vertex that its planes
  /// place beyond its cell, into a cell whose vertex lies on those planes as
  /// well, welds to that vertex, and a vertex free along an edge takes a
  /// corner only within half a grid spacing of its cell. Without it, as for
  /// the distances of Reinitialize, a vertex free along an edge takes a
  /// corner within the vertex margin, and vertices beyond their cells stay:
  /// slivers along sharp edges then face anywhere, but over many rebuilds of
  /// AdvectHermite the band keeps sharp edges more accurately.
  void setAvoidFolds(bool passedAvoidFolds) { avoidFolds = passedAvoidFolds; }

  void apply() {
    if (levelSet == nullptr) {
      Logger::getInstance()
          .addError("No level set was passed to ToDualContouringMesh.")
          .print();
      return;
    }
    if (mesh == nullptr) {
      Logger::getInstance()
          .addError("No mesh was passed to ToDualContouringMesh.")
          .print();
      return;
    }

    mesh->clear();
    const T gridDelta = levelSet->getGrid().getGridDelta();
    auto &pointData = levelSet->getPointData();
    const auto *closest = pointData.getVectorData(
        FromSurfaceMesh<T, D>::closestPointsLabel, true);
    const auto *normals = pointData.getVectorData(
        FromSurfaceMesh<T, D>::surfaceNormalsLabel, true);
    if (closest == nullptr || normals == nullptr ||
        closest->size() != levelSet->getNumberOfPoints() ||
        normals->size() != levelSet->getNumberOfPoints()) {
      closest = nullptr;
      normals = nullptr;
    }
    const auto *priority = pointData.getScalarData(priorityLabel, true);
    // nodes whose plane is a tangent plane of a convex rounded part, as
    // AdvectHermite marks them
    const auto *smoothData = pointData.getScalarData("HermiteSmooth", true);
    if (smoothData != nullptr &&
        smoothData->size() != levelSet->getNumberOfPoints())
      smoothData = nullptr;
    if (closest == nullptr || priority == nullptr ||
        priority->size() != levelSet->getNumberOfPoints())
      priority = nullptr;

    const auto &grid = levelSet->getGrid();

    // the least-squares point of the planes around the mass point, with a
    // truncated pseudo-inverse, kept within the cell plus the margin, and its
    // largest distance to the planes
    auto solve = [&](const std::vector<std::pair<VectorD, T>> &planes,
                     const VectorD &mass, const hrleIndex &idx,
                     VectorD &vertex) {
      MatrixD ATA = MatrixD::Zero();
      VectorD ATb = VectorD::Zero();
      for (const auto &[n, offset] : planes) {
        ATA += n * n.transpose();
        ATb += n * offset;
      }
      Eigen::SelfAdjointEigenSolver<MatrixD> eigen(ATA);
      const VectorD lambda = eigen.eigenvalues();
      const T largest = lambda.cwiseAbs().maxCoeff();
      VectorD inverse = VectorD::Zero();
      for (int i = 0; i < D; ++i)
        if (largest > T(0) &&
            std::abs(lambda(i)) > singularValueThreshold * largest)
          inverse(i) = T(1) / lambda(i);
      const MatrixD V = eigen.eigenvectors();
      vertex =
          mass + V * inverse.asDiagonal() * V.transpose() * (ATb - ATA * mass);
      for (int i = 0; i < D; ++i)
        vertex(i) = std::clamp(vertex(i), idx[i] - vertexMargin,
                               idx[i] + T(1) + vertexMargin);
      T residual = T(0);
      for (const auto &[n, offset] : planes)
        residual = std::max(residual, std::abs(n.dot(vertex) - offset));
      return residual;
    };

    auto insertVertex = [&](const VectorD &vertex) {
      Vec3D<T> node{T(0), T(0), T(0)};
      for (int i = 0; i < D; ++i)
        node[i] = vertex(i) * gridDelta;
      mesh->minimumExtent = Min(mesh->minimumExtent, node);
      mesh->maximumExtent = Max(mesh->maximumExtent, node);
      return mesh->insertNextNode(node);
    };

    // a cell whose vertex a cell across a periodic or reflective boundary
    // copies keeps one vertex
    auto copiedAcross = [&](const hrleIndex &idx) {
      for (int i = 0; i < D; ++i)
        if ((grid.isBoundaryPeriodic(i) || grid.isBoundaryReflective(i)) &&
            (idx[i] <= grid.getMinGridPoint(i) ||
             idx[i] + 1 >= grid.getMaxGridPoint(i)))
          return true;
      return false;
    };

    // a narrow face through a 3D cell, as the face a silhouette edge sweeps
    // in directional deposition while it is narrower than a cell: around the
    // surface in the cell the plane of that face appears twice, between the
    // planes on either side of it, where at a corner every plane appears
    // once. The two edges of the face then get a vertex each, from the
    // planes on one side with the plane of the face. The crossing edges on
    // the one side take the first, those on the face and on the other side
    // the second (bit a * D + d of takesSecond for the edge from corner a
    // along d). Taken where the two vertices lie on their planes within a
    // tenth of the distance of the single vertex to its planes.
    auto splitCell = [&](const std::vector<int> &edges,
                         const std::vector<int> &planeOf,
                         const std::vector<VectorD> &points,
                         const std::vector<std::pair<VectorD, T>> &planes,
                         const hrleIndex &idx, T residual, VectorD &first,
                         VectorD &second, std::uint32_t &takesSecond) {
      const int m = int(edges.size());
      // two crossing edges on a face of the cell follow each other around
      // the surface; four on a face leave its course open
      std::array<int, (1 << D) * D> slot;
      slot.fill(-1);
      for (int k = 0; k < m; ++k)
        slot[edges[k]] = k;
      std::vector<std::array<int, 2>> link(m, {-1, -1});
      for (int k = 0; k < D; ++k)
        for (int s = 0; s < 2; ++s) {
          std::array<int, 4> onFace;
          int found = 0;
          for (int d = 0; d < D; ++d) {
            if (d == k)
              continue;
            const int o = 3 - k - d;
            for (int t = 0; t < 2; ++t) {
              const int e = slot[((s << k) | (t << o)) * D + d];
              if (e >= 0)
                onFace[found++] = e;
            }
          }
          if (found == 4)
            return false;
          if (found != 2)
            continue;
          for (int q = 0; q < 2; ++q) {
            auto &l = link[onFace[q]];
            if (l[0] < 0)
              l[0] = onFace[1 - q];
            else if (l[1] < 0)
              l[1] = onFace[1 - q];
            else
              return false;
          }
        }
      // the crossings around the surface, one loop through all of them
      std::vector<int> order;
      for (int previous = -1, current = 0;;) {
        if (link[current][1] < 0)
          return false;
        order.push_back(current);
        const int next = (link[current][0] != previous) ? link[current][0]
                                                        : link[current][1];
        previous = current;
        current = next;
        if (current == 0 || int(order.size()) > m)
          break;
      }
      if (int(order.size()) != m)
        return false;
      // runs of crossings on one plane, from a change of plane on
      int start = -1;
      for (int i = 0; i < m && start < 0; ++i)
        if (planeOf[order[i]] != planeOf[order[(i + m - 1) % m]])
          start = i;
      if (start < 0)
        return false;
      std::vector<int> around(m), runPlane, runBegin;
      for (int i = 0; i < m; ++i)
        around[i] = order[(start + i) % m];
      for (int i = 0; i < m; ++i)
        if (i == 0 || planeOf[around[i]] != planeOf[around[i - 1]]) {
          runPlane.push_back(planeOf[around[i]]);
          runBegin.push_back(i);
        }
      // the plane of the face in two runs, every other plane in one
      const int runs = int(runPlane.size());
      int face = -1, firstRun = -1, secondRun = -1;
      for (int r = 0; r < runs; ++r) {
        int before = -1;
        for (int q = 0; q < r; ++q)
          if (runPlane[q] == runPlane[r]) {
            if (before >= 0)
              return false;
            before = q;
          }
        if (before < 0)
          continue;
        if (face >= 0)
          return false;
        face = runPlane[r];
        firstRun = before;
        secondRun = r;
      }
      if (face < 0 || secondRun - firstRun < 2 ||
          runs - secondRun + firstRun < 2)
        return false;
      // the face meets the planes on either side at edges as sharp as the
      // truncation of the quadric keeps them; planes closer to it, as the
      // sampled support planes of a fan around one direction, are one face
      const T sharp =
          (T(1) - singularValueThreshold) / (T(1) + singularValueThreshold);
      for (int k = 0; k < m; ++k)
        if (planeOf[k] != face &&
            planes[planeOf[k]].first.dot(planes[face].first) >= sharp)
          return false;
      // the side between the two runs of the face, and the other side with
      // the face itself
      std::vector<int> onFirst{face}, onSecond{face};
      VectorD massFirst = VectorD::Zero(), massSecond = VectorD::Zero();
      int countFirst = 0, countSecond = 0;
      takesSecond = 0;
      for (int i = 0; i < m; ++i) {
        const int k = around[i];
        const bool firstSide =
            i >= runBegin[firstRun + 1] && i < runBegin[secondRun];
        auto &on = firstSide ? onFirst : onSecond;
        if (std::find(on.begin(), on.end(), planeOf[k]) == on.end())
          on.push_back(planeOf[k]);
        if (firstSide) {
          massFirst += points[k];
          ++countFirst;
        } else {
          massSecond += points[k];
          ++countSecond;
          takesSecond |= std::uint32_t(1) << edges[k];
        }
      }
      std::vector<std::pair<VectorD, T>> planesFirst, planesSecond;
      for (const int p : onFirst)
        planesFirst.push_back(planes[p]);
      for (const int p : onSecond)
        planesSecond.push_back(planes[p]);
      const T residualFirst =
          solve(planesFirst, massFirst / T(countFirst), idx, first);
      const T residualSecond =
          solve(planesSecond, massSecond / T(countSecond), idx, second);
      return std::max(residualFirst, residualSecond) <= T(0.1) * residual;
    };

    // one vertex per cell crossed by the surface, two where a narrow face
    // crosses it; bit a * D + d of the masks for the edge from corner a along
    // d
    struct CellVertices {
      unsigned first = 0, second = 0;
      std::uint32_t crossing = 0, takesSecond = 0;
      // placed halfway between its planes' point and its mass point
      bool halfway = false;
    };
    std::map<hrleIndex, CellVertices> cellVertex;
    // in 2D the cell and mass point of the vertex of every cell whose planes
    // lie within 35 degrees of each other, for the repair of folds below
    std::vector<std::pair<hrleIndex, VectorD>> nearlyParallel;
    std::vector<char> hasNearlyParallel;
    // the normals of the crossing edges of the cells with two vertices, for
    // the polygons that close the mesh between them
    std::map<hrleIndex, std::array<std::array<T, 3>, (1 << D) * D>>
        splitNormals;
    // the distinct planes of every cell, for the cells beside a corner below
    std::map<hrleIndex, std::vector<std::pair<VectorD, T>>> cellPlanes;
    for (ConstSparseCellIterator cellIt(levelSet->getDomain());
         !cellIt.isFinished(); cellIt.next()) {
      unsigned inside = 0;
      for (int c = 0; c < (1 << D); ++c)
        if (cellIt.getCorner(c).getValue() < T(0))
          inside |= 1u << c;
      if (inside == 0 || inside == (1u << (1 << D)) - 1)
        continue;

      VectorD mass = VectorD::Zero();
      int count = 0;
      // each distinct plane enters once: a face crossing several edges of
      // the cell is one constraint, and counting it per edge lets the
      // truncation drop the direction of a plane that crosses fewer edges
      std::vector<std::pair<VectorD, T>> planes;
      // the crossing edges with their Hermite points, normals and planes
      std::vector<int> edges, planeOf;
      std::vector<VectorD> points, edgeNormals;
      CellVertices cell;
      bool allSmooth = true;
      for (int a = 0; a < (1 << D); ++a) {
        for (int d = 0; d < D; ++d) {
          if ((a >> d) & 1)
            continue;
          const int b = a | (1 << d);
          if (((inside >> a) & 1) == ((inside >> b) & 1))
            continue;
          VectorD x, n;
          int source;
          hermiteData(cellIt, a, b, closest, normals, priority, gridDelta, x, n,
                      source);
          // whether the plane is one of those
          if (source < 0 || smoothData == nullptr ||
              (*smoothData)[cellIt.getCorner(source).getPointId()] < T(0.5))
            allSmooth = false;
          mass += x;
          ++count;
          const T offset = n.dot(x);
          int plane = -1;
          for (std::size_t k = 0; k < planes.size() && plane < 0; ++k) {
            const auto &[m, o] = planes[k];
            if ((m - n).norm() < std::max(T(1e-9), roundOff) &&
                std::abs(o - offset) < std::max(T(1e-9), roundOff) *
                                           std::max(T(1), std::abs(offset)))
              plane = int(k);
          }
          if (plane < 0) {
            planes.emplace_back(n, offset);
            plane = int(planes.size()) - 1;
          }
          cell.crossing |= std::uint32_t(1) << (a * D + d);
          edges.push_back(a * D + d);
          planeOf.push_back(n.isZero() ? -1 : plane);
          points.push_back(x);
          edgeNormals.push_back(n);
        }
      }
      mass /= T(count);

      const hrleIndex &idx = cellIt.getIndices();
      VectorD vertex, second;
      const T residual = solve(planes, mass, idx, vertex);
      // on a convex rounded part, as the arc a deposition rounds an edge
      // into, the tangent planes of the crossings meet outside it by about
      // R <theta^2> / 2 and the mass point of the crossings lies inside by as
      // much: halfway between them lies on the surface to second order. Only
      // where every plane of the cell is such a tangent plane and the planes
      // lie within 45 degrees of each other
      if (allSmooth && planes.size() > 1) {
        bool close = true;
        for (std::size_t a = 0; a < planes.size() && close; ++a)
          for (std::size_t b = a + 1; b < planes.size() && close; ++b)
            close = planes[a].first.dot(planes[b].first) >
                    T(0.70710678118654752); // cos 45°
        if (close) {
          vertex = (vertex + mass) / T(2);
          cell.halfway = true;
        }
      }
      // a 2D cell whose two planes meet outside it, as at the end of a cavity
      // below a mask, where the plane of the interface meets the arc's
      // tangent at a crossing a little before the end beyond it, takes the
      // point where one of them meets the plane that a corner of the cell off
      // its crossing edges stores, where that point lies in the cell, as the
      // plane of the end that a node of the mask past it holds
      if constexpr (D == 2) {
        bool outside = false;
        for (int i = 0; i < D; ++i)
          outside = outside || vertex(i) < idx[i] - T(0.05) ||
                    vertex(i) > idx[i] + T(1.05);
        if (outside && planes.size() == 2 && closest != nullptr &&
            !copiedAcross(idx)) {
          unsigned onCrossing = 0;
          for (const int e : edges)
            onCrossing |= (1u << (e / D)) | (1u << ((e / D) | (1 << (e % D))));
          bool taken = false;
          for (int c = 0; c < (1 << D) && !taken; ++c) {
            if ((onCrossing >> c) & 1)
              continue;
            const auto &corner = cellIt.getCorner(c);
            if (!corner.isDefined())
              continue;
            VectorD cp, nrm;
            for (int i = 0; i < D; ++i) {
              cp(i) = (*closest)[corner.getPointId()][i] / gridDelta;
              nrm(i) = (*normals)[corner.getPointId()][i];
            }
            for (const auto &[m, o] : planes) {
              const T det = nrm(0) * m(1) - nrm(1) * m(0);
              if (std::abs(det) < T(0.1))
                continue;
              const T op = nrm.dot(cp);
              VectorD q;
              q(0) = (op * m(1) - nrm(1) * o) / det;
              q(1) = (nrm(0) * o - op * m(0)) / det;
              bool in = true;
              for (int i = 0; i < D; ++i)
                in = in && q(i) >= idx[i] - T(0.05) && q(i) <= idx[i] + T(1.05);
              if (in) {
                vertex = q;
                taken = true;
                break;
              }
            }
          }
        }
      }
      bool split = false;
      if constexpr (D == 3)
        split =
            residual > T(1e-3) && !copiedAcross(idx) &&
            std::find(planeOf.begin(), planeOf.end(), -1) == planeOf.end() &&
            splitCell(edges, planeOf, points, planes, idx, residual, vertex,
                      second, cell.takesSecond);
      cell.first = insertVertex(vertex);
      cell.second = cell.first;
      if constexpr (D == 2) {
        bool close = planes.size() > 1;
        for (std::size_t a = 0; a < planes.size(); ++a)
          for (std::size_t b = a + 1; b < planes.size(); ++b)
            close = close && planes[a].first.dot(planes[b].first) >
                                 T(0.81915204428899178); // cos 35°
        if (close) {
          nearlyParallel.resize(cell.first + 1);
          hasNearlyParallel.resize(cell.first + 1, 0);
          nearlyParallel[cell.first] = {idx, mass};
          hasNearlyParallel[cell.first] = 1;
        }
      }
      if (split) {
        cell.second = insertVertex(second);
        auto &stored = splitNormals[idx];
        for (std::size_t k = 0; k < edges.size(); ++k)
          for (int i = 0; i < D; ++i)
            stored[edges[k]][i] = edgeNormals[k](i);
      } else {
        cell.takesSecond = 0;
      }
      cellVertex[idx] = cell;
      cellPlanes[idx] = planes;
    }

    // a cell whose planes leave its vertex free along a line, as a 3D cell
    // beside a corner whose crossing on a face takes the support plane that
    // a fan of the corner stores, or a 2D cell beside a corner whose two
    // planes the quadric takes for one, takes the corner that the planes of
    // the cells around it fix, where its own planes pass through that point:
    // the vertex of a cell there whose planes fix it alone (the two weld),
    // or else the point. Along an edge (3D) or a straight line (2D) the
    // planes around fix no point, and its vertices stay. Cells whose
    // neighbours reach across a periodic or reflective boundary keep theirs,
    // since the planes there are those of another frame, and so do cells of
    // estimated planes
    if (closest != nullptr) {
      auto nearBoundary = [&](const hrleIndex &idx) {
        for (int i = 0; i < D; ++i)
          if ((grid.isBoundaryPeriodic(i) || grid.isBoundaryReflective(i)) &&
              (idx[i] <= grid.getMinGridPoint(i) + 1 ||
               idx[i] + 2 >= grid.getMaxGridPoint(i)))
            return true;
        return false;
      };
      // the number of directions the planes fix
      auto rankOf = [&](const std::vector<std::pair<VectorD, T>> &planes) {
        MatrixD ATA = MatrixD::Zero();
        for (const auto &[n, offset] : planes)
          ATA += n * n.transpose();
        Eigen::SelfAdjointEigenSolver<MatrixD> eigen(ATA);
        const VectorD lambda = eigen.eigenvalues();
        const T largest = lambda.cwiseAbs().maxCoeff();
        int rank = 0;
        for (int i = 0; i < D; ++i)
          rank += largest > T(0) &&
                  std::abs(lambda(i)) > singularValueThreshold * largest;
        return rank;
      };
      std::vector<std::pair<unsigned, Vec3D<T>>> moves;
      for (const auto &[idx, cell] : cellVertex) {
        if (cell.second != cell.first || nearBoundary(idx))
          continue;
        const auto &own = cellPlanes[idx];
        if (rankOf(own) != D - 1)
          continue;
        // in 2D only two planes more than half a degree apart, which the
        // quadric takes for one up to some 11 degrees, as the wall and a
        // support plane of a fan of the trench bottom 9 degrees off it: a
        // straight line, also one held twice or at two offsets close
        // together, as a floor, leaves nothing to fix
        if constexpr (D == 2) {
          bool apart = false;
          for (std::size_t a = 0; a < own.size() && !apart; ++a)
            for (std::size_t b = a + 1; b < own.size() && !apart; ++b)
              apart = own[a].first.norm() > T(0.5) &&
                      own[b].first.norm() > T(0.5) &&
                      own[a].first.dot(own[b].first) <
                          T(0.99996192306417131); // cos 0.5°
          if (!apart)
            continue;
        }
        std::vector<std::pair<VectorD, T>> around = own;
        std::vector<hrleIndex> others;
        for (int code = 0; code < (D == 2 ? 9 : 27); ++code) {
          hrleIndex other = idx;
          int rest = code;
          for (int i = 0; i < D; ++i) {
            other[i] += rest % 3 - 1;
            rest /= 3;
          }
          if (other == idx)
            continue;
          auto found = cellPlanes.find(other);
          if (found == cellPlanes.end())
            continue;
          others.push_back(other);
          around.insert(around.end(), found->second.begin(),
                        found->second.end());
        }
        // a point the planes around fix well, as at a corner whose faces
        // meet at more than about 13 degrees
        MatrixD ATA = MatrixD::Zero();
        VectorD ATb = VectorD::Zero();
        for (const auto &[n, offset] : around) {
          ATA += n * n.transpose();
          ATb += n * offset;
        }
        Eigen::SelfAdjointEigenSolver<MatrixD> eigen(ATA);
        const VectorD lambda = eigen.eigenvalues();
        if (!(lambda(0) > T(0.05) * lambda(D - 1)))
          continue;
        VectorD point =
            eigen.eigenvectors() *
            (eigen.eigenvectors().transpose() * ATb).cwiseQuotient(lambda);
        T ownResidual = T(0), residual = T(0);
        for (const auto &[n, offset] : own)
          ownResidual = std::max(ownResidual, std::abs(n.dot(point) - offset));
        for (const auto &[n, offset] : around)
          residual = std::max(residual, std::abs(n.dot(point) - offset));
        // in 2D a corner that all the planes around fix exactly, as the
        // corner of a trench bottom before noise reaches it; planes meeting
        // only roughly, as below the foot of a mask wall beside a curved
        // undercut, leave the vertex where it is
        if (ownResidual > ((D == 2) ? T(1e-3) : T(0.01)) ||
            residual > ((D == 2) ? T(1e-3) : T(0.05)))
          continue;
        // avoiding folds in 3D, a corner within half a grid spacing of the
        // cell: one further off, past the vertices of the cells along the
        // edge to it, folds the quads around the cell
        const T reach = (D == 3 && avoidFolds) ? T(0.5) : vertexMargin;
        bool within = true;
        for (int i = 0; i < D; ++i)
          within = within && point(i) >= idx[i] - reach &&
                   point(i) <= idx[i] + T(1) + reach;
        if (!within)
          continue;
        // the vertex of the cell holding the corner, so that the two weld
        Vec3D<T> moved{T(0), T(0), T(0)};
        for (int i = 0; i < D; ++i)
          moved[i] = point(i) * gridDelta;
        for (const auto &other : others) {
          const auto &neighbour = cellVertex.at(other);
          if (neighbour.second != neighbour.first ||
              rankOf(cellPlanes.at(other)) != D)
            continue;
          const auto &node = mesh->getNodes()[neighbour.first];
          T gap = T(0);
          for (int i = 0; i < D; ++i)
            gap = std::max(gap, std::abs(node[i] / gridDelta - point(i)));
          if (gap <= T(0.02)) {
            moved = node;
            break;
          }
        }
        moves.push_back({cell.first, moved});
      }
      for (const auto &[id, moved] : moves) {
        mesh->nodes[id] = moved;
        mesh->minimumExtent = Min(mesh->minimumExtent, moved);
        mesh->maximumExtent = Max(mesh->maximumExtent, moved);
      }

      // avoiding folds in 3D, a vertex that its planes place beyond its
      // cell, into a cell whose vertex lies in it and on those planes as
      // well, where the planes of that cell fix as many directions, takes
      // that vertex (the two weld), as a cell beside a sharp edge whose
      // crossing near the edge takes a support plane of the edge: two
      // vertices on the edge from cells on either side of a face of the grid
      // leave a quad along the edge without area, whose triangles face
      // anywhere
      if (D == 3 && avoidFolds) {
        const T tolerance = T(1e-6);
        auto cellOf = [&](const Vec3D<T> &node, hrleIndex &into) {
          bool outside = false;
          for (int i = 0; i < D; ++i) {
            const T x = node[i] / gridDelta;
            into[i] = static_cast<viennahrle::IndexType>(std::floor(x));
            outside = outside || x - T(into[i]) < tolerance ||
                      x - T(into[i]) > T(1) - tolerance;
          }
          return !outside;
        };
        moves.clear();
        for (const auto &[idx, cell] : cellVertex) {
          if (cell.second != cell.first || cell.halfway || nearBoundary(idx))
            continue;
          const auto &node = mesh->getNodes()[cell.first];
          hrleIndex into = idx;
          if (!cellOf(node, into) || into == idx)
            continue;
          auto found = cellVertex.find(into);
          if (found == cellVertex.end() ||
              found->second.second != found->second.first ||
              found->second.halfway)
            continue;
          const auto &target = mesh->getNodes()[found->second.first];
          hrleIndex holding = into;
          if (!cellOf(target, holding) || holding != into ||
              rankOf(cellPlanes[into]) < rankOf(cellPlanes[idx]))
            continue;
          T residual = T(0);
          for (const auto &[n, offset] : cellPlanes[idx]) {
            T distance = -offset;
            for (int i = 0; i < D; ++i)
              distance += n(i) * target[i] / gridDelta;
            residual = std::max(residual, std::abs(distance));
          }
          if (residual <= T(0.01))
            moves.push_back({cell.first, target});
        }
        for (const auto &[id, moved] : moves)
          mesh->nodes[id] = moved;
      }
    }

    // the vertex of a cell for its crossing edge from corner a along d (edge
    // a * D + d); for a cell across a periodic boundary a copy of the wrapped
    // cell's vertex moved by the period, for a cell across a reflective
    // boundary a copy of its mirror cell's vertex mirrored at the boundary
    std::map<std::pair<hrleIndex, hrleIndex>, unsigned> copiedVertex;
    constexpr viennahrle::IndexType mirrorLow = 1 << 20, mirrorHigh = 1 << 21;
    // the copies and the vertices they copy, which the repair of folds below
    // leaves where they are
    struct Copy {
      unsigned id, source;
    };
    std::vector<Copy> copies;
    auto vertexOf = [&](const hrleIndex &cell, int edge, unsigned &id) {
      hrleIndex source = cell, change;
      bool copied = false;
      for (int i = 0; i < D; ++i) {
        change[i] = 0;
        if (grid.isBoundaryPeriodic(i)) {
          source[i] = grid.globalIndex2LocalIndex(i, cell[i]);
          change[i] = cell[i] - source[i];
        } else if (grid.isBoundaryReflective(i)) {
          const auto low = grid.getMinGridPoint(i),
                     high = grid.getMaxGridPoint(i);
          if (cell[i] < low) {
            source[i] = 2 * low - cell[i] - 1;
            change[i] = mirrorLow;
          } else if (cell[i] >= high) {
            source[i] = 2 * high - cell[i] - 1;
            change[i] = mirrorHigh;
          }
        }
        copied = copied || change[i] != 0;
      }
      auto it = cellVertex.find(source);
      if (it == cellVertex.end())
        return false;
      if (!copied) {
        id = ((it->second.takesSecond >> edge) & 1) ? it->second.second
                                                    : it->second.first;
        return true;
      }
      const auto key = std::make_pair(source, change);
      auto found = copiedVertex.find(key);
      if (found != copiedVertex.end()) {
        id = found->second;
        return true;
      }
      Vec3D<T> node = mesh->getNodes()[it->second.first];
      for (int i = 0; i < D; ++i) {
        if (change[i] == mirrorLow)
          node[i] = 2 * grid.getMinGridPoint(i) * gridDelta - node[i];
        else if (change[i] == mirrorHigh)
          node[i] = 2 * grid.getMaxGridPoint(i) * gridDelta - node[i];
        else
          node[i] += change[i] * gridDelta;
      }
      id = mesh->insertNextNode(node);
      mesh->minimumExtent = Min(mesh->minimumExtent, node);
      mesh->maximumExtent = Max(mesh->maximumExtent, node);
      copiedVertex[key] = id;
      copies.push_back({id, it->second.first});
      return true;
    };

    // a quad (3D) or line (2D) for the crossing edge from corner start along
    // d, from the cells around it, with the normal of the edge's tangent
    // plane; outward means away from the negative side
    std::vector<std::array<unsigned, 2 * (D - 1)>> faces;
    std::vector<VectorD> faceNormals;
    VectorD edgeNormal = VectorD::Zero();
    auto emit = [&](const hrleIndex &start, int d, bool insideLow) {
      if constexpr (D == 3) {
        const int u = (d + 1) % 3, v = (d + 2) % 3;
        // cells around the edge, counterclockwise seen from +d, and the edge
        // in each of them
        const int offsets[4][2] = {{-1, -1}, {0, -1}, {0, 0}, {-1, 0}};
        std::array<unsigned, 4> quad;
        for (int q = 0; q < 4; ++q) {
          hrleIndex cell = start;
          cell[u] += offsets[q][0];
          cell[v] += offsets[q][1];
          const int a = (-offsets[q][0] << u) | (-offsets[q][1] << v);
          if (!vertexOf(cell, a * D + d, quad[q]))
            return;
        }
        if (!insideLow)
          std::reverse(quad.begin(), quad.end());
        faces.push_back(quad);
      } else {
        // lines have the outside on their left
        hrleIndex below = start;
        below[1 - d] -= 1;
        unsigned low, high;
        if (!vertexOf(below, (1 << (1 - d)) * D + d, low) ||
            !vertexOf(start, d, high))
          return;
        const bool lowToHigh = (d == 1) == insideLow;
        if (lowToHigh)
          faces.push_back({low, high});
        else
          faces.push_back({high, low});
      }
      faceNormals.push_back(edgeNormal);
    };

    // every crossing edge once, from the cell whose lowest corner starts it;
    // edges in an upper reflective boundary start at a corner without such a
    // cell and are taken from the cell below the boundary
    for (ConstSparseCellIterator cellIt(levelSet->getDomain());
         !cellIt.isFinished(); cellIt.next()) {
      const hrleIndex idx = cellIt.getIndices();
      const bool insideLow = cellIt.getCorner(0).getValue() < T(0);
      VectorD x;
      int source;
      for (int d = 0; d < D; ++d)
        if (insideLow != (cellIt.getCorner(1 << d).getValue() < T(0))) {
          hermiteData(cellIt, 0, 1 << d, closest, normals, priority, gridDelta,
                      x, edgeNormal, source);
          emit(idx, d, insideLow);
        }
      // upper reflective boundaries the cell touches; each set of them
      // gives the edges that start at the cell corner on all of them
      int upper = 0;
      for (int i = 0; i < D; ++i)
        if (grid.isBoundaryReflective(i) &&
            idx[i] + 1 == grid.getMaxGridPoint(i))
          upper |= 1 << i;
      for (int set = upper; set > 0; set = (set - 1) & upper) {
        hrleIndex start = idx;
        for (int i = 0; i < D; ++i)
          start[i] += (set >> i) & 1;
        const bool inside = cellIt.getCorner(set).getValue() < T(0);
        for (int d = 0; d < D; ++d)
          if (!((set >> d) & 1) &&
              inside != (cellIt.getCorner(set | (1 << d)).getValue() < T(0))) {
            hermiteData(cellIt, set, set | (1 << d), closest, normals, priority,
                        gridDelta, x, edgeNormal, source);
            emit(start, d, inside);
          }
      }
    }

    // where a cell with two vertices goes over from the one to the other
    // between the two crossing edges on a face of the cell, the quads of the
    // two edges leave a gap. It closes with the vertex of the cell across
    // that face, or with its two vertices where that cell goes over as well
    // (that polygon from the lower of the two cells), facing out as the
    // planes of the two edges
    if constexpr (D == 3)
      for (const auto &[idx, outward] : splitNormals) {
        const auto &cell = cellVertex.at(idx);
        auto vertexAt = [&](unsigned id) {
          const auto &node = mesh->getNodes()[id];
          return VectorD(node[0], node[1], node[2]);
        };
        for (int k = 0; k < D; ++k)
          for (int s = 0; s < 2; ++s) {
            std::array<int, 2> pair;
            int found = 0;
            for (int d = 0; d < D; ++d) {
              if (d == k)
                continue;
              const int o = 3 - k - d;
              for (int t = 0; t < 2; ++t) {
                const int e = ((s << k) | (t << o)) * D + d;
                if ((cell.crossing >> e) & 1) {
                  if (found < 2)
                    pair[found] = e;
                  ++found;
                }
              }
            }
            if (found != 2)
              continue;
            const bool second0 = (cell.takesSecond >> pair[0]) & 1;
            if (second0 == bool((cell.takesSecond >> pair[1]) & 1))
              continue;
            // the two edges in the cell across, whose corners lie on the
            // other side along k
            hrleIndex across = idx;
            across[k] += s ? 1 : -1;
            auto acrossEdge = [&](int e) {
              return ((e / D) ^ (1 << k)) * D + e % D;
            };
            unsigned w0, w1;
            if (!vertexOf(across, acrossEdge(pair[0]), w0) ||
                !vertexOf(across, acrossEdge(pair[1]), w1) ||
                (w0 != w1 && across < idx))
              continue;
            std::array<unsigned, 4> polygon{second0 ? cell.second : cell.first,
                                            second0 ? cell.first : cell.second,
                                            w1, w0};
            VectorD n;
            for (int i = 0; i < D; ++i)
              n(i) = outward[pair[0]][i] + outward[pair[1]][i];
            const VectorD p0 = vertexAt(polygon[0]), p1 = vertexAt(polygon[1]),
                          p2 = vertexAt(polygon[2]), p3 = vertexAt(polygon[3]);
            const VectorD area =
                (p1 - p0).cross(p2 - p0) + (p2 - p0).cross(p3 - p0);
            if (area.dot(n) < T(0))
              std::reverse(polygon.begin(), polygon.end());
            faces.push_back(polygon);
            faceNormals.push_back(n);
          }
      }

    // corners of a quad or line that coincide up to rounding, as where
    // neighbouring cells place their vertices at the same point of a sharp
    // edge or corner, become one node; each set of welded nodes keeps its
    // lowest id
    auto &nodes = mesh->nodes;
    std::vector<unsigned> root(nodes.size());
    std::iota(root.begin(), root.end(), 0u);
    auto find = [&](unsigned i) {
      while (root[i] != i)
        i = root[i] = root[root[i]];
      return i;
    };
    const T tolerance = T(100) * std::numeric_limits<T>::epsilon();
    for (const auto &face : faces)
      for (std::size_t i = 0; i < face.size(); ++i)
        for (std::size_t j = i + 1; j < face.size(); ++j) {
          const auto &p = nodes[face[i]], &q = nodes[face[j]];
          T scale = gridDelta, gap = T(0);
          for (int k = 0; k < D; ++k) {
            scale = std::max({scale, std::abs(p[k]), std::abs(q[k])});
            gap = std::max(gap, std::abs(p[k] - q[k]));
          }
          if (gap <= tolerance * scale) {
            const unsigned a = find(face[i]), b = find(face[j]);
            root[std::max(a, b)] = std::min(a, b);
          }
        }
    std::vector<unsigned> newId(nodes.size());
    std::vector<Vec3D<T>> kept;
    for (unsigned i = 0; i < nodes.size(); ++i) {
      const unsigned r = find(i);
      if (r == i) {
        newId[i] = kept.size();
        kept.push_back(nodes[i]);
      } else {
        newId[i] = newId[r];
      }
    }
    nodes = std::move(kept);

    // folds in 2D: a line whose outside faces against the tangent plane of
    // its crossing edge, more than about 96 degrees off, where the vertex of
    // a cell lies beyond that of the cell next to it. The quadric of nearly
    // parallel planes, as two sampled planes of a noisy trench floor next to
    // its corner, puts its point far out, up to the margin beyond the cell,
    // where a sharp corner needs it, and a rebuild from such a mesh takes the
    // fold for a corner of the surface and etches from it. A vertex of a
    // folded line outside its cell whose planes lie within 35 degrees of each
    // other goes to the mass point of its crossing points, in a few passes,
    // unless it stands for several cells, as the welded vertices of a sharp
    // corner do, or is copied across a boundary. In 3D a triangle on the far
    // face of a sharp concave edge faces against the tangent plane of a
    // crossing edge on the near face, so that the same test there takes sharp
    // features for folds
    if constexpr (D == 2) {
      std::vector<int> count(nodes.size(), 0);
      std::vector<long> from(nodes.size(), -1);
      for (unsigned i = 0; i < newId.size(); ++i) {
        ++count[newId[i]];
        from[newId[i]] = i;
      }
      std::vector<char> copied(nodes.size(), 0);
      for (const auto &c : copies) {
        copied[newId[c.id]] = 1;
        copied[newId[c.source]] = 1;
      }
      auto folded = [&](std::size_t k) {
        const unsigned a = newId[faces[k][0]], b = newId[faces[k][1]];
        if (a == b)
          return false;
        const T tx = nodes[b][0] - nodes[a][0], ty = nodes[b][1] - nodes[a][1];
        const VectorD &n = faceNormals[k];
        return -ty * n(0) + tx * n(1) <
               -T(0.1) * std::sqrt(tx * tx + ty * ty) * n.norm();
      };
      for (int pass = 0; pass < 4; ++pass) {
        bool moved = false;
        for (std::size_t k = 0; k < faces.size(); ++k) {
          if (!folded(k))
            continue;
          for (const unsigned original : faces[k]) {
            const unsigned id = newId[original];
            const long v = from[id];
            if (count[id] != 1 || copied[id] || v < 0 ||
                std::size_t(v) >= hasNearlyParallel.size() ||
                !hasNearlyParallel[v])
              continue;
            const auto &[idx, mass] = nearlyParallel[v];
            bool outside = false;
            for (int i = 0; i < D; ++i)
              outside = outside || nodes[id][i] < idx[i] * gridDelta ||
                        nodes[id][i] > (idx[i] + 1) * gridDelta;
            if (!outside)
              continue;
            for (int i = 0; i < D; ++i)
              nodes[id][i] = mass(i) * gridDelta;
            moved = true;
          }
        }
        if (!moved)
          break;
      }
    }

    for (std::size_t k = 0; k < faces.size(); ++k) {
      const auto &face = faces[k];
      if constexpr (D == 3) {
        // the corners of the quad without repeats: a quad with two welded
        // neighbouring corners is a triangle, and one with welded opposite
        // corners or two corners left encloses nothing
        std::array<unsigned, 4> corner;
        int count = 0;
        for (const unsigned id : face)
          if (count == 0 || corner[count - 1] != newId[id])
            corner[count++] = newId[id];
        if (count > 1 && corner[count - 1] == corner[0])
          --count;
        if (count == 3) {
          mesh->insertNextTriangle({corner[0], corner[1], corner[2]});
          continue;
        }
        if (count < 4 || corner[0] == corner[2] || corner[1] == corner[3])
          continue;
        // split along the diagonal 0-2 unless that leaves a triangle without
        // area, as where three vertices lie on one sharp edge, or one that
        // faces against the tangent plane of the crossing edge, as where
        // vertices of neighbouring cells cross over, and the diagonal 1-3
        // does not
        const VectorD &n = faceNormals[k];
        auto poor = [&](unsigned a, unsigned b, unsigned c) {
          const auto &pa = nodes[a], &pb = nodes[b], &pc = nodes[c];
          const Vec3D<T> ab{pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]},
              ac{pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]},
              bc{pc[0] - pb[0], pc[1] - pb[1], pc[2] - pb[2]};
          const T longest = std::max({Norm2(ab), Norm2(ac), Norm2(bc)});
          const auto normal = CrossProduct(ab, ac);
          return Norm2(normal) <= std::max(T(1e-20), roundOff * roundOff) *
                                      longest * longest ||
                 normal[0] * n(0) + normal[1] * n(1) + normal[2] * n(2) < T(0);
        };
        if ((poor(corner[0], corner[1], corner[2]) ||
             poor(corner[0], corner[2], corner[3])) &&
            !poor(corner[0], corner[1], corner[3]) &&
            !poor(corner[1], corner[2], corner[3])) {
          mesh->insertNextTriangle({corner[0], corner[1], corner[3]});
          mesh->insertNextTriangle({corner[1], corner[2], corner[3]});
        } else {
          mesh->insertNextTriangle({corner[0], corner[1], corner[2]});
          mesh->insertNextTriangle({corner[0], corner[2], corner[3]});
        }
      } else {
        const unsigned a = newId[face[0]], b = newId[face[1]];
        if (a == b)
          continue;
        mesh->insertNextLine({a, b});
      }
    }
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(ToDualContouringMesh)

} // namespace viennals
