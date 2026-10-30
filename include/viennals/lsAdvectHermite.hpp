#pragma once

#include <lsPreCompileMacros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

#include <hrleSparseIterator.hpp>
#include <hrleSparseStarIterator.hpp>

#include <lsBooleanOperation.hpp>
#include <lsCheck.hpp>
#include <lsDomain.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMaterialStack.hpp>
#include <lsMesh.hpp>
#include <lsReinitialize.hpp>
#include <lsToDualContouringMesh.hpp>
#include <lsVelocityField.hpp>

namespace viennals {

using namespace viennacore;

/// Advection of a Euclidean level set carrying Hermite data, the closest
/// points and normals stored by the Euclidean FromSurfaceMesh.
///
/// On the side the front moves into, and on static parts, every band node
/// moves its surface point along its characteristic (Huygens' construction).
/// Behind the front it moves its closest face plane with the speed of that
/// face. At edges and corners the fan is the support function over directions
/// sampled in the cone of the feature, planes through the feature point offset
/// by t V(D). Pieces of the front moving at different speeds combine through
/// the envelope of neighbouring planes: the largest signed distance at convex
/// junctions, the smallest at concave ones, with up to three pieces combined
/// exactly at saddle corners. Once the front has moved a given distance, the
/// band is rebuilt exactly from the dual contouring mesh of the planes. Sharp
/// edges and corners stay sharp where the motion keeps them sharp, and fans
/// round them where it rounds them.
///
/// A static convex feature, as the corner of a mask, acts on the gas nodes
/// around it with its exact distance, its support function over the cone of
/// normals of the two faces meeting there (in 3D at an edge), and a node that
/// takes it over becomes a fan of that feature with its convexity. A face
/// meets a fan with the fan's convexity where it is one of the fan's faces,
/// and a convex feature behind a face's plane at the rebuild meets that face
/// convexly, also once the face has moved and cut part of the feature off.
///
/// A convex edge whose two faces grow at different speeds, as a silhouette
/// edge of directional deposition, creates a face between the moved face and
/// the other one. Its fan acts on every node around it with its support
/// function over the whole cone, offset by the time over which the speeds
/// have acted, fitted to the offsets its planes accumulated. This is exact
/// for velocities that change only their rate, as between calls. A node
/// beside the new face takes the plane of that face and becomes a fan of the
/// edge.
///
/// The level set must be Euclidean with closest points and normals, as created
/// by FromSurfaceMesh or MakeGeometry. Level sets below it in a material stack
/// (as in Advect) give the material at each surface point. Where etching
/// reaches a lower material, the step is split at the interface, and where it
/// lays an interface open from below, as under a mask, the interface plane
/// comes from the lower level set. A fan takes the material of its region. A
/// neighbouring face acts on a node where the node's foot on it lies on the
/// face's material, so that the wall of a mask does not reach into the
/// undercut below it, and a fan acts on the gas and on its own material, so
/// that a pit etched into the substrate does not cut into the mask, in either
/// order of the stack. A plane taken over from a neighbour keeps the
/// neighbour's point on it, which gives its material. Behind a junction of
/// materials the interface splits the cone of a fan on the solid side: each
/// material moves over its own part at its own speed, as the substrate below
/// the foot of a mask wall on a slope while the mask beside it stays. Where a
/// face of the cone is of another material, as a wall of a mask leaning over
/// the substrate surface, the cone of the fan's material reaches past it up
/// to the interface and meets it convexly, while the interface itself meets
/// the fan concavely. Scalar and vector velocities. The ends of a fan's cone
/// keep the offsets they accumulated, and in between the fan takes the
/// current velocities.
///
/// After apply() the point data hold the current planes and their priorities
/// for ToDualContouringMesh and the band state, from which a further apply()
/// continues as long as the level set keeps the values it left. Once the
/// whole advection time is done, a band that is no valid level set, since
/// values far from the surface lag behind edges created since the last
/// rebuild, is rebuilt (setRebuildAtEnd). A band narrower than the rebuilt
/// one, as left by Reduce or Prune, is widened first.
template <class T, int D> class AdvectHermite {
  using Vec = Vec3D<T>;
  using hrleIndex = viennahrle::Index<D>;
  static constexpr int numberOfNeighbours = (D == 2) ? 8 : 26;
  static constexpr int slots = 3;
  static constexpr T junction = T(0.98480775301220806);      // cos 10°
  static constexpr T coneTolerance = T(0.99984769515639127); // cos 1°
  static constexpr T bendLimit = T(0.1);
  // accuracy of the planes of the rebuilt mesh beside edges and corners: a
  // distance (grid spacings) plus a slope, and a tilt (sin 2 degrees)
  static constexpr T rebuiltLength = T(0.02);
  static constexpr T rebuiltSlope = T(0.01);
  static constexpr T rebuiltTilt = T(0.034899496702500969);

  // material stack, the last level set is the surface that is advected
  std::vector<SmartPointer<Domain<T, D>>> levelSets;
  SmartPointer<Domain<T, D>> levelSet = nullptr;
  // the lower level sets, for the material at a point
  MaterialStack<T, D> stack;
  SmartPointer<VelocityField<T>> velocities = nullptr;
  T advectionTime = 0.;
  T timeStepRatio = 0.5;  // largest move per step, grid spacings
  T rebuildDistance = 2.; // front motion between rebuilds, grid spacings
  int bandWidth = 12;     // level set width of the rebuilt band
  bool singleStep = false;
  bool rebuildAtEnd = true;
  T advectedTime = 0.;
  unsigned numberOfTimeSteps = 0;
  unsigned numberOfRebuilds = 0;

  // tolerances of the precision: the fixed values with double, scaled to
  // the rounding of float; for unit vectors, for values and lengths in grid
  // spacings, and for times and motions relative to their size
  static constexpr T epsilon = std::numeric_limits<T>::epsilon();
  static constexpr T unitTolerance = std::max(T(1e-9), T(64) * epsilon);
  static constexpr T valueTolerance = std::max(T(1e-9), T(1024) * epsilon);
  static constexpr T smallLength = std::max(T(1e-6), T(1024) * epsilon);
  static constexpr T tinyLength = std::max(T(1e-12), T(1024) * epsilon);
  // a node whose direction from its closest point leaves the normal stored
  // there by more than this (0.08 degrees in double) is a fan of the edge or
  // corner it lies on; float keeps 2.5 degrees, against the rounding of the
  // directions of nodes near the surface
  static constexpr T fanCosine = 1 - std::max(T(1e-6), T(8192) * epsilon);

  // band nodes, indexed by point id; positions in grid spacings
  T gridDelta = 1.;
  std::vector<std::array<long, numberOfNeighbours>> neighbours;
  // period shift of every neighbour across periodic boundaries, as an index
  // into shiftTable
  std::vector<std::array<unsigned char, numberOfNeighbours>> neighbourShift;
  std::vector<Vec> shiftTable;
  std::vector<Vec> X, HA, FA, NRM, NRM0, C0, W, Dir, PL, NL;
  // a point of every node's piece at the rebuild, kept with the piece when
  // another node takes it over: a face's plane then, a fan's feature point
  std::vector<Vec> RW;
  // band nodes sorted by grid index, for the node at a point
  std::vector<std::pair<hrleIndex, long>> nodeOrder;
  std::vector<T> S, side0;
  std::vector<char> good, fan, moved;
  // whether a fan's feature is convex (its fan on the gas side), kept when
  // another node takes the fan over
  std::vector<char> convexFan;
  // cone of every fan whose feature is a 2D corner or a 3D edge: the
  // normals of the two faces meeting there, its direction lying between them
  std::vector<Vec> arcA, arcB;
  std::vector<char> arcKnown;
  // solid nodes that took a moving plane of the etched material around them
  // from exposeInterfaces, as the plane of the end of a cavity below a mask
  // past which they lie; dual contouring prefers it (priority 3) for the
  // crossing edges of their cells it cuts
  std::vector<char> boundPlane;
  // nodes that kept the tangent plane of a convex rounded part at a rebuild
  // (keepAdvected), as on the arc a deposition rounds an edge into
  std::vector<char> smooth;
  T sinceRebuild = 0.; // time since the last rebuild
  T motionSinceRebuild =
      0.; // front motion since the last rebuild, grid spacings

  static T dot(const Vec &a, const Vec &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }
  static T norm(const Vec &a) { return std::sqrt(dot(a, a)); }
  static Vec scaled(const Vec &a, T s) {
    return Vec{a[0] * s, a[1] * s, a[2] * s};
  }
  static Vec sum(const Vec &a, const Vec &b) {
    return Vec{a[0] + b[0], a[1] + b[1], a[2] + b[2]};
  }
  static Vec difference(const Vec &a, const Vec &b) {
    return Vec{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  }
  static Vec cross(const Vec &a, const Vec &b) {
    return Vec{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
               a[0] * b[1] - a[1] * b[0]};
  }
  // whether the unit vector u lies on the shorter arc from a to b
  static bool within(const Vec &u, const Vec &a, const Vec &b) {
    const Vec z = cross(a, b);
    return dot(cross(a, u), z) >= -unitTolerance &&
           dot(cross(u, b), z) >= -unitTolerance && dot(u, sum(a, b)) > 0;
  }
  static T sign(T v) { return T((v > 0) - (v < 0)); }

  // neighbour offsets ordered faces, edges, corners
  static std::vector<hrleIndex> offsets() {
    std::vector<hrleIndex> result;
    const int dz0 = (D == 3) ? -1 : 0, dz1 = (D == 3) ? 1 : 0;
    for (int dz = dz0; dz <= dz1; ++dz)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dx = -1; dx <= 1; ++dx) {
          if (dx == 0 && dy == 0 && dz == 0)
            continue;
          hrleIndex o;
          o[0] = dx;
          o[1] = dy;
          if constexpr (D == 3)
            o[2] = dz;
          result.push_back(o);
        }
    std::stable_sort(result.begin(), result.end(),
                     [](const hrleIndex &a, const hrleIndex &b) {
                       int la = 0, lb = 0;
                       for (int i = 0; i < D; ++i) {
                         la += std::abs(a[i]);
                         lb += std::abs(b[i]);
                       }
                       return la < lb;
                     });
    return result;
  }

  // band node at the grid point nearest to p, -1 where there is none
  long nodeAt(const Vec &p) const {
    const auto &grid = levelSet->getGrid();
    hrleIndex q;
    for (int i = 0; i < D; ++i) {
      q[i] = static_cast<viennahrle::IndexType>(std::lround(p[i]));
      if (grid.isBoundaryPeriodic(i))
        q[i] = grid.globalIndex2LocalIndex(i, q[i]);
    }
    auto found = std::lower_bound(
        nodeOrder.begin(), nodeOrder.end(), q,
        [](const auto &a, const hrleIndex &b) { return a.first < b; });
    return (found != nodeOrder.end() && !(q < found->first)) ? found->second
                                                             : -1;
  }

  // band node at the grid point q, -1 where there is none, and the shift
  // into the frame of q of its geometry across periodic boundaries
  long nodeAtIndex(hrleIndex q, Vec &shift) const {
    const auto &grid = levelSet->getGrid();
    shift = Vec{};
    for (int i = 0; i < D; ++i)
      if (grid.isBoundaryPeriodic(i)) {
        const auto local = grid.globalIndex2LocalIndex(i, q[i]);
        shift[i] = T(q[i] - local);
        q[i] = local;
      }
    auto found = std::lower_bound(
        nodeOrder.begin(), nodeOrder.end(), q,
        [](const auto &a, const hrleIndex &b) { return a.first < b; });
    return (found != nodeOrder.end() && !(q < found->first)) ? found->second
                                                             : -1;
  }

  // the cone of every fan: among the nodes around its feature point, those
  // whose closest point at the rebuild lay on a face whose plane contains the
  // feature point give the normals of the faces meeting there; with exactly
  // two, the fan's direction between them (and across the edge in 3D), the
  // feature is a 2D corner or a 3D edge, otherwise (3D corners) unknown
  void featureCones() {
    const long n = long(S.size());
    arcA.assign(n, Vec{});
    arcB.assign(n, Vec{});
    arcKnown.assign(n, 0);
    constexpr int span = 4;
    constexpr int boxSize = (D == 2) ? span * span : span * span * span;
#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      if (!fan[i])
        continue;
      const Vec &c = C0[i];
      hrleIndex base;
      for (int k = 0; k < D; ++k)
        base[k] = static_cast<viennahrle::IndexType>(std::floor(c[k])) - 1;
      std::array<Vec, 3> normals;
      int found = 0;
      for (int code = 0; code < boxSize && found < 3; ++code) {
        hrleIndex q = base;
        int rest = code;
        for (int k = 0; k < D; ++k) {
          q[k] += rest % span;
          rest /= span;
        }
        Vec shift;
        const long m = nodeAtIndex(q, shift);
        if (m < 0 || fan[m])
          continue;
        const Vec r = difference(X[m], C0[m]);
        const T length = norm(r);
        if (length > smallLength &&
            std::abs(dot(r, NRM0[m])) < T(0.999) * length)
          continue;
        // to the accuracy of the rebuilt mesh, whose faces beside an edge
        // pass a few hundredths of a grid spacing beside the feature point
        // and tilt by a degree or two
        const Vec toFeature = difference(c, sum(C0[m], shift));
        if (std::abs(dot(toFeature, NRM0[m])) >
            rebuiltLength + rebuiltSlope * norm(toFeature))
          continue;
        bool seen = false;
        for (int v = 0; v < found; ++v)
          seen = seen ||
                 dot(normals[v], NRM0[m]) >= 1 - rebuiltTilt * rebuiltTilt / 2;
        if (!seen)
          normals[found++] = NRM0[m];
      }
      if (found != 2 || dot(normals[0], normals[1]) < T(-1) + unitTolerance)
        continue;
      const Vec z = cross(normals[0], normals[1]);
      if (D == 3 && std::abs(dot(Dir[i], z)) > rebuiltTilt * norm(z))
        continue;
      // or past one of them, into the part of the cone of its material
      // beyond a face of another material
      if (!within(Dir[i], normals[0], normals[1])) {
        if (levelSets.size() < 2 || convexFan[i])
          continue;
        Vec a = normals[0], b = normals[1], past, end;
        if (materialCone(C0[i], nodeMaterial(i, HA[i], Dir[i]), a, b, past,
                         end) != 2 ||
            !within(Dir[i], a, b))
          continue;
      }
      arcA[i] = normals[0];
      arcB[i] = normals[1];
      arcKnown[i] = 1;
    }
  }

  // normal speed in grid spacings per unit time of the surface point p (grid
  // spacings) with unit normal n; velocity fields keeping data per point (as
  // visibilities) get the node nearest to p, the point the speed is for
  T speedOf(const Vec &p, const Vec &n, unsigned long id, int material) const {
    const Vec coordinate = scaled(p, gridDelta);
    const long near = nodeAt(p);
    const unsigned long point = (near >= 0) ? (unsigned long)near : id;
    const T scalar =
        velocities->getScalarVelocity(coordinate, material, n, point);
    const Vec vector =
        velocities->getVectorVelocity(coordinate, material, n, point);
    return (scalar + dot(vector, n)) / gridDelta;
  }

  // the material of a surface point is that just behind it, against its
  // normal, so that both sides of an interface are told apart
  int materialBehind(const Vec &p, const Vec &n) const {
    return material(difference(p, scaled(n, T(0.01))));
  }

  T speed(const Vec &p, const Vec &n, unsigned long id) const {
    return speedOf(p, n, id, materialBehind(p, n));
  }

  // the material of node id's piece through p along n. A solid node whose
  // surface point is an edge or corner (a fan on the solid side) lies behind
  // a reflex junction, as below the foot of a mask wall beside the opening,
  // and the front that reaches it is that of its own material: just behind
  // the junction may lie a sliver of another, as where a mask reaches a
  // little below the substrate surface
  int nodeMaterial(long id, const Vec &p, const Vec &n) const {
    if (levelSets.size() > 1 && fan[id] && S[id] < 0)
      return material(X[id]);
    return materialBehind(p, n);
  }

  // the part of the cone from a to b of a fan on the solid side whose
  // material, behind the feature point c against the direction, is m. Behind
  // a junction of materials, as at the foot of a mask wall on a slope, the
  // interface splits the solid: the directions of the other material belong
  // to its own fan, which moves at its own speed (1). Where one face of the
  // cone is of the other material, the fan's material reaches past that face
  // up to the interface, as below the mask beside the foot of its wall when
  // the wall leans over the substrate surface (2, past is the face's normal)
  int materialCone(const Vec &c, int m, Vec &a, Vec &b, Vec &past,
                   Vec &end) const {
    // a direction belongs to m where m lies behind the feature point against
    // it, just behind or a grid spacing behind: a sliver of another material
    // there, as a mask reaching a little below the substrate surface, splits
    // nothing
    constexpr T probe = T(0.1);
    auto inside = [&](const Vec &d) {
      return material(difference(c, scaled(d, probe))) == m ||
             material(difference(c, d)) == m;
    };
    // the other material, just beyond the cut
    auto behind = [&](const Vec &d) {
      const int near = material(difference(c, scaled(d, probe)));
      return (near != m) ? near : material(difference(c, d));
    };
    const bool inA = inside(a), inB = inside(b);
    if (inA != inB) {
      auto direction = [&](T f) {
        const Vec d = sum(scaled(a, 1 - f), scaled(b, f));
        return scaled(d, 1 / norm(d));
      };
      T lo = 0, hi = 1;
      for (int iteration = 0; iteration < 40; ++iteration) {
        const T mid = (lo + hi) / 2;
        if (inside(direction(mid)) == inA)
          lo = mid;
        else
          hi = mid;
      }
      // an interface along the end of the cone, as below a mask on a flat
      // substrate, splits nothing, also within a degree of it: the part left
      // out moves the front by less than the elapsed motion times sin 1°
      Vec cut = inA ? direction(lo) : direction(hi);
      if (!alongInterface(c, m, behind(inA ? direction(hi) : direction(lo)), a,
                          b, cut) ||
          dot(cut, inA ? b : a) > coneTolerance || !within(cut, a, b))
        return 0;
      (inA ? b : a) = cut;
      end = cut;
      return 1;
    }
    if (!inA)
      return 0;
    // the face with normal e lies along the part of the other normal o
    // normal to e
    auto along = [&](const Vec &e, const Vec &o) {
      const Vec t = difference(o, scaled(e, dot(o, e)));
      return scaled(t, 1 / norm(t));
    };
    const bool otherA =
                   materialBehind(sum(c, scaled(along(a, b), probe)), a) != m,
               otherB =
                   materialBehind(sum(c, scaled(along(b, a), probe)), b) != m;
    if (otherA == otherB)
      return 0;
    Vec &e = otherB ? b : a;
    const Vec t = otherB ? along(b, a) : along(a, b);
    // directions turning from e away from the other face
    const Vec start = e;
    auto turned = [&](T phi) {
      return sum(scaled(start, std::cos(phi)), scaled(t, -std::sin(phi)));
    };
    const T quarter = T(1.5707963267948966);
    if (inside(turned(quarter)))
      return 0;
    T lo = 0, hi = quarter;
    for (int iteration = 0; iteration < 40; ++iteration) {
      const T mid = (lo + hi) / 2;
      if (inside(turned(mid)))
        lo = mid;
      else
        hi = mid;
    }
    Vec cut = turned(lo);
    if (!alongInterface(c, m, behind(turned(hi)), a, b, cut) ||
        dot(cut, start) > coneTolerance || dot(cut, t) > 0)
      return 0;
    past = start;
    e = cut;
    end = e;
    return 2;
  }

  // the direction of the interface, in the plane of the cone from a to b,
  // where the probes behind the feature point c found m to end along cut and
  // the material other to begin: the interface a grid spacing behind c,
  // normal to the gradient there of the level set between the two. A step
  // of the interface right at c, as a mask reaching a little below the
  // substrate surface, leaves it as it lies further out
  bool alongInterface(const Vec &c, int m, int other, const Vec &a,
                      const Vec &b, Vec &cut) const {
    const Vec q = difference(c, cut);
    const std::size_t l = std::size_t(std::min(m, other));
    if (other == m || l + 1 >= levelSets.size())
      return false;
    Vec gradient{};
    constexpr T h = T(0.05);
    for (int k = 0; k < D; ++k) {
      Vec step{};
      step[k] = h;
      gradient[k] =
          (lowerValue(l, sum(q, step)) - lowerValue(l, difference(q, step))) /
          (2 * h);
    }
    // in the plane of the cone
    if constexpr (D == 3) {
      Vec z = cross(a, b);
      z = scaled(z, 1 / norm(z));
      gradient = difference(gradient, scaled(z, dot(gradient, z)));
    }
    const T length = norm(gradient);
    if (length < tinyLength)
      return false;
    gradient = scaled(gradient, 1 / length);
    // the direction normal to it nearest the probed one
    Vec along = difference(cut, scaled(gradient, dot(cut, gradient)));
    const T lengthAlong = norm(along);
    if (lengthAlong < T(0.5))
      return false;
    cut = scaled(along, 1 / lengthAlong);
    return true;
  }

  // multilinear value of lower level set l at p (grid spacings)
  T lowerValue(std::size_t l, const Vec &p) const { return stack.value(l, p); }

  // whether lower level set l holds a value at the grid node p
  bool lowerDefined(std::size_t l, const Vec &p) const {
    return stack.defined(l, p);
  }

  // material at p: the lowest level set of the stack that contains it
  // strictly; a point on an interface belongs to the material above
  int material(const Vec &p) const { return stack.material(p); }

  // normal speed of a surface point p moving along n during dt; where it
  // etches into another material, the step is split where it crosses the
  // boundary of the lower one of the two and the rest of it taken with the
  // speed of the new material
  T splitAtInterface(const Vec &p, const Vec &n, T v, T dt,
                     unsigned long id) const {
    const int m0 = nodeMaterial(long(id), p, n);
    const Vec p1 = sum(p, scaled(n, dt * v));
    const int m1 = materialBehind(p1, n);
    if (m1 == m0)
      return v;
    const std::size_t l = std::size_t(std::min(m0, m1));
    const T a = lowerValue(l, p), b = lowerValue(l, p1);
    const T f = (a * b < 0) ? std::clamp(a / (a - b), T(0), T(1)) : T(1);
    const Vec crossing = sum(p, scaled(n, f * dt * v));
    return f * v + (1 - f) * speedOf(crossing, n, id, m1);
  }

  // +1 where two pieces meet convexly, -1 concavely: the normal turns with
  // the path between the witnesses on convex bends; pieces meeting on one
  // edge or corner (fans) count as convex outside the material and concave
  // inside
  static int classify(const Vec &na, const Vec &wa, const Vec &nb,
                      const Vec &wb, T side) {
    const Vec dn = difference(nb, na), dp = difference(wb, wa);
    const T b = dot(dn, dp) / std::max(norm(dn) * norm(dp), T(1e-12));
    if (std::abs(b) <= bendLimit)
      return int(sign(side));
    return (b > bendLimit) ? 1 : -1;
  }

  // shift into the frame of node i of the geometry of its neighbour k
  const Vec &shiftOf(long i, int k) const {
    return shiftTable[neighbourShift[i][k]];
  }

  // band nodes, their values and the neighbour table from the level set;
  // neighbours across periodic boundaries are wrapped and their shift kept
  void collectNodes() {
    auto &domain = levelSet->getDomain();
    const auto &grid = levelSet->getGrid();
    gridDelta = grid.getGridDelta();
    int combinations = 1;
    for (int i = 0; i < D; ++i)
      combinations *= 3;
    shiftTable.assign(combinations, Vec{});
    for (int code = 0; code < combinations; ++code) {
      int rest = code;
      for (int i = 0; i < D; ++i) {
        shiftTable[code][i] =
            T(rest % 3 - 1) * T(grid.getMaxIndex(i) - grid.getMinIndex(i));
        rest /= 3;
      }
    }
    const auto n = levelSet->getNumberOfPoints();
    std::vector<std::pair<hrleIndex, long>> sorted;
    sorted.reserve(n);
    S.assign(n, T(0));
    X.assign(n, Vec{});
    for (viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType> it(
             domain);
         !it.isFinished(); it.next()) {
      if (!it.isDefined())
        continue;
      const auto id = it.getPointId();
      const auto &idx = it.getStartIndices();
      S[id] = it.getValue();
      for (int i = 0; i < D; ++i)
        X[id][i] = T(idx[i]);
      sorted.emplace_back(idx, long(id));
    }
    std::sort(sorted.begin(), sorted.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });
    nodeOrder = sorted;
    const auto offs = offsets();
    neighbours.assign(n, {});
    neighbourShift.assign(n, {});
#pragma omp parallel for
    for (long p = 0; p < long(sorted.size()); ++p) {
      const auto &[idx, id] = sorted[p];
      for (int k = 0; k < numberOfNeighbours; ++k) {
        hrleIndex q = idx;
        int code = 0, power = 1;
        for (int i = 0; i < D; ++i) {
          q[i] += offs[k][i];
          int wrap = 0;
          if (grid.isBoundaryPeriodic(i)) {
            const auto local = grid.globalIndex2LocalIndex(i, q[i]);
            wrap = (q[i] - local) / (grid.getMaxIndex(i) - grid.getMinIndex(i));
            q[i] = local;
          }
          code += (wrap + 1) * power;
          power *= 3;
        }
        neighbourShift[id][k] = static_cast<unsigned char>(code);
        auto found = std::lower_bound(
            sorted.begin(), sorted.end(), q,
            [](const auto &a, const hrleIndex &b) { return a.first < b; });
        neighbours[id][k] = (found != sorted.end() && !(q < found->first) &&
                             !(found->first < q))
                                ? found->second
                                : -1;
      }
    }
  }

  // Hermite data from the closest points and normals of the level set
  bool initializeFromClosestPoints() {
    const auto n = S.size();
    auto &pointData = levelSet->getPointData();
    const auto *closest = pointData.getVectorData(
        FromSurfaceMesh<T, D>::closestPointsLabel, true);
    const auto *normals = pointData.getVectorData(
        FromSurfaceMesh<T, D>::surfaceNormalsLabel, true);
    if (closest == nullptr || normals == nullptr || closest->size() != n ||
        normals->size() != n) {
      VIENNACORE_LOG_ERROR("AdvectHermite: the level set has no closest points "
                           "and normals (Euclidean FromSurfaceMesh).");
      return false;
    }
    HA.assign(n, Vec{});
    NRM.assign(n, Vec{});
    Dir.assign(n, Vec{});
    side0.assign(n, T(0));
    good.assign(n, 0);
    fan.assign(n, 0);
    moved.assign(n, 0);
    convexFan.assign(n, 0);
    boundPlane.assign(n, 0);
    smooth.assign(n, 0);
    if (const auto *kept = pointData.getScalarData(smoothLabel, true);
        kept != nullptr && kept->size() == n)
      for (std::size_t id = 0; id < n; ++id)
        smooth[id] = (*kept)[id] > T(0.5);
#pragma omp parallel for
    for (long id = 0; id < long(n); ++id) {
      Vec cp{}, nrm{};
      for (int i = 0; i < 3; ++i) {
        cp[i] = (*closest)[id][i] / gridDelta;
        nrm[i] = (*normals)[id][i];
      }
      const Vec r = difference(X[id], cp);
      const T length = norm(r);
      const T cosine = dot(r, nrm) / std::max(length, T(1e-12));
      // a folded triangle of the rebuilt mesh can give a node whose closest
      // point lies inside it a normal against its own side
      if (std::abs(cosine) > T(0.999) && sign(cosine) != sign(S[id]))
        nrm = scaled(nrm, T(-1));
      const Vec d =
          (length > smallLength) ? scaled(r, sign(S[id]) / length) : nrm;
      // the material behind the closest point along the node's own direction,
      // also for its face, whose normal at a junction of materials lies along
      // the interface. A solid node whose closest point is an edge or corner
      // lies behind a reflex junction, as below the foot of a mask wall beside
      // the opening, and the front that reaches it is that of its own
      // material: just behind the junction may lie a sliver of another, as
      // where a mask reaches a little below the substrate surface
      const bool behindCorner =
          S[id] < 0 && length > smallLength && std::abs(cosine) < T(0.999);
      const int m = behindCorner ? material(X[id]) : materialBehind(cp, d);
      const T v = speedOf(cp, d, id, m);
      HA[id] = cp;
      NRM[id] = nrm;
      Dir[id] = d;
      side0[id] = sign(S[id]);
      convexFan[id] = S[id] > 0;
      // Huygens planes on the side the front moves into and on static parts,
      // the closest face planes behind the front; a node static only along
      // its own direction whose stored face moves (an edge equidistant from
      // a moving and a static face) carries that face's moved plane
      const bool faceMoves = (v == 0) && speedOf(cp, nrm, id, m) != 0;
      good[id] = ((S[id] * v >= 0) || (v == 0)) && !faceMoves;
      fan[id] =
          good[id] && length > smallLength && std::abs(cosine) < fanCosine;
    }
    FA = HA;
    W = HA;
    C0 = HA;
    PL = HA;
    NRM0 = NRM;
    NL = NRM;
    RW = HA;
    sinceRebuild = 0.;
    motionSinceRebuild = 0.;
    featureCones();
    return true;
  }

  // the band state written by the last apply(), if the level set still holds
  // the values it left
  bool restoreState() {
    const auto n = S.size();
    auto &pointData = levelSet->getPointData();
    std::array<const typename PointData<T>::VectorDataType *, 12> vectors{};
    const char *vectorLabels[12] = {anchorLabel,
                                    facePointLabel,
                                    faceNormalLabel,
                                    rebuildNormalLabel,
                                    featurePointLabel,
                                    witnessLabel,
                                    directionLabel,
                                    FromSurfaceMesh<T, D>::closestPointsLabel,
                                    FromSurfaceMesh<T, D>::surfaceNormalsLabel,
                                    rebuildWitnessLabel,
                                    coneALabel,
                                    coneBLabel};
    for (int k = 0; k < 12; ++k) {
      vectors[k] = pointData.getVectorData(vectorLabels[k], true);
      if (vectors[k] == nullptr || vectors[k]->size() != n)
        return false;
    }
    const auto *values = pointData.getScalarData(valueLabel, true);
    const auto *sides = pointData.getScalarData(sideLabel, true);
    const auto *flags = pointData.getScalarData(flagsLabel, true);
    const auto *times = pointData.getScalarData(timeLabel, true);
    if (values == nullptr || sides == nullptr || flags == nullptr ||
        times == nullptr || values->size() != n || sides->size() != n ||
        flags->size() != n || times->size() < 2)
      return false;
    for (std::size_t i = 0; i < n; ++i)
      if (std::abs((*values)[i] - S[i]) >
          valueTolerance * std::max(T(1), std::abs(S[i])))
        return false;
    std::array<std::vector<Vec> *, 12> targets = {
        &HA, &FA, &NRM, &NRM0, &C0, &W, &Dir, &PL, &NL, &RW, &arcA, &arcB};
    for (int k = 0; k < 12; ++k) {
      targets[k]->assign(n, Vec{});
      // points in grid spacings, normals as they are
      const T scale = (k == 2 || k == 3 || k == 6 || k == 8 || k >= 10)
                          ? T(1)
                          : 1 / gridDelta;
      for (std::size_t i = 0; i < n; ++i)
        for (int c = 0; c < 3; ++c)
          (*targets[k])[i][c] = (*vectors[k])[i][c] * scale;
    }
    side0.assign(n, T(0));
    good.assign(n, 0);
    fan.assign(n, 0);
    moved.assign(n, 0);
    convexFan.assign(n, 0);
    arcKnown.assign(n, 0);
    boundPlane.assign(n, 0);
    smooth.assign(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
      side0[i] = (*sides)[i];
      const int bits = int((*flags)[i]);
      good[i] = bits & 1;
      fan[i] = (bits >> 1) & 1;
      moved[i] = (bits >> 2) & 1;
      convexFan[i] = (bits >> 3) & 1;
      arcKnown[i] = (bits >> 4) & 1;
      boundPlane[i] = (bits >> 5) & 1;
      smooth[i] = (bits >> 6) & 1;
    }
    sinceRebuild = (*times)[0];
    motionSinceRebuild = (*times)[1];
    // the cones as the rebuild found them and the steps passed them on, since
    // the faces around a feature, which featureCones reads, can have become
    // fans of it since then, as nodes that took the plane of a new face
    return true;
  }

  // sample the directions between dx and dy for node i: support planes
  // through the feature point c offset by the elapsed time t times V(D) with
  // the speeds of the fan's material m, corrected linearly so that the two
  // ends keep the offsets ox and oy they accumulated (exact for velocities
  // that keep their shape over t); keep the largest signed distance where the
  // feature is convex, the smallest where it is concave; the plane's point is
  // its support point on the moved fan
  void arc(long i, const Vec dx, T ox, const Vec &dy, T oy, const Vec &c, int m,
           T t, bool convex, bool split, const Vec &ma, const Vec &mb,
           bool refine, T &s, Vec &p, Vec &nn) const {
    static constexpr T samples[] = {0.125, 0.25, 0.375, 0.5,
                                    0.625, 0.75, 0.875, 1.0};
    const T cx = ox - t * speedOf(c, dx, i, m);
    const T cy = oy - t * speedOf(c, dy, i, m);
    // the samples, and the ends of a split cone between dx and dy
    std::array<T, 10> fractions{};
    int count = 0;
    for (T f : samples)
      fractions[count++] = f;
    if (split)
      for (const Vec &e : {ma, mb}) {
        const T ex = norm(cross(dx, e)), ey = norm(cross(dy, e));
        if (ex + ey > tinyLength && within(e, dx, dy))
          fractions[count++] = ex / (ex + ey);
      }
    // the signed distance to the support plane at fraction f, false where
    // the direction is not one of the fan's
    auto planeAt = [&](T f, Vec &d, T &offset, T &value) {
      d = sum(scaled(dx, 1 - f), scaled(dy, f));
      const T length = norm(d);
      if (length < tinyLength)
        return false;
      d = scaled(d, 1 / length);
      // directions of its own material only
      if (split && !within(d, ma, mb))
        return false;
      offset = t * speedOf(c, d, i, m) + (1 - f) * cx + f * cy;
      value = dot(difference(X[i], c), d) - offset;
      return true;
    };
    T best = -1;
    for (int q = 0; q < count; ++q) {
      const T f = fractions[q];
      Vec d;
      T offset, value;
      if (!planeAt(f, d, offset, value))
        continue;
      if (convex ? value > s + tinyLength : value < s - tinyLength) {
        s = value;
        nn = d;
        p = sum(c, scaled(d, offset));
        best = f;
      }
    }
    // a growing fan's largest value lies where the speed has its kink, the
    // direction of the face it creates, which a sample hits only by chance:
    // a golden section search between the samples beside the best one finds
    // it, so that the planes of the new face lie in one plane
    if (!refine || !convex || best < 0)
      return;
    const T g = T(0.61803398874989485);
    T lo = std::max(best - T(0.125), T(0)),
      hi = std::min(best + T(0.125), T(1));
    T x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo);
    Vec d1, d2;
    T o1 = 0, o2 = 0, f1 = -std::numeric_limits<T>::max(),
      f2 = -std::numeric_limits<T>::max();
    if (!planeAt(x1, d1, o1, f1))
      f1 = -std::numeric_limits<T>::max();
    if (!planeAt(x2, d2, o2, f2))
      f2 = -std::numeric_limits<T>::max();
    for (int iteration = 0; iteration < 48; ++iteration) {
      if (f1 < f2) {
        lo = x1;
        x1 = x2;
        f1 = f2;
        d1 = d2;
        o1 = o2;
        x2 = lo + g * (hi - lo);
        if (!planeAt(x2, d2, o2, f2))
          f2 = -std::numeric_limits<T>::max();
      } else {
        hi = x2;
        x2 = x1;
        f2 = f1;
        d2 = d1;
        o2 = o1;
        x1 = hi - g * (hi - lo);
        if (!planeAt(x1, d1, o1, f1))
          f1 = -std::numeric_limits<T>::max();
      }
    }
    const bool first = f1 >= f2;
    const T refined = first ? f1 : f2;
    if (refined >= s) {
      s = refined;
      nn = first ? d1 : d2;
      p = sum(c, scaled(nn, first ? o1 : o2));
    }
  }

  void step(T dt, const std::vector<T> &Vd, const std::vector<T> &Vn) {
    const long n = long(S.size());
    std::vector<Vec> P(n), N(n), Wm(n);
    std::vector<T> Vp(n), Sown(n);
    std::vector<char> moving(n);
#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      moving[i] = good[i] && Vd[i] != 0;
      P[i] = good[i] ? sum(HA[i], scaled(Dir[i], dt * Vd[i]))
                     : sum(FA[i], scaled(NRM[i], dt * Vn[i]));
      N[i] = (good[i] && (moved[i] || Vd[i] != 0)) ? Dir[i] : NRM[i];
      Vp[i] = good[i] ? Vd[i] : Vn[i];
      // moving nodes take the value of their own characteristic, since a
      // fan node may hold the value of a sampled support plane
      Sown[i] =
          moving[i] ? dot(difference(X[i], P[i]), N[i]) : S[i] - dt * Vp[i];
      // witnesses move with their piece; a fan's stays at its feature point
      Wm[i] = fan[i] ? W[i] : sum(W[i], scaled(N[i], dt * Vp[i]));
    }
    sinceRebuild += dt;

    // in a material stack, the material of every piece; a fan's is that of
    // its region, since the ends of its cone may lie on an interface, and on
    // the solid side, behind a reflex junction, that of the node itself
    const bool stacked = levelSets.size() > 1;
    std::vector<int> piece(n, 0), region(n, 0);
    if (stacked) {
#pragma omp parallel for
      for (long i = 0; i < n; ++i) {
        region[i] = material(X[i]);
        piece[i] = fan[i] ? nodeMaterial(i, HA[i], Dir[i])
                          : materialBehind(P[i], N[i]);
      }
    }

    // the part of the cone of a fan on the solid side of its own material
    std::vector<Vec> coneA(arcA), coneB(arcB), past(n, Vec{}), cut(n, Vec{});
    std::vector<char> split(n, 0), extended(n, 0);
    if (stacked) {
#pragma omp parallel for
      for (long i = 0; i < n; ++i)
        if (fan[i] && arcKnown[i] && !convexFan[i]) {
          const int kind = materialCone(C0[i], piece[i], coneA[i], coneB[i],
                                        past[i], cut[i]);
          split[i] = kind != 0;
          extended[i] = kind == 2;
        }
    }

    // a convex edge whose two faces grow at different speeds, as a
    // silhouette edge of directional deposition, creates a face between the
    // moved one and the other. Its fan samples its support function with the
    // time over which the speeds have acted since the rebuild, fitted to the
    // offsets that the fan's own plane and the planes of its faces on
    // neighbouring nodes accumulated over their speeds now: exact where the
    // velocity keeps its shape and only its rate changes, as between calls.
    // Planes much slower than the faster face, as the plane of the new face,
    // leave the time open, and a fan without other planes takes the time of
    // a fan of the same edge next to it
    std::vector<char> arcGrows(n, 0), timed(n, 0);
    std::vector<T> elapsed(n, sinceRebuild);
#pragma omp parallel for
    for (long j = 0; j < n; ++j) {
      if (!fan[j] || !arcKnown[j] || !convexFan[j])
        continue;
      const T va = speedOf(C0[j], arcA[j], j, piece[j]);
      const T vb = speedOf(C0[j], arcB[j], j, piece[j]);
      arcGrows[j] = va >= 0 && vb >= 0 && va != vb;
      if (!arcGrows[j])
        continue;
      const T slow = T(0.1) * std::max(va, vb);
      T ov = 0, vv = 0;
      const T own = speedOf(C0[j], N[j], j, piece[j]);
      if (own > slow) {
        ov += own * dot(difference(P[j], C0[j]), N[j]);
        vv += own * own;
      }
      for (int q = 0; q < numberOfNeighbours; ++q) {
        const long m = neighbours[j][q];
        if (m < 0 || fan[m] || !good[m] || piece[m] != piece[j])
          continue;
        const Vec &sh = shiftOf(j, q);
        if ((dot(NRM0[m], arcA[j]) > 1 - unitTolerance ||
             dot(NRM0[m], arcB[j]) > 1 - unitTolerance) &&
            dot(N[m], NRM0[m]) > 1 - unitTolerance &&
            std::abs(dot(difference(sum(C0[m], sh), C0[j]), NRM0[m])) <=
                smallLength) {
          const T v = speedOf(C0[j], NRM0[m], j, piece[j]);
          if (v > slow) {
            ov += v * dot(difference(sum(P[m], sh), C0[j]), NRM0[m]);
            vv += v * v;
          }
        }
      }
      if (vv > 0) {
        elapsed[j] = ov / vv;
        timed[j] = 1;
      }
    }
    std::vector<long> open, from;
    for (long j = 0; j < n; ++j)
      if (arcGrows[j] && !timed[j])
        open.push_back(j);
    while (!open.empty()) {
      from.assign(open.size(), -1);
#pragma omp parallel for
      for (long o = 0; o < long(open.size()); ++o) {
        const long j = open[o];
        for (int q = 0; q < numberOfNeighbours; ++q) {
          const long m = neighbours[j][q];
          if (m < 0 || !arcGrows[m] || !timed[m])
            continue;
          const Vec dc = difference(sum(C0[m], shiftOf(j, q)), C0[j]);
          const bool sameArc = (dot(arcA[m], arcA[j]) > 1 - unitTolerance &&
                                dot(arcB[m], arcB[j]) > 1 - unitTolerance) ||
                               (dot(arcA[m], arcB[j]) > 1 - unitTolerance &&
                                dot(arcB[m], arcA[j]) > 1 - unitTolerance);
          if (sameArc && std::abs(dot(dc, arcA[j])) <= smallLength &&
              std::abs(dot(dc, arcB[j])) <= smallLength) {
            from[o] = m;
            break;
          }
        }
      }
      std::vector<long> still;
      for (std::size_t o = 0; o < open.size(); ++o)
        if (from[o] < 0)
          still.push_back(open[o]);
        else
          elapsed[open[o]] = elapsed[from[o]];
      for (std::size_t o = 0; o < open.size(); ++o)
        if (from[o] >= 0)
          timed[open[o]] = 1;
      if (still.size() == open.size())
        break;
      open.swap(still);
    }

    // support function of the fans: towards fan neighbours of the same
    // feature and towards the normals of faces whose plane at the rebuild
    // contains the feature point; neighbour by neighbour, each pass reading
    // the directions the previous passes left
    std::vector<T> Sf = Sown;
    std::vector<Vec> Pf = P, Nf = N, previous, previousP;
    for (int k = 0; k < numberOfNeighbours; ++k) {
      previous = Nf;
      previousP = Pf;
#pragma omp parallel for
      for (long i = 0; i < n; ++i) {
        const long j = neighbours[i][k];
        if (!fan[i] || !moving[i] || j < 0)
          continue;
        const bool convex = convexFan[i];
        const T t = elapsed[i];
        const Vec &sh = shiftOf(i, k);
        const Vec dc = difference(sum(C0[j], sh), C0[i]);
        const bool same = fan[j] && std::abs(dot(dc, Nf[i])) <= smallLength &&
                          std::abs(dot(dc, previous[j])) <= smallLength &&
                          dot(Nf[i], previous[j]) < 1 - unitTolerance;
        const bool face = good[j] && !fan[j] &&
                          std::abs(dot(dc, NRM0[j])) <= smallLength &&
                          dot(Nf[i], NRM0[j]) < 1 - unitTolerance;
        // a fan of the same feature of another material moves at another
        // speed, as the static mask behind the foot of its wall
        if (same)
          arc(i, Nf[i], dot(difference(Pf[i], C0[i]), Nf[i]), previous[j],
              (piece[j] == piece[i])
                  ? dot(difference(sum(previousP[j], sh), C0[i]), previous[j])
                  : t * speedOf(C0[i], previous[j], i, piece[i]),
              C0[i], piece[i], t, convex, split[i], coneA[i], coneB[i],
              arcGrows[i], Sf[i], Pf[i], Nf[i]);
        if (face) {
          // the face's accumulated offset while the neighbour still carries
          // that face of the fan's material; one overtaken since holds another
          // plane, a face of another material moves at another speed
          const T oy =
              (dot(N[j], NRM0[j]) > 1 - unitTolerance && piece[j] == piece[i])
                  ? dot(difference(sum(P[j], sh), C0[i]), NRM0[j])
                  : t * speedOf(C0[i], NRM0[j], i, piece[i]);
          arc(i, Nf[i], dot(difference(Pf[i], C0[i]), Nf[i]), NRM0[j], oy,
              C0[i], piece[i], t, convex, split[i], coneA[i], coneB[i],
              arcGrows[i], Sf[i], Pf[i], Nf[i]);
        }
      }
    }

    // speed of the plane each node offers: a moving fan's plane is the support
    // plane in its sampled direction and moves with the speed along that
    // direction, as the static plane of a face swept by a moving neighbour
    std::vector<T> Vplane(Vp);
#pragma omp parallel for
    for (long i = 0; i < n; ++i)
      if (fan[i] && moving[i])
        Vplane[i] = speedOf(C0[i], Nf[i], i, piece[i]);

    // envelope of neighbouring pieces moving at different speeds. A static
    // convex fan whose cone is known acts on the gas nodes around it with its
    // support function there, the distance to the feature point (in 3D to the
    // edge) along the node's direction from it clamped to the cone; a node
    // that takes it becomes a fan of the same feature. Its support plane in
    // the fan's own direction alone underestimates the distance the further
    // the node lies from that direction, as above the middle of a trench
    // beside the corners of its mask. A convex fan whose faces grow at
    // different speeds acts on every node around it with its best support
    // plane over the whole cone at its elapsed time, which beside the face it
    // creates is the plane of that face; a node that takes it becomes a fan of
    // the same feature, and the plane moves with the speed in its direction.
    // Other pieces act with their planes.
    std::vector<T> Snew(Sf);
    std::vector<Vec> Pnew(Pf), Nnew(Nf), Wnew(Wm), Cnew(C0), Anew(arcA),
        Bnew(arcB), RWnew(RW);
    // whether a face is one of the faces of fan f: its plane at the rebuild
    // contains the feature point c and its normal ends the fan's cone
    auto faceOf = [&](const Vec &normal, const Vec &rebuildPoint, const Vec &c,
                      long f) {
      return (dot(normal, arcA[f]) > 1 - unitTolerance ||
              dot(normal, arcB[f]) > 1 - unitTolerance) &&
             std::abs(dot(difference(c, rebuildPoint), normal)) <= smallLength;
    };
    std::vector<char> over(n, 0), fanNew(fan), knownNew(arcKnown),
        convexNew(convexFan);
#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      std::array<Vec, slots> sN, sP, sW0, sWm, sC, sRW;
      std::array<T, slots> sV, sS;
      std::array<long, slots> sFan;
      std::array<int, slots> sPiece{};
      std::array<char, slots> sSupport;
      // the distance of each slot's node to the surface at the step's start
      std::array<T, slots> sDepth;
      const Vec &Ni = Nf[i];
      int count = 0;
      // the growing fans first, whose planes the slots must hold
      for (int pass = 0; pass < 2 * numberOfNeighbours && count < slots;
           ++pass) {
        const int k = pass % numberOfNeighbours;
        const long j = neighbours[i][k];
        if (j < 0 || (pass < numberOfNeighbours) != bool(arcGrows[j]))
          continue;
        const Vec &sh = shiftOf(i, k);
        const Vec &Nj = Nf[j];
        const Vec Pj = sum(Pf[j], sh);
        // a face acts where the node's foot on it lies on its material, as
        // beside the wall of a mask and not in the undercut below it; a fan
        // on the gas and on its own material, as the pit in a substrate does
        // not cut into the mask
        if (stacked) {
          bool valid;
          if (fan[j]) {
            valid = Sf[i] > 0 || region[i] == piece[j];
          } else {
            const Vec foot =
                difference(X[i], scaled(Nj, dot(difference(X[i], Pj), Nj)));
            valid = materialBehind(foot, Nj) == piece[j];
          }
          if (!valid)
            continue;
        }
        T value = dot(difference(X[i], Pj), Nj);
        Vec point = Pj, normal = Nj, witness = sum(W[j], sh),
            movedWitness = sum(Wm[j], sh), rebuildPoint = sum(RW[j], sh),
            feature = sum(C0[j], sh);
        bool support = false;
        T supportSpeed = 0;
        if (arcGrows[j]) {
          const Vec c = sum(C0[j], sh);
          // signed distance to the support plane at f along the cone
          auto planeAt = [&](T f, Vec &u, T &offset) {
            u = sum(scaled(arcA[j], 1 - f), scaled(arcB[j], f));
            const T length = norm(u);
            if (length < tinyLength)
              return -std::numeric_limits<T>::max();
            u = scaled(u, 1 / length);
            offset = elapsed[j] * speedOf(c, u, j, piece[j]);
            return dot(difference(X[i], c), u) - offset;
          };
          constexpr int samples = 32;
          T best = -std::numeric_limits<T>::max(), bestOffset = 0;
          Vec bestU = Nj;
          int bestS = 0;
          for (int s = 0; s <= samples; ++s) {
            Vec u;
            T offset = 0;
            const T v = planeAt(T(s) / samples, u, offset);
            if (v > best) {
              best = v;
              bestU = u;
              bestOffset = offset;
              bestS = s;
            }
          }
          // a golden section search between the samples beside the best one
          // finds the direction where the speed has its kink, the normal of
          // the new face, which a sample hits only by chance: a face off by a
          // fraction of a degree moves, and its rebuilt mesh is no plane
          const T g = T(0.61803398874989485);
          T lo = T(std::max(bestS - 1, 0)) / samples,
            hi = T(std::min(bestS + 1, samples)) / samples;
          T x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo), o1 = 0, o2 = 0;
          Vec u1, u2;
          T f1 = planeAt(x1, u1, o1), f2 = planeAt(x2, u2, o2);
          for (int iteration = 0; iteration < 48; ++iteration) {
            if (f1 < f2) {
              lo = x1;
              x1 = x2;
              f1 = f2;
              u1 = u2;
              o1 = o2;
              x2 = lo + g * (hi - lo);
              f2 = planeAt(x2, u2, o2);
            } else {
              hi = x2;
              x2 = x1;
              f2 = f1;
              u2 = u1;
              o2 = o1;
              x1 = hi - g * (hi - lo);
              f1 = planeAt(x1, u1, o1);
            }
          }
          if (std::max(f1, f2) > best) {
            best = std::max(f1, f2);
            bestU = (f1 >= f2) ? u1 : u2;
            bestOffset = (f1 >= f2) ? o1 : o2;
          }
          value = best;
          normal = bestU;
          point = sum(c, scaled(bestU, bestOffset));
          witness = movedWitness = rebuildPoint = feature = c;
          supportSpeed = speedOf(c, bestU, j, piece[j]);
          support = true;
        } else if (fan[j] && arcKnown[j] && convexFan[j] && Sf[i] > 0) {
          Vec c = sum(C0[j], sh);
          Vec r = difference(X[i], c);
          if constexpr (D == 3) {
            Vec e = cross(arcA[j], arcB[j]);
            e = scaled(e, 1 / norm(e));
            const T along = dot(r, e);
            c = sum(c, scaled(e, along));
            r = difference(r, scaled(e, along));
          }
          const T length = norm(r);
          if (length > smallLength) {
            Vec u = scaled(r, 1 / length);
            if (!within(u, arcA[j], arcB[j]))
              u = (dot(u, arcA[j]) >= dot(u, arcB[j])) ? arcA[j] : arcB[j];
            // static features only, whose support function needs no offsets
            // in other directions than the fan's own
            if (std::abs(dot(difference(Pj, c), Nj)) <= tinyLength &&
                speedOf(c, Nj, i, piece[j]) == 0 &&
                speedOf(c, u, i, piece[j]) == 0) {
              value = dot(r, u);
              point = c;
              normal = u;
              witness = movedWitness = rebuildPoint = feature = c;
              support = true;
            }
          }
        }
        // a moving concave fan, as at a corner of an etched trench bottom,
        // acts where its arc is the nearest part of the surface: on nodes
        // whose direction from its feature, reversed, lies in its cone.
        // Elsewhere its plane would reach nodes behind the faces next to it,
        // as below the floor beside the corner. Convex fans offer their
        // planes beyond the cone, for the faces the motion creates.
        bool beyond = false;
        if (fan[j] && arcKnown[j] && !convexFan[j]) {
          const Vec c = sum(C0[j], sh);
          Vec r = difference(X[i], c);
          if constexpr (D == 3) {
            Vec e = cross(arcA[j], arcB[j]);
            e = scaled(e, 1 / norm(e));
            r = difference(r, scaled(e, dot(r, e)));
          }
          const T length = norm(r);
          const Vec u = scaled(r, T(-1) / std::max(length, tinyLength));
          if (length > smallLength && !within(u, coneA[j], coneB[j]))
            continue;
          // where the cone reaches past a face of another material, the fan
          // acts on the nodes there with its support plane in their own
          // direction, offset as the fan's own plane: beside the face, as
          // the static wall of a mask, the two meet convexly, also where
          // their normals lie close
          if (extended[j] && length > smallLength &&
              !within(u, arcA[j], arcB[j])) {
            const T t = elapsed[j];
            const T offset = t * speedOf(c, u, j, piece[j]) +
                             dot(difference(Pj, c), Nj) -
                             t * speedOf(c, Nj, j, piece[j]);
            value = dot(difference(X[i], c), u) - offset;
            normal = u;
            point = sum(c, scaled(u, offset));
            witness = movedWitness = rebuildPoint = feature = c;
            supportSpeed = speedOf(c, u, j, piece[j]);
            support = true;
            beyond = true;
          }
        }
        if (!beyond && dot(Ni, normal) >= junction)
          continue;
        bool distinct = true;
        for (int q = 0; q < count; ++q)
          distinct = distinct && dot(sN[q], normal) < junction;
        if (!distinct)
          continue;
        sN[count] = normal;
        sP[count] = point;
        sW0[count] = witness;
        sWm[count] = movedWitness;
        sC[count] = feature;
        sRW[count] = rebuildPoint;
        sFan[count] = fan[j] ? j : -1;
        sPiece[count] = piece[j];
        sSupport[count] = support;
        sDepth[count] = std::abs(S[j]);
        sV[count] = support ? supportSpeed : Vplane[j];
        sS[count] = value;
        ++count;
      }
      if (count == 0)
        continue;
      const T side = side0[i];
      const long ownFan = fan[i] ? i : -1;
      // a fan and a face meet with the fan's convexity where the face is one
      // of its faces or one of its support planes, through its feature point
      // at the rebuild with the normal in its cone, as a plane of a growing
      // fan that a node took in a direction between the faces. They meet
      // convexly where a convex feature lies behind the face's plane at the
      // rebuild and concavely where it lies in front of it. This holds also
      // once a face moving away from a static feature has cut part of it
      // off, where the witnesses make the junction look concave. Other pairs
      // are classified from their witnesses.
      // (to the accuracy of the planes of the rebuilt mesh, which pass a
      // few hundredths of a grid spacing beside the feature point and tilt
      // by up to two degrees beside edges and corners, so that such a plane
      // of a face, as the moving face of a growing edge, does not meet the
      // fan concavely)
      auto supportOf = [&](const Vec &normal, const Vec &rebuildPoint,
                           const Vec &c, long f) {
        const Vec r = difference(c, rebuildPoint);
        if (std::abs(dot(r, normal)) > rebuiltLength + rebuiltSlope * norm(r))
          return false;
        if constexpr (D == 3) {
          const Vec z = cross(arcA[f], arcB[f]);
          if (std::abs(dot(normal, z)) > rebuiltTilt * norm(z))
            return false;
        }
        return within(normal, arcA[f], arcB[f]) ||
               dot(normal, arcA[f]) >= 1 - rebuiltTilt * rebuiltTilt / 2 ||
               dot(normal, arcB[f]) >= 1 - rebuiltTilt * rebuiltTilt / 2;
      };
      auto fanAndFace = [&](long f, const Vec &c, const Vec &normal,
                            const Vec &rebuildPoint, int faceMaterial,
                            const Vec &point) {
        // faces of another material through the feature point of a fan
        // whose cone an interface splits: the interface itself, as the
        // underside of a mask that an undercut lays open, bounds the fan's
        // front concavely, and the face that the cone of the fan's material
        // reaches past, as the wall of the mask, convexly
        if (split[f] && faceMaterial != piece[f] &&
            std::abs(dot(difference(c, rebuildPoint), normal)) <= smallLength) {
          if (extended[f] && dot(normal, past[f]) > 1 - unitTolerance)
            return 1;
          if (std::abs(dot(normal, cut[f])) < T(1e-3))
            return -1;
        }
        if (arcKnown[f] && (faceOf(normal, rebuildPoint, c, f) ||
                            supportOf(normal, rebuildPoint, c, f)))
          return convexFan[f] ? 1 : -1;
        // a concave fan, as the corner of an etched trench bottom, is a
        // corner of the convex region on the other side: a plane facing into
        // its cone that has moved past the feature point, as the floor going
        // down beside the corner, bounds that region there and meets the fan
        // concavely
        if (arcKnown[f] && !convexFan[f] && within(normal, arcA[f], arcB[f]) &&
            dot(difference(c, point), normal) > smallLength)
          return -1;
        if (convexFan[f]) {
          const T d = dot(difference(c, rebuildPoint), normal);
          if (d < -smallLength)
            return 1;
          if (d > smallLength)
            return -1;
        }
        return 0;
      };
      auto relation = [&](const Vec &na, const Vec &wa, const Vec &ra, long fa,
                          const Vec &ca, int ma, const Vec &pa, const Vec &nb,
                          const Vec &wb, const Vec &rb, long fb, const Vec &cb,
                          int mb, const Vec &pb) {
        int r = 0;
        if (fa >= 0 && fb < 0)
          r = fanAndFace(fa, ca, nb, rb, mb, pb);
        else if (fb >= 0 && fa < 0)
          r = fanAndFace(fb, cb, na, ra, ma, pa);
        return (r != 0) ? r : classify(na, wa, nb, wb, side);
      };
      std::array<int, slots> withOwn{};
      std::array<std::array<int, slots>, slots> between{};
      for (int k = 0; k < count; ++k) {
        withOwn[k] =
            relation(Ni, W[i], RW[i], ownFan, C0[i], piece[i], Pf[i], sN[k],
                     sW0[k], sRW[k], sFan[k], sC[k], sPiece[k], sP[k]);
        for (int q = 0; q < count; ++q)
          if (q != k)
            between[k][q] = relation(sN[k], sW0[k], sRW[k], sFan[k], sC[k],
                                     sPiece[k], sP[k], sN[q], sW0[q], sRW[q],
                                     sFan[q], sC[q], sPiece[q], sP[q]);
      }
      // a convex partner raises the value after the minimum with partners
      // concave to it, a concave partner lowers it after the maximum with
      // partners convex to it: the exact combination of three half-spaces.
      // Within one step a neighbouring piece changes the value by at most dt
      // times the speed difference, larger changes come from folds of the
      // mesh or from planes that lag behind, as deep in the solid behind the
      // corner of an etched trench bottom: the piece that would change it by
      // more drops out, and the others combine without it, as the wall
      // beside that corner, which keeps a node beyond the end of the arc's
      // plane in the solid
      std::array<char, slots> dropped{};
      T value = Sf[i];
      int active = -1;
      for (int attempt = 0; attempt < count; ++attempt) {
        T up = -std::numeric_limits<T>::max(),
          down = std::numeric_limits<T>::max();
        int upK = -1, downK = -1;
        for (int k = 0; k < count; ++k) {
          if (dropped[k])
            continue;
          T raiseV = sS[k], lowerV = sS[k];
          int raiseK = k, lowerK = k;
          for (int q = 0; q < count; ++q) {
            if (q == k || dropped[q])
              continue;
            if (between[k][q] == -1 && sS[q] < raiseV) {
              raiseV = sS[q];
              raiseK = q;
            }
            if (between[k][q] == 1 && sS[q] > lowerV) {
              lowerV = sS[q];
              lowerK = q;
            }
          }
          if (withOwn[k] == 1 && raiseV > up) {
            up = raiseV;
            upK = raiseK;
          }
          if (withOwn[k] == -1 && lowerV < down) {
            down = lowerV;
            downK = lowerK;
          }
        }
        const bool raised = up > Sf[i] + valueTolerance;
        value = raised ? up : Sf[i];
        const bool lowered = down < value - valueTolerance;
        if (lowered)
          value = down;
        active = lowered ? downK : (raised ? upK : -1);
        if (active < 0)
          break;
        T budget = 0.;
        for (int k = 0; k < count; ++k)
          if (!dropped[k])
            budget = std::max(budget, dt * std::abs(sV[k] - Vplane[i]));
        // where the support plane of a growing fan gives the value, it may
        // have run ahead of the node over the whole time since the rebuild,
        // as beside a face the motion creates, whose plane the node's own
        // plane passes
        if (sSupport[active] && sFan[active] >= 0 && arcGrows[sFan[active]])
          budget = std::max(budget, (elapsed[sFan[active]] + dt) *
                                        std::abs(sV[active] - Vplane[i]));
        // where its own plane would carry a solid node into the gas, a piece
        // of a node at most a grid spacing deeper that keeps it solid holds
        // it there beyond the budget, as the arc at the corner of an etched
        // trench bottom, which reaches a node beside the wall a step after
        // the node above took it, once the floor has passed: values lag
        // behind only deeper in the solid
        const bool kept = S[i] < 0 && Sf[i] >= 0 && value < 0 &&
                          sDepth[active] <= std::abs(S[i]) + T(1);
        if (std::abs(value - Sf[i]) <= budget + smallLength || kept)
          break;
        dropped[active] = 1;
        active = -1;
      }
      if (active < 0)
        continue;
      // the plane keeps the neighbour's point on it, a point of the piece
      // itself, which gives its material and speed; the node's own foot point
      // can lie where the piece does not exist, as in another material
      Snew[i] = value;
      Nnew[i] = sN[active];
      Pnew[i] = sP[active];
      Wnew[i] = sWm[active];
      RWnew[i] = sRW[active];
      over[i] = 1;
      const long source = sSupport[active] ? sFan[active] : -1;
      fanNew[i] = source >= 0;
      if (source >= 0) {
        Cnew[i] = sC[active];
        Anew[i] = arcA[source];
        Bnew[i] = arcB[source];
        knownNew[i] = 1;
        convexNew[i] = convexFan[source];
      }
    }

#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      if (over[i])
        HA[i] = Pnew[i];
      else if (good[i])
        HA[i] = sum(HA[i], scaled(Dir[i], dt * Vd[i]));
      if (!good[i] || over[i])
        FA[i] = Pnew[i];
      if (over[i]) {
        NRM[i] = Nnew[i];
        Dir[i] = Nnew[i];
        C0[i] = Cnew[i];
        boundPlane[i] = 0;
        smooth[i] = 0;
      }
      W[i] = Wnew[i];
      RW[i] = RWnew[i];
      PL[i] = Pnew[i];
      NL[i] = Nnew[i];
      moved[i] = moved[i] || moving[i] || over[i];
      fan[i] = fanNew[i];
      arcA[i] = Anew[i];
      arcB[i] = Bnew[i];
      arcKnown[i] = knownNew[i];
      convexFan[i] = convexNew[i];
      good[i] = good[i] || over[i];
      S[i] = Snew[i];
    }
  }

  // interfaces of the material stack laid open since the last rebuild, as
  // when etching undercuts a mask: across the interface of a lower level set
  // from a solid node, a node etched since then takes the smaller and the
  // solid node the larger signed distance to the plane of the interface,
  // from the lower level set (exact for planar interfaces), and both carry
  // that plane; only where the material on the solid side does not etch,
  // since a front etching through the interface is not held there. Beyond
  // the end of the cavity, as past the tip of an undercut, the interface is
  // no surface: a solid node there lies no closer to the surface than the
  // plane where the front meets the interface, from the fans around it whose
  // cone the interface splits, and carries that plane
  void exposeInterfaces() {
    const long n = long(S.size());
    for (std::size_t l = 0; l + 1 < levelSets.size(); ++l) {
      // from the state before the pass, applied after it
      std::vector<char> take(n, 0), takeBound(n, 0);
      std::vector<T> newValue(n, T(0));
      std::vector<Vec> newPoint(n), newNormal(n);
#pragma omp parallel for
      for (long i = 0; i < n; ++i) {
        const bool etched = S[i] > 0 && side0[i] < 0, solid = S[i] < 0;
        if (!etched && !solid)
          continue;
        const T li = lowerValue(l, X[i]);
        if (li == T(0))
          continue;
        bool exposed = false;
        for (int k = 0; k < numberOfNeighbours && !exposed; ++k) {
          const long j = neighbours[i][k];
          if (j < 0)
            continue;
          const T lj = lowerValue(l, sum(X[j], shiftOf(i, k)));
          if (li * lj >= T(0))
            continue;
          exposed = etched ? S[j] < 0 : (S[j] > 0 && side0[j] < 0);
        }
        if (!exposed || !lowerDefined(l, X[i]))
          continue;
        // the interface plane: from the Hermite data of a Euclidean lower
        // level set, otherwise from differences within its band, whose values
        // further out only give the side
        Vec gradient{};
        Vec planePoint, planeDirection;
        const bool hermite =
            stack.nodePlane(l, X[i], planePoint, planeDirection);
        if (hermite) {
          gradient = planeDirection;
        } else {
          for (int d = 0; d < D; ++d) {
            Vec step{};
            step[d] = 1;
            const Vec forward = sum(X[i], step),
                      backward = difference(X[i], step);
            const bool f = lowerDefined(l, forward),
                       b = lowerDefined(l, backward);
            if (f == b)
              gradient[d] =
                  (lowerValue(l, forward) - lowerValue(l, backward)) / 2;
            else
              gradient[d] = f ? lowerValue(l, forward) - li
                              : li - lowerValue(l, backward);
          }
        }
        const T length = norm(gradient);
        if (length < tinyLength)
          continue;
        // the normal points to the gas, on the side of the etched node
        const T towardsGas = etched ? sign(li) : -sign(li);
        const Vec normal = scaled(gradient, towardsGas / length);
        const T distance = towardsGas * li / length;
        T value = distance;
        bool bounded = false;
        Vec point = difference(X[i], scaled(normal, distance)),
            planeNormal = normal;
        if (solid) {
          // the moving planes of the other material around the node, faces
          // and the ends of the cones of fans on the solid side
          const int own = material(X[i]);
          // the interface counts where the node's foot on it lies in the
          // cavity, on the gas side of every plane of the front
          const Vec interfaceFoot = point;
          T bound = std::numeric_limits<T>::max(), footSide = bound;
          Vec boundNormal{}, endNormal{};
          auto consider = [&](const Vec &p, const Vec &d) {
            const T g = dot(difference(X[i], p), d);
            const T side = dot(difference(interfaceFoot, p), d);
            if (side < footSide) {
              footSide = side;
              endNormal = d;
            }
            if (g < bound) {
              bound = g;
              boundNormal = d;
            }
          };
          for (int k = 0; k < numberOfNeighbours; ++k) {
            const long j = neighbours[i][k];
            if (j < 0)
              continue;
            const Vec &sh = shiftOf(i, k);
            const int m = fan[j] ? nodeMaterial(j, HA[j], Dir[j])
                                 : materialBehind(PL[j], NL[j]);
            if (m == own || speedOf(PL[j], NL[j], j, m) == 0)
              continue;
            consider(sum(PL[j], sh), NL[j]);
            if (!fan[j] || !arcKnown[j] || convexFan[j])
              continue;
            const Vec &c = C0[j];
            const T correction = dot(difference(PL[j], c), NL[j]) -
                                 sinceRebuild * speedOf(c, NL[j], j, m);
            Vec a = arcA[j], b = arcB[j], past, end;
            const bool split = materialCone(c, m, a, b, past, end) != 0;
            for (const Vec &d : {arcA[j], arcB[j], end}) {
              if (&d == &end && !split)
                continue;
              const T offset = sinceRebuild * speedOf(c, d, j, m) + correction;
              consider(sum(sum(c, scaled(d, offset)), sh), d);
            }
          }
          // beyond the end of the cavity the node lies nearest to the edge
          // where the plane of the end meets the interface
          const Vec along =
              difference(endNormal, scaled(normal, dot(endNormal, normal)));
          const bool beyond = footSide < -smallLength && norm(along) > T(0.5);
          if (beyond || bound < value) {
            bounded = true;
            if (beyond) {
              const T past = -footSide / norm(along);
              value = -std::sqrt(distance * distance + past * past);
              bound =
                  dot(difference(X[i], interfaceFoot), endNormal) + footSide;
              boundNormal = endNormal;
            } else {
              value = bound;
            }
            planeNormal = boundNormal;
            // its point is the node's foot on the plane, moved along the plane
            // to the solid side of the interface, which gives its material
            Vec foot = difference(X[i], scaled(planeNormal, bound));
            Vec up = sum(scaled(normal, T(-1)),
                         scaled(planeNormal, dot(normal, planeNormal)));
            const T lengthUp = norm(up);
            if (lengthUp > T(0.5))
              foot = sum(foot, scaled(up, T(0.1) / lengthUp));
            point = foot;
          }
        }
        if (etched ? value >= S[i] : value <= S[i])
          continue;
        const Vec foot = difference(X[i], scaled(normal, distance));
        if (speedOf(foot, normal, i, materialBehind(foot, normal)) < T(0))
          continue;
        take[i] = 1;
        takeBound[i] = solid && bounded;
        newValue[i] = value;
        newPoint[i] = point;
        newNormal[i] = planeNormal;
      }
      for (long i = 0; i < n; ++i) {
        if (!take[i])
          continue;
        S[i] = newValue[i];
        boundPlane[i] = takeBound[i];
        HA[i] = FA[i] = W[i] = PL[i] = RW[i] = newPoint[i];
        NRM[i] = Dir[i] = NL[i] = newNormal[i];
        good[i] = moved[i] = 1;
        fan[i] = 0;
      }
    }
  }

  // values into the level set, current planes and priorities into its point
  // data for ToDualContouringMesh
  void writeBack() {
    auto &domain = levelSet->getDomain();
    for (viennahrle::SparseIterator<typename Domain<T, D>::DomainType> it(
             domain);
         !it.isFinished(); it.next())
      if (it.isDefined())
        it.getDefinedValue() = S[it.getPointId()];
    typename PointData<T>::VectorDataType points(S.size()), normals(S.size());
    typename PointData<T>::ScalarDataType priority(S.size()),
        smoothData(S.size());
    const std::vector<char> rounded = roundedNodes();
    for (std::size_t i = 0; i < S.size(); ++i) {
      for (int k = 0; k < 3; ++k) {
        points[i][k] = PL[i][k] * gridDelta;
        normals[i][k] = NL[i][k];
      }
      // all planes count equally, the one cutting an edge nearest to the
      // interpolated crossing wins; priorities for the planes behind the front
      // pick tilted fan planes on faces next to collapsed fans
      priority[i] = boundPlane[i] ? T(3) : T(2);
      smoothData[i] = rounded[i] ? T(1) : T(0);
    }
    auto &pointData = levelSet->getPointData();
    pointData.insertReplaceVectorData(
        std::move(points), FromSurfaceMesh<T, D>::closestPointsLabel);
    pointData.insertReplaceVectorData(
        std::move(normals), FromSurfaceMesh<T, D>::surfaceNormalsLabel);
    pointData.insertReplaceScalarData(
        std::move(priority), ToDualContouringMesh<T, D>::priorityLabel);
    pointData.insertReplaceScalarData(std::move(smoothData), smoothLabel);

    // the band state, so that a further apply() continues from it
    const std::array<const std::vector<Vec> *, 10> sources = {
        &HA, &FA, &NRM, &NRM0, &C0, &W, &Dir, &RW, &arcA, &arcB};
    const char *labels[10] = {
        anchorLabel,       facePointLabel, faceNormalLabel, rebuildNormalLabel,
        featurePointLabel, witnessLabel,   directionLabel,  rebuildWitnessLabel,
        coneALabel,        coneBLabel};
    for (int k = 0; k < 10; ++k) {
      const T scale = (k == 2 || k == 3 || k == 6 || k >= 8) ? T(1) : gridDelta;
      typename PointData<T>::VectorDataType data(S.size());
      for (std::size_t i = 0; i < S.size(); ++i)
        for (int c = 0; c < 3; ++c)
          data[i][c] = (*sources[k])[i][c] * scale;
      pointData.insertReplaceVectorData(std::move(data), labels[k]);
    }
    typename PointData<T>::ScalarDataType values(S), sides(side0),
        flags(S.size()), times = {sinceRebuild, motionSinceRebuild};
    for (std::size_t i = 0; i < S.size(); ++i)
      flags[i] = T(int(good[i]) + 2 * int(fan[i]) + 4 * int(moved[i]) +
                   8 * int(convexFan[i]) + 16 * int(arcKnown[i]) +
                   32 * int(boundPlane[i]) + 64 * int(smooth[i]));
    pointData.insertReplaceScalarData(std::move(values), valueLabel);
    pointData.insertReplaceScalarData(std::move(sides), sideLabel);
    pointData.insertReplaceScalarData(std::move(flags), flagsLabel);
    pointData.insertReplaceScalarData(std::move(times), timeLabel);
  }

  // a value further from the surface than a neighbour's by more than a grid
  // spacing lags behind it, as a node far out past a rounded part that the
  // last rebuild placed by the dual contouring mesh beside a neighbour that
  // kept its exact value: a distance changes by at most the distance moved,
  // and such a value takes that bound. Only beyond two grid spacings, where
  // no cell crossed by the surface has its corners. Returns whether any
  // value changed
  bool boundFarValues() {
    bool changed = false;
    for (bool again = true; again;) {
      again = false;
      std::vector<std::pair<std::size_t, T>> bounded;
      for (viennahrle::ConstSparseStarIterator<
               typename Domain<T, D>::DomainType, 1>
               it(levelSet->getDomain());
           !it.isFinished(); it.next()) {
        const auto &center = it.getCenter();
        if (!center.isDefined())
          continue;
        const std::size_t id = center.getPointId();
        const T value = S[id];
        if (std::abs(value) <= T(2))
          continue;
        T bound = std::abs(value);
        for (int k = 0; k < 2 * D; ++k) {
          const auto &neighbour = it.getNeighbor(k);
          if (!neighbour.isDefined())
            continue;
          const T other = S[neighbour.getPointId()];
          if ((other < T(0)) == (value < T(0)))
            bound = std::min(bound, std::abs(other) + T(1));
        }
        if (bound < std::abs(value))
          bounded.push_back({id, value < T(0) ? -bound : bound});
      }
      for (const auto &[id, value] : bounded)
        S[id] = value;
      again = !bounded.empty();
      changed = changed || again;
    }
    return changed;
  }

  // exact band from the dual contouring mesh of the current planes
  void rebuild() {
    writeBack();
    const auto before = nodeOrder;
    const std::vector<T> valueBefore = S;
    const std::vector<Vec> normalBefore = NL;
    const std::vector<char> rounded = roundedNodes();
    Reinitialize<T, D>(levelSet, bandWidth).apply();
    keepAdvected(before, valueBefore, normalBefore, rounded);
    ++numberOfRebuilds;
  }

  // nodes whose plane is a tangent plane of a convex rounded part: convex
  // fans whose speed, the same in every direction of their cone, at the node
  // and at its plane, moves the surface outwards, as a deposition rounds an
  // edge into an arc, and nodes that kept such a plane at a rebuild. A rate
  // that varies in space, as a flux or a noisy view factor, or with the
  // direction, as directional deposition, rounds nothing exactly
  std::vector<char> roundedNodes() const {
    const long n = long(S.size());
    std::vector<char> rounded(n, 0);
#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      if (long(smooth.size()) == n && smooth[i]) {
        rounded[i] = 1;
        continue;
      }
      if (!fan[i] || !good[i] || !convexFan[i])
        continue;
      const int m = nodeMaterial(i, HA[i], Dir[i]);
      const T v = speedOf(C0[i], NL[i], i, m);
      if (v <= T(0))
        continue;
      auto same = [&](const Vec &p, const Vec &d) {
        return std::abs(speedOf(p, d, i, m) - v) <= T(1e-6) * v;
      };
      bool uniform = same(C0[i], NRM0[i]) && same(C0[i], Dir[i]) &&
                     same(X[i], NL[i]) && same(PL[i], NL[i]);
      if (arcKnown[i])
        uniform = uniform && same(C0[i], arcA[i]) && same(C0[i], arcB[i]);
      rounded[i] = uniform;
    }
    return rounded;
  }

  // the rebuilt band takes its values from the dual contouring mesh, whose
  // vertices on a convex curved part, as the arc a deposition rounds an edge
  // into, lie where tangent planes meet, outside the arc: the rounding grows
  // on from there, and every rebuild adds to it. A node of such a part keeps
  // its advected value and tangent plane, exact under a uniform rate, where
  // they agree with the rebuilt ones to a quarter of a grid spacing, and is
  // marked for dual contouring (smoothLabel)
  void keepAdvected(const std::vector<std::pair<hrleIndex, long>> &before,
                    const std::vector<T> &valueBefore,
                    const std::vector<Vec> &normalBefore,
                    const std::vector<char> &rounded) {
    auto &pointData = levelSet->getPointData();
    auto *closest = pointData.getVectorData(
        FromSurfaceMesh<T, D>::closestPointsLabel, true);
    auto *normals = pointData.getVectorData(
        FromSurfaceMesh<T, D>::surfaceNormalsLabel, true);
    if (closest == nullptr || normals == nullptr ||
        closest->size() != levelSet->getNumberOfPoints() ||
        normals->size() != levelSet->getNumberOfPoints())
      return;
    typename PointData<T>::ScalarDataType kept(levelSet->getNumberOfPoints(),
                                               T(0));
    for (viennahrle::SparseIterator<typename Domain<T, D>::DomainType> it(
             levelSet->getDomain());
         !it.isFinished(); it.next()) {
      if (!it.isDefined())
        continue;
      const hrleIndex &idx = it.getStartIndices();
      auto found = std::lower_bound(
          before.begin(), before.end(), idx,
          [](const auto &a, const hrleIndex &b) { return a.first < b; });
      if (found == before.end() || found->first < idx || idx < found->first)
        continue;
      const long j = found->second;
      if (!rounded[j] || norm(normalBefore[j]) < T(0.5) ||
          std::abs(valueBefore[j] - it.getValue()) > T(0.25))
        continue;
      const auto id = it.getPointId();
      Vec x{};
      for (int i = 0; i < D; ++i)
        x[i] = T(idx[i]);
      it.getDefinedValue() = valueBefore[j];
      const Vec foot = difference(x, scaled(normalBefore[j], valueBefore[j]));
      for (int i = 0; i < 3; ++i) {
        (*closest)[id][i] = foot[i] * gridDelta;
        (*normals)[id][i] = normalBefore[j][i];
      }
      kept[id] = T(1);
    }
    pointData.insertReplaceScalarData(std::move(kept), smoothLabel);
  }

public:
  // point data of the band state between apply() calls
  static constexpr char anchorLabel[] = "HermiteAnchor";
  static constexpr char facePointLabel[] = "HermiteFacePoint";
  static constexpr char faceNormalLabel[] = "HermiteFaceNormal";
  static constexpr char rebuildNormalLabel[] = "HermiteRebuildNormal";
  static constexpr char featurePointLabel[] = "HermiteFeaturePoint";
  static constexpr char witnessLabel[] = "HermiteWitness";
  static constexpr char rebuildWitnessLabel[] = "HermiteRebuildWitness";
  static constexpr char directionLabel[] = "HermiteDirection";
  static constexpr char valueLabel[] = "HermiteValue";
  static constexpr char sideLabel[] = "HermiteSide";
  static constexpr char flagsLabel[] = "HermiteFlags";
  // nodes on a convex rounded part (roundedNodes), whose cells dual
  // contouring places halfway between the planes' point and the mass point
  static constexpr char smoothLabel[] = "HermiteSmooth";
  static constexpr char coneALabel[] = "HermiteConeA";
  static constexpr char coneBLabel[] = "HermiteConeB";
  static constexpr char timeLabel[] = "HermiteTime";

  AdvectHermite() = default;

  AdvectHermite(SmartPointer<Domain<T, D>> passedLevelSet,
                SmartPointer<VelocityField<T>> passedVelocities)
      : levelSets{passedLevelSet}, velocities(passedVelocities) {}

  /// Material stack with the advected surface last, as in Advect.
  AdvectHermite(std::vector<SmartPointer<Domain<T, D>>> passedLevelSets,
                SmartPointer<VelocityField<T>> passedVelocities)
      : levelSets(passedLevelSets), velocities(passedVelocities) {}

  /// The level set to advect, as the only one of the stack.
  void setLevelSet(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSets = {passedLevelSet};
  }

  /// Puts a level set on top of the material stack; the last one is the
  /// Euclidean surface that is advected, the ones below give the material at
  /// each surface point and are intersected with it afterwards.
  void insertNextLevelSet(SmartPointer<Domain<T, D>> passedLevelSet) {
    levelSets.push_back(passedLevelSet);
  }

  void clearLevelSets() { levelSets.clear(); }

  void setVelocityField(SmartPointer<VelocityField<T>> passedVelocities) {
    velocities = passedVelocities;
  }

  /// Time to advect for.
  void setAdvectionTime(T time) { advectionTime = time; }

  /// Largest move of the front per time step in grid spacings, default 0.5.
  void setTimeStepRatio(T ratio) { timeStepRatio = ratio; }

  /// Front motion in grid spacings after which the band is rebuilt from the
  /// dual contouring mesh, default 2. It must stay well inside the band.
  void setRebuildDistance(T distance) { rebuildDistance = distance; }

  /// Level set width of the rebuilt band, default 12 (6 grid spacings on
  /// either side).
  void setBandWidth(int width) { bandWidth = width; }

  /// Whether apply() stops after one time step, as Advect's single step for
  /// velocities that are updated between steps. Default false.
  void setSingleStep(bool passedSingleStep) { singleStep = passedSingleStep; }

  /// Whether apply() rebuilds the band from the dual contouring mesh of its
  /// planes once the whole advection time is done if the band is no valid
  /// level set (Check), default true. Between rebuilds, values further from
  /// the surface can lag behind edges the motion created, as at the tip of
  /// an undercut; the rebuild leaves the surface as it is. Without it, a
  /// further apply() continues from the band state.
  void setRebuildAtEnd(bool passedRebuildAtEnd) {
    rebuildAtEnd = passedRebuildAtEnd;
  }

  T getAdvectedTime() const { return advectedTime; }

  unsigned getNumberOfTimeSteps() const { return numberOfTimeSteps; }

  unsigned getNumberOfRebuilds() const { return numberOfRebuilds; }

  void apply() {
    advectedTime = 0.;
    numberOfTimeSteps = 0;
    numberOfRebuilds = 0;
    if (levelSets.empty() || levelSets.back() == nullptr ||
        velocities == nullptr) {
      VIENNACORE_LOG_ERROR("AdvectHermite: no level set or velocity field.");
      return;
    }
    levelSet = levelSets.back();
    stack = MaterialStack<T, D>(std::vector<SmartPointer<Domain<T, D>>>(
        levelSets.begin(), levelSets.end() - 1));
    if (levelSet->getNormalization() != NormalizationEnum::EUCLIDEAN) {
      VIENNACORE_LOG_ERROR("AdvectHermite: the level set is not Euclidean.");
      return;
    }
    bool rebuiltLast = false;
    collectNodes();
    if (!restoreState()) {
      // a band narrower than the rebuilt one, as after Reduce or Prune, is
      // widened from its planes first
      if (levelSet->getLevelSetWidth() < bandWidth) {
        Reinitialize<T, D>(levelSet, bandWidth).apply();
        collectNodes();
      }
      if (!initializeFromClosestPoints())
        return;
    }

    while (advectedTime < advectionTime) {
      // rebuild once the front has moved far enough, also where a previous
      // apply() left a band that is due
      if (motionSinceRebuild >= rebuildDistance - tinyLength) {
        rebuild();
        collectNodes();
        if (!initializeFromClosestPoints())
          return;
      }
      const long n = long(S.size());
      std::vector<T> Vd(n, T(0)), Vn(n, T(0));
      T maxSpeed = 0.;
#pragma omp parallel for reduction(max : maxSpeed)
      for (long i = 0; i < n; ++i) {
        Vd[i] = speedOf(HA[i], Dir[i], i, nodeMaterial(i, HA[i], Dir[i]));
        Vn[i] = good[i] ? Vd[i] : speed(FA[i], NRM[i], i);
        maxSpeed = std::max(maxSpeed, std::abs(good[i] ? Vd[i] : Vn[i]));
      }
      const T remaining = advectionTime - advectedTime;
      const T dt = (maxSpeed > 0)
                       ? std::min(timeStepRatio / maxSpeed, remaining)
                       : remaining;
      if (levelSets.size() > 1) {
#pragma omp parallel for
        for (long i = 0; i < n; ++i) {
          T &v = good[i] ? Vd[i] : Vn[i];
          if (v < 0)
            v = splitAtInterface(good[i] ? HA[i] : FA[i],
                                 good[i] ? Dir[i] : NRM[i], v, dt, i);
        }
      }
      if (maxSpeed > 0)
        step(dt, Vd, Vn);
      if (levelSets.size() > 1)
        exposeInterfaces();
      advectedTime += dt;
      ++numberOfTimeSteps;
      motionSinceRebuild += maxSpeed * dt;
      // rebuild right after the step that reaches the distance, so that
      // between calls the band and its node ids are those of the next step,
      // as quantities computed on the level set in between (visibilities)
      // expect; the last step of a call leaves the rebuilt band as it is, and
      // the next call classifies its nodes with the velocities given for them
      const bool last =
          singleStep ||
          advectedTime >=
              advectionTime - tinyLength * std::max(T(1), advectionTime);
      if (motionSinceRebuild >= rebuildDistance - tinyLength) {
        rebuild();
        if (last) {
          rebuiltLast = true;
          break;
        }
        collectNodes();
        if (!initializeFromClosestPoints())
          return;
      }
      if (singleStep)
        break;
    }
    const bool done =
        advectedTime >=
        advectionTime - tinyLength * std::max(T(1), advectionTime);
    if (!rebuiltLast) {
      writeBack();
      // values further from the surface can lag behind edges the motion
      // created since the last rebuild, as at the tip of an undercut; a band
      // that is no valid level set then is rebuilt
      if (done && rebuildAtEnd) {
        Check<T, D> check(levelSet);
        check.apply();
        if (!check.isValid() && boundFarValues()) {
          writeBack();
          check.apply();
        }
        if (!check.isValid()) {
          Reinitialize<T, D>(levelSet, bandWidth).apply();
          ++numberOfRebuilds;
        }
      }
    }
    // the lower level sets end at the surface, as in Advect, once the whole
    // advection time is done: between single steps they keep the interfaces
    // and materials of the stack, which deviations of the surface from a
    // static material would otherwise erode step by step
    if (done)
      for (std::size_t l = 0; l + 1 < levelSets.size(); ++l)
        BooleanOperation<T, D>(levelSets[l], levelSet,
                               BooleanOperationEnum::INTERSECT)
            .apply();
  }
};

// add all template specialisations for this class
PRECOMPILE_PRECISION_DIMENSION(AdvectHermite)

} // namespace viennals
