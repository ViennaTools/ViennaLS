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
/// by t V(D), with t the time over which V gives the fan its own offset. Pieces
/// of the front moving at different speeds combine through the envelope of
/// neighbouring planes: the largest signed distance at convex junctions, the
/// smallest at concave ones, with up to three pieces combined exactly at saddle
/// corners. Once the front has moved a given distance, the band moves with it:
/// the dual contouring mesh of the planes gives the nodes of the new band, the
/// nodes of the old one keep their states, and a node new to the band takes the
/// plane of an old node near it, or the fan of a rounded part whose cone holds
/// its direction, or beside another fan a fan of the same feature with the
/// mesh's plane in its own direction. Sharp edges and corners stay sharp where
/// the motion keeps
/// them sharp, and fans round them where it rounds them. A node that takes the
/// plane of a fan becomes a fan of the same feature, and pieces of one fan meet
/// with its convexity.
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
/// edge. A concave fan whose faces a velocity varying in space has moved
/// away from its feature point, as the corner at the foot of a trench that
/// slides down the wall, hands the nodes beside its faces to those faces,
/// which move with the speeds at their own points.
///
/// The level set must be Euclidean with closest points and normals, as created
/// by FromSurfaceMesh or MakeGeometry. Level sets below it in a material stack
/// (as in Advect) give the material at each surface point. Where etching
/// reaches a lower material, the step is split at the interface, and where it
/// lays an interface open from below, as under a mask, the interface plane
/// comes from the lower level set. A neighbouring face acts on a node where
/// the node's foot on it lies on the face's material, so that the wall of a
/// mask does not reach into the undercut below it, and a fan acts on the gas
/// and on its own material, so that a pit etched into the substrate does not
/// cut into the mask, in either order of the stack. Every plane carries its
/// material, a face's as a fan's, which a node taking it takes over and which
/// changes only where the face's own motion crosses an interface; a fan's node
/// in the solid takes the material at its own position. Behind a junction of
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
/// values far from the surface lag behind edges the motion created, is
/// repaired from the mesh of its planes (setRebuildAtEnd). A band narrower
/// than the rebuilt one, as left by Reduce or Prune, is widened first.
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
  // cone of every fan whose feature is a 3D corner: the normals of the three
  // faces meeting there, in arcA, arcB and arcC
  std::vector<Vec> arcC;
  std::vector<char> cornerKnown;
  // solid nodes that took a moving plane of the etched material around them
  // from exposeInterfaces, as the plane of the end of a cavity below a mask
  // past which they lie; dual contouring prefers it (priority 3) for the
  // crossing edges of their cells it cuts
  std::vector<char> boundPlane;
  // nodes that a previous apply() left on a rounded part (smoothLabel), as on
  // the arc a deposition rounds an edge into
  std::vector<char> smooth;
  // the rounded parts of the band a rebuild made, from those of the old band:
  // the data a velocity field keeps per point (as visibilities) are those of
  // the old band until the next step, and roundedNodes reads no speeds of
  // the new one before it
  std::vector<char> rebuiltRounded;
  bool rebuiltBand = false;
  // speeds known only for the surface's own normals (concaveFansToFaces, the
  // fitted turn)
  bool surfaceSpeedsOnly = false;
  // time over which each node's planes have moved since they were set, at
  // the rebuild that found the node or with the fan it took over: the time
  // of a fan's offsets at its feature point
  std::vector<T> since;
  // material of every node's piece, set with the piece and carried with it
  // when another node takes it over; it changes where the piece's own motion
  // crosses an interface of the material stack
  std::vector<int> pieceMaterial;
  // side of every node when it entered the band, which tells the nodes
  // etched since
  std::vector<T> entrySide;
  // nodes whose plane exposeInterfaces set: the interface or the plane of
  // the end of a cavity
  std::vector<char> exposedState;
  // the interface plane each such node found, its normal towards the gas,
  // which the solid nodes further from the cavity take on, and the lower
  // level set it belongs to
  std::vector<Vec> interfacePoint, interfaceNormal;
  std::vector<int> interfaceLevel;
  T sinceRebuild = 0.; // time since the last rebuild
  T motionSinceRebuild =
      0.; // front motion since the last rebuild, grid spacings

  // the two faces featureCones found through a fan's feature point, also
  // where its cone stays unknown, for concaveFansToFaces
  std::vector<char> facePair;
  std::vector<Vec> pairA, pairB;

  // whether the unit vector u lies on the shorter arc from a to b
  static bool within(const Vec &u, const Vec &a, const Vec &b) {
    const Vec z = CrossProduct(a, b);
    return DotProduct(CrossProduct(a, u), z) >= -unitTolerance &&
           DotProduct(CrossProduct(u, b), z) >= -unitTolerance &&
           DotProduct(u, a + b) > 0;
  }
  static T sign(T v) { return T((v > 0) - (v < 0)); }
  // whether the unit vector u lies in the cone of the normals a, b and c
  static bool inCornerCone(const Vec &u, const Vec &a, const Vec &b,
                           const Vec &c) {
    const T det = DotProduct(a, CrossProduct(b, c));
    if (std::abs(det) < T(1e-6))
      return false;
    return DotProduct(u, CrossProduct(b, c)) / det >= -unitTolerance &&
           DotProduct(a, CrossProduct(u, c)) / det >= -unitTolerance &&
           DotProduct(a, CrossProduct(b, u)) / det >= -unitTolerance;
  }

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
  // feature is a 2D corner or a 3D edge, with three in 3D a corner, otherwise
  // unknown
  void featureCones() {
    const long n = long(S.size());
    arcA.assign(n, Vec{});
    arcB.assign(n, Vec{});
    arcKnown.assign(n, 0);
    arcC.assign(n, Vec{});
    cornerKnown.assign(n, 0);
    facePair.assign(n, 0);
    pairA.assign(n, Vec{});
    pairB.assign(n, Vec{});
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
        const Vec r = X[m] - C0[m];
        const T length = Norm(r);
        if (length > smallLength &&
            std::abs(DotProduct(r, NRM0[m])) < T(0.999) * length)
          continue;
        // to the accuracy of the rebuilt mesh, whose faces beside an edge
        // pass a few hundredths of a grid spacing beside the feature point
        // and tilt by a degree or two
        const Vec toFeature = c - (C0[m] + shift);
        if (std::abs(DotProduct(toFeature, NRM0[m])) >
            rebuiltLength + rebuiltSlope * Norm(toFeature))
          continue;
        bool seen = false;
        for (int v = 0; v < found; ++v)
          seen = seen || DotProduct(normals[v], NRM0[m]) >=
                             1 - rebuiltTilt * rebuiltTilt / 2;
        if (!seen)
          normals[found++] = NRM0[m];
      }
      // a 3D corner of three faces
      if (D == 3 && found == 3) {
        arcA[i] = normals[0];
        arcB[i] = normals[1];
        arcC[i] = normals[2];
        cornerKnown[i] = 1;
        continue;
      }
      if (found != 2 ||
          DotProduct(normals[0], normals[1]) < T(-1) + unitTolerance)
        continue;
      facePair[i] = 1;
      pairA[i] = normals[0];
      pairB[i] = normals[1];
      const Vec z = CrossProduct(normals[0], normals[1]);
      if (D == 3 && std::abs(DotProduct(Dir[i], z)) > rebuiltTilt * Norm(z))
        continue;
      // or past one of them, into the part of the cone of its material
      // beyond a face of another material
      if (!within(Dir[i], normals[0], normals[1])) {
        if (levelSets.size() < 2 || convexFan[i])
          continue;
        Vec a = normals[0], b = normals[1], past, end;
        if (materialCone(C0[i], nodeMaterial(i), a, b, past, end) != 2 ||
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
  // visibilities) get the node nearest to the point a grid spacing in front
  // of p along n, which sees what the surface point sees: the node nearest to
  // p itself can lie behind the surface or beside it, as below the edge of a
  // face that sees the source, where it is shadowed
  T speedOf(const Vec &p, const Vec &n, unsigned long id, int material) const {
    const Vec coordinate = p * gridDelta;
    const long near = nodeAt(p + n);
    const unsigned long point = (near >= 0) ? (unsigned long)near : id;
    const T scalar =
        velocities->getScalarVelocity(coordinate, material, n, point);
    const Vec vector =
        velocities->getVectorVelocity(coordinate, material, n, point);
    return (scalar + DotProduct(vector, n)) / gridDelta;
  }

  // whether node m, a face with normal n now, still carries the face it held
  // at the rebuild: its plane has that face's normal and its piece's point at
  // the rebuild lies on the face, so that a node that took the plane of a
  // parallel face, as the floor below the top of a mask, does not stand for
  // that face
  bool carriesFace(long m, const Vec &n) const {
    return DotProduct(n, NRM0[m]) > 1 - unitTolerance &&
           std::abs(DotProduct(RW[m] - C0[m], NRM0[m])) <= rebuiltLength;
  }

  // the material of a surface point is that just behind it, against its
  // normal, so that both sides of an interface are told apart
  int materialBehind(const Vec &p, const Vec &n) const {
    return material(p - n * T(0.01));
  }

  // the material of node id's piece: the one its plane carries, set where
  // the plane was created and passed on with it, also a fan's, whose surface
  // point can lie on an interface, as the end of the cone of the foot of a
  // mask wall along the mask's underside, where the material just behind it
  // is that above the interface. A solid node whose surface point is an edge
  // or corner (a fan on the solid side) lies behind a reflex junction, as
  // below the foot of a mask wall beside the opening, and the front that
  // reaches it is that of its own material: just behind the junction may lie
  // a sliver of another, as where a mask reaches a little below the
  // substrate surface
  int nodeMaterial(long id) const {
    if (levelSets.size() > 1 && fan[id] && S[id] < 0)
      return material(X[id]);
    return pieceMaterial[id];
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
      return material(c - d * probe) == m || material(c - d) == m;
    };
    // the other material, just beyond the cut
    auto behind = [&](const Vec &d) {
      const int near = material(c - d * probe);
      return (near != m) ? near : material(c - d);
    };
    const bool inA = inside(a), inB = inside(b);
    if (inA != inB) {
      auto direction = [&](T f) {
        const Vec d = a * (1 - f) + b * f;
        return d * (1 / Norm(d));
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
          DotProduct(cut, inA ? b : a) > coneTolerance || !within(cut, a, b))
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
      const Vec t = o - e * DotProduct(o, e);
      return t * (1 / Norm(t));
    };
    const bool otherA = materialBehind(c + along(a, b) * probe, a) != m,
               otherB = materialBehind(c + along(b, a) * probe, b) != m;
    if (otherA == otherB)
      return 0;
    Vec &e = otherB ? b : a;
    const Vec t = otherB ? along(b, a) : along(a, b);
    // directions turning from e away from the other face
    const Vec start = e;
    auto turned = [&](T phi) {
      return start * std::cos(phi) + t * -std::sin(phi);
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
        DotProduct(cut, start) > coneTolerance || DotProduct(cut, t) > 0)
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
    const Vec q = c - cut;
    const std::size_t l = std::size_t(std::min(m, other));
    if (other == m || l + 1 >= levelSets.size())
      return false;
    Vec gradient{};
    constexpr T h = T(0.05);
    for (int k = 0; k < D; ++k) {
      Vec step{};
      step[k] = h;
      gradient[k] =
          (lowerValue(l, q + step) - lowerValue(l, q - step)) / (2 * h);
    }
    // in the plane of the cone
    if constexpr (D == 3) {
      Vec z = CrossProduct(a, b);
      z = z * (1 / Norm(z));
      gradient = gradient - z * DotProduct(gradient, z);
    }
    const T length = Norm(gradient);
    if (length < tinyLength)
      return false;
    gradient = gradient * (1 / length);
    // the direction normal to it nearest the probed one
    Vec along = cut - gradient * DotProduct(cut, gradient);
    const T lengthAlong = Norm(along);
    if (lengthAlong < T(0.5))
      return false;
    cut = along * (1 / lengthAlong);
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
    const int m0 = nodeMaterial(long(id));
    const Vec p1 = p + n * (dt * v);
    const int m1 = materialBehind(p1, n);
    if (m1 == m0)
      return v;
    const std::size_t l = std::size_t(std::min(m0, m1));
    const T a = lowerValue(l, p), b = lowerValue(l, p1);
    const T f = (a * b < 0) ? std::clamp(a / (a - b), T(0), T(1)) : T(1);
    const Vec crossing = p + n * (f * dt * v);
    return f * v + (1 - f) * speedOf(crossing, n, id, m1);
  }

  // +1 where two pieces meet convexly, -1 concavely: the normal turns with
  // the path between the witnesses on convex bends; pieces meeting on one
  // edge or corner (fans) count as convex outside the material and concave
  // inside
  static int classify(const Vec &na, const Vec &wa, const Vec &nb,
                      const Vec &wb, T side) {
    const Vec dn = nb - na, dp = wb - wa;
    T b = DotProduct(dn, dp) / std::max(Norm(dn) * Norm(dp), T(1e-12));
    // witnesses apart along the line where the two planes meet leave the bend
    // undecided over the whole path: its part across that line decides
    if constexpr (D == 3)
      if (std::abs(b) <= bendLimit) {
        const Vec e = CrossProduct(na, nb);
        const T e2 = DotProduct(e, e);
        if (e2 > T(1e-4)) {
          const Vec across = dp - e * (DotProduct(dp, e) / e2);
          if (Norm(across) > smallLength)
            b = DotProduct(dn, across) / (Norm(dn) * Norm(across));
        }
      }
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
    since.assign(n, T(0));
    pieceMaterial.assign(n, 0);
    exposedState.assign(n, 0);
    interfacePoint.assign(n, Vec{});
    interfaceNormal.assign(n, Vec{});
    interfaceLevel.assign(n, -1);
    entrySide.assign(n, T(0));
    for (std::size_t id = 0; id < n; ++id)
      entrySide[id] = sign(S[id]);
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
      const Vec r = X[id] - cp;
      const T length = Norm(r);
      const T cosine = DotProduct(r, nrm) / std::max(length, T(1e-12));
      // a folded triangle of the rebuilt mesh can give a node whose closest
      // point lies inside it a normal against its own side
      if (std::abs(cosine) > T(0.999) && sign(cosine) != sign(S[id]))
        nrm = nrm * T(-1);
      const Vec d = (length > smallLength) ? r * (sign(S[id]) / length) : nrm;
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
      pieceMaterial[id] = m;
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
    concaveFansToFaces();
    return true;
  }

  // two directions along a face with the normal n0, orthonormal with it
  std::array<Vec, D - 1> faceTangents(const Vec &n0) const {
    std::array<Vec, D - 1> tangents;
    if constexpr (D == 2) {
      tangents[0] = Vec{-n0[1], n0[0], T(0)};
    } else {
      int k = 0;
      for (int c = 1; c < 3; ++c)
        if (std::abs(n0[c]) < std::abs(n0[k]))
          k = c;
      Vec e{};
      e[k] = 1;
      Vec t1 = CrossProduct(n0, e);
      t1 = t1 * (1 / Norm(t1));
      tangents[0] = t1;
      tangents[1] = CrossProduct(n0, t1);
    }
    return tangents;
  }

  // under speeds known only for the surface's own normals, as fluxes from
  // ray tracing, a concave fan that moves would need the speeds of the
  // directions between its faces at its feature point, which the velocity
  // field does not define: such a corner is the meeting point of its faces
  // instead, of two at an edge and of three at a 3D corner, as at the foot
  // of a curved mask wall. A fan whose cone stays unknown with two faces
  // found through its feature point is such an edge as well. At the foot of
  // a curved mask wall, the face of the wall found beside the feature point
  // makes the edge look tilted against the fan's direction. Kept as a fan,
  // it would act on the floor nodes beyond its cone and hold them to its
  // plane, which moves with the speed at the mask's foot. Each node of the
  // fan takes the plane of the face that gives it the smallest value,
  // through the feature point, and moves with the speed at its own surface
  // point; its witness lies on the face with the next smallest value a grid
  // spacing from the corner, so that the faces meet concavely
  void concaveFansToFaces() {
    if (!surfaceSpeedsOnly)
      return;
    const long n = long(S.size());
    std::vector<char> take(n, 0);
    std::vector<Vec> normal(n), witness(n);
    std::vector<T> value(n);
#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      const bool corner =
          !arcKnown[i] && cornerKnown.size() == S.size() && cornerKnown[i];
      const bool pair =
          !arcKnown[i] && !corner && facePair.size() == S.size() && facePair[i];
      if (!fan[i] || convexFan[i] || !good[i] ||
          !(arcKnown[i] || corner || pair))
        continue;
      const Vec &faceA = pair ? pairA[i] : arcA[i];
      const Vec &faceB = pair ? pairB[i] : arcB[i];
      const int m = nodeMaterial(i);
      if (speedOf(C0[i], faceA, i, m) == 0 &&
          speedOf(C0[i], faceB, i, m) == 0 &&
          (!corner || speedOf(C0[i], arcC[i], i, m) == 0))
        continue;
      // the faces through the feature point, the nearest first
      const Vec *faces[3] = {&faceA, &faceB, &arcC[i]};
      T values[3] = {DotProduct(X[i] - C0[i], faceA),
                     DotProduct(X[i] - C0[i], faceB),
                     corner ? DotProduct(X[i] - C0[i], arcC[i])
                            : std::numeric_limits<T>::max()};
      int order[3] = {0, 1, 2};
      std::sort(order, order + 3,
                [&](int p, int q) { return values[p] < values[q]; });
      const T a = values[order[0]], b = values[order[1]];
      const Vec &own = *faces[order[0]];
      const Vec &other = *faces[order[1]];
      Vec along = other - own * DotProduct(own, other);
      const T length = Norm(along);
      if (length < smallLength)
        continue;
      along = along * (1 / length);
      take[i] = 1;
      normal[i] = own;
      value[i] = std::min(a, b);
      witness[i] = C0[i] + along;
    }
    for (long i = 0; i < n; ++i) {
      if (!take[i])
        continue;
      S[i] = value[i];
      const Vec foot = X[i] - normal[i] * S[i];
      HA[i] = FA[i] = PL[i] = C0[i] = RW[i] = foot;
      W[i] = witness[i];
      NRM[i] = NRM0[i] = Dir[i] = NL[i] = normal[i];
      fan[i] = 0;
      arcKnown[i] = 0;
      cornerKnown[i] = 0;
      since[i] = 0;
      moved[i] = 1;
    }
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
    const auto *sinces = pointData.getScalarData(sinceLabel, true);
    const auto *materials = pointData.getScalarData(materialLabel, true);
    const auto *entries = pointData.getScalarData(entryLabel, true);
    if (values == nullptr || sides == nullptr || flags == nullptr ||
        times == nullptr || sinces == nullptr || values->size() != n ||
        sides->size() != n || flags->size() != n || sinces->size() != n ||
        times->size() < 2 || materials == nullptr || materials->size() != n ||
        entries == nullptr || entries->size() != n)
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
    exposedState.assign(n, 0);
    pieceMaterial.assign(n, 0);
    interfacePoint.assign(n, Vec{});
    interfaceNormal.assign(n, Vec{});
    interfaceLevel.assign(n, -1);
    entrySide.assign(entries->begin(), entries->end());
    for (std::size_t i = 0; i < n; ++i) {
      side0[i] = (*sides)[i];
      pieceMaterial[i] = int(std::lround((*materials)[i]));
      const int bits = int((*flags)[i]);
      good[i] = bits & 1;
      fan[i] = (bits >> 1) & 1;
      moved[i] = (bits >> 2) & 1;
      convexFan[i] = (bits >> 3) & 1;
      arcKnown[i] = (bits >> 4) & 1;
      boundPlane[i] = (bits >> 5) & 1;
      smooth[i] = (bits >> 6) & 1;
      exposedState[i] = (bits >> 7) & 1;
    }
    arcC.assign(n, Vec{});
    cornerKnown.assign(n, 0);
    if (const auto *cones = pointData.getVectorData(coneCLabel, true);
        cones != nullptr && cones->size() == n)
      for (std::size_t i = 0; i < n; ++i) {
        for (int c = 0; c < 3; ++c)
          arcC[i][c] = (*cones)[i][c];
        cornerKnown[i] = (int((*flags)[i]) >> 8) & 1;
      }
    since.assign(sinces->begin(), sinces->end());
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
  // its support point on the moved fan. Only the directions of the fan's cone
  // count: between ma and mb, the faces of an edge or the part of its own
  // material where an interface splits it, and in 3D across the edge (cone
  // 1), or in the cone of the three faces ma, mb and mc of a 3D corner (cone
  // 2). A direction outside it, as the normal of another face that a
  // neighbour carries, belongs to another feature, and a plane of the fan in
  // that direction would cut into the material
  void arc(long i, const Vec dx, T ox, const Vec &dy, T oy, const Vec &c, int m,
           T t, bool convex, int cone, const Vec &ma, const Vec &mb,
           const Vec &mc, bool refine, T &s, Vec &p, Vec &nn) const {
    static constexpr T samples[] = {0.125, 0.25, 0.375, 0.5,
                                    0.625, 0.75, 0.875, 1.0};
    const T cx = ox - t * speedOf(c, dx, i, m);
    const T cy = oy - t * speedOf(c, dy, i, m);
    // the samples, and the ends of the cone between dx and dy
    std::array<T, 10> fractions{};
    int count = 0;
    for (T f : samples)
      fractions[count++] = f;
    if (cone == 1)
      for (const Vec &e : {ma, mb}) {
        const T ex = Norm(CrossProduct(dx, e)), ey = Norm(CrossProduct(dy, e));
        if (ex + ey > tinyLength && within(e, dx, dy))
          fractions[count++] = ex / (ex + ey);
      }
    // the signed distance to the support plane at fraction f, false where
    // the direction is not one of the fan's
    auto planeAt = [&](T f, Vec &d, T &offset, T &value) {
      d = dx * (1 - f) + dy * f;
      const T length = Norm(d);
      if (length < tinyLength)
        return false;
      d = d * (1 / length);
      if (cone == 1) {
        if (!within(d, ma, mb))
          return false;
        if constexpr (D == 3) {
          const Vec z = CrossProduct(ma, mb);
          if (std::abs(DotProduct(d, z)) > rebuiltTilt * Norm(z))
            return false;
        }
      }
      if (cone == 2 && !inCornerCone(d, ma, mb, mc))
        return false;
      offset = t * speedOf(c, d, i, m) + (1 - f) * cx + f * cy;
      value = DotProduct(X[i] - c, d) - offset;
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
        p = c + d * offset;
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
      p = c + nn * (first ? o1 : o2);
    }
  }

  void step(T dt, const std::vector<T> &Vd, const std::vector<T> &Vn) {
    const long n = long(S.size());
    rebuiltBand = false;
    std::vector<Vec> P(n), N(n), Wm(n);
    std::vector<T> Vp(n), Sown(n);
    std::vector<char> moving(n);
#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      moving[i] = good[i] && Vd[i] != 0;
      P[i] = good[i] ? HA[i] + Dir[i] * (dt * Vd[i])
                     : FA[i] + NRM[i] * (dt * Vn[i]);
      N[i] = (good[i] && (moved[i] || Vd[i] != 0)) ? Dir[i] : NRM[i];
      Vp[i] = good[i] ? Vd[i] : Vn[i];
      // moving nodes take the value of their own characteristic, since a
      // fan node may hold the value of a sampled support plane
      Sown[i] = moving[i] ? DotProduct(X[i] - P[i], N[i]) : S[i] - dt * Vp[i];
      // witnesses move with their piece; a fan's stays at its feature point
      Wm[i] = fan[i] ? W[i] : W[i] + N[i] * (dt * Vp[i]);
    }
    // faces turn with the rate along them: along a characteristic of the
    // level set equation the normal turns as dn/dt = -(I - n n^T) grad V,
    // with the gradient of the rate along the face at the face's own normal,
    // from central differences a grid spacing to either side of the face's
    // point; it vanishes for a rate of the normal and the material alone.
    // The point moves along the normal halfway through the turn, with the
    // rate halfway along its path, and the plane takes the turned normal
    std::vector<char> turned(n, 0);
    const std::vector<Vec> startN = N;
#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      if (fan[i])
        continue;
      const Vec a = good[i] ? HA[i] : FA[i];
      const Vec n0 = N[i];
      const int m = nodeMaterial(i);
      const std::array<Vec, D - 1> tangents = faceTangents(n0);
      Vec g{};
      for (const Vec &t : tangents) {
        const T vp = speedOf(a + t, n0, i, m), vm = speedOf(a - t, n0, i, m);
        g = g + t * ((vp - vm) / T(2));
      }
      // a rate of the normal and the material alone keeps the face as it is
      if (Norm(g) == 0)
        continue;
      const Vec delta = g * (-dt);
      Vec half = n0 + delta * T(0.5), end = n0 + delta;
      half = half * (1 / Norm(half));
      end = end * (1 / Norm(end));
      const T v0 = speedOf(a, n0, i, m),
              vHalf = speedOf(a + n0 * (dt * Vp[i] / 2), half, i, m);
      const T v = Vp[i] + (vHalf - v0);
      P[i] = a + half * (dt * v);
      N[i] = end;
      Vp[i] = v;
      Sown[i] = DotProduct(X[i] - P[i], end);
      Wm[i] = W[i] + half * (dt * v);
      turned[i] = 1;
    }
    // speeds known only at the surface's own normals put into the gradient
    // the change of the rate with the normal times the change of the normal.
    // The turned normal is then that of a least-squares plane through the
    // moved points of the face nodes two grid spacings to either side along
    // the face (normals within 20 degrees, the same material, points within
    // 1.5 grid spacings of the face and 0.75 of their place along it), with a
    // slope only along the directions in which these points spread (variance
    // of their positions of at least a quarter grid spacing squared) and the
    // gradient's slope along the others; with fewer than three points the
    // face keeps the gradient's turn. On a line of points, as around a curved
    // corner where the normals of the others differ by more than 20 degrees,
    // the plane would take any slope across the line
    if (surfaceSpeedsOnly) {
      constexpr int R = 2;
      const T fitSpread = T(0.25);
      std::vector<Vec> fitted(N);
#pragma omp parallel for
      for (long i = 0; i < n; ++i) {
        if (!turned[i])
          continue;
        const Vec n0 = startN[i];
        const int m = nodeMaterial(i);
        const std::array<Vec, D - 1> tangents = faceTangents(n0);
        // normal equations of h = c0 + sum c_k s_k over the points
        constexpr int U = D; // unknowns: c0 and D - 1 slopes
        T A[U][U] = {}, b[U] = {};
        int count = 0;
        std::vector<long> used;
        auto add = [&](const Vec &q) {
          const Vec r = q - P[i];
          T row[U];
          row[0] = 1;
          for (int k = 0; k < D - 1; ++k)
            row[k + 1] = DotProduct(r, tangents[k]);
          const T h = DotProduct(r, n0);
          for (int a = 0; a < U; ++a) {
            b[a] += row[a] * h;
            for (int c = 0; c < U; ++c)
              A[a][c] += row[a] * row[c];
          }
          ++count;
        };
        add(P[i]);
        used.push_back(i);
        for (int k1 = -R; k1 <= R; ++k1)
          for (int k2 = (D == 3 ? -R : 0); k2 <= (D == 3 ? R : 0); ++k2) {
            if (k1 == 0 && k2 == 0)
              continue;
            Vec q = P[i] + tangents[0] * T(k1);
            if constexpr (D == 3)
              q = q + tangents[1] * T(k2);
            const long j = nodeAt(q);
            if (j < 0 || fan[j] || nodeMaterial(j) != m)
              continue;
            const Vec pj = P[j], nj = startN[j];
            if (std::find(used.begin(), used.end(), j) != used.end())
              continue;
            if (DotProduct(nj, n0) < T(0.93969262078590838)) // 20 deg
              continue;
            const Vec r = pj - P[i];
            if (std::abs(DotProduct(r, n0)) > T(1.5) ||
                std::abs(DotProduct(r, tangents[0]) - T(k1)) > T(0.75))
              continue;
            if constexpr (D == 3)
              if (std::abs(DotProduct(r, tangents[1]) - T(k2)) > T(0.75))
                continue;
            used.push_back(j);
            add(pj);
          }
        if (count < 3)
          continue;
        // the slopes of the gradient's turn (fitted holds it here)
        const Vec ncd = fitted[i];
        const T along = DotProduct(ncd, n0);
        const T cd1 =
            along > T(0) ? -DotProduct(ncd, tangents[0]) / along : T(0);
        const T cd2 = (D == 3 && along > T(0))
                          ? -DotProduct(ncd, tangents[D - 2]) / along
                          : T(0);
        // the slopes along the principal directions of the points' spread
        Vec nn = n0;
        const T cnt = A[0][0];
        if constexpr (D == 3) {
          const T c11 = A[1][1] - A[0][1] * A[0][1] / cnt;
          const T c22 = A[2][2] - A[0][2] * A[0][2] / cnt;
          const T c12 = A[1][2] - A[0][1] * A[0][2] / cnt;
          const T g1 = b[1] - A[0][1] * b[0] / cnt;
          const T g2 = b[2] - A[0][2] * b[0] / cnt;
          const T mean = T(0.5) * (c11 + c22);
          const T dev =
              std::sqrt(T(0.25) * (c11 - c22) * (c11 - c22) + c12 * c12);
          const T lam[2] = {mean + dev, mean - dev};
          for (int e = 0; e < 2; ++e) {
            T ex, ey;
            if (std::abs(c12) > T(1e-12) * (std::abs(c11) + std::abs(c22))) {
              ex = lam[e] - c22;
              ey = c12;
            } else {
              const bool first = (c11 >= c22) == (e == 0);
              ex = first ? T(1) : T(0);
              ey = first ? T(0) : T(1);
            }
            const T len = std::sqrt(ex * ex + ey * ey);
            ex /= len;
            ey /= len;
            const T slope = lam[e] >= fitSpread * cnt
                                ? (ex * g1 + ey * g2) / lam[e]
                                : ex * cd1 + ey * cd2;
            nn = nn - (tangents[0] * ex + tangents[D - 2] * ey) * slope;
          }
        } else {
          const T c11 = A[1][1] - A[0][1] * A[0][1] / cnt;
          const T g1 = b[1] - A[0][1] * b[0] / cnt;
          nn = nn - tangents[0] * (c11 >= fitSpread * cnt ? g1 / c11 : cd1);
        }
        fitted[i] = nn * (1 / Norm(nn));
      }
#pragma omp parallel for
      for (long i = 0; i < n; ++i)
        if (turned[i]) {
          N[i] = fitted[i];
          Sown[i] = DotProduct(X[i] - P[i], N[i]);
        }
    }
    sinceRebuild += dt;
    for (auto &t : since)
      t += dt;

    // in a material stack, the material of every piece: a face's is the one
    // it carries; a fan's that of its region, since the ends of its cone may
    // lie on an interface, and on the solid side, behind a reflex junction,
    // that of the node itself
    const bool stacked = levelSets.size() > 1;
    std::vector<int> piece(n, 0), region(n, 0);
    if (stacked) {
#pragma omp parallel for
      for (long i = 0; i < n; ++i) {
        region[i] = material(X[i]);
        piece[i] = nodeMaterial(i);
      }
    }
    // the material of every piece after its own motion in this step: a
    // face whose point crosses an interface is of the material beyond it
    std::vector<int> pieceAfter(piece);
    if (stacked) {
#pragma omp parallel for
      for (long i = 0; i < n; ++i) {
        if (fan[i])
          continue;
        const Vec &start = good[i] ? HA[i] : FA[i];
        const int before = materialBehind(start, N[i]),
                  after = materialBehind(P[i], N[i]);
        if (after != before)
          pieceAfter[i] = after;
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
    std::vector<T> elapsed(since);
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
        ov += own * DotProduct(P[j] - C0[j], N[j]);
        vv += own * own;
      }
      for (int q = 0; q < numberOfNeighbours; ++q) {
        const long m = neighbours[j][q];
        if (m < 0 || fan[m] || !good[m] || piece[m] != piece[j])
          continue;
        const Vec &sh = shiftOf(j, q);
        if ((DotProduct(NRM0[m], arcA[j]) > 1 - unitTolerance ||
             DotProduct(NRM0[m], arcB[j]) > 1 - unitTolerance) &&
            carriesFace(m, N[m]) &&
            std::abs(DotProduct(C0[m] + sh - C0[j], NRM0[m])) <= smallLength) {
          const T v = speedOf(C0[j], NRM0[m], j, piece[j]);
          if (v > slow) {
            ov += v * DotProduct(P[m] + sh - C0[j], NRM0[m]);
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
          const Vec dc = C0[m] + shiftOf(j, q) - C0[j];
          const bool sameArc =
              (DotProduct(arcA[m], arcA[j]) > 1 - unitTolerance &&
               DotProduct(arcB[m], arcB[j]) > 1 - unitTolerance) ||
              (DotProduct(arcA[m], arcB[j]) > 1 - unitTolerance &&
               DotProduct(arcB[m], arcA[j]) > 1 - unitTolerance);
          if (sameArc && std::abs(DotProduct(dc, arcA[j])) <= smallLength &&
              std::abs(DotProduct(dc, arcB[j])) <= smallLength) {
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

    // the time of the offsets of every other moving fan: the time over which
    // the speed at its feature point in its own direction gives the offset
    // its own plane accumulated, so that the offsets of directions no
    // neighbour carries scale with its own. Like the elapsed time of a
    // growing fan it is exact where the velocity changes only its rate, and
    // for a rate alike in every direction it does not depend on the point
    // the rate is read at. Where no time gives the fan its own offset, as
    // where a flux read at a feature point that the front has left far
    // behind is zero, the fan extrapolates no offsets (time -1), and so where
    // the time exceeds ten times the time since its planes were set: a speed
    // below a tenth of the mean that built the offset, as a flux that the
    // walls of a deepening trench shadow at the foot of its mask, would scale
    // the offsets of other directions by an arbitrary time. A concave fan
    // whose own direction moves at less than half the speed of the faster
    // end of its cone takes the time since its planes were set: the error of
    // the slow offset, as of the wall beside the corner at the foot of a
    // trench, would reach the fast directions, as the floor, multiplied by
    // the ratio of the speeds
    std::vector<T> fanTime(elapsed);
#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      if (!fan[i] || !moving[i] || arcGrows[i])
        continue;
      const T v = speedOf(C0[i], N[i], i, piece[i]);
      const T own = DotProduct(P[i] - C0[i], N[i]);
      if (!convexFan[i] && arcKnown[i] && since[i] > 0) {
        const T va = std::abs(speedOf(C0[i], arcA[i], i, piece[i])),
                vb = std::abs(speedOf(C0[i], arcB[i], i, piece[i]));
        if (std::abs(v) < T(0.5) * std::max(va, vb)) {
          fanTime[i] = since[i];
          continue;
        }
      }
      if (v != 0 && own / v > 0 && own / v <= 10 * since[i])
        fanTime[i] = own / v;
      else if (std::abs(own) > smallLength)
        fanTime[i] = -1;
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
        // a fan without a time samples only between offsets that its
        // neighbours carry, along the linear blend of the two
        const bool timed = fanTime[i] >= 0;
        const T t = timed ? fanTime[i] : T(0);
        // the fan's cone: of an edge or 2D corner (the part of its material
        // where an interface splits it), or of a 3D corner
        const int cone = arcKnown[i] ? 1 : (cornerKnown[i] ? 2 : 0);
        const Vec &sh = shiftOf(i, k);
        const Vec dc = C0[j] + sh - C0[i];
        const bool same =
            fan[j] && std::abs(DotProduct(dc, Nf[i])) <= smallLength &&
            std::abs(DotProduct(dc, previous[j])) <= smallLength &&
            DotProduct(Nf[i], previous[j]) < 1 - unitTolerance;
        const bool face = good[j] && !fan[j] &&
                          std::abs(DotProduct(dc, NRM0[j])) <= smallLength &&
                          DotProduct(Nf[i], NRM0[j]) < 1 - unitTolerance;
        // a fan of the same feature of another material moves at another
        // speed, as the static mask behind the foot of its wall
        if (same && (timed || piece[j] == piece[i]))
          arc(i, Nf[i], DotProduct(Pf[i] - C0[i], Nf[i]), previous[j],
              (piece[j] == piece[i])
                  ? DotProduct(previousP[j] + sh - C0[i], previous[j])
                  : t * speedOf(C0[i], previous[j], i, piece[i]),
              C0[i], piece[i], t, convex, cone, coneA[i], coneB[i], arcC[i],
              arcGrows[i], Sf[i], Pf[i], Nf[i]);
        // the face's accumulated offset while the neighbour still carries
        // that face of the fan's material; one overtaken since holds another
        // plane, a face of another material moves at another speed
        const bool carried =
            face && carriesFace(j, N[j]) && piece[j] == piece[i];
        // a face of the fan's material that holds a face of a concave fan's
        // cone, which the speed at its own points has moved away from the
        // feature point (facesFromConcaveFans), gives that end its plane
        const bool endFace = !face && !convex && arcKnown[i] && good[j] &&
                             !fan[j] && piece[j] == piece[i] &&
                             (DotProduct(N[j], arcA[i]) > 1 - unitTolerance ||
                              DotProduct(N[j], arcB[i]) > 1 - unitTolerance) &&
                             DotProduct(Nf[i], N[j]) < 1 - unitTolerance;
        if (endFace)
          arc(i, Nf[i], DotProduct(Pf[i] - C0[i], Nf[i]), N[j],
              DotProduct(P[j] + sh - C0[i], N[j]), C0[i], piece[i], t, convex,
              cone, coneA[i], coneB[i], arcC[i], arcGrows[i], Sf[i], Pf[i],
              Nf[i]);
        if (face && (carried || timed)) {
          const T oy = carried ? DotProduct(P[j] + sh - C0[i], NRM0[j])
                               : t * speedOf(C0[i], NRM0[j], i, piece[i]);
          arc(i, Nf[i], DotProduct(Pf[i] - C0[i], Nf[i]), NRM0[j], oy, C0[i],
              piece[i], t, convex, cone, coneA[i], coneB[i], arcC[i],
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
        Bnew(arcB), RWnew(RW), CCnew(arcC);
    std::vector<char> cornerNew(cornerKnown);
    // whether a face is one of the faces of fan f: its plane at the rebuild
    // contains the feature point c and its normal ends the fan's cone
    auto faceOf = [&](const Vec &normal, const Vec &rebuildPoint, const Vec &c,
                      long f) {
      return (DotProduct(normal, arcA[f]) > 1 - unitTolerance ||
              DotProduct(normal, arcB[f]) > 1 - unitTolerance) &&
             std::abs(DotProduct(c - rebuildPoint, normal)) <= smallLength;
    };
    std::vector<char> over(n, 0), fanNew(fan), knownNew(arcKnown),
        convexNew(convexFan);
    std::vector<T> sinceNew(since);
    std::vector<int> materialNew(pieceAfter);
#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      std::array<Vec, slots> sN, sP, sW0, sWm, sC, sRW;
      std::array<T, slots> sV, sS;
      std::array<long, slots> sFan;
      std::array<int, slots> sPiece{}, sAfter{};
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
        const Vec Pj = Pf[j] + sh;
        // a face acts where the node's foot on it lies on its material, as
        // beside the wall of a mask and not in the undercut below it; a fan
        // on the gas and on its own material, as the pit in a substrate does
        // not cut into the mask
        if (stacked) {
          bool valid;
          if (fan[j]) {
            valid = Sf[i] > 0 || region[i] == piece[j];
          } else {
            const Vec foot = X[i] - Nj * DotProduct(X[i] - Pj, Nj);
            valid = materialBehind(foot, Nj) == piece[j];
          }
          if (!valid)
            continue;
        }
        T value = DotProduct(X[i] - Pj, Nj);
        Vec point = Pj, normal = Nj, witness = W[j] + sh,
            movedWitness = Wm[j] + sh, rebuildPoint = RW[j] + sh,
            feature = C0[j] + sh;
        bool support = false;
        T supportSpeed = 0;
        if (arcGrows[j]) {
          const Vec c = C0[j] + sh;
          // signed distance to the support plane at f along the cone
          auto planeAt = [&](T f, Vec &u, T &offset) {
            u = arcA[j] * (1 - f) + arcB[j] * f;
            const T length = Norm(u);
            if (length < tinyLength)
              return -std::numeric_limits<T>::max();
            u = u * (1 / length);
            offset = elapsed[j] * speedOf(c, u, j, piece[j]);
            return DotProduct(X[i] - c, u) - offset;
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
          point = c + bestU * bestOffset;
          witness = movedWitness = rebuildPoint = feature = c;
          supportSpeed = speedOf(c, bestU, j, piece[j]);
          support = true;
        } else if (fan[j] && arcKnown[j] && convexFan[j] && Sf[i] > 0) {
          Vec c = C0[j] + sh;
          Vec r = X[i] - c;
          if constexpr (D == 3) {
            Vec e = CrossProduct(arcA[j], arcB[j]);
            e = e * (1 / Norm(e));
            r = r - e * DotProduct(r, e);
          }
          const T length = Norm(r);
          if (length > smallLength) {
            Vec u = r * (1 / length);
            if (!within(u, arcA[j], arcB[j]))
              u = (DotProduct(u, arcA[j]) >= DotProduct(u, arcB[j])) ? arcA[j]
                                                                     : arcB[j];
            // static features only, whose support function needs no offsets
            // in other directions than the fan's own
            if (std::abs(DotProduct(Pj - c, Nj)) <= tinyLength &&
                speedOf(c, Nj, i, piece[j]) == 0 &&
                speedOf(c, u, i, piece[j]) == 0) {
              value = DotProduct(r, u);
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
          const Vec c = C0[j] + sh;
          Vec r = X[i] - c;
          if constexpr (D == 3) {
            Vec e = CrossProduct(arcA[j], arcB[j]);
            e = e * (1 / Norm(e));
            r = r - e * DotProduct(r, e);
          }
          const T length = Norm(r);
          const Vec u = r * (T(-1) / std::max(length, tinyLength));
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
                             DotProduct(Pj - c, Nj) -
                             t * speedOf(c, Nj, j, piece[j]);
            value = DotProduct(X[i] - c, u) - offset;
            normal = u;
            point = c + u * offset;
            witness = movedWitness = rebuildPoint = feature = c;
            supportSpeed = speedOf(c, u, j, piece[j]);
            support = true;
            beyond = true;
          }
        }
        // planes within the junction angle of the node's own act only where
        // their normal makes them move at another speed at the node's own
        // point, as a plane of a bevel beside the face a motion creates; a
        // rate that varies in space alone, as a noisy flux, makes no junction
        // between the planes of one face
        if (!beyond && DotProduct(Ni, normal) >= junction) {
          const T va = speedOf(Pf[i], normal, i, piece[j]),
                  vb = speedOf(Pf[i], Ni, i, piece[i]);
          if (std::abs(va - vb) <=
              T(1e-3) * std::max(std::abs(va), std::abs(vb)))
            continue;
        }
        bool distinct = true;
        for (int q = 0; q < count; ++q)
          distinct = distinct && DotProduct(sN[q], normal) < junction;
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
        sAfter[count] = support ? piece[j] : pieceAfter[j];
        sSupport[count] = support;
        sDepth[count] = S[j];
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
        const Vec r = c - rebuildPoint;
        if (std::abs(DotProduct(r, normal)) >
            rebuiltLength + rebuiltSlope * Norm(r))
          return false;
        if constexpr (D == 3) {
          const Vec z = CrossProduct(arcA[f], arcB[f]);
          if (std::abs(DotProduct(normal, z)) > rebuiltTilt * Norm(z))
            return false;
        }
        return within(normal, arcA[f], arcB[f]) ||
               DotProduct(normal, arcA[f]) >=
                   1 - rebuiltTilt * rebuiltTilt / 2 ||
               DotProduct(normal, arcB[f]) >= 1 - rebuiltTilt * rebuiltTilt / 2;
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
            std::abs(DotProduct(c - rebuildPoint, normal)) <= smallLength) {
          if (extended[f] && DotProduct(normal, past[f]) > 1 - unitTolerance)
            return 1;
          if (std::abs(DotProduct(normal, cut[f])) < T(1e-3))
            return -1;
        }
        // an interface through the feature point of a concave fan whose
        // front, swept in the directions against its cone, lies wholly on
        // the interface's gas side, as the underside of a mask that an
        // undercut lays open beside the arc around the mask's foot on a flat
        // substrate, where the interface lies along the end of the cone
        if (arcKnown[f] && !convexFan[f] && faceMaterial != piece[f] &&
            std::abs(DotProduct(c - rebuildPoint, normal)) <= smallLength &&
            DotProduct(normal, arcA[f]) <= unitTolerance &&
            DotProduct(normal, arcB[f]) <= unitTolerance)
          return -1;
        if (arcKnown[f] && (faceOf(normal, rebuildPoint, c, f) ||
                            supportOf(normal, rebuildPoint, c, f)))
          return convexFan[f] ? 1 : -1;
        // a concave fan, as the corner of an etched trench bottom, is a
        // corner of the convex region on the other side: a plane facing into
        // its cone that has moved past the feature point, as the floor going
        // down beside the corner, bounds that region there and meets the fan
        // concavely
        if (arcKnown[f] && !convexFan[f] && within(normal, arcA[f], arcB[f]) &&
            DotProduct(c - point, normal) > smallLength)
          return -1;
        if (convexFan[f]) {
          const T d = DotProduct(c - rebuildPoint, normal);
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
        // pieces of one fan meet with its convexity, also on the side the
        // front has passed since the rebuild
        else if (fa >= 0 && fb >= 0 && Norm(ca - cb) <= smallLength &&
                 convexFan[fa] == convexFan[fb])
          r = convexFan[fa] ? 1 : -1;
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
        // the budget of the piece that decides
        T budget = dt * std::abs(sV[active] - Vplane[i]);
        // where the support plane of a growing fan gives the value, it may
        // have run ahead of the node over the whole time since the rebuild,
        // as beside a face the motion creates, whose plane the node's own
        // plane passes
        if (sSupport[active] && sFan[active] >= 0 && arcGrows[sFan[active]])
          budget = std::max(budget, (elapsed[sFan[active]] + dt) *
                                        std::abs(sV[active] - Vplane[i]));
        // where its own plane would carry a solid node into the gas, a piece
        // of a solid node at most a grid spacing deeper that keeps it solid
        // holds it there beyond the budget, as the arc at the corner of an
        // etched trench bottom, which reaches a node beside the wall a step
        // after the node above took it, once the floor has passed: values lag
        // behind only deeper in the solid
        const bool kept = S[i] < 0 && Sf[i] >= 0 && value < 0 &&
                          sDepth[active] < 0 &&
                          -sDepth[active] <= std::abs(S[i]) + T(1);
        if (std::abs(value - Sf[i]) <= budget + smallLength || kept)
          break;
        dropped[active] = 1;
        active = -1;
      }
      if (active < 0)
        continue;
      // the plane keeps the neighbour's point on it, a point of the piece
      // itself, where its speed is taken, and the piece's material
      Snew[i] = value;
      Nnew[i] = sN[active];
      Pnew[i] = sP[active];
      Wnew[i] = sWm[active];
      RWnew[i] = sRW[active];
      materialNew[i] = sAfter[active];
      over[i] = 1;
      // a node that takes the plane of a fan, its own or a support plane,
      // becomes a fan of the same feature
      const long source = sFan[active];
      fanNew[i] = source >= 0;
      if (source >= 0) {
        Cnew[i] = sC[active];
        Anew[i] = arcA[source];
        Bnew[i] = arcB[source];
        knownNew[i] = arcKnown[source];
        convexNew[i] = convexFan[source];
        sinceNew[i] = since[source];
        CCnew[i] = arcC[source];
        cornerNew[i] = cornerKnown[source];
      }
    }

#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      if (over[i])
        HA[i] = Pnew[i];
      else if (good[i])
        HA[i] = P[i];
      if (!good[i] || over[i])
        FA[i] = Pnew[i];
      if (over[i]) {
        NRM[i] = Nnew[i];
        Dir[i] = Nnew[i];
        C0[i] = Cnew[i];
        since[i] = sinceNew[i];
        boundPlane[i] = 0;
        smooth[i] = 0;
        exposedState[i] = 0;
      } else if (turned[i]) {
        NRM[i] = N[i];
        Dir[i] = N[i];
      }
      if (stacked)
        pieceMaterial[i] = materialNew[i];
      W[i] = Wnew[i];
      RW[i] = RWnew[i];
      PL[i] = Pnew[i];
      NL[i] = Nnew[i];
      moved[i] = moved[i] || moving[i] || over[i];
      fan[i] = fanNew[i];
      arcA[i] = Anew[i];
      arcB[i] = Bnew[i];
      arcKnown[i] = knownNew[i];
      arcC[i] = CCnew[i];
      cornerKnown[i] = cornerNew[i];
      convexFan[i] = convexNew[i];
      good[i] = good[i] || over[i];
      S[i] = Snew[i];
    }
  }

  // interfaces of the material stack that etching lays open, as when it
  // undercuts a mask: across the interface of a lower level set from a solid
  // node, a node etched since it entered the band takes the smaller and the
  // solid node the larger signed distance to the plane of the interface,
  // from the lower level set (exact for planar interfaces), and both carry
  // that plane; only where the material on the solid side does not etch,
  // since a front etching through the interface is not held there, and for
  // the etched node only where its foot on the interface lies on that
  // material. Beyond the end of the cavity, as past the tip of an undercut,
  // the interface is no surface: a solid node there lies no closer to the
  // surface than the edge where the plane of the end meets the interface,
  // from the planes of the other material around it and the ends of the
  // cones of its fans. Pass by pass, solid nodes further from the cavity
  // beside one that took the interface take it on, with the plane that one
  // found where the lower level set holds no value there
  void exposeInterfaces() {
    const long n = long(S.size());
    for (int pass = 0; pass < 2 * bandWidth; ++pass) {
      long taken = 0;
      for (std::size_t l = 0; l + 1 < levelSets.size(); ++l) {
        // from the state before the pass, applied after it
        std::vector<char> take(n, 0), takeBound(n, 0);
        std::vector<T> newValue(n, T(0));
        std::vector<Vec> newPoint(n), newNormal(n);
        std::vector<int> newMaterial(n, 0);
        std::vector<Vec> newInterfacePoint(n), newInterfaceNormal(n);
#pragma omp parallel for
        for (long i = 0; i < n; ++i) {
          const bool etched = S[i] > 0 && entrySide[i] < 0, solid = S[i] < 0;
          if (!etched && !solid)
            continue;
          const T li = lowerValue(l, X[i]);
          if (li == T(0))
            continue;
          bool exposed = false;
          long carrier = -1;
          Vec carrierShift{};
          // the solid nodes across the interface from an etched node
          std::array<long, numberOfNeighbours> across;
          std::array<int, numberOfNeighbours> acrossK;
          int acrossCount = 0;
          for (int k = 0; k < numberOfNeighbours && (!exposed || etched); ++k) {
            const long j = neighbours[i][k];
            if (j < 0)
              continue;
            // a solid node beside one whose plane the interface gave
            if (solid && exposedState[j] && S[j] < 0 &&
                interfaceLevel[j] == int(l) &&
                material(X[j]) == material(X[i])) {
              exposed = true;
              carrier = j;
              carrierShift = shiftOf(i, k);
              continue;
            }
            const T lj = lowerValue(l, X[j] + shiftOf(i, k));
            if (li * lj >= T(0))
              continue;
            const bool here =
                etched ? S[j] < 0 : (S[j] > 0 && entrySide[j] < 0);
            exposed = exposed || here;
            if (here && etched) {
              across[acrossCount] = j;
              acrossK[acrossCount++] = k;
            }
          }
          if (!exposed)
            continue;
          Vec normal{};
          T distance = 0;
          if (lowerDefined(l, X[i])) {
            // the interface plane: from the Hermite data of a Euclidean lower
            // level set, otherwise from differences within its band, whose
            // values further out only give the side
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
                const Vec forward = X[i] + step, backward = X[i] - step;
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
            const T length = Norm(gradient);
            if (length < tinyLength)
              continue;
            // the normal points to the gas, on the side of the etched node
            const T towardsGas = etched ? sign(li) : -sign(li);
            normal = gradient * (towardsGas / length);
            distance = towardsGas * li / length;
          } else if (carrier >= 0) {
            // beyond the band of the lower level set, the plane a neighbour
            // nearer the cavity found
            normal = interfaceNormal[carrier];
            distance = DotProduct(
                X[i] - (interfacePoint[carrier] + carrierShift), normal);
            if (distance >= 0)
              continue;
          } else {
            continue;
          }
          T value = distance;
          bool bounded = false;
          Vec point = X[i] - normal * distance, planeNormal = normal;
          // an etched node's foot on the interface lies on the solid material
          // across it: on the solid side of the plane of one of the solid nodes
          // there, as below a mask and not beside its wall, where the interface
          // continues past the mask into the opening
          if (etched) {
            bool onSolid = false;
            for (int q = 0; q < acrossCount && !onSolid; ++q) {
              const long j = across[q];
              onSolid = DotProduct(point - (PL[j] + shiftOf(i, acrossK[q])),
                                   NL[j]) <= smallLength;
            }
            if (!onSolid)
              continue;
          }
          if (solid) {
            // the planes of the other material around the node, faces and the
            // ends of the cones of fans on the solid side; a static one, as the
            // wall that directional etching cuts below the edge of a mask,
            // bounds the cavity as well
            const int own = material(X[i]);
            // the interface counts where the node's foot on it lies in the
            // cavity, on the gas side of every plane of the front
            const Vec interfaceFoot = point;
            T bound = std::numeric_limits<T>::max(), footSide = bound;
            Vec boundNormal{}, endNormal{};
            auto consider = [&](const Vec &p, const Vec &d) {
              const T g = DotProduct(X[i] - p, d);
              const T side = DotProduct(interfaceFoot - p, d);
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
              const int m = nodeMaterial(j);
              if (m == own)
                continue;
              consider(PL[j] + sh, NL[j]);
              if (!fan[j] || !arcKnown[j] || convexFan[j])
                continue;
              const Vec &c = C0[j];
              const T correction = DotProduct(PL[j] - c, NL[j]) -
                                   since[j] * speedOf(c, NL[j], j, m);
              Vec a = arcA[j], b = arcB[j], past, end;
              const bool split = materialCone(c, m, a, b, past, end) != 0;
              auto plane = [&](const Vec &d) {
                const T offset = since[j] * speedOf(c, d, j, m) + correction;
                consider(c + d * offset + sh, d);
              };
              plane(arcA[j]);
              plane(arcB[j]);
              // the end of the cone, which materialCone sets only where the
              // interface splits it
              if (split)
                plane(end);
            }
            // beyond the end of the cavity the node lies nearest to the edge
            // where the plane of the end meets the interface
            const Vec along =
                endNormal - normal * DotProduct(endNormal, normal);
            const bool beyond = footSide < -smallLength && Norm(along) > T(0.5);
            if (beyond || bound < value) {
              bounded = true;
              if (beyond) {
                const T past = -footSide / Norm(along);
                value = -std::sqrt(distance * distance + past * past);
                bound = DotProduct(X[i] - interfaceFoot, endNormal) + footSide;
                boundNormal = endNormal;
              } else {
                value = bound;
              }
              planeNormal = boundNormal;
              // its point is the node's foot on the plane
              point = X[i] - planeNormal * bound;
            }
          }
          if (etched ? value >= S[i] : value <= S[i])
            continue;
          const Vec foot = X[i] - normal * distance;
          const int solidSide = materialBehind(foot, normal);
          if (speedOf(foot, normal, i, solidSide) < T(0))
            continue;
          take[i] = 1;
          takeBound[i] = solid && bounded;
          newValue[i] = value;
          newPoint[i] = point;
          newNormal[i] = planeNormal;
          // the interface is of the material on its solid side, which does not
          // etch; the plane of the end of a cavity that holds a solid node is
          // that node's material, static, and every step sets it anew
          newMaterial[i] = solid ? material(X[i]) : solidSide;
          newInterfacePoint[i] = foot;
          newInterfaceNormal[i] = normal;
        }
        for (long i = 0; i < n; ++i) {
          if (!take[i])
            continue;
          ++taken;
          S[i] = newValue[i];
          boundPlane[i] = takeBound[i];
          HA[i] = FA[i] = W[i] = PL[i] = RW[i] = newPoint[i];
          NRM[i] = Dir[i] = NL[i] = newNormal[i];
          good[i] = moved[i] = 1;
          fan[i] = 0;
          exposedState[i] = 1;
          pieceMaterial[i] = newMaterial[i];
          interfacePoint[i] = newInterfacePoint[i];
          interfaceNormal[i] = newInterfaceNormal[i];
          interfaceLevel[i] = int(l);
        }
      }
      if (taken == 0)
        break;
    }
  }

  // a band that is no valid level set at the end, where two neighbours along
  // an axis lie further apart than Check allows: of each such pair the node
  // further from the value the mesh of the planes gives it takes that value,
  // its planes moved to it. Returns whether any value changed
  bool repairFromMesh() {
    auto copy = SmartPointer<Domain<T, D>>::New(levelSet);
    Reinitialize<T, D>(copy, bandWidth).apply();
    std::vector<std::pair<hrleIndex, T>> mesh;
    for (viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType> it(
             copy->getDomain());
         !it.isFinished(); it.next())
      if (it.isDefined())
        mesh.emplace_back(it.getStartIndices(), it.getValue());
    std::sort(mesh.begin(), mesh.end(),
              [](const auto &a, const auto &b) { return a.first < b.first; });
    const long n = long(S.size());
    std::vector<T> target(n, std::numeric_limits<T>::quiet_NaN());
    for (long i = 0; i < n; ++i) {
      hrleIndex idx;
      for (int k = 0; k < D; ++k)
        idx[k] = static_cast<viennahrle::IndexType>(std::lround(X[i][k]));
      auto found = std::lower_bound(
          mesh.begin(), mesh.end(), idx,
          [](const auto &a, const hrleIndex &b) { return a.first < b; });
      if (found != mesh.end() && !(idx < found->first) && !(found->first < idx))
        target[i] = found->second;
    }
    // Check's rounding of a value to a grid point
    auto status = [](T v) {
      const int x = static_cast<int>(v);
      if (v >= 0)
        return (v <= T(x) + T(0.5)) ? x : x + 1;
      return (v >= T(x) - T(0.5)) ? x : x - 1;
    };
    bool changed = false;
    for (int pass = 0; pass < 3; ++pass) {
      std::vector<char> take(n, 0);
      for (long i = 0; i < n; ++i)
        for (int k = 0; k < 2 * D; ++k) {
          const long m = neighbours[i][k];
          if (m < 0 || std::abs(status(S[i]) - status(S[m])) <= 1)
            continue;
          const T gi =
              std::isnan(target[i]) ? T(-1) : std::abs(S[i] - target[i]);
          const T gm =
              std::isnan(target[m]) ? T(-1) : std::abs(S[m] - target[m]);
          if (gi >= gm && gi >= 0)
            take[i] = 1;
          else if (gm >= 0)
            take[m] = 1;
        }
      bool any = false;
      for (long i = 0; i < n; ++i) {
        if (!take[i])
          continue;
        const T delta = target[i] - S[i];
        S[i] = target[i];
        PL[i] = PL[i] - NL[i] * delta;
        HA[i] = HA[i] - Dir[i] * delta;
        FA[i] = FA[i] - NRM[i] * delta;
        any = changed = true;
      }
      if (!any)
        break;
    }
    return changed;
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
                   32 * int(boundPlane[i]) + 64 * int(smooth[i]) +
                   128 * int(exposedState[i]) +
                   256 * int(cornerKnown.size() == S.size() && cornerKnown[i]));
    pointData.insertReplaceScalarData(std::move(values), valueLabel);
    pointData.insertReplaceScalarData(std::move(sides), sideLabel);
    pointData.insertReplaceScalarData(std::move(flags), flagsLabel);
    pointData.insertReplaceScalarData(std::move(times), timeLabel);
    pointData.insertReplaceScalarData(
        typename PointData<T>::ScalarDataType(since), sinceLabel);
    typename PointData<T>::ScalarDataType materials(S.size());
    for (std::size_t i = 0; i < S.size(); ++i)
      materials[i] = T(pieceMaterial[i]);
    pointData.insertReplaceScalarData(std::move(materials), materialLabel);
    pointData.insertReplaceScalarData(
        typename PointData<T>::ScalarDataType(entrySide), entryLabel);
    if (arcC.size() == S.size()) {
      typename PointData<T>::VectorDataType cones(S.size());
      for (std::size_t i = 0; i < S.size(); ++i)
        for (int c = 0; c < 3; ++c)
          cones[i][c] = arcC[i][c];
      pointData.insertReplaceVectorData(std::move(cones), coneCLabel);
    }
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

  // the band moves with the surface: the dual contouring mesh of the planes
  // gives the nodes of the new band, every node of the old one keeps its
  // value and state, and a node new to the band, at most the rebuild
  // distance from the old ones, takes the state of an old node near it
  void rebuild() {
    writeBack();
    const auto before = nodeOrder;
    const std::vector<char> oRounded = roundedNodes();
    const std::vector<Vec> oldX = X;
    const std::vector<T> oS = S, oSide = side0, oSince = since,
                         oEntry = entrySide;
    const std::vector<Vec> oHA = HA, oFA = FA, oNRM = NRM, oNRM0 = NRM0,
                           oC0 = C0, oW = W, oDir = Dir, oPL = PL, oNL = NL,
                           oRW = RW, oA = arcA, oB = arcB, oC3 = arcC,
                           oIP = interfacePoint, oIN = interfaceNormal;
    const std::vector<char> oGood = good, oFan = fan, oMoved = moved,
                            oConvex = convexFan, oKnown = arcKnown,
                            oCorner = cornerKnown, oBound = boundPlane,
                            oSmooth = smooth, oExposed = exposedState;
    const std::vector<int> oMaterial = pieceMaterial, oIL = interfaceLevel;
    const T clock = sinceRebuild;
    Reinitialize<T, D>(levelSet, bandWidth).apply();
    ++numberOfRebuilds;
    collectNodes();
    const long n = long(S.size());
    // the values the mesh gives the new band
    const std::vector<T> mesh = S;
    // a node keeps the rounded part of its piece, a node beside a rounded
    // fan that takes it is on one
    std::vector<char> roundedNew(n, 0);
    for (auto *v : {&HA, &FA, &NRM, &NRM0, &C0, &W, &Dir, &PL, &NL, &RW, &arcA,
                    &arcB, &arcC, &interfacePoint, &interfaceNormal})
      v->assign(n, Vec{});
    for (auto *v : {&side0, &entrySide, &since})
      v->assign(n, T(0));
    for (auto *v : {&good, &fan, &moved, &convexFan, &arcKnown, &cornerKnown,
                    &boundPlane, &smooth, &exposedState})
      v->assign(n, 0);
    pieceMaterial.assign(n, 0);
    interfaceLevel.assign(n, -1);
    sinceRebuild = clock;
    // offsets to the old nodes around a new one, nearest first, as far as
    // the front moved since the last rebuild: the step that reaches the
    // rebuild distance can pass it, and a new node then lies further from
    // the old band, as on the cap of a rounded corner, whose new nodes would
    // otherwise find only the fans of the edges beside it
    std::vector<hrleIndex> around;
    {
      const int r = static_cast<int>(std::ceil(
          std::max(rebuildDistance, motionSinceRebuild) - smallLength));
      const int dz0 = (D == 3) ? -r : 0, dz1 = (D == 3) ? r : 0;
      for (int dz = dz0; dz <= dz1; ++dz)
        for (int dy = -r; dy <= r; ++dy)
          for (int dx = -r; dx <= r; ++dx) {
            hrleIndex o;
            o[0] = dx;
            o[1] = dy;
            if constexpr (D == 3)
              o[2] = dz;
            around.push_back(o);
          }
      std::stable_sort(around.begin(), around.end(),
                       [](const hrleIndex &a, const hrleIndex &b) {
                         int la = 0, lb = 0;
                         for (int i = 0; i < D; ++i) {
                           la += a[i] * a[i];
                           lb += b[i] * b[i];
                         }
                         return la < lb;
                       });
    }
    auto oldAt = [&](const hrleIndex &q) {
      auto found = std::lower_bound(
          before.begin(), before.end(), q,
          [](const auto &a, const hrleIndex &b) { return a.first < b; });
      return (found != before.end() && !(q < found->first) &&
              !(found->first < q))
                 ? found->second
                 : long(-1);
    };
    auto &pointData = levelSet->getPointData();
    const auto *closest = pointData.getVectorData(
        FromSurfaceMesh<T, D>::closestPointsLabel, true);
    const auto *normals = pointData.getVectorData(
        FromSurfaceMesh<T, D>::surfaceNormalsLabel, true);
#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      hrleIndex idx;
      for (int k = 0; k < D; ++k)
        idx[k] = static_cast<viennahrle::IndexType>(std::lround(X[i][k]));
      // the node itself in the old band; else, among the old nodes around
      // it, the plane that gives the node the value nearest to the mesh's,
      // unless the node's direction from the feature of a rounded convex fan
      // lies in the fan's cone, a 3D corner's before an edge's: a rounded
      // part grows alike in every direction of its cone, and its support
      // there gives the node's value exactly. Its radius is the old node's
      // distance from the feature less its signed value, also for an old
      // node inside the rounded part
      long j = oldAt(idx);
      const bool same = j >= 0;
      bool asFan = false;
      Vec fanPoint{}, fanR{};
      T fanValue = 0;
      int fanRank = 0;
      if (!same) {
        T best = std::numeric_limits<T>::max();
        for (const auto &o : around) {
          hrleIndex q = idx;
          for (int k = 0; k < D; ++k)
            q[k] += o[k];
          const long m = oldAt(q);
          if (m < 0)
            continue;
          const T v = DotProduct(oNL[m], X[i] - oPL[m]);
          if (!asFan && std::abs(v - S[i]) < best) {
            best = std::abs(v - S[i]);
            j = m;
          }
          if (!oFan[m] || !oGood[m] || !oRounded[m] || !oConvex[m] ||
              fanRank == 2)
            continue;
          const Vec c = oC0[m];
          if (D == 3 && oCorner[m]) {
            const Vec r = X[i] - c, rj = oldX[m] - c;
            const T length = Norm(r);
            if (length > smallLength &&
                inCornerCone(r * (1 / length), oA[m], oB[m], oC3[m])) {
              const T vf = length - (Norm(rj) - oS[m]);
              if (vf > 0) {
                fanRank = 2;
                j = m;
                asFan = true;
                fanPoint = c;
                fanR = r;
                fanValue = vf;
              }
            }
            continue;
          }
          if (!oKnown[m] || fanRank == 1)
            continue;
          Vec feature = c, r = X[i] - c, rj = oldX[m] - c;
          if constexpr (D == 3) {
            Vec e = CrossProduct(oA[m], oB[m]);
            const T le = Norm(e);
            if (le < tinyLength)
              continue;
            e = e * (1 / le);
            const T along = DotProduct(r, e);
            feature = feature + e * along;
            r = r - e * along;
            rj = rj - e * DotProduct(rj, e);
          }
          const T length = Norm(r);
          if (length < smallLength || !within(r * (1 / length), oA[m], oB[m]))
            continue;
          const T vf = length - (Norm(rj) - oS[m]);
          if (vf > 0) {
            fanRank = 1;
            j = m;
            asFan = true;
            fanPoint = feature;
            fanR = r;
            fanValue = vf;
          }
        }
      }
      // the mesh's closest point and the node's direction from it
      Vec cp{}, nm{};
      if (closest != nullptr && normals != nullptr)
        for (int k = 0; k < 3; ++k) {
          cp[k] = (*closest)[i][k] / gridDelta;
          nm[k] = (*normals)[i][k];
        }
      const Vec r = X[i] - cp;
      const T length = Norm(r);
      const Vec d = (length > smallLength) ? r * (sign(S[i]) / length) : nm;
      // the plane of a concave fan, as the cavity an etch rounds below the
      // foot of a mask, is its support plane in the direction the fan sampled
      // at that node, which in the node's own direction lies on the gas side
      // of the cavity: a node new to the band beside it becomes a fan of the
      // same feature with the mesh's plane in its own direction, where that
      // lies in the fan's cone, and else takes the mesh's plane. Beside a
      // convex fan the plane of the old node lies on the solid side, and the
      // fan's support over its cone, as the time fitted to the planes of a
      // growing edge, gives the node its value in the step
      bool meshFan = false;
      if (!same && !asFan && j >= 0 && oFan[j] && !oConvex[j]) {
        if (D == 3 && oCorner[j]) {
          meshFan = inCornerCone(d, oA[j], oB[j], oC3[j]);
        } else if (oKnown[j]) {
          const Vec z = CrossProduct(oA[j], oB[j]);
          meshFan =
              within(d, oA[j], oB[j]) &&
              (D == 2 || std::abs(DotProduct(d, z)) <= rebuiltTilt * Norm(z));
        }
        if (!meshFan)
          j = -1;
      }
      if (j < 0) {
        // no old node near, or a concave fan's outside its cone: the mesh's
        // closest point, the material as at the initialization, at its own
        // position for a solid node behind a corner of the mesh
        HA[i] = FA[i] = C0[i] = W[i] = PL[i] = RW[i] = cp;
        NRM[i] = NRM0[i] = NL[i] = Dir[i] = d;
        side0[i] = entrySide[i] = sign(S[i]);
        pieceMaterial[i] = (S[i] < 0 && length > smallLength &&
                            std::abs(DotProduct(r, nm)) < T(0.999) * length)
                               ? material(X[i])
                               : materialBehind(cp, d);
        good[i] = 1;
        convexFan[i] = S[i] > 0;
        continue;
      }
      HA[i] = oHA[j];
      FA[i] = oFA[j];
      NRM[i] = oNRM[j];
      NRM0[i] = oNRM0[j];
      C0[i] = oC0[j];
      W[i] = oW[j];
      Dir[i] = oDir[j];
      PL[i] = oPL[j];
      NL[i] = oNL[j];
      RW[i] = oRW[j];
      arcA[i] = oA[j];
      arcB[i] = oB[j];
      arcC[i] = oC3[j];
      side0[i] = oSide[j];
      good[i] = oGood[j];
      fan[i] = oFan[j];
      moved[i] = oMoved[j];
      convexFan[i] = oConvex[j];
      arcKnown[i] = oKnown[j];
      cornerKnown[i] = oCorner[j];
      boundPlane[i] = oBound[j];
      smooth[i] = oSmooth[j];
      since[i] = oSince[j];
      pieceMaterial[i] = oMaterial[j];
      if (same) {
        roundedNew[i] = oRounded[j];
        S[i] = oS[j];
        entrySide[i] = oEntry[j];
        exposedState[i] = oExposed[j];
        interfacePoint[i] = oIP[j];
        interfaceNormal[i] = oIN[j];
        interfaceLevel[i] = oIL[j];
        continue;
      }
      entrySide[i] = sign(S[i]);
      // a node new to the band beside a rounded fan: a fan of the same
      // feature, in its own direction from it
      if (asFan) {
        roundedNew[i] = 1;
        S[i] = fanValue;
        const Vec d = fanR * (1 / Norm(fanR));
        C0[i] = RW[i] = W[i] = fanPoint;
        Dir[i] = NL[i] = d;
        HA[i] = FA[i] = PL[i] = X[i] - d * S[i];
        side0[i] = sign(S[i]);
        moved[i] = 1;
        boundPlane[i] = 0;
        continue;
      }
      if (meshFan) {
        HA[i] = FA[i] = PL[i] = cp;
        NRM[i] = NRM0[i] = Dir[i] = NL[i] = d;
        side0[i] = sign(S[i]);
        moved[i] = 1;
        boundPlane[i] = 0;
        smooth[i] = 0;
        continue;
      }
      // a node new to the band: the plane of the old node
      S[i] = DotProduct(oNL[j], X[i] - oPL[j]);
      HA[i] = X[i] - oNL[j] * S[i];
      FA[i] = PL[i] = oPL[j];
      NRM[i] = Dir[i] = NL[i] = oNL[j];
      side0[i] = sign(S[i]);
      fan[i] = 0;
      arcKnown[i] = 0;
      cornerKnown[i] = 0;
      moved[i] = 1;
      boundPlane[i] = 0;
      smooth[i] = 0;
    }
    rebuiltRounded.swap(roundedNew);
    rebuiltBand = true;
    // the witness of a face node back at its foot on its plane: a witness
    // moves with its plane along the normal, while a face that a directional
    // velocity translates moves along it, so that over many rebuilds the
    // witnesses would leave the part of the piece around the node
    if (closest != nullptr && closest->size() == std::size_t(n))
#pragma omp parallel for
      for (long i = 0; i < n; ++i)
        if (!fan[i])
          for (int k = 0; k < 3; ++k)
            W[i][k] = (*closest)[i][k] / gridDelta;
    // a node whose kept value the mesh does not reproduce to the accuracy of
    // its planes, as a node beside a face the motion created that no node
    // held, takes the plane of the neighbour that gives it the mesh's value
    // to that accuracy; not on fans and rounded parts, whose mesh lies up to
    // a tenth of a grid spacing off. A node that takes the plane of a
    // concave fan becomes a fan of the same feature, as in the envelope: as
    // a face, a plane of the corner at the foot of a trench in the node's
    // direction would move with its own speed and cut a notch into the wall
    // once the corner has slid past
    {
      const std::vector<char> rounded = roundedNodes();
      for (int pass = 0; pass < 3; ++pass) {
        std::vector<long> from(n, -1);
        std::vector<T> to(n, T(0));
#pragma omp parallel for
        for (long i = 0; i < n; ++i) {
          if (std::abs(S[i] - mesh[i]) <= rebuiltLength || fan[i] || rounded[i])
            continue;
          T best = rebuiltLength;
          for (int k = 0; k < numberOfNeighbours; ++k) {
            const long m = neighbours[i][k];
            if (m < 0)
              continue;
            const T v = DotProduct(NL[m], X[i] - (PL[m] + shiftOf(i, k)));
            if (std::abs(v - mesh[i]) < best) {
              best = std::abs(v - mesh[i]);
              from[i] = m;
              to[i] = v;
            }
          }
        }
        long changed = 0;
        const auto NRM0b = NRM0, C00 = C0, W0 = W, RW0 = RW, NL0 = NL, PL0 = PL,
                   A0 = arcA, B0 = arcB;
        const auto good0 = good, convex0 = convexFan, fan0 = fan,
                   known0 = arcKnown, corner0 = cornerKnown;
        const auto C30 = arcC;
        const auto since0 = since;
        const auto material0 = pieceMaterial;
        for (long i = 0; i < n; ++i) {
          const long m = from[i];
          if (m < 0)
            continue;
          ++changed;
          S[i] = to[i];
          NRM0[i] = NRM0b[m];
          C0[i] = C00[m];
          W[i] = W0[m];
          RW[i] = RW0[m];
          FA[i] = PL[i] = PL0[m];
          NRM[i] = Dir[i] = NL[i] = NL0[m];
          HA[i] = X[i] - NL0[m] * S[i];
          arcA[i] = A0[m];
          arcB[i] = B0[m];
          good[i] = good0[m];
          convexFan[i] = convex0[m];
          since[i] = since0[m];
          pieceMaterial[i] = material0[m];
          exposedState[i] = 0;
          side0[i] = sign(S[i]);
          const bool concave = fan0[m] && !convex0[m];
          fan[i] = concave;
          arcKnown[i] = concave && known0[m];
          cornerKnown[i] = concave && corner0[m];
          if (concave)
            arcC[i] = C30[m];
          moved[i] = 1;
          boundPlane[i] = 0;
          smooth[i] = 0;
        }
        if (changed == 0)
          break;
      }
    }
    // the side at the rebuild is the side now: the side of undecided
    // junctions counts from this rebuild
    for (long i = 0; i < n; ++i)
      side0[i] = sign(S[i]);
    motionSinceRebuild = 0.;
    writeBack();
  }

  // nodes whose plane is a tangent plane of a rounded part: convex fans
  // whose speed, the same in every direction of their cone, at the node and
  // at its plane, moves the surface outwards, as a deposition rounds an edge
  // into an arc, concave fans whose speed moves it inwards, as the cavity an
  // etching rounds around a static corner, and nodes that a previous apply()
  // left on such a part. A rate that varies in space, as a flux or a noisy
  // view factor, or with the direction, as directional deposition, rounds
  // nothing exactly
  std::vector<char> roundedNodes() const {
    if (rebuiltBand && rebuiltRounded.size() == S.size())
      return rebuiltRounded;
    const long n = long(S.size());
    std::vector<char> rounded(n, 0);
#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      if (long(smooth.size()) == n && smooth[i]) {
        rounded[i] = 1;
        continue;
      }
      if (!fan[i] || !good[i])
        continue;
      const int m = nodeMaterial(i);
      const T v = speedOf(C0[i], NL[i], i, m);
      if (convexFan[i] ? v <= T(0) : v >= T(0))
        continue;
      auto same = [&](const Vec &p, const Vec &d) {
        return std::abs(speedOf(p, d, i, m) - v) <= T(1e-6) * std::abs(v);
      };
      bool uniform = same(C0[i], NRM0[i]) && same(C0[i], Dir[i]) &&
                     same(X[i], NL[i]) && same(PL[i], NL[i]);
      if (arcKnown[i])
        uniform = uniform && same(C0[i], arcA[i]) && same(C0[i], arcB[i]);
      rounded[i] = uniform;
    }
    return rounded;
  }

  // a convex fan takes its offsets as its elapsed time times the speeds at
  // its feature point now, exact only while the velocity keeps its shape. A
  // fan that rounds its edge, as the arc a deposition grows over the edge of
  // a mask, keeps its feature point under the growing film over rebuilds; a
  // speed that changes its shape there, as a flux that a closing opening
  // shadows or a rate whose isotropic part stops, collapses the arc. A node
  // of such a fan whose plane has moved more than a grid spacing from the
  // feature point carries that plane as a face from then on, through its
  // foot, and moves with the speed at its own surface point. The fan rounds
  // where the middle of its arc moves slower than the edge, (va + vb) / |a +
  // b| from the speeds of the arc's ends; a fan whose planes pass through
  // one moving corner, as the silhouette edge of a directional deposition
  // sweeping out a face, stays a fan, and so do fans of a rate alike in every
  // direction (roundedNodes), which round exactly
  void facesFromRoundingFans() {
    const long n = long(S.size());
    const std::vector<char> rounded = roundedNodes();
#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      if (!fan[i] || !convexFan[i] || !good[i] || !arcKnown[i] || rounded[i] ||
          DotProduct(PL[i] - C0[i], NL[i]) <= T(1))
        continue;
      const Vec middle = arcA[i] + arcB[i];
      const T length = Norm(middle);
      if (length < smallLength)
        continue;
      const int m = nodeMaterial(i);
      const T va = speedOf(C0[i], arcA[i], i, m),
              vb = speedOf(C0[i], arcB[i], i, m);
      if (!(speedOf(C0[i], middle * (1 / length), i, m) <
            (va + vb) / length -
                T(1e-3) * std::max(std::abs(va), std::abs(vb))))
        continue;
      const Vec d = NL[i], foot = X[i] - NL[i] * S[i];
      HA[i] = FA[i] = PL[i] = C0[i] = W[i] = RW[i] = foot;
      NRM[i] = NRM0[i] = Dir[i] = NL[i] = d;
      fan[i] = 0;
      arcKnown[i] = 0;
      cornerKnown[i] = 0;
      since[i] = 0;
      moved[i] = 1;
    }
  }

  // a concave fan takes its offsets as its time times the speeds at its
  // feature point, exact for a velocity alike in space, as for the cavity an
  // etch rounds below the foot of a mask. Its feature point stays where its
  // faces met, while a velocity that varies in space moves each face with
  // the speed at its own points, as the floor of a trench that slides the
  // corner at its foot down the wall. A node beside one of the faces holds
  // the plane of that face, the end of the cone in its direction. Once the
  // ends of the cone have moved more than a grid spacing in the time since
  // the fan's planes were set, such a node carries that plane as a face,
  // through its foot, and moves with the speed at its own surface point. A
  // node behind the plane of the other end, that of a face neighbour holding
  // it or else the fan's own, lies behind the moved corner and stays a fan,
  // and so do the fans of rounded parts (roundedNodes)
  void facesFromConcaveFans() {
    const long n = long(S.size());
    const std::vector<char> rounded = roundedNodes();
    std::vector<char> face(n, 0);
#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      if (!fan[i] || convexFan[i] || !good[i] || !arcKnown[i] || rounded[i])
        continue;
      const bool atA = DotProduct(NL[i], arcA[i]) >= 1 - unitTolerance;
      if (!atA && DotProduct(NL[i], arcB[i]) < 1 - unitTolerance)
        continue;
      const int m = nodeMaterial(i);
      const T va = std::abs(speedOf(C0[i], arcA[i], i, m)),
              vb = std::abs(speedOf(C0[i], arcB[i], i, m));
      if (since[i] * std::max(va, vb) <= T(1))
        continue;
      const Vec &other = atA ? arcB[i] : arcA[i];
      Vec otherPoint = C0[i] + other * (since[i] * speedOf(C0[i], other, i, m));
      for (int k = 0; k < numberOfNeighbours; ++k) {
        const long j = neighbours[i][k];
        if (j >= 0 && !fan[j] && good[j] &&
            pieceMaterial[j] == pieceMaterial[i] &&
            DotProduct(NL[j], other) >= 1 - unitTolerance) {
          otherPoint = PL[j] + shiftOf(i, k);
          break;
        }
      }
      face[i] = DotProduct(X[i] - otherPoint, other) >= -smallLength;
    }
#pragma omp parallel for
    for (long i = 0; i < n; ++i) {
      if (!face[i])
        continue;
      const Vec d = NL[i], foot = X[i] - NL[i] * S[i];
      HA[i] = FA[i] = PL[i] = C0[i] = W[i] = RW[i] = foot;
      NRM[i] = NRM0[i] = Dir[i] = NL[i] = d;
      fan[i] = 0;
      arcKnown[i] = 0;
      cornerKnown[i] = 0;
      since[i] = 0;
      moved[i] = 1;
    }
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
  // nodes on a rounded part (roundedNodes), whose cells dual contouring
  // places halfway between the planes' point and the mass point
  static constexpr char smoothLabel[] = "HermiteSmooth";
  static constexpr char coneALabel[] = "HermiteConeA";
  static constexpr char coneBLabel[] = "HermiteConeB";
  static constexpr char timeLabel[] = "HermiteTime";
  static constexpr char sinceLabel[] = "HermiteSince";
  static constexpr char materialLabel[] = "HermiteMaterial";
  static constexpr char coneCLabel[] = "HermiteConeC";
  static constexpr char entryLabel[] = "HermiteEntrySide";

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

  /// Front motion in grid spacings after which the band moves with the
  /// surface, its nodes from the dual contouring mesh of the planes, default
  /// 2. It must stay well inside the band.
  void setRebuildDistance(T distance) { rebuildDistance = distance; }

  /// Level set width of the rebuilt band, default 12 (6 grid spacings on
  /// either side).
  void setBandWidth(int width) { bandWidth = width; }

  /// Whether apply() stops after one time step, as Advect's single step for
  /// velocities that are updated between steps. Default false.
  void setSingleStep(bool passedSingleStep) { singleStep = passedSingleStep; }

  /// Whether the velocity field knows the speeds only for the normals the
  /// surface has at its points, as fluxes from ray tracing: a concave corner
  /// that moves is then the meeting point of its faces, and a face turns to
  /// the normal of the moved points around it. Default false.
  void setSurfaceSpeedsOnly(bool passed) { surfaceSpeedsOnly = passed; }

  /// Whether apply() repairs the band from the dual contouring mesh of its
  /// planes once the whole advection time is done if the band is no valid
  /// level set (Check), default true. Values further from the surface can lag
  /// behind edges the motion created, as at the tip of an undercut; the
  /// repair leaves the surface as it is. Without it, a further apply()
  /// continues from the band state.
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
    // the velocities hold the data of the band as this call finds it
    rebuiltBand = false;
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
      if (motionSinceRebuild >= rebuildDistance - tinyLength)
        rebuild();
      const long n = long(S.size());
      // the point of each face's plane, where its speed is asked and from
      // which it moves and turns, at the node's own foot on that plane: a
      // node that took a plane over from a neighbour would otherwise keep
      // the neighbour's point, which can lie grid spacings away, past the end
      // of the face. Fans keep their feature point, nodes holding an
      // interface plane theirs, and so does a node whose foot lies on another
      // material
#pragma omp parallel for
      for (long i = 0; i < n; ++i) {
        if (fan[i] || exposedState[i] || boundPlane[i])
          continue;
        Vec &a = good[i] ? HA[i] : FA[i];
        const Vec &nn = good[i] ? Dir[i] : NRM[i];
        const Vec foot = X[i] - nn * DotProduct(X[i] - a, nn);
        if (Norm(foot - a) <= smallLength)
          continue;
        if (levelSets.size() > 1 && materialBehind(foot, nn) != nodeMaterial(i))
          continue;
        a = foot;
      }
      std::vector<T> Vd(n, T(0)), Vn(n, T(0));
      T maxSpeed = 0.;
#pragma omp parallel for reduction(max : maxSpeed)
      for (long i = 0; i < n; ++i) {
        Vd[i] = speedOf(HA[i], Dir[i], i, nodeMaterial(i));
        Vn[i] = good[i] ? Vd[i] : speedOf(FA[i], NRM[i], i, nodeMaterial(i));
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
      if (maxSpeed > 0) {
        step(dt, Vd, Vn);
        facesFromRoundingFans();
        facesFromConcaveFans();
      }
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
      }
      if (singleStep)
        break;
    }
    const bool done =
        advectedTime >=
        advectionTime - tinyLength * std::max(T(1), advectionTime);
    if (!rebuiltLast)
      writeBack();
    // values further from the surface can lag behind edges the motion
    // created, as at the tip of an undercut, also right after a rebuild,
    // which keeps the values of the old band. Of a band that is no valid
    // level set, the nodes of the neighbour pairs that Check finds too far
    // apart take the values of the mesh of the planes; only where that does
    // not suffice, the band is rebuilt from the mesh
    if (done && rebuildAtEnd) {
      Check<T, D> check(levelSet);
      check.apply();
      if (!check.isValid() && boundFarValues()) {
        writeBack();
        check.apply();
      }
      if (!check.isValid() && repairFromMesh()) {
        writeBack();
        check.apply();
      }
      if (!check.isValid()) {
        Reinitialize<T, D>(levelSet, bandWidth).apply();
        ++numberOfRebuilds;
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
