#pragma once

#include <lsDomain.hpp>

#include <hrleSparseIterator.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <vector>

namespace viennals {

using namespace viennacore;

namespace detail {

// Key of a grid index.  The per-node maps (pressure, stress and
// concentration of the oxidation solvers, the mask carry velocity) are keyed
// by this value alone, so it must be injective: the coordinates are packed
// into 64 bits, 32 per axis in 2-D and 21 in 3-D (|index| < 2^20).  The
// hash_combine it replaces collided for 9 in 10 indices of a 5 nm trench
// domain ((x, y) with (x - 1, y + 67)), and nodes read each other's stress.
template <int D>
inline std::size_t gridIndexHash(const viennahrle::Index<D> &index) {
  static_assert(D == 2 || D == 3);
  constexpr unsigned bits = D == 3 ? 21 : 32;
  constexpr std::uint64_t mask = (std::uint64_t(1) << bits) - 1;
  constexpr std::int64_t offset = std::int64_t(1) << (bits - 1);
  std::uint64_t key = 0;
  for (unsigned i = 0; i < static_cast<unsigned>(D); ++i)
    key = (key << bits) | (static_cast<std::uint64_t>(
                               static_cast<std::int64_t>(index[i]) + offset) &
                           mask);
  return static_cast<std::size_t>(key);
}

// Hash functor for viennahrle::Index<D> (the key above).
template <int D> struct IndexTypeHasher {
  std::size_t operator()(const viennahrle::Index<D> &idx) const {
    return gridIndexHash<D>(idx);
  }
};

template <class T> inline Vec3D<T> vecScaled(const Vec3D<T> &source, T factor) {
  Vec3D<T> result{0., 0., 0.};
  for (unsigned i = 0; i < 3; ++i)
    result[i] = source[i] * factor;
  return result;
}

template <class T>
inline Vec3D<T> vecAdd(const Vec3D<T> &a, const Vec3D<T> &b) {
  Vec3D<T> result{0., 0., 0.};
  for (unsigned i = 0; i < 3; ++i)
    result[i] = a[i] + b[i];
  return result;
}

template <class T>
inline Vec3D<T> vecSubtract(const Vec3D<T> &a, const Vec3D<T> &b) {
  Vec3D<T> result{0., 0., 0.};
  for (unsigned i = 0; i < 3; ++i)
    result[i] = a[i] - b[i];
  return result;
}

template <class T>
inline void vecAddTo(Vec3D<T> &target, const Vec3D<T> &source) {
  for (unsigned i = 0; i < 3; ++i)
    target[i] += source[i];
}

/// Clamp HRLE far-field sentinels (±DBL_MAX) to ±1 before differencing to
/// prevent DBL_MAX² overflow that silently returns the zero vector.
template <class T> inline T clampLevelSetPhi(T v) {
  return v > T(1) ? T(1) : (v < T(-1) ? T(-1) : v);
}

template <class T>
inline T levelSetCrossingDistance(T insidePhi, T outsidePhi,
                                  T minBoundaryFraction, T gridDelta) {
  const T denom = std::abs(insidePhi) + std::abs(outsidePhi);
  if (denom <= std::numeric_limits<T>::epsilon())
    return gridDelta;
  return std::max(minBoundaryFraction * gridDelta,
                  gridDelta * std::abs(insidePhi) / denom);
}

} // namespace detail

/// Common Cartesian-grid infrastructure shared by the oxidation solver
/// classes (diffusion, deformation). Provides the node lookup
/// table, grid extents, and the core grid utility methods that operate on them.
template <class T, int D> class OxidationSolverBase {
protected:
  using IndexType = viennahrle::Index<D>;
  using ConstSparseIterator =
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>;

  static constexpr std::size_t noNode = std::numeric_limits<std::size_t>::max();
  std::vector<std::size_t> nodeLookupFlat;
  IndexType minIndex{};
  IndexType maxIndex{};
  std::array<std::size_t, D> extents{};
  std::array<std::size_t, D> strides{};
  T gridDelta = 1.;
  // Box sides on a reflective domain boundary: the mirror plane through the
  // outermost grid points.
  std::array<bool, D> mirrorLow_{}, mirrorHigh_{};

  template <class GridType> void setMirrorSides(const GridType &grid) {
    for (unsigned i = 0; i < D; ++i) {
      const bool reflective = grid.getBoundaryConditions(i) ==
                              viennahrle::BoundaryType::REFLECTIVE_BOUNDARY;
      mirrorLow_[i] = reflective && minIndex[i] == grid.getMinGridPoint(i);
      mirrorHigh_[i] = reflective && maxIndex[i] == grid.getMaxGridPoint(i);
    }
  }

  // The stencil side (direction, offset) of the grid point at index: beyond
  // a mirror side the neighbour is the mirror image of the opposite one, so
  // the side is -offset; else offset.
  int stencilOffset(const IndexType &index, unsigned direction,
                    int offset) const {
    IndexType neighbor = index;
    neighbor[direction] += offset;
    if (inBounds(neighbor) ||
        !(offset < 0 ? mirrorLow_[direction] : mirrorHigh_[direction]))
      return offset;
    IndexType opposite = index;
    opposite[direction] -= offset;
    return inBounds(opposite) ? -offset : offset;
  }

  struct GridBounds {
    IndexType min{};
    IndexType max{};
    bool foundDefinedPoints = false;
  };

  bool crosses(T a, T b) const {
    // Use a small tolerance consistent with isInsideOxide so that grid points
    // that land exactly on an interface (phi ≈ 0 from floating-point rounding
    // in GeometricAdvect) are treated as surface points in boundary detection.
    constexpr T eps = T(1e-9);
    return (a <= eps && b >= -eps) || (a >= -eps && b <= eps);
  }

  T valueAt(ConstSparseIterator &it, const IndexType &index) const {
    it.goToIndices(index);
    return it.getValue();
  }

  // Gradient direction of a level set at a grid point.  Per axis the central
  // difference where both neighbours are defined, one-sided (with the point
  // itself) where one is.  Clamping phi to one grid unit before differencing
  // -- the earlier form, kept only where the point itself is outside the
  // narrow band -- clipped the steeper component next to the interface and
  // turned the normals of tilted fronts towards the diagonal: a 30-degree
  // front by 2.6 degrees at phi = 0.3 and 5.4 at 0.6, at any grid spacing.
  // onGrid(idx, ok) maps a neighbour onto the grid (e.g. a reflection) and
  // clears ok if it is off the grid.  Zero if degenerate.
  template <class OnGrid>
  Vec3D<T> definedLevelSetGradient(ConstSparseIterator &it,
                                   const IndexType &index,
                                   OnGrid onGrid) const {
    // value at idx (the point itself if idx is off the grid), and whether
    // it is a defined value
    auto sample = [&](const IndexType &idx, T &value) {
      bool ok = true;
      const IndexType q = onGrid(idx, ok);
      it.goToIndices(ok ? q : index);
      value = it.getValue();
      return ok && it.isDefined() && std::isfinite(value);
    };
    T centre = T(0);
    const bool centreDefined = sample(index, centre);
    Vec3D<T> gradient{T(0), T(0), T(0)};
    for (unsigned d = 0; d < D; ++d) {
      IndexType plus = index, minus = index;
      plus[d] += 1;
      minus[d] -= 1;
      T vp = T(0), vm = T(0);
      const bool okP = sample(plus, vp), okM = sample(minus, vm);
      if (okP && okM)
        gradient[d] = T(0.5) * (vp - vm);
      else if (centreDefined && okP)
        gradient[d] = vp - centre;
      else if (centreDefined && okM)
        gradient[d] = centre - vm;
      else
        gradient[d] = T(0.5) * (detail::clampLevelSetPhi(vp) -
                                detail::clampLevelSetPhi(vm));
    }
    return gradient;
  }

  bool inBounds(const IndexType &index) const {
    for (unsigned i = 0; i < D; ++i)
      if (index[i] < minIndex[i] || index[i] > maxIndex[i])
        return false;
    return true;
  }

  void initNodeLookup() {
    std::size_t total = 1;
    for (unsigned i = 0; i < D; ++i)
      total *= extents[i];
    nodeLookupFlat.assign(total, noNode);
  }

  std::size_t lookupNode(const IndexType &index) const {
    if (!inBounds(index))
      return noNode;
    return nodeLookupFlat[linearIndex(index)];
  }

  std::size_t linearIndex(const IndexType &index) const {
    if (!inBounds(index))
      return std::numeric_limits<std::size_t>::max();
    std::size_t result = 0;
    for (unsigned i = 0; i < D; ++i)
      result += static_cast<std::size_t>(index[i] - minIndex[i]) * strides[i];
    return result;
  }

  bool increment(IndexType &index) const {
    for (unsigned i = 0; i < D; ++i) {
      if (index[i] < maxIndex[i]) {
        ++index[i];
        return true;
      }
      index[i] = minIndex[i];
    }
    return false;
  }

  // Searches with expanding radius shells so that RK2 Stage 2 can find oxide
  // nodes even after the surface has advanced several grid cells beyond the
  // band.
  std::size_t findNearbyNode(const IndexType &index) const {
    for (int radius = 1; radius <= 4; ++radius) {
      T bestDistance2 = std::numeric_limits<T>::max();
      std::size_t bestNode = std::numeric_limits<std::size_t>::max();

      IndexType offset{};
      offset.fill(-radius);
      while (true) {
        IndexType candidate = index;
        T distance2 = 0.;
        for (unsigned i = 0; i < D; ++i) {
          candidate[i] += offset[i];
          distance2 += static_cast<T>(offset[i] * offset[i]);
        }

        if (distance2 > 0 && inBounds(candidate)) {
          const std::size_t foundId = nodeLookupFlat[linearIndex(candidate)];
          if (foundId != noNode && distance2 < bestDistance2) {
            bestDistance2 = distance2;
            bestNode = foundId;
          }
        }

        unsigned dim = 0;
        for (; dim < D; ++dim) {
          if (offset[dim] < radius) {
            ++offset[dim];
            break;
          }
          offset[dim] = -radius;
        }
        if (dim == D)
          break;
      }

      if (bestNode != std::numeric_limits<std::size_t>::max())
        return bestNode;
    }
    return std::numeric_limits<std::size_t>::max();
  }

  // All nodes at the smallest distance from index (the first shell, up to 4
  // cells, holding a node that passes accept(j)).  Averaging them instead of
  // taking the first one found keeps a result free of the scan order (-y
  // first, then -x), which moved the free surface of a stress-free oxide
  // cylinder by up to 4 nm in y in 0.3 hr.
  template <class Accept>
  std::vector<std::size_t> findNearbyNodes(const IndexType &index,
                                           Accept accept) const {
    std::vector<std::size_t> best;
    for (int radius = 1; radius <= 4 && best.empty(); ++radius) {
      long bestDistance2 = std::numeric_limits<long>::max();
      IndexType offset{};
      offset.fill(-radius);
      while (true) {
        IndexType candidate = index;
        long distance2 = 0;
        for (unsigned i = 0; i < D; ++i) {
          candidate[i] += offset[i];
          distance2 += static_cast<long>(offset[i]) * offset[i];
        }
        if (distance2 > 0 && distance2 <= bestDistance2 &&
            inBounds(candidate)) {
          const std::size_t foundId = nodeLookupFlat[linearIndex(candidate)];
          if (foundId != noNode && accept(foundId)) {
            if (distance2 < bestDistance2) {
              best.clear();
              bestDistance2 = distance2;
            }
            best.push_back(foundId);
          }
        }
        unsigned dim = 0;
        for (; dim < D; ++dim) {
          if (offset[dim] < radius) {
            ++offset[dim];
            break;
          }
          offset[dim] = -radius;
        }
        if (dim == D)
          break;
      }
    }
    return best;
  }

  bool initializeGridFromInterfaces(
      SmartPointer<Domain<T, D>> reactionInterface,
      SmartPointer<Domain<T, D>> ambientInterface,
      SmartPointer<Domain<T, D>> maskInterface, bool useRequestedBounds,
      const IndexType &requestedMinIndex, const IndexType &requestedMaxIndex,
      std::size_t maxGridPoints, const std::string &solverName) {
    auto &reactionGrid = reactionInterface->getGrid();
    auto &ambientGrid = ambientInterface->getGrid();
    gridDelta = reactionGrid.getGridDelta();

    if (std::abs(gridDelta - ambientGrid.getGridDelta()) >
            std::numeric_limits<T>::epsilon() ||
        (maskInterface != nullptr &&
         std::abs(gridDelta - maskInterface->getGrid().getGridDelta()) >
             std::numeric_limits<T>::epsilon())) {
      Logger::getInstance()
          .addError(solverName + ": Interface grid deltas must match.")
          .print();
      return false;
    }

    IndexType gridMin = reactionGrid.getMinGridPoint();
    IndexType gridMax = reactionGrid.getMaxGridPoint();
    for (unsigned i = 0; i < D; ++i) {
      gridMin[i] = std::max(gridMin[i], ambientGrid.getMinGridPoint(i));
      gridMax[i] = std::min(gridMax[i], ambientGrid.getMaxGridPoint(i));
      if (maskInterface != nullptr) {
        gridMin[i] =
            std::max(gridMin[i], maskInterface->getGrid().getMinGridPoint(i));
        gridMax[i] =
            std::min(gridMax[i], maskInterface->getGrid().getMaxGridPoint(i));
      }
    }

    auto bounds = definedPointBounds(reactionInterface);
    mergeDefinedPointBounds(bounds, ambientInterface);
    if (maskInterface != nullptr)
      mergeDefinedPointBounds(bounds, maskInterface);

    minIndex = bounds.foundDefinedPoints ? bounds.min : gridMin;
    maxIndex = bounds.foundDefinedPoints ? bounds.max : gridMax;
    applyBoundsPaddingAndClamp(minIndex, maxIndex, gridMin, gridMax);

    std::size_t numGridPoints = 1;
    for (unsigned i = 0; i < D; ++i) {
      if (useRequestedBounds) {
        minIndex[i] = std::max(minIndex[i], requestedMinIndex[i]);
        maxIndex[i] = std::min(maxIndex[i], requestedMaxIndex[i]);
      }
      if (maxIndex[i] < minIndex[i]) {
        Logger::getInstance()
            .addError(solverName + ": Cartesian solve region is empty.")
            .print();
        return false;
      }
      extents[i] = static_cast<std::size_t>(maxIndex[i] - minIndex[i] + 1);
      numGridPoints *= extents[i];
      strides[i] = (i == 0) ? 1 : strides[i - 1] * extents[i - 1];
    }

    if (numGridPoints > maxGridPoints) {
      Logger::getInstance()
          .addError(solverName +
                    ": Cartesian solve region exceeds maxGridPoints.")
          .print();
      return false;
    }
    setMirrorSides(reactionGrid);
    return true;
  }

  bool initializeGridFromMask(SmartPointer<Domain<T, D>> maskInterface,
                              bool useRequestedBounds,
                              const IndexType &requestedMinIndex,
                              const IndexType &requestedMaxIndex,
                              std::size_t maxGridPoints,
                              const std::string &solverName) {
    auto &grid = maskInterface->getGrid();
    gridDelta = grid.getGridDelta();

    auto bounds = definedPointBounds(maskInterface);
    minIndex = bounds.foundDefinedPoints ? bounds.min : grid.getMinGridPoint();
    maxIndex = bounds.foundDefinedPoints ? bounds.max : grid.getMaxGridPoint();
    applyBoundsPaddingAndClamp(minIndex, maxIndex, grid.getMinGridPoint(),
                               grid.getMaxGridPoint());

    std::size_t numGridPoints = 1;
    for (unsigned i = 0; i < D; ++i) {
      if (useRequestedBounds) {
        minIndex[i] = std::max(minIndex[i], requestedMinIndex[i]);
        maxIndex[i] = std::min(maxIndex[i], requestedMaxIndex[i]);
      }
      if (maxIndex[i] < minIndex[i]) {
        Logger::getInstance()
            .addError(solverName + ": Cartesian solve region is empty.")
            .print();
        return false;
      }
      extents[i] = static_cast<std::size_t>(maxIndex[i] - minIndex[i] + 1);
      numGridPoints *= extents[i];
      strides[i] = (i == 0) ? 1 : strides[i - 1] * extents[i - 1];
    }

    if (numGridPoints > maxGridPoints) {
      Logger::getInstance()
          .addError(solverName +
                    ": Cartesian solve region exceeds maxGridPoints.")
          .print();
      return false;
    }
    setMirrorSides(grid);
    return true;
  }

private:
  GridBounds definedPointBounds(SmartPointer<Domain<T, D>> levelSet) const {
    GridBounds bounds;
    ConstSparseIterator it(levelSet->getDomain());
    for (; !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;

      const auto &index = it.getStartIndices();
      if (!bounds.foundDefinedPoints) {
        bounds.min = index;
        bounds.max = index;
        bounds.foundDefinedPoints = true;
      } else {
        for (unsigned i = 0; i < D; ++i) {
          bounds.min[i] = std::min(bounds.min[i], index[i]);
          bounds.max[i] = std::max(bounds.max[i], index[i]);
        }
      }
    }
    return bounds;
  }

  void mergeDefinedPointBounds(GridBounds &target,
                               SmartPointer<Domain<T, D>> levelSet) const {
    const auto source = definedPointBounds(levelSet);
    if (!source.foundDefinedPoints)
      return;
    if (!target.foundDefinedPoints) {
      target = source;
      return;
    }
    for (unsigned i = 0; i < D; ++i) {
      target.min[i] = std::min(target.min[i], source.min[i]);
      target.max[i] = std::max(target.max[i], source.max[i]);
    }
  }

  static void applyBoundsPaddingAndClamp(IndexType &minIndex,
                                         IndexType &maxIndex,
                                         const IndexType &gridMin,
                                         const IndexType &gridMax) {
    constexpr int padding = 4;
    for (unsigned i = 0; i < D; ++i) {
      // Clamp to the physical grid before adding/subtracting padding.
      // lsInterior on an INFINITE_BOUNDARY level-set can leave far-field
      // sentinel entries (near std::numeric_limits<IndexType>::min/max) in
      // the HRLE structure.  definedPointBounds() visits them and records
      // the extreme index.  Subtracting padding from that sentinel directly
      // overflows, turning the tiny sentinel into a large positive value and
      // making extents[i] enormous.  Clamping first makes the arithmetic safe.
      minIndex[i] =
          std::max(gridMin[i], std::max(gridMin[i], minIndex[i]) - padding);
      maxIndex[i] =
          std::min(gridMax[i], std::min(gridMax[i], maxIndex[i]) + padding);
    }
  }
};

} // namespace viennals
