#pragma once

#include <lsDomain.hpp>

#include <hrleSparseIterator.hpp>

#include <algorithm>
#include <cmath>
#include <array>
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <unordered_map>

namespace viennals {

using namespace viennacore;

namespace detail {

template <int D>
inline std::size_t gridIndexHash(const viennahrle::Index<D> &index) {
  std::size_t seed = 0;
  for (unsigned i = 0; i < static_cast<unsigned>(D); ++i) {
    seed ^= std::hash<long long>{}(static_cast<long long>(index[i])) +
             0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}

template <class T>
inline Vec3D<T> vecScaled(const Vec3D<T> &source, T factor) {
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
template <class T>
inline T clampLevelSetPhi(T v) {
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

/// Common Cartesian-grid infrastructure shared by the three oxidation solver
/// classes (diffusion, deformation, mask bending). Provides the node lookup
/// table, grid extents, and the core grid utility methods that operate on them.
template <class T, int D>
class OxidationSolverBase {
protected:
  using IndexType = viennahrle::Index<D>;
  using ConstSparseIterator =
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>;

  std::unordered_map<std::size_t, std::size_t> nodeLookup;
  IndexType minIndex{};
  IndexType maxIndex{};
  std::array<std::size_t, D> extents{};
  std::array<std::size_t, D> strides{};
  T gridDelta = 1.;

  bool crosses(T a, T b) const {
    return (a <= 0. && b >= 0.) || (a >= 0. && b <= 0.);
  }

  T valueAt(ConstSparseIterator &it, const IndexType &index) const {
    it.goToIndices(index);
    return it.getValue();
  }

  bool inBounds(const IndexType &index) const {
    for (unsigned i = 0; i < D; ++i)
      if (index[i] < minIndex[i] || index[i] > maxIndex[i])
        return false;
    return true;
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
  // nodes even after the surface has advanced several grid cells beyond the band.
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
          const auto found = nodeLookup.find(linearIndex(candidate));
          if (found != nodeLookup.end() && distance2 < bestDistance2) {
            bestDistance2 = distance2;
            bestNode = found->second;
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

    minIndex = reactionGrid.getMinGridPoint();
    maxIndex = reactionGrid.getMaxGridPoint();
    std::size_t numGridPoints = 1;
    for (unsigned i = 0; i < D; ++i) {
      minIndex[i] = std::max(minIndex[i], ambientGrid.getMinGridPoint(i));
      maxIndex[i] = std::min(maxIndex[i], ambientGrid.getMaxGridPoint(i));
      if (maskInterface != nullptr) {
        minIndex[i] =
            std::max(minIndex[i], maskInterface->getGrid().getMinGridPoint(i));
        maxIndex[i] =
            std::min(maxIndex[i], maskInterface->getGrid().getMaxGridPoint(i));
      }
      if (useRequestedBounds) {
        minIndex[i] = std::max(minIndex[i], requestedMinIndex[i]);
        maxIndex[i] = std::min(maxIndex[i], requestedMaxIndex[i]);
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

    minIndex = grid.getMinGridPoint();
    maxIndex = grid.getMaxGridPoint();
    std::size_t numGridPoints = 1;
    for (unsigned i = 0; i < D; ++i) {
      if (useRequestedBounds) {
        minIndex[i] = std::max(minIndex[i], requestedMinIndex[i]);
        maxIndex[i] = std::min(maxIndex[i], requestedMaxIndex[i]);
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
    return true;
  }
};

} // namespace viennals
