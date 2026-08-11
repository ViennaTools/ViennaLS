#pragma once

// Generic level-set resample engine (see repo TODO.md, "unified design").
//
// Produces a copy of a level set on a NEW grid spacing over the same physical
// bounds and boundary conditions.  Two strategies:
//
//   SurfaceMesh (DEFAULT): ToSurfaceMesh -> FromSurfaceMesh round trip — the
//     dev dual-grid technique.  Measured EXACT on analytic tests (plane and
//     sphere, refine 2x and coarsen 3x: max error <= 0.02 nm): the
//     triangulation nodes are the sub-cell zero crossings, so nothing is
//     quantised away.
//
//   Interpolate (EXPERIMENTAL, do not use): direct multilinear interpolation
//     of stored band values.  Measured BROKEN for curved surfaces (sphere
//     r=50nm: ~60 nm max error): values beyond the first band layer are
//     Expand-propagated approximations, not Euclidean distances, so
//     interpolating them physically is invalid.  Kept only as a record;
//     making it work requires true redistancing.
//
// Usage modes (policy lives with the caller, not here):
//   PERSISTENT regrid: resample the WHOLE level-set stack with one Resample
//     per LS onto the SAME newDelta, then replace the simulation state.
//   TRANSIENT solver refinement: resample only the subset of level sets a
//     field solver consumes; use the copies for the solve, discard them.
//     A solver's subset must share one newDelta so the cut-cell boundary
//     distances see consistent relative sub-cell interface positions.
//
// Field data is NOT carried: caches keyed by quantised physical coordinate
// survive a resample by convention; grid-index-keyed data must be re-keyed
// by the caller (see TODO.md).

#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMesh.hpp>
#include <lsPreCompileMacros.hpp>
#include <lsToSurfaceMesh.hpp>

#include <hrleSparseIterator.hpp>

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace viennals {

using namespace viennacore;

enum class ResampleStrategyEnum { INTERPOLATE, SURFACE_MESH };

template <class T, int D> class Resample {
  using DomainType = SmartPointer<Domain<T, D>>;
  using IndexType = viennahrle::Index<D>;
  using ConstSparseIterator =
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>;

  DomainType input_ = nullptr;
  DomainType output_ = nullptr;
  double newDelta_ = 0.;
  ResampleStrategyEnum strategy_ = ResampleStrategyEnum::SURFACE_MESH;
  // Width (in new-grid layers) of the produced narrow band.
  int expandWidth_ = 3;

  // |phi| (grid units) beyond which a stored value is treated as background
  // rather than a usable distance.  Defined narrow-band values live within
  // ~half the band width; anything larger is a run/background artefact.
  static constexpr T usableBand_ = T(4);

public:
  Resample() = default;
  Resample(DomainType input, double newDelta)
      : input_(input), newDelta_(newDelta) {}

  void setLevelSet(DomainType input) { input_ = input; }
  void setNewGridDelta(double delta) { newDelta_ = delta; }
  void setStrategy(ResampleStrategyEnum s) { strategy_ = s; }
  void setExpandWidth(int layers) { expandWidth_ = std::max(2, layers); }

  DomainType getResult() const { return output_; }

  void apply() {
    if (input_ == nullptr || newDelta_ <= 0.) {
      Logger::getInstance()
          .addError("Resample: input level set or newDelta not set.")
          .print();
      return;
    }

    output_ = makeEmptyTarget();
    if (strategy_ == ResampleStrategyEnum::SURFACE_MESH) {
      auto mesh = SmartPointer<Mesh<T>>::New();
      ToSurfaceMesh<T, D>(input_, mesh).apply();
      FromSurfaceMesh<T, D>(output_, mesh).apply();
      Expand<T, D>(output_, 2 * expandWidth_ + 1).apply();
      return;
    }
    applyInterpolate();
  }

private:
  DomainType makeEmptyTarget() const {
    const auto &grid = input_->getGrid();
    const double oldDelta = grid.getGridDelta();
    double bounds[2 * D];
    BoundaryConditionEnum bcs[D];
    for (int i = 0; i < D; ++i) {
      bounds[2 * i] = grid.getMinBounds(i) * oldDelta;
      bounds[2 * i + 1] = grid.getMaxBounds(i) * oldDelta;
      bcs[i] = grid.getBoundaryConditions(i);
    }
    return Domain<T, D>::New(bounds, bcs, newDelta_);
  }

  void applyInterpolate() {
    const double oldDelta = input_->getGrid().getGridDelta();

    // Coarsening needs source values out to ~1 new cell around the surface;
    // widen a DEEP COPY of the source so the input is never mutated.
    // Always widen a deep copy of the source: interpolation stencils need
    // complete 2^D corners even diagonally off the band, for refining too.
    const int neededSourceLayers = std::max(
        8, 2 * (static_cast<int>(std::ceil(1.5 * newDelta_ / oldDelta)) + 2));
    DomainType source = Domain<T, D>::New(input_);
    Expand<T, D>(source, neededSourceLayers).apply();

    // Gather source values into a flat lookup (grid index -> phi in units of
    // oldDelta) and remember the near-surface points as candidate seeds.
    std::unordered_map<std::size_t, T> phiOld;
    std::vector<IndexType> seeds;
    for (ConstSparseIterator it(source->getDomain()); !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const auto idx = it.getStartIndices();
      const T value = it.getValue();
      if (std::abs(value) > usableBand_)
        continue;
      phiOld[hashIndex(idx)] = value;
      if (std::abs(value) <= T(1))
        seeds.push_back(idx);
    }
    if (seeds.empty()) {
      Logger::getInstance()
          .addWarning("Resample: source level set has no near-surface points.")
          .print();
      return;
    }

    // Candidate new-grid points: everything within ~1.5 new cells of any
    // near-surface source point.
    const double reach = std::max(oldDelta, newDelta_) * 1.5;
    const int span = static_cast<int>(std::ceil(reach / newDelta_));
    std::unordered_set<std::size_t> candidateKeys;
    std::vector<IndexType> candidates;
    for (const auto &seed : seeds) {
      IndexType base;
      for (int d = 0; d < D; ++d)
        base[d] = static_cast<viennahrle::IndexType>(
            std::round(seed[d] * oldDelta / newDelta_));
      IndexType offset;
      offset.fill(-span);
      while (true) {
        IndexType cand = base;
        for (int d = 0; d < D; ++d)
          cand[d] += offset[d];
        if (candidateKeys.insert(hashIndex(cand)).second)
          candidates.push_back(cand);
        int d = 0;
        for (; d < D; ++d) {
          if (offset[d] < span) {
            ++offset[d];
            break;
          }
          offset[d] = -span;
        }
        if (d == D)
          break;
      }
    }

    // Interpolate phi (physical units) at each candidate's coordinate from
    // the source lattice; keep points whose full multilinear stencil is
    // available and whose new-grid distance lies in the narrow band.
    typename Domain<T, D>::PointValueVectorType points;
    points.reserve(candidates.size());
    for (const auto &cand : candidates) {
      double coord[D];
      for (int d = 0; d < D; ++d)
        coord[d] = cand[d] * newDelta_;
      T phiPhys;
      if (!interpolatePhys(phiOld, oldDelta, coord, phiPhys))
        continue;
      const T phiNew = phiPhys / static_cast<T>(newDelta_);
      if (std::abs(phiNew) <= T(1.2))
        points.emplace_back(cand, phiNew);
    }
    if (points.empty()) {
      Logger::getInstance()
          .addWarning("Resample: interpolation produced no band points.")
          .print();
      return;
    }

    output_->insertPoints(points);
    output_->getDomain().segment();
    output_->finalize(2);
    Expand<T, D>(output_, 2 * expandWidth_ + 1).apply();
  }

  static std::size_t hashIndex(const IndexType &idx) {
    std::size_t seed = 0;
    for (int d = 0; d < D; ++d) {
      const auto v = static_cast<std::size_t>(
          static_cast<long long>(idx[d]) + (1ll << 30));
      seed ^= v + 0x9e3779b97f4a7c15ull + (seed << 6) + (seed >> 2);
    }
    return seed;
  }

  /// Multilinear interpolation of phi in PHYSICAL units at `coord` from the
  /// source lattice.  Fails (returns false) unless all 2^D stencil corners
  /// hold usable values — candidates outside the source band are dropped
  /// rather than guessed.
  bool interpolatePhys(const std::unordered_map<std::size_t, T> &phiOld,
                       double oldDelta, const double *coord, T &out) const {
    long long base[D];
    double frac[D];
    for (int d = 0; d < D; ++d) {
      const double c = coord[d] / oldDelta;
      const double f = std::floor(c);
      base[d] = static_cast<long long>(f);
      frac[d] = c - f;
    }
    T acc = T(0);
    for (int corner = 0; corner < (1 << D); ++corner) {
      IndexType idx;
      double w = 1.;
      for (int d = 0; d < D; ++d) {
        const int bit = (corner >> d) & 1;
        idx[d] = static_cast<viennahrle::IndexType>(base[d] + bit);
        w *= bit ? frac[d] : 1. - frac[d];
      }
      const auto found = phiOld.find(hashIndex(idx));
      if (found == phiOld.end())
        return false;
      acc += static_cast<T>(w) * found->second * static_cast<T>(oldDelta);
    }
    out = acc;
    return true;
  }
};

PRECOMPILE_PRECISION_DIMENSION(Resample)

} // namespace viennals
