#pragma once

// Generic level-set resample engine.
//
// Produces a copy of a level set on a NEW grid spacing over the same physical
// bounds and boundary conditions (a reflective edge moves inwards by less than
// one new cell when the extent is not a multiple of the new spacing), by a
// ToSurfaceMesh -> FromSurfaceMesh round trip.
// Measured EXACT on analytic tests (plane and sphere, refine 2x and coarsen
// 3x: max error <= 0.02 nm): the triangulation nodes are the sub-cell zero
// crossings, so nothing is quantised away.  (Multilinear interpolation of the
// stored band values is not an alternative: values beyond the first band
// layer are Expand-propagated approximations, not Euclidean distances --
// ~60 nm max error on a r = 50 nm sphere.)
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
// by the caller.

#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsMesh.hpp>
#include <lsPreCompileMacros.hpp>
#include <lsToSurfaceMesh.hpp>

#include <algorithm>
#include <cmath>

namespace viennals {

using namespace viennacore;

template <class T, int D> class Resample {
  using DomainType = SmartPointer<Domain<T, D>>;

  DomainType input_ = nullptr;
  DomainType output_ = nullptr;
  double newDelta_ = 0.;
  // Width (in new-grid layers) of the produced narrow band.
  int expandWidth_ = 3;

public:
  Resample() = default;
  Resample(DomainType input, double newDelta)
      : input_(input), newDelta_(newDelta) {}

  void setLevelSet(DomainType input) { input_ = input; }
  void setNewGridDelta(double delta) { newDelta_ = delta; }
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
    auto mesh = SmartPointer<Mesh<T>>::New();
    ToSurfaceMesh<T, D>(input_, mesh).apply();
    FromSurfaceMesh<T, D>(output_, mesh).apply();
    Expand<T, D>(output_, 2 * expandWidth_ + 1).apply();
  }

private:
  // The new grid covers the old one's extent.  A reflective edge must lie
  // INSIDE that extent: Domain rounds bounds outwards, and when the extent
  // was not a multiple of the new spacing the edge column fell beyond the
  // old grid, where the source has no surface.  That column was filled with
  // the undefined sentinel (+-1.8e308) as defined values, the Si normal there
  // became NaN and the oxide diffusion matrix non-finite: a LOCOS run whose
  // regrid left such a column aborted.  The edge now moves inwards by less
  // than one new cell.  (A periodic extent must be a multiple of the new
  // spacing; it is left as it was.)
  DomainType makeEmptyTarget() const {
    const auto &grid = input_->getGrid();
    const double oldDelta = grid.getGridDelta();
    double bounds[2 * D];
    BoundaryConditionEnum bcs[D];
    for (int i = 0; i < D; ++i) {
      bcs[i] = grid.getBoundaryConditions(i);
      double low = grid.getMinBounds(i) * oldDelta;
      double high = grid.getMaxBounds(i) * oldDelta;
      const bool lowReflective =
          bcs[i] == BoundaryConditionEnum::REFLECTIVE_BOUNDARY ||
          bcs[i] == BoundaryConditionEnum::POS_INFINITE_BOUNDARY;
      const bool highReflective =
          bcs[i] == BoundaryConditionEnum::REFLECTIVE_BOUNDARY ||
          bcs[i] == BoundaryConditionEnum::NEG_INFINITE_BOUNDARY;
      const double lowIndex = std::ceil(low / newDelta_ - 1e-9);
      const double highIndex = std::floor(high / newDelta_ + 1e-9);
      const bool moved =
          (lowReflective &&
           std::abs(lowIndex * newDelta_ - low) > 1e-9 * newDelta_) ||
          (highReflective &&
           std::abs(highIndex * newDelta_ - high) > 1e-9 * newDelta_);
      if (moved)
        Logger::getInstance()
            .addWarning("Resample: the extent in direction " +
                        std::to_string(i) +
                        " is not a multiple of the new grid spacing; the "
                        "reflective boundaries move from " +
                        std::to_string(low) + " / " + std::to_string(high) +
                        " to " + std::to_string(lowIndex * newDelta_) + " / " +
                        std::to_string(highIndex * newDelta_) + ".")
            .print();
      // Domain floors the lower bound and ceils the upper one
      if (lowReflective)
        low = (lowIndex + 0.25) * newDelta_;
      if (highReflective)
        high = (highIndex - 0.25) * newDelta_;
      bounds[2 * i] = low;
      bounds[2 * i + 1] = high;
    }
    return Domain<T, D>::New(bounds, bcs, newDelta_);
  }
};

PRECOMPILE_PRECISION_DIMENSION(Resample)

} // namespace viennals
