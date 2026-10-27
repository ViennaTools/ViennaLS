#pragma once

#include <lsInterior.hpp>
#include <lsOxidationDeformation.hpp>

#include <cmath>
#include <unordered_map>
#include <vector>

namespace viennals {

/// Nitride (Si₃N₄) mask of the LOCOS process.  The grid points inside the
/// mask join the oxide mechanics solve (OxidationDeformation) as a
/// linear-elastic body with stress memory: oxide and nitride share nodes and
/// interface faces, so velocity and traction are continuous across their
/// interface by construction.
struct OxidationMaskParameters {
  // Young's modulus (Pa) and Poisson ratio of the nitride.
  // Si₃N₄ at 1000 °C: ~250–270 GPa.
  double youngModulus = 270e9;
  double poissonRatio = 0.27;
  // Creep of the deviatoric stress (Maxwell): viscosity in Pa hr (0 =
  // linear-elastic mask) and its Eyring activation volume in m^3 (0 =
  // linear viscosity).  Si3N4: OxidationPresets::siliconNitrideMask.
  double viscosity = 0.;
  double shearActivationVolume = 0.;
};

/// Carries a mechanically transparent (barrier-only) mask along with the oxide
/// it rests on.
///
/// The mask's contact face already gets the right velocity by sampling the
/// oxide field at its own coordinates — that face lies on the oxide.  The rest
/// of the mask does not: those points sit off the oxide entirely, where the
/// field is undefined and evaluates to zero.  Advecting with the raw field
/// therefore moves the contact face and leaves the rest behind, deforming the
/// mask.
///
/// The contact velocity is extended through the mask by relaxation over the
/// mask's OWN points.  Interior fills the mask body first, so every point the
/// sweep touches is guaranteed to be defined — no value is ever read from
/// another level set outside its narrow band, which is what made isolated grid
/// columns fall back to a wrong velocity and gouge the oxide.
///
/// Each sweep is synchronous and averages over ALL known neighbours, so the
/// result does not depend on visit order and cannot shear where the contact
/// velocity varies laterally.  Direction never appears: propagation follows the
/// mask's own connectivity, so this is equally valid for a mask lying flat, one
/// standing against a sidewall, and general 3-D geometry.
template <class T, int D>
class OxidationMaskCarry final : public VelocityField<T> {
  using IndexType = viennahrle::Index<D>;
  using ConstSparseIterator =
      viennahrle::ConstSparseIterator<typename Domain<T, D>::DomainType>;

  std::unordered_map<std::size_t, Vec3D<T>> velocity_;
  // The oxide MECHANICS field.  The carry reads only its resolved Stokes
  // velocity (getResolvedVectorVelocity) — never the generic VelocityField
  // interface, whose fallback derives a speed from oxidant concentration,
  // which has nothing to do with mask motion.
  SmartPointer<OxidationDeformation<T, D>> oxideMechanics_ = nullptr;
  T gridDelta_ = 1.;

public:
  OxidationMaskCarry() = default;

  OxidationMaskCarry(SmartPointer<OxidationDeformation<T, D>> oxideMechanics,
                     SmartPointer<Domain<T, D>> maskInterface)
      : oxideMechanics_(oxideMechanics) {
    if (oxideMechanics == nullptr || maskInterface == nullptr)
      return;
    gridDelta_ = maskInterface->getGrid().getGridDelta();

    // Own Interior-filled copy: the caller's mask may be a bare narrow band at
    // this point, and the sweep needs the whole body defined.
    auto body = Domain<T, D>::New(maskInterface);
    Interior<T, D>(body).apply();

    std::vector<IndexType> index;
    std::vector<Vec3D<T>> value;
    std::vector<uint8_t> known;
    std::unordered_map<std::size_t, std::size_t> lookup;
    const T eps = std::numeric_limits<T>::epsilon();

    // Seed: a mask point that sees a non-zero RESOLVED Stokes velocity is on
    // (or within extension range of) the oxide contact.
    for (ConstSparseIterator it(body->getDomain()); !it.isFinished(); ++it) {
      if (!it.isDefined())
        continue;
      const auto idx = it.getStartIndices();
      Vec3D<T> coordinate{};
      for (unsigned d = 0; d < D; ++d)
        coordinate[d] = static_cast<T>(idx[d]) * gridDelta_;
      auto sampled =
          oxideMechanics->getResolvedVectorVelocity(coordinate);
      T magnitude = T(0);
      for (unsigned d = 0; d < D; ++d)
        magnitude += sampled[d] * sampled[d];
      const bool onOxide = magnitude > eps;
      lookup[detail::gridIndexHash<D>(idx)] = index.size();
      index.push_back(idx);
      value.push_back(onOxide ? sampled : Vec3D<T>{T(0), T(0), T(0)});
      known.push_back(onOxide ? uint8_t(1) : uint8_t(0));
    }

    // Relaxation sweeps until the whole body has a velocity.
    for (bool changed = true; changed;) {
      changed = false;
      auto nextValue = value;
      auto nextKnown = known;
      for (std::size_t i = 0; i < index.size(); ++i) {
        if (known[i])
          continue;
        Vec3D<T> sum{T(0), T(0), T(0)};
        unsigned count = 0;
        for (unsigned d = 0; d < D; ++d) {
          for (int offset = -1; offset <= 1; offset += 2) {
            auto neighbor = index[i];
            neighbor[d] += offset;
            const auto found = lookup.find(detail::gridIndexHash<D>(neighbor));
            if (found == lookup.end() || !known[found->second])
              continue;
            for (unsigned k = 0; k < D; ++k)
              sum[k] += value[found->second][k];
            ++count;
          }
        }
        if (count == 0)
          continue;
        for (unsigned k = 0; k < D; ++k)
          nextValue[i][k] = sum[k] / static_cast<T>(count);
        nextKnown[i] = 1;
        changed = true;
      }
      value.swap(nextValue);
      known.swap(nextKnown);
    }

    for (std::size_t i = 0; i < index.size(); ++i)
      if (known[i])
        velocity_[detail::gridIndexHash<D>(index[i])] = value[i];
  }

  template <class... Args> static auto New(Args &&...args) {
    return SmartPointer<OxidationMaskCarry>::New(std::forward<Args>(args)...);
  }

  Vec3D<T> getVectorVelocity(const Vec3D<T> &coordinate, int material,
                             const Vec3D<T> &normalVector,
                             unsigned long pointId) final {
    // On the contact, return EXACTLY what the oxide's own free surface gets.
    //
    // This is the kinematic tie that keeps the two surfaces together: the
    // shared surface advects on a single value.  Returning the grid-snapped
    // cached value here instead would advect the mask with a slightly
    // different velocity than the oxide at the same place; the contact then
    // slips by a fraction of a cell every substep and the mask peels away.
    if (oxideMechanics_ != nullptr) {
      const auto onOxide =
          oxideMechanics_->getResolvedVectorVelocity(coordinate);
      T magnitude = T(0);
      for (unsigned d = 0; d < D; ++d)
        magnitude += onOxide[d] * onOxide[d];
      if (magnitude > std::numeric_limits<T>::epsilon())
        return onOxide;
    }

    // Off the oxide (mask interior and outer surface): use the extended value.
    IndexType index;
    for (unsigned d = 0; d < D; ++d)
      index[d] = std::llround(coordinate[d] / gridDelta_);
    const auto found = velocity_.find(detail::gridIndexHash<D>(index));
    if (found == velocity_.end())
      return {0., 0., 0.};
    return found->second;
  }

  T getScalarVelocity(const Vec3D<T> &, int, const Vec3D<T> &,
                      unsigned long) final {
    // Carried, not grown: a normal growth term would change mask thickness.
    return 0.;
  }
};

} // namespace viennals
