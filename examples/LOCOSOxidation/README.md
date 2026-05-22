# LOCOSOxidation

Demonstrates full LOCOS (Local Oxidation of Silicon) simulation with a
Si₃N₄ nitride mask. Three level sets represent the Si/SiO₂ reaction
interface, the SiO₂/ambient free surface, and the nitride mask. The mask
blocks oxidant access on the left half of the structure, leaving the right
half as an open oxidation window. The growing oxide bends the mask at its
edge, producing the characteristic bird's beak geometry.

The simulation uses the `LOCOSOxidation<T,D>` wrapper class from
`lsLOCOSOxidation.hpp`, which orchestrates the complete per-timestep workflow
described below.

## Geometry

Three level sets are created on the same Cartesian grid:

| Level Set | Description |
|---|---|
| `siInterface` | Si/SiO₂ reaction interface. Flat plane at y = 0. |
| `ambientInterface` | SiO₂/ambient free surface. Created by geometric offset (spherical distribution) of `siInterface` by `padOxideThickness = 0.15 μm`. Represents the initial pad oxide surface. |
| `maskInterface` | Si₃N₄ nitride mask. Box from x = −4 to x = 0 (left half), y = 0.15 μm (pad oxide top) to y = 0.35 μm. The mask bottom is placed at `padOxideThickness − maskContactEpsilon` with `maskContactEpsilon = 1×10⁻⁶ μm` so the mask sits flat on the pad oxide while Cartesian stencils unambiguously resolve the mask boundary. |

The simulation domain is:
```
x ∈ [-4, 4] um   REFLECTIVE boundaries
y ∈ [-1, 2] um   INFINITE boundaries
gridDelta = 0.05 um
```

The mask covers `x ∈ [-4, 0]` and leaves `x ∈ [0, 4]` as the open oxidation
window.

## Physical Model

The LOCOS model has three coupled physics solves per time step. The first two
are the same as in StepOxidation; the third is specific to the mask.

### 1. Coupled Diffusion and Deformation Solve

`OxidationCoupledModel` drives `OxidationDiffusionVelocityField` and
`OxidationDeformationVelocityField` to convergence, passing the deformation
pressure back into the diffusion reaction rate at each iteration. See the
StepOxidation README for the full description of those two solvers.

The mask modifies both solves:

**Diffusion:** Nodes inside the mask (where `maskSign · φ_mask >= 0`) are
excluded from the oxide solve domain. Any stencil edge that would cross into
the mask is replaced by a mask boundary condition:

```
-D ∂C/∂n = h_m · (C_m - C)
```

With `maskTransferCoefficient = 0`, this becomes a zero-flux Neumann boundary
(`D grad(C) · n = 0`). The nitride is therefore a perfect oxidant block; no
concentration enters the oxide from the mask side.

The sign convention for the diffusion and deformation solves is `maskSign = -1`:
nodes where `φ_mask < 0` (inside the nitride box) are treated as being inside
the mask and excluded from the oxide solve region.

**Deformation:** Mask crossings in the mechanics stencils use `maskVelocityScale`
and `maskNormalStiffness` to define the contact boundary condition for the
growing oxide at the mask base:

- `maskVelocityScale = 0.15`: the oxide velocity at the mask contact is
  partially constrained toward zero.
- `maskNormalStiffness = 2×10⁹ Pa/μm`: elastic normal resistance at the mask
  contact surface, analogous to `substrateNormalStiffness` at the Si side.

### 2. Mask Bending Solve (`OxidationMaskBendingVelocityField`)

The nitride mask is a thin elastic plate. The oxide exerts a distributed
normal load on the mask's lower face; this bends the mask upward near its
edge. The bending solve determines a Cartesian vector velocity field inside
the mask body, which is used to advect the mask level set.

**Solve domain:** All Cartesian grid nodes at indices inside the mask level set
(where `maskSign · φ_mask >= 0`). Here `maskSign = +1`, so nodes where
`φ_mask >= 0` are interior mask nodes. The solve bounds are set to bracket the
mask geometry with a one-cell margin:

```
xMin: toIndex(-xExtent) = -80          (left domain boundary)
xMax: toIndex(maskEdge)  =   0          (mask right edge)
yMin: toIndex(padOxideThickness) - 1 =  2  (one row below mask bottom)
yMax: toIndex(padOxideThickness + maskThickness) + 1 = 8  (one row above mask top)
```

**Contact nodes:** A mask node is a contact node if it touches the oxide/mask
interface — specifically, if the axis-aligned face on the side where
the mask-normal has a downward component (outward mask normal pointing toward
the oxide) has a level-set crossing. Additionally, all nodes in the bottom-most
occupied grid row are treated as contact nodes. Contact nodes receive a
Dirichlet velocity from the oxide.

**Contact velocity:** The velocity imposed on each contact node is sampled from
the oxide deformation field one grid step below the contact index (i.e., in
the oxide just outside the mask boundary):

```
v_contact = elasticCompliance · velocityScale · v_oxide(x_contact − Δy ê_y)
```

where `v_oxide` is the solved Cartesian deformation velocity from
`OxidationDeformationVelocityField`.

**Elastic compliance factor:** A thicker or stiffer mask bends less. The
compliance is modeled as the ratio of a reference bending rigidity to the
actual bending rigidity, clamped to [0, 1]:

```
D(t) = E t³ / (12 (1 - ν²))           (bending rigidity of thickness t)
D_ref = E t_ref³ / (12 (1 - ν²))      (reference rigidity at t_ref)

elasticCompliance = clamp(D_ref / D(t), 0, 1)
                  = clamp((t_ref / t)³, 0, 1)
```

With `thickness = referenceThickness`, compliance = 1 (full contact). A mask
twice as thick has compliance = 0.125 (resists deflection strongly).

The Lamé parameters are derived from the mask elastic constants:
```
μ = E / (2(1 + ν))                     (shear modulus)
λ = E ν / ((1 + ν)(1 − 2ν))           (first Lamé parameter)
```

**Interior node update:** Interior nodes satisfy the quasi-static linear
elasticity equation:

```
μ ∇²v + (λ + μ) ∇(∇ · v) = 0
```

This is solved by Jacobi-style point relaxation, where each update has two
contributions:

1. **Laplacian average** (Gauss-Seidel-like step for `μ ∇²v`):
   ```
   v_new = (1 / count) · Σ_neighbors v_neighbor
   ```
   For a neighbor outside the mask (not a mask node), the contact velocity is
   substituted at the oxide-contact face; for an out-of-bounds or air-face
   neighbor the current node velocity is mirrored (Neumann).

2. **Grad-div correction** (contribution from `(λ + μ) ∇(∇ · v)`):
   ```
   v_new += gradDivWeight · h² / (2D) · ∇(∇ · v)
   ```
   where `gradDivWeight = (λ + μ) / max(λ + 2μ, ε)`. The central-difference
   divergence `∇ · v` is computed at each neighbor, and its gradient gives the
   correction. This term enforces volumetric compatibility and is what
   distinguishes elasticity from a pure Laplace smoothing.

The relaxed update is:
```
v = relaxation · v_new + (1 - relaxation) · v_old
```

Convergence: `max|Δv_component| < tolerance` over all interior nodes.

**Velocity field:** After convergence the solved nodal velocities are stored.
Queries outside the mask node set fall back to the nearest-node velocity. This
vector field is used directly to advect the mask level set.

### 3. Constrained Ambient Velocity Field

`OxidationConstrainedAmbientVelocityField` adapts the ambient interface
velocity depending on whether a point is under the mask or in the open window:

```
if maskSign · φ_mask(x) >= 0:     (point is inside the mask)
    return mask bending vector velocity   (zero scalar growth)
else:
    return oxide deformation velocity     (free-surface scalar growth)
```

The `isMaskContact` test is performed for every velocity query at advection
time. The mask contact sign for this field is `maskSign = -1`, so any ambient
point where `φ_mask < 0` (inside the nitride) is classified as a mask contact
and receives the mask bending velocity instead of the oxide expansion velocity.

This ensures that:
- Ambient points in the open window move by the kinematic oxide expansion split.
- Ambient points under the nitride move coherently with the mask body at the
  bending velocity, and receive zero additional scalar free-surface growth
  (the oxide cannot grow through the mask).

### Mandatory Boolean Clips

Before and after the three level-set advections, a boolean RELATIVE_COMPLEMENT
operation clips the ambient interface against the mask:

```
ambientInterface = ambientInterface \ maskInterface
```

This operation removes any part of the ambient interface that falls inside the
mask volume. It is required at two points:

**Pre-advection clip:** The ambient and mask level sets must not overlap at
the start of the advection step. The constrained ambient velocity field queries
`φ_mask` for each ambient surface point; if the two surfaces overlap the query
is ambiguous and the velocity field misbehaves.

**Post-advection clip:** The ambient interface can drift slightly into the mask
volume during the advection time step, especially near the bird's beak where
the mask edge moves. The post-advection clip corrects any such penetration
before the geometry is used in the next step.

These clips are structural — `LOCOSOxidation::apply()` always executes both —
not optional cleanup.

## LOCOSOxidation Wrapper Class

`LOCOSOxidation<T,D>` in `lsLOCOSOxidation.hpp` encapsulates the complete
per-step workflow. It is the recommended entry point for LOCOS simulation.

### API

```cpp
auto locos = ls::LOCOSOxidation<double, 2>::New(siInterface, ambientInterface,
                                                maskInterface);
locos->setOxidationParameters(oxParams);
locos->setDeformationParameters(defParams);
locos->setCouplingParameters(couplingParams);
locos->setMaskParameters(maskParams);
locos->setSolveBounds(diffMinIndex, diffMaxIndex);
locos->setMaskBendingBounds(maskMinIndex, maskMaxIndex);
locos->setSpatialScheme(ls::SpatialSchemeEnum::ENGQUIST_OSHER_1ST_ORDER);
locos->setTemporalScheme(ls::TemporalSchemeEnum::RUNGE_KUTTA_2ND_ORDER);
locos->apply(advectionTime);
```

After `apply()`, the internal velocity fields are accessible for diagnostics:

```cpp
locos->getDiffusionField()->getConcentration(pt);
locos->getDeformationField()->getPressure(pt);
locos->getMaskBendingField()->getNumberOfContactNodes();
```

### Workflow Inside `apply(advectionTime)`

```
1. Create OxidationDiffusionVelocityField   (maskSign = -1)
2. Create OxidationDeformationVelocityField (maskSign = -1)
3. OxidationCoupledModel::apply()            → coupled diffusion + deformation solve
4. Create OxidationMaskBendingVelocityField (maskSign = +1)
5. OxidationMaskBendingVelocityField::apply() → elastic mask bending solve
6. Create OxidationConstrainedAmbientVelocityField (maskSign = -1)
7. BooleanOperation(ambientInterface, maskInterface, RELATIVE_COMPLEMENT)  ← pre-clip
8. Advect ambientInterface with constrainedAmbient velocity
9. Advect siInterface with diffusion velocity
10. Advect maskInterface with mask bending velocity
11. BooleanOperation(ambientInterface, maskInterface, RELATIVE_COMPLEMENT) ← post-clip
```

### Sign Convention Summary

| Velocity Field | `maskSign` | Meaning |
|---|---|---|
| `OxidationDiffusionVelocityField` | −1 | Nodes with `φ_mask < 0` (inside nitride) excluded from oxide solve |
| `OxidationDeformationVelocityField` | −1 | Same exclusion for mechanics |
| `OxidationMaskBendingVelocityField` | +1 | Nodes with `φ_mask > 0` (outside nitride) are not mask nodes; solve is inside the nitride |
| `OxidationConstrainedAmbientVelocityField` | −1 | Ambient points with `φ_mask < 0` are under the mask and get bending velocity |

The sign inversion between the bending field (`+1`) and the others (`−1`)
arises because the bending solve domain is the mask interior, whereas the
diffusion/deformation solve domain is the oxide interior. The mask interior
and oxide interior are on opposite sides of the mask level set's zero
isosurface.

## Parameters

### Oxidation (`OxidationParameters`)

| Parameter | Value | Notes |
|---|---|---|
| `diffusionCoefficient` | 0.157 μm²/hr | Wet oxidation ⟨100⟩ Si at 1000 °C |
| `reactionRate` | 0.74 μm/hr | Deal-Grove B/A |
| `transferCoefficient` | 100 μm/hr | Large → C ≈ C* at ambient surface |
| `equilibriumConcentration` | 1 | C*/N normalized |
| `oxidantMoleculeDensity` | 1 | Normalized |
| `expansionCoefficient` | 2.27 | SiO₂/Si volume ratio |
| `velocitySign` | −1 | Si consumed |
| `stressCouplingCoefficient` | 1×10⁻¹⁵ Pa⁻¹ | Weak Arrhenius pressure correction on k |
| `diffusionStressCouplingCoefficient` | 0 | Stress-dependent D (0 = off) |
| `reactionRateRatio111` | 1 | (111)/(100) rate ratio (1 = isotropic) |
| `crystalAxis` | {0,1,0} | (100) wafer normal direction |
| `maskTransferCoefficient` | 0 | Nitride is a perfect oxidant block |
| `maskConcentration` | 0 | Zero oxidant inside nitride |

### Deformation (`OxidationDeformationParameters`)

| Parameter | Value | Notes |
|---|---|---|
| `viscosity` | 1×10⁷ Pa·hr | Oxide viscosity |
| `bulkModulus` | 7.5×10⁸ Pa | Pressure ← divergence coupling |
| `shearModulus` | 3×10¹⁰ Pa | Maxwell deviatoric relaxation |
| `stressTimeStep` | 0.35 hr | Maxwell relaxation time step |
| `freeSurfaceTractionScale` | 1 | Traction-free pressure at ambient |
| `substrateNormalStiffness` | 1×10⁹ Pa/μm | Elastic Si substrate |
| `maskVelocityScale` | 0.15 | Partial constraint at mask contact |
| `maskNormalStiffness` | 2×10⁹ Pa/μm | Elastic resistance at mask base |
| `pressureGradientScale` | 0.001 | Scales pressure gradient in Stokes RHS |
| `mechanicsIterations` | 2 | Pressure/velocity outer iterations |
| `pressureIterations` | 500 | Inner pressure Jacobi iterations |
| `stokesIterations` | 100 | Inner Stokes Jacobi iterations |

### Coupling (`OxidationCouplingParameters`)

| Parameter | Value | Notes |
|---|---|---|
| `maxIterations` | 8 | Outer diffusion/deformation coupling iterations |
| `tolerance` | 1×10⁻⁶ | Relative pressure change threshold |
| `relaxation` | 1 | No under-relaxation in pressure feedback |

### Mask Bending (`OxidationMaskParameters`)

| Parameter | Value | Notes |
|---|---|---|
| `youngModulus` | 2.5×10¹¹ Pa | LPCVD Si₃N₄ Young's modulus |
| `poissonRatio` | 0.27 | Si₃N₄ Poisson's ratio |
| `thickness` | 0.2 μm | Physical mask thickness |
| `referenceThickness` | 0.2 μm | Set equal to thickness → compliance = 1 (full contact) |
| `velocityScale` | 0.35 | Scales oxide contact velocity imposed on mask nodes |
| `pressureVelocityScale` | 0 | Pressure-normal term disabled (mask path uses contact velocity) |
| `maxVelocity` | 0.25 μm/hr | Clamp on mask bending velocity magnitude |
| `relaxation` | 0.9 | Bending solve under-relaxation |
| `tolerance` | 5×10⁻⁶ | Bending solve convergence threshold |

## Diagnostics

After `apply()` the example samples key quantities to verify correct behavior:

**Oxidation suppression check:**
```
openConcentration    ≈ 0.72   (high oxidant under open window)
maskedConcentration  ≈ 1.4e-4 (oxidant almost fully blocked under nitride)
openSiliconSpeed     ≈ 0.235 μm/hr
maskedSiliconSpeed   ≈ 4.6e-5 μm/hr
suppression ratio    ≈ 1.9e-4  (< 0.05 → mask sanity check passes)
```

**Mask bending:**
```
Mask elasticity nodes: 247   (all nodes inside the nitride)
Contact nodes:          82   (nodes at the oxide/mask interface)
Mask bottom velocity:  near (0, 0) far from edge (mask interior is rigid there)
Mask contact-node velocity at x ≈ -0.05: (0.001, 0.043) μm/hr (upward push)
```

The contact-node velocity is largest near the mask edge (x ≈ 0) because that
is where the oxide deformation velocity is strongest (the oxide is growing
rapidly in the open window and pushing against the mask corner). It is near
zero under the mask center (x ≈ −1.5) because the oxidant is blocked there
and the oxide barely moves.

**Diagnostics CSV** (`locos_oxidation_diagnostics.csv`):

| Column | Description |
|---|---|
| `x`, `y` | Node coordinates (μm) |
| `concentration` | Oxidant concentration at that node |
| `velocity_x`, `velocity_y` | Oxide deformation velocity (μm/hr) |
| `pressure` | Mechanical pressure (Pa) |
| `strain_trace` | Volumetric strain rate `div(v)` (hr⁻¹) |
| `von_mises_stress` | Von Mises equivalent stress (Pa) |

## Building and Running

```bash
cmake --build build --target LOCOSOxidation
./build/examples/LOCOSOxidation/LOCOSOxidation
```

## Output Files

| File | Contents |
|---|---|
| `locos_si_initial.vtk` | Si/SiO₂ interface before oxidation |
| `locos_ambient_initial.vtk` | Pad oxide/ambient surface before oxidation |
| `locos_mask.vtk` | Nitride mask (unchanging geometry reference) |
| `locos_si_after.vtk` | Si/SiO₂ interface after 0.35 hr |
| `locos_ambient_after.vtk` | SiO₂/ambient surface after 0.35 hr (bird's beak visible) |
| `locos_mask_after.vtk` | Nitride mask after 0.35 hr (slight upward bending at edge) |
| `locos_oxidation_diagnostics.csv` | Cartesian-grid concentration and mechanics fields |
