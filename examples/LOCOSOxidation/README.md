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

LOCOS has four coupled physics solves per time step: oxidant diffusion, oxide
deformation (Stokes flow), pressure–concentration coupling, and mask bending
elasticity.

### Oxidant Diffusion — Deal-Grove

Oxidant (O₂ or H₂O) diffuses through the growing SiO₂ film and reacts at the
Si/SiO₂ interface. Under the steady-state assumption, the oxidant concentration
`C` inside the oxide satisfies:

```
∇ · (D ∇C) = 0
```

where `D` is `diffusionCoefficient`.

**Reaction boundary condition** at the Si/SiO₂ interface:

```
-D ∂C/∂n = k_eff · C
```

where `k_eff` is the stress-modulated effective reaction rate (see below).

**Gas-transfer boundary condition** at the SiO₂/ambient interface:

```
-D ∂C/∂n = h · (C* - C)
```

where `h` is `transferCoefficient` and `C*` is `equilibriumConcentration`.
With `h → ∞` this reduces to a Dirichlet condition `C = C*` at the free surface.

**Mask boundary condition** at any stencil edge that exits through the nitride:

```
-D ∂C/∂n = h_m · (C_m - C)
```

With `maskTransferCoefficient = 0` this is a zero-flux Neumann condition
(`D ∇C · n = 0`). The nitride is a perfect oxidant block; no concentration
enters the oxide from the mask side.

**Stencil assembly:** For a node with sub-grid boundary distances `d₋`, `d₊`
along one axis the second-derivative coefficients are:

```
α₊ = 2D / (d₊ · (d₋ + d₊))
α₋ = 2D / (d₋ · (d₋ + d₊))
```

Each boundary returns `(nodeCoefficient, constant)` such that
`C_boundary = nodeCoefficient · C₀ + constant`:

| Boundary | `nodeCoefficient` | `constant` |
|---|---|---|
| Interior neighbor | 0 | `C_neighbor` |
| Reaction (dist `d`, g=D/d) | `g/(g + k_eff)` | 0 |
| Ambient (dist `d`, g=D/d) | `g/(g + h)` | `h·C*/(g+h)` |
| Mask (`h_m = 0`) | 1 | 0 (zero-flux Neumann) |
| Out-of-bounds | 1 | 0 (zero-flux Neumann) |

The updated concentration is `C_new = rhs / diag`, then relaxed:
`C = relax · C_new + (1−relax) · C_old`.

### Reaction and Expansion Velocities

The Si/SiO₂ interface recedes at speed:

```
v_Si = velocitySign · k_eff · C / (N · γ)
```

where `N` is `oxidantMoleculeDensity` and `γ` is `expansionCoefficient`.
With `velocitySign = −1`, positive speed moves the interface into silicon.

Volume expansion produces an outward velocity at the free surface. For an
ideally flat oxide the Si fraction of total displacement is `1/γ` and the
ambient fraction is `(γ−1)/γ`. For `γ = 2.27`:

```
Si fraction:      1/γ     ≈ 0.441
Ambient fraction: (γ−1)/γ ≈ 0.559
```

The local expansion velocity fed to the deformation solver at each Si/SiO₂
crossing is:

```
v_exp = ((γ − 1) / γ) · k_eff · C / N
```

directed along the outward reaction-interface normal.

### Stress-Coupled Reaction Rate

Compressive pressure in the oxide lowers the reaction rate through an
Arrhenius-like factor (Sutardja and Oldham, 1988):

```
k_eff = k · clamp(exp(−α · (p − p_ref)), f_min, f_max)
```

where `k` is `reactionRate`, `α` is `stressCouplingCoefficient`, `p_ref` is
`referencePressure`, and `f_min`/`f_max` are the rate-modulation clamps.

### Stress-Coupled Diffusion Coefficient

The same Arrhenius mechanism optionally applies to the diffusivity:

```
D_eff = D · clamp(exp(−β · (p − p_ref)), f_min_D, f_max_D)
```

where `β = diffusionStressCouplingCoefficient`. With `β = 0` (default) the
diffusivity is constant.

### Crystal-Orientation Reaction Rate

Silicon oxidizes at different rates on (100), (110), and (111) faces. The
orientation factor in terms of the Si surface normal `n̂` and the wafer
crystallographic axis `ê₁₀₀` is:

```
k(n̂) = k_eff · [1 + (r₁₁₁ − 1) · (1 − (n̂ · ê₁₀₀)²)]
```

where `r₁₁₁ = reactionRateRatio111`. With `r₁₁₁ = 1` (default) the factor
is identically 1 for all orientations (isotropic). The parameter `crystalAxis`
sets `ê₁₀₀`; the default `{0, 1, 0}` corresponds to a (100) wafer.

### Oxide Deformation — Quasi-Static Stokes Flow

The growing oxide is treated as a viscous material in the quasi-static limit.
The solve has three nested stages.

#### Stage 1 — Harmonic Extension (Predictor Velocity)

At each Si/SiO₂ crossing, the expansion velocity `v_exp` is set as a Dirichlet
boundary condition directed along the local reaction-interface normal. This is
then harmonically extended through the oxide band by iteratively averaging each
interior node over its Cartesian neighbors until convergence (`harmonicIterations`,
`tolerance`). The result serves as the predictor velocity.

At mask crossings the oxide velocity uses the current mask velocity from the
oxide/mask interface solve. If no mask velocity field is attached, the fallback
is a stationary no-slip mask boundary.

#### Stage 2 — Pressure Solve

From the current velocity field, the divergence `div(v)` is computed at each
node with central-difference stencils using sub-grid boundary distances. The
pressure Poisson equation

```
∇²p = −K · div(v)
```

is solved with `K = bulkModulus`. Boundary values:

- **Free surface (SiO₂/ambient):** Dirichlet from the traction-free condition.
  With `freeSurfaceTractionScale = 1`:
  ```
  p_surface ≈ n · s_dev · n
  ```
  where `s_dev` is the deviatoric stress from the previous mechanics
  iteration (zero on the first iteration).

- **Reaction interface (Si/SiO₂):** Elastic substrate support:
  ```
  p_reaction = substrateNormalStiffness · Δt · (v · n)
  ```
  This models the silicon resisting normal displacement. Setting
  `substrateNormalStiffness = 0` gives a zero-flux Neumann boundary.

- **Mask contact:** Elastic mask support:
  ```
  p_mask = maskNormalStiffness · Δt · (v · n)
  ```
  This models the nitride mask resisting normal penetration from the oxide.
  Setting `maskNormalStiffness = 0` gives a free slip (Neumann) boundary.

The pressure Poisson equation is solved by point-Jacobi iteration
(`pressureIterations`, `pressureTolerance`).

#### Stage 3 — Stokes Velocity Update

The quasi-static Stokes momentum equation is:

```
η · ∇²v = ∇p − ∇ · s_dev
```

where `η = viscosity`. The right-hand side forcing `∇p − ∇ · s_dev` is
computed with central-difference stencils at each node. Boundary values:

- **Reaction interface:** Dirichlet `v = v_exp` (expansion velocity from Stage 1).
- **Free surface:** Ghost velocity from the traction-free condition:
  ```
  v_ghost = 2 · v_surface_analytical − v_node
  ```
  giving a second-order one-sided estimate of the traction gradient.
- **Mask contact:** Dirichlet velocity from the current mask mechanics solve.
  In the standalone deformation solver, without an attached mask velocity
  field, this becomes a stationary no-slip mask boundary.

The Stokes velocity equation is solved by point-Jacobi iteration
(`stokesIterations`, `stokesTolerance`).

#### Mechanics Outer Loop

Stages 2 and 3 are repeated for `mechanicsIterations` outer iterations until
the relative change in both pressure and velocity falls below
`mechanicsTolerance`. The final pressure and velocity fields are stored in the
node array.

### Maxwell Viscoelastic Deviatoric Stress

After each velocity update the symmetric strain-rate tensor is computed:

```
D_ij = 0.5 · (∂v_i/∂x_j + ∂v_j/∂x_i)
```

The deviatoric stress evolves with a Maxwell relaxation law:

```
s_new = exp(−Δt/τ) · s_old + (1 − exp(−Δt/τ)) · 2η · dev(D)
```

where `Δt = stressTimeStep` and `τ = viscosity / shearModulus` (or
`stressRelaxationTime` if specified). Without a shear modulus, `τ = ∞` and
`s_dev = 0` (purely viscous). The full Cauchy stress tensor is

```
σ = −p I + s_dev
```

### Pressure–Concentration Coupling Loop

`OxidationCoupledModel` iterates the two one-way solves into a coupled loop:

```
for up to couplingParams.maxIterations:
    diffusionField->apply()        // solve C using current k_eff(p)
    deformationField->apply()      // solve v, p using current C
    for each deformation node:
        relaxed_p = relax * p_new + (1-relax) * p_old
        diffusionField->setPressure(index, relaxed_p)
    residual = max|Δp| / max|p|
    if residual < couplingParams.tolerance: break
diffusionField->apply()            // one final solve at converged state
deformationField->apply()
```

On each pass, `stressCouplingCoefficient` causes the deformation pressure to
modulate the local reaction rate seen by the diffusion solver, closing the
feedback loop. With the small coupling used here, convergence is fast.

### Mask Bending — Quasi-Static Linear Elasticity

The nitride mask is treated as a viscous body (velocity formulation) with
an effective Si₃N₄ creep viscosity. This is identical in mathematical form to the oxide
Stokes solve but without a pressure equation — the mask is modeled as
incompressible-like, so only the vector displacement-rate field is solved.
The governing equation inside the mask is:

```
μ ∇²v + (λ + μ) ∇(∇ · v) = 0
```

where the Lamé viscosity parameters are derived from `maskViscosity` and
Poisson's ratio:

```
μ_v = maskViscosity / (2(1 + ν))
λ_v = maskViscosity · ν / ((1 + ν)(1 − 2ν))
```

#### Solve Domain

All Cartesian grid nodes inside the mask level set are included. For the LOCOS
wrapper, the mask follows the usual ViennaLS solid convention: negative level
set values are inside the nitride. The solve bounds are clipped to a
user-specified box that brackets the mask geometry.

#### Contact Nodes

A mask node is a contact node when at least one of its Cartesian neighbor
faces exits the mask through the oxide side. The oxide side is determined by
the outward mask normal: the face in direction `−n̂_mask` (toward the oxide)
is the contact face. Additionally, if the ambient interface is provided, only
nodes that are inside the oxide band (below the ambient surface) are classified
as contact nodes — this prevents mask nodes above the oxide from being driven.

The contact side is detected from the signed gradient of `φ_mask` using central
differences. HRLE far-field sentinel values are clamped before differencing so
that the contact classifier remains stable at the sparse-domain boundary.

#### Contact Boundary Condition — Traction from Oxide Stress

Contact nodes are not given a prescribed velocity. Instead they participate in
the same interior Jacobi solve as all other mask nodes, but their out-of-mask
neighbors (in the oxide) are replaced by ghost velocities derived from the
oxide traction.

**Traction at the contact face:**

The face outward normal (from the mask into the oxide) at a contact face in
direction `d` with offset `s` is `n̂ = s·ê_d`. The oxide full Cauchy stress
tensor at the adjacent oxide node is

```
σ_oxide = s_dev − p·I          (row-major: σ[3i+j] = σ_ij)
```

The traction acting on the mask from the oxide is:

```
t_i = Σ_j σ_oxide_ij · n̂_j
```

For a purely compressive oxide (positive p, negligible deviatoric), `t` points
into the mask (upward for a mask above the oxide), correctly representing the
oxide pushing the mask upward.

**Ghost velocity:**

For a first-order Neumann BC the ghost node just outside the contact face
satisfies:

```
(λ + 2μ) · (v_ghost_normal − v_node_normal) / h = t_n      (normal)
        μ · (v_ghost_tangential − v_node_tangential) / h = t_t   (tangential)
```

Rearranging:

```
v_ghost = v_node + h/(λ+2μ) · t_n · n̂  +  h/μ · t_tangential
```

where `t_n = t · n̂` and `t_tangential = t − t_n · n̂`. This ghost velocity is
substituted in place of the actual (out-of-mask) neighbor in the Laplacian
average for the Jacobi update.

The contact is unilateral by default. If the normal traction is tensile
(`t_n >= 0` with `n̂` pointing from mask to oxide), the ghost velocity falls back
to the current mask-node velocity and the oxide does not pull the mask.

Contact can also open. The oxide-side ghost node must be inside the oxide band,
or within `contactGapTolerance` level-set grid units of the oxide/ambient
surface. If the mask and oxide are separated by a larger gas gap, no traction is
applied on that face.

**Physical interpretation:**

- Larger oxide pressure → larger `t_n` → larger ghost offset → contact node
  accelerates in the normal direction.
- Larger `maskViscosity` → larger `μ` and `λ` → smaller ghost offset →
  contact node responds less to the same traction. Stiffer mask bends slower.
- Mask thickness enters naturally through the solve domain geometry: a thicker
  mask has more Cartesian cells in the vertical direction, so the traction-driven
  contact velocity decays more before reaching the top face. No explicit
  compliance scaling factor is needed.
- Both the pressure (`−p·I`) and deviatoric stress (`s_dev`) components of
  `σ_oxide` contribute to the traction. The deviatoric part encodes shear
  stresses at the mask corner (bird's beak region) that the purely scalar
  pressure misses.

#### Interior Node Update

Interior (non-contact) nodes are updated by Jacobi relaxation of the Lamé
equation. Each iteration has two contributions:

**Laplacian average** (resolves `μ ∇²v`):
```
v_avg = (1 / count) · Σ_neighbors v_neighbor
```
For a neighbor outside the mask on the oxide-contact face, the contact velocity
is substituted. For an out-of-bounds or air-face neighbor, the current node
velocity is mirrored (zero-flux Neumann).

The pure traction problem has a rigid translation mode, so the remote lateral
mask boundary is anchored explicitly. The example uses
`anchorMode = MIN_BOUNDARY` and `anchorDirection = 0`, representing the covered
nitride continuing away from the LOCOS window.

**Grad-div correction** (resolves `(λ + μ) ∇(∇ · v)`):
```
v_update += gradDivWeight · h² / (2D) · ∇(∇ · v)
```
where `gradDivWeight = (λ + μ) / max(λ + 2μ, ε)`. This term enforces
volumetric compatibility and distinguishes elasticity from pure Laplace
smoothing.

The relaxed update is:
```
v = relaxation · v_update + (1 − relaxation) · v_old
```

Convergence: `max|Δv_component| < tolerance` over all interior nodes.

### Free-Surface Velocity — Local Projection

`getScalarVelocity` on the deformation field returns the local free-surface
normal speed. For a query point on the ambient interface the code searches
inward along the interface normal until it finds the nearest Si/SiO₂ crossing,
then evaluates the expansion velocity `v_exp` at that crossing. This gives a
spatially varying speed that correctly reflects the local geometry without a
global average.

### Velocity Fields Returned to the Advector

`OxidationDiffusionVelocityField` returns `getScalarVelocity`:
```
v_Si = velocitySign · k_eff · C / (N · γ)
```

`OxidationDeformationVelocityField` returns:
- `getScalarVelocity = 0` (scalar channel unused; reserved for mask-contact points)
- `getVectorVelocity = V(x)` (full solved Stokes vector field)

The ambient interface is advected by the full Stokes vector field, capturing
lateral flow and shear effects. This is especially important at step corners
where a net horizontal displacement component from the Stokes flow would be
zero under a purely normal-speed advection.

`OxidationMaskBendingVelocityField` returns:
- `getVectorVelocity = V_mask(x)` (solved elastic displacement-rate field)

### Constrained Ambient Velocity Field

`OxidationConstrainedAmbientVelocityField` adapts the ambient interface
velocity depending on whether a query point is under the mask or in the open
window:

```
if x is inside the mask:
    return mask bending vector velocity   (zero scalar growth)
else:
    return oxide deformation velocity     (free-surface expansion)
```

Any ambient-interface point that overlaps the nitride volume moves with the
mask body and receives no additional free-surface scalar growth, so the oxide
cannot grow through the mask.

### Mandatory Boolean Clips

Before and after the three level-set advections, a boolean
RELATIVE_COMPLEMENT clips the ambient interface against the mask:

```
ambientInterface = ambientInterface \ maskInterface
```

**Pre-advection clip:** The ambient and mask level sets must not overlap when
the constrained velocity field queries `φ_mask`. Overlap makes the query
ambiguous and corrupts the velocity field.

**Post-advection clip:** The ambient interface can drift slightly into the mask
volume during the advection step, especially near the bird's beak where the
mask edge moves. The post-advection clip corrects any such penetration before
the geometry is used in the next step.

Both clips are structural — `LOCOSOxidation::apply()` always executes both.

## LOCOSOxidation Wrapper Class

`LOCOSOxidation<T,D>` in `lsLOCOSOxidation.hpp` encapsulates the complete
per-step workflow.

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

The oxide/mask interface fixed-point solve uses model defaults of six coupling
iterations and a mask-velocity residual tolerance of `2.e-2` unless overridden
with `setMaskCouplingIterations()` or `setMaskCouplingTolerance()`.

After `apply()`, the internal velocity fields are accessible for diagnostics:

```cpp
locos->getDiffusionField()->getConcentration(pt);
locos->getDeformationField()->getPressure(pt);
locos->getMaskBendingField()->getNumberOfContactNodes();
locos->getMaskCouplingResidual();
```

### Workflow Inside `apply(advectionTime)`

```
1. Create OxidationDiffusionVelocityField
2. Create OxidationDeformationVelocityField
3. OxidationCoupledModel::apply()            → coupled diffusion + deformation solve
4. Create OxidationMaskBendingVelocityField
5. OxidationMaskBendingVelocityField::apply() → elastic mask bending solve
6. Iterate the oxide/mask interface solve:
   - oxide mechanics uses the current mask velocity as its mask boundary
   - mask mechanics uses the updated oxide stress as contact traction
   - stop when the contact-velocity fixed-point residual is below tolerance
8. Create OxidationConstrainedAmbientVelocityField
9. BooleanOperation(ambientInterface, maskInterface, RELATIVE_COMPLEMENT)  ← pre-clip
10. Advect ambientInterface with constrainedAmbient velocity
11. Advect siInterface with diffusion velocity
12. Advect maskInterface with mask bending velocity
13. BooleanOperation(ambientInterface, maskInterface, RELATIVE_COMPLEMENT) ← post-clip
```

### Mask Sign Convention

`LOCOSOxidation<T,D>` assumes the mask level set uses the usual ViennaLS solid
convention: `φ_mask < 0` is inside the nitride. That convention is handled
inside the wrapper. User code does not need to pass a mask sign to the
diffusion, deformation, constrained ambient, or mask bending fields.

The same mask interior is used in different ways by the internal solvers:
diffusion and oxide deformation exclude nitride nodes from the oxide solve,
mask bending solves on those nitride nodes, and constrained ambient advection
moves any ambient surface under the nitride with the mask velocity.

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

The Deal-Grove rate constants relate to the above as:
```
B   = 2D · C*/N = 0.314 μm²/hr   (parabolic rate constant)
B/A = k · C*/N  = 0.74  μm/hr    (linear rate constant)
```

### Deformation (`OxidationDeformationParameters`)

| Parameter | Value | Notes |
|---|---|---|
| `viscosity` | 1×10¹⁰ Pa·hr | Effective oxide viscosity |
| `bulkModulus` | 7.5×10⁸ Pa | Pressure ← divergence coupling |
| `shearModulus` | 3×10¹⁰ Pa | Maxwell deviatoric relaxation |
| `stressTimeStep` | 0.35 hr | Maxwell relaxation time step |
| `freeSurfaceTractionScale` | 1 | Traction-free pressure at ambient surface |
| `substrateNormalStiffness` | 1×10⁹ Pa/μm | Elastic Si substrate resistance |
| `maskNormalStiffness` | 2×10⁹ Pa/μm | Elastic normal resistance at mask base |
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
| `maskViscosity` | 5×10¹¹ Pa·hr | Effective creep viscosity of the Si₃N₄ mask stack |
| `poissonRatio` | 0.27 | Si₃N₄ Poisson's ratio; sets λ/μ ratio in Lamé viscosity |
| `maxVelocity` | 0.25 μm/hr | Clamp on mask bending velocity magnitude |
| `unilateralContact` | true | Oxide can push the mask but not pull it |
| `contactGapTolerance` | 0 | Contact opens as soon as the oxide-side ghost node leaves the oxide |
| `anchorMode` | `MIN_BOUNDARY` | Remote lateral mask boundary is fixed |
| `anchorDirection` | 0 | x-direction anchoring for the 2D LOCOS window |
| `relaxation` | 0.9 | Bending solve under-relaxation |
| `tolerance` | 5×10⁻⁶ | Bending solve convergence threshold |

The Lamé viscosity parameters are derived as:
```
μ_v = maskViscosity / (2(1+ν))
λ_v = maskViscosity · ν / ((1+ν)(1−2ν))
```

These play the same role as `viscosity` does in the oxide Stokes solve. The
example uses a large effective value so the stress-driven response stays below
the `maxVelocity` clamp and the lateral bending profile remains visible.

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
Mask elasticity nodes: 324   (all nodes inside the nitride)
Contact nodes:          81   (nodes at the oxide/mask interface)
Interface iterations:   6, residual ≈ 1.4e-2
Mask bottom velocity:   upward, no downward pull from oxide contact
Mask lateral samples:   anchored far side ≈ 0, larger upward velocity near edge
```

The oxide/mask interface solve is a fixed-point mechanics solve, not a one-way
post-process: the oxide boundary uses the current mask velocity, the mask
boundary uses the resulting oxide traction, and the residual measures the
relative contact-velocity change between interface iterates. Tighter tolerances
are possible but require additional full oxide/mask solves. The contact-node
velocity is largest near the mask edge (x ≈ 0) because that is where the oxide
deformation velocity has the largest outward normal component relative to the
mask face. The remote covered side is anchored to remove the rigid Neumann mode,
so the solved field bends laterally instead of translating the whole mask.

**Open-window volume conservation:**
```
Si consumed area:        ≈ 0.203 μm²
Ambient lift area:       ≈ 0.277 μm²
Expected ambient lift:   ≈ 0.258 μm²
lift/Si ratio:           ≈ 1.36   (expected γ−1 = 1.27)
relative error:          ≈ 7.4%
```

This diagnostic samples vertical level-set crossings in the open window and
checks the integral relation `ambient lift ≈ (γ−1) · silicon consumption`. It is
a coarse grid-column check, not an exact volume integral, but it catches
nonphysical Si/oxide/ambient motion immediately.

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
