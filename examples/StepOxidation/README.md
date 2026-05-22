# StepOxidation

Demonstrates the coupled oxidant diffusion and oxide deformation workflow on a
two-dimensional silicon step geometry, using `lsOxidationModel.hpp`. This is
the simplest configuration: two level sets, no mask, pressure-coupled
diffusion and deformation, one advection per interface.

## Physical Model

### Oxidant Diffusion — Deal-Grove

Oxidant (O₂ or H₂O) diffuses through the growing SiO₂ film and reacts at the
Si/SiO₂ interface. Under the steady-state assumption (diffusion time scale much
shorter than oxide growth time scale), the oxidant concentration `C` inside the
oxide satisfies:

```
∇ · (D ∇C) = 0
```

where `D` is `diffusionCoefficient`.

**Reaction boundary condition** at the Si/SiO₂ interface: oxidant is consumed
by a first-order surface reaction,

```
-D ∂C/∂n = k_eff · C
```

where `k_eff` is the effective reaction rate (see stress coupling below).

**Gas-transfer boundary condition** at the SiO₂/ambient interface: oxidant
enters the oxide from the ambient gas at a rate proportional to the deficit
from the equilibrium concentration `C*`,

```
-D ∂C/∂n = h · (C* - C)
```

where `h` is `transferCoefficient` and `C*` is `equilibriumConcentration`.
With `h → ∞` this becomes a Dirichlet condition `C = C*` at the free surface.

### Reaction and Expansion Velocities

The Si/SiO₂ interface moves inward (silicon is consumed) at speed

```
v_Si = velocitySign · k_eff · C / (N · γ)
```

where `N` is `oxidantMoleculeDensity` and `γ` is `expansionCoefficient`. With
`velocitySign = -1`, positive speed moves the interface into silicon.

Volume expansion produces an outward displacement at the free surface. For an
ideally flat oxide, the Si fraction of total interface displacement is `1/γ`
and the ambient fraction is `(γ-1)/γ`. For `γ = 2.27` (SiO₂):

```
Si fraction:          1/γ     ≈ 0.441
Oxide/ambient fraction: (γ-1)/γ ≈ 0.559
```

The local expansion velocity fed to the deformation solver at each Si/SiO₂
boundary crossing is

```
v_exp = ((γ - 1) / γ) · k_eff · C / N
```

directed along the outward reaction-interface normal.

### Stress-Coupled Reaction Rate

Compressive pressure in the oxide lowers the reaction rate through an
Arrhenius-like factor (Sutardja and Oldham, 1988):

```
k_eff = k · clamp(exp(-α · (p - p_ref)), f_min, f_max)
```

where:
- `k`     = `reactionRate` (zero-stress rate constant)
- `α`     = `stressCouplingCoefficient` (1/Pa; positive means compression slows growth)
- `p_ref` = `referencePressure`
- `f_min` / `f_max` = `minStressRateFactor` / `maxStressRateFactor`

Pressure `p` at each grid node is supplied by the deformation solver and fed
back between iterations by `OxidationCoupledModel`.

### Stress-Coupled Diffusion Coefficient

The same Arrhenius mechanism applies to the diffusivity (Massoud and Plummer, 1987):

```
D_eff = D · clamp(exp(-β · (p - p_ref)), f_min_D, f_max_D)
```

where `β = diffusionStressCouplingCoefficient`. When `β = 0` (default) the
diffusivity is constant. A non-zero `β` models the effect of compressive oxide
stress reducing the free volume available for oxidant diffusion.

### Crystal-Orientation Reaction Rate

Silicon oxidizes at different rates on (100), (110), and (111) faces because the
surface atom density and bond configuration differ (Irene, 1978). The orientation
factor is expressed in terms of the Si surface normal `n̂` and the wafer
crystallographic axis `ê₁₀₀`:

```
k(n̂) = k_eff · [1 + (r₁₁₁ − 1) · (1 − (n̂ · ê₁₀₀)²)]
```

where `r₁₁₁ = reactionRateRatio111` is the (111)/(100) rate ratio (≈ 1.7 for
dry oxidation at 1000 °C). The unit Si normal `n̂` is computed from the
gradient of φ_Si at each oxide-side boundary node. With `r₁₁₁ = 1` (default)
the factor is identically 1 for all faces (isotropic). The parameter
`crystalAxis` sets `ê₁₀₀`; the default `{0, 1, 0}` corresponds to a (100)
wafer with the surface normal along the simulation y-axis.

### Oxide Deformation

The growing oxide is treated as a viscous material in the quasi-static limit.
The solve has three nested stages.

#### Stage 1 — Harmonic Extension (Predictor Velocity)

The boundary velocity at every Si/SiO₂ crossing is set to `v_exp` along the
local reaction-interface normal. This Dirichlet condition is then harmonically
extended through the oxide band by iteratively averaging each interior node
over its Cartesian neighbors until convergence (`harmonicIterations`,
`tolerance`). The resulting field serves as the predictor velocity. At mask
crossings (if a mask is present), a scaled no-slip condition
`maskVelocityScale · v_interior` is applied (typically `maskVelocityScale = 0`).

#### Stage 2 — Pressure Solve

From the current velocity field the divergence `div(v)` is computed at each
node with central-difference stencils using sub-grid distances. The pressure
Poisson equation

```
∇²p = -K · div(v)
```

is solved with `K = bulkModulus`. Boundary values:

- **Free surface (SiO₂/ambient):** Dirichlet from the approximate traction-free
  condition. With `freeSurfaceTractionScale = 1`:
  ```
  p_surface ≈ n · s_dev · n
  ```
  where `s_dev` is the deviatoric stress from the previous mechanics
  iteration. On the first iteration `s_dev = 0` and `p_surface = 0`.

- **Reaction interface (Si/SiO₂):** Optionally includes an elastic substrate
  support:
  ```
  p_reaction = substrateNormalStiffness · Δt · (v · n)
  ```
  This models the Si resisting normal displacement at rate `v · n`. Setting
  `substrateNormalStiffness = 0` gives a zero-flux pressure Neumann boundary.

The pressure Poisson equation is solved by point-Jacobi iteration
(`pressureIterations`, `pressureTolerance`).

#### Stage 3 — Stokes Velocity Update

The quasi-static Stokes momentum equation is

```
η · ∇²v = pressureGradientScale · (∇p - ∇ · s_dev)
```

where `η = viscosity`. The right-hand side forcing `∇p - ∇ · s_dev` is
computed with central-difference stencils at each node. Boundary values:

- **Reaction interface:** Dirichlet `v = v_exp` (expansion velocity from Stage 1).
- **Free surface:** Ghost velocity from the traction-free condition:
  ```
  v_ghost = 2 · v_surface_analytical - v_node
  ```
  giving a second-order one-sided estimate of the traction gradient.
- **Mask contact:** No-slip or scaled: `maskVelocityScale · v_interior`.

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
s_new = exp(-Δt/τ) · s_old + (1 - exp(-Δt/τ)) · 2η · dev(D)
```

where `Δt = stressTimeStep` and `τ = viscosity / shearModulus` (or
`stressRelaxationTime` if specified). Without a shear modulus, `τ = ∞` and
`s_dev = 0` at all times (purely viscous). The full Cauchy stress tensor is

```
σ = -p I + s_dev
```

### Pressure-Concentration Coupling Loop (`OxidationCoupledModel`)

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

## Numerical Implementation

### Grid and Node Selection

Both velocity fields operate on the same Cartesian grid clipped to
`setSolveBounds`. A node at grid index `(i, j)` is included in the oxide solve
when it satisfies both:

```
reactionSign · φ_Si(i,j)  >= 0   (inside the Si/SiO₂ surface)
ambientSign  · φ_amb(i,j) >= 0   (inside the SiO₂/ambient surface)
```

With defaults (`reactionSign = +1`, `ambientSign = -1`), this selects nodes
above the Si interface and below the ambient interface — the oxide band.

All selected nodes are stored in a flat `std::vector<Node>` and addressed by a
hash map from Cartesian index to node position. The hash uses an FNV-style
mixing on each coordinate.

### Embedded Boundary Stencils

When a Cartesian edge from a node to its neighbor exits the oxide, the
level-set zero crossing along that edge gives a sub-grid boundary point at
distance

```
d = h · |φ_inside| / (|φ_inside| + |φ_outside|)
```

clamped below by `minBoundaryDistance * h` to avoid singular stencils.
When multiple interfaces are crossed (e.g., a reaction crossing and an ambient
crossing on the same edge), the nearest one wins.

### Diffusion Stencil Assembly

For a node with neighbors at sub-grid distances `d₋` and `d₊` along one axis,
the second-derivative coefficients are:

```
α₊ = 2D / (d₊ · (d₋ + d₊))
α₋ = 2D / (d₋ · (d₋ + d₊))
```

Each boundary type returns a `StencilSide(distance, nodeCoefficient, constant)`,
which encodes the boundary value in the form `C_boundary = nodeCoefficient · C₀ + constant`:

| Boundary | `nodeCoefficient` | `constant` |
|---|---|---|
| Interior neighbor | 0 | `C_neighbor` (treated as known from last iteration) |
| Reaction (dist `d`) | `g/(g + k_eff)`, `g = D/d` | 0 |
| Ambient (dist `d`) | `g/(g + h)`, `g = D/d` | `h · C* / (g + h)` |
| Mask (dist `d`, `h_m > 0`) | `g/(g + h_m)`, `g = D/d` | `h_m · C_m / (g + h_m)` |
| Mask (`h_m = 0`) | 1 | 0 (zero-flux Neumann) |
| Out-of-bounds | 1 | 0 (zero-flux Neumann) |

The stencil contributions update the system as:

```
rhs   += α · constant
diag  += α · (1 - nodeCoefficient)
```

The updated concentration is `C_new = rhs / diag`, then relaxed:
`C = relax · C_new + (1-relax) · C_old`.

Convergence: `max|C_new - C_old| < tolerance`.

### Free-Surface Velocity (Local Projection)

`getScalarVelocity` on the deformation field returns the local free-surface
normal speed. For a query point on the ambient interface, the code searches
inward along the interface normal until it finds the nearest Si/SiO₂ crossing,
then evaluates the expansion velocity `v_exp` at that crossing. This avoids
a global average and gives a spatially varying speed that correctly reflects
the local geometry.

### Velocity Fields Returned to the Advector

`OxidationDiffusionVelocityField` returns `getScalarVelocity`:
```
v_Si = velocitySign · k_eff · C / (N · γ)
```

`OxidationDeformationVelocityField` returns:
- `getScalarVelocity = 0` (unused; the scalar channel is reserved for mask-contact points)
- `getVectorVelocity = V(x)` (full solved Stokes vector field)

The ambient interface is therefore advected by the full Stokes vector field.
This captures lateral flow and shear effects that a scalar local-projection
kinematic split would miss — in particular, the surface at step corners
experiences a net horizontal displacement component from the Stokes flow that
would be zero under a purely normal-speed advection.

## Geometry

The silicon step is constructed from a flat plane at `y = 0` (left side) with
a box `x ∈ [0, 5] um, y ∈ [0, 1] um` added by a boolean UNION. The initial
oxide/ambient interface is generated by geometrically advancing the Si/SiO₂
interface by `oxideThickness = 0.2 um` using a spherical distribution (isotropic
offset), which correctly handles the step corners.

```
Domain:   x ∈ [-5, 5] um   REFLECTIVE boundaries
          y ∈ [-2, 4] um   INFINITE boundaries
Grid:     gridDelta = 0.05 um
Bounds:   {-100, -40} to {100, 80} grid indices
```

## Parameters

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
| `minStressRateFactor` | 0.25 | Floor on k rate modulation |
| `maxStressRateFactor` | 4 | Cap on k rate modulation |
| `diffusionStressCouplingCoefficient` | 0 | Stress-dependent D (0 = off) |
| `reactionRateRatio111` | 1 | (111)/(100) rate ratio (1 = isotropic) |
| `crystalAxis` | {0,1,0} | (100) wafer normal direction |
| `advectionTime` | 0.1 hr | |
| `viscosity` | 1×10⁷ Pa·hr | Oxide viscosity |
| `bulkModulus` | 7.5×10⁸ Pa | Pressure ← divergence coupling |
| `shearModulus` | 3×10¹⁰ Pa | Maxwell deviatoric relaxation |
| `freeSurfaceTractionScale` | 1 | Traction-free pressure at ambient |
| `substrateNormalStiffness` | 1×10⁹ Pa/μm | Elastic Si substrate |
| `pressureGradientScale` | 0.001 | Scales pressure gradient in Stokes RHS |
| `mechanicsIterations` | 2 | Pressure/velocity outer iterations |
| `pressureIterations` | 500 | Inner pressure Jacobi iterations |
| `stokesIterations` | 100 | Inner Stokes Jacobi iterations |

The Deal-Grove rate constants relate to the above as:
```
B   = 2D · C*/N = 0.314 μm²/hr   (parabolic rate constant)
B/A = k · C*/N  = 0.74  μm/hr    (linear rate constant)
```

## Running

```bash
cmake --build build --target StepOxidation
./build/examples/StepOxidation/StepOxidation
```

## Output Files

| File | Contents |
|---|---|
| `step_oxidation_si_initial.vtk` | Si/SiO₂ interface before oxidation |
| `step_oxidation_ambient_initial.vtk` | SiO₂/ambient interface before oxidation |
| `step_oxidation_si_after.vtk` | Si/SiO₂ interface after 0.1 hr |
| `step_oxidation_ambient_after.vtk` | SiO₂/ambient interface after 0.1 hr |
| `step_oxidation_deformation.csv` | Cartesian-grid mechanics diagnostics (x, y, velocity, pressure, strain, stress) |
