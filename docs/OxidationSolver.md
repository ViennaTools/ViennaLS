# Thermal Oxidation Solver

This document describes the physics model, numerical algorithms, and implementation
of the coupled oxidation solver in ViennaLS/ViennaPS. The solver computes thermal
silicon oxidation by solving a coupled system of a diffusion PDE, a Stokes-flow
mechanical deformation PDE, and (in LOCOS mode) a mask bending PDE, all on a shared
Cartesian grid extracted from the ViennaLS level-set narrow bands.

---

## 1. Physical Model: Deal-Grove Oxidation

Silicon oxidation follows the Deal-Grove linear-parabolic model (B. E. Deal and A. S.
Grove, *J. Appl. Phys.* **36**, 3770, 1965). Three sequential steps limit the overall
rate:

1. **Oxidant transport** from the gas phase through the existing oxide to the Si/SiO2
   interface.
2. **Surface reaction** of the oxidant with silicon at the Si/SiO2 interface.
3. **Volume expansion**: consuming 1 µm³ of Si produces 2.27 µm³ of SiO2, so both
   the Si/SiO2 interface recedes and the SiO2/gas (ambient) interface rises.

The two rate constants are:
- **B** (µm²/hr) — parabolic rate constant, proportional to the oxidant
  diffusivity D through the oxide: `B = 2 D C* / N`, where C* is the
  equilibrium oxidant concentration at the ambient surface and N is the
  oxidant molecule density in the oxide.
- **B/A** (µm/hr) — linear rate constant, equal to the surface reaction rate k_s
  scaled by the concentration ratio: `B/A = k_s C* / N`.

Both constants are thermally activated (Arrhenius):

```
B    = B₀    · exp(−E_B   / k_B T)
B/A  = BoA₀  · exp(−E_BoA / k_B T)
```

The temperature-dependent Arrhenius table is implemented in `psOxidation.hpp`
(`dealGroveRow()`). Separate regimes apply for dry oxidation above and below 950 °C
(Massoud, Plummer & Irene, *J. Electrochem. Soc.* **132**, 2685, 1985).

### Orientation Dependence

The linear rate B/A depends on the exposed crystal face. The solver applies a
continuous correction per reaction-boundary face:

```
k(n̂) = k_base · [1 + (r − 1) · (1 − (n̂ · â)²)]
```

where n̂ is the outward normal at the Si/SiO2 crossing, â is the wafer-normal
crystal axis (Y by convention), and r = `reactionRateRatio111` is the ratio
of the perpendicular-face rate to the bulk-orientation rate. For Si(100): r = 1.45
(perpendicular faces are (110)-like); for Si(111): r = 1/1.68; for Si(110): r = 1/1.45.
This captures the Bird's-Beak anisotropy observed in LOCOS without requiring explicit
crystal-plane tracking.

### Stress-Dependent Kinetics

Compressive stress in the growing oxide modifies both the reaction rate and the
diffusivity via an activation-volume coupling (Kao et al., 1987):

```
k_eff(p) = k · exp(−(p − p_ref) · V_k / (k_B T))
D_eff(p) = D · exp(−(p − p_ref) · V_D / (k_B T))
```

where p is the local hydrostatic pressure in Pa, V_k and V_D are the stress-coupling
activation volumes (m³), and p_ref is the reference pressure (zero by default). With
the default `V_k = 1.76×10⁻³⁵ m³`, compressive pressures of ~1 GPa reduce the
interface reaction rate by ~20% at 1000 °C, self-limiting the oxide growth in
confined geometries.

---

## 2. Level-Set Representation

The solver tracks two (or three) signed-distance level-set functions:

| Level set | Sign convention | Physical meaning |
|---|---|---|
| `reactionInterface` (Si) | positive = solid Si | Si/SiO2 boundary |
| `ambientInterface` (SiO2) | positive = outside oxide | SiO2/gas boundary |
| `maskInterface` (Si3N4) | positive = outside mask | Mask/oxide boundary (LOCOS only) |

A point is inside the oxide when:
- `reactionSign · φ_reaction ≥ 0` (on the Si side of the Si/SiO2 interface or
  beyond it into the oxide), **and**
- `ambientSign · φ_ambient ≥ 0` (below the ambient surface).

ViennaLS uses HRLE (hierarchical run-length encoding) to represent the narrow bands
efficiently. Both level sets are iterated simultaneously to find oxide-interior
grid points.

### Floating-Point Robustness

`GeometricAdvect` (used to seed the initial thin oxide) computes interface positions
from `k * gridDelta` arithmetic. For non-zero k this is not exactly representable in
IEEE 754 double, so the ambient surface can land at `φ = +8.87×10⁻¹⁶` instead of
exactly 0.0. Without a tolerance, `isInsideOxide` and `crosses()` both fail silently:
no oxide nodes are created at the ambient surface row, no AMBIENT boundary condition
is imposed, and the diffusion field gets no source — leaving zero oxide growth.

The fix (eps = 1×10⁻⁹) in `isInsideOxide()` and `crosses()` handles this across the
board: it is large enough to absorb the ~4×DBL_EPSILON residual from `GeometricAdvect`
and small enough to never merge two grid-cell-separated surfaces (gridDelta ≥ 0.001 µm
in practice).

---

## 3. Cartesian Solve Grid

Before each time step the solver builds a structured Cartesian grid over the oxide
region. The grid is extracted from the union of the HRLE narrow bands:

1. **Scan defined points** in both level sets; record the index bounding box with
   ±4-cell padding (`definedPointBounds` in `OxidationSolverBase`).
2. **Allocate a flat look-up array** (`nodeLookupFlat`) of size
   `(maxI - minI + 1) × (maxJ - minJ + 1)` (times `(maxK - minK + 1)` in 3D)
   initialised to `noNode`.
3. **Walk the bounding box**; for each grid index call `isInsideOxide(φ_reaction, φ_ambient)`.
   Every point that passes is registered as a solve node and its linear index is written
   into `nodeLookupFlat`. Faces between two oxide nodes are skipped; only faces that
   cross a level-set surface get classified as boundaries.

The maximum node count is guarded by `maxGridPoints` (default 5 million). The same
index range is used by both the diffusion and deformation solvers, giving them an
identical stencil neighbourhood and simplifying pressure ↔ concentration coupling.

---

## 4. Sub-Grid Boundary Geometry

Finite-difference operators at boundary faces need the distance from the node centre
to the physical interface, not just whether a crossing exists. The solver uses linear
interpolation of the signed-distance values:

```
d = gridDelta · |φ_inside| / (|φ_inside| + |φ_outside|)
```

clamped to `[minBoundaryDistance × gridDelta, gridDelta]`. The minimum distance guard
(default 1×10⁻⁶ · gridDelta) prevents division-by-zero when a surface coincides with
a grid node; the gridDelta upper cap prevents negative distances.

This sub-grid position enters every boundary coefficient:

- **AMBIENT face** (Robin BC): the transfer coefficient term `h·C*·(Δx/d)` and the
  diffusion-over-distance term `D·(Δx/d)` use d directly, so the Robin stiffness
  scales correctly with the actual oxidant path length regardless of where the surface
  falls within the grid cell.
- **REACTION face** (Neumann/Robin BC): the local reaction rate `k_eff(p)` is
  multiplied by the surface concentration C at the sub-grid crossing, and the
  resulting flux `k_eff · C / N` drives the Si surface velocity.
- **MASK face** (LOCOS): the traction BC requires the oxide velocity at the actual
  interface position; the deformation solver interpolates linearly using d.

---

## 5. Diffusion Solver (`lsOxidationDiffusion.hpp`)

The steady-state diffusion equation (oxidant transport is fast compared to interface
motion) is solved on the Cartesian oxide grid:

```
∇ · (D_eff(p) ∇C) = 0   inside oxide
```

with boundary conditions at every face that crosses a level-set surface:

| Boundary type | Condition |
|---|---|
| **AMBIENT** (SiO2/gas) | Robin (mixed): `−D ∂C/∂n = h(C − C*)` where h = transfer coefficient |
| **REACTION** (Si/SiO2) | Robin: `−D ∂C/∂n = k_eff(p) · C / N` |
| **MASK** (LOCOS) | Robin or Dirichlet: `C = C_mask` (small non-zero value) |

In the limit h → ∞ the ambient BC becomes Dirichlet C = C* (instantaneous gas-phase
replenishment), which is the classical Deal-Grove assumption. The default transfer
coefficient of 100 µm/hr approximates this well.

The linear system is assembled face-by-face as a sparse matrix in CSR format and solved
with **BiCGSTAB** with Jacobi preconditioning. An optional GPU back-end
(`lsOxidationBiCGSTABInterface.hpp`) offloads the solve to CUDA if the node count
exceeds the threshold (default 20 000 nodes) and a CUDA device is available; the
CPU path is always available as a fallback.

The Si surface velocity driven by diffusion is:

```
v_Si = −k_eff(p) · C(x_Si) / (N · ρ)
```

where C(x_Si) is the concentration at the sub-grid reaction-surface crossing and
ρ normalises by the atomic density. The expansion coefficient β = 2.27 splits this
into the Si recession rate (inward, 1/β) and the ambient growth rate (outward,
(β−1)/β).

---

## 6. Deformation Solver (`lsOxidationDeformation.hpp`)

The oxide grows as a viscous incompressible fluid driven by the volume source at the
Si/SiO2 interface. The deformation solver finds the oxide velocity field that
is consistent with both the source and the mechanics. The governing equations are the
Stokes flow problem with an incompressibility constraint modified by the expansion:

```
−∇p + η ∇²u = 0          (momentum)
∇ · u = ṡ                 (modified continuity: local volume source ṡ at reaction boundary)
```

supplemented by viscoelastic stress tracking:

```
σ = σ_elastic + σ_viscous
  = G · ε_dev + η · ε̇_dev + (K − 2G/3) · tr(ε) · I
```

where G = shear modulus (~3×10¹⁰ Pa), K = bulk modulus, and η = viscosity
(~10¹⁰ Pa·hr at 1000 °C, Arrhenius-scaled). The stress relaxes over a
characteristic time τ = η/G.

### SIMPLE Algorithm

The incompressible Stokes system is solved with the **SIMPLE** (Semi-Implicit Method
for Pressure-Linked Equations) pressure-velocity coupling:

1. **Predictor step**: solve for ũ (velocity) ignoring the pressure-gradient
   correction, using the current pressure field pⁿ.
2. **Pressure correction**: solve the Poisson equation
   `∇²p' = ∇ · ũ / (α_p Δt)` for the pressure correction p'.
3. **Velocity correction**: `u = ũ − α_v Δt ∇p'` where α_v and α_p are
   under-relaxation factors (defaults: 0.7 and 0.5).
4. **Update pressure**: `p^{n+1} = p^n + p'`.
5. **Iterate** until `||∇·u|| / ||u|| < stokesTolerance`.

The outer SIMPLE loop runs until convergence of both the velocity divergence and the
pressure field (controlled by `mechanicsIterations` and `mechanicsTolerance`). Typical
settings for stable LOCOS runs: `mechanicsIterations=200`, `mechanicsTolerance=1e-2`,
`pressureTolerance=1e-3`, `stokesTolerance=1e-3`. The tolerance hierarchy matters:
`mechanicsTolerance ≥ 5 × max(pressureTolerance, stokesTolerance)` — the mechanics
residual cannot go below the inner solver noise floor.

The ambient interface velocity is read directly from the solved velocity field u(x)
evaluated at the sub-grid ambient surface crossing. The Si surface velocity is the
scalar projection of u onto the Si surface normal.

---

## 7. Pressure–Concentration Coupling

The diffusion coefficient D_eff and reaction rate k_eff both depend on the local
oxide pressure p (from the deformation solver), while the velocity field depends on
the concentration gradient (through the volume-source boundary condition). The
coupling loop in `OxidationModel` iterates:

```
For i = 1, 2, ..., maxIterations:
    1. Solve diffusion with current {p_i} → get C_i
    2. Solve deformation with current {C_i} → get p_i+1
    3. Compute Aitken Δ² blended pressure p̃_i+1
    4. Feed p̃_i+1 back into diffusion as new pressure field
    5. Check relative change ‖p̃_i+1 − p̃_i‖ / ‖p̃_i+1‖ < tolerance
```

**Aitken Δ² acceleration** computes the mixing parameter θ ∈ [0.1, 1.0] that
minimises the squared residual norm at each coupling iteration, updating it as
`θ = −⟨F_{k-1}, ΔF⟩ / ‖ΔF‖²`. Clamping θ ≤ 1 disables extrapolation: the
exponential dependence of rate constants on pressure makes over-shooting unstable.
The coupling typically converges in 3–5 iterations for standard geometries.

---

## 8. CFL-Limited Time Stepping

Interface advection with explicit Euler or Runge-Kutta requires the CFL condition:

```
Δt ≤ α × gridDelta / max_velocity
```

where α < 0.5 (default 0.499). The time-step controller in `psOxidation.hpp` and
`lsOxidation.hpp` implements a two-pass predictor-corrector:

1. **Predictor**: solve fields at the trial Δt (first step: estimated from B/A).
2. **Compute max velocity** across all interface points.
3. **CFL check**: if `v_max × Δt / gridDelta > α`, set `Δt_actual = α × gridDelta / v_max`
   and re-solve fields at the smaller step.
4. **Accept** if the corrected velocity is consistent (i.e. re-checking with the
   newly computed v_max still satisfies CFL).

The outer loop in `psOxidation.hpp` accumulates steps until the total time reaches
the requested oxidation time, growing each step up to `maxStepGrowth = 2×` the
previous accepted step to avoid oscillations. If a solve fails to converge, the step
is halved (up to 16 retries) before throwing.

---

## 9. Interior Fill and Warm-Start Persistence

After each advection step ViennaLS's `lsAdvect` leaves only a narrow band around the
new surface position; interior data (concentration, pressure, stress) are discarded.
The solver restores them in two phases:

1. **Before advection** (`diffusionField->writePersistentFields()`): write solved
   concentration, pressure, and viscoelastic stress tensors back into the oxide
   level-set's `pointData` arrays as named scalar/vector fields (`OxConcentration`,
   `OxPressure`, `OxStressR0..R2`).
2. **After advection + interior fill** (`lsInterior`):
   a. `lsInterior` floods the interior of the advected level-set with the interior
      value (−1 for SiO2), expanding the HRLE back to the full oxide domain.
   b. `writePersistentFields()` is called again — now over the full interior HRLE —
      so the next step's `buildNodes()` finds saved values at every interior node and
      can warm-start the BiCGSTAB solve with the previous solution instead of zero.

This warm-start consistently halves the BiCGSTAB iteration count on subsequent steps.

**Important**: `lsInterior::apply()` intentionally does **not** call
`getDomain().segment()` at the end. Calling `segment()` on the freshly filled HRLE
invalidates the iterator state used by the downstream solvers and triggers a
multi-segment-domain bug in multi-threaded builds. The single-segment domain left by
`lsInterior` works correctly with all downstream iterators.

---

## 10. LOCOS Mode

When a Si3N4 mask level set is present (`maskInterface`), the solver switches to LOCOS
(Local Oxidation of Silicon) mode. LOCOS has three additional physics components:

### 10.1 Mask Bending Solver (`lsOxidationMask.hpp`)

The nitride mask deforms under the traction exerted by the growing oxide. The mask
is treated as a viscoelastic sheet; the governing equation is the 2D/3D biharmonic
(plate-bending) equation with a Neumann boundary condition at oxide-contact faces:

```
σ_mask · n̂ = σ_oxide · n̂     (traction continuity at mask/oxide interface)
```

The traction σ_oxide · n̂ is computed from the Stokes-flow stress tensor at the
mask surface and applied as a Neumann load on the mask boundary. The mask solve uses
a multigrid GMRES solver with a Gauss-Seidel (SOR) smoother.

**Contact modes** (`contactMode`):
- **0 — kinematic**: mask velocity at contact faces equals the oxide velocity
  (legacy Dirichlet BC; no traction computation).
- **1 — one-way (default)**: oxide traction drives mask; mask stiffness resists.
  Unilateral contact (compression-only) is optional.
- **2 — elastic**: two-way: elastic equilibrium of the mask, with oxide contact
  prescribing displacement v_oxide × Δt. The outer coupling loop feeds the
  resulting mask displacement back into the next oxide solve.

**Arrhenius creep viscosity** scales the mask stiffness with temperature:
```
η_mask(T) = η_ref · exp(E_a / R × (1/T − 1/T_ref))
```
with E_a ≈ 3.86×10⁵ J/mol (~4 eV, Senez et al., *IEEE Trans. Electron Devices*,
1994) and η_ref = 5×10¹¹ Pa·hr at 1000 °C.

### 10.2 Constrained Ambient Advection

Under the nitride mask the oxide is prevented from growing upward; instead it flows
laterally (the Bird's Beak). `OxidationConstrainedAmbient` enforces this by:

1. Accepting the full Stokes velocity at oxide faces away from the mask.
2. Projecting the velocity onto the mask-tangent plane at oxide faces in contact with
   the mask, eliminating the normal component that would push the oxide through the
   nitride.
3. Clipping the ambient level set to the mask exterior before advection
   (`BooleanOperation: RELATIVE_COMPLEMENT`).

### 10.3 Mask–Oxide Coupling Loop

After the first coupled diffusion+deformation solve and mask bending solve, the
solver iterates the pair until the mask velocity field stabilises:

```
For j = 1, ..., maskCouplingIterations:
    1. Re-solve coupled diffusion + deformation
    2. Re-solve mask bending with updated oxide traction
    3. Compute residual: ‖v_mask^{j} − v_mask^{j−1}‖ / ‖v_mask^{j}‖
    4. Converge if residual < maskCouplingTolerance
         OR mask displacement < tolerance × gridDelta
```

An adaptive relaxation scale (`adaptiveMaskRelaxationScale`) halves the outer
relaxation when the mask/oxide fixed-point oscillates, and reverts to the minimum
time step as a last resort before freezing the mask for one step.

---

## 11. Level-Set Advection

After all field solves converge, the three interfaces are advanced using `lsAdvect`:

- **Ambient interface**: advected with the deformation velocity field (Stokes
  velocity at the SiO2/gas surface). In LOCOS mode: `OxidationConstrainedAmbient`
  is used instead.
- **Si/SiO2 interface**: advected with the diffusion velocity field (proportional
  to the local flux at the Si surface).
- **Mask interface** (LOCOS): advected with the mask bending velocity.

The default integration scheme is first-order Engquist-Osher for spatial
derivatives, forward Euler for time (controllable via `setSpatialScheme` /
`setTemporalScheme`).

---

## 12. File Map

| File | Role |
|---|---|
| `lsOxidationSolverBase.hpp` | Cartesian grid: node lookup, bounds, `crosses()`, `isInsideOxide()` |
| `lsOxidationDiffusion.hpp` | Diffusion PDE, BiCGSTAB solve, Si surface velocity |
| `lsOxidationDeformation.hpp` | Stokes flow, SIMPLE pressure-velocity coupling, oxide velocity |
| `lsOxidationMask.hpp` | Mask bending, multigrid, contact mechanics, traction BC |
| `lsOxidationModel.hpp` | Pressure–concentration coupling loop (Aitken Δ²) |
| `lsOxidation.hpp` | Per-step orchestrator: CFL stepping, advection, interior fill |
| `lsOxidationPresets.hpp` | Named parameter presets (wet/dry at 1000 °C, SiN mask) |
| `lsOxidationBiCGSTABInterface.hpp` | CPU/GPU BiCGSTAB dispatcher |
| `psOxidation.hpp` (ViennaPS) | High-level process model: Arrhenius table, domain management |
