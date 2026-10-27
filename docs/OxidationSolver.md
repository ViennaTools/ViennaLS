# Thermal Oxidation Solver

This document describes the physics model, numerical algorithms, and implementation
of the coupled oxidation solver in ViennaLS/ViennaPS. The solver computes thermal
silicon oxidation by solving a coupled system of a diffusion PDE and a Stokes-flow
mechanical deformation PDE (which, in LOCOS mode, includes the elastic nitride mask),
both on a shared Cartesian grid extracted from the ViennaLS level-set narrow bands.

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
crystal axis (the last axis: y in 2-D, z in 3-D; ViennaPS sets it, the ViennaLS
default is the 2-D one), and r = `reactionRateRatio111` is the ratio
of the perpendicular-face rate to the bulk-orientation rate. For Si(100): r = 1.45
(perpendicular faces are (110)-like); for Si(111): r = 1/1.68; for Si(110): r = 1/1.45.
This captures the Bird's-Beak anisotropy observed in LOCOS without requiring explicit
crystal-plane tracking.

### Stress-Dependent Kinetics

Stress in the growing oxide modifies the reaction rate and the diffusivity via
activation-volume couplings (Kao et al., IEEE TED 35, 25 (1988); the `Vr`/`Vd` of
SUPREM-IV and TCAD):

```
k_eff = k · exp( (σ_nn + p_ref) · V_k / (k_B T))
D_eff = D · exp(−(p − p_ref)   · V_D / (k_B T))
```

The reaction responds to the stress NORMAL to the Si/SiO₂ interface,
σ_nn = n·σ·n with σ the oxide's Cauchy stress (tension positive) and n the Si-surface
normal: the new oxide must expand against it, so normal compression (σ_nn < 0)
retards the reaction.  The diffusivity responds to the hydrostatic pressure p
(compression positive).  The two are different stresses: at a convex Si surface the
interface is in normal compression while the oxide pressure there can be tensile.
V_k and V_D are activation volumes (m³), p_ref is the reference pressure (zero by
default).  Both laws read the oxide stress 2 and 3 cells into the oxide along the
front normal, extrapolated linearly to the front (σ_nn) or to the node (p, at the
nodes next to the front): the MAC stress of the first node layer is O(1) wrong next
to the front (on a Kao cylinder the pressure 0.2–0.4 GPa too high, σ_rr 0.1–0.2 GPa
too compressive, at any grid spacing). The σ_nn values are then smoothed along the
interface over about two cells; in a film too thin to sample (< ~4 cells) the median
over the interface node and its oxide face neighbours takes their place. Each node's
σ_nn belongs to its anchor on the front (the node moved back along its normal by its
front distance); the reaction law is evaluated at every interface point (§4) with
σ_nn interpolated from the anchors (Gaussian weights, one cell wide). Check: a
2-D Kao cylinder reproduces the exact 1-D radial model to 0.5 nm of consumed Si,
with each law alone and with both (test `OxidationKaoCylinder`).  The reaction factor is clamped to [0.1, 1]: stress only retards the
reaction.  Kao's interface stress is compressive for convex and concave Si alike, and
a tension-enhanced rate deepens any dent in a convex front (the oxide's tensile hoop
stress pulls a dent open), which made Kao pillars dent by ~13 nm in 0.1 h.  The presets
keep `V_k = 1.76×10⁻³⁵ m³` and `V_D = 0`, i.e. the coupling is effectively off.

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
- **REACTION face** (Neumann/Robin BC): the reaction rate `k_eff` at the
  crossing (its interface point, below) is
  multiplied by the surface concentration C at the sub-grid crossing, and the
  resulting flux `k_eff · C / N` drives the Si surface velocity.  The condition is
  imposed per crossing axis arm, where only the axis component of the normal flux is
  available: arm a uses `D ∂C/∂x_a = |n_a| k_eff C`, with n the Si-surface normal.
  Imposing the full `k_eff` on every crossing arm would sum to a normal sink
  `k_eff (|n_x| + |n_y| + |n_z|)` — exact on an axis-aligned front but √2 too strong
  at 45°, which made tilted fronts oxidize too slowly (Kao cylinders consumed 3–7%
  too little Si).  With the scaling, the interface concentration reconstructed on
  any arm is `C_node (D/d) / (D/d + k_eff)` for the normal distance d.
- **Interface points**: the crossings of the REACTION faces are the Si front's
  interface points, and everything the front does is evaluated there, never per
  grid node. Each carries its normal (the level-set gradients of the node and
  its neighbour, interpolated to the crossing), its reaction rate `k_eff` (the
  stress law with σ_nn at the crossing, the orientation factor with the
  crossing's normal) and its front speed `k_eff C_i / (N α)`, with C_i the
  crossing's own interface concentration. The level-set points of the front
  (§11) and the mechanics ghosts (§6.3) take the speed and normal at their
  point of the front, averaged over the crossings within two cells with
  Gaussian weights 0.75 cells wide times |n_a| (each crossing's share of the
  front area: the extra crossings at a stair step do not count twice, and
  grazing arms, whose interface concentration is the least accurate, hardly
  count). Per-node rates and speeds copied to nearby points (the earlier
  form) made the front speed depend on where the front sat in its cell; with
  the stress feedback that stair-step pattern grew into a corrugated front
  with ±1–3 GPa stress columns in a dry oxide (Al₂O₃-masked LOCOS, 1.2 nm
  grid, 2 h).
- **Deformation faces**: on the staggered mechanics grid (§6.3) a velocity
  face whose midpoint lies in Si takes a ghost value, linear through the
  crossing, where the oxide moves with the expansion velocity; the crossing's
  fractional distance is floored at `minMechanicsBoundaryDistance` (0.2). At the
  free surface the stress is extrapolated to its traction-free (purely
  tangential) value at the sub-grid crossing. Oxide/nitride faces (LOCOS) are
  interior faces of the shared mesh.

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

Every face carries the harmonic mean of its two nodes' D, which keeps the flux
conservative where D varies with the pressure (and continuous across the mask/oxide
jump). With each node's own D on its faces the flux leaving a node differed from the
flux entering its neighbour: a Kao cylinder captured 73 % of the tension enhancement
of the diffusivity, and a compressed LOCOS bird's beak let too much oxidant through
(Senez W = 0.6 µm, 0.25 h: H3/T_ox 0.37 → 0.28 with the front sampling of §1, which
alone gives 0.36).

The linear system is assembled face-by-face as a sparse matrix in CSR format. On the CPU
(the default back-end) it is solved with **BiCGSTAB preconditioned by smoothed-
aggregation AMG** with an ILU0 smoother (AMGCL, as the 3-D mechanics), to a residual of
`tolerance` x b_norm in the max norm, where b_norm is a typical row (the median diagonal)
times the concentration scale (max |b_i / A_ii|). It was max|b|, which the rows of nodes
on the ambient surface set: their diagonal grows as 1/distance (2e10 for a node exactly
on it), and a dry Kao pillar whose surface passed through grid points stopped at its
initial guess (no oxidant in the oxide). The former Jacobi-preconditioned BiCGSTAB broke
down on LOCOS bird's beaks (the pressure-dependent D spans five decades there and the
cut-cell rows make the matrix strongly non-symmetric) after thousands of iterations;
AMG takes 3–4 on those systems. An optional GPU back-end
(`lsOxidationBiCGSTABInterface.hpp`, BiCGSTAB with a Jacobi or ILU0 preconditioner)
offloads the solve to CUDA when `GpuMode::Gpu` is selected via `setGpuMode()` and a CUDA
device is available.

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

The oxide grows as a creeping Maxwell viscoelastic fluid driven by the Si→SiO2
volume expansion, which enters as a velocity boundary condition at the Si/SiO2
interface (expansion speed (β−1)/β · k_eff C/N along the interface normal). The
deformation solver finds the oxide velocity field consistent with that boundary
motion, the traction-free oxide/gas surface and, in LOCOS mode, the nitride mask
(§10). The governing equations are quasi-static incompressible Stokes flow:

```
−∇p + ∇ · s = 0   (momentum)
∇ · v = 0         (continuity)
```

with the Cauchy stress σ = −p·I + s, compression positive in p, and s the
Maxwell deviatoric stress of §6.2: G = shear modulus (~3×10¹⁰ Pa), η = viscosity
(~10¹⁰ Pa·hr at 1000 °C, Arrhenius-scaled), relaxation time τ = η/G (or
`stressRelaxationTime`). There is no volumetric elastic law: the oxide is
incompressible.

### 6.1 Stress Tensor Storage Convention

Throughout the solver `nodes[i].stressTensor` stores the **Cauchy stress** σ defined as:

```
σ_{rc} = s_{rc} − δ_{rc} · p
```

where s is the deviatoric stress at the end of the substep and p =
`nodes[i].pressure` is the pressure. This means `stressTensor` is **not** purely
deviatoric; the diagonal entries include −p.

### 6.2 Maxwell Relaxation

Over one substep Δt the deviatoric stress follows the exact Maxwell update for a
constant strain rate,

```
s_end = d · s_start + (1 − d) · 2η · ε̇',     d = exp(−Δt/τ)
```

where ε̇' is the deviatoric strain rate and s_start the stored stress of the
material, which stays fixed for the whole substep (it does not depend on the
number of mechanics iterations). The momentum equation treats the instantaneous
part implicitly, with the modulus η_impl = (1 − d)·η (about G·Δt for Δt ≪ τ, η
for a purely viscous oxide), and the decayed memory d·s_start explicitly:

```
∇ · (2 η_impl ε̇') − ∇p = −∇ · (d · s_start)
```

After the step s_end is carried with the oxide (semi-Lagrangian, the value at
x − v·Δt) and is the next substep's s_start. Oxide whose upstream point x − v·Δt
lies in the Si formed during the step and starts stress-free; Si grid points store
no history. (Interpolating over the oxide nodes there instead handed the first
oxide layer its own old stress every step, so it never took in new material: on a
corrugated front the stress next to the front grew as 1/h, tensile over the
fast-growing columns, and the Kao law turned that into growth of the corrugation.) The normal components live at the
grid points and are stored on the ambient level set (`OxStressR0..R2`); the
shear components live on the corners of the staggered grid (§6.3), are
carried there and stored on the level sets between solves (`OxCornerShear`,
`NitCornerShear`, keyed by the corner's lower grid point); the node fields hold
their mean for output and the stress feedback of §7.

**Stress-dependent viscosity** (`shearActivationVolume` V_c > 0; the presets use
300 Å³, the common TCAD default):

```
η(τ) = η₀ · x / sinh(x),     x = τ · V_c / (2 k T),     τ = √J₂ of the deviatoric stress
```

(Sutardja & Oldham, IEEE Trans. Electron Devices 36, 2415 (1989)). Below the
critical shear stress 2kT/V_c (0.117 GPa at 1000 °C) the viscosity stays near
η₀; above it, it falls exponentially, and the relaxation time η/G with it. A
linear viscous incompressible oxide has a pressure singularity at sharp concave
corners that grows as the grid is refined; the stress-dependent viscosity
bounds it. The viscosity lives where the stress is evaluated: at the grid
points (normal stress) and on the corners of the staggered grid (shear stress).
Each solves η = η₀ f(|s_end(η)|) for its current strain rate (bisection in
log η; unique root), and the viscosity iteration repeats the velocity–pressure
solve until viscosities and flow agree. The momentum equation carries the full
∇·(2ηε̇') with the local viscosity. The viscosity factor η/η₀ is written as
`OxViscosity`, the pressure (Pa, compression positive) as
`OxMechanicsPressure`.

### 6.3 Staggered (MAC) Discretisation

The mechanics uses a marker-and-cell grid on the level-set grid points (Harlow
& Welch, *Phys. Fluids* **8**, 2182, 1965):

| Unknown | Location |
|---|---|
| pressure p | oxide grid points |
| velocity component v_d | midpoint of the grid edge from a to a + e_d (a "face") |
| shear stress σ_ce | corner (2-D) / edge (3-D) a + (e_c + e_e)/2 |

A pressure difference between two neighbouring grid points drives the face
velocity between them directly. On a collocated grid (every unknown at the grid
points, central differences) an alternating ±p pattern has zero discrete
gradient at every point and hides from the momentum equation; a thick, stiff
nitride over a thin pad oxide drives it to ±10–30 GPa. The staggered grid has no
such mode.

- **Momentum**: one row per face, the force balance
  [σ_dd(a + e_d) − σ_dd(a)]/h + Σ_e [σ_de(upper corner) − σ_de(lower corner)]/h = 0.
- **Continuity**: one row per oxide grid point, Σ_d [v_d(+face) − v_d(−face)]/h = 0.
  The nitride has no pressure unknown.

The interfaces come from the signed-distance functions and do not follow the
grid:

- **Si/SiO₂ front.** φ_R at a face midpoint is the mean of its two grid points.
  A face with its midpoint on the oxide side is an unknown; a face in Si takes
  the ghost value v_g = (1 − 1/θ) v_p + v_front/θ, linear through the adjacent
  face v_p and the front velocity at the crossing, with θ = φ_p/(φ_p − φ_g)
  floored at `minMechanicsBoundaryDistance`. v_front is the expansion velocity
  (α − 1) × the front speed along the front normal, both taken at the crossing
  itself from the interface points (§4); where none is near, the mean of the
  face's oxide nodes' expansion velocity. The normals are level-set gradients
  from central differences of the defined values, one-sided where a neighbour
  lies outside the narrow band. (Clamping φ to one grid unit before differencing
  clipped the steeper component and turned the normals of tilted fronts
  towards the diagonal by up to 7°, at any grid spacing: the oxide around a Si
  disk flowed towards the grid diagonals, tangential velocity 7–10 % of the
  radial one at the front.) Nitride resting on bare Si is welded
  (v_front = 0).
- **Free surface (oxide/gas, nitride/gas).** Every face next to a solid grid
  point is an unknown. The normal stress at a gas point is the solid point's,
  extrapolated linearly to its traction-free value at the sub-grid crossing:
  σ_gas = (1 − 1/θ) σ_s,dd + (1/θ) (P σ_s P)_dd, with P = I − n nᵀ and n the
  surface normal from the defined level-set values. Only the tangential stress
  survives a free surface; extrapolating every component to zero strips it
  from the top cell row of each step of a sloped surface (a bending nitride
  then hinges). A corner with two or more gas points lies on the surface and
  carries the traction-free shear (P σ P)_ce of its solid points, plus 10⁻³ of
  the material modulus so that a one-cell feature of a level set stays
  attached; a gas face in its stencil takes the value of its partner face. A
  corner with one gas point (a step of the surface) keeps its full shear and
  needs no gas face.
- **Oxide/nitride.** One velocity field for both materials. The corner modulus
  is the harmonic mean over the corner's solid points (η_impl for oxide, μΔt for
  nitride), so velocity and traction are continuous across the interface.
- **Domain boundary.** A reflective level-set boundary is a mirror plane
  through the outermost grid points: the normal velocity is odd, the tangential
  velocity even, and a shear stress reflected across it changes sign.

The node velocity that advects the surfaces is the mean of the two faces of
each component; the node strain rate takes the normal rates at the point and the
mean of the surrounding corners' shear rates.

Check: a film on a tilted substrate with a uniform in-plane pre-stress is an
exact equilibrium (v = 0). An oxide film is reproduced to round-off for slopes
0–1 in 2-D. At a tilted oxide/nitride interface (a stress jump) the staircase
leaves a first-order velocity error and an O(1) pressure error in the interface
cells. The oxide around a Si disk (isotropic kinetics, linear viscosity) flows
radially: tangential velocity < 2 % of the radial one next to the front, the
radial velocity uniform to ±1.5 % over the angle, the free surface round to
1.5 nm and centred to 0.001 nm (test `OxidationIsotropy`).

### 6.4 Solution

The coupled velocity–pressure system is solved exactly:

- **2-D**: sparse LU (Eigen `SparseLU`, COLAMD ordering) of the row- and
  column-equilibrated system, with two refinement passes. Ordering and symbolic
  analysis are reused while the node set is unchanged.
- **3-D**: AMGCL FGMRES(80) to 10⁻⁹, preconditioned by a block-triangular
  Schur-complement pressure correction (one smoothed-aggregation V-cycle with
  ILU0 for the velocities, SPAI0 for the approximate Schur complement),
  warm-started from the previous iterate. The aggregation's near nullspace is
  one translation per velocity component: with AMGCL's default, the constant
  vector, an aggregate's components move together, and 3-D LOCOS needed about
  8× the iterations. A direct solve does not scale in 3-D; in 2-D the LU is
  about 2× faster than AMGCL on the stiff-nitride LOCOS stack.

A closed oxide (no gas contact and no nitride) pins the pressure to
`ambientPressure` at its first grid point.

**Viscosity iteration** (stress-dependent viscosity only): at most
`mechanicsIterations` solves. After each one every grid point and corner takes
its consistent viscosity for the current strain rate; the next viscosities are
the Anderson-accelerated (history 5) fixed point in log η, clamped to
[10⁻⁶ η₀, η₀]. The iteration stops when the largest change of a face velocity,
relative to the largest face velocity, is below `mechanicsTolerance`
(default 10⁻⁴). A linear oxide needs one solve. A re-solve on the same node
set (the next iteration of the stress coupling of §7, a CFL re-solve) starts
from the last converged viscosities.

A solve is accepted when the viscosity iteration converged and the residual
of the linear solve (componentwise backward error for the LU, relative residual
for AMGCL) is below 10⁻⁶. One INFO line per solve reports the
iterations, the continuity residual (Σ|∇·v| h^D over the oxide grid points,
relative to the oxide volume created per unit time at the Si/SiO₂ interface),
|∇·v|, the net volume rate and the pressure range.

The mechanics runs on the CPU (OpenMP assembly); `GpuMode` applies to the
diffusion solve (§5).

---

## 7. Pressure–Concentration Coupling

The diffusion coefficient D_eff depends on the oxide pressure p and the reaction rate
k_eff on the interface normal stress σ_nn (both from the deformation solver), while
the velocity field depends on the concentration gradient (through the volume-source
boundary condition). The coupling loop in `OxidationModel` iterates on the feedback
vector s = {p, and with V_k ≠ 0 the Cauchy stress σ} per node:

```
For i = 1, 2, ..., maxIterations:
    1. Solve diffusion with current {s_i} → get C_i
    2. Solve deformation with current {C_i} → get s_i+1 (one viscosity pass, below)
    3. Compute Aitken Δ² blended feedback s̃_i+1
    4. Feed s̃_i+1 back into diffusion (pressure and stress)
    5. Check relative change ‖s̃_i+1 − s̃_i‖ / max(‖s̃_i+1‖, k_B T / V) < tolerance
```

The stresses only enter through the Kao factors exp(±s V / k_B T), so a change
small against k_B T / V (V the larger activation volume) does not move the rates:
the denominator is floored there, and without stress feedback (V = 0) the loop
accepts the first iterate.

**Aitken Δ² acceleration** computes the mixing parameter θ ∈ [0.1, 1.0] that
minimises the squared residual norm at each coupling iteration, updating it as
`θ = −⟨F_{k-1}, ΔF⟩ / ‖ΔF‖²`. Clamping θ ≤ 1 disables extrapolation: the
exponential dependence of rate constants on pressure makes over-shooting unstable.
With stress kinetics on, a coupled solve takes about 10 iterations (8-24
on the Senez LOCOS stack).

**One viscosity pass per coupling iteration.** While the feedback iterates,
a deformation solve does one linear solve of the stress-dependent viscosity
iteration and then updates the viscosity for the next coupling iteration, so
the viscosities converge together with the feedback. When the coupling has
converged, one more deformation solve runs the viscosity iteration to
`mechanicsTolerance` (warm-started) and gives the velocities the interfaces
are advected with. Without stress feedback the loop ends after its first
solve, which runs the full iteration. With the warm start of §6.4, a
stress-on run on the Senez LOCOS stack (W = 0.65 µm, 0.3 h) is 3.1x faster
than restarting and converging the viscosity iteration in every coupling
iteration (7906 -> 1256 linear solves), the interfaces within 0.8 nm.
A coupling iteration also skips the walk over the ambient interface that
estimates the gas-side velocities for the CFL step and the advection's
dissipation: only the converged solve, whose velocities are advected, needs
it. It took ~30 % of a stress-on LOCOS run (SiN stack, 0.05 h: 72.5 -> 51.1 s,
results bit-identical).
In 2-D a coupling iteration's linear solve also reuses the last LU
factorisation: within a coupled solve the node set, and so the pattern, is
fixed, and a coupling iteration changes only the viscosity and the front
velocities. GMRES preconditioned by that factorisation, warm-started from the
last solution, reaches a relative residual of 10⁻¹⁰ in 3-9 iterations; a
fresh factorisation follows when it does not converge in 30. Stress on /
stress off: SiN stack 0.05 h 1.92-2.10 -> 1.59-1.72x (two pairs), Senez stack
W = 0.60 µm 0.3 h 1.25 -> 1.00x, interfaces identical.

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

1. **Before advection** (`diffusionField->writePersistentFields()`,
   `deformationField->writeFieldsToLevelSet()`): write solved concentration,
   pressure, velocity and the end-of-step viscoelastic stress back into the oxide
   level-set's `pointData` arrays as named scalar/vector fields (`OxConcentration`,
   `OxPressure`, `OxVelocity`, `OxStressR0..R2`; the nitride stress goes to the
   mask level set as `NitStressR0..R2`).
2. **After advection + interior fill** (`lsInterior`):
   a. `lsInterior` floods the interior of the advected level-set with the interior
      value (−1 for SiO2), expanding the HRLE back to the full oxide domain.
   b. `writePersistentFields()` is called again — now over the full interior HRLE —
      so the next step's `buildNodes()` finds saved values at every interior node and
      can warm-start the diffusion solve with the previous solution instead of zero.

The oxide fill starts at the gas surface, adds one layer per cycle and stops at the
Si, which is its guide (`setGuide`); with a guide it runs until no point is added.
Without one `lsInterior` stops after 100 layers (the interior of an open surface,
e.g. a half-space with an INFINITE boundary, never ends). The 100-layer stop used to
apply to the oxide fill too: in a film over ~100 cells the oxide next to the Si
front kept no stress history and no warm start.

The fill keeps the level set's width: it adds points inside only, and the band around
the surface is unchanged. It used to reset the width to 1 (the default of the domain
each cycle copies from). A filled oxide then reported one layer, and anything that
widens a level set narrower than two layers rebuilt it as a narrow band, dropping
the filled points with their stress history: `ToSurfaceMesh` does so in place, so
meshing the oxide between steps reset the stress of all but its outer layers (a dry
Kao hole, meshed every 10 steps, then consumed 24 % more Si than the 1-D model).
Operations that narrow the oxide level set between steps (`Reduce`, `Prune`) still
drop the history.

**Important**: `lsInterior::apply()` intentionally does **not** call
`getDomain().segment()` at the end. Calling `segment()` on the freshly filled HRLE
invalidates the iterator state used by the downstream solvers and triggers a
multi-segment-domain bug in multi-threaded builds. The single-segment domain left by
`lsInterior` works correctly with all downstream iterators.

---

## 10. LOCOS Mode

When a Si3N4 mask level set is present (`maskInterface`), the solver switches to LOCOS
(Local Oxidation of Silicon) mode. The mask is a diffusion barrier for the oxidant
(§5) and, by default, part of the mechanics.

### 10.1 Monolithic Oxide + Nitride Mechanics

The grid points inside the mask (above the Si surface) join the deformation solve
of §6 as nodes of a nitride body (`OxidationMaskParameters`: `youngModulus`,
`poissonRatio`, and for creep `viscosity`, `shearActivationVolume`). Oxide and
nitride share nodes and interface faces, so velocity and traction are continuous
across their interface by construction;
there is no separate mask solve and no coupling iteration.

- **Nitride**: a Maxwell body like the oxide, but compressible (Senez et al.,
  IEEE TED 43, 720 (1996)): the mean stress is elastic, the deviator s
  relaxes over τ = η_N/μ,
  σ_end = (m_start + K Δt tr(ε̇)) I + d s_start + 2 (1 − d) η_N ε̇',
  d = exp(−Δt/τ), K = λ + 2μ/3. For Δt ≪ τ this is incremental linear
  elasticity, σ_end = σ_start + Δt · (λ tr(ε̇) I + 2μ ε̇), which is also the
  nitride without a viscosity. With an activation volume the viscosity falls
  with the nitride's shear stress like the oxide's (Eyring, §6.2), iterated
  together with it (§6.4). The Si₃N₄ preset (`OxidationPresets::siliconNitrideMask`)
  takes η_N and V_N from Senez's Table II (5×10¹⁵ P and 170 Å³ at 1000 °C:
  τ ≈ 1.3 h at E 270 GPa, critical shear 2kT/V_N 0.21 GPa), interpolated in
  temperature. The step-start stress σ_start is explicit. The nitride has no
  pressure unknown. Its normal stresses are carried
  with the nitride (semi-Lagrangian) and stored on the mask level set
  (`NitStressR0..R2`), its shear stresses on the corners (§6.3).
- **Supports**: the nitride is bonded to the oxide (one velocity field) and
  welded to bare Si; at a reflective domain boundary (the symmetry plane) it has
  no normal velocity and no shear, like the oxide. There is no other clamp.
- **Solve**: one coupled system for both materials (§6.4).

The oxide surface is advected with the oxide nodes' velocity and the mask with the
nitride nodes' velocity (`OxidationMaterialVelocity`), so the two surfaces never
mix velocities across the oxide/nitride corner.

The nitride surface is tracked in a Lagrangian way (`OxidationMaskTracker`). The
nitride moves almost rigidly, by up to hundreds of nm, and level-set advection
rounds its convex corners by a numerical diffusion that grows with the distance
travelled, for any scheme order: a nitride edge lifted 0.1 um on a 5 nm grid had
its side wall 2.5 nm above the bottom corner moved in by 13 nm (Engquist-Osher
1st order) or 6.5 nm (WENO5 + Runge-Kutta 3), and a 50 nm strip lost its whole
contact with the oxide. After each advection the mask is re-sampled as
phi(x) = phi_0(x - u(x)): phi_0 is the mask at the start of the process call,
u the nitride displacement accumulated from the velocity the mask was advected
with (D-linear interpolation of phi_0). The advection still moves the narrow
band and carries the point data. phi_0 is sampled once per process call, which
costs ~1 nm at a corner over a few calls (0.9 nm for the 0.1 um lift above).

### 10.2 Barrier-Only Mask

With `setMaskBarrierOnly(true)` the mask is withheld from the mechanics: the oxide
expands as at a free surface and the mask is carried along on top
(`OxidationMaskCarry` extends the oxide velocity through the mask body). This is
exact for a full-width mask (uniform growth translates the mask) and is the mode
for flat barrier-coating validation.

### 10.3 Boolean Clips

Before and after advection the ambient level set is clipped to the mask exterior
(`BooleanOperation: RELATIVE_COMPLEMENT`), so the oxide never overlaps the mask.

---

## 11. Level-Set Advection

After all field solves converge, the interfaces are advanced using `lsAdvect`:

- **Ambient interface**: advected with the deformation velocity field (Stokes
  velocity at the SiO2/gas surface; in LOCOS mode the oxide nodes' velocity).
  A gas-side grid point, which has no node, takes the mean velocity of the
  nearest oxide nodes (all those at the smallest distance, up to 4 cells);
  the first one found in scan order shifted the free surface of a Si disk's
  oxide by 4 nm in y.
- **Si/SiO2 interface**: advected with the diffusion velocity field: each
  level-set point takes the front speed at its point of the front,
  x − φ ∇φ/|∇φ|, from the interface points (§4).
- **Mask interface** (LOCOS): advected with the nitride nodes' velocity and
  re-sampled from its shape at the start of the process call and the nitride
  displacement (§10.1), or carried with the oxide in barrier-only mode.

The default integration scheme is first-order Engquist-Osher for spatial
derivatives, forward Euler for time (controllable via `setSpatialScheme` /
`setTemporalScheme`).

---

## 12. File Map

| File | Role |
|---|---|
| `lsOxidationSolverBase.hpp` | Cartesian grid: node lookup, bounds, `crosses()`, `isInsideOxide()` |
| `lsOxidationDiffusion.hpp` | Diffusion PDE, AMG-preconditioned BiCGSTAB solve (AMGCL), Si surface velocity |
| `lsOxidationDeformation.hpp` | Maxwell Stokes flow on a staggered grid, coupled velocity–pressure solve (Eigen sparse LU in 2-D, AMGCL in 3-D), monolithic nitride, oxide/nitride velocities |
| `lsOxidationBiCGSTABInterface.hpp` | C++ (g++) interface to GPU BiCGSTAB (diffusion); safe for .cpp inclusion |
| `lsOxidationBiCGSTAB.cuh` | GPU BiCGSTAB kernels (`solveBiCGSTABInPlace`); nvcc-only |
| `lsOxidationMask.hpp` | Nitride parameters; barrier-only mask carry |
| `lsOxidationModel.hpp` | Pressure–concentration coupling loop (Aitken Δ²) |
| `lsOxidation.hpp` | Per-step orchestrator: CFL stepping, advection, interior fill |
| `lsOxidationPresets.hpp` | Named parameter presets (wet/dry at 1000 °C, SiN mask) |
| `psOxidation.hpp` (ViennaPS) | High-level process model: Arrhenius table, domain management |
