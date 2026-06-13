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

### 6.1 Stress Tensor Storage Convention

Throughout the solver — both CPU and GPU — `nodes[i].stressTensor` (CPU) and
`d_stress[k*n+i]` (GPU) store the **Cauchy stress** σ defined as:

```
σ_{rc} = τ_{rc} − δ_{rc} · p
```

where τ is the true deviatoric (traceless) stress and p = `nodes[i].pressure` is the
isotropic pressure. This means `stressTensor` is **not** purely deviatoric; the
diagonal entries include −p.

To recover the true deviatoric τ for traction computations:
```
τ_{rc} = stressTensor[r][c] + δ_{rc} · pressure    (CPU: deviatoricStressAt())
       = d_stress[(3r+c)*n+i] + (r==c ? pressure : 0)   (GPU inline)
```

This convention is important for boundary conditions: the free-surface (AMBIENT) ghost
velocity uses `p·n_comp − τ_comp·n`, where τ·n = Σ_j σ_{comp,j}·n_j + p·n_comp.

### 6.2 Maxwell Relaxation

At each SIMPLE iteration the deviatoric stress τ is updated with Maxwell relaxation:

```
τ^{new}_{rc} = e^{−Δt/τ_r} · τ^{prev}_{rc} + (1 − e^{−Δt/τ_r}) · 2η · ε̇^{dev}_{rc}
```

where τ_r = η/G is the stress relaxation time and ε̇^{dev} is the deviatoric part of
the strain rate tensor (traceless, symmetric).

The CPU stores τ history in `deviatoricStressHistory` (an `unordered_map<IndexType,
array<T,9>>`). The GPU stores it in the flat device array `d_stressHist[9*n]`
(component-major: `(3r+c)*n+i`). Both are updated **every SIMPLE iteration**, so the
history reflects the most recent velocity field, not just the beginning of each CFL step.

### 6.3 SIMPLE Algorithm

The incompressible Stokes system is solved with the **SIMPLE** (Semi-Implicit Method
for Pressure-Linked Equations) pressure-velocity coupling. The algorithm is
unconditionally stable for steady Stokes and avoids the spectral-radius > 1 problem
that makes a naive Gauss-Seidel p→v→p iteration diverge on thin geometries.

The three steps per iteration are:

```
1. Momentum predictor:  A_v · v* = vBC − (∇p^k − ∇·τ(v^k)) / η
2. Pressure update:     A_p · p^{k+1} = pBC(τ(v*)) + K · ∇·v*
3. Velocity correction: v^{k+1} = v* − ∇(p^{k+1}−p^k) / (η · diag(A_v)[i])
```

Step 3 ensures the corrected velocity is consistent with the new pressure without
re-solving the full momentum equation.

**Convergence criterion** (both CPU and GPU):
```
residual = max(|Δv|_max / |v|_max,  |Δp|_max / |p|_max)  <  mechanicsTolerance (1e-2)
```
Up to `mechanicsIterations` (200) outer iterations are allowed.

**Inner solver tolerances**: `stokesTolerance=1e-3`, `pressureTolerance=1e-3`.
The hierarchy `mechanicsTolerance ≥ 5 × max(stokesTolerance, pressureTolerance)`
must hold — the outer residual cannot go below the inner solver noise floor.

**Pressure relaxation**: the solved p^{k+1} is blended with the previous pressure
before applying the velocity correction:
```
p^{k+1}_relaxed = (1 − β) · p^k + β · p^{k+1}_solved,   β = pressureRelaxation (0.5)
```

**Velocity correction diagonal**: two distinct diagonal arrays are maintained:
- `stokesDiag` / `d_diagV` — diagonal of the Stokes operator A_v, computed from
  REACTION+MASK+interior faces only (OOB and AMBIENT contributions cancel in the
  matvec). Used as the preconditioner diagonal in the Stokes BiCGSTAB solve.
- `corrDiag` / `d_diagV_corr` — full diagonal including OOB and AMBIENT face terms,
  returned by `computeVelocityDiagonals()`. Used exclusively for the velocity
  correction `v^{k+1} = v* − ∇δp / (η · corrDiag_i)`.

Using `corrDiag` for the correction (instead of `stokesDiag`) is critical: at
ambient-adjacent nodes, `corrDiag > stokesDiag`, producing a smaller correction
magnitude that avoids SIMPLE oscillation.

### 6.4 CPU SIMPLE Loop

The CPU implementation lives in `solveMechanics()` → CPU branch
(lines ~1764–1847 of `lsOxidationDeformation.hpp`).

**Pre-loop setup** (once per `solveMechanics()` call):
```
diagV = computeVelocityDiagonals()   // corrDiag, geometry-fixed for all iterations
```

**Per-iteration** (up to `mechanicsIterations`):
```
previousVelocity = collectVelocities()     // save v^k
previousPressure = collectPressures()      // save p^k
computeDiagnostics()                       // div(v^k), vonMises — uses OLD stress
computeStressTensors()                     // Maxwell update: d_stress = τ−p·I, history updated
solveStokesVelocity()                      // v* in nodes[i].velocity; uses p^k, τ(v^k)
solvePressure()                            // p^{k+1}; ambient BC uses τ freshly from v*
applySimpleVelocityCorrection(p^k, diagV) // v^{k+1} = v* − ∇δp / (η·corrDiag)
residual = max(maxVelChange, maxPresChange)
```

**After convergence** (or max iterations):
```
computeDiagnostics()     // refresh div(v), vonMises with final velocity
computeStressTensors()   // refresh stressTensor with final velocity
```

**Key detail — ambient BC timing**: `solvePressure()` calls
`freeSurfacePressureBoundary()` which calls `currentBoundaryDeviatoricStress()`,
computing τ freshly from the current (post-Stokes) velocity v*. This means the
ambient pressure BC for the pressure RHS uses τ(v*), not τ(v^k). This distinction
matters for correctness and is reproduced by the GPU's step 6b.

**Ghost velocity at AMBIENT faces** (`freeSurfaceVelocityBoundary()`): uses
`deviatoricStressAt()` which recovers the true deviatoric τ by adding back p·δ to
the stored σ. The normalTraction formula is:
```
normalTraction[comp] = p · n[comp] − τ_comp · n̂
                     = p · n[comp] − Σ_j τ_{comp,j} · n_j
```

### 6.5 GPU SIMPLE Loop

The GPU SIMPLE loop (`solveMechanicsGpu()`, ~lines 1523–1716, compiled only under
`VIENNALS_GPU_BICGSTAB`) offloads the entire outer SIMPLE iteration to the GPU. It is
activated when `gpuSimpleIsValid(gpuSimpleBufs_)` is true — set up during
`buildNodes()` when `tryGpu && !useIlu0` (Jacobi preconditioner required; ILU0
requires per-component re-factorization that is not implemented for the SIMPLE loop).

**PCIe traffic**: only 4×8 = 32 bytes cross the bus per SIMPLE iteration (the
convergence metric). All assembly (RHS, stress, strain, ambient BC) is done entirely
on-device.

**Pre-SIMPLE setup** (CPU side, before the iteration loop):

1. Pack `nodes[]` into flat device-layout buffers (component-major `[c*n+i]`):
   - `h_vel[D*n]`, `h_pressure[n]`, `h_stress[9*n]` from `nodes[id]`
   - `h_stressHist[9*n]` from `deviatoricStressHistory` map (zero for unknown nodes)
2. Upload state: `gpuSimpleUploadState()`
3. Compute and upload **solid vBC** — the REACTION + MASK ghost-velocity contribution
   to the Stokes RHS. This is precomputed CPU-side and held fixed for all iterations
   because REACTION/MASK velocities do not depend on the current velocity or pressure.
   AMBIENT face contributions are **omitted** here and recomputed on-device each
   iteration from the current stress state.
4. Upload per-node **reaction velocities** (`h_reactionVel[D*n]`) for use in the
   strain-trace and strain-rate kernels (REACTION ghost velocity = oxidation velocity).
5. Upload geometry-fixed pressure and Stokes coefficient arrays
   (`gpuUploadSolverArrays` for both pressure and Stokes BiCGSTAB buffers).

**Per-iteration** (`gpuSimpleRunIterationImpl()` in `lsOxidationSimpleGpu.cuh`):

```
Step 1:  simpleSaveStateKernel            d_prevVel[D*n] ← d_vel; d_prevPressure[n] ← d_pressure
Step 2:  simpleComputeStrainTraceKernel   div(v^k), using OLD d_stress (matches CPU computeDiagnostics)
Step 3:  simpleComputeStressTensorsKernel Maxwell update → new d_stress = τ−p·I, d_stressHist = τ,
                                           d_strainRate, d_vonMises
Step 4:  simpleComputeAmbientBPKernel     d_ambientBP = p_amb + n^T·τ(v^k)·n (uses NEW d_stress)
Step 5:  simpleAssembleStokesRhsKernel    b_stokes[D*n] = solidVelBC + ambientVelBC − (∇p − ∇·τ)/η
Step 6:  Stokes BiCGSTAB (per component)  v* in d_vel, using stokesDiag (d_diagV)
Step 6b: simpleRecomputeAmbientBPFromVStarKernel
          Recompute d_ambientBP = p_amb + n^T·τ(v*)·n (τ from v*, matching CPU solvePressure)
          d_stress is NOT updated — still needed for div(v*) ghost vel in step 7
Step 7:  simpleComputeStrainTraceKernel   div(v*) with updated d_vel — for pressure RHS
Step 8:  simpleAssemblePressureRhsKernel  b_press = pBC(τ(v*)) + K·div(v*)
Step 9:  simpleSetPressureGuessKernel     warm-start: touchesAmbient → ambBP, else current p
         Pressure BiCGSTAB → p^{k+1} in pressBufs->d_x
Step 10: simpleApplyPressureRelaxationKernel  d_pressure = (1−β)·d_prevPressure + β·p^{k+1}
Step 11: simpleVelocityCorrectionKernel   d_vel −= ∇δp / (η · d_diagV_corr)  [corrDiag]
Step 12: simpleMaxChangeKernel + D2H      residual = max(|Δv|/|v|, |Δp|/|p|)
```

**Post-SIMPLE** (CPU side, after the loop):
- Download final state: `gpuSimpleDownloadState()` → `h_vel`, `h_pressure`, `h_stress`, `h_stressHist`
- Unpack into `nodes[]` and rebuild `deviatoricStressHistory` map
- Run `computeDiagnostics()` + `computeStressTensors()` on CPU to refresh derived
  quantities (`strainTrace`, `vonMisesStress`, etc.) for downstream LOCOS/mask solvers

**GpuSimpleBuffers device layout** (all arrays in `viennals::gpu::GpuSimpleBuffers`):

| Array | Size | Content |
|---|---|---|
| `d_vel` | D×n doubles | velocity, component-major `[c*n+i]` |
| `d_pressure` | n doubles | isotropic pressure |
| `d_stress` | 9×n doubles | Cauchy stress σ = τ−p·I, `[(3r+c)*n+i]` |
| `d_stressHist` | 9×n doubles | previous deviatoric τ (Maxwell history) |
| `d_strainRate` | 9×n doubles | strain rate tensor ε̇ |
| `d_strainTrace` | n doubles | div(v) = tr(ε̇) |
| `d_vonMises` | n doubles | von Mises stress |
| `d_prevVel` | D×n doubles | v^k saved at start of each iteration |
| `d_prevPressure` | n doubles | p^k saved at start of each iteration |
| `d_ambientBP` | n doubles | free-surface pressure BC value |
| `d_stokesRhs` | D×n doubles | assembled Stokes RHS |
| `d_pressRhs` | n doubles | assembled pressure RHS |
| `d_nb` | 2D×n uint32 | face-major neighbor node IDs (0xFFFFFFFF = boundary) |
| `d_bcType` | 2D×n uint8 | 0=NONE, 1=REACTION, 2=AMBIENT, 3=MASK |
| `d_bcDist` | 2D×n doubles | sub-grid interface crossing distances |
| `d_diagV` | D×n doubles | Stokes diagonal (excl. OOB/AMBIENT) — for BiCGSTAB |
| `d_diagV_corr` | D×n doubles | full correction diagonal (incl. OOB/AMBIENT) — for velocity correction |
| `d_ambNormal` | D×n doubles | ambient level-set gradient normal |
| `d_ambStokesCoeff` | 2D×n doubles | ±2/dSum for AMBIENT face Stokes vBC |
| `d_ambPressFaceCoeff` | n doubles | Σ 2/(dist·dSum) for AMBIENT pressure faces |
| `d_ambPressNeighCoeff` | 2D×n doubles | 2/(gδ·dSum) for ambient-touching pressure neighbors |
| `d_solidVelBC` | D×n doubles | REACTION+MASK Stokes RHS contribution (constant per step) |
| `d_reactionVel` | D×n doubles | oxidation velocity at each node for strain kernels |
| `d_touchesAmbient` | n uint8 | 1 if node carries Dirichlet pressure identity row |
| `d_ull4` | 4 uint64 | IEEE-754 bit-trick accumulators for maxVelChange/Vel/PresChange/Pres |

**Step 6b — why it matters**: The CPU assembles the pressure RHS after calling
`solvePressure()`, which in turn calls `freeSurfacePressureBoundary()` →
`currentBoundaryDeviatoricStress(v*)` — computing τ fresh from v* (post-Stokes).
Step 4 had computed d_ambientBP from τ(v^k) (pre-Stokes). Without step 6b, the GPU's
pressure RHS ambient BC would lag by one iteration, accumulating a systematic offset
in `d_stressHist` and eventually causing pressure divergence.

**Ghost velocity at AMBIENT faces** (`gpuSimpleFaceVelBC()`): `d_stress` stores σ = τ−p·I,
so:
```
devTraction = Σ_j σ_{comp,j} · n_j = τ_comp·n − p·n_comp
normalTraction = p·n_comp − τ_comp·n = −devTraction
ghostVel[comp] = interiorVel[comp] + offset · dist · normalTraction · n[dir] / η
```
Using `-devTraction` directly avoids adding and then subtracting p·n_comp (which
would introduce a double-pressure error when p ≠ 0).

### 6.6 GPU Activation Logic

GPU support is conditionally compiled under `VIENNALS_GPU_BICGSTAB`. When enabled,
`setupDeformationGpuBuffers()` is called from `buildNodes()` and allocates/initialises
the following CUDA handles:

| Handle | Purpose |
|---|---|
| `gpuPressBufs_` | BiCGSTAB buffers for the pressure Poisson solve |
| `gpuStokesBufs_` | BiCGSTAB buffers for the Stokes velocity solve (per component) |
| `gpuHarmonicBufs_` | BiCGSTAB buffers for the harmonic extension solve |
| `gpuSimpleBufs_` | Full GPU SIMPLE loop device state (`GpuSimpleBuffers`) |

**Activation conditions**:
```
tryGpu  = (gpuMode_ == GpuMode::Gpu)           // forced
        || (gpuMode_ == GpuMode::Auto && n >= kGpuThreshold)   // kGpuThreshold = 20000

gpuPressBufs_, gpuStokesBufs_ allocated when: tryGpu
gpuSimpleBufs_ allocated when:               tryGpu && !useIlu0
                                              (ILU0 not supported for GPU SIMPLE loop)
```

`GpuMode::Gpu` forces GPU for any node count; `GpuMode::Auto` uses CPU for small
problems and switches to GPU at ≥ 20 000 nodes. When `gpuSimpleBufs_` is valid, the
**entire** SIMPLE outer loop runs on-device (zero per-iteration CPU work). When only
`gpuPressBufs_`/`gpuStokesBufs_` are valid (e.g. when ILU0 is requested), the CPU
assembles the RHS each iteration and calls the GPU only for the linear solve.

**Preconditioner options** (`GpuPreconditioner`):
- `Jacobi` (default) — diagonal scaling, supports both individual solves and GPU SIMPLE loop
- `ILU0` — incomplete LU(0) factorization via CUSPARSE, better convergence for
  ill-conditioned systems; supported only for individual Stokes/pressure solves,
  **not** for the GPU SIMPLE loop (per-component re-factorization not implemented)

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
| `lsOxidationDeformation.hpp` | Stokes flow, SIMPLE pressure-velocity coupling, oxide velocity; CPU + GPU dispatch |
| `lsOxidationSimpleGpu.cuh` | GPU SIMPLE loop: CUDA kernels + `GpuSimpleBuffers`; nvcc-only |
| `lsOxidationBiCGSTABInterface.hpp` | C++ (g++) interface to GPU BiCGSTAB; safe for .cpp inclusion |
| `lsOxidationBiCGSTAB.cuh` | GPU BiCGSTAB kernels (`solveBiCGSTABInPlace`); nvcc-only |
| `lsOxidationMask.hpp` | Mask bending, multigrid, contact mechanics, traction BC |
| `lsOxidationModel.hpp` | Pressure–concentration coupling loop (Aitken Δ²) |
| `lsOxidation.hpp` | Per-step orchestrator: CFL stepping, advection, interior fill |
| `lsOxidationPresets.hpp` | Named parameter presets (wet/dry at 1000 °C, SiN mask) |
| `psOxidation.hpp` (ViennaPS) | High-level process model: Arrhenius table, domain management |
