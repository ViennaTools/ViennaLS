# StepOxidation Example

This example demonstrates a two-level-set silicon oxidation workflow using
`lsOxidationModel.hpp`. It is based on the Cartesian finite-difference idea used
in the Suvorov, Hoessinger, Djuric, and Ljepojevic oxidation model: the moving
interfaces stay represented as level sets, while the oxidant concentration and
oxide deformation quantities are solved on the Cartesian grid that carries those
level sets.

## Geometry

The simulation uses two level sets:

- `siInterface`: the Si/SiO2 reaction interface.
- `ambientInterface`: the SiO2/ambient free surface.

The oxide is the band between these two surfaces. In the default sign
convention used by the model, grid nodes belong to the oxide when

```text
reactionPhi >= 0
ambientPhi  <= 0
```

The initial silicon shape is a 2D step. The bottom of the wafer is not generated
as a surface; the silicon level set is built from an infinite plane plus a raised
right-hand block. The initial oxide/ambient surface is then created by
geometrically offsetting the Si/SiO2 interface by `oxideThickness`, which avoids
overlapping sidewalls at the step.

## Oxidant Diffusion Solve

`OxidationDiffusionVelocityField` solves a steady diffusion problem inside the
oxide band:

```text
div(D grad C) = 0
```

where:

- `C` is oxidant concentration in the oxide.
- `D` is `diffusionCoefficient`.

Only Cartesian grid nodes inside the oxide are solved. If a finite-difference
edge leaves the oxide, the code inserts an axis-aligned cross point at the
nearest crossed level-set interface and applies the appropriate boundary
condition there. The distance from the grid node to that cross point is kept as
a sub-grid distance.

At the Si/SiO2 reaction interface, a Robin reaction boundary is used:

```text
-D dC/dn = k C
```

where `k` is `reactionRate`.

At the oxide/ambient interface, a gas-transfer Robin boundary is used:

```text
-D dC/dn = h (C* - C)
```

where:

- `h` is `transferCoefficient`.
- `C*` is `equilibriumConcentration`.

The Cartesian update is assembled as a nonuniform three-point stencil along each
axis. For one axis with distances `d-` and `d+` from the current node to the
minus and plus stencil points, the second derivative is discretized as:

```text
2 / (d- + d+) * ((C+ - C0) / d+ - (C0 - C-) / d-)
```

For ordinary neighboring oxide grid nodes, `d` is the grid spacing and the
neighbor concentration is the unknown from the previous iteration. For a
cross-point, the Robin condition is used to eliminate the boundary
concentration. For example, at the reaction interface:

```text
g = D / d
Cb = g C0 / (g + k)
```

and at the ambient interface:

```text
g = D / d
Cb = (g C0 + h C*) / (g + h)
```

The resulting `Cb = a C0 + b` form is inserted into the nonuniform stencil
implicitly, so the local diagonal coefficient and right-hand side both reflect
the actual sub-grid cross-point distance. This is the part of the example that
matches the paper's idea of adding extra points on the Cartesian axes.

The diffusion field also supports an optional mask level set through
`setMaskInterface(mask, sign)`. Grid nodes inside the mask are removed from the
oxide solve region. If a stencil crosses into the mask, the mask boundary is used
instead of the ambient boundary. The default `maskTransferCoefficient = 0` makes
the mask oxidant-blocking:

```text
D grad(C) . n = 0
```

Setting `maskTransferCoefficient > 0` gives an imperfect-blocking Robin boundary
toward `maskConcentration`, which can represent a leaky cap or another material
in the stack.

The Si/SiO2 interface velocity returned by the diffusion field is

```text
v_Si = sign * k C / (N gamma)
```

where:

- `N` is `oxidantMoleculeDensity`.
- `gamma` is `expansionCoefficient`.
- `sign` is `velocitySign`; this example uses `-1` so silicon is consumed.

For the example parameters, `C*/N` is normalized to 1. The values are chosen so
that the Deal-Grove coefficients approximately correspond to wet oxidation of
`<100>` silicon at 1000 C:

```text
D = 0.157 um^2/hr
k = 0.74  um/hr
h = 100   um/hr
gamma = 2.27
```

With this normalization, `B = 2D` and `B/A ~= k` when gas transfer is large.

## Oxide Deformation Solve

`OxidationDeformationVelocityField` propagates the volume expansion generated at
the Si/SiO2 interface through the oxide. The local expansion speed associated
with silicon consumption is

```text
v_exp = ((gamma - 1) / gamma) * k C / N
```

The mechanics solve starts with a component-wise harmonic extension of this
reaction-boundary velocity through the oxide band. This is used as a predictor
velocity. The model then performs a fixed-point quasi-static Stokes update on
the same Cartesian oxide grid.

First, it solves a pressure equation from the current velocity divergence:

```text
laplace(p) ~= -bulkModulus * div(v)
```

The pressure Laplacian is evaluated with nonuniform finite-difference stencils.
If a Cartesian axis leaves the oxide before reaching the next grid point, the
level-set crossing is used as a sub-grid boundary point. Very small cut cells are
clamped by `minMechanicsBoundaryDistance` to avoid a sliver cell dominating the
mechanical solve.

At the oxide/ambient surface, the pressure boundary is chosen from the normal
traction condition:

```text
p_surface = p_ambient + freeSurfaceTractionScale * n . s . n
```

where `s` is the current deviatoric stress estimate. With
`freeSurfaceTractionScale = 1`, this is the scalar normal part of a traction-free
boundary. At the Si/SiO2 side, the pressure boundary can include an elastic
silicon support:

```text
p_reaction = p_local + substrateNormalStiffness * dt * (v_reaction . n)
```

Setting `substrateNormalStiffness = 0` recovers the older no-flux pressure
treatment.

If a mask level set is attached to the deformation field, mask crossings are
treated as constrained material-stack boundaries rather than oxide/ambient free
surfaces. The default `maskVelocityScale = 0` is a no-slip mask contact for the
oxide velocity. `maskNormalStiffness` adds elastic normal support, and
`maskPressure` is the reference pressure used when that support is disabled.

The velocity is then relaxed with a finite-difference momentum equation:

```text
viscosity * laplace(v) = pressureGradientScale * (grad(p) - div(s))
strain_trace = div(v)
```

Here `s` is the current deviatoric stress estimate. The velocity Laplacian,
pressure gradient, velocity divergence, and strain-rate tensor use the same
sub-grid boundary distances. Reaction-boundary grid crossings use the oxidation
expansion velocity as a Dirichlet condition. Ambient boundary crossings use a
ghost velocity from the same approximate traction-free condition, so the free
surface responds to both pressure and deviatoric stress rather than receiving
only an after-the-fact pressure-gradient correction.

This pressure/velocity update is repeated for `mechanicsIterations` until the
relative pressure and velocity changes fall below `mechanicsTolerance`. This is
still a compact finite-difference approximation, but pressure, velocity, and
deviatoric stress now participate in the mechanics balance.

`pressureIterations` controls the pressure solve, and `stokesIterations` controls
the velocity momentum relaxation. They are separate because the harmonic
predictor and Stokes correction have different convergence costs.

After each velocity update, the model computes the symmetric strain-rate tensor:

```text
D_ij = 0.5 * (dv_i/dx_j + dv_j/dx_i)
```

The deviatoric part is evolved with a Maxwell-style relaxation law:

```text
s_new = exp(-dt/tau) * s_old
        + (1 - exp(-dt/tau)) * 2 * viscosity * dev(D)
```

where `tau` is either `stressRelaxationTime` or `viscosity / shearModulus` if a
shear modulus is provided. The full Cauchy stress tensor is then:

```text
sigma = -p I + s_new
```

The example writes selected tensor components and von Mises stress for
inspection. This is still a compact finite-difference approximation rather than
a fully calibrated viscoelastic process model, but the stress tensor is now a
state variable with relaxation memory.

## Pressure-Coupled Iteration

The example now uses `OxidationCoupledModel` to iterate the one-way fields into a
first pressure-coupled loop:

```text
solve oxidant diffusion
solve oxide deformation and pressure
feed pressure back into the diffusion reaction rate
repeat until the relative pressure change is small
```

The stress feedback enters the reaction boundary through an effective reaction
rate:

```text
k_eff = k * clamp(exp(-alpha * (p - p_ref)), minFactor, maxFactor)
```

where:

- `alpha` is `stressCouplingCoefficient`.
- `p_ref` is `referencePressure`.
- `minFactor` and `maxFactor` bound the correction.

With positive `alpha`, positive pressure reduces the local reaction rate. The
coefficient in this example is deliberately small because the pressure model is
still a compact finite-difference approximation rather than a calibrated
viscoelastic stress model.

The example writes these diagnostics to:

```text
step_oxidation_deformation.csv
```

with columns:

```text
x,y,velocity_x,velocity_y,pressure,strain_trace,stress_xx,stress_xy,stress_yy,von_mises_stress
```

## Free-Surface Motion and the 0.44:0.56 Split

Silicon oxidation expands the material volume. For an expansion coefficient
`gamma`, the ideal flat-interface displacement split is:

```text
Si/SiO2 inward fraction      = 1 / gamma
SiO2/ambient outward fraction = (gamma - 1) / gamma
```

For `gamma = 2.27`, this gives approximately:

```text
Si/SiO2      = 0.4405
SiO2/ambient = 0.5595
```

The example prints this benchmark and also samples the flat top of the step to
verify that the local velocity field has the same split.

The free surface uses `freeSurfaceVelocityScale = 1.0`. The harmonic vector
deformation field is still solved and written for diagnostics, but the example
sets:

```text
vectorVelocityScale = 0.0
```

This avoids double-counting the free-surface motion during advection. In other
words, the ambient level set is advected by the kinematic volume-expansion split,
while the vector field remains available for pressure and strain inspection.

To evaluate the local free-surface speed on a non-planar shape, the deformation
field projects the oxide/ambient point back to the nearest Si/SiO2 reaction
crossing along both normal directions and samples the reaction speed just inside
the oxide. This makes the free-surface speed local to the corresponding reaction
interface instead of using a global average.

## Output

Running the example writes:

```text
step_oxidation_si_initial.vtk
step_oxidation_ambient_initial.vtk
step_oxidation_si_after.vtk
step_oxidation_ambient_after.vtk
step_oxidation_deformation.csv
```

The VTK files contain the initial and final interfaces. The CSV contains the
deformation velocity, pressure, and strain-trace diagnostics on the Cartesian
solve grid.

## Running

From the `ViennaLS` build directory:

```bash
cmake --build build --target StepOxidation
./build/examples/StepOxidation/StepOxidation
```

The example currently uses `advectionTime = 0.1 hr`. This keeps the flat-region
displacement benchmark close to the imposed local velocity split on the current
grid. Larger one-shot advection times can introduce visible level-set
discretization error, especially near the step corner.
