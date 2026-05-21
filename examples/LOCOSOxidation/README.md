# LOCOS Oxidation Example

This example exercises the oxidation mask boundary conditions. It creates a
flat silicon surface with a thin pad oxide and a finite-thickness nitride-like
mask covering the left half of the structure. The mask bottom lies flat on the
pad-oxide/ambient interface, with only a tiny sub-grid contact offset used
internally to make the cut-cell boundary unambiguous. The right half is an open
oxidation window.

The parameters are intentionally set as a demonstration dose rather than a
calibrated process recipe: the oxidation time and mask pressure-to-velocity
coupling are large enough that mask lift/bending is visible in the output VTK
files.

The mask is attached to both oxidation fields:

```cpp
oxidationVelocity->setMaskInterface(maskInterface, -1);
deformationVelocity->setMaskInterface(maskInterface, -1);
```

For diffusion, `maskTransferCoefficient = 0` makes the mask oxidant-blocking.
For oxide mechanics, `maskVelocityScale` and `maskNormalStiffness` define the
contact condition seen by the growing oxide. After oxidizing the pad oxide, the
example advects the mask with `OxidationMaskBendingVelocityField`: it builds a
Cartesian vector-displacement solve inside the SiN mask level set, imposes the
oxide contact velocity and pressure on the oxide-facing mask boundary, and
relaxes the displacement-rate field through the mask volume with linear-elastic
Lamé parameters derived from the mask Young's modulus and Poisson ratio. The
whole mask level set is then advected with that vector field, so the top and
bottom surfaces move coherently. A thicker or stiffer SiN mask bends less;
stronger oxidation pressure near the mask edge bends it more.

Run from the build directory:

```bash
./examples/LOCOSOxidation/LOCOSOxidation
```

The example prints oxidant concentration and silicon consumption speed in the
open window and under the mask. The masked values should be suppressed relative
to the open-window values. It also prints the number of mask elasticity nodes,
the number of oxide-contact boundary nodes, and sample velocities at the bottom
and top of the SiN mask to verify that the full mask body is being moved by a
vector field. The masked oxide/mask contact reports zero scalar free-surface
growth; its motion comes from the mask vector velocity, so the oxide cannot grow
through the mask and the expansion is redirected toward exposed ambient regions.

Outputs:

- `locos_si_initial.vtk`
- `locos_ambient_initial.vtk`
- `locos_mask.vtk`
- `locos_si_after.vtk`
- `locos_ambient_after.vtk`
- `locos_mask_after.vtk`
- `locos_oxidation_diagnostics.csv`
