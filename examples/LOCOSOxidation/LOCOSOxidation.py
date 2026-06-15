#!/usr/bin/env python3
"""
LOCOS (Local Oxidation of Silicon) example — ViennaLS low-level API.

Geometry (2-D cross-section):
    · Si substrate at y = 0
    · Pad SiO2 of thickness padOxideThickness grown by geometric advect
    · Si3N4 mask box covering x ∈ [−xExtent, maskEdge], sitting on pad oxide

This script uses the ViennaLS Oxidation solver directly (lsOxidation),
bypassing ViennaPS.  It mirrors LOCOSOxidation.cpp exactly.

Usage:
    python LOCOSOxidation.py

All lengths are in micrometers, time in hours, pressure in Pa.
"""

import math
import sys
import csv

import viennals as vls

vls.setDimension(2)

# Pass "gpu" as first argument to use the GPU BiCGSTAB solver, e.g.:
#   python LOCOSOxidation.py gpu
# Requires a GPU build (-DVIENNALS_USE_GPU=ON). CPU is the default.
USE_GPU = len(sys.argv) > 1 and sys.argv[1].lower() == "gpu"

# ── Geometry constants ────────────────────────────────────────────────────────
GRID_DELTA          = 0.01   # um
X_EXTENT            = 1.0    # um  (REFLECTIVE at x = 0)
Y_MIN               = -1.0   # um
Y_MAX               = 2.0    # um
PAD_OXIDE_THICKNESS = 0.03   # um
MASK_THICKNESS      = 0.05   # um
MASK_EDGE           = 0.0    # um  (open window: x > 0)
ADVECTION_TIME      = 0.1    # hr
TIME_STEP           = 0.01   # hr per outer step
MASK_CONTACT_EPS    = 1.e-6  # um


def to_index(x: float) -> int:
    return int(round(x / GRID_DELTA))


def write_surface(level_set, file_name: str) -> None:
    if level_set.getLevelSetWidth() < 2:
        vls.Expand(level_set, 2).apply()
    mesh = vls.Mesh()
    surf = vls.ToSurfaceMesh(level_set, mesh)
    surf.setSharpCorners(False)
    surf.apply()
    vls.VTKWriter(mesh, file_name).apply()


# ── Domain setup ──────────────────────────────────────────────────────────────
bounds = [-X_EXTENT, X_EXTENT, Y_MIN, Y_MAX]
bcs    = [vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
          vls.BoundaryConditionEnum.INFINITE_BOUNDARY]

# Si substrate at y = 0.
si_interface = vls.Domain(bounds, bcs, GRID_DELTA)
vls.MakeGeometry(
    si_interface,
    vls.Plane([0., 0.], [0., 1.]),
).apply()

# Pad SiO2: geometric offset of Si by padOxideThickness.
ambient_interface = vls.Domain(si_interface)
vls.GeometricAdvect(
    ambient_interface,
    vls.SphereDistribution(PAD_OXIDE_THICKNESS),
).apply()

# Capture pre-oxidation state for volume conservation check.
si_initial      = vls.Domain(si_interface)
ambient_initial = vls.Domain(ambient_interface)

# Si3N4 mask box: covers x ∈ [−xExtent, maskEdge], sitting on the pad oxide.
mask_interface = vls.Domain(bounds, bcs, GRID_DELTA)
mask_geom = vls.MakeGeometry(
    mask_interface,
    vls.Box(
        [-X_EXTENT,  PAD_OXIDE_THICKNESS - MASK_CONTACT_EPS],
        [MASK_EDGE,  PAD_OXIDE_THICKNESS + MASK_THICKNESS],
    ),
)
mask_geom.setIgnoreBoundaryConditions([False, True, False])
mask_geom.apply()

write_surface(si_interface,     "locos_si_initial.vtp")
write_surface(ambient_interface, "locos_ambient_initial.vtp")
write_surface(mask_interface,   "locos_mask.vtp")

# ── Oxidation parameters ──────────────────────────────────────────────────────
ox_params = vls.OxidationPresets.wet1000CDealGrove100()
ox_params.velocitySign              = -1.
ox_params.maskTransferCoefficient   = 0.
ox_params.maskConcentration         = 0.
ox_params.maxIterations             = 10000
ox_params.tolerance                 = 1.e-7

def_params = vls.OxidationPresets.oxideMechanics1000C(TIME_STEP)

coupling_params = vls.OxidationCouplingParameters()
coupling_params.maxIterations = 100
coupling_params.tolerance     = 2.e-2
coupling_params.relaxation    = 1.

mask_params = vls.OxidationPresets.siliconNitrideMask1000C()

# ── Solve bounds ──────────────────────────────────────────────────────────────
diff_min = [to_index(-X_EXTENT), to_index(Y_MIN)]
diff_max = [to_index(X_EXTENT),  to_index(Y_MAX)]

mask_min = [to_index(-X_EXTENT),                              to_index(PAD_OXIDE_THICKNESS) - 1]
mask_max = [to_index(MASK_EDGE), to_index(PAD_OXIDE_THICKNESS + MASK_THICKNESS) + 1]

# ── LOCOS solver ──────────────────────────────────────────────────────────────
locos = vls.Oxidation(si_interface, ambient_interface, mask_interface)
locos.setOxidationParameters(ox_params)
locos.setDeformationParameters(def_params)
locos.setCouplingParameters(coupling_params)
locos.setMaskParameters(mask_params)
locos.setSolveBounds(diff_min, diff_max)
locos.setMaskBendingBounds(mask_min, mask_max)
locos.setMaskCouplingIterations(30)
locos.setMaskCouplingTolerance(1.e-2)
if USE_GPU:
    locos.setGpuMode(vls.GpuMode.Gpu)

# ── Time-stepping loop ────────────────────────────────────────────────────────
elapsed  = 0.0
time_eps = ADVECTION_TIME * 1.e-8
while ADVECTION_TIME - elapsed > time_eps:
    dt       = min(TIME_STEP, ADVECTION_TIME - elapsed)
    elapsed += locos.applyCFLLimited(dt, 0.499)
    print(f"t = {elapsed:.6f} hr")

# ── Diagnostics ───────────────────────────────────────────────────────────────
diffusion   = locos.getDiffusionField()
deformation = locos.getDeformationField()
mask_bending = locos.getMaskBendingField()

open_pt   = [0.5,  GRID_DELTA, 0.]
masked_pt = [-0.5, GRID_DELTA, 0.]
normal_up = [0., 1., 0.]

open_conc           = diffusion.getConcentration(open_pt)
masked_conc         = diffusion.getConcentration(masked_pt)
open_si_speed       = abs(diffusion.getScalarVelocity(open_pt,   0, normal_up, 0))
masked_si_speed     = abs(diffusion.getScalarVelocity(masked_pt, 0, normal_up, 0))
suppression_ratio   = (masked_si_speed / open_si_speed) if open_si_speed > 0. else 0.

print(f"Diffusion   nodes: {diffusion.getNumberOfSolutionNodes()}, "
      f"iters: {diffusion.getIterations()}, "
      f"residual: {diffusion.getResidual():.3e}")
print(f"Deformation nodes: {deformation.getNumberOfSolutionNodes()}, "
      f"iters: {deformation.getIterations()}, "
      f"residual: {deformation.getResidual():.3e}")
print(f"Average oxide expansion speed: {deformation.avgExpansionSpeed():.6f} um/hr")
print(f"Open-window concentration: {open_conc:.6f}, "
      f"masked concentration: {masked_conc:.6f}")
print(f"Open-window Si speed: {open_si_speed:.6f} um/hr, "
      f"masked Si speed: {masked_si_speed:.6f} um/hr, "
      f"suppression ratio: {suppression_ratio:.4f}")

if not math.isfinite(open_si_speed) or not math.isfinite(masked_si_speed) \
        or open_si_speed <= 0. or suppression_ratio > 0.05:
    print("LOCOS mask sanity check FAILED: masked oxidation is not sufficiently suppressed.",
          file=sys.stderr)
    sys.exit(1)

# Write diagnostics CSV.
with open("locos_oxidation_diagnostics.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["x", "y", "concentration", "velocity_x", "velocity_y",
                     "pressure", "strain_trace", "von_mises_stress"])
    for j in range(diff_min[1], diff_max[1] + 1):
        for i in range(diff_min[0], diff_max[0] + 1):
            coord = [i * GRID_DELTA, j * GRID_DELTA, 0.]
            vel   = deformation.getVelocity(coord)
            writer.writerow([
                coord[0], coord[1],
                diffusion.getConcentration(coord),
                vel[0], vel[1],
                deformation.getPressure(coord),
                deformation.getStrainTrace(coord),
                deformation.getVonMisesStress(coord),
            ])

# Mask bending samples.
normal_dn  = [0., -1., 0.]
mask_pts = {
    "bottom":       ([-0.5, PAD_OXIDE_THICKNESS,                         0.], normal_dn),
    "top":          ([-0.5, PAD_OXIDE_THICKNESS + MASK_THICKNESS,        0.], normal_up),
    "contact_node": ([-GRID_DELTA, PAD_OXIDE_THICKNESS,                  0.], normal_dn),
    "anchor":       ([-X_EXTENT + GRID_DELTA,
                      PAD_OXIDE_THICKNESS + MASK_THICKNESS * 0.5,        0.], normal_up),
    "mid":          ([-0.5, PAD_OXIDE_THICKNESS + MASK_THICKNESS * 0.5,  0.], normal_up),
    "edge":         ([-GRID_DELTA, PAD_OXIDE_THICKNESS + MASK_THICKNESS * 0.5, 0.], normal_up),
}

vels = {name: mask_bending.getVectorVelocity(pt, 0, n, 0)
        for name, (pt, n) in mask_pts.items()}
mask_bottom_pressure = deformation.getPressure(mask_pts["bottom"][0])

print(f"Mask elasticity nodes: {mask_bending.getNumberOfSolutionNodes()}, "
      f"contact nodes: {mask_bending.getNumberOfContactNodes()}, "
      f"iters: {mask_bending.getIterations()}, "
      f"residual: {mask_bending.getResidual():.3e}")
print(f"Oxide/mask interface iters: {locos.getMaskCouplingIterations()}, "
      f"residual: {locos.getMaskCouplingResidual():.3e}")
v = vels["bottom"]
print(f"Mask bottom velocity: ({v[0]:.6f}, {v[1]:.6f}) um/hr, "
      f"top velocity: ({vels['top'][0]:.6f}, {vels['top'][1]:.6f}) um/hr, "
      f"oxide pressure: {mask_bottom_pressure:.3e} Pa")
vc = vels["contact_node"]
print(f"Mask contact-node velocity: ({vc[0]:.6f}, {vc[1]:.6f}) um/hr")
va, vm, ve = vels["anchor"], vels["mid"], vels["edge"]
print(f"Mask lateral bending: anchor ({va[0]:.6f}, {va[1]:.6f}), "
      f"mid ({vm[0]:.6f}, {vm[1]:.6f}), "
      f"edge ({ve[0]:.6f}, {ve[1]:.6f}) um/hr")

# ── Volume conservation ───────────────────────────────────────────────────────
conservation = vls.computeLOCOSOpenWindowConservation(
    si_initial, si_interface,
    ambient_initial, ambient_interface,
    0.1, 0.9,
    to_index(Y_MIN), to_index(Y_MAX),
    ox_params.expansionCoefficient,
)
print(f"Open-window conservation samples: {conservation.samples}, "
      f"Si consumed area: {conservation.siliconRecession:.6f} um^2, "
      f"ambient lift area: {conservation.ambientLift:.6f} um^2, "
      f"expected ambient lift: {conservation.expectedAmbientLift:.6f} um^2, "
      f"lift/Si ratio: {conservation.ambientLiftRatio:.4f} "
      f"(expected {ox_params.expansionCoefficient - 1.:.4f}), "
      f"relative error: {conservation.relativeError:.4f}")

# ── Final surface meshes ──────────────────────────────────────────────────────
write_surface(si_interface,     "locos_si_after.vtp")
write_surface(ambient_interface, "locos_ambient_after.vtp")
write_surface(mask_interface,   "locos_mask_after.vtp")

print("Wrote locos_si_initial.vtp, locos_ambient_initial.vtp, locos_mask.vtp, "
      "locos_si_after.vtp, locos_ambient_after.vtp, locos_mask_after.vtp, "
      "and locos_oxidation_diagnostics.csv")
