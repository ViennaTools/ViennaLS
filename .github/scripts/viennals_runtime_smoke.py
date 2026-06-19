from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import textwrap


SMOKE_TEST = r"""
from pathlib import Path
import tempfile

import viennals as vls

print(vls.__file__)
print(vls.d2)
print(vls.d3)

tmp = Path(tempfile.mkdtemp(prefix="viennals-runtime-smoke-"))

# Exercise direct VTK mesh writing from the installed wheel.
mesh = vls.Mesh()
mesh.insertNextNode([0.0, 0.0, 0.0])
mesh.insertNextNode([1.0, 0.0, 0.0])
mesh.insertNextNode([0.0, 1.0, 0.0])
mesh.insertNextTriangle([0, 1, 2])
mesh_file = tmp / "mesh.vtp"
vls.VTKWriter(mesh, str(mesh_file)).apply()

if not mesh_file.exists() or mesh_file.stat().st_size == 0:
    raise RuntimeError(f"VTKWriter did not create {mesh_file}")

# Exercise a small level-set-to-surface path before writing. This hits more of
# the bundled VTK runtime than importing the module or writing a hand-built mesh.
d2 = vls.d2
bounds = [-1.0, 1.0, -1.0, 1.0]
bcs = [
    vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
    vls.BoundaryConditionEnum.INFINITE_BOUNDARY,
]
domain = d2.Domain(bounds, bcs, 0.1)
d2.MakeGeometry(domain, d2.Sphere([0.0, 0.0], 0.5)).apply()

surface = vls.Mesh()
d2.ToSurfaceMesh(domain, surface).apply()
surface_file = tmp / "surface.vtp"
vls.VTKWriter(surface, str(surface_file)).apply()

if not surface_file.exists() or surface_file.stat().st_size == 0:
    raise RuntimeError(f"VTKWriter did not create {surface_file}")

# Exercise the higher-level VTK visualization writer when the wheel was built
# with VTK support.
d3 = vls.d3
if hasattr(d3, "WriteVisualizationMesh"):
    vis_domain = d3.Domain(0.25)
    d3.MakeGeometry(vis_domain, d3.Sphere([0.0, 0.0, 0.0], 1.0)).apply()
    vis_writer = d3.WriteVisualizationMesh(vis_domain)
    vis_writer.setFileName(str(tmp / "visualization"))
    vis_writer.apply()
"""


def main() -> int:
    env = os.environ.copy()
    # The test should catch runtime OpenMP warnings, not inherit a developer's
    # local suppression setting.
    env.pop("KMP_WARNINGS", None)

    result = subprocess.run(
        [sys.executable, "-c", SMOKE_TEST],
        capture_output=True,
        cwd=tempfile.mkdtemp(prefix="viennals-runtime-smoke-run-"),
        env=env,
        text=True,
    )

    if result.stdout:
        print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="", file=sys.stderr)

    combined_output = result.stdout + result.stderr
    if "omp_set_nested" in combined_output:
        print(
            textwrap.dedent(
                """
                ViennaLS runtime smoke test detected deprecated OpenMP nested
                parallelism output. This usually means a runtime path still
                called omp_set_nested instead of omp_set_max_active_levels.
                """
            ).strip(),
            file=sys.stderr,
        )
        return 1

    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
