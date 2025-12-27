import viennals as vls
import argparse

vls.Logger.setLogLevel(vls.LogLevel.INFO)


# 1. Define the Velocity Field
# We inherit from viennals.VelocityField to define custom logic in Python
class EtchingField(vls.VelocityField):
    def getScalarVelocity(self, coordinate, material, normalVector, pointId):
        if material == 2:
            return -0.1
        if material == 1:
            return -1.0
        return 0.0

    def getVectorVelocity(self, coordinate, material, normalVector, pointId):
        return [0.0] * len(coordinate)


# 1. Define the Velocity Field
# We inherit from viennals.VelocityField to define custom logic in Python
class DepositionField(vls.VelocityField):
    def getScalarVelocity(self, coordinate, material, normalVector, pointId):
        return 0.75

    def getVectorVelocity(self, coordinate, material, normalVector, pointId):
        return [0.0] * len(coordinate)


def main():
    # 1. Parse Arguments
    parser = argparse.ArgumentParser(
        description="Run Square Etch simulation in 2D or 3D."
    )
    parser.add_argument(
        "-D",
        "--dim",
        type=int,
        default=2,
        choices=[2, 3],
        help="Dimension of the simulation (2 or 3). Default is 2.",
    )
    args = parser.parse_args()

    DIMENSION = args.dim
    vls.setDimension(DIMENSION)
    vls.setNumThreads(8)

    extent = 30.0
    gridDelta = 0.47

    # Define bounds and boundary conditions
    trenchBottom = -2.0
    bounds = [-extent, extent, -extent, extent]
    origin = [0.0, 0.0]
    if DIMENSION == 3:
        bounds.append(-extent)
        bounds.append(extent)
        origin.append(0.0)
    boundaryCons = []
    for i in range(DIMENSION - 1):
        boundaryCons.append(vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY)
    boundaryCons.append(vls.BoundaryConditionEnum.INFINITE_BOUNDARY)

    planeNormal = list(origin)
    downNormal = list(origin)
    planeNormal[DIMENSION - 1] = 1.0
    downNormal[DIMENSION - 1] = -1.0

    # Create the Substrate Domain
    substrate = vls.Domain(bounds, boundaryCons, gridDelta)

    # Create initial flat geometry (Plane at y/z=0)
    plane = vls.Plane(origin, planeNormal)
    vls.MakeGeometry(substrate, plane).apply()

    # --------------------------------------
    # 3. Trench Geometry
    # --------------------------------------
    trench = vls.Domain(bounds, boundaryCons, gridDelta)

    # Define Box Corners based on dimension
    minCorner = list(origin)
    maxCorner = list(origin)
    originMask = list(origin)
    for i in range(DIMENSION - 1):
        minCorner[i] = -extent / 1.5
        maxCorner[i] = extent / 1.5
    minCorner[DIMENSION - 1] = trenchBottom
    maxCorner[DIMENSION - 1] = 1.0
    originMask[DIMENSION - 1] = trenchBottom + 1e-9

    box = vls.Box(minCorner, maxCorner)
    vls.MakeGeometry(trench, box).apply()

    # Subtract trench from substrate (Relative Complement)
    vls.BooleanOperation(
        substrate, trench, vls.BooleanOperationEnum.RELATIVE_COMPLEMENT
    ).apply()

    # 4. Create the Mask Layer
    # The mask prevents etching at the bottom of the trench in specific configurations,
    # or acts as a hard stop/masking material.
    mask = vls.Domain(bounds, boundaryCons, gridDelta)

    vls.MakeGeometry(mask, vls.Plane(originMask, downNormal)).apply()

    # Intersect with substrate geometry
    vls.BooleanOperation(mask, substrate, vls.BooleanOperationEnum.INTERSECT).apply()

    # 5. Export Initial State
    print(f"Running in {DIMENSION}D.")
    print("Extracting initial meshes...")
    mesh = vls.Mesh()

    # Save substrate
    vls.ToSurfaceMesh(substrate, mesh).apply()
    vls.VTKWriter(mesh, f"substrate-{DIMENSION}D.vtp").apply()

    # Save mask
    vls.ToSurfaceMesh(mask, mesh).apply()
    vls.VTKWriter(mesh, f"mask-{DIMENSION}D.vtp").apply()

    print("Creating new layer...")
    polymer = vls.Domain(substrate)

    # 6. Setup Advection
    deposition = DepositionField()
    etching = EtchingField()

    print("Advecting")
    advectionKernel = vls.Advect()

    # the level set to be advected has to be inserted last
    # the other could be taken as a mask layer for advection
    advectionKernel.insertNextLevelSet(mask)
    advectionKernel.insertNextLevelSet(substrate)
    advectionKernel.insertNextLevelSet(polymer)

    advectionKernel.setSaveAdvectionVelocities(True)
    # advectionKernel.setVelocityField(etching)

    advectionKernel.setVelocityField(deposition)
    advectionKernel.apply()

    vls.ToSurfaceMesh(polymer, mesh).apply()
    vls.VTKWriter(mesh, f"newLayer-{DIMENSION}D.vtp").apply()

    advectionKernel.setSaveAdvectionVelocities(True)
    advectionKernel.setVelocityField(etching)

    # Use default spatial discretization scheme (Lax Friedrichs 1st order) as in the C++ else branch
    # advectionKernel.setSpatialScheme(viennals.SpatialSchemeEnum.LAX_FRIEDRICHS_1ST_ORDER)

    # 7. Time Loop
    finalTime = 10.0
    counter = 1
    currentTime = 0.0

    # The advect kernel calculates stable time steps automatically.
    # We call apply() repeatedly until we reach finalTime.

    # Note: In the C++ example, the loop structure is:
    # for (double time = 0.; time < finalTime; time += advectionKernel.getAdvectedTime())
    # This implies one step is taken, then time is updated.

    # Save initial state
    vls.ToSurfaceMesh(polymer, mesh).apply()
    vls.VTKWriter(mesh, f"numerical-{DIMENSION}D-0.vtp").apply()

    while currentTime < finalTime:
        advectionKernel.apply()

        stepTime = advectionKernel.getAdvectedTime()
        currentTime += stepTime

        print(
            f"Advection step: {counter}, time: {stepTime:.4f} (Total: {currentTime:.4f})"
        )

        # Export result
        vls.ToSurfaceMesh(polymer, mesh).apply()
        vls.VTKWriter(mesh, f"numerical-{DIMENSION}D-{counter}.vtp").apply()

        counter += 1

    print(f"\nNumber of Advection steps taken: {counter}")

    # Final export
    vls.ToSurfaceMesh(polymer, mesh).apply()
    vls.VTKWriter(mesh, f"final-{DIMENSION}D.vtp").apply()


if __name__ == "__main__":
    main()
