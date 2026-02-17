import viennals as vls


# @example AirGapDeposition.py
#  2D Example showing how to use the library for topography
#  simulation, by creating a trench geometry. A layer of a different material is
#  then grown directionally on top.


class velocityField(vls.VelocityField):
    # coord and normalVec are lists with 3 elements
    # in 2D coord[2] and normalVec[2] are zero
    # getScalarVelocity must return a scalar
    def getScalarVelocity(self, coord, material, normal, pointId):
        # velocity is proportional to the normal vector
        velocity = 0.0
        for i in range(len(normal)):
            velocity += abs(normal[i])
        return velocity

    def getVectorVelocity(self, coord, material, normal, pointId):
        return (0.0, 0.0, 0.0)


def run_simulation(kernel, new_layer, total_time, output_interval, name):
    passed_time = 0.0
    step_counter = 0
    mesh = vls.Mesh()

    while passed_time < total_time:
        dt = output_interval
        if passed_time + dt > total_time:
            dt = total_time - passed_time

        kernel.setAdvectionTime(dt)
        kernel.apply()
        passed_time += kernel.getAdvectedTime()

        print(
            f"\r{name} Advection time: {passed_time:04.1f} / {total_time:04.1f}s",
            end="",
            flush=True,
        )

        vls.ToSurfaceMesh(new_layer, mesh, minNodeDistFactor=0.02).apply()
        writer = vls.VTKWriter(mesh, f"trench_{name}_{step_counter}.vtp")
        writer.addMetaData("time", passed_time)
        writer.apply()

        vls.ToMesh(new_layer, mesh).apply()
        vls.VTKWriter(mesh, f"trench_{name}_{step_counter}.vtu").apply()

        step_counter += 1
    print()
    return passed_time


vls.setNumThreads(8)

extent = 30.0
gridDelta = 0.5

bounds = (-extent, extent, -extent, extent)
boundaryCons = (
    vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
    vls.BoundaryConditionEnum.INFINITE_BOUNDARY,
)

# create level set
substrate = vls.Domain(bounds, boundaryCons, gridDelta)

# create plane
origin = (0.0, 0.0)
planeNormal = (0.0, 1.0)

vls.MakeGeometry(substrate, vls.Plane(origin, planeNormal)).apply()

print("Extracting plane")
mesh = vls.Mesh()
vls.ToSurfaceMesh(substrate, mesh, minNodeDistFactor=0.02).apply()
vls.VTKWriter(mesh, "plane.vtp").apply()

# create layer used for booling
print("Creating box...")
trench = vls.Domain(bounds, boundaryCons, gridDelta)
minCorner = (-extent / 6.0, -25.0)
maxCorner = (extent / 6.0, 1.0)
vls.MakeGeometry(trench, vls.Box(minCorner, maxCorner)).apply()

print("Extracting box")
vls.ToMesh(trench, mesh).apply()
vls.VTKWriter(mesh, "box.vtu").apply()
surfaceMeshBox = vls.ToSurfaceMesh(trench, mesh, minNodeDistFactor=0.02)
surfaceMeshBox.setSharpCorners(True)
surfaceMeshBox.apply()
vls.VTKWriter(mesh, "box.vtp").apply()

# Create trench geometry
print("Booling trench")
vls.BooleanOperation(
    substrate, trench, vls.BooleanOperationEnum.RELATIVE_COMPLEMENT
).apply()

print("Extracting trench")
vls.ToMesh(substrate, mesh).apply()
vls.VTKWriter(mesh, "trench.vtu").apply()
surfaceMeshTrench = vls.ToSurfaceMesh(substrate, mesh, minNodeDistFactor=0.02)
surfaceMeshTrench.setSharpCorners(True)
surfaceMeshTrench.apply()
vls.VTKWriter(mesh, "trench.vtp").apply()

# Now grow new material

# Create copies for FE, RK2, and RK simulations
print("Creating new layers...")
substrateFE = vls.Domain(substrate)
newLayerFE = vls.Domain(substrateFE)

substrateRK2 = vls.Domain(substrate)
newLayerRK2 = vls.Domain(substrateRK2)

substrateRK = vls.Domain(substrate)
newLayerRK = vls.Domain(substrateRK)

velocities = velocityField()

totalSimulationTime = 16.5
outputInterval = 0.5

print("Advecting")

# FE Kernel
advectionKernelFE = vls.Advect()
advectionKernelFE.insertNextLevelSet(substrateFE)
advectionKernelFE.insertNextLevelSet(newLayerFE)
advectionKernelFE.setVelocityField(velocities)
advectionKernelFE.setIgnoreVoids(True)
advectionKernelFE.setTemporalScheme(vls.TemporalSchemeEnum.FORWARD_EULER)

passedTimeFE = run_simulation(
    advectionKernelFE, newLayerFE, totalSimulationTime, outputInterval, "FE"
)

# RK2 Kernel
advectionKernelRK2 = vls.Advect()
advectionKernelRK2.insertNextLevelSet(substrateRK2)
advectionKernelRK2.insertNextLevelSet(newLayerRK2)
advectionKernelRK2.setVelocityField(velocities)
advectionKernelRK2.setIgnoreVoids(True)
advectionKernelRK2.setTemporalScheme(vls.TemporalSchemeEnum.RUNGE_KUTTA_2ND_ORDER)

passedTimeRK2 = run_simulation(
    advectionKernelRK2, newLayerRK2, totalSimulationTime, outputInterval, "RK2"
)

# RK3 Kernel
advectionKernelRK = vls.Advect()
advectionKernelRK.insertNextLevelSet(substrateRK)
advectionKernelRK.insertNextLevelSet(newLayerRK)
advectionKernelRK.setVelocityField(velocities)
advectionKernelRK.setIgnoreVoids(True)
advectionKernelRK.setTemporalScheme(vls.TemporalSchemeEnum.RUNGE_KUTTA_3RD_ORDER)

passedTimeRK = run_simulation(
    advectionKernelRK, newLayerRK, totalSimulationTime, outputInterval, "RK3"
)

print(f"Time passed FE: {passedTimeFE}")
print(f"Time passed RK2: {passedTimeRK2}")
print(f"Time passed RK3: {passedTimeRK}")

# FE Output
writer = vls.WriteVisualizationMesh()
writer.insertNextLevelSet(substrateFE)
writer.insertNextLevelSet(newLayerFE)
writer.addMetaData("time", passedTimeFE)
writer.setFileName("airgap_FE")
writer.setExtractHullMesh(True)
writer.apply()

multiMeshKernel = vls.ToMultiSurfaceMesh(minNodeDistFactor=0.02)
multiMeshKernel.insertNextLevelSet(substrateFE)
multiMeshKernel.insertNextLevelSet(newLayerFE)
multiMesh = vls.Mesh()
multiMeshKernel.setMesh(multiMesh)
multiMeshKernel.setSharpCorners(True)
multiMeshKernel.apply()
vls.VTKWriter(multiMesh, "multimesh_FE.vtp").apply()

# RK2 Output
writer = vls.WriteVisualizationMesh()
writer.insertNextLevelSet(substrateRK2)
writer.insertNextLevelSet(newLayerRK2)
writer.addMetaData("time", passedTimeRK2)
writer.setFileName("airgap_RK2")
writer.setExtractHullMesh(True)
writer.apply()

multiMeshKernel = vls.ToMultiSurfaceMesh(minNodeDistFactor=0.02)
multiMeshKernel.insertNextLevelSet(substrateRK2)
multiMeshKernel.insertNextLevelSet(newLayerRK2)
multiMesh = vls.Mesh()
multiMeshKernel.setMesh(multiMesh)
multiMeshKernel.setSharpCorners(True)
multiMeshKernel.apply()
vls.VTKWriter(multiMesh, "multimesh_RK2.vtp").apply()

# RK3 Output
writer = vls.WriteVisualizationMesh()
writer.insertNextLevelSet(substrateRK)
writer.insertNextLevelSet(newLayerRK)
writer.addMetaData("time", passedTimeRK)
writer.setFileName("airgap_RK3")
writer.setExtractHullMesh(True)
writer.apply()

multiMeshKernel = vls.ToMultiSurfaceMesh(minNodeDistFactor=0.02)
multiMeshKernel.insertNextLevelSet(substrateRK)
multiMeshKernel.insertNextLevelSet(newLayerRK)
multiMesh = vls.Mesh()
multiMeshKernel.setMesh(multiMesh)
multiMeshKernel.setSharpCorners(True)
multiMeshKernel.apply()
vls.VTKWriter(multiMesh, "multimesh_RK3.vtp").apply()
