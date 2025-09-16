import viennals as vls

# @example GeometricAdvection.py
#  3D Example showing how to use the library for topography
# emulation, by creating a trench geometry. A uniform
# layer of a different material is then grown on top. It is
# the same example as Deposition but emulates the deposition
# rather than simulating a slow growth.

vls.setDimension(3)

extent = 30
gridDelta = 0.5

bounds = (-extent, extent, -extent, extent, -extent, extent)
boundaryCons = (0, 0, 1)  # 0 = reflective, 1 = infinite, 2 = periodic

# create level set
substrate = vls.Domain(bounds, boundaryCons, gridDelta)

# create plane
origin = (0, 0, 0)
planeNormal = (0, 0, 1)

vls.MakeGeometry(substrate, vls.Plane(origin, planeNormal)).apply()

# create layer used for booling
print("Creating box...")
trench = vls.Domain(bounds, boundaryCons, gridDelta)
minCorner = (-extent - 1, -extent / 4.0, -15.0)
maxCorner = (extent + 1, extent / 4.0, 1.0)
vls.MakeGeometry(trench, vls.Box(minCorner, maxCorner)).apply()

# Create trench geometry
print("Booling trench")
vls.BooleanOperation(
    substrate, trench, vls.BooleanOperationEnum.RELATIVE_COMPLEMENT
).apply()

mesh = vls.Mesh()
vls.ToSurfaceMesh(substrate, mesh).apply()
vls.VTKWriter(mesh, "trench-initial.vtp").apply()

# Now grow new material

# create new levelset for new material, which will be grown
# since it has to wrap around the substrate, just copy it
print("Creating new layer...")
newLayer = vls.Domain(substrate)

print("Advecting")
# Advect the level set
dist = vls.SphereDistribution(4.0, gridDelta)
vls.GeometricAdvect(newLayer, dist).apply()

vls.ToSurfaceMesh(newLayer, mesh).apply()
vls.VTKWriter(mesh, "trench-final.vtp").apply()
