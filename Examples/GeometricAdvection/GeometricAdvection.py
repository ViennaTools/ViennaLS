import viennaLS3d as vls

## @example GeometricAdvection.py
#  3D Example showing how to use the library for topography
# emulation, by creating a trench geometry. A uniform
# layer of a different material is then grown on top. It is
# the same example as Deposition but emulates the deposition
# rather than simulating a slow growth.

extent = 30
gridDelta = 0.5

bounds = (-extent, extent, -extent, extent, -extent, extent)
boundaryCons = (0, 0, 1) # 0 = reflective, 1 = infinite, 2 = periodic

# create level set
substrate = vls.lsDomain(bounds, boundaryCons, gridDelta)

# create plane
origin = (0,0,0)
planeNormal = (0,0,1)

vls.lsMakeGeometry(substrate, vls.lsPlane(origin, planeNormal)).apply()

# create layer used for booling
print("Creating box...")
trench = vls.lsDomain(bounds, boundaryCons, gridDelta)
minCorner = (-extent - 1, -extent / 4., -15.)
maxCorner = (extent + 1, extent / 4., 1.)
vls.lsMakeGeometry(trench, vls.lsBox(minCorner, maxCorner)).apply()

# Create trench geometry
print("Booling trench")
vls.lsBooleanOperation(substrate, trench, vls.lsBooleanOperationEnum.RELATIVE_COMPLEMENT).apply()

mesh = vls.lsMesh()
vls.lsToSurfaceMesh(substrate, mesh).apply()
vls.lsVTKWriter(mesh, "trench-initial.vtk").apply()

# Now grow new material

# create new levelset for new material, which will be grown
# since it has to wrap around the substrate, just copy it
print("Creating new layer...")
newLayer = vls.lsDomain(substrate)

print("Advecting")
# Advect the level set
dist = vls.lsSphereDistribution(4.0, gridDelta)
vls.lsGeometricAdvect(newLayer, dist).apply()

vls.lsToSurfaceMesh(newLayer, mesh).apply()
vls.lsVTKWriter(mesh, "trench-final.vtk").apply()
