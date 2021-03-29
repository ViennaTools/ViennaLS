import viennaLS3d as vls

## @example Deposition.py
#  3D Example showing how to use the library for topography
#  simulation, by creating a trench geometry. A uniform
#  layer of a different material is then grown on top.

class velocityField(vls.lsVelocityField):
  # coord and normalVec are lists with 3 elements
  # in 2D coord[2] and normalVec[2] are zero
  # getScalarVelocity must return a scalar
  def getScalarVelocity(self, coord, material, normal,  pointId):
    # some arbitrary velocity function of your liking
    # (try changing it and see what happens :)
    velocity = 1
    return velocity

  def getVectorVelocity(self, coord, material, normal,  pointId):
    return (0,0,0)

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

# Now grow new material

# create new levelset for new material, which will be grown
# since it has to wrap around the substrate, just copy it
print("Creating new layer...")
newLayer = vls.lsDomain(substrate)

velocities = velocityField()

print("Advecting")
advectionKernel = vls.lsAdvect()

# the level set to be advected has to be inserted last
# the other could be taken as a mask layer for advection
advectionKernel.insertNextLevelSet(substrate)
advectionKernel.insertNextLevelSet(newLayer)

advectionKernel.setVelocityField(velocities)

# Advect the level set
counter = 1
passedTime = 0

mesh = vls.lsMesh()
while(passedTime < 4):
  advectionKernel.apply()
  passedTime += advectionKernel.getAdvectedTime()

  vls.lsToSurfaceMesh(newLayer, mesh).apply()
  vls.lsVTKWriter(mesh, "trench-{}.vtk".format(counter)).apply()

  vls.lsToMesh(newLayer, mesh).apply()
  vls.lsVTKWriter(mesh, "LS-{}.vtk".format(counter)).apply()

  counter = counter + 1

print("Time passed during advection: {}".format(passedTime))
