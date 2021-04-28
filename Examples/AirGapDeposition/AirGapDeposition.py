import viennaLS2d as vls

## @example AirGapDeposition.py
#  Example showing how to use the library for topography
#  simulation, by creating a trench geometry. A layer of a different material is
#  then grown directionally on top.

class velocityField(vls.lsVelocityField):
  # coord and normalVec are lists with 3 elements
  # in 2D coord[2] and normalVec[2] are zero
  # getScalarVelocity must return a scalar
  def getScalarVelocity(self, coord, material, normal,  pointId):
    return abs(normal[0]) + abs(normal[1])

  def getVectorVelocity(self, coord, material, normal,  pointId):
    return (0,0,0)

extent = 30
gridDelta = 0.5

bounds = (-extent, extent, -extent, extent)
boundaryCons = (0, 1, 0) # 0 = reflective, 1 = infinite, 2 = periodic

# create level set
substrate = vls.lsDomain(bounds, boundaryCons, gridDelta)

# create plane
origin = (0,0,0)
planeNormal = (0,1,0)

vls.lsMakeGeometry(substrate, vls.lsPlane(origin, planeNormal)).apply()

print("Extracting")
mesh = vls.lsMesh()
vls.lsToSurfaceMesh(substrate, mesh).apply()
vls.lsVTKWriter(mesh, "plane.vtk").apply()

# create layer used for booling
print("Creating box...")
trench = vls.lsDomain(bounds, boundaryCons, gridDelta)
minCorner = (-extent / 6., -25.)
maxCorner = (extent / 6., 1.)
vls.lsMakeGeometry(trench, vls.lsBox(minCorner, maxCorner)).apply()

print("Extracting")
vls.lsToMesh(trench, mesh).apply()
vls.lsVTKWriter(mesh, "box.vtk").apply()

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
advectionKernel.setIgnoreVoids(True)

# Now advect the level set 50 times, outputting every
# advection step. Save the physical time that
# passed during the advection.
passedTime = 0
numberOfSteps = 60
for i in range(numberOfSteps):
  advectionKernel.apply()
  passedTime += advectionKernel.getAdvectedTime()

  print("Advection step {} / {}".format(i, numberOfSteps))

  vls.lsToSurfaceMesh(newLayer, mesh).apply()
  vls.lsVTKWriter(mesh, "trench{}.vtk".format(i)).apply()

print("Time passed during advection: {}".format(passedTime))
