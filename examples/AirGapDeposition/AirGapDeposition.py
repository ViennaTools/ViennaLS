import viennals.d2 as vls
from viennals.common import VelocityField, Mesh, VTKWriter, BooleanOperationEnum

# @example AirGapDeposition.py
#  Example showing how to use the library for topography
#  simulation, by creating a trench geometry. A layer of a different material is
#  then grown directionally on top.


class velocityField(VelocityField):
    # coord and normalVec are lists with 3 elements
    # in 2D coord[2] and normalVec[2] are zero
    # getScalarVelocity must return a scalar
    def getScalarVelocity(self, coord, material, normal, pointId):
        return abs(normal[0]) + abs(normal[1])

    def getVectorVelocity(self, coord, material, normal, pointId):
        return (0.0, 0.0, 0.0)


extent = 30
gridDelta = 0.5

bounds = (-extent, extent, -extent, extent)
boundaryCons = (0, 1, 0)  # 0 = reflective, 1 = infinite, 2 = periodic

# create level set
substrate = vls.Domain(bounds, boundaryCons, gridDelta)

# create plane
origin = (0, 0, 0)
planeNormal = (0, 1, 0)

vls.MakeGeometry(substrate, vls.Plane(origin, planeNormal)).apply()

print("Extracting")
mesh = Mesh()
vls.ToSurfaceMesh(substrate, mesh).apply()
VTKWriter(mesh, "plane.vtp").apply()

# create layer used for booling
print("Creating box...")
trench = vls.Domain(bounds, boundaryCons, gridDelta)
minCorner = (-extent / 6.0, -25.0)
maxCorner = (extent / 6.0, 1.0)
vls.MakeGeometry(trench, vls.Box(minCorner, maxCorner)).apply()

print("Extracting")
vls.ToMesh(trench, mesh).apply()
VTKWriter(mesh, "box.vtp").apply()

# Create trench geometry
print("Booling trench")
vls.BooleanOperation(
    substrate, trench, BooleanOperationEnum.RELATIVE_COMPLEMENT
).apply()

# Now grow new material

# create new levelset for new material, which will be grown
# since it has to wrap around the substrate, just copy it
print("Creating new layer...")
newLayer = vls.Domain(substrate)

velocities = velocityField()

print("Advecting")
advectionKernel = vls.Advect()

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

    vls.ToSurfaceMesh(newLayer, mesh).apply()
    writer = VTKWriter(mesh, "trench{}.vtp".format(i))
    writer.addMetaData("time", passedTime)
    writer.apply()

print("Time passed during advection: {}".format(passedTime))
