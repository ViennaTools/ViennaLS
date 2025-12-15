import viennals as vls
vls.setDimension(3)

# Define a constant velocity field
class ConstantVelocity(vls.VelocityField):
    def __init__(self, vel):
        super().__init__()
        self.velocity = vel

    def getScalarVelocity(self, coord, material, normal, pointId):
        return 0.0

    def getVectorVelocity(self, coord, material, normal, pointId):
        return self.velocity

def main():
    # Define grid and domain bounds
    gridDelta = 0.1
    bounds = [-5.0, 5.0, -5.0, 5.0, -5.0, 5.0]
    boundaryCons = [vls.BoundaryConditionEnum.INFINITE_BOUNDARY] * 3

    # Create initial level set (Sphere)
    sphere = vls.Domain(bounds, boundaryCons, gridDelta)
    origin = [0.0, 0.0, 0.0]
    radius = 1.5
    vls.MakeGeometry(sphere, vls.Sphere(origin, radius)).apply()

    # Create copies for Forward Euler and RK3
    sphereFE = vls.Domain(sphere)
    sphereRK3 = vls.Domain(sphere)

    # Define constant velocity field (moving in x-direction)
    vel = [1.0, 0.0, 0.0]
    velocityField = ConstantVelocity(vel)

    # Setup Advection: Forward Euler
    advectFE = vls.AdvectForwardEuler()
    advectFE.insertNextLevelSet(sphereFE)
    advectFE.setVelocityField(velocityField)
    advectFE.setAdvectionTime(2.0)
    advectFE.setIntegrationScheme(vls.IntegrationSchemeEnum.ENGQUIST_OSHER_1ST_ORDER)

    # Setup Advection: Runge-Kutta 3
    advectRK3 = vls.AdvectRungeKutta3()
    advectRK3.insertNextLevelSet(sphereRK3)
    advectRK3.setVelocityField(velocityField)
    advectRK3.setAdvectionTime(2.0)
    advectRK3.setIntegrationScheme(vls.IntegrationSchemeEnum.ENGQUIST_OSHER_1ST_ORDER)

    # Run Advection
    print("Running Forward Euler Advection...")
    advectFE.apply()
    
    checkFE = vls.Check(sphereFE)
    checkFE.apply()

    meshFE = vls.Mesh()
    vls.ToSurfaceMesh(sphereFE, meshFE).apply()
    vls.VTKWriter(meshFE, "sphereFE.vtp").apply()

    print("Running Runge-Kutta 3 Advection...")
    advectRK3.apply()
    
    checkRK3 = vls.Check(sphereRK3)
    checkRK3.apply()

    meshRK3 = vls.Mesh()
    vls.ToSurfaceMesh(sphereRK3, meshRK3).apply()
    vls.VTKWriter(meshRK3, "sphereRK3.vtp").apply()

if __name__ == "__main__":
    main()
