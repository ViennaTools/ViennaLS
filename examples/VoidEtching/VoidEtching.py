import viennals

viennals.setDimension(3)

# Implement own velocity field
class VelocityField(viennals.VelocityField):
    def getScalarVelocity(self, coord, material, normal, pointId):
        # isotropic etch rate
        return -1.0

    def getVectorVelocity(self, coord, material, normal, pointId):
        return [0.0, 0.0, 0.0]

def main():
    # Simulation parameters
    extent = 30.0
    gridDelta = 1.0
    
    # Bounds and boundary conditions
    bounds = [-extent, extent, -extent, extent, -extent, extent]
    boundaryCons = [
        viennals.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
        viennals.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
        viennals.BoundaryConditionEnum.INFINITE_BOUNDARY
    ]
    
    # Create substrate domain
    substrate = viennals.Domain(bounds, boundaryCons, gridDelta)

    
    # Create initial plane geometry
    origin = [0.0, 0.0, 0.0]
    planeNormal = [0.0, 0.0, 1.0]
    plane = viennals.Plane(origin, planeNormal)
    viennals.MakeGeometry(substrate, plane).apply()
    
    # Create spheres used for booling
    print("Creating spheres...")
    sphere = viennals.Domain(bounds, boundaryCons, gridDelta)
    
    # Define boolean operation (Relative Complement: substrate - sphere)
    boolOp = viennals.BooleanOperation(
        substrate, sphere, viennals.BooleanOperationEnum.RELATIVE_COMPLEMENT
    )
    
    # Sphere 1
    origin = [-12.0, -5.0, -15.0]
    radius = 10.0
    viennals.MakeGeometry(sphere, viennals.Sphere(origin, radius)).apply()
    boolOp.apply()
    
    # Sphere 2
    origin = [-7.0, -30.0, -20.0]
    radius = 8.0
    viennals.MakeGeometry(sphere, viennals.Sphere(origin, radius)).apply()
    boolOp.apply()
    
    # Sphere 3
    origin = [5.0, 15.0, -2.0]
    radius = 8.0
    viennals.MakeGeometry(sphere, viennals.Sphere(origin, radius)).apply()
    boolOp.apply()
    
    # Sphere 4
    origin = [2.0, 8.0, -27.0]
    radius = 8.0
    viennals.MakeGeometry(sphere, viennals.Sphere(origin, radius)).apply()
    boolOp.apply()
    
    # Now etch the substrate isotropically
    velocities = VelocityField()
    
    print("Advecting")
    
    advectionKernel = viennals.Advect()
    advectionKernel.insertNextLevelSet(substrate)
    advectionKernel.setVelocityField(velocities)
    advectionKernel.setIgnoreVoids(True)
    
    passedTime = 0.0
    numberOfSteps = 50
    
    for i in range(numberOfSteps):
        print(f"\rAdvection step {i} / {numberOfSteps}", end="", flush=True)
        
        mesh = viennals.Mesh()
        viennals.ToSurfaceMesh(substrate, mesh).apply()
        viennals.VTKWriter(mesh, f"void-{i}.vtp").apply()
        
        advectionKernel.apply()
        passedTime += advectionKernel.getAdvectedTime()
        
    print()
    print(f"Time passed during advection: {passedTime}")

if __name__ == "__main__":
    main()
