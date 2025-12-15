import viennals as vls
import time

# Dimension
D = 3
vls.setNumThreads(8)

class EpitaxyVelocity(vls.VelocityField):
    def __init__(self, velocities):
        super().__init__()
        self.velocities = velocities
        self.R111 = 0.5
        self.R100 = 1.0
        # D is 3 in this context
        self.low = 0.5773502691896257 if D > 2 else 0.7071067811865476
        self.high = 1.0

    def getScalarVelocity(self, coord, material, normal, point_id):
        vel = max(abs(normal[0]), abs(normal[2]))
        factor = (self.R100 - self.R111) / (self.high - self.low)
        vel = (vel - self.low) * factor + self.R111

        if abs(normal[0]) < abs(normal[2]):
            vel *= 2.0

        mat_vel = 0.0
        if material < len(self.velocities):
            mat_vel = self.velocities[material]

        return vel * mat_vel

    def getVectorVelocity(self, coord, material, normal, point_id):
        return (0., 0., 0.)

def write_surface(domain, filename):
    mesh = vls.Mesh()
    vls.ToSurfaceMesh(domain, mesh).apply()
    vls.VTKWriter(mesh, filename).apply()

def main():
    # Set dimension to 3D
    vls.setDimension(D)

    # Parameters
    grid_delta = 0.03
    fin_width = 0.5
    fin_height = 0.2

    bounds = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0]
    boundary_conditions = [
        vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
        vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
        vls.BoundaryConditionEnum.INFINITE_BOUNDARY
    ]

    # Create Geometry
    mask = vls.Domain(bounds, boundary_conditions, grid_delta)
    plane = vls.Plane([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
    vls.MakeGeometry(mask, plane).apply()

    substrate = vls.Domain(bounds, boundary_conditions, grid_delta)
    min_point = [-fin_width / 2, -fin_width / 2, 0.0]
    max_point = [fin_width / 2, fin_width / 2, fin_height]
    box = vls.Box(min_point, max_point)
    vls.MakeGeometry(substrate, box).apply()

    vls.BooleanOperation(substrate, mask, vls.BooleanOperationEnum.UNION).apply()

    write_surface(mask, "mask.vtp")
    write_surface(substrate, "substrate.vtp")

    level_sets = [mask, substrate]
    vls.PrepareStencilLocalLaxFriedrichs(level_sets, [False, True])

    # Advection
    velocities = [0.0, -0.5]
    velocity_field = EpitaxyVelocity(velocities)

    advection = vls.Advect()
    for ls in level_sets:
        advection.insertNextLevelSet(ls)
    advection.setVelocityField(velocity_field)
    advection.setIntegrationScheme(vls.IntegrationSchemeEnum.STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER)
    advection.setAdvectionTime(0.5)

    start_time = time.time()
    advection.apply()
    end_time = time.time()
    print(f"Epitaxy took {end_time - start_time}s")

    vls.FinalizeStencilLocalLaxFriedrichs(level_sets)

    write_surface(substrate, "epitaxy.vtp")

if __name__ == "__main__":
    main()
