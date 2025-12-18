import viennals as vls

vls.setDimension(3)

# Implement own velocity field
class VelocityField(vls.VelocityField):
    def __init__(self):
        super().__init__()

    def getScalarVelocity(self, coord, material, normal, pointId):
        # Some arbitrary velocity function
        velocity = 1.0 + (2.3 if normal[0] > 0 else 0.5) * abs(normal[0] * normal[0])
        return velocity

    def getVectorVelocity(self, coord, material, normal, pointId):
        return [0.0, 0.0, 0.0]

def main():
    # Set number of threads
    vls.setNumThreads(8)

    gridDelta = 0.6999999

    discretizationSchemes = [
        vls.DiscretizationSchemeEnum.ENGQUIST_OSHER_1ST_ORDER,
        vls.DiscretizationSchemeEnum.ENGQUIST_OSHER_2ND_ORDER,
        vls.DiscretizationSchemeEnum.LAX_FRIEDRICHS_1ST_ORDER,
        vls.DiscretizationSchemeEnum.LAX_FRIEDRICHS_2ND_ORDER,
        vls.DiscretizationSchemeEnum.LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER,
        vls.DiscretizationSchemeEnum.LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER,
        vls.DiscretizationSchemeEnum.LOCAL_LAX_FRIEDRICHS_1ST_ORDER,
        vls.DiscretizationSchemeEnum.LOCAL_LAX_FRIEDRICHS_2ND_ORDER,
        vls.DiscretizationSchemeEnum.WENO_5TH_ORDER
    ]

    for scheme in discretizationSchemes:
        sphere1 = vls.Domain(gridDelta)

        origin = [5.0, 0.0, 0.0]
        radius = 7.3

        vls.MakeGeometry(sphere1, vls.Sphere(origin, radius)).apply()

        # Instantiate velocities
        velocities = VelocityField()

        advectionKernel = vls.AdvectRungeKutta3()
        advectionKernel.insertNextLevelSet(sphere1)
        advectionKernel.setVelocityField(velocities)
        advectionKernel.setDiscretizationScheme(scheme)
        advectionKernel.setSaveAdvectionVelocities(True)

        time = 0.0
        i = 0
        while time < 1.0 and i < 100:
            advectionKernel.apply()
            time += advectionKernel.getAdvectedTime()

            fileName = f"{scheme}_{i}.vtp"
            mesh = vls.Mesh()
            vls.ToMesh(sphere1, mesh).apply()
            vls.VTKWriter(mesh, "points_" + fileName).apply()
            vls.ToSurfaceMesh(sphere1, mesh).apply()
            vls.VTKWriter(mesh, "surface_" + fileName).apply()
            i += 1
        print(f"Done scheme {scheme}")

if __name__ == "__main__":
    main()
