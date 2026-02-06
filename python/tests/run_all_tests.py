import unittest
import os
import viennals as vls
import viennals.d2 as d2
import viennals.d3 as d3

class TestVelocityField(vls.VelocityField):
    def __init__(self):
        super().__init__()

    def getScalarVelocity(self, coord, material, normal, pointId):
        return 1.0

    def getVectorVelocity(self, coord, material, normal, pointId):
        return [0.0, 0.0, 0.0]

class WrappingTest(unittest.TestCase):
    def test_enums_and_constants(self):
        print("\n[TEST] test_enums_and_constants")
        # Verify Enums exist
        self.assertTrue(hasattr(vls, "SpatialSchemeEnum"))
        self.assertTrue(hasattr(vls, "TemporalSchemeEnum"))
        self.assertTrue(hasattr(vls, "BooleanOperationEnum"))
        self.assertTrue(hasattr(vls, "CurvatureEnum"))
        self.assertTrue(hasattr(vls, "FeatureDetectionEnum"))
        self.assertTrue(hasattr(vls, "BoundaryConditionEnum"))
        self.assertTrue(hasattr(vls, "FileFormatEnum"))
        self.assertTrue(hasattr(vls, "VoidTopSurfaceEnum"))
        self.assertTrue(hasattr(vls, "TransformEnum"))
        
        # Verify Logger
        logger = vls.Logger.getInstance()
        logger.setLogLevel(vls.LogLevel.INFO)
        logger.addInfo("Testing Logger")
        print("  > Done test_enums_and_constants")

    def test_mesh_and_data(self):
        print("\n[TEST] test_mesh_and_data")
        # PointData
        print("  > Testing PointData scalar insertion")
        pd = vls.PointData()
        pd.insertNextScalarData([1.0, 2.0], "Scalars")
        self.assertEqual(pd.getScalarDataSize(), 1)
        
        # Mesh
        print("  > Testing Mesh insertion")
        mesh = vls.Mesh()
        mesh.insertNextNode([0.0, 0.0, 0.0])
        mesh.insertNextNode([1.0, 0.0, 0.0])
        mesh.insertNextNode([0.0, 1.0, 0.0])
        mesh.insertNextTriangle([0, 1, 2])
        self.assertEqual(len(mesh.getNodes()), 3)
        self.assertEqual(len(mesh.getTriangles()), 1)
        
        # MaterialMap
        print("  > Testing MaterialMap")
        mm = vls.MaterialMap()
        mm.insertNextMaterial(0)
        self.assertEqual(mm.getNumberOfMaterials(), 1)
        
        # Vector Data
        print("  > Testing PointData vector insertion")
        pd.insertNextVectorData([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], "Vectors")
        self.assertEqual(pd.getVectorDataSize(), 1)
        
        print("  > Done test_mesh_and_data")

    def test_algorithms(self):
        print("\n[TEST] test_algorithms (2D & 3D)")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            gridDelta = 0.2
            bounds = [-10.0, 10.0] * dim
            dom = module.Domain(bounds, [vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY] * dim, gridDelta)
            
            # Geometries
            origin = [0.0] * dim
            radius = 5.0
            print("    >> MakeGeometry")
            sphere_geom = module.Sphere(origin, radius)
            module.MakeGeometry(dom, sphere_geom).apply()
            self.assertGreater(dom.getNumberOfPoints(), 0)
            
            print("    >> Advect")
            # Advect
            vel = TestVelocityField()
            adv = module.Advect()
            adv.insertNextLevelSet(dom)
            adv.setVelocityField(vel)
            adv.setAdvectionTime(0.1)
            adv.setSpatialScheme(vls.SpatialSchemeEnum.ENGQUIST_OSHER_1ST_ORDER)
            adv.setTemporalScheme(vls.TemporalSchemeEnum.FORWARD_EULER)
            adv.apply()
            
            print("    >> GeometricAdvect")
            # Geometric Advect & Distributions
            dist = module.SphereDistribution(radius)
            gadv = module.GeometricAdvect(dom, dist)
            gadv.apply()
            
            # Check other distributions instantiation
            print("    >> Checking Distribution instantiation")
            _ = module.BoxDistribution([1.0] * 3) # BoxDistribution always takes 3D bounds in constructor? No, usually dim.
            # Checking pyWrap.hpp: BoxDistribution<T, D> takes std::array<T, 3> in constructor?
            # py::class_<BoxDistribution<T, D> ... .def(py::init(... std::array<T, 3>))
            # It seems BoxDistribution constructor in Python bindings expects 3 args regardless of D?
            # Let's check pyWrap.hpp again.
            # .def(py::init(&SmartPointer<BoxDistribution<T, D>>::template New<const std::array<T, 3>>))
            # Yes, it seems hardcoded to 3 in the binding signature for BoxDistribution.
            _ = module.BoxDistribution([1.0, 1.0, 1.0])
            _ = module.CustomSphereDistribution([1.0, 2.0])
            
            print("    >> Boolean")
            # Boolean
            dom2 = module.Domain(dom)
            bool_op = module.BooleanOperation(dom, dom2, vls.BooleanOperationEnum.UNION)
            bool_op.apply()
            
            print("    >> CompareVolume")
            # CompareVolume
            dom_vol1 = module.Domain(gridDelta)
            module.MakeGeometry(dom_vol1, module.Sphere([0.0] * dim, 5.0)).apply()
            dom_vol2 = module.Domain(gridDelta)
            module.MakeGeometry(dom_vol2, module.Sphere([0.0] * dim, 6.0)).apply()
            cv = module.CompareVolume(dom_vol1, dom_vol2)
            cv.apply()
            self.assertGreater(cv.getVolumeMismatch(), 0.0)

            print("    >> Calc Curvatures")
            # Calc Curvatures
            curv = module.CalculateCurvatures(dom)
            curv.apply()
            
            print("    >> Calc Normals")
            # Calc Normals
            norms = module.CalculateNormalVectors(dom)
            norms.apply()
            
            print("    >> Calc Visibilities")
            # Calc Visibilities
            # CalculateVisibilities expects Vec3D (3 components) even in 2D
            module.CalculateVisibilities(dom, [0.0, 0.0, 1.0], "Visibility").apply()
            
            print("    >> Check")
            # Check
            module.Check(dom).apply()
            
            print("    >> Prune")
            # Prune
            module.Prune(dom).apply()
            
            print("    >> Expand")
            # Expand
            module.Expand(dom, 1).apply()
            
            print("    >> Reduce")
            # Reduce
            module.Reduce(dom, 1).apply()
            
            print("    >> MarkVoidPoints")
            # Mark Void Points
            module.MarkVoidPoints(dom).apply()
            
            # Expand to ensure sufficient width for feature detection (curvature)
            module.Prune(dom).apply()
            module.Expand(dom, 5).apply()

            print("    >> DetectFeatures")
            # Detect Features
            module.DetectFeatures(dom).apply()
            
            print("    >> Mesh Conversions")
            # Mesh Conversions
            # Expand back to a safe width for meshing
            module.Expand(dom, 3).apply()

            mesh = vls.Mesh()
            print("      >>> ToMesh")
            module.ToMesh(dom, mesh).apply()
            print("      >>> ToSurfaceMesh")
            module.ToSurfaceMesh(dom, mesh).apply()
            print("      >>> ToDiskMesh")
            module.ToDiskMesh(dom, mesh).apply()
            
            print("    >> IO")
            # IO
            fname = f"test_{dim}d.vtp"
            vls.VTKWriter(mesh, fname).apply()
            if os.path.exists(fname):
                os.remove(fname)
            
            print("    >> RemoveStrayPoints")
            # Remove Stray Points
            module.RemoveStrayPoints(dom).apply()

        print("  > Done test_algorithms")

    def test_geometries(self):
        print("\n[TEST] test_geometries")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            # Initialize domain with bounds to ensure MakeGeometry(Plane) works efficiently
            bounds = [-10.0, 10.0] * dim
            boundaryCons = [vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY] * dim
            # Make one boundary infinite for variety
            boundaryCons[-1] = vls.BoundaryConditionEnum.INFINITE_BOUNDARY
            dom = module.Domain(bounds, boundaryCons, 0.2)
            
            # Plane
            print("    >> Plane")
            normal = [0.0] * dim
            normal[-1] = 1.0
            plane = module.Plane([0.0] * dim, normal)
            module.MakeGeometry(dom, plane).apply()
            self.assertGreater(dom.getNumberOfPoints(), 0)
            
            # Box
            print("    >> Box")
            box = module.Box([-1.0] * dim, [1.0] * dim)
            module.MakeGeometry(dom, box).apply()
            
            # Cylinder
            print("    >> Cylinder")
            # Cylinder constructor takes 3D vectors for origin and axis even in 2D?
            # pyWrap.hpp: .def(py::init(&SmartPointer<Cylinder<T, D>>::template New<const std::vector<T> & ...>))
            # It uses std::vector, so it should adapt to D.
            axis = [0.0] * dim
            axis[-1] = 1.0
            cyl = module.Cylinder([0.0] * dim, axis, 2.0, 1.0, 1.0)
            module.MakeGeometry(dom, cyl).apply()

        print("  > Done test_geometries")

    def test_io_and_transforms(self):
        print("\n[TEST] test_io_and_transforms")
        # Create a mesh
        mesh = vls.Mesh()
        mesh.insertNextNode([0.0, 0.0, 0.0])
        mesh.insertNextNode([1.0, 0.0, 0.0])
        mesh.insertNextNode([0.0, 1.0, 0.0])
        mesh.insertNextTriangle([0, 1, 2])
        
        # Transform
        print("  > TransformMesh")
        trans = vls.TransformMesh(mesh, vls.TransformEnum.TRANSLATION, [1.0, 0.0, 0.0], 0.0)
        trans.apply()
        nodes = mesh.getNodes()
        self.assertEqual(nodes[0][0], 1.0)
        
        # VTK Writer
        print("  > VTKWriter")
        vls.VTKWriter(mesh, "test_io.vtp").apply()
        
        # VTK Reader
        print("  > VTKReader")
        mesh2 = vls.Mesh()
        vls.VTKReader(mesh2, vls.FileFormatEnum.VTP, "test_io.vtp").apply()
        self.assertEqual(len(mesh2.getNodes()), 3)
        
        if os.path.exists("test_io.vtp"):
            os.remove("test_io.vtp")
        print("  > Done test_io_and_transforms")

    def test_mesh_conversion_roundtrip(self):
        print("\n[TEST] test_mesh_conversion_roundtrip")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            dom = module.Domain(0.2)
            module.MakeGeometry(dom, module.Sphere([0.0]*dim, 5.0)).apply()
            
            mesh = vls.Mesh()
            module.ToSurfaceMesh(dom, mesh).apply()
            
            dom2 = module.Domain(0.2)
            module.FromSurfaceMesh(dom2, mesh).apply()
            self.assertGreater(dom2.getNumberOfPoints(), 0)
        print("  > Done test_mesh_conversion_roundtrip")

    def test_native_io(self):
        print("\n[TEST] test_native_io")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            dom = module.Domain(0.2)
            module.MakeGeometry(dom, module.Sphere([0.0]*dim, 5.0)).apply()
            
            filename = f"test_native_{dim}d.lvst"
            print("    >> Writer")
            module.Writer(dom, filename).apply()
            
            dom_in = module.Domain(0.2)
            print("    >> Reader")
            module.Reader(dom_in, filename).apply()
            self.assertGreater(dom_in.getNumberOfPoints(), 0)
            
            if os.path.exists(filename):
                os.remove(filename)
        print("  > Done test_native_io")

    def test_advanced_meshing(self):
        print("\n[TEST] test_advanced_meshing")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            dom1 = module.Domain(0.2)
            module.MakeGeometry(dom1, module.Sphere([0.0]*dim, 5.0)).apply()
            
            dom2 = module.Domain(0.2)
            center = [0.0]*dim
            center[0] = 2.0
            module.MakeGeometry(dom2, module.Sphere(center, 5.0)).apply()
            
            # ToMultiSurfaceMesh
            print("    >> ToMultiSurfaceMesh")
            mesh_multi = vls.Mesh()
            module.ToMultiSurfaceMesh([dom1, dom2], mesh_multi).apply()
            self.assertGreater(len(mesh_multi.getNodes()), 0)
            
            # ToVoxelMesh
            print("    >> ToVoxelMesh")
            mesh_voxel = vls.Mesh()
            module.ToVoxelMesh([dom1, dom2], mesh_voxel).apply()
            self.assertGreater(len(mesh_voxel.getNodes()), 0)
        print("  > Done test_advanced_meshing")

    def test_stencil_functions(self):
        print("\n[TEST] test_stencil_functions")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            dom = module.Domain(0.2)
            module.MakeGeometry(dom, module.Sphere([0.0]*dim, 5.0)).apply()
            module.Expand(dom, 5).apply()

            # Create a second domain to form a valid multi-layer system
            dom2 = module.Domain(0.2)
            module.MakeGeometry(dom2, module.Sphere([0.0]*dim, 10.0)).apply()
            module.Expand(dom2, 5).apply()

            # PrepareStencilLocalLaxFriedrichs
            print("    >> PrepareStencilLocalLaxFriedrichs")
            module.PrepareStencilLocalLaxFriedrichs([dom, dom2], [False, True])
            
            # Finalize
            print("    >> FinalizeStencilLocalLaxFriedrichs")
            module.FinalizeStencilLocalLaxFriedrichs([dom, dom2])
        print("  > Done test_stencil_functions")

    def test_geometric_advect_box(self):
        print("\n[TEST] test_geometric_advect_box")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            dom = module.Domain(0.2)
            module.MakeGeometry(dom, module.Sphere([0.0]*dim, 5.0)).apply()
            
            print("    >> BoxDistribution")
            # BoxDistribution constructor in Python bindings expects 3 args regardless of D
            dist = module.BoxDistribution([1.0, 1.0, 1.0])
            
            print("    >> GeometricAdvect")
            module.GeometricAdvect(dom, dist).apply()
            self.assertGreater(dom.getNumberOfPoints(), 0)
        print("  > Done test_geometric_advect_box")

    def test_from_mesh_volume(self):
        print("\n[TEST] test_from_mesh_volume")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            mesh = vls.Mesh()
            if dim == 3:
                # Create a simple tetrahedral mesh
                mesh.insertNextNode([0.0, 0.0, 0.0])
                mesh.insertNextNode([5.0, 0.0, 0.0])
                mesh.insertNextNode([0.0, 5.0, 0.0])
                mesh.insertNextNode([0.0, 0.0, 5.0])
                mesh.insertNextTetra([0, 1, 2, 3])
                ls_values = [-1.0, 1.0, 1.0, 1.0]
            else:
                # Create a simple triangle mesh for 2D volume
                mesh.insertNextNode([0.0, 0.0, 0.0])
                mesh.insertNextNode([5.0, 0.0, 0.0])
                mesh.insertNextNode([0.0, 5.0, 0.0])
                mesh.insertNextTriangle([0, 1, 2])
                ls_values = [-1.0, 1.0, 1.0]
            
            # FromMesh requires "LSValues" in cellData corresponding to nodes
            mesh.getPointData().insertNextScalarData(ls_values, "LSValues")

            dom = module.Domain(0.2)
            # FromMesh should handle volume meshes if implemented in the wrapper logic
            print("    >> FromMesh (Volume)")
            module.FromMesh(dom, mesh).apply()
            self.assertGreater(dom.getNumberOfPoints(), 0)
        print("  > Done test_from_mesh_volume")

    def test_boolean_variants(self):
        print("\n[TEST] test_boolean_variants")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            dom1 = module.Domain(0.2)
            module.MakeGeometry(dom1, module.Sphere([0.0]*dim, 5.0)).apply()
            
            dom2 = module.Domain(0.2)
            center = [0.0]*dim
            center[0] = 2.0
            module.MakeGeometry(dom2, module.Sphere(center, 5.0)).apply()
            
            # Intersect
            print("    >> Intersect")
            dom_int = module.Domain(dom1)
            module.BooleanOperation(dom_int, dom2, vls.BooleanOperationEnum.INTERSECT).apply()
            self.assertGreater(dom_int.getNumberOfPoints(), 0)
            
            # Relative Complement (Difference)
            print("    >> Relative Complement")
            dom_diff = module.Domain(dom1)
            module.BooleanOperation(dom_diff, dom2, vls.BooleanOperationEnum.RELATIVE_COMPLEMENT).apply()
            self.assertGreater(dom_diff.getNumberOfPoints(), 0)
            
            # Invert
            print("    >> Invert")
            dom_inv = module.Domain(dom1)
            module.BooleanOperation(dom_inv, vls.BooleanOperationEnum.INVERT).apply()
            # Inverted domain is infinite, but points should exist
            self.assertGreater(dom_inv.getNumberOfPoints(), 0)
        print("  > Done test_boolean_variants")

    def test_domain_manipulation(self):
        print("\n[TEST] test_domain_manipulation")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            dom = module.Domain(0.2)
            module.MakeGeometry(dom, module.Sphere([0.0]*dim, 5.0)).apply()
            
            # Deep Copy
            print("    >> deepCopy")
            dom_copy = module.Domain(0.2)
            dom_copy.deepCopy(dom)
            self.assertEqual(dom.getNumberOfPoints(), dom_copy.getNumberOfPoints())
            
            # Level Set Width
            print("    >> setLevelSetWidth")
            dom.setLevelSetWidth(3)
            self.assertEqual(dom.getLevelSetWidth(), 3)
        print("  > Done test_domain_manipulation")

    def test_vtk_metadata(self):
        print("\n[TEST] test_vtk_metadata")
        mesh = vls.Mesh()
        mesh.insertNextNode([0.0, 0.0, 0.0])
        
        writer = vls.VTKWriter(mesh, "test_meta.vtp")
        writer.addMetaData("Value", 42.0)
        writer.apply()
        
        if os.path.exists("test_meta.vtp"):
            os.remove("test_meta.vtp")
        print("  > Done test_vtk_metadata")

    def test_write_visualization_mesh(self):
        print("\n[TEST] test_write_visualization_mesh")
        if hasattr(d3, "WriteVisualizationMesh"):
            dom = d3.Domain(0.2)
            d3.MakeGeometry(dom, d3.Sphere([0.0, 0.0, 0.0], 5.0)).apply()
            
            wvm = d3.WriteVisualizationMesh(dom)
            wvm.setFileName("test_vis")
            wvm.apply()
        else:
            print("  > WriteVisualizationMesh not available (VIENNALS_USE_VTK off?)")
        print("  > Done test_write_visualization_mesh")

    def test_comparisons(self):
        print("\n[TEST] test_comparisons")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            gridDelta = 0.5
            dom = module.Domain(gridDelta)
            
            # Make Geometry
            module.MakeGeometry(dom, module.Sphere([0.0]*dim, 5.0)).apply()
            
            dom2 = module.Domain(dom)
            # Shift dom2 slightly
            center = [0.0]*dim
            center[0] = 0.1
            module.MakeGeometry(dom2, module.Sphere(center, 5.0)).apply()
            
            print("    >> Comparisons")
            # Comparisons
            module.CompareVolume(dom, dom2).apply()
            module.CompareChamfer(dom, dom2).apply()
            
            comp_crit = module.CompareCriticalDimensions(dom, dom2)
            comp_crit.addXRange(-1.0, 1.0, True)
            comp_crit.apply()
            
            module.CompareNarrowBand(dom, dom2).apply()
            
            print("    >> Sparse Field")
            # Sparse Field (needs expanded/reduced)
            dom_exp = module.Domain(dom)
            module.Expand(dom_exp, 55).apply()
            module.CompareSparseField(dom_exp, dom).apply()
        print("  > Done test_comparisons")

    def test_cross_dimension(self):
        print("\n[TEST] test_cross_dimension")
        # Slice 3D -> 2D
        print("  > Creating 3D domain")
        dom3 = d3.Domain(0.2)
        d3.MakeGeometry(dom3, d3.Sphere([0.0, 0.0, 0.0], 5.0)).apply()
        
        dom2 = d2.Domain(0.2)
        
        # Slice at z=0
        print("  > Slice")
        s = vls.Slice(dom3, dom2, 2, 0.0)
        s.apply()
        self.assertGreater(dom2.getNumberOfPoints(), 0)
        
        # Extrude 2D -> 3D
        print("  > Extrude")
        dom3_out = d3.Domain(0.2)
        extent = [-1.0, 1.0]
        bcs = [vls.BoundaryConditionEnum.INFINITE_BOUNDARY] * 3
        
        e = vls.Extrude(dom2, dom3_out, extent, 2, bcs)
        e.apply()
        self.assertGreater(dom3_out.getNumberOfPoints(), 0)
        print("  > Done test_cross_dimension")

    def test_point_cloud_convex_hull(self):
        print("\n[TEST] test_point_cloud_convex_hull")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            # PointCloud
            print("    >> Creating PointCloud")
            # Create points on axes
            points = [[0.0]*dim] + [[1.0 if i == j else 0.0 for i in range(dim)] for j in range(dim)]
            pc = module.PointCloud(points)
            
            # ConvexHull
            print("    >> ConvexHull")
            mesh = vls.Mesh()
            hull = module.ConvexHull(mesh, pc)
            hull.apply()
            self.assertGreater(len(mesh.getNodes()), 0)
        print("  > Done test_point_cloud_convex_hull")

    def test_advection_callback(self):
        print("\n[TEST] test_advection_callback")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            print("    >> Setting up advection")
            dom = module.Domain(0.2)
            module.MakeGeometry(dom, module.Sphere([0.0]*dim, 5.0)).apply()
            
            vel = TestVelocityField()
            adv = module.Advect()
            adv.insertNextLevelSet(dom)
            adv.setVelocityField(vel)
            adv.setAdvectionTime(0.1)
            adv.setTemporalScheme(vls.TemporalSchemeEnum.RUNGE_KUTTA_2ND_ORDER)
            
            # Use a mutable container to capture callback execution
            callback_data = {'called': False}
            def callback(d):
                callback_data['called'] = True
                return True
                
            print("    >> Running advection with callback")
            adv.setVelocityUpdateCallback(callback)
            adv.apply()
            
            self.assertTrue(callback_data['called'])
        print("  > Done test_advection_callback")

    def test_make_geometry_boundary(self):
        print("\n[TEST] test_make_geometry_boundary")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            print("    >> Creating domain with bounds")
            bounds = [-5.0, 5.0] * dim
            boundaryCons = [vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY] * dim
            dom = module.Domain(bounds, boundaryCons, 1.0)
            
            print("    >> MakeGeometry with IgnoreBoundaryConditions")
            maker = module.MakeGeometry(dom, module.Sphere([0.0]*dim, 4.0))
            maker.setIgnoreBoundaryConditions(True)
            maker.apply()
            
            self.assertGreater(dom.getNumberOfPoints(), 0)
        print("  > Done test_make_geometry_boundary")

    def test_mesh_operations(self):
        print("\n[TEST] test_mesh_operations")
        print("  > Inserting duplicate nodes")
        mesh = vls.Mesh()
        # Insert duplicate nodes
        mesh.insertNextNode([0.0, 0.0, 0.0])
        mesh.insertNextNode([0.0, 0.0, 0.0])
        self.assertEqual(len(mesh.getNodes()), 2)
        
        print("  > removeDuplicateNodes")
        mesh.removeDuplicateNodes()
        self.assertEqual(len(mesh.getNodes()), 1)
        
        # Append
        print("  > append")
        mesh2 = vls.Mesh()
        mesh2.insertNextNode([1.0, 1.0, 1.0])
        
        mesh.append(mesh2)
        self.assertEqual(len(mesh.getNodes()), 2)
        print("  > Done test_mesh_operations")

    def test_advect_configuration(self):
        print("\n[TEST] test_advect_configuration")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            print("    >> Creating domain")
            dom = module.Domain(0.2)
            module.MakeGeometry(dom, module.Sphere([0.0]*dim, 5.0)).apply()
            
            vel = TestVelocityField()
            adv = module.Advect()
            adv.insertNextLevelSet(dom)
            adv.setVelocityField(vel)
            
            # Test Single Step
            print("    >> Testing Single Step")
            adv.setAdvectionTime(10.0) # Large time that would normally require multiple steps
            adv.setSingleStep(True)
            adv.apply()
            self.assertEqual(adv.getNumberOfTimeSteps(), 1)
        print("  > Done test_advect_configuration")

    def test_mark_void_points_advanced(self):
        print("\n[TEST] test_mark_void_points_advanced")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            print("    >> Creating disjoint spheres")
            dom = module.Domain(0.2)
            # Create two disjoint spheres
            module.MakeGeometry(dom, module.Sphere([0.0]*dim, 3.0)).apply()
            
            dom2 = module.Domain(0.2)
            center = [0.0]*dim
            center[0] = 10.0
            module.MakeGeometry(dom2, module.Sphere(center, 3.0)).apply()
            
            # Union them to create a domain with two separate components
            print("    >> Union")
            module.BooleanOperation(dom, dom2, vls.BooleanOperationEnum.UNION).apply()
            
            print("    >> MarkVoidPoints")
            marker = module.MarkVoidPoints(dom)
            marker.setSaveComponentIds(True)
            marker.apply()
            
            # Should find 3 connected components (2 spheres + 1 background)
            self.assertEqual(marker.getNumberOfComponents(), 3)
        print("  > Done test_mark_void_points_advanced")

    def test_adaptive_time_stepping(self):
        print("\n[TEST] test_adaptive_time_stepping")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            print("    >> Setting up advection")
            dom = module.Domain(0.2)
            module.MakeGeometry(dom, module.Sphere([0.0]*dim, 5.0)).apply()
            
            vel = TestVelocityField()
            adv = module.Advect()
            adv.insertNextLevelSet(dom)
            adv.setVelocityField(vel)
            adv.setAdvectionTime(0.1)
            
            # Enable adaptive time stepping
            print("    >> Enabling Adaptive Time Stepping")
            adv.setAdaptiveTimeStepping(True, 10)
            adv.apply()
            
            self.assertGreater(dom.getNumberOfPoints(), 0)
        print("  > Done test_adaptive_time_stepping")

    def test_geometric_distribution_functions(self):
        print("\n[TEST] test_geometric_distribution_functions")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            print("    >> SphereDistribution")
            # Sphere Distribution
            dist = module.SphereDistribution(5.0)
            # Point inside (0,0,0) relative to center (0,0,0)
            print("    >> Checking isInside")
            # isInside takes 3D points in bindings?
            # pyWrap.hpp: .def("isInside", &GeometricAdvectDistribution<T, D>::isInside)
            # GeometricAdvectDistribution uses std::array<T, 3> in pyWrap.hpp for PylsGeometricAdvectDistribution trampoline
            # but the base class is templated on D.
            # However, the trampoline class PylsGeometricAdvectDistribution hardcodes vectorType to std::array<CoordType, 3>.
            # This suggests that even for 2D, the distribution functions might expect 3D coordinates in the Python binding layer if they go through the trampoline.
            # But SphereDistribution is a direct binding.
            # Let's assume it takes 3D points because of the trampoline definition in pyWrap.hpp which seems to force 3D for the virtual methods.
            # Wait, SphereDistribution binding doesn't use the trampoline unless it inherits from it in Python, which it doesn't here.
            # But the C++ class SphereDistribution<T, D> uses VectorType<T, D>.
            # Let's try passing 3D points to be safe, as extra dimensions are usually ignored in 2D logic if passed as std::array<T, 3> to a function expecting std::array<T, 2>? No, pybind11 would complain.
            # Actually, looking at pyWrap.hpp, PylsGeometricAdvectDistribution is defined with `typedef std::array<viennahrle::CoordType, 3> vectorType;`
            # This looks like a bug/limitation in the bindings for 2D if it enforces 3D.
            # However, let's try using `[0.0]*3` which is safe for 3D and might work for 2D if the binding expects 3 args.
            # If the binding expects 2 args for 2D, `[0.0]*3` will fail.
            # Let's check `GeometricAdvectDistribution` binding in `pyWrap.hpp`.
            # It uses `PylsGeometricAdvectDistribution<D>` as trampoline.
            # `PylsGeometricAdvectDistribution` has `typedef std::array<viennahrle::CoordType, 3> vectorType;` regardless of D.
            # This strongly suggests that the Python side expects 3-element lists/tuples even for 2D distributions.
            
            pt_zero = [0.0] * 3
            pt_out = [0.0] * 3
            pt_out[0] = 6.0
            
            self.assertTrue(dist.isInside(pt_zero, pt_zero, 1e-6))
            self.assertFalse(dist.isInside(pt_zero, pt_out, 1e-6))
            
            # Signed Distance
            print("    >> Checking getSignedDistance")
            dist_val = dist.getSignedDistance(pt_zero, pt_zero, 0)
            self.assertAlmostEqual(dist_val, -5.0, delta=1e-5)
        print("  > Done test_geometric_distribution_functions")

    def test_advect_additional_options(self):
        print("\n[TEST] test_advect_additional_options")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            print("    >> Setting up advection")
            dom = module.Domain(0.2)
            module.MakeGeometry(dom, module.Sphere([0.0]*dim, 5.0)).apply()
            
            vel = TestVelocityField()
            adv = module.Advect()
            adv.insertNextLevelSet(dom)
            adv.setVelocityField(vel)
            adv.setAdvectionTime(0.1)
            
            print("    >> Setting IgnoreVoids")
            # Test Ignore Voids
            adv.setIgnoreVoids(True)
            
            print("    >> Setting DissipationAlpha")
            # Test Dissipation Alpha
            adv.setDissipationAlpha(0.5)
            
            print("    >> Setting CheckDissipation")
            # Test Check Dissipation
            adv.setCheckDissipation(False)
            
            adv.apply()
            
            self.assertGreater(dom.getNumberOfPoints(), 0)
        print("  > Done test_advect_additional_options")

    def test_additional_options(self):
        print("\n[TEST] test_additional_options")
        for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")
            print("    >> Creating domain with bounds")
            bounds = [-10.0, 10.0] * dim
            boundaryCons = [vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY] * dim
            dom = module.Domain(bounds, boundaryCons, 0.2)
            module.MakeGeometry(dom, module.Sphere([0.0]*dim, 5.0)).apply()
            
            # ToMesh options
            print("    >> ToMesh options (OnlyDefined, OnlyActive)")
            mesh = vls.Mesh()
            tomesh = module.ToMesh(dom, mesh)
            tomesh.setOnlyDefined(True)
            tomesh.setOnlyActive(True)
            tomesh.apply()
            self.assertGreater(len(mesh.getNodes()), 0)
            
            # DetectFeatures options
            print("    >> DetectFeatures options (Threshold, Method)")
            df = module.DetectFeatures(dom)
            df.setDetectionThreshold(10.0)
            df.setDetectionMethod(vls.FeatureDetectionEnum.NORMALS_ANGLE)
            df.apply()
            
            # MarkVoidPoints options
            print("    >> MarkVoidPoints options (Reverse, LargestSurface)")
            mvp = module.MarkVoidPoints(dom)
            mvp.setReverseVoidDetection(True)
            mvp.setDetectLargestSurface(True)
            mvp.apply()
            
            # RemoveStrayPoints options
            print("    >> RemoveStrayPoints options (VoidTopSurface)")
            rsp = module.RemoveStrayPoints(dom)
            rsp.setVoidTopSurface(vls.VoidTopSurfaceEnum.LARGEST)
            rsp.apply()
            
        print("  > Done test_additional_options")

    def test_corner_preservation(self):
        print("\n[TEST] test_corner_preservation")
        gridDelta = 0.0485
        tolerance = 0.03 * gridDelta

        for module, dim in [(d2, 2)]:
        # for module, dim in [(d2, 2), (d3, 3)]:
            print(f"  > Testing {dim}D")

            # --- Box Test ---
            print("    >> Box Test")
            bounds = [-3.0, 3.0] * dim
            boundaryCons = [vls.BoundaryConditionEnum.INFINITE_BOUNDARY] * dim
            dom = module.Domain(bounds, boundaryCons, gridDelta)

            min_corner = [-1.0] * dim
            max_corner = [1.0] * dim
            box = module.Box(min_corner, max_corner)
            module.MakeGeometry(dom, box).apply()

            mesh = vls.Mesh()
            tsm = module.ToSurfaceMesh(dom, mesh)
            tsm.setSharpCorners(True)
            tsm.apply()

            nodes = mesh.getNodes()
            # Generate expected corners
            expected_corners = []
            num_corners = 1 << dim
            for i in range(num_corners):
                corner = []
                for j in range(dim):
                    val = 1.0 if ((i >> j) & 1) else -1.0
                    corner.append(val)
                expected_corners.append(corner)

            # Verify
            for ex in expected_corners:
                found = False
                for n in nodes:
                    dist_sq = sum([(n[k] - ex[k]) ** 2 for k in range(dim)])
                    if dist_sq < tolerance:
                        found = True
                        break
                if not found:
                    self.fail(f"Box Corner {ex} not found in {dim}D mesh")

            # --- Box Cavity Test ---
            print("    >> Box Cavity Test")
            # Boundary conditions: Reflective for D-1, Infinite for last
            bc_cavity = [vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY] * (dim - 1)
            bc_cavity.append(vls.BoundaryConditionEnum.INFINITE_BOUNDARY)

            substrate = module.Domain(bounds, bc_cavity, gridDelta)

            origin = [0.0] * dim
            origin[dim - 1] = 0.025
            normal = [0.0] * dim
            normal[dim - 1] = 1.0
            plane = module.Plane(origin, normal)
            module.MakeGeometry(substrate, plane).apply()

            box_domain = module.Domain(bounds, bc_cavity, gridDelta)
            module.MakeGeometry(box_domain, box).apply()  # Reuse box geometry (-1 to 1)

            module.BooleanOperation(
                substrate, box_domain, vls.BooleanOperationEnum.RELATIVE_COMPLEMENT
            ).apply()

            mesh_cavity = vls.Mesh()
            tsm_cavity = module.ToSurfaceMesh(substrate, mesh_cavity)
            tsm_cavity.setSharpCorners(True)
            tsm_cavity.apply()

            nodes_cavity = mesh_cavity.getNodes()

            # Expected corners for cavity
            expected_corners_cavity = []
            for i in range(num_corners):
                corner = []
                for j in range(dim - 1):
                    val = 1.0 if ((i >> j) & 1) else -1.0
                    corner.append(val)
                # Last dimension
                val_last = 0.025 if ((i >> (dim - 1)) & 1) else -1.0
                corner.append(val_last)
                expected_corners_cavity.append(corner)

            for ex in expected_corners_cavity:
                found = False
                for n in nodes_cavity:
                    dist_sq = sum([(n[k] - ex[k]) ** 2 for k in range(dim)])
                    if dist_sq < tolerance:
                        found = True
                        break
                if not found:
                    self.fail(f"Cavity Corner {ex} not found in {dim}D mesh")

        print("  > Done test_corner_preservation")

def run_test_method(method_name):
    """Run a single test method in a separate process."""
    suite = unittest.TestSuite()
    suite.addTest(WrappingTest(method_name))
    # buffer=True captures stdout/stderr to prevent interleaving
    runner = unittest.TextTestRunner(verbosity=2, buffer=True)
    result = runner.run(suite)
    return result.wasSuccessful()

if __name__ == "__main__":
    vls.setNumThreads(8)
    unittest.main()
