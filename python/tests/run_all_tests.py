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

    def test_3d_algorithms(self):
        print("\n[TEST] test_3d_algorithms")
        # Use explicit d3 module to be safe
        gridDelta = 0.2
        bounds = [-10.0, 10.0, -10.0, 10.0, -10.0, 10.0]
        dom = d3.Domain(bounds, [vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY] * 3, gridDelta)
        
        # Geometries
        origin = [0.0, 0.0, 0.0]
        radius = 5.0
        print("  > MakeGeometry")
        sphere_geom = d3.Sphere(origin, radius)
        d3.MakeGeometry(dom, sphere_geom).apply()
        self.assertGreater(dom.getNumberOfPoints(), 0)
        
        print("  > Advect")
        # Advect
        vel = TestVelocityField()
        adv = d3.Advect()
        adv.insertNextLevelSet(dom)
        adv.setVelocityField(vel)
        adv.setAdvectionTime(0.1)
        adv.setSpatialScheme(vls.SpatialSchemeEnum.ENGQUIST_OSHER_1ST_ORDER)
        adv.setTemporalScheme(vls.TemporalSchemeEnum.FORWARD_EULER)
        adv.apply()
        
        print("  > GeometricAdvect")
        # Geometric Advect & Distributions
        dist = d3.SphereDistribution(radius)
        gadv = d3.GeometricAdvect(dom, dist)
        gadv.apply()
        
        # Check other distributions instantiation
        print("  > Checking Distribution instantiation")
        _ = d3.BoxDistribution([1.0, 1.0, 1.0])
        _ = d3.CustomSphereDistribution([1.0, 2.0])
        
        print("  > Boolean")
        # Boolean
        dom2 = d3.Domain(dom)
        bool_op = d3.BooleanOperation(dom, dom2, vls.BooleanOperationEnum.UNION)
        bool_op.apply()
        
        print("  > Calc Curvatures")
        # Calc Curvatures
        curv = d3.CalculateCurvatures(dom)
        curv.apply()
        
        print("  > Calc Normals")
        # Calc Normals
        norms = d3.CalculateNormalVectors(dom)
        norms.apply()
        
        print("  > Calc Visibilities")
        # Calc Visibilities
        d3.CalculateVisibilities(dom, [0.0, 0.0, 1.0], "Visibility").apply()
        
        print("  > Check")
        # Check
        d3.Check(dom).apply()
        
        print("  > Prune")
        # Prune
        d3.Prune(dom).apply()
        
        print("  > Expand")
        # Expand
        d3.Expand(dom, 1).apply()
        
        print("  > Reduce")
        # Reduce
        d3.Reduce(dom, 1).apply()
        
        print("  > MarkVoidPoints")
        # Mark Void Points
        d3.MarkVoidPoints(dom).apply()
        
        # Expand to ensure sufficient width for feature detection (curvature)
        d3.Prune(dom).apply()
        d3.Expand(dom, 5).apply()

        print("  > DetectFeatures")
        # Detect Features
        d3.DetectFeatures(dom).apply()
        
        print("  > Mesh Conversions")
        # Mesh Conversions
        # Expand back to a safe width for meshing
        d3.Expand(dom, 3).apply()
        print("    >> ToMesh")

        mesh = vls.Mesh()
        print("    >> ToMesh")
        d3.ToMesh(dom, mesh).apply()
        print("    >> ToSurfaceMesh")
        d3.ToSurfaceMesh(dom, mesh).apply()
        print("    >> ToDiskMesh")
        d3.ToDiskMesh(dom, mesh).apply()
        
        print("  > IO")
        # IO
        vls.VTKWriter(mesh, "test_3d.vtp").apply()
        if os.path.exists("test_3d.vtp"):
            os.remove("test_3d.vtp")
        
        print("  > RemoveStrayPoints")
        # Remove Stray Points
        d3.RemoveStrayPoints(dom).apply()

        print("  > Done test_3d_algorithms")

    def test_geometries(self):
        print("\n[TEST] test_geometries")
        # Initialize domain with bounds to ensure MakeGeometry(Plane) works efficiently
        bounds = [-10.0, 10.0, -10.0, 10.0, -10.0, 10.0]
        boundaryCons = [vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
                        vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY,
                        vls.BoundaryConditionEnum.INFINITE_BOUNDARY]
        dom = d3.Domain(bounds, boundaryCons, 0.2)
        
        # Plane
        print("  > Plane")
        plane = d3.Plane([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
        d3.MakeGeometry(dom, plane).apply()
        self.assertGreater(dom.getNumberOfPoints(), 0)
        
        # Box
        print("  > Box")
        box = d3.Box([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])
        d3.MakeGeometry(dom, box).apply()
        
        # Cylinder
        print("  > Cylinder")
        cyl = d3.Cylinder([0.0,0.0,0.0], [0.0,0.0,1.0], 2.0, 1.0, 1.0)
        d3.MakeGeometry(dom, cyl).apply()
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
        dom = d3.Domain(0.2)
        d3.MakeGeometry(dom, d3.Sphere([0.0,0.0,0.0], 5.0)).apply()
        
        mesh = vls.Mesh()
        d3.ToSurfaceMesh(dom, mesh).apply()
        
        dom2 = d3.Domain(0.2)
        d3.FromSurfaceMesh(dom2, mesh).apply()
        self.assertGreater(dom2.getNumberOfPoints(), 0)
        print("  > Done test_mesh_conversion_roundtrip")

    def test_native_io(self):
        print("\n[TEST] test_native_io")
        dom = d3.Domain(0.2)
        d3.MakeGeometry(dom, d3.Sphere([0.0, 0.0, 0.0], 5.0)).apply()
        
        filename = "test_native.lvst"
        print("  > Writer")
        d3.Writer(dom, filename).apply()
        
        dom_in = d3.Domain(0.2)
        print("  > Reader")
        d3.Reader(dom_in, filename).apply()
        self.assertGreater(dom_in.getNumberOfPoints(), 0)
        
        if os.path.exists(filename):
            os.remove(filename)
        print("  > Done test_native_io")

    def test_advanced_meshing(self):
        print("\n[TEST] test_advanced_meshing")
        dom1 = d3.Domain(0.2)
        d3.MakeGeometry(dom1, d3.Sphere([0.0, 0.0, 0.0], 5.0)).apply()
        
        dom2 = d3.Domain(0.2)
        d3.MakeGeometry(dom2, d3.Sphere([2.0, 0.0, 0.0], 5.0)).apply()
        
        # ToMultiSurfaceMesh
        print("  > ToMultiSurfaceMesh")
        mesh_multi = vls.Mesh()
        d3.ToMultiSurfaceMesh([dom1, dom2], mesh_multi).apply()
        self.assertGreater(len(mesh_multi.getNodes()), 0)
        
        # ToVoxelMesh
        print("  > ToVoxelMesh")
        mesh_voxel = vls.Mesh()
        d3.ToVoxelMesh([dom1, dom2], mesh_voxel).apply()
        self.assertGreater(len(mesh_voxel.getNodes()), 0)
        print("  > Done test_advanced_meshing")

    def test_stencil_functions(self):
        print("\n[TEST] test_stencil_functions")
        dom = d3.Domain(0.2)
        d3.MakeGeometry(dom, d3.Sphere([0.0, 0.0, 0.0], 5.0)).apply()
        d3.Expand(dom, 5).apply()

        # Create a second domain to form a valid multi-layer system
        dom2 = d3.Domain(0.2)
        d3.MakeGeometry(dom2, d3.Sphere([0.0, 0.0, 0.0], 10.0)).apply()
        d3.Expand(dom2, 5).apply()

        # PrepareStencilLocalLaxFriedrichs
        print("  > PrepareStencilLocalLaxFriedrichs")
        d3.PrepareStencilLocalLaxFriedrichs([dom, dom2], [False, True])
        
        # Finalize
        print("  > FinalizeStencilLocalLaxFriedrichs")
        d3.FinalizeStencilLocalLaxFriedrichs([dom, dom2])
        print("  > Done test_stencil_functions")

    def test_geometric_advect_box(self):
        print("\n[TEST] test_geometric_advect_box")
        dom = d3.Domain(0.2)
        d3.MakeGeometry(dom, d3.Sphere([0.0, 0.0, 0.0], 5.0)).apply()
        
        print("  > BoxDistribution")
        dist = d3.BoxDistribution([1.0, 1.0, 1.0])
        
        print("  > GeometricAdvect")
        d3.GeometricAdvect(dom, dist).apply()
        self.assertGreater(dom.getNumberOfPoints(), 0)
        print("  > Done test_geometric_advect_box")

    def test_from_mesh_volume(self):
        print("\n[TEST] test_from_mesh_volume")
        # Create a simple tetrahedral mesh
        mesh = vls.Mesh()
        mesh.insertNextNode([0.0, 0.0, 0.0])
        mesh.insertNextNode([5.0, 0.0, 0.0])
        mesh.insertNextNode([0.0, 5.0, 0.0])
        mesh.insertNextNode([0.0, 0.0, 5.0])
        mesh.insertNextTetra([0, 1, 2, 3])
        
        # FromMesh requires "LSValues" in cellData corresponding to nodes
        ls_values = [-1.0, 1.0, 1.0, 1.0]
        mesh.getPointData().insertNextScalarData(ls_values, "LSValues")

        dom = d3.Domain(0.2)
        # FromMesh should handle volume meshes if implemented in the wrapper logic
        print("  > FromMesh (Volume)")
        d3.FromMesh(dom, mesh).apply()
        self.assertGreater(dom.getNumberOfPoints(), 0)
        print("  > Done test_from_mesh_volume")

    def test_boolean_variants(self):
        print("\n[TEST] test_boolean_variants")
        dom1 = d3.Domain(0.2)
        d3.MakeGeometry(dom1, d3.Sphere([0.0, 0.0, 0.0], 5.0)).apply()
        
        dom2 = d3.Domain(0.2)
        d3.MakeGeometry(dom2, d3.Sphere([2.0, 0.0, 0.0], 5.0)).apply()
        
        # Intersect
        print("  > Intersect")
        dom_int = d3.Domain(dom1)
        d3.BooleanOperation(dom_int, dom2, vls.BooleanOperationEnum.INTERSECT).apply()
        self.assertGreater(dom_int.getNumberOfPoints(), 0)
        
        # Relative Complement (Difference)
        print("  > Relative Complement")
        dom_diff = d3.Domain(dom1)
        d3.BooleanOperation(dom_diff, dom2, vls.BooleanOperationEnum.RELATIVE_COMPLEMENT).apply()
        self.assertGreater(dom_diff.getNumberOfPoints(), 0)
        
        # Invert
        print("  > Invert")
        dom_inv = d3.Domain(dom1)
        d3.BooleanOperation(dom_inv, vls.BooleanOperationEnum.INVERT).apply()
        # Inverted domain is infinite, but points should exist
        self.assertGreater(dom_inv.getNumberOfPoints(), 0)
        print("  > Done test_boolean_variants")

    def test_domain_manipulation(self):
        print("\n[TEST] test_domain_manipulation")
        dom = d3.Domain(0.2)
        d3.MakeGeometry(dom, d3.Sphere([0.0, 0.0, 0.0], 5.0)).apply()
        
        # Deep Copy
        print("  > deepCopy")
        dom_copy = d3.Domain(0.2)
        dom_copy.deepCopy(dom)
        self.assertEqual(dom.getNumberOfPoints(), dom_copy.getNumberOfPoints())
        
        # Level Set Width
        print("  > setLevelSetWidth")
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

    def test_2d_algorithms(self):
        print("\n[TEST] test_2d_algorithms")
        # Use explicit d2 module
        gridDelta = 0.2
        gridDelta = 0.5
        print("  > Creating 2D domains")
        dom = d2.Domain(gridDelta)
        
        # Make Geometry
        d2.MakeGeometry(dom, d2.Sphere([0.0, 0.0], 5.0)).apply()
        
        dom2 = d2.Domain(dom)
        # Shift dom2 slightly
        d2.MakeGeometry(dom2, d2.Sphere([0.1, 0.0], 5.0)).apply()
        
        print("  > Comparisons")
        # Comparisons
        d2.CompareArea(dom, dom2).apply()
        d2.CompareChamfer(dom, dom2).apply()
        
        comp_crit = d2.CompareCriticalDimensions(dom, dom2)
        comp_crit.addXRange(-1.0, 1.0, True)
        comp_crit.apply()
        
        d2.CompareNarrowBand(dom, dom2).apply()
        
        print("  > Sparse Field")
        # Sparse Field (needs expanded/reduced)
        dom_exp = d2.Domain(dom)
        d2.Expand(dom_exp, 55).apply()
        d2.CompareSparseField(dom_exp, dom).apply()
        print("  > Done test_2d_algorithms")

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
        # PointCloud
        print("  > Creating PointCloud")
        points = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        pc = d3.PointCloud(points)
        
        # ConvexHull
        print("  > ConvexHull")
        mesh = vls.Mesh()
        hull = d3.ConvexHull(mesh, pc)
        hull.apply()
        self.assertGreater(len(mesh.getNodes()), 0)
        print("  > Done test_point_cloud_convex_hull")

    def test_advection_callback(self):
        print("\n[TEST] test_advection_callback")
        print("  > Setting up advection")
        dom = d3.Domain(0.2)
        d3.MakeGeometry(dom, d3.Sphere([0.0, 0.0, 0.0], 5.0)).apply()
        
        vel = TestVelocityField()
        adv = d3.Advect()
        adv.insertNextLevelSet(dom)
        adv.setVelocityField(vel)
        adv.setAdvectionTime(0.1)
        adv.setTemporalScheme(vls.TemporalSchemeEnum.RUNGE_KUTTA_2ND_ORDER)
        
        # Use a mutable container to capture callback execution
        callback_data = {'called': False}
        def callback(d):
            callback_data['called'] = True
            return True
            
        print("  > Running advection with callback")
        adv.setVelocityUpdateCallback(callback)
        adv.apply()
        
        self.assertTrue(callback_data['called'])
        print("  > Done test_advection_callback")

    def test_make_geometry_boundary(self):
        print("\n[TEST] test_make_geometry_boundary")
        print("  > Creating domain with bounds")
        bounds = [-5.0, 5.0, -5.0, 5.0, -5.0, 5.0]
        boundaryCons = [vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY] * 3
        dom = d3.Domain(bounds, boundaryCons, 1.0)
        
        print("  > MakeGeometry with IgnoreBoundaryConditions")
        maker = d3.MakeGeometry(dom, d3.Sphere([0.0, 0.0, 0.0], 4.0))
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
        print("  > Creating domain")
        dom = d3.Domain(0.2)
        d3.MakeGeometry(dom, d3.Sphere([0.0, 0.0, 0.0], 5.0)).apply()
        
        vel = TestVelocityField()
        adv = d3.Advect()
        adv.insertNextLevelSet(dom)
        adv.setVelocityField(vel)
        
        # Test Single Step
        print("  > Testing Single Step")
        adv.setAdvectionTime(10.0) # Large time that would normally require multiple steps
        adv.setSingleStep(True)
        adv.apply()
        self.assertEqual(adv.getNumberOfTimeSteps(), 1)
        print("  > Done test_advect_configuration")

    def test_mark_void_points_advanced(self):
        print("\n[TEST] test_mark_void_points_advanced")
        print("  > Creating disjoint spheres")
        dom = d3.Domain(0.2)
        # Create two disjoint spheres
        d3.MakeGeometry(dom, d3.Sphere([0.0, 0.0, 0.0], 3.0)).apply()
        
        dom2 = d3.Domain(0.2)
        d3.MakeGeometry(dom2, d3.Sphere([10.0, 0.0, 0.0], 3.0)).apply()
        
        # Union them to create a domain with two separate components
        print("  > Union")
        d3.BooleanOperation(dom, dom2, vls.BooleanOperationEnum.UNION).apply()
        
        print("  > MarkVoidPoints")
        marker = d3.MarkVoidPoints(dom)
        marker.setSaveComponentIds(True)
        marker.apply()
        
        # Should find 3 connected components (2 spheres + 1 background)
        self.assertEqual(marker.getNumberOfComponents(), 3)
        print("  > Done test_mark_void_points_advanced")

    def test_adaptive_time_stepping(self):
        print("\n[TEST] test_adaptive_time_stepping")
        print("  > Setting up advection")
        dom = d3.Domain(0.2)
        d3.MakeGeometry(dom, d3.Sphere([0.0, 0.0, 0.0], 5.0)).apply()
        
        vel = TestVelocityField()
        adv = d3.Advect()
        adv.insertNextLevelSet(dom)
        adv.setVelocityField(vel)
        adv.setAdvectionTime(0.1)
        
        # Enable adaptive time stepping
        print("  > Enabling Adaptive Time Stepping")
        adv.setAdaptiveTimeStepping(True, 10)
        adv.apply()
        
        self.assertGreater(dom.getNumberOfPoints(), 0)
        print("  > Done test_adaptive_time_stepping")

    def test_geometric_distribution_functions(self):
        print("\n[TEST] test_geometric_distribution_functions")
        print("  > SphereDistribution")
        # Sphere Distribution
        dist = d3.SphereDistribution(5.0)
        # Point inside (0,0,0) relative to center (0,0,0)
        print("  > Checking isInside")
        self.assertTrue(dist.isInside([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1e-6))
        # Point outside (6,0,0) relative to center (0,0,0)
        self.assertFalse(dist.isInside([0.0, 0.0, 0.0], [6.0, 0.0, 0.0], 1e-6))
        
        # Signed Distance
        print("  > Checking getSignedDistance")
        # At 0,0,0 distance should be -5.0
        dist_val = dist.getSignedDistance([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0)
        self.assertAlmostEqual(dist_val, -5.0, delta=1e-5)
        print("  > Done test_geometric_distribution_functions")

    def test_advect_additional_options(self):
        print("\n[TEST] test_advect_additional_options")
        print("  > Setting up advection")
        dom = d3.Domain(0.2)
        d3.MakeGeometry(dom, d3.Sphere([0.0, 0.0, 0.0], 5.0)).apply()
        
        vel = TestVelocityField()
        adv = d3.Advect()
        adv.insertNextLevelSet(dom)
        adv.setVelocityField(vel)
        adv.setAdvectionTime(0.1)
        
        print("  > Setting IgnoreVoids")
        # Test Ignore Voids
        adv.setIgnoreVoids(True)
        
        print("  > Setting DissipationAlpha")
        # Test Dissipation Alpha
        adv.setDissipationAlpha(0.5)
        
        print("  > Setting CheckDissipation")
        # Test Check Dissipation
        adv.setCheckDissipation(False)
        
        adv.apply()
        
        self.assertGreater(dom.getNumberOfPoints(), 0)
        print("  > Done test_advect_additional_options")

    def test_additional_options(self):
        print("\n[TEST] test_additional_options")
        print("  > Creating domain with bounds")
        bounds = [-10.0, 10.0, -10.0, 10.0, -10.0, 10.0]
        boundaryCons = [vls.BoundaryConditionEnum.REFLECTIVE_BOUNDARY] * 3
        dom = d3.Domain(bounds, boundaryCons, 0.2)
        d3.MakeGeometry(dom, d3.Sphere([0.0, 0.0, 0.0], 5.0)).apply()
        
        # ToMesh options
        print("  > ToMesh options (OnlyDefined, OnlyActive)")
        mesh = vls.Mesh()
        tomesh = d3.ToMesh(dom, mesh)
        tomesh.setOnlyDefined(True)
        tomesh.setOnlyActive(True)
        tomesh.apply()
        self.assertGreater(len(mesh.getNodes()), 0)
        
        # DetectFeatures options
        print("  > DetectFeatures options (Threshold, Method)")
        df = d3.DetectFeatures(dom)
        df.setDetectionThreshold(10.0)
        df.setDetectionMethod(vls.FeatureDetectionEnum.NORMALS_ANGLE)
        df.apply()
        
        # MarkVoidPoints options
        print("  > MarkVoidPoints options (Reverse, LargestSurface)")
        mvp = d3.MarkVoidPoints(dom)
        mvp.setReverseVoidDetection(True)
        mvp.setDetectLargestSurface(True)
        mvp.apply()
        
        # RemoveStrayPoints options
        print("  > RemoveStrayPoints options (VoidTopSurface)")
        rsp = d3.RemoveStrayPoints(dom)
        rsp.setVoidTopSurface(vls.VoidTopSurfaceEnum.LARGEST)
        rsp.apply()
        
        print("  > Done test_additional_options")

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
