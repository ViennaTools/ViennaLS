# ViennaLS Python API

This document lists all functions, classes, and enums available in the **ViennaLS** Python package.  
Functions are grouped into **dimension-bound** and **dimension-independent** components.

---

## **Dimension-Bound Functions**

These functions are tied to either the **2D** or **3D** version of the library.  
They must be imported from either `viennals.d2` or `viennals.d3`.

### **Core Classes and Operations**
- `Domain`
- `Advect`
- `BooleanOperation`
- `CalculateCurvatures`
- `CalculateNormalVectors`
- `CalculateVisibilities`
- `Check`
- `PointCloud`
- `ConvexHull`
- `DetectFeatures`
- `GeometricAdvect`
- `GeometricAdvectDistributions`
- `SphereDistribution`
- `BoxDistribution`
- `Expand`
- `FromSurfaceMesh`
- `FromVolumeMesh`
- `FromMesh`
- `Sphere`
- `Plane`
- `Box`
- `Cylinder`
- `MakeGeometry`
- `MarkVoidPoints`
- `Prune`
- `Reader`
- `Reduce`
- `RemoveStrayPoints`
- `ToDiskMesh`
- `ToMesh`
- `ToMultiSurfaceMesh`
- `ToSurfaceMesh`
- `ToVoxelMesh`
- `Writer`
- `WriteVisualizationMesh`

### **2D-Only Functions**
- `CompareArea`
- `CompareNarrowBand`
- `CompareSparseField`

---

## **Dimension-Independent Components**

These are imported directly from the **root** `viennals` module.

### **Enums**
- `LogLevel`
- `DiscretizationSchemeEnum`
- `BooleanOperationEnum`
- `CurvatureEnum`
- `FeatureDetectionEnum`
- `BoundaryConditionEnum`
- `FileFormatEnum`
- `VoidTopSurfaceEnum`
- `TransformEnum`

### **Mesh and I/O Classes**
- `Mesh`
- `PointData`
- `TransformMesh`
- `VTKReader`
- `VTKWriter`

### **Utilities**
- `Logger`
- `MaterialMap`
- `VelocityField`

---

## **Transformations Between Dimensions**

These functions enable conversion between 2D and 3D representations and are also imported from the root `viennals` module:
- `Slice`
- `Extrude`

---

## **Example Imports**

```python
# 2D functions and domain-specific operations
import viennals.d2 as vls

domain = vls.Domain()
vls.Advect(domain, ...)

# Common functions and enums
import viennals

mesh = viennals.Mesh()
enum = viennals.BooleanOperationEnum.UNION
