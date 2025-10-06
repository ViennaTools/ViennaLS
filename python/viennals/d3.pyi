"""
3D bindings
"""

from __future__ import annotations
import collections.abc
import typing
import viennals._core

__all__: list[str] = [
    "Advect",
    "BooleanOperation",
    "Box",
    "BoxDistribution",
    "CalculateCurvatures",
    "CalculateNormalVectors",
    "CalculateVisibilities",
    "Check",
    "ConvexHull",
    "CustomSphereDistribution",
    "Cylinder",
    "DetectFeatures",
    "Domain",
    "Expand",
    "FinalizeStencilLocalLaxFriedrichs",
    "FromMesh",
    "FromSurfaceMesh",
    "FromVolumeMesh",
    "GeometricAdvect",
    "GeometricAdvectDistribution",
    "MakeGeometry",
    "MarkVoidPoints",
    "Plane",
    "PointCloud",
    "PrepareStencilLocalLaxFriedrichs",
    "Prune",
    "Reader",
    "Reduce",
    "RemoveStrayPoints",
    "Sphere",
    "SphereDistribution",
    "StencilLocalLaxFriedrichsScalar",
    "ToDiskMesh",
    "ToMesh",
    "ToSurfaceMesh",
    "ToVoxelMesh",
    "WriteVisualizationMesh",
    "Writer",
    "hrleGrid",
]

class Advect:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: viennals._core.VelocityField) -> None: ...
    def apply(self) -> None:
        """
        Perform advection.
        """

    def applyIntegration(self, arg0: typing.SupportsFloat) -> None:
        """
        Apply the integration scheme and calculate rates and maximum time step, but it do **not** move the surface.
        """

    def getAdvectedTime(self) -> float:
        """
        Get the time passed during advection.
        """

    def getCalculateNormalVectors(self) -> bool:
        """
        Get whether normal vectors are computed during advection.
        """

    def getCurrentTimeStep(self) -> float:
        """
        Get the current time step.
        """

    def getNumberOfTimeSteps(self) -> int:
        """
        Get how many advection steps were performed after the last apply() call.
        """

    def getTimeStepRatio(self) -> float:
        """
        Get the time step ratio used for advection.
        """

    def insertNextLevelSet(self, arg0: Domain) -> None:
        """
        Insert next level set to use for advection.
        """

    def prepareLS(self) -> None:
        """
        Prepare the level-set.
        """

    def setAdvectionTime(self, arg0: typing.SupportsFloat) -> None:
        """
        Set the time until when the level set should be advected.
        """

    def setCalculateNormalVectors(self, arg0: bool) -> None:
        """
        Set whether normal vectors are needed for the supplied velocity field.
        """

    def setDissipationAlpha(self, arg0: typing.SupportsFloat) -> None:
        """
        Set the dissipation value to use for Lax Friedrichs integration.
        """

    def setIgnoreVoids(self, arg0: bool) -> None:
        """
        Set whether voids in the geometry should be ignored during advection or not.
        """

    def setIntegrationScheme(self, arg0: viennals._core.IntegrationSchemeEnum) -> None:
        """
        Set the integration scheme to use during advection.
        """

    def setTimeStepRatio(self, arg0: typing.SupportsFloat) -> None:
        """
        Set the maximum time step size relative to grid size. Advection is only stable for <0.5.
        """

    def setVelocityField(self, arg0: viennals._core.VelocityField) -> None:
        """
        Set the velocity to use for advection.
        """

class BooleanOperation:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: Domain) -> None: ...
    @typing.overload
    def __init__(
        self, arg0: Domain, arg1: viennals._core.BooleanOperationEnum
    ) -> None: ...
    @typing.overload
    def __init__(
        self, arg0: Domain, arg1: Domain, arg2: viennals._core.BooleanOperationEnum
    ) -> None: ...
    def apply(self) -> None:
        """
        Perform the boolean operation.
        """

    def setBooleanOperation(self, arg0: viennals._core.BooleanOperationEnum) -> None:
        """
        Set which type of boolean operation should be performed.
        """

    def setLevelset(self, arg0: Domain) -> None:
        """
        Set levelset on which the boolean operation should be performed.
        """

    def setSecondLevelSet(self, arg0: Domain) -> None:
        """
        Set second levelset for boolean operation.
        """

class Box:
    def __init__(
        self,
        minPoint: collections.abc.Sequence[typing.SupportsFloat],
        maxPoint: collections.abc.Sequence[typing.SupportsFloat],
    ) -> None: ...

class BoxDistribution(GeometricAdvectDistribution):
    def __init__(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
    ) -> None: ...
    def getBounds(self) -> typing.Annotated[list[float], "FixedSize(6)"]:
        """
        Get the cartesian bounds of the distribution.
        """

    def getSignedDistance(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg1: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg2: typing.SupportsInt,
    ) -> float:
        """
        Get the signed distance of the passed point to the surface of the distribution.
        """

    def isInside(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg1: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg2: typing.SupportsFloat,
    ) -> bool:
        """
        Check whether passed point is inside the distribution.
        """

class CalculateCurvatures:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: viennals._core.CurvatureEnum) -> None: ...
    def apply(self) -> None:
        """
        Perform normal vector calculation.
        """

    def setCurvatureType(self, arg0: viennals._core.CurvatureEnum) -> None:
        """
        Set which method to use for calculation: Defaults to mean curvature.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset for which to calculate the curvatures.
        """

    def setMaxValue(self, arg0: typing.SupportsFloat) -> None:
        """
        Curvatures will be calculated for all LS values < maxValue.
        """

class CalculateNormalVectors:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    def apply(self) -> None:
        """
        Perform normal vector calculation.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset for which to calculate normal vectors.
        """

class CalculateVisibilities:
    def __init__(
        self,
        arg0: Domain,
        arg1: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg2: str,
    ) -> None: ...
    def apply(self) -> None: ...

class Check:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    def apply(self) -> None:
        """
        Perform check.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset for which to calculate normal vectors.
        """

class ConvexHull:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: viennals._core.Mesh, arg1: PointCloud) -> None: ...
    def apply(self) -> None:
        """
        Generate Hull.
        """

    def setMesh(self, arg0: viennals._core.Mesh) -> None:
        """
        Set mesh object where the generated mesh should be stored.
        """

    def setPointCloud(self, arg0: PointCloud) -> None:
        """
        Set point cloud used to generate mesh.
        """

class CustomSphereDistribution(GeometricAdvectDistribution):
    def __init__(
        self, arg0: collections.abc.Sequence[typing.SupportsFloat]
    ) -> None: ...
    def getBounds(self) -> typing.Annotated[list[float], "FixedSize(6)"]:
        """
        Get the cartesian bounds of the distribution.
        """

    def getSignedDistance(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg1: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg2: typing.SupportsInt,
    ) -> float:
        """
        Get the signed distance of the passed point to the surface of the distribution.
        """

    def isInside(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg1: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg2: typing.SupportsFloat,
    ) -> bool:
        """
        Check whether passed point is inside the distribution.
        """

class Cylinder:
    def __init__(
        self,
        origin: collections.abc.Sequence[typing.SupportsFloat],
        axisDirection: collections.abc.Sequence[typing.SupportsFloat],
        height: typing.SupportsFloat,
        radius: typing.SupportsFloat,
        topRadius: typing.SupportsFloat = 0.0,
    ) -> None: ...

class DetectFeatures:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: typing.SupportsFloat) -> None: ...
    @typing.overload
    def __init__(
        self,
        arg0: Domain,
        arg1: typing.SupportsFloat,
        arg2: viennals._core.FeatureDetectionEnum,
    ) -> None: ...
    def apply(self) -> None:
        """
        Detect features.
        """

    def setDetectionMethod(self, arg0: viennals._core.FeatureDetectionEnum) -> None:
        """
        Set which method to use to detect features. Defaults to Curvature.
        """

    def setDetectionThreshold(self, arg0: typing.SupportsFloat) -> None:
        """
        Set the curvature value above which a point is considered a feature.
        """

class Domain:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, gridDelta: typing.SupportsFloat = 1.0) -> None: ...
    @typing.overload
    def __init__(
        self,
        bounds: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(6)"
        ],
        boundaryConditions: typing.Annotated[
            collections.abc.Sequence[viennals._core.BoundaryConditionEnum],
            "FixedSize(3)",
        ],
        gridDelta: typing.SupportsFloat = 1.0,
    ) -> None: ...
    @typing.overload
    def __init__(
        self,
        bounds: collections.abc.Sequence[typing.SupportsFloat],
        boundaryConditions: collections.abc.Sequence[typing.SupportsInt],
        gridDelta: typing.SupportsFloat = 1.0,
    ) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    @typing.overload
    def __init__(self, arg0: hrleGrid) -> None: ...
    def clearMetaData(self) -> None:
        """
        Clear all metadata stored in the level set.
        """

    def deepCopy(self, arg0: Domain) -> None:
        """
        Copy lsDomain in this lsDomain.
        """

    def getLevelSetWidth(self) -> int:
        """
        Get the number of layers of level set points around the explicit surface.
        """

    def getNumberOfPoints(self) -> int:
        """
        Get the number of defined level set values.
        """

    def getNumberOfSegments(self) -> int:
        """
        Get the number of segments, the level set structure is divided into.
        """

    def print(self, stream: ... = ...) -> None: ...
    def setLevelSetWidth(self, arg0: typing.SupportsInt) -> None:
        """
        Set the number of layers of level set points which should be stored around the explicit surface.
        """

class Expand:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: typing.SupportsInt) -> None: ...
    def apply(self) -> None:
        """
        Perform expansion.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset to expand.
        """

    def setWidth(self, arg0: typing.SupportsInt) -> None:
        """
        Set the width to expand to.
        """

class FromMesh:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: viennals._core.Mesh) -> None: ...
    def apply(self) -> None: ...
    def setMesh(self, arg0: viennals._core.Mesh) -> None:
        """
        Set the mesh to read from.
        """

    def setSortPointList(self, arg0: bool) -> None: ...

class FromSurfaceMesh:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: viennals._core.Mesh) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: viennals._core.Mesh, arg2: bool) -> None: ...
    def apply(self) -> None:
        """
        Construct a levelset from a surface mesh.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset to read into.
        """

    def setMesh(self, arg0: viennals._core.Mesh) -> None:
        """
        Set the mesh to read from.
        """

    @typing.overload
    def setRemoveBoundaryTriangles(self, arg0: bool) -> None:
        """
        Set whether to include mesh elements outside of the simulation domain.
        """

    @typing.overload
    def setRemoveBoundaryTriangles(
        self, arg0: typing.Annotated[collections.abc.Sequence[bool], "FixedSize(3)"]
    ) -> None:
        """
        Set whether to include mesh elements outside of the simulation domain.
        """

class FromVolumeMesh:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: hrleGrid, arg1: viennals._core.Mesh) -> None: ...
    @typing.overload
    def __init__(
        self, arg0: hrleGrid, arg1: viennals._core.Mesh, arg2: bool
    ) -> None: ...
    def apply(self) -> None:
        """
        Construct a levelset from a volume mesh.
        """

    def setGrid(self, arg0: hrleGrid) -> None:
        """
        Set the grid used to read in the level sets.
        """

    def setMesh(self, arg0: viennals._core.Mesh) -> None:
        """
        Set the mesh to read from.
        """

    def setRemoveBoundaryTriangles(self, arg0: bool) -> None:
        """
        Set whether to include mesh elements outside of the simulation domain.
        """

class GeometricAdvect:
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @staticmethod
    def setAdvectionDistribution(*args, **kwargs) -> None:
        """
        Set advection distribution to use as kernel for the fast advection.
        """

    @typing.overload
    def __init__(self) -> None: ...
    def apply(self) -> None:
        """
        Perform advection.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset to advect.
        """

class GeometricAdvectDistribution:
    def __init__(self) -> None: ...
    def finalize(self) -> None:
        """
        Finalize the distribution after use with the level set.
        """

    def getBounds(self) -> typing.Annotated[list[float], "FixedSize(6)"]:
        """
        Get the cartesian bounds of the distribution.
        """

    def getSignedDistance(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg1: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg2: typing.SupportsInt,
    ) -> float:
        """
        Get the signed distance of the passed point to the surface of the distribution.
        """

    def isInside(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg1: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg2: typing.SupportsFloat,
    ) -> bool:
        """
        Check whether passed point is inside the distribution.
        """

    def prepare(self, arg0: Domain) -> None:
        """
        Prepare the distribution for use with the passed level set.
        """

class MakeGeometry:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: Sphere) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: Plane) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: Box) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: Cylinder) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: PointCloud) -> None: ...
    def apply(self) -> None:
        """
        Generate the geometry.
        """

    @typing.overload
    def setGeometry(self, arg0: Sphere) -> None: ...
    @typing.overload
    def setGeometry(self, arg0: Plane) -> None: ...
    @typing.overload
    def setGeometry(self, arg0: Box) -> None: ...
    @typing.overload
    def setGeometry(self, arg0: Cylinder) -> None: ...
    @typing.overload
    def setGeometry(self, arg0: PointCloud) -> None: ...
    @typing.overload
    def setIgnoreBoundaryConditions(self, arg0: bool) -> None: ...
    @typing.overload
    def setIgnoreBoundaryConditions(
        self, arg0: typing.Annotated[collections.abc.Sequence[bool], "FixedSize(3)"]
    ) -> None: ...
    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set the levelset in which to create the geometry.
        """

class MarkVoidPoints:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: bool) -> None: ...
    def apply(self) -> None:
        """
        Mark void points.
        """

    def setDetectLargestSurface(self, arg0: bool) -> None:
        """
        Set that the top surface should be the one with the most connected LS points.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set the levelset to mark void points in.
        """

    def setReverseVoidDetection(self, arg0: bool) -> None:
        """
        Reverse the logic of detecting the top surface.
        """

    def setSaveComponentsId(self, arg0: bool) -> None:
        """
        Save the connectivity information of all LS points in the pointData of the level set.
        """

    def setVoidTopSurface(self, arg0: viennals._core.VoidTopSurfaceEnum) -> None:
        """
        Set the logic by which to choose the surface which is non-void. All other connected surfaces will then be marked as void points.
        """

class Plane:
    def __init__(
        self,
        origin: collections.abc.Sequence[typing.SupportsFloat],
        normal: collections.abc.Sequence[typing.SupportsFloat],
    ) -> None: ...

class PointCloud:
    def __init__(
        self,
        arg0: collections.abc.Sequence[
            typing.Annotated[
                collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
            ]
        ],
    ) -> None: ...
    def insertNextPoint(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
    ) -> None: ...

class Prune:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    def apply(self) -> None:
        """
        Perform pruning operation.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset to prune.
        """

class Reader:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: str) -> None: ...
    def apply(self) -> None:
        """
        Write to file.
        """

    def setFileName(self, arg0: str) -> None:
        """
        Set the filename for the output file.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset to write to file.
        """

class Reduce:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: typing.SupportsInt) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: typing.SupportsInt, arg2: bool) -> None: ...
    def apply(self) -> None:
        """
        Perform reduction.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset to reduce.
        """

    def setNoNewSegment(self, arg0: bool) -> None:
        """
        Set whether the levelset should be segmented anew (balanced across cores) after reduction.
        """

    def setWidth(self, arg0: typing.SupportsInt) -> None:
        """
        Set the width to reduce to.
        """

class RemoveStrayPoints:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    def apply(self) -> None:
        """
        Remove stray points.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset for stray point removal.
        """

    def setVoidTopSurface(self, arg0: viennals._core.VoidTopSurfaceEnum) -> None:
        """
        Set the logic by which to choose the surface which should be kept. All other LS values will be marked as stray points and removed.
        """

class Sphere:
    def __init__(
        self,
        origin: collections.abc.Sequence[typing.SupportsFloat],
        radius: typing.SupportsFloat,
    ) -> None: ...

class SphereDistribution(GeometricAdvectDistribution):
    def __init__(self, arg0: typing.SupportsFloat) -> None: ...
    def getBounds(self) -> typing.Annotated[list[float], "FixedSize(6)"]:
        """
        Get the cartesian bounds of the distribution.
        """

    def getSignedDistance(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg1: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg2: typing.SupportsInt,
    ) -> float:
        """
        Get the signed distance of the passed point to the surface of the distribution.
        """

    def isInside(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg1: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg2: typing.SupportsFloat,
    ) -> bool:
        """
        Check whether passed point is inside the distribution.
        """

class StencilLocalLaxFriedrichsScalar:
    @staticmethod
    def setMaxDissipation(maxDissipation: typing.SupportsFloat) -> None: ...

class ToDiskMesh:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: viennals._core.Mesh) -> None: ...
    def apply(self) -> None:
        """
        Convert the levelset to a surface mesh.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset to mesh.
        """

    def setMesh(self, arg0: viennals._core.Mesh) -> None:
        """
        Set the mesh to generate.
        """

class ToMesh:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: viennals._core.Mesh) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: viennals._core.Mesh, arg2: bool) -> None: ...
    @typing.overload
    def __init__(
        self, arg0: Domain, arg1: viennals._core.Mesh, arg2: bool, arg3: bool
    ) -> None: ...
    def apply(self) -> None:
        """
        Convert the levelset to a surface mesh.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset to mesh.
        """

    def setMesh(self, arg0: viennals._core.Mesh) -> None:
        """
        Set the mesh to generate.
        """

    def setOnlyActive(self, arg0: bool) -> None:
        """
        Set whether only level set points <0.5 should be output.
        """

    def setOnlyDefined(self, arg0: bool) -> None:
        """
        Set whether only defined points should be output to the mesh.
        """

class ToSurfaceMesh:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: viennals._core.Mesh) -> None: ...
    def apply(self) -> None:
        """
        Convert the levelset to a surface mesh.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset to mesh.
        """

    def setMesh(self, arg0: viennals._core.Mesh) -> None:
        """
        Set the mesh to generate.
        """

class ToVoxelMesh:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: viennals._core.Mesh) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: viennals._core.Mesh) -> None: ...
    @typing.overload
    def __init__(
        self, arg0: collections.abc.Sequence[Domain], arg1: viennals._core.Mesh
    ) -> None: ...
    def apply(self) -> None:
        """
        Convert the levelset to a surface mesh.
        """

    def insertNextLevelSet(self, arg0: Domain) -> None:
        """
        Insert next level set to output in the mesh.
        """

    def setMesh(self, arg0: viennals._core.Mesh) -> None:
        """
        Set the mesh to generate.
        """

class WriteVisualizationMesh:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    @typing.overload
    def addMetaData(self, arg0: str, arg1: typing.SupportsFloat) -> None:
        """
        Add a single metadata entry to the file.
        """

    @typing.overload
    def addMetaData(
        self, arg0: str, arg1: collections.abc.Sequence[typing.SupportsFloat]
    ) -> None:
        """
        Add a single metadata entry to the file.
        """

    @typing.overload
    def addMetaData(
        self,
        arg0: collections.abc.Mapping[
            str, collections.abc.Sequence[typing.SupportsFloat]
        ],
    ) -> None:
        """
        Add metadata to the file.
        """

    def apply(self) -> None:
        """
        Make and write mesh.
        """

    def insertNextLevelSet(self, arg0: Domain) -> None:
        """
        Insert next level set to convert. Bigger level sets wrapping smaller ones, should be inserted last.
        """

    def setExtractHullMesh(self, arg0: bool) -> None:
        """
        Whether to extract a hull mesh. Defaults to false.
        """

    def setExtractVolumeMesh(self, arg0: bool) -> None:
        """
        Whether to extract a tetra volume mesh. Defaults to true.
        """

    def setFileName(self, arg0: str) -> None:
        """
        Set Name of File to write.
        """

    def setMetaData(
        self,
        arg0: collections.abc.Mapping[
            str, collections.abc.Sequence[typing.SupportsFloat]
        ],
    ) -> None:
        """
        Set the metadata to be written to the file.
        """

class Writer:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: str) -> None: ...
    def apply(self) -> None:
        """
        Write to file.
        """

    def setFileName(self, arg0: str) -> None:
        """
        Set the filename for the output file.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset to write to file.
        """

class hrleGrid:
    pass

def FinalizeStencilLocalLaxFriedrichs(
    levelSets: collections.abc.Sequence[Domain],
) -> None: ...
def PrepareStencilLocalLaxFriedrichs(
    levelSets: collections.abc.Sequence[Domain], isDepo: collections.abc.Sequence[bool]
) -> None: ...
