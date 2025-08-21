"""
ViennaLS is a header-only C++ level set library developed for high performance topography simulations. The main design goals are simplicity and efficiency, tailored towards scientific simulations. ViennaLS can also be used for visualization applications, although this is not the main design target.
"""

from __future__ import annotations
import collections.abc
import typing

__all__: list[str] = [
    "Advect",
    "BooleanOperation",
    "BooleanOperationEnum",
    "BoundaryConditionEnum",
    "Box",
    "BoxDistribution",
    "CalculateCurvatures",
    "CalculateNormalVectors",
    "CalculateVisibilities",
    "Check",
    "ConvexHull",
    "CurvatureEnum",
    "Cylinder",
    "DetectFeatures",
    "Domain",
    "Expand",
    "FeatureDetectionEnum",
    "FileFormatEnum",
    "FromMesh",
    "FromSurfaceMesh",
    "FromVolumeMesh",
    "GeometricAdvect",
    "GeometricAdvectDistribution",
    "IntegrationSchemeEnum",
    "LogLevel",
    "Logger",
    "MakeGeometry",
    "MarkVoidPoints",
    "MaterialMap",
    "Mesh",
    "Plane",
    "PointCloud",
    "PointData",
    "Prune",
    "Reader",
    "Reduce",
    "RemoveStrayPoints",
    "Slice",
    "Sphere",
    "SphereDistribution",
    "StencilLocalLaxFriedrichsScalar",
    "ToDiskMesh",
    "ToMesh",
    "ToSurfaceMesh",
    "ToVoxelMesh",
    "TransformEnum",
    "TransformMesh",
    "VTKReader",
    "VTKWriter",
    "VelocityField",
    "VoidTopSurfaceEnum",
    "WriteVisualizationMesh",
    "Writer",
    "hrleGrid",
    "setNumThreads",
    "version",
]

class Advect:
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @staticmethod
    def insertNextLevelSet(*args, **kwargs) -> None:
        """
        Insert next level set to use for advection.
        """

    @typing.overload
    def __init__(self) -> None: ...
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

    def setIntegrationScheme(self, arg0: ...) -> None:
        """
        Set the integration scheme to use during advection.
        """

    def setTimeStepRatio(self, arg0: typing.SupportsFloat) -> None:
        """
        Set the maximum time step size relative to grid size. Advection is only stable for <0.5.
        """

    def setVelocityField(self, arg0: ...) -> None:
        """
        Set the velocity to use for advection.
        """

class BooleanOperation:
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @staticmethod
    def setLevelset(*args, **kwargs) -> None:
        """
        Set levelset on which the boolean operation should be performed.
        """

    @staticmethod
    def setSecondLevelSet(*args, **kwargs) -> None:
        """
        Set second levelset for boolean operation.
        """

    @typing.overload
    def __init__(self) -> None: ...
    def apply(self) -> None:
        """
        Perform the boolean operation.
        """

    def setBooleanOperation(self, arg0: ...) -> None:
        """
        Set which type of boolean operation should be performed.
        """

class BooleanOperationEnum:
    """
    Members:

      INTERSECT

      UNION

      RELATIVE_COMPLEMENT

      INVERT
    """

    INTERSECT: typing.ClassVar[
        BooleanOperationEnum
    ]  # value = <BooleanOperationEnum.INTERSECT: 0>
    INVERT: typing.ClassVar[
        BooleanOperationEnum
    ]  # value = <BooleanOperationEnum.INVERT: 3>
    RELATIVE_COMPLEMENT: typing.ClassVar[
        BooleanOperationEnum
    ]  # value = <BooleanOperationEnum.RELATIVE_COMPLEMENT: 2>
    UNION: typing.ClassVar[
        BooleanOperationEnum
    ]  # value = <BooleanOperationEnum.UNION: 1>
    __members__: typing.ClassVar[
        dict[str, BooleanOperationEnum]
    ]  # value = {'INTERSECT': <BooleanOperationEnum.INTERSECT: 0>, 'UNION': <BooleanOperationEnum.UNION: 1>, 'RELATIVE_COMPLEMENT': <BooleanOperationEnum.RELATIVE_COMPLEMENT: 2>, 'INVERT': <BooleanOperationEnum.INVERT: 3>}
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: typing.SupportsInt) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class BoundaryConditionEnum:
    """
    Members:

      REFLECTIVE_BOUNDARY

      INFINITE_BOUNDARY

      PERIODIC_BOUNDARY

      POS_INFINITE_BOUNDARY

      NEG_INFINITE_BOUNDARY
    """

    INFINITE_BOUNDARY: typing.ClassVar[
        BoundaryConditionEnum
    ]  # value = <BoundaryConditionEnum.INFINITE_BOUNDARY: 1>
    NEG_INFINITE_BOUNDARY: typing.ClassVar[
        BoundaryConditionEnum
    ]  # value = <BoundaryConditionEnum.NEG_INFINITE_BOUNDARY: 4>
    PERIODIC_BOUNDARY: typing.ClassVar[
        BoundaryConditionEnum
    ]  # value = <BoundaryConditionEnum.PERIODIC_BOUNDARY: 2>
    POS_INFINITE_BOUNDARY: typing.ClassVar[
        BoundaryConditionEnum
    ]  # value = <BoundaryConditionEnum.POS_INFINITE_BOUNDARY: 3>
    REFLECTIVE_BOUNDARY: typing.ClassVar[
        BoundaryConditionEnum
    ]  # value = <BoundaryConditionEnum.REFLECTIVE_BOUNDARY: 0>
    __members__: typing.ClassVar[
        dict[str, BoundaryConditionEnum]
    ]  # value = {'REFLECTIVE_BOUNDARY': <BoundaryConditionEnum.REFLECTIVE_BOUNDARY: 0>, 'INFINITE_BOUNDARY': <BoundaryConditionEnum.INFINITE_BOUNDARY: 1>, 'PERIODIC_BOUNDARY': <BoundaryConditionEnum.PERIODIC_BOUNDARY: 2>, 'POS_INFINITE_BOUNDARY': <BoundaryConditionEnum.POS_INFINITE_BOUNDARY: 3>, 'NEG_INFINITE_BOUNDARY': <BoundaryConditionEnum.NEG_INFINITE_BOUNDARY: 4>}
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: typing.SupportsInt) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

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
        arg1: typing.SupportsFloat,
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
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @staticmethod
    def setLevelSet(*args, **kwargs) -> None:
        """
        Set levelset for which to calculate the curvatures.
        """

    @typing.overload
    def __init__(self) -> None: ...
    def apply(self) -> None:
        """
        Perform normal vector calculation.
        """

    def setCurvatureType(self, arg0: ...) -> None:
        """
        Set which method to use for calculation: Defaults to mean curvature.
        """

    def setMaxValue(self, arg0: typing.SupportsFloat) -> None:
        """
        Curvatures will be calculated for all LS values < maxValue.
        """

class CalculateNormalVectors:
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @staticmethod
    def setLevelSet(*args, **kwargs) -> None:
        """
        Set levelset for which to calculate normal vectors.
        """

    @typing.overload
    def __init__(self) -> None: ...
    def apply(self) -> None:
        """
        Perform normal vector calculation.
        """

class CalculateVisibilities:
    @staticmethod
    def __init__(*args, **kwargs) -> None: ...
    def apply(self) -> None: ...

class Check:
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @staticmethod
    def setLevelSet(*args, **kwargs) -> None:
        """
        Set levelset for which to calculate normal vectors.
        """

    @typing.overload
    def __init__(self) -> None: ...
    def apply(self) -> None:
        """
        Perform check.
        """

class ConvexHull:
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @staticmethod
    def setPointCloud(*args, **kwargs) -> None:
        """
        Set point cloud used to generate mesh.
        """

    @typing.overload
    def __init__(self) -> None: ...
    def apply(self) -> None:
        """
        Generate Hull.
        """

    def setMesh(self, arg0: ...) -> None:
        """
        Set mesh object where the generated mesh should be stored.
        """

class CurvatureEnum:
    """
    Members:

      MEAN_CURVATURE

      GAUSSIAN_CURVATURE

      MEAN_AND_GAUSSIAN_CURVATURE
    """

    GAUSSIAN_CURVATURE: typing.ClassVar[
        CurvatureEnum
    ]  # value = <CurvatureEnum.GAUSSIAN_CURVATURE: 1>
    MEAN_AND_GAUSSIAN_CURVATURE: typing.ClassVar[
        CurvatureEnum
    ]  # value = <CurvatureEnum.MEAN_AND_GAUSSIAN_CURVATURE: 2>
    MEAN_CURVATURE: typing.ClassVar[
        CurvatureEnum
    ]  # value = <CurvatureEnum.MEAN_CURVATURE: 0>
    __members__: typing.ClassVar[
        dict[str, CurvatureEnum]
    ]  # value = {'MEAN_CURVATURE': <CurvatureEnum.MEAN_CURVATURE: 0>, 'GAUSSIAN_CURVATURE': <CurvatureEnum.GAUSSIAN_CURVATURE: 1>, 'MEAN_AND_GAUSSIAN_CURVATURE': <CurvatureEnum.MEAN_AND_GAUSSIAN_CURVATURE: 2>}
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: typing.SupportsInt) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

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
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @typing.overload
    def __init__(self) -> None: ...
    def apply(self) -> None:
        """
        Detect features.
        """

    def setDetectionMethod(self, arg0: ...) -> None:
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
            collections.abc.Sequence[...], "FixedSize(3)"
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
    def __init__(self, arg0: ...) -> None: ...
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

class FeatureDetectionEnum:
    """
    Members:

      CURVATURE

      NORMALS_ANGLE
    """

    CURVATURE: typing.ClassVar[
        FeatureDetectionEnum
    ]  # value = <FeatureDetectionEnum.CURVATURE: 0>
    NORMALS_ANGLE: typing.ClassVar[
        FeatureDetectionEnum
    ]  # value = <FeatureDetectionEnum.NORMALS_ANGLE: 1>
    __members__: typing.ClassVar[
        dict[str, FeatureDetectionEnum]
    ]  # value = {'CURVATURE': <FeatureDetectionEnum.CURVATURE: 0>, 'NORMALS_ANGLE': <FeatureDetectionEnum.NORMALS_ANGLE: 1>}
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: typing.SupportsInt) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class FileFormatEnum:
    """
    Members:

      VTK_LEGACY

      VTP

      VTU
    """

    VTK_LEGACY: typing.ClassVar[
        FileFormatEnum
    ]  # value = <FileFormatEnum.VTK_LEGACY: 0>
    VTP: typing.ClassVar[FileFormatEnum]  # value = <FileFormatEnum.VTP: 1>
    VTU: typing.ClassVar[FileFormatEnum]  # value = <FileFormatEnum.VTU: 2>
    __members__: typing.ClassVar[
        dict[str, FileFormatEnum]
    ]  # value = {'VTK_LEGACY': <FileFormatEnum.VTK_LEGACY: 0>, 'VTP': <FileFormatEnum.VTP: 1>, 'VTU': <FileFormatEnum.VTU: 2>}
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: typing.SupportsInt) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class FromMesh:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: ...) -> None: ...
    def apply(self) -> None: ...
    def setMesh(self, arg0: ...) -> None:
        """
        Set the mesh to read from.
        """

    def setSortPointList(self, arg0: bool) -> None: ...

class FromSurfaceMesh:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: ...) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: ..., arg2: bool) -> None: ...
    def apply(self) -> None:
        """
        Construct a levelset from a surface mesh.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset to read into.
        """

    def setMesh(self, arg0: ...) -> None:
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
    def __init__(self, arg0: ..., arg1: ...) -> None: ...
    @typing.overload
    def __init__(self, arg0: ..., arg1: ..., arg2: bool) -> None: ...
    def apply(self) -> None:
        """
        Construct a levelset from a volume mesh.
        """

    def setGrid(self, arg0: ...) -> None:
        """
        Set the grid used to read in the level sets.
        """

    def setMesh(self, arg0: ...) -> None:
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
    @typing.overload
    def __init__(self) -> None: ...
    def apply(self) -> None:
        """
        Perform advection.
        """

    def setAdvectionDistribution(self, arg0: PylsGeometricAdvectDistribution) -> None:
        """
        Set advection distribution to use as kernel for the fast advection.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset to advect.
        """

class GeometricAdvectDistribution:
    def __init__(self) -> None: ...
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

class IntegrationSchemeEnum:
    """
    Members:

      ENGQUIST_OSHER_1ST_ORDER

      ENGQUIST_OSHER_2ND_ORDER

      LAX_FRIEDRICHS_1ST_ORDER

      LAX_FRIEDRICHS_2ND_ORDER

      LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER

      LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER

      LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER

      LOCAL_LAX_FRIEDRICHS_1ST_ORDER

      LOCAL_LAX_FRIEDRICHS_2ND_ORDER

      STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER
    """

    ENGQUIST_OSHER_1ST_ORDER: typing.ClassVar[
        IntegrationSchemeEnum
    ]  # value = <IntegrationSchemeEnum.ENGQUIST_OSHER_1ST_ORDER: 0>
    ENGQUIST_OSHER_2ND_ORDER: typing.ClassVar[
        IntegrationSchemeEnum
    ]  # value = <IntegrationSchemeEnum.ENGQUIST_OSHER_2ND_ORDER: 1>
    LAX_FRIEDRICHS_1ST_ORDER: typing.ClassVar[
        IntegrationSchemeEnum
    ]  # value = <IntegrationSchemeEnum.LAX_FRIEDRICHS_1ST_ORDER: 2>
    LAX_FRIEDRICHS_2ND_ORDER: typing.ClassVar[
        IntegrationSchemeEnum
    ]  # value = <IntegrationSchemeEnum.LAX_FRIEDRICHS_2ND_ORDER: 3>
    LOCAL_LAX_FRIEDRICHS_1ST_ORDER: typing.ClassVar[
        IntegrationSchemeEnum
    ]  # value = <IntegrationSchemeEnum.LOCAL_LAX_FRIEDRICHS_1ST_ORDER: 7>
    LOCAL_LAX_FRIEDRICHS_2ND_ORDER: typing.ClassVar[
        IntegrationSchemeEnum
    ]  # value = <IntegrationSchemeEnum.LOCAL_LAX_FRIEDRICHS_2ND_ORDER: 8>
    LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER: typing.ClassVar[
        IntegrationSchemeEnum
    ]  # value = <IntegrationSchemeEnum.LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER: 4>
    LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER: typing.ClassVar[
        IntegrationSchemeEnum
    ]  # value = <IntegrationSchemeEnum.LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER: 5>
    LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER: typing.ClassVar[
        IntegrationSchemeEnum
    ]  # value = <IntegrationSchemeEnum.LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER: 6>
    STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER: typing.ClassVar[
        IntegrationSchemeEnum
    ]  # value = <IntegrationSchemeEnum.STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER: 9>
    __members__: typing.ClassVar[
        dict[str, IntegrationSchemeEnum]
    ]  # value = {'ENGQUIST_OSHER_1ST_ORDER': <IntegrationSchemeEnum.ENGQUIST_OSHER_1ST_ORDER: 0>, 'ENGQUIST_OSHER_2ND_ORDER': <IntegrationSchemeEnum.ENGQUIST_OSHER_2ND_ORDER: 1>, 'LAX_FRIEDRICHS_1ST_ORDER': <IntegrationSchemeEnum.LAX_FRIEDRICHS_1ST_ORDER: 2>, 'LAX_FRIEDRICHS_2ND_ORDER': <IntegrationSchemeEnum.LAX_FRIEDRICHS_2ND_ORDER: 3>, 'LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER': <IntegrationSchemeEnum.LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER: 4>, 'LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER': <IntegrationSchemeEnum.LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER: 5>, 'LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER': <IntegrationSchemeEnum.LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER: 6>, 'LOCAL_LAX_FRIEDRICHS_1ST_ORDER': <IntegrationSchemeEnum.LOCAL_LAX_FRIEDRICHS_1ST_ORDER: 7>, 'LOCAL_LAX_FRIEDRICHS_2ND_ORDER': <IntegrationSchemeEnum.LOCAL_LAX_FRIEDRICHS_2ND_ORDER: 8>, 'STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER': <IntegrationSchemeEnum.STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER: 9>}
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: typing.SupportsInt) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class LogLevel:
    """
    Members:

      ERROR

      WARNING

      INFO

      INTERMEDIATE

      TIMING

      DEBUG
    """

    DEBUG: typing.ClassVar[LogLevel]  # value = <LogLevel.DEBUG: 5>
    ERROR: typing.ClassVar[LogLevel]  # value = <LogLevel.ERROR: 0>
    INFO: typing.ClassVar[LogLevel]  # value = <LogLevel.INFO: 2>
    INTERMEDIATE: typing.ClassVar[LogLevel]  # value = <LogLevel.INTERMEDIATE: 3>
    TIMING: typing.ClassVar[LogLevel]  # value = <LogLevel.TIMING: 4>
    WARNING: typing.ClassVar[LogLevel]  # value = <LogLevel.WARNING: 1>
    __members__: typing.ClassVar[
        dict[str, LogLevel]
    ]  # value = {'ERROR': <LogLevel.ERROR: 0>, 'WARNING': <LogLevel.WARNING: 1>, 'INFO': <LogLevel.INFO: 2>, 'INTERMEDIATE': <LogLevel.INTERMEDIATE: 3>, 'TIMING': <LogLevel.TIMING: 4>, 'DEBUG': <LogLevel.DEBUG: 5>}
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: typing.SupportsInt) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Logger:
    @staticmethod
    def appendToLogFile(arg0: str) -> bool: ...
    @staticmethod
    def closeLogFile() -> None: ...
    @staticmethod
    def getInstance() -> Logger: ...
    @staticmethod
    def getLogLevel() -> int: ...
    @staticmethod
    def setLogFile(arg0: str) -> bool: ...
    @staticmethod
    def setLogLevel(arg0: LogLevel) -> None: ...
    def addDebug(self, arg0: str) -> Logger: ...
    def addError(self, s: str, shouldAbort: bool = True) -> Logger: ...
    def addInfo(self, arg0: str) -> Logger: ...
    @typing.overload
    def addTiming(self, arg0: str, arg1: typing.SupportsFloat) -> Logger: ...
    @typing.overload
    def addTiming(
        self, arg0: str, arg1: typing.SupportsFloat, arg2: typing.SupportsFloat
    ) -> Logger: ...
    def addWarning(self, arg0: str) -> Logger: ...
    def print(self) -> None: ...

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

    def setVoidTopSurface(self, arg0: ...) -> None:
        """
        Set the logic by which to choose the surface which is non-void. All other connected surfaces will then be marked as void points.
        """

class MaterialMap:
    def __init__(self) -> None: ...
    def getMaterialId(self, arg0: typing.SupportsInt) -> int: ...
    def getNumberOfLayers(self) -> int:
        """
        Get the number of level-sets in the material map.
        """

    def getNumberOfMaterials(self) -> int: ...
    def insertNextMaterial(self, arg0: typing.SupportsInt) -> None:
        """
        Insert a new material into the map.
        """

    def setMaterialId(
        self, arg0: typing.SupportsInt, arg1: typing.SupportsInt
    ) -> None: ...

class Mesh:
    def __init__(self) -> None: ...
    def append(self, arg0: Mesh) -> None:
        """
        Append another mesh to this mesh.
        """

    def clear(self) -> None:
        """
        Clear all data in the mesh.
        """

    def getCellData(self) -> PointData:
        """
        Return a reference to the cell data of the mesh.
        """

    def getHexas(self) -> list[typing.Annotated[list[int], "FixedSize(8)"]]:
        """
        Get a list of hexahedrons of the mesh.
        """

    def getLines(self) -> list[typing.Annotated[list[int], "FixedSize(2)"]]:
        """
        Get a list of lines of the mesh.
        """

    def getNodes(self) -> list[typing.Annotated[list[float], "FixedSize(3)"]]:
        """
        Get all nodes of the mesh as a list.
        """

    def getPointData(self) -> PointData:
        """
        Return a reference to the point data of the mesh.
        """

    def getTetras(self) -> list[typing.Annotated[list[int], "FixedSize(4)"]]:
        """
        Get a list of tetrahedrons of the mesh.
        """

    def getTriangles(self) -> list[typing.Annotated[list[int], "FixedSize(3)"]]:
        """
        Get a list of triangles of the mesh.
        """

    def getVerticies(self) -> list[typing.Annotated[list[int], "FixedSize(1)"]]:
        """
        Get a list of verticies of the mesh.
        """

    def insertNextHexa(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsInt], "FixedSize(8)"
        ],
    ) -> int:
        """
        Insert a hexahedron in the mesh.
        """

    def insertNextLine(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsInt], "FixedSize(2)"
        ],
    ) -> int:
        """
        Insert a line in the mesh.
        """

    def insertNextNode(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
    ) -> int:
        """
        Insert a node in the mesh.
        """

    def insertNextTetra(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsInt], "FixedSize(4)"
        ],
    ) -> int:
        """
        Insert a tetrahedron in the mesh.
        """

    def insertNextTriangle(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsInt], "FixedSize(3)"
        ],
    ) -> int:
        """
        Insert a triangle in the mesh.
        """

    def insertNextVertex(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsInt], "FixedSize(1)"
        ],
    ) -> int:
        """
        Insert a vertex in the mesh.
        """

    def print(self) -> None:
        """
        Print basic information about the mesh.
        """

    def removeDuplicateNodes(self) -> None:
        """
        Remove nodes which occur twice in the mesh, and replace their IDs in the mesh elements.
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

class PointData:
    def __init__(self) -> None: ...
    @typing.overload
    def getScalarData(self, arg0: typing.SupportsInt) -> list[float]: ...
    @typing.overload
    def getScalarData(self, arg0: str, arg1: bool) -> list[float]: ...
    def getScalarDataLabel(self, arg0: typing.SupportsInt) -> str: ...
    def getScalarDataSize(self) -> int: ...
    @typing.overload
    def getVectorData(
        self, arg0: typing.SupportsInt
    ) -> list[typing.Annotated[list[float], "FixedSize(3)"]]: ...
    @typing.overload
    def getVectorData(
        self, arg0: str, arg1: bool
    ) -> list[typing.Annotated[list[float], "FixedSize(3)"]]: ...
    def getVectorDataLabel(self, arg0: typing.SupportsInt) -> str: ...
    def getVectorDataSize(self) -> int: ...
    def insertNextScalarData(
        self,
        scalars: collections.abc.Sequence[typing.SupportsFloat],
        label: str = "Scalars",
    ) -> None: ...
    def insertNextVectorData(
        self,
        vectors: collections.abc.Sequence[
            typing.Annotated[
                collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
            ]
        ],
        label: str = "Vectors",
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

    def setVoidTopSurface(self, arg0: VoidTopSurfaceEnum) -> None:
        """
        Set the logic by which to choose the surface which should be kept. All other LS values will be marked as stray points and removed.
        """

class Slice:
    @staticmethod
    @typing.overload
    def __init__(*args, **kwargs) -> None: ...
    @staticmethod
    def setSliceLevelSet(*args, **kwargs) -> None:
        """
        Set the 2D level set where the extracted slice will be stored.
        """

    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(
        self, arg0: Domain, arg1: typing.SupportsInt, arg2: typing.SupportsFloat
    ) -> None: ...
    def apply(self) -> None:
        """
        Extract the 2D slice from the 3D domain.
        """

    def getSliceLevelSet(self) -> ...:
        """
        Get the 2D slice level set after extraction.
        """

    def setReflectX(self, arg0: bool) -> None:
        """
        Set whether to reflect all x-coordinates in the resulting slice.
        """

    def setSliceDimension(self, arg0: typing.SupportsInt) -> None:
        """
        Set the dimension along which to slice (0=x, 1=y, 2=z).
        """

    def setSlicePosition(self, arg0: typing.SupportsFloat) -> None:
        """
        Set the position along the slice dimension where to extract the slice.
        """

    def setSourceLevelSet(self, arg0: Domain) -> None:
        """
        Set the 3D source level set from which to extract the slice.
        """

    def setWritePath(self, arg0: str, arg1: bool) -> None:
        """
        Set the path where the slice should be written to.
        """

class Sphere:
    def __init__(
        self,
        origin: collections.abc.Sequence[typing.SupportsFloat],
        radius: typing.SupportsFloat,
    ) -> None: ...

class SphereDistribution(GeometricAdvectDistribution):
    def __init__(
        self, arg0: typing.SupportsFloat, arg1: typing.SupportsFloat
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

class StencilLocalLaxFriedrichsScalar:
    @staticmethod
    def setMaxDissipation(maxDissipation: typing.SupportsFloat) -> None: ...

class ToDiskMesh:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: Mesh) -> None: ...
    def apply(self) -> None:
        """
        Convert the levelset to a surface mesh.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset to mesh.
        """

    def setMesh(self, arg0: Mesh) -> None:
        """
        Set the mesh to generate.
        """

class ToMesh:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: Mesh) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: Mesh, arg2: bool) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: Mesh, arg2: bool, arg3: bool) -> None: ...
    def apply(self) -> None:
        """
        Convert the levelset to a surface mesh.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset to mesh.
        """

    def setMesh(self, arg0: Mesh) -> None:
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
    def __init__(self, arg0: Domain, arg1: Mesh) -> None: ...
    def apply(self) -> None:
        """
        Convert the levelset to a surface mesh.
        """

    def setLevelSet(self, arg0: Domain) -> None:
        """
        Set levelset to mesh.
        """

    def setMesh(self, arg0: Mesh) -> None:
        """
        Set the mesh to generate.
        """

class ToVoxelMesh:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Mesh) -> None: ...
    @typing.overload
    def __init__(self, arg0: Domain, arg1: Mesh) -> None: ...
    @typing.overload
    def __init__(self, arg0: collections.abc.Sequence[Domain], arg1: Mesh) -> None: ...
    def apply(self) -> None:
        """
        Convert the levelset to a surface mesh.
        """

    def insertNextLevelSet(self, arg0: Domain) -> None:
        """
        Insert next level set to output in the mesh.
        """

    def setMesh(self, arg0: Mesh) -> None:
        """
        Set the mesh to generate.
        """

class TransformEnum:
    """
    Members:

      TRANSLATION

      ROTATION

      SCALE
    """

    ROTATION: typing.ClassVar[TransformEnum]  # value = <TransformEnum.ROTATION: 1>
    SCALE: typing.ClassVar[TransformEnum]  # value = <TransformEnum.SCALE: 2>
    TRANSLATION: typing.ClassVar[
        TransformEnum
    ]  # value = <TransformEnum.TRANSLATION: 0>
    __members__: typing.ClassVar[
        dict[str, TransformEnum]
    ]  # value = {'TRANSLATION': <TransformEnum.TRANSLATION: 0>, 'ROTATION': <TransformEnum.ROTATION: 1>, 'SCALE': <TransformEnum.SCALE: 2>}
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: typing.SupportsInt) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class TransformMesh:
    def __init__(
        self,
        mesh: Mesh,
        transform: TransformEnum = ...,
        transformVector: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ] = [0.0, 0.0, 0.0],
        angle: typing.SupportsFloat = 0.0,
    ) -> None: ...
    def apply(self) -> None:
        """
        Apply the transformation.
        """

class VTKReader:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Mesh) -> None: ...
    @typing.overload
    def __init__(self, arg0: Mesh, arg1: str) -> None: ...
    @typing.overload
    def __init__(self, arg0: Mesh, arg1: FileFormatEnum, arg2: str) -> None: ...
    def apply(self) -> None:
        """
        Read the mesh.
        """

    def getMetaData(self) -> dict[str, list[float]]:
        """
        Get the metadata from the file.
        """

    def setFileFormat(self, arg0: FileFormatEnum) -> None:
        """
        Set the file format of the file to be read.
        """

    def setFileName(self, arg0: str) -> None:
        """
        Set the name of the input file.
        """

    def setMesh(self, arg0: Mesh) -> None:
        """
        Set the mesh to read into.
        """

class VTKWriter:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: Mesh) -> None: ...
    @typing.overload
    def __init__(self, arg0: Mesh, arg1: str) -> None: ...
    @typing.overload
    def __init__(self, arg0: Mesh, arg1: FileFormatEnum, arg2: str) -> None: ...
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
        Write the mesh.
        """

    def setFileFormat(self, arg0: FileFormatEnum) -> None:
        """
        Set the file format, the mesh should be written to.
        """

    def setFileName(self, arg0: str) -> None:
        """
        Set the name of the output file.
        """

    def setMesh(self, arg0: Mesh) -> None:
        """
        Set the mesh to output.
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

class VelocityField:
    def __init__(self) -> None: ...
    def getDissipationAlpha(
        self,
        arg0: typing.SupportsInt,
        arg1: typing.SupportsInt,
        arg2: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
    ) -> float:
        """
        Return the analytical dissipation alpha value if the lsLocalLaxFriedrichsAnalytical scheme is used for advection.
        """

    def getScalarVelocity(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg1: typing.SupportsInt,
        arg2: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg3: typing.SupportsInt,
    ) -> float:
        """
        Return the scalar velocity for a point of material at coordinate with normal vector normal.
        """

    def getVectorVelocity(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg1: typing.SupportsInt,
        arg2: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(3)"
        ],
        arg3: typing.SupportsInt,
    ) -> typing.Annotated[list[float], "FixedSize(3)"]:
        """
        Return the vector velocity for a point of material at coordinate with normal vector normal.
        """

class VoidTopSurfaceEnum:
    """
    Members:

      LEX_LOWEST

      LEX_HIGHEST

      LARGEST

      SMALLEST
    """

    LARGEST: typing.ClassVar[
        VoidTopSurfaceEnum
    ]  # value = <VoidTopSurfaceEnum.LARGEST: 2>
    LEX_HIGHEST: typing.ClassVar[
        VoidTopSurfaceEnum
    ]  # value = <VoidTopSurfaceEnum.LEX_HIGHEST: 1>
    LEX_LOWEST: typing.ClassVar[
        VoidTopSurfaceEnum
    ]  # value = <VoidTopSurfaceEnum.LEX_LOWEST: 0>
    SMALLEST: typing.ClassVar[
        VoidTopSurfaceEnum
    ]  # value = <VoidTopSurfaceEnum.SMALLEST: 3>
    __members__: typing.ClassVar[
        dict[str, VoidTopSurfaceEnum]
    ]  # value = {'LEX_LOWEST': <VoidTopSurfaceEnum.LEX_LOWEST: 0>, 'LEX_HIGHEST': <VoidTopSurfaceEnum.LEX_HIGHEST: 1>, 'LARGEST': <VoidTopSurfaceEnum.LARGEST: 2>, 'SMALLEST': <VoidTopSurfaceEnum.SMALLEST: 3>}
    def __eq__(self, other: typing.Any) -> bool: ...
    def __getstate__(self) -> int: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __init__(self, value: typing.SupportsInt) -> None: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: typing.Any) -> bool: ...
    def __repr__(self) -> str: ...
    def __setstate__(self, state: typing.SupportsInt) -> None: ...
    def __str__(self) -> str: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

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

def setNumThreads(arg0: typing.SupportsInt) -> None: ...

__version__: str = "4.5.0"
version: str = "4.5.0"
