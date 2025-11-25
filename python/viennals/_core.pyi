"""
ViennaLS is a header-only C++ level set library developed for high performance topography simulations. The main design goals are simplicity and efficiency, tailored towards scientific simulations. ViennaLS can also be used for visualization applications, although this is not the main design target.
"""

from __future__ import annotations
import collections.abc
import enum
import typing
import viennals.d2
from viennals import d2
import viennals.d3
from viennals import d3

__all__: list[str] = [
    "BooleanOperationEnum",
    "BoundaryConditionEnum",
    "CurvatureEnum",
    "Extrude",
    "FeatureDetectionEnum",
    "FileFormatEnum",
    "IntegrationSchemeEnum",
    "LogLevel",
    "Logger",
    "MaterialMap",
    "Mesh",
    "PointData",
    "Slice",
    "TransformEnum",
    "TransformMesh",
    "VTKReader",
    "VTKWriter",
    "VelocityField",
    "VoidTopSurfaceEnum",
    "d2",
    "d3",
    "setNumThreads",
    "version",
]

class BooleanOperationEnum(enum.IntEnum):
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
    @classmethod
    def __new__(cls, value): ...
    def __format__(self, format_spec):
        """
        Convert to a string according to format_spec.
        """

class BoundaryConditionEnum(enum.IntEnum):
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
    @classmethod
    def __new__(cls, value): ...
    def __format__(self, format_spec):
        """
        Convert to a string according to format_spec.
        """

class CurvatureEnum(enum.IntEnum):
    GAUSSIAN_CURVATURE: typing.ClassVar[
        CurvatureEnum
    ]  # value = <CurvatureEnum.GAUSSIAN_CURVATURE: 1>
    MEAN_AND_GAUSSIAN_CURVATURE: typing.ClassVar[
        CurvatureEnum
    ]  # value = <CurvatureEnum.MEAN_AND_GAUSSIAN_CURVATURE: 2>
    MEAN_CURVATURE: typing.ClassVar[
        CurvatureEnum
    ]  # value = <CurvatureEnum.MEAN_CURVATURE: 0>
    @classmethod
    def __new__(cls, value): ...
    def __format__(self, format_spec):
        """
        Convert to a string according to format_spec.
        """

class Extrude:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(
        self,
        arg0: viennals.d2.Domain,
        arg1: viennals.d3.Domain,
        arg2: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(2)"
        ],
        arg3: typing.SupportsInt,
        arg4: typing.Annotated[
            collections.abc.Sequence[BoundaryConditionEnum], "FixedSize(3)"
        ],
    ) -> None: ...
    def apply(self) -> None:
        """
        Perform extrusion.
        """

    @typing.overload
    def setBoundaryConditions(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[BoundaryConditionEnum], "FixedSize(3)"
        ],
    ) -> None:
        """
        Set the boundary conditions in the 3D extruded domain.
        """

    @typing.overload
    def setBoundaryConditions(self, arg0: BoundaryConditionEnum) -> None:
        """
        Set the boundary conditions in the 3D extruded domain.
        """

    def setExtent(
        self,
        arg0: typing.Annotated[
            collections.abc.Sequence[typing.SupportsFloat], "FixedSize(2)"
        ],
    ) -> None:
        """
        Set the extent in the extruded dimension
        """

    def setExtrusionAxis(self, arg0: typing.SupportsInt) -> None:
        """
        Set the axis in which to extrude (0=x, 1=y, 2=z).
        """

    def setInputLevelSet(self, arg0: viennals.d2.Domain) -> None:
        """
        Set 2D input Level Set
        """

    def setOutputLevelSet(self, arg0: viennals.d3.Domain) -> None:
        """
        Set 3D output Level Set
        """

class FeatureDetectionEnum(enum.IntEnum):
    CURVATURE: typing.ClassVar[
        FeatureDetectionEnum
    ]  # value = <FeatureDetectionEnum.CURVATURE: 0>
    NORMALS_ANGLE: typing.ClassVar[
        FeatureDetectionEnum
    ]  # value = <FeatureDetectionEnum.NORMALS_ANGLE: 1>
    @classmethod
    def __new__(cls, value): ...
    def __format__(self, format_spec):
        """
        Convert to a string according to format_spec.
        """

class FileFormatEnum(enum.IntEnum):
    VTK_LEGACY: typing.ClassVar[
        FileFormatEnum
    ]  # value = <FileFormatEnum.VTK_LEGACY: 0>
    VTP: typing.ClassVar[FileFormatEnum]  # value = <FileFormatEnum.VTP: 1>
    VTU: typing.ClassVar[FileFormatEnum]  # value = <FileFormatEnum.VTU: 2>
    @classmethod
    def __new__(cls, value): ...
    def __format__(self, format_spec):
        """
        Convert to a string according to format_spec.
        """

class IntegrationSchemeEnum(enum.IntEnum):
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
    @classmethod
    def __new__(cls, value): ...
    def __format__(self, format_spec):
        """
        Convert to a string according to format_spec.
        """

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

class MaterialMap:
    def __init__(self) -> None: ...
    def clear(self) -> None: ...
    def getMaterialId(self, arg0: typing.SupportsInt) -> int: ...
    def getMaterialMap(self) -> list[int]:
        """
        Get the material map.
        """

    def getMaterials(self) -> set[int]:
        """
        Get a list of all materials in the map.
        """

    def getNumberOfLayers(self) -> int:
        """
        Get the number of level-sets in the material map.
        """

    def getNumberOfMaterials(self) -> int: ...
    def hasMaterial(self, arg0: typing.SupportsInt) -> bool: ...
    def insertNextMaterial(self, arg0: typing.SupportsInt) -> None:
        """
        Insert a new material into the map.
        """

    def isValidIndex(self, arg0: typing.SupportsInt) -> bool: ...
    def reserve(self, arg0: typing.SupportsInt) -> None: ...
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

class Slice:
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(
        self,
        arg0: viennals.d3.Domain,
        arg1: viennals.d2.Domain,
        arg2: typing.SupportsInt,
        arg3: typing.SupportsFloat,
    ) -> None: ...
    @typing.overload
    def __init__(
        self,
        arg0: viennals.d3.Domain,
        arg1: typing.SupportsInt,
        arg2: typing.SupportsFloat,
    ) -> None: ...
    def apply(self) -> None:
        """
        Extract the 2D slice from the 3D domain.
        """

    def getSliceLevelSet(self) -> viennals.d2.Domain:
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

    def setSliceLevelSet(self, arg0: viennals.d2.Domain) -> None:
        """
        Set the 2D level set where the extracted slice will be stored.
        """

    def setSlicePosition(self, arg0: typing.SupportsFloat) -> None:
        """
        Set the position along the slice dimension where to extract the slice.
        """

    def setSourceLevelSet(self, arg0: viennals.d3.Domain) -> None:
        """
        Set the 3D source level set from which to extract the slice.
        """

    def setWritePath(self, arg0: str, arg1: bool) -> None:
        """
        Set the path where the slice should be written to.
        """

class TransformEnum(enum.IntEnum):
    ROTATION: typing.ClassVar[TransformEnum]  # value = <TransformEnum.ROTATION: 1>
    SCALE: typing.ClassVar[TransformEnum]  # value = <TransformEnum.SCALE: 2>
    TRANSLATION: typing.ClassVar[
        TransformEnum
    ]  # value = <TransformEnum.TRANSLATION: 0>
    @classmethod
    def __new__(cls, value): ...
    def __format__(self, format_spec):
        """
        Convert to a string according to format_spec.
        """

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

class VoidTopSurfaceEnum(enum.IntEnum):
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
    @classmethod
    def __new__(cls, value): ...
    def __format__(self, format_spec):
        """
        Convert to a string according to format_spec.
        """

def setNumThreads(arg0: typing.SupportsInt) -> None: ...

__version__: str = "5.2.0"
version: str = "5.2.0"
