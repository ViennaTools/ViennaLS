from typing import ClassVar, overload

__version__: str
version: str

class Advect:
    def __init__(self) -> None: ...
    def apply(self) -> None: ...
    def applyIntegration(self, arg0: float) -> None: ...
    def getAdvectedTime(self) -> float: ...
    def getCalculateNormalVectors(self) -> bool: ...
    def getCurrentTimeStep(self) -> float: ...
    def getNumberOfTimeSteps(self) -> int: ...
    def getTimeStepRatio(self) -> float: ...
    def insertNextLevelSet(self, *args, **kwargs): ...
    def prepareLS(self) -> None: ...
    def setAdvectionTime(self, arg0: float) -> None: ...
    def setCalculateNormalVectors(self, arg0: bool) -> None: ...
    def setDissipationAlpha(self, arg0: float) -> None: ...
    def setIgnoreVoids(self, arg0: bool) -> None: ...
    def setIntegrationScheme(self, arg0) -> None: ...
    def setTimeStepRatio(self, arg0: float) -> None: ...
    def setVelocityField(self, arg0) -> None: ...

class BooleanOperation:
    def __init__(self) -> None: ...
    def apply(self) -> None: ...
    def setBooleanOperation(self, arg0) -> None: ...
    def setLevelset(self, *args, **kwargs): ...
    def setSecondLevelSet(self, *args, **kwargs): ...

class BooleanOperationEnum:
    __members__: ClassVar[dict] = ...  # read-only
    INTERSECT: ClassVar[BooleanOperationEnum] = ...
    INVERT: ClassVar[BooleanOperationEnum] = ...
    RELATIVE_COMPLEMENT: ClassVar[BooleanOperationEnum] = ...
    UNION: ClassVar[BooleanOperationEnum] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class BoundaryConditionEnum:
    __members__: ClassVar[dict] = ...  # read-only
    INFINITE_BOUNDARY: ClassVar[BoundaryConditionEnum] = ...
    NEG_INFINITE_BOUNDARY: ClassVar[BoundaryConditionEnum] = ...
    PERIODIC_BOUNDARY: ClassVar[BoundaryConditionEnum] = ...
    POS_INFINITE_BOUNDARY: ClassVar[BoundaryConditionEnum] = ...
    REFLECTIVE_BOUNDARY: ClassVar[BoundaryConditionEnum] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Box:
    def __init__(self, minPoint: list[float], maxPoint: list[float]) -> None: ...

class BoxDistribution(GeometricAdvectDistribution):
    def __init__(self, arg0, arg1: float) -> None: ...
    def getBounds(self, *args, **kwargs): ...
    def getSignedDistance(self, arg0, arg1, arg2: int) -> float: ...
    def isInside(self, arg0, arg1, arg2: float) -> bool: ...

class CalculateCurvatures:
    def __init__(self) -> None: ...
    def apply(self) -> None: ...
    def setCurvatureType(self, arg0) -> None: ...
    def setLevelSet(self, *args, **kwargs): ...
    def setMaxValue(self, arg0: float) -> None: ...

class CalculateNormalVectors:
    def __init__(self) -> None: ...
    def apply(self) -> None: ...
    def setLevelSet(self, *args, **kwargs): ...

class CalculateVisibilities:
    def __init__(self, *args, **kwargs) -> None: ...
    def apply(self) -> None: ...

class Check:
    def __init__(self) -> None: ...
    def apply(self) -> None: ...
    def setLevelSet(self, *args, **kwargs): ...

class ConvexHull:
    def __init__(self) -> None: ...
    def apply(self) -> None: ...
    def setMesh(self, arg0) -> None: ...
    def setPointCloud(self, *args, **kwargs): ...

class CurvatureEnum:
    __members__: ClassVar[dict] = ...  # read-only
    GAUSSIAN_CURVATURE: ClassVar[CurvatureEnum] = ...
    MEAN_AND_GAUSSIAN_CURVATURE: ClassVar[CurvatureEnum] = ...
    MEAN_CURVATURE: ClassVar[CurvatureEnum] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Cylinder:
    def __init__(
        self,
        origin: list[float],
        axisDirection: list[float],
        height: float,
        radius: float,
        topRadius: float = ...,
    ) -> None: ...

class DetectFeatures:
    def __init__(self) -> None: ...
    def apply(self) -> None: ...
    def setDetectionMethod(self, arg0) -> None: ...
    def setDetectionThreshold(self, arg0: float) -> None: ...

class Domain:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, gridDelta: float = ...) -> None: ...
    @overload
    def __init__(self, bounds, boundaryConditions, gridDelta: float = ...) -> None: ...
    @overload
    def __init__(
        self, bounds: list[float], boundaryConditions: list[int], gridDelta: float = ...
    ) -> None: ...
    @overload
    def __init__(self, arg0: Domain) -> None: ...
    @overload
    def __init__(self, arg0) -> None: ...
    def clearMetaData(self) -> None: ...
    def deepCopy(self, arg0: Domain) -> None: ...
    def getLevelSetWidth(self) -> int: ...
    def getNumberOfPoints(self) -> int: ...
    def getNumberOfSegments(self) -> int: ...
    def print(self, stream: object = ...) -> None: ...
    def setLevelSetWidth(self, arg0: int) -> None: ...

class Expand:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Domain) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: int) -> None: ...
    def apply(self) -> None: ...
    def setLevelSet(self, arg0: Domain) -> None: ...
    def setWidth(self, arg0: int) -> None: ...

class FeatureDetectionEnum:
    __members__: ClassVar[dict] = ...  # read-only
    CURVATURE: ClassVar[FeatureDetectionEnum] = ...
    NORMALS_ANGLE: ClassVar[FeatureDetectionEnum] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class FileFormatEnum:
    __members__: ClassVar[dict] = ...  # read-only
    VTK_LEGACY: ClassVar[FileFormatEnum] = ...
    VTP: ClassVar[FileFormatEnum] = ...
    VTU: ClassVar[FileFormatEnum] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class FromMesh:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1) -> None: ...
    def apply(self) -> None: ...
    def setMesh(self, arg0) -> None: ...
    def setSortPointList(self, arg0: bool) -> None: ...

class FromSurfaceMesh:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1, arg2: bool) -> None: ...
    def apply(self) -> None: ...
    def setLevelSet(self, arg0: Domain) -> None: ...
    def setMesh(self, arg0) -> None: ...
    @overload
    def setRemoveBoundaryTriangles(self, arg0: bool) -> None: ...
    @overload
    def setRemoveBoundaryTriangles(self, arg0) -> None: ...

class FromVolumeMesh:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0, arg1) -> None: ...
    @overload
    def __init__(self, arg0, arg1, arg2: bool) -> None: ...
    def apply(self) -> None: ...
    def setGrid(self, arg0) -> None: ...
    def setMesh(self, arg0) -> None: ...
    def setRemoveBoundaryTriangles(self, arg0: bool) -> None: ...

class GeometricAdvect:
    def __init__(self) -> None: ...
    def apply(self) -> None: ...
    def setAdvectionDistribution(
        self, arg0: PylsGeometricAdvectDistribution
    ) -> None: ...
    def setLevelSet(self, arg0: Domain) -> None: ...

class GeometricAdvectDistribution:
    def __init__(self) -> None: ...
    def getBounds(self, *args, **kwargs): ...
    def getSignedDistance(self, arg0, arg1, arg2: int) -> float: ...
    def isInside(self, arg0, arg1, arg2: float) -> bool: ...

class IntegrationSchemeEnum:
    __members__: ClassVar[dict] = ...  # read-only
    ENGQUIST_OSHER_1ST_ORDER: ClassVar[IntegrationSchemeEnum] = ...
    ENGQUIST_OSHER_2ND_ORDER: ClassVar[IntegrationSchemeEnum] = ...
    LAX_FRIEDRICHS_1ST_ORDER: ClassVar[IntegrationSchemeEnum] = ...
    LAX_FRIEDRICHS_2ND_ORDER: ClassVar[IntegrationSchemeEnum] = ...
    LOCAL_LAX_FRIEDRICHS_1ST_ORDER: ClassVar[IntegrationSchemeEnum] = ...
    LOCAL_LAX_FRIEDRICHS_2ND_ORDER: ClassVar[IntegrationSchemeEnum] = ...
    LOCAL_LAX_FRIEDRICHS_ANALYTICAL_1ST_ORDER: ClassVar[IntegrationSchemeEnum] = ...
    LOCAL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER: ClassVar[IntegrationSchemeEnum] = ...
    LOCAL_LOCAL_LAX_FRIEDRICHS_2ND_ORDER: ClassVar[IntegrationSchemeEnum] = ...
    STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER: ClassVar[IntegrationSchemeEnum] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class LogLevel:
    __members__: ClassVar[dict] = ...  # read-only
    DEBUG: ClassVar[LogLevel] = ...
    ERROR: ClassVar[LogLevel] = ...
    INFO: ClassVar[LogLevel] = ...
    INTERMEDIATE: ClassVar[LogLevel] = ...
    TIMING: ClassVar[LogLevel] = ...
    WARNING: ClassVar[LogLevel] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class Logger:
    def __init__(self, *args, **kwargs) -> None: ...
    def addDebug(self, arg0: str) -> Logger: ...
    def addError(self, s: str, shouldAbort: bool = ...) -> Logger: ...
    def addInfo(self, arg0: str) -> Logger: ...
    @overload
    def addTiming(self, arg0: str, arg1: float) -> Logger: ...
    @overload
    def addTiming(self, arg0: str, arg1: float, arg2: float) -> Logger: ...
    def addWarning(self, arg0: str) -> Logger: ...
    @staticmethod
    def getInstance() -> Logger: ...
    @staticmethod
    def getLogLevel() -> int: ...
    def print(self) -> None: ...
    @staticmethod
    def setLogLevel(arg0: LogLevel) -> None: ...

class MakeGeometry:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Domain) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: Sphere) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: Plane) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: Box) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: Cylinder) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: PointCloud) -> None: ...
    def apply(self) -> None: ...
    @overload
    def setGeometry(self, arg0: Sphere) -> None: ...
    @overload
    def setGeometry(self, arg0: Plane) -> None: ...
    @overload
    def setGeometry(self, arg0: Box) -> None: ...
    @overload
    def setGeometry(self, arg0: Cylinder) -> None: ...
    @overload
    def setGeometry(self, arg0: PointCloud) -> None: ...
    @overload
    def setIgnoreBoundaryConditions(self, arg0: bool) -> None: ...
    @overload
    def setIgnoreBoundaryConditions(self, arg0) -> None: ...
    def setLevelSet(self, arg0: Domain) -> None: ...

class MarkVoidPoints:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Domain) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: bool) -> None: ...
    def apply(self) -> None: ...
    def setDetectLargestSurface(self, arg0: bool) -> None: ...
    def setLevelSet(self, arg0: Domain) -> None: ...
    def setReverseVoidDetection(self, arg0: bool) -> None: ...
    def setSaveComponentsId(self, arg0: bool) -> None: ...
    def setVoidTopSurface(self, arg0) -> None: ...

class MaterialMap:
    def __init__(self) -> None: ...
    def getMaterialId(self, arg0: int) -> int: ...
    def getNumberOfLayers(self) -> int: ...
    def getNumberOfMaterials(self) -> int: ...
    def insertNextMaterial(self, arg0: int) -> None: ...
    def setMaterialId(self, arg0: int, arg1: int) -> None: ...

class Mesh:
    def __init__(self) -> None: ...
    def append(self, arg0: Mesh) -> None: ...
    def clear(self) -> None: ...
    def getCellData(self) -> PointData: ...
    def getHexas(self, *args, **kwargs): ...
    def getLines(self, *args, **kwargs): ...
    def getNodes(self, *args, **kwargs): ...
    def getPointData(self) -> PointData: ...
    def getTetras(self, *args, **kwargs): ...
    def getTriangles(self, *args, **kwargs): ...
    def getVerticies(self, *args, **kwargs): ...
    def insertNextHexa(self, arg0) -> int: ...
    def insertNextLine(self, arg0) -> int: ...
    def insertNextNode(self, arg0) -> int: ...
    def insertNextTetra(self, arg0) -> int: ...
    def insertNextTriangle(self, arg0) -> int: ...
    def insertNextVertex(self, arg0) -> int: ...
    def print(self) -> None: ...
    def removeDuplicateNodes(self) -> None: ...

class Plane:
    def __init__(self, origin: list[float], normal: list[float]) -> None: ...

class PointCloud:
    def __init__(self, arg0) -> None: ...
    def insertNextPoint(self, arg0) -> None: ...

class PointData:
    def __init__(self) -> None: ...
    @overload
    def getScalarData(self, arg0: int) -> list[float]: ...
    @overload
    def getScalarData(self, arg0: str, arg1: bool) -> list[float]: ...
    def getScalarDataLabel(self, arg0: int) -> str: ...
    def getScalarDataSize(self) -> int: ...
    def getVectorData(self, *args, **kwargs): ...
    def getVectorDataLabel(self, arg0: int) -> str: ...
    def getVectorDataSize(self) -> int: ...
    def insertNextScalarData(self, scalars: list[float], label: str = ...) -> None: ...
    def insertNextVectorData(self, vectors, label: str = ...) -> None: ...

class Prune:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Domain) -> None: ...
    def apply(self) -> None: ...
    def setLevelSet(self, arg0: Domain) -> None: ...

class Reader:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Domain) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: str) -> None: ...
    def apply(self) -> None: ...
    def setFileName(self, arg0: str) -> None: ...
    def setLevelSet(self, arg0: Domain) -> None: ...

class Reduce:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Domain) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: int) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: int, arg2: bool) -> None: ...
    def apply(self) -> None: ...
    def setLevelSet(self, arg0: Domain) -> None: ...
    def setNoNewSegment(self, arg0: bool) -> None: ...
    def setWidth(self, arg0: int) -> None: ...

class RemoveStrayPoints:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Domain) -> None: ...
    def apply(self) -> None: ...
    def setLevelSet(self, arg0: Domain) -> None: ...
    def setVoidTopSurface(self, arg0: VoidTopSurfaceEnum) -> None: ...

class Slice:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: int, arg2: float) -> None: ...
    def apply(self) -> None: ...
    def getSliceLevelSet(self, *args, **kwargs): ...
    def setReflectX(self, arg0: bool) -> None: ...
    def setSliceDimension(self, arg0: int) -> None: ...
    def setSliceLevelSet(self, *args, **kwargs): ...
    def setSlicePosition(self, arg0: float) -> None: ...
    def setSourceLevelSet(self, arg0: Domain) -> None: ...
    def setWritePath(self, arg0: str, arg1: bool) -> None: ...

class Sphere:
    def __init__(self, origin: list[float], radius: float) -> None: ...

class SphereDistribution(GeometricAdvectDistribution):
    def __init__(self, arg0: float, arg1: float) -> None: ...
    def getBounds(self, *args, **kwargs): ...
    def getSignedDistance(self, arg0, arg1, arg2: int) -> float: ...
    def isInside(self, arg0, arg1, arg2: float) -> bool: ...

class ToDiskMesh:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: Mesh) -> None: ...
    def apply(self) -> None: ...
    def setLevelSet(self, arg0: Domain) -> None: ...
    def setMesh(self, arg0: Mesh) -> None: ...

class ToMesh:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: Mesh) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: Mesh, arg2: bool) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: Mesh, arg2: bool, arg3: bool) -> None: ...
    def apply(self) -> None: ...
    def setLevelSet(self, arg0: Domain) -> None: ...
    def setMesh(self, arg0: Mesh) -> None: ...
    def setOnlyActive(self, arg0: bool) -> None: ...
    def setOnlyDefined(self, arg0: bool) -> None: ...

class ToSurfaceMesh:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: Mesh) -> None: ...
    def apply(self) -> None: ...
    def setLevelSet(self, arg0: Domain) -> None: ...
    def setMesh(self, arg0: Mesh) -> None: ...

class ToVoxelMesh:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Mesh) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: Mesh) -> None: ...
    @overload
    def __init__(self, arg0: list[Domain], arg1: Mesh) -> None: ...
    def apply(self) -> None: ...
    def insertNextLevelSet(self, arg0: Domain) -> None: ...
    def setMesh(self, arg0: Mesh) -> None: ...

class TransformEnum:
    __members__: ClassVar[dict] = ...  # read-only
    ROTATION: ClassVar[TransformEnum] = ...
    SCALE: ClassVar[TransformEnum] = ...
    TRANSLATION: ClassVar[TransformEnum] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class TransformMesh:
    def __init__(
        self,
        mesh: Mesh,
        transform: TransformEnum = ...,
        transformVector=...,
        angle: float = ...,
    ) -> None: ...
    def apply(self) -> None: ...

class VTKReader:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Mesh) -> None: ...
    @overload
    def __init__(self, arg0: Mesh, arg1: str) -> None: ...
    @overload
    def __init__(self, arg0: Mesh, arg1: FileFormatEnum, arg2: str) -> None: ...
    def apply(self) -> None: ...
    def getMetaData(self) -> dict[str, list[float]]: ...
    def setFileFormat(self, arg0: FileFormatEnum) -> None: ...
    def setFileName(self, arg0: str) -> None: ...
    def setMesh(self, arg0: Mesh) -> None: ...

class VTKWriter:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Mesh) -> None: ...
    @overload
    def __init__(self, arg0: Mesh, arg1: str) -> None: ...
    @overload
    def __init__(self, arg0: Mesh, arg1: FileFormatEnum, arg2: str) -> None: ...
    @overload
    def addMetaData(self, arg0: str, arg1: float) -> None: ...
    @overload
    def addMetaData(self, arg0: str, arg1: list[float]) -> None: ...
    @overload
    def addMetaData(self, arg0: dict[str, list[float]]) -> None: ...
    def apply(self) -> None: ...
    def setFileFormat(self, arg0: FileFormatEnum) -> None: ...
    def setFileName(self, arg0: str) -> None: ...
    def setMesh(self, arg0: Mesh) -> None: ...
    def setMetaData(self, arg0: dict[str, list[float]]) -> None: ...

class VelocityField:
    def __init__(self) -> None: ...
    def getDissipationAlpha(self, arg0: int, arg1: int, arg2) -> float: ...
    def getScalarVelocity(self, arg0, arg1: int, arg2, arg3: int) -> float: ...
    def getVectorVelocity(self, *args, **kwargs): ...

class VoidTopSurfaceEnum:
    __members__: ClassVar[dict] = ...  # read-only
    LARGEST: ClassVar[VoidTopSurfaceEnum] = ...
    LEX_HIGHEST: ClassVar[VoidTopSurfaceEnum] = ...
    LEX_LOWEST: ClassVar[VoidTopSurfaceEnum] = ...
    SMALLEST: ClassVar[VoidTopSurfaceEnum] = ...
    __entries: ClassVar[dict] = ...
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    def __ne__(self, other: object) -> bool: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

class WriteVisualizationMesh:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Domain) -> None: ...
    @overload
    def addMetaData(self, arg0: str, arg1: float) -> None: ...
    @overload
    def addMetaData(self, arg0: str, arg1: list[float]) -> None: ...
    @overload
    def addMetaData(self, arg0: dict[str, list[float]]) -> None: ...
    def apply(self) -> None: ...
    def insertNextLevelSet(self, arg0: Domain) -> None: ...
    def setExtractHullMesh(self, arg0: bool) -> None: ...
    def setExtractVolumeMesh(self, arg0: bool) -> None: ...
    def setFileName(self, arg0: str) -> None: ...
    def setMetaData(self, arg0: dict[str, list[float]]) -> None: ...

class Writer:
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, arg0: Domain) -> None: ...
    @overload
    def __init__(self, arg0: Domain, arg1: str) -> None: ...
    def apply(self) -> None: ...
    def setFileName(self, arg0: str) -> None: ...
    def setLevelSet(self, arg0: Domain) -> None: ...

class hrleGrid:
    def __init__(self, *args, **kwargs) -> None: ...

def setNumThreads(arg0: int) -> None: ...
