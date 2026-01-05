"""

ViennaLS
========

ViennaLS is a level set library developed for high performance
topography simulations. The main design goals are simplicity and efficiency,
tailored towards scientific simulations. ViennaLS can also be used for
visualisation applications, although this is not the main design target.
"""
from __future__ import annotations
import sys as _sys
from viennals._core import BooleanOperationEnum
from viennals._core import BoundaryConditionEnum
from viennals._core import CurvatureEnum
from viennals._core import Extrude
from viennals._core import FeatureDetectionEnum
from viennals._core import FileFormatEnum
from viennals._core import LogLevel
from viennals._core import Logger
from viennals._core import MaterialMap
from viennals._core import Mesh
from viennals._core import PointData
from viennals._core import Slice
from viennals._core import SpatialSchemeEnum as IntegrationSchemeEnum
from viennals._core import SpatialSchemeEnum
from viennals._core import TemporalSchemeEnum
from viennals._core import TransformEnum
from viennals._core import TransformMesh
from viennals._core import VTKReader
from viennals._core import VTKRenderWindow
from viennals._core import VTKWriter
from viennals._core import VelocityField
from viennals._core import VoidTopSurfaceEnum
from viennals._core import setNumThreads
from viennals.d2 import Advect
from viennals.d2 import BooleanOperation
from viennals.d2 import Box
from viennals.d2 import BoxDistribution
from viennals.d2 import CalculateCurvatures
from viennals.d2 import CalculateNormalVectors
from viennals.d2 import CalculateVisibilities
from viennals.d2 import Check
from viennals.d2 import CompareArea
from viennals.d2 import CompareChamfer
from viennals.d2 import CompareCriticalDimensions
from viennals.d2 import CompareNarrowBand
from viennals.d2 import CompareSparseField
from viennals.d2 import ConvexHull
from viennals.d2 import CustomSphereDistribution
from viennals.d2 import Cylinder
from viennals.d2 import DetectFeatures
from viennals.d2 import Domain
from viennals.d2 import Expand
from viennals.d2 import FinalizeStencilLocalLaxFriedrichs
from viennals.d2 import FromMesh
from viennals.d2 import FromSurfaceMesh
from viennals.d2 import FromVolumeMesh
from viennals.d2 import GeometricAdvect
from viennals.d2 import GeometricAdvectDistribution
from viennals.d2 import MakeGeometry
from viennals.d2 import MarkVoidPoints
from viennals.d2 import Plane
from viennals.d2 import PointCloud
from viennals.d2 import PrepareStencilLocalLaxFriedrichs
from viennals.d2 import Prune
from viennals.d2 import Reader
from viennals.d2 import Reduce
from viennals.d2 import RemoveStrayPoints
from viennals.d2 import Sphere
from viennals.d2 import SphereDistribution
from viennals.d2 import StencilLocalLaxFriedrichsScalar
from viennals.d2 import ToDiskMesh
from viennals.d2 import ToMesh
from viennals.d2 import ToMultiSurfaceMesh
from viennals.d2 import ToSurfaceMesh
from viennals.d2 import ToVoxelMesh
from viennals.d2 import WriteVisualizationMesh
from viennals.d2 import Writer
from viennals.d2 import hrleGrid
from . import _core
from . import d2
from . import d3
__all__: list[str] = ['Advect', 'BooleanOperation', 'BooleanOperationEnum', 'BoundaryConditionEnum', 'Box', 'BoxDistribution', 'CalculateCurvatures', 'CalculateNormalVectors', 'CalculateVisibilities', 'Check', 'CompareArea', 'CompareChamfer', 'CompareCriticalDimensions', 'CompareNarrowBand', 'CompareSparseField', 'ConvexHull', 'CurvatureEnum', 'CustomSphereDistribution', 'Cylinder', 'DetectFeatures', 'Domain', 'Expand', 'Extrude', 'FeatureDetectionEnum', 'FileFormatEnum', 'FinalizeStencilLocalLaxFriedrichs', 'FromMesh', 'FromSurfaceMesh', 'FromVolumeMesh', 'GeometricAdvect', 'GeometricAdvectDistribution', 'IntegrationSchemeEnum', 'LogLevel', 'Logger', 'MakeGeometry', 'MarkVoidPoints', 'MaterialMap', 'Mesh', 'PROXY_DIM', 'Plane', 'PointCloud', 'PointData', 'PrepareStencilLocalLaxFriedrichs', 'Prune', 'Reader', 'Reduce', 'RemoveStrayPoints', 'Slice', 'SpatialSchemeEnum', 'Sphere', 'SphereDistribution', 'StencilLocalLaxFriedrichsScalar', 'TemporalSchemeEnum', 'ToDiskMesh', 'ToMesh', 'ToMultiSurfaceMesh', 'ToSurfaceMesh', 'ToVoxelMesh', 'TransformEnum', 'TransformMesh', 'VTKReader', 'VTKRenderWindow', 'VTKWriter', 'VelocityField', 'VoidTopSurfaceEnum', 'WriteVisualizationMesh', 'Writer', 'd2', 'd3', 'getDimension', 'hrleGrid', 'setDimension', 'setNumThreads', 'version']
def __dir__():
    ...
def __getattr__(name):
    ...
def _windows_dll_path():
    ...
def getDimension() -> int:
    """
    Get the current dimension of the simulation.
    
        Returns
        -------
        int
            The currently set dimension (2 or 3).
        
    """
def setDimension(d: int):
    """
    Set the dimension of the simulation (2 or 3).
    
        Parameters
        ----------
        d: int
            Dimension of the simulation (2 or 3).
        
    """
PROXY_DIM: int = 2
__version__: str = '5.4.0'
version: str = '5.4.0'
_C = _core
