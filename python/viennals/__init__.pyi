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
from viennals._core import IntegrationSchemeEnum
from viennals._core import LogLevel
from viennals._core import Logger
from viennals._core import MaterialMap
from viennals._core import Mesh
from viennals._core import PointData
from viennals._core import Slice
from viennals._core import TransformEnum
from viennals._core import TransformMesh
from viennals._core import VTKReader
from viennals._core import VTKWriter
from viennals._core import VelocityField
from viennals._core import VoidTopSurfaceEnum
from viennals._core import setNumThreads
from . import _core
from . import d2
from . import d3

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

def __dir__(): ...
def __getattr__(name): ...
def _windows_dll_path(): ...

__version__: str = ""
version: str = ""
_C = _core
