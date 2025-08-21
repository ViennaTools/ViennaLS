"""

ViennaLS
========

ViennaLS is a level set library developed for high performance
topography simulations. The main design goals are simplicity and efficiency,
tailored towards scientific simulations. ViennaLS can also be used for
visualisation applications, although this is not the main design target.
"""

from __future__ import annotations
import sys as sys
from viennals.viennals import BooleanOperationEnum
from viennals.viennals import BoundaryConditionEnum
from viennals.viennals import CurvatureEnum
from viennals.viennals import Extrude
from viennals.viennals import FeatureDetectionEnum
from viennals.viennals import FileFormatEnum
from viennals.viennals import IntegrationSchemeEnum
from viennals.viennals import LogLevel
from viennals.viennals import Logger
from viennals.viennals import MaterialMap
from viennals.viennals import Mesh
from viennals.viennals import PointData
from viennals.viennals import Slice
from viennals.viennals import TransformEnum
from viennals.viennals import TransformMesh
from viennals.viennals import VTKReader
from viennals.viennals import VTKWriter
from viennals.viennals import VelocityField
from viennals.viennals import VoidTopSurfaceEnum
from viennals.viennals import d2
from viennals.viennals import d3
from viennals.viennals import setNumThreads
from . import viennals

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
    "sys",
    "version",
    "viennals",
]

def _windows_dll_path(): ...

version: str = "5.0.0"
