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
from viennals3d.viennals3d import Advect
from viennals3d.viennals3d import BooleanOperation
from viennals3d.viennals3d import BooleanOperationEnum
from viennals3d.viennals3d import BoundaryConditionEnum
from viennals3d.viennals3d import Box
from viennals3d.viennals3d import BoxDistribution
from viennals3d.viennals3d import CalculateCurvatures
from viennals3d.viennals3d import CalculateNormalVectors
from viennals3d.viennals3d import CalculateVisibilities
from viennals3d.viennals3d import Check
from viennals3d.viennals3d import ConvexHull
from viennals3d.viennals3d import CurvatureEnum
from viennals3d.viennals3d import Cylinder
from viennals3d.viennals3d import DetectFeatures
from viennals3d.viennals3d import Domain
from viennals3d.viennals3d import Expand
from viennals3d.viennals3d import FeatureDetectionEnum
from viennals3d.viennals3d import FileFormatEnum
from viennals3d.viennals3d import FromMesh
from viennals3d.viennals3d import FromSurfaceMesh
from viennals3d.viennals3d import FromVolumeMesh
from viennals3d.viennals3d import GeometricAdvect
from viennals3d.viennals3d import GeometricAdvectDistribution
from viennals3d.viennals3d import IntegrationSchemeEnum
from viennals3d.viennals3d import LogLevel
from viennals3d.viennals3d import Logger
from viennals3d.viennals3d import MakeGeometry
from viennals3d.viennals3d import MarkVoidPoints
from viennals3d.viennals3d import MaterialMap
from viennals3d.viennals3d import Mesh
from viennals3d.viennals3d import Plane
from viennals3d.viennals3d import PointCloud
from viennals3d.viennals3d import PointData
from viennals3d.viennals3d import Prune
from viennals3d.viennals3d import Reader
from viennals3d.viennals3d import Reduce
from viennals3d.viennals3d import RemoveStrayPoints
from viennals3d.viennals3d import Slice
from viennals3d.viennals3d import Sphere
from viennals3d.viennals3d import SphereDistribution
from viennals3d.viennals3d import StencilLocalLaxFriedrichsScalar
from viennals3d.viennals3d import ToDiskMesh
from viennals3d.viennals3d import ToMesh
from viennals3d.viennals3d import ToSurfaceMesh
from viennals3d.viennals3d import ToVoxelMesh
from viennals3d.viennals3d import TransformEnum
from viennals3d.viennals3d import TransformMesh
from viennals3d.viennals3d import VTKReader
from viennals3d.viennals3d import VTKWriter
from viennals3d.viennals3d import VelocityField
from viennals3d.viennals3d import VoidTopSurfaceEnum
from viennals3d.viennals3d import WriteVisualizationMesh
from viennals3d.viennals3d import Writer
from viennals3d.viennals3d import hrleGrid
from viennals3d.viennals3d import setNumThreads
from . import viennals3d

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
    "sys",
    "version",
    "viennals3d",
]

def _windows_dll_path(): ...

version: str = "4.5.0"
