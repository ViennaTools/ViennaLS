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
from viennals2d.viennals2d import Advect
from viennals2d.viennals2d import BooleanOperation
from viennals2d.viennals2d import BooleanOperationEnum
from viennals2d.viennals2d import BoundaryConditionEnum
from viennals2d.viennals2d import Box
from viennals2d.viennals2d import BoxDistribution
from viennals2d.viennals2d import CalculateCurvatures
from viennals2d.viennals2d import CalculateNormalVectors
from viennals2d.viennals2d import CalculateVisibilities
from viennals2d.viennals2d import Check
from viennals2d.viennals2d import CompareArea
from viennals2d.viennals2d import CompareNarrowBand
from viennals2d.viennals2d import CompareSparseField
from viennals2d.viennals2d import ConvexHull
from viennals2d.viennals2d import CurvatureEnum
from viennals2d.viennals2d import Cylinder
from viennals2d.viennals2d import DetectFeatures
from viennals2d.viennals2d import Domain
from viennals2d.viennals2d import Expand
from viennals2d.viennals2d import FeatureDetectionEnum
from viennals2d.viennals2d import FileFormatEnum
from viennals2d.viennals2d import FromMesh
from viennals2d.viennals2d import FromSurfaceMesh
from viennals2d.viennals2d import FromVolumeMesh
from viennals2d.viennals2d import GeometricAdvect
from viennals2d.viennals2d import GeometricAdvectDistribution
from viennals2d.viennals2d import IntegrationSchemeEnum
from viennals2d.viennals2d import LogLevel
from viennals2d.viennals2d import Logger
from viennals2d.viennals2d import MakeGeometry
from viennals2d.viennals2d import MarkVoidPoints
from viennals2d.viennals2d import MaterialMap
from viennals2d.viennals2d import Mesh
from viennals2d.viennals2d import Plane
from viennals2d.viennals2d import PointCloud
from viennals2d.viennals2d import PointData
from viennals2d.viennals2d import Prune
from viennals2d.viennals2d import Reader
from viennals2d.viennals2d import Reduce
from viennals2d.viennals2d import RemoveStrayPoints
from viennals2d.viennals2d import Slice
from viennals2d.viennals2d import Sphere
from viennals2d.viennals2d import SphereDistribution
from viennals2d.viennals2d import StencilLocalLaxFriedrichsScalar
from viennals2d.viennals2d import ToDiskMesh
from viennals2d.viennals2d import ToMesh
from viennals2d.viennals2d import ToSurfaceMesh
from viennals2d.viennals2d import ToVoxelMesh
from viennals2d.viennals2d import TransformEnum
from viennals2d.viennals2d import TransformMesh
from viennals2d.viennals2d import VTKReader
from viennals2d.viennals2d import VTKWriter
from viennals2d.viennals2d import VelocityField
from viennals2d.viennals2d import VoidTopSurfaceEnum
from viennals2d.viennals2d import WriteVisualizationMesh
from viennals2d.viennals2d import Writer
from viennals2d.viennals2d import hrleGrid
from viennals2d.viennals2d import setNumThreads
from . import viennals2d

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
    "CompareArea",
    "CompareNarrowBand",
    "CompareSparseField",
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
    "viennals2d",
]

def _windows_dll_path(): ...

version: str = "4.5.0"
