/*
  This header creates all the template specializations
  in order to build the shared library for later linking.
*/

// add all headers which require template specialisation
#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsCalculateCurvatures.hpp>
#include <lsCalculateNormalVectors.hpp>
#include <lsCheck.hpp>
#include <lsCompareArea.hpp>
#include <lsCompareNarrowBand.hpp>
#include <lsCompareSparseField.hpp>
#include <lsDetectFeatures.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromSurfaceMesh.hpp>
#include <lsFromVolumeMesh.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsReader.hpp>
#include <lsReduce.hpp>
#include <lsSlice.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsWriteVisualizationMesh.hpp>
#include <lsWriter.hpp>

namespace viennals {

// now call the specialize macro to precompile them
PRECOMPILE_SPECIALIZE_PRECISION(PointData)
PRECOMPILE_SPECIALIZE_PRECISION(Mesh)
PRECOMPILE_SPECIALIZE(Advect)
PRECOMPILE_SPECIALIZE(BooleanOperation)
PRECOMPILE_SPECIALIZE(CalculateCurvatures)
PRECOMPILE_SPECIALIZE(CalculateNormalVectors)
PRECOMPILE_SPECIALIZE(Check)
PRECOMPILE_SPECIALIZE(ConvexHull)
PRECOMPILE_SPECIALIZE(CompareArea)
PRECOMPILE_SPECIALIZE(CompareNarrowBand)
PRECOMPILE_SPECIALIZE(CompareSparseField)
PRECOMPILE_SPECIALIZE(Domain)
PRECOMPILE_SPECIALIZE(Expand)
PRECOMPILE_SPECIALIZE(GeometricAdvect)
PRECOMPILE_SPECIALIZE(DetectFeatures)
PRECOMPILE_SPECIALIZE(FromMesh)
PRECOMPILE_SPECIALIZE(FromSurfaceMesh)
PRECOMPILE_SPECIALIZE(FromVolumeMesh)
PRECOMPILE_SPECIALIZE(Sphere)
PRECOMPILE_SPECIALIZE(Plane)
PRECOMPILE_SPECIALIZE(Box)
PRECOMPILE_SPECIALIZE(PointCloud)
PRECOMPILE_SPECIALIZE(MakeGeometry)
PRECOMPILE_SPECIALIZE(Prune)
PRECOMPILE_SPECIALIZE(Reader)
PRECOMPILE_SPECIALIZE(Reduce)
PRECOMPILE_SPECIALIZE_PRECISION(Slice)
PRECOMPILE_SPECIALIZE(ToDiskMesh)
PRECOMPILE_SPECIALIZE(ToMesh)
PRECOMPILE_SPECIALIZE(ToSurfaceMesh)
PRECOMPILE_SPECIALIZE(ToVoxelMesh)
PRECOMPILE_SPECIALIZE(Writer)
#ifdef VIENNALS_USE_VTK
PRECOMPILE_SPECIALIZE(WriteVisualizationMesh)
#endif

} // namespace viennals