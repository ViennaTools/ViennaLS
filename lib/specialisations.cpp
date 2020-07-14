/*
  This header creates all the template specializations
  in order to build the shared library for later linking.
*/

// add all headers which require template specialisation
#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsCalculateNormalVectors.hpp>
#include <lsCheck.hpp>
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
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsToVoxelMesh.hpp>
#include <lsWriter.hpp>

// now call the specialize macro to precompile them
PRECOMPILE_SPECIALIZE(lsAdvect)
PRECOMPILE_SPECIALIZE(lsBooleanOperation)
PRECOMPILE_SPECIALIZE(lsCalculateNormalVectors)
PRECOMPILE_SPECIALIZE(lsCheck)
PRECOMPILE_SPECIALIZE(lsConvexHull)
PRECOMPILE_SPECIALIZE(lsDomain)
PRECOMPILE_SPECIALIZE(lsExpand)
PRECOMPILE_SPECIALIZE(lsGeometricAdvect)
PRECOMPILE_SPECIALIZE(lsFromMesh)
PRECOMPILE_SPECIALIZE(lsFromSurfaceMesh)
PRECOMPILE_SPECIALIZE(lsFromVolumeMesh)
PRECOMPILE_SPECIALIZE(lsSphere)
PRECOMPILE_SPECIALIZE(lsPlane)
PRECOMPILE_SPECIALIZE(lsBox)
PRECOMPILE_SPECIALIZE(lsPointCloud)
PRECOMPILE_SPECIALIZE(lsMakeGeometry)
PRECOMPILE_SPECIALIZE(lsPrune)
PRECOMPILE_SPECIALIZE(lsReader)
PRECOMPILE_SPECIALIZE(lsReduce)
PRECOMPILE_SPECIALIZE(lsToDiskMesh)
PRECOMPILE_SPECIALIZE(lsToMesh)
PRECOMPILE_SPECIALIZE(lsToSurfaceMesh)
PRECOMPILE_SPECIALIZE(lsToVoxelMesh)
PRECOMPILE_SPECIALIZE(lsWriter)
