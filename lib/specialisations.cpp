/*
  This header creates all the template specializations
  in order to build the shared library for later linking.
*/

// add all headers which require template specialisation
#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsCalculateNormalVectors.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsFromExplicitMesh.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsReduce.hpp>
#include <lsToExplicitMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToVoxelMesh.hpp>

// now call the specialize macro to precompile them
PRECOMPILE_SPECIALIZE(lsAdvect)
PRECOMPILE_SPECIALIZE(lsBooleanOperation)
PRECOMPILE_SPECIALIZE(lsCalculateNormalVectors)
PRECOMPILE_SPECIALIZE(lsDomain)
PRECOMPILE_SPECIALIZE(lsExpand)
PRECOMPILE_SPECIALIZE(lsFromExplicitMesh)
PRECOMPILE_SPECIALIZE(lsMakeGeometry)
PRECOMPILE_SPECIALIZE(lsPrune)
PRECOMPILE_SPECIALIZE(lsReduce)
PRECOMPILE_SPECIALIZE(lsToExplicitMesh)
PRECOMPILE_SPECIALIZE(lsToMesh)
PRECOMPILE_SPECIALIZE(lsToVoxelMesh)
