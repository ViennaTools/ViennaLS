/*
  This file declares the precompiled specialisations for the lsDomain class
*/

#include <lsToVoxelMesh_template.hpp>

// SPECIALISATION CONVENIENCE TYPEDEFS
typedef lsToVoxelMesh<double, 2> lsToVoxelMesh_double_2;
typedef lsToVoxelMesh<double, 3> lsToVoxelMesh_double_3;

#ifndef VIENNALS_USE_HEADER_ONLY

// Specialise template to declare all members which only depend on class
// template type
extern template class lsToVoxelMesh<double, 2>;
extern template class lsToVoxelMesh<double, 3>;

// Specialise members which depend on a different (nested) templated type

#endif // VIENNA_LS_REBUILD_TEMPLATES
