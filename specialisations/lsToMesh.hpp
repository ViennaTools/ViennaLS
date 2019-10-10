/*
  This file declares the precompiled specialisations for the lsDomain class
*/

#include <lsToMesh_template.hpp>

// SPECIALISATION CONVENIENCE TYPEDEFS
typedef lsToMesh<double, 2> lsToMesh_double_2;
typedef lsToMesh<double, 3> lsToMesh_double_3;

#ifndef VIENNALS_USE_HEADER_ONLY

// Specialise template to declare all members which only depend on class
// template type
extern template class lsToMesh<double, 2>;
extern template class lsToMesh<double, 3>;

// Specialise members which depend on a different (nested) templated type

#endif // VIENNA_LS_REBUILD_TEMPLATES
