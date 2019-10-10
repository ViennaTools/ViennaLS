/*
  This file declares the precompiled specialisations for the lsDomain class
*/

#include <lsToExplicitMesh_template.hpp>

// SPECIALISATION CONVENIENCE TYPEDEFS
typedef lsToExplicitMesh<double, 2> lsToExplicitMesh_double_2;
typedef lsToExplicitMesh<double, 3> lsToExplicitMesh_double_3;

#ifndef VIENNALS_USE_HEADER_ONLY

// Specialise template to declare all members which only depend on class
// template type
extern template class lsToExplicitMesh<double, 2>;
extern template class lsToExplicitMesh<double, 3>;

// Specialise members which depend on a different (nested) templated type

#endif // VIENNA_LS_REBUILD_TEMPLATES
