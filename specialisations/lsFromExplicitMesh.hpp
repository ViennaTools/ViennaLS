/*
  This file declares the precompiled specialisations for the lsDomain class
*/

#include <lsFromExplicitMesh_template.hpp>

// SPECIALISATION CONVENIENCE TYPEDEFS
typedef lsFromExplicitMesh<double, 2> lsFromExplicitMesh_double_2;
typedef lsFromExplicitMesh<double, 3> lsFromExplicitMesh_double_3;

#ifndef VIENNALS_USE_HEADER_ONLY

// Specialise template to declare all members which only depend on class
// template type
extern template class lsFromExplicitMesh<double, 2>;
extern template class lsFromExplicitMesh<double, 3>;

// Specialise members which depend on a different (nested) templated type

#endif // VIENNA_LS_REBUILD_TEMPLATES
