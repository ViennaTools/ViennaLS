/*
  Build file for shared library instantiating precompiled
  template specialisations
*/
#include <lsToExplicitMesh_template.hpp>

// Specialise template to declare all members which only depend on class
// template type
template class lsToExplicitMesh<double, 2>;
template class lsToExplicitMesh<double, 3>;

// Specialise members which depend on a different (nested) templated type
