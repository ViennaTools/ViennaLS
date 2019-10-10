/*
  Build file for shared library instantiating precompiled
  template specialisations
*/
#include <lsExpand_template.hpp>

// Specialise template to declare all members which only depend on class
// template type
template class lsExpand<double, 2>;
template class lsExpand<double, 3>;

// Specialise members which depend on a different (nested) templated type
