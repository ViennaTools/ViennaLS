/*
  Build file for shared library instantiating precompiled
  template specialisations
*/
#include <lsDomain_template.hpp>

// Specialise template to declare all members which only depend on class
// template type
template class lsDomain<double, 2>;
template class lsDomain<double, 3>;

// Specialise members which depend on a different (nested) templated type
