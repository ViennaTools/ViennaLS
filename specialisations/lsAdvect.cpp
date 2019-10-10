/*
  Build file for shared library instantiating precompiled
  template specialisations
*/
#include <lsAdvect_template.hpp>

// Specialise template to declare all members which only depend on class
// template type
template class lsAdvect<double, 2>;
template class lsAdvect<double, 3>;

// Specialise members which depend on a different (nested) templated type
