/*
  This file declares the precompiled specialisations for the lsDomain class
*/

#include <lsDomain_template.hpp>

// SPECIALISATION CONVENIENCE TYPEDEFS
typedef lsDomain<double, 2> lsDomain_double_2;
typedef lsDomain<double, 3> lsDomain_double_3;

#ifndef VIENNALS_USE_HEADER_ONLY

// Specialise template to declare all members which only depend on class
// template type
extern template class lsDomain<double, 2>;
extern template class lsDomain<double, 3>;

// Specialise members which depend on a different (nested) templated type

#endif // VIENNA_LS_REBUILD_TEMPLATES
