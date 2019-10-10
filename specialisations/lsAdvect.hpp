/*
  This file declares the precompiled specialisations for the lsDomain class
*/

#include <lsAdvect_template.hpp>

// SPECIALISATION CONVENIENCE TYPEDEFS
typedef lsAdvect<double, 2> lsAdvect_double_2;
typedef lsAdvect<double, 3> lsAdvect_double_3;

#ifndef VIENNALS_USE_HEADER_ONLY

// Specialise template to declare all members which only depend on class
// template type
extern template class lsAdvect<double, 2>;
extern template class lsAdvect<double, 3>;

// Specialise members which depend on a different (nested) templated type

#endif // VIENNA_LS_REBUILD_TEMPLATES
