/*
  This file declares the precompiled specialisations for the lsDomain class
*/

#include <lsBooleanOperation_template.hpp>

// SPECIALISATION CONVENIENCE TYPEDEFS
typedef lsBooleanOperation<double, 2> lsBooleanOperation_double_2;
typedef lsBooleanOperation<double, 3> lsBooleanOperation_double_3;

#ifndef VIENNALS_USE_HEADER_ONLY

// Specialise template to declare all members which only depend on class
// template type
extern template class lsBooleanOperation<double, 2>;
extern template class lsBooleanOperation<double, 3>;

// Specialise members which depend on a different (nested) templated type

#endif // VIENNA_LS_REBUILD_TEMPLATES
