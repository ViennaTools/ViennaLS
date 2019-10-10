/*
  This file declares the precompiled specialisations for the lsDomain class
*/

#include <lsReduce_template.hpp>

// SPECIALISATION CONVENIENCE TYPEDEFS
typedef lsReduce<double, 2> lsReduce_double_2;
typedef lsReduce<double, 3> lsReduce_double_3;

#ifndef VIENNALS_USE_HEADER_ONLY

// Specialise template to declare all members which only depend on class
// template type
extern template class lsReduce<double, 2>;
extern template class lsReduce<double, 3>;

// Specialise members which depend on a different (nested) templated type

#endif // VIENNA_LS_REBUILD_TEMPLATES
