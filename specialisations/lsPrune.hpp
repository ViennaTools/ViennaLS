/*
  This file declares the precompiled specialisations for the lsDomain class
*/

#include <lsPrune_template.hpp>

// SPECIALISATION CONVENIENCE TYPEDEFS
typedef lsPrune<double, 2> lsPrune_double_2;
typedef lsPrune<double, 3> lsPrune_double_3;

#ifndef VIENNALS_USE_HEADER_ONLY

// Specialise template to declare all members which only depend on class
// template type
extern template class lsPrune<double, 2>;
extern template class lsPrune<double, 3>;

// Specialise members which depend on a different (nested) templated type

#endif // VIENNA_LS_REBUILD_TEMPLATES
