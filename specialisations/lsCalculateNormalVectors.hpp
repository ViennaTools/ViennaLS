/*
  This file declares the precompiled specialisations for the lsDomain class
*/

#include <lsCalculateNormalVectors_template.hpp>

// SPECIALISATION CONVENIENCE TYPEDEFS
typedef lsCalculateNormalVectors<double, 2> lsCalculateNormalVectors_double_2;
typedef lsCalculateNormalVectors<double, 3> lsCalculateNormalVectors_double_3;

#ifndef VIENNALS_USE_HEADER_ONLY

// Specialise template to declare all members which only depend on class
// template type
extern template class lsCalculateNormalVectors<double, 2>;
extern template class lsCalculateNormalVectors<double, 3>;

// Specialise members which depend on a different (nested) templated type

#endif // VIENNA_LS_REBUILD_TEMPLATES
