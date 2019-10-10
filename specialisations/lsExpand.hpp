/*
  This file declares the precompiled specialisations for the lsDomain class
*/

#include <lsExpand_template.hpp>

// SPECIALISATION CONVENIENCE TYPEDEFS
typedef lsExpand<double, 2> lsExpand_double_2;
typedef lsExpand<double, 3> lsExpand_double_3;

#ifndef VIENNALS_USE_HEADER_ONLY

// Specialise template to declare all members which only depend on class
// template type
extern template class lsExpand<double, 2>;
extern template class lsExpand<double, 3>;

// Specialise members which depend on a different (nested) templated type

#endif // VIENNA_LS_REBUILD_TEMPLATES
