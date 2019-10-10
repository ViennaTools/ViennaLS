/*
  This file declares the precompiled specialisations for the lsDomain class
*/

#include <lsMakeGeometry_template.hpp>

// SPECIALISATION CONVENIENCE TYPEDEFS
typedef lsMakeGeometry<double, 2> lsMakeGeometry_double_2;
typedef lsMakeGeometry<double, 3> lsMakeGeometry_double_3;

#ifndef VIENNALS_USE_HEADER_ONLY

// Specialise template to declare all members which only depend on class
// template type
extern template class lsMakeGeometry<double, 2>;
extern template class lsMakeGeometry<double, 3>;

// Specialise members which depend on a different (nested) templated type
extern template void
lsMakeGeometry<double, 2>::makeSphere<double *>(double *, double, int);
extern template void
lsMakeGeometry<double, 3>::makeSphere<double *>(double *, double, int);
extern template void
    lsMakeGeometry<double, 2>::makeSphere<hrleVectorType<hrleIndexType, 2>>(
        hrleVectorType<hrleIndexType, 2>, double, int);
extern template void
    lsMakeGeometry<double, 3>::makeSphere<hrleVectorType<hrleIndexType, 3>>(
        hrleVectorType<hrleIndexType, 3>, double, int);

#endif // VIENNA_LS_REBUILD_TEMPLATES
