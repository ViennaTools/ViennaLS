/*
  Build file for shared library instantiating precompiled
  template specialisations
*/
#include <lsMakeGeometry_template.hpp>

// Specialise template to declare all members which only depend on class
// template type
template class lsMakeGeometry<double, 2>;
template class lsMakeGeometry<double, 3>;

// Specialise members which depend on a different (nested) templated type
template void lsMakeGeometry<double, 2>::makeSphere<double *>(double *, double,
                                                              int);
template void lsMakeGeometry<double, 3>::makeSphere<double *>(double *, double,
                                                              int);
template void
    lsMakeGeometry<double, 2>::makeSphere<hrleVectorType<hrleIndexType, 2>>(
        hrleVectorType<hrleIndexType, 2>, double, int);
template void
    lsMakeGeometry<double, 3>::makeSphere<hrleVectorType<hrleIndexType, 3>>(
        hrleVectorType<hrleIndexType, 3>, double, int);
