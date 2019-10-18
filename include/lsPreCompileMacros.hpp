#ifndef LS_PRE_COMPILE_MACROS_HPP
#define LS_PRE_COMPILE_MACROS_HPP

#ifndef VIENNALS_USE_HEADER_ONLY

#define PRECOMPILE_PRECISION_DIMENSION(className) \
  typedef className<double, 2> className##_double_2;  \
  typedef className<double, 3> className##_double_3;  \
  extern template class className<double, 2>;  \
  extern template class className<double, 3>;

#else

// do nothing if we use header only
#define PRECOMPILE_PRECISION_DIMENSION(className) \
  typedef className<double, 2> className##_double_2;  \
  typedef className<double, 3> className##_double_3;

#endif

#define PRECOMPILE_SPECIALIZE(className)  \
  template class className<double, 2>;  \
  template class className<double, 3>;

#endif // LS_PRE_COMPILE_MACROS_HPP
