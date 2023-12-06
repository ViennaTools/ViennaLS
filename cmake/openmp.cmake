if(APPLE)
  find_package(OpenMP COMPONENTS CXX)

  if(NOT OpenMP_CXX_FOUND)
    execute_process(
      COMMAND brew --prefix libomp
      OUTPUT_VARIABLE HOMEBREW_LIBOMP_PREFIX
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${HOMEBREW_LIBOMP_PREFIX}/include")
    set(OpenMP_omp_LIBRARY ${HOMEBREW_LIBOMP_PREFIX}/lib/libomp.dylib)
  endif()
endif()
