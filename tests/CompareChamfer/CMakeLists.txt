project(CompareChamfer LANGUAGES CXX)

add_executable(${PROJECT_NAME} "${PROJECT_NAME}.cpp")
target_link_libraries(${PROJECT_NAME} PRIVATE ViennaLS)

# Add OpenMP support for MSVC
if(MSVC)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-openmp:llvm" COMPILER_SUPPORTS_OPENMP_LLVM)
  if(COMPILER_SUPPORTS_OPENMP_LLVM)
    target_compile_options(${PROJECT_NAME} PRIVATE -openmp:llvm)
    target_link_options(${PROJECT_NAME} PRIVATE -openmp:llvm)
  endif()
endif()

add_dependencies(ViennaLS_Tests ${PROJECT_NAME})
add_test(NAME ${PROJECT_NAME} COMMAND $<TARGET_FILE:${PROJECT_NAME}>)
