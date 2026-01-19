# Add local libtorch directory to search path
list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../../libtorch")

find_package(Torch QUIET)

if(Torch_FOUND)
  add_executable(NeuralAdvection NeuralAdvection.cpp)
  target_link_libraries(NeuralAdvection PRIVATE ViennaTools::ViennaLS ${TORCH_LIBRARIES})
  set_property(TARGET NeuralAdvection PROPERTY CXX_STANDARD 17)
  add_custom_command(TARGET NeuralAdvection POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy_if_different
                     ${CMAKE_CURRENT_SOURCE_DIR}/sdf_super_res.pt
                     $<TARGET_FILE_DIR:NeuralAdvection>)
else()
  message(STATUS "[ViennaLS] LibTorch not found. Skipping NeuralAdvection example.")
endif()
