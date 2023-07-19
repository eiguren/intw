


include(IntwTesting)

if(USE_MPI_ON_TESTS)
  set(MPI_CMD "mpiexec -np 4")
else()
  set(MPI_CMD "")
endif()

add_library(test_modules STATIC module/test_module.f90)
target_include_directories( test_modules INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(test_modules PUBLIC intw_modules)


add_subdirectory(test_dummy)

add_subdirectory(test_IO)

add_subdirectory(test_mmn)

add_subdirectory(test_matrix-elements)

