

include(IntwTesting)

if(USE_MPI_ON_TESTS)
  set(MPI_CMD "mpiexec -np 4")
else()
  set(MPI_CMD "")
endif()

add_library(test_modules STATIC module/test_module.f90)
target_link_libraries(test_modules PUBLIC intw_modules)



# Adding a test that uses a python script:
add_test(
  NAME test_dummy/simple_test.py
  COMMAND simple_test.py
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_dummy
)

# NAME could be whatever.

# COMMAND must point to the script itself (relative to WORKING_DIRECTORY), and
# could also contain arguments.

# WORKING_DIRECTORY specifies the working directory in which to execute the test.
# I think that it is convenient to be the same directory where the script is located,
# as all input and output files would be found there.


# Here I add a simple bash script that uses command line arguments:
add_test(
  NAME test_dummy/simple_test.bash
  COMMAND simple_test.bash 3 5
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_dummy
)


# Adding a test written in fortran is done with the add_fortran_test_executable function,
# which is defined in above.
# In this case there are some particularities:
# 1- The first argument of add_fortran_test_executable must be the directory where
#    the test is located, relative to the test directory.
# 2- If the file that contains the test is test_dummy/test.f90, then, the name of
#    the function written there must be test_dummy_test, i.e. the same path, but
#    replacing "/" with "_", and without .f90.
# 3- The NAME of the test in this case will be test_dummy/test.f90, and there is no
#    option to change it. Thus, the name of the test file should be meaningful.
add_fortran_test_executable(
  test_dummy
test_dummy/simple_test.f90)

# If more than one test are present in the same directory, it is enough to add it
# in an additional line:
add_fortran_test_executable(
  test_IO
  test_IO/read_input.f90
  test_IO/read_data.f90
)


add_test(
  NAME test_matrix-elements/Fe-O
  COMMAND bash -c "if [ \"\${INTW_RUN_LONG_TESTS:-}\" -ne 0 ]; then ./run-test.bash; else exit 127; fi"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/Fe-O
)
set_tests_properties(test_matrix-elements/Fe-O PROPERTIES SKIP_RETURN_CODE 127)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/Fe-O/matrix-element-test.bash.in"
	"${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/Fe-O/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/Fe-O/QE-6.7.0/Makefile.in"
"${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/Fe-O/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/Cu-fullmesh
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/Cu-fullmesh
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/Cu-fullmesh/matrix-element-test.bash.in"
	"${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/Cu-fullmesh/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/Cu-fullmesh/QE-6.7.0/Makefile.in"
	"${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/Cu-fullmesh/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/Si-fullmesh-off-diagonal
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/Si-fullmesh-off-diagonal
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/Si-fullmesh-off-diagonal/matrix-element-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/Si-fullmesh-off-diagonal/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/Si-fullmesh-off-diagonal/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/Si-fullmesh-off-diagonal/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/GaAs-fullmesh
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/GaAs-fullmesh
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/GaAs-fullmesh/matrix-element-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/GaAs-fullmesh/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/GaAs-fullmesh/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/test_matrix-elements/GaAs-fullmesh/QE-6.7.0/Makefile" @ONLY)
