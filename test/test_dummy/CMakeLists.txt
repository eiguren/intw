

# Adding a test that uses a python script:
add_test(
  NAME test_dummy/simple_test.py
  COMMAND simple_test.py
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
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
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)


# Adding a test written in fortran is done with the add_fortran_test_executable function,
# which is defined in above.
# In this case there are some particularities:
# 1- The first argument of add_fortran_test_executable must be the directory where
#    the test is located, relative to the test directory, i.e. the current directory.
# 2- If the file that contains the test is filename.f90, then, the name of
#    the function written there must be filename_test.
# 3- The NAME of the test in this case will be test_dummy/simpe_test.f90, and there
#    is no option to change it. Thus, the folder and file names should be meaningful.
# 4- If more than one test are present in the same directory, it is enough to add it
#    in an additional line (see test_IO as an example).
add_fortran_test_executable(
  test_dummy
  simple_test.f90)

