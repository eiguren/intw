

## JLB: w902intw
add_test(
  NAME test_w902intw/Si-isolated
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Si-isolated
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-isolated/w902intw-isolated-test.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-isolated/QE-6.7.0/w902intw-isolated-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-isolated/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-isolated/QE-6.7.0/Makefile" @ONLY)


## Setup dependencies
#set_tests_properties(test_w902intw/Si-isolated PROPERTIES FIXTURES_SETUP test_w902intw_Si-isolated_DEPENDENCIES) # Define test_w902intw_Si-isolated_DEPENDENCIES
#set_tests_properties(test_w902intw/Si-isolated PROPERTIES FIXTURES_REQUIRED test_w902intw_Si-isolated_DEPENDENCIES)
