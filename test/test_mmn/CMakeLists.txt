

## JLB: mmn
add_test(
  NAME test_mmn/Si-fullmesh
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh/mmn-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh/QE-6.7.0/mmn-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_mmn/Si-sym
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Si-sym
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/mmn-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/QE-6.7.0/mmn-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_mmn/GaAs-soc
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc/mmn-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc/QE-6.7.0/mmn-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc/QE-6.7.0/Makefile" @ONLY)

# Setup dependencies
set_tests_properties(test_mmn/Si-fullmesh PROPERTIES FIXTURES_SETUP test_mmn_Si-sym_DEPENDENCIES) # Define test_mmn_Si-sym_DEPENDENCIES
set_tests_properties(test_mmn/Si-sym PROPERTIES FIXTURES_REQUIRED test_mmn_Si-sym_DEPENDENCIES)
