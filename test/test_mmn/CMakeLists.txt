

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
  NAME test_mmn/Pb-sym
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Pb-sym
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Pb-sym/mmn-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Pb-sym/QE-6.7.0/mmn-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Pb-sym/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Pb-sym/QE-6.7.0/Makefile" @ONLY)
