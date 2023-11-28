

add_test(
  NAME test_matrix-elements/Cu-fullmesh
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Cu-fullmesh
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
	"${CMAKE_CURRENT_SOURCE_DIR}/Cu-fullmesh/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Cu-fullmesh/QE-6.7.0/Makefile.in"
	"${CMAKE_CURRENT_SOURCE_DIR}/Cu-fullmesh/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/Si-fullmesh-off-diagonal
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh-off-diagonal
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh-off-diagonal/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh-off-diagonal/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh-off-diagonal/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/Si-sym
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Si-sym
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/GaAs-sym
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GaAs-sym
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-sym/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-sym/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-sym/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/Cu-qpoint-fullmesh
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-fullmesh
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-fullmesh/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-fullmesh/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-fullmesh/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/Cu-qpoint-sym
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-sym
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-sym/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-sym/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-sym/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/Si-qpoint-sym
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Si-qpoint-sym
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-qpoint-sym/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-qpoint-sym/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-qpoint-sym/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/GaAs-qpoint-sym
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-sym
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-sym/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-sym/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-sym/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/GaAs-soc-fullmesh
  COMMAND bash -c "if [ \"\${INTW_RUN_LONG_TESTS:-}\" -ne 0 ]; then ./run-test.bash; else exit 127; fi"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-fullmesh
)
set_tests_properties(test_matrix-elements/GaAs-soc-fullmesh PROPERTIES SKIP_RETURN_CODE 127)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-fullmesh/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-fullmesh/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-fullmesh/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/GaAs-soc-sym
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-sym
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-sym/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-sym/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-sym/QE-6.7.0/Makefile" @ONLY)

add_test(
  NAME test_matrix-elements/GaAs-qpoint-soc-sym
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-soc-sym
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-soc-sym/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-soc-sym/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-soc-sym/QE-6.7.0/Makefile" @ONLY)
