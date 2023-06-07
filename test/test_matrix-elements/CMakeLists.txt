
add_test(
  NAME test_matrix-elements/Fe-O
  COMMAND bash -c "if [ \"\${INTW_RUN_LONG_TESTS:-}\" -ne 0 ]; then ./run-test.bash; else exit 127; fi"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Fe-O
)
set_tests_properties(test_matrix-elements/Fe-O PROPERTIES SKIP_RETURN_CODE 127)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Fe-O/matrix-element-test.bash.in"
	"${CMAKE_CURRENT_SOURCE_DIR}/Fe-O/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Fe-O/QE-6.7.0/Makefile.in"
"${CMAKE_CURRENT_SOURCE_DIR}/Fe-O/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/Cu-fullmesh
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Cu-fullmesh
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Cu-fullmesh/matrix-element-test.bash.in"
	"${CMAKE_CURRENT_SOURCE_DIR}/Cu-fullmesh/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Cu-fullmesh/QE-6.7.0/Makefile.in"
	"${CMAKE_CURRENT_SOURCE_DIR}/Cu-fullmesh/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/Si-fullmesh-off-diagonal
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh-off-diagonal
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh-off-diagonal/matrix-element-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh-off-diagonal/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh-off-diagonal/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh-off-diagonal/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/Si-sym
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Si-sym
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/matrix-element-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/QE-6.7.0/Makefile" @ONLY)


add_test(
  NAME test_matrix-elements/GaAs-sym
  COMMAND run-test.bash
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GaAs-sym
)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-sym/matrix-element-test.bash.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-sym/QE-6.7.0/matrix-element-test.bash" @ONLY)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-sym/QE-6.7.0/Makefile.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-sym/QE-6.7.0/Makefile" @ONLY)
