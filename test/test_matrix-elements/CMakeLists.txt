
if(USE_QE)

  add_test(
    NAME test_matrix-elements/Cu-fullmesh
    COMMAND run-test.bash
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Cu-fullmesh
  )
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/Cu-fullmesh/qe/matrix-element-test.bash" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Cu-fullmesh/qe/Makefile.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/Cu-fullmesh/qe/Makefile" @ONLY)


  add_test(
    NAME test_matrix-elements/Si-fullmesh-off-diagonal
    COMMAND run-test.bash
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh-off-diagonal
  )
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh-off-diagonal/qe/matrix-element-test.bash" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh-off-diagonal/qe/Makefile.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh-off-diagonal/qe/Makefile" @ONLY)


  add_test(
    NAME test_matrix-elements/Si-sym
    COMMAND run-test.bash
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Si-sym
  )
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/qe/matrix-element-test.bash" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/qe/Makefile.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/qe/Makefile" @ONLY)


  add_test(
    NAME test_matrix-elements/GaAs-sym
    COMMAND run-test.bash
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GaAs-sym
  )
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-sym/qe/matrix-element-test.bash" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-sym/qe/Makefile.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-sym/qe/Makefile" @ONLY)


  add_test(
    NAME test_matrix-elements/Cu-qpoint-fullmesh
    COMMAND run-test.bash
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-fullmesh
  )
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-fullmesh/qe/matrix-element-test.bash" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-fullmesh/qe/Makefile.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-fullmesh/qe/Makefile" @ONLY)


  add_test(
    NAME test_matrix-elements/Cu-qpoint-sym
    COMMAND run-test.bash
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-sym
  )
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-sym/qe/matrix-element-test.bash" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-sym/qe/Makefile.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/Cu-qpoint-sym/qe/Makefile" @ONLY)


  add_test(
    NAME test_matrix-elements/Si-qpoint-sym
    COMMAND run-test.bash
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Si-qpoint-sym
  )
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/Si-qpoint-sym/qe/matrix-element-test.bash" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-qpoint-sym/qe/Makefile.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/Si-qpoint-sym/qe/Makefile" @ONLY)


  add_test(
    NAME test_matrix-elements/GaAs-qpoint-sym
    COMMAND run-test.bash
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-sym
  )
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-sym/qe/matrix-element-test.bash" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-sym/qe/Makefile.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-sym/qe/Makefile" @ONLY)


  add_test(
    NAME test_matrix-elements/GaAs-soc-fullmesh
    COMMAND bash -c "if [ \"\${INTW_RUN_LONG_TESTS:-}\" -ne 0 ]; then ./run-test.bash; else exit 127; fi"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-fullmesh
  )
  set_tests_properties(test_matrix-elements/GaAs-soc-fullmesh PROPERTIES SKIP_RETURN_CODE 127)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-fullmesh/qe/matrix-element-test.bash" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-fullmesh/qe/Makefile.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-fullmesh/qe/Makefile" @ONLY)


  add_test(
    NAME test_matrix-elements/GaAs-soc-sym
    COMMAND run-test.bash
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-sym
  )
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-sym/qe/matrix-element-test.bash" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-sym/qe/Makefile.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc-sym/qe/Makefile" @ONLY)

  add_test(
    NAME test_matrix-elements/GaAs-qpoint-soc-sym
    COMMAND run-test.bash
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-soc-sym
  )
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/matrix-element-test.bash.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-soc-sym/qe/matrix-element-test.bash" @ONLY)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-soc-sym/qe/Makefile.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-qpoint-soc-sym/qe/Makefile" @ONLY)

else()
  message(STATUS "test_matrix-elements is disabled")
  set(DISABLED_TESTS "${DISABLED_TESTS} test_matrix-elements" CACHE INTERNAL "")
endif()