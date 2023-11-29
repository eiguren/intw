
if(USE_WANNIER AND (USE_QE OR USE_SIESTA))

  add_test(
    NAME test_mmn/Si-fullmesh
    COMMAND run-test.bash
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh
  )
  if(USE_QE)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh/mmn-test.bash.in"
      "${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh/QE-6.7.0/mmn-test.bash" @ONLY)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh/QE-6.7.0/Makefile.in"
      "${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh/QE-6.7.0/Makefile" @ONLY)
  endif()
  if(USE_SIESTA)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh/mmn-test.bash.in"
      "${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh/siesta/mmn-test.bash" @ONLY)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh/siesta/Makefile.in"
      "${CMAKE_CURRENT_SOURCE_DIR}/Si-fullmesh/siesta/Makefile" @ONLY)
  endif()


  add_test(
    NAME test_mmn/Si-sym
    COMMAND run-test.bash
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Si-sym
  )
  # Setup dependencies
  set_tests_properties(test_mmn/Si-fullmesh PROPERTIES FIXTURES_SETUP test_mmn_Si-sym_DEPENDENCIES) # Define test_mmn_Si-sym_DEPENDENCIES
  set_tests_properties(test_mmn/Si-sym PROPERTIES FIXTURES_REQUIRED test_mmn_Si-sym_DEPENDENCIES)

  if(USE_QE)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/mmn-test.bash.in"
      "${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/QE-6.7.0/mmn-test.bash" @ONLY)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/QE-6.7.0/Makefile.in"
      "${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/QE-6.7.0/Makefile" @ONLY)
  endif()
  if(USE_SIESTA)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/mmn-test.bash.in"
      "${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/siesta/mmn-test.bash" @ONLY)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/siesta/Makefile.in"
      "${CMAKE_CURRENT_SOURCE_DIR}/Si-sym/siesta/Makefile" @ONLY)
  endif()


  if(USE_QE)
    add_test(
      NAME test_mmn/GaAs-soc
      COMMAND run-test.bash
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc
    )
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc/mmn-test.bash.in"
      "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc/QE-6.7.0/mmn-test.bash" @ONLY)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc/QE-6.7.0/Makefile.in"
      "${CMAKE_CURRENT_SOURCE_DIR}/GaAs-soc/QE-6.7.0/Makefile" @ONLY)
  endif()

else()
  message(STATUS "test_mmn is disabled")
  set(DISABLED_TESTS "${DISABLED_TESTS} test_mmn" CACHE INTERNAL "")
endif()
