
set(SIESTA2PH_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION} CACHE PATH "Path to the siesta2ph src directory." FORCE)

set(SIESTA2PH_SRC_FILES ${SIESTA2PH_SRC_DIR}/siesta2ph_io.f90
                        ${SIESTA2PH_SRC_DIR}/siesta2ph_symmetry.f90
                        ${SIESTA2PH_SRC_DIR}/siesta2ph_utils.f90
                        ${SIESTA2PH_SRC_DIR}/siesta2ph_linalg.F90
                        ${SIESTA2PH_SRC_DIR}/siesta2ph_system.f90
                        ${SIESTA2PH_SRC_DIR}/spglib_f08.f90
                        )


# Add some debugging flags
string(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
if(BT STREQUAL "DEBUG")
  set(SIESTA2PH_FFLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
  set(SIESTA2PH_FPPFLAGS -DDEBUG)
endif()


if("${SIESTA_CMAKE}" STREQUAL "SIESTA_CMAKE-NOTFOUND")

  set(SIESTA2PH_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${SIESTA_VERSION} CACHE PATH "Path to the siesta2ph build directory." FORCE)


  set(SIESTA2PH_OBJ_FILES ${SIESTA2PH_BUILD_DIR}/siesta2ph_io.o
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_linalg.o
                          ${SIESTA2PH_BUILD_DIR}/spglib_f08.o
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_symmetry.o
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_system.o
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_utils.o
                          )
  set(SIESTA2PH_MOD_FILES ${SIESTA2PH_BUILD_DIR}/siesta2ph_io.mod
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_linalg.mod
                          ${SIESTA2PH_BUILD_DIR}/spglib_f08.mod
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_symmetry.mod
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_system.mod
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_utils.mod
                          )


  list(JOIN SPGLIB_LIBRARIES " " SPGLIB_LIBRARIES_STR)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/Makefile.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/${SIESTA_VERSION}/Makefile" @ONLY)

  add_custom_command(OUTPUT ${SIESTA2PH_BUILD_DIR}/siesta2ph.x
    COMMAND make -C ${SIESTA2PH_BUILD_DIR} siesta2ph.x
    COMMENT "Building siesta2ph.x"
    DEPENDS ${SIESTA2PH_SRC_FILES} ${SIESTA2PH_SRC_DIR}/siesta2ph.F90
    BYPRODUCTS ${SIESTA2PH_OBJ_FILES} ${SIESTA2PH_MOD_FILES} ${SIESTA2PH_BUILD_DIR}/siesta2ph.o
    )

  add_custom_command(OUTPUT ${SIESTA2PH_BUILD_DIR}/siesta2fc.x
    COMMAND make -C ${SIESTA2PH_BUILD_DIR} siesta2fc.x
    COMMENT "Building siesta2fc.x"
    DEPENDS ${SIESTA2PH_SRC_FILES} ${SIESTA2PH_SRC_DIR}/siesta2fc.F90
    BYPRODUCTS ${SIESTA2PH_OBJ_FILES} ${SIESTA2PH_MOD_FILES} ${SIESTA2PH_BUILD_DIR}/siesta2fc.o
    )

  add_custom_command(OUTPUT ${SIESTA2PH_BUILD_DIR}/siesta2dv.x
    COMMAND make -C ${SIESTA2PH_BUILD_DIR} siesta2dv.x
    COMMENT "Building siesta2dv.x"
    DEPENDS ${SIESTA2PH_SRC_FILES} ${SIESTA2PH_SRC_DIR}/siesta2dv.F90
    BYPRODUCTS ${SIESTA2PH_OBJ_FILES} ${SIESTA2PH_MOD_FILES} ${SIESTA2PH_BUILD_DIR}/siesta2dv.o
    )

  add_custom_target(siesta2ph.x ALL
    DEPENDS ${SIESTA2PH_BUILD_DIR}/siesta2ph.x
    )

  add_custom_target(siesta2fc.x ALL
    DEPENDS ${SIESTA2PH_BUILD_DIR}/siesta2fc.x
    )

  add_custom_target(siesta2dv.x ALL
    DEPENDS ${SIESTA2PH_BUILD_DIR}/siesta2dv.x
    )

else()

  set(SIESTA2PH_BUILD_DIR ${SIESTA_BUILD_DIR}/Util/Contrib/siesta2ph CACHE PATH "Path to the siesta2ph build directory." FORCE)

  set(SIESTA2PH_OBJ_FILES ${SIESTA2PH_BUILD_DIR}/siesta2ph_io.o
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_linalg.o
                          ${SIESTA2PH_BUILD_DIR}/spglib_f08.o
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_symmetry.o
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_system.o
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_utils.o
                          )
  set(SIESTA2PH_MOD_FILES ${SIESTA2PH_BUILD_DIR}/siesta2ph_io.mod
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_linalg.mod
                          ${SIESTA2PH_BUILD_DIR}/spglib_f08.mod
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_symmetry.mod
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_system.mod
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_utils.mod
                          )


  list(APPEND SIESTA2PH_SRC_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/siesta2ph.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/siesta2fc.F90
    ${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/siesta2dv.F90
    )

  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/CMakeLists.txt.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/CMakeLists.txt" @ONLY)


  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp
    COMMAND mkdir -p "${SIESTA_HOME}/Util/Contrib/siesta2ph"
    COMMAND cp ${SIESTA2PH_SRC_FILES} "${SIESTA_HOME}/Util/Contrib/siesta2ph/"
    COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/CMakeLists.txt "${SIESTA_HOME}/Util/Contrib/siesta2ph/"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp"
    COMMENT "Copying siesta2ph to Util/Contrib"
    DEPENDS ${SIESTA2PH_SRC_FILES}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/CMakeLists.txt
    )

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp
    COMMAND bash -c "if ! grep -q siesta2ph ${SIESTA_HOME}/Util/CMakeLists.txt; then echo \"add_subdirectory(\\\"Contrib/siesta2ph\\\")\" >> ${SIESTA_HOME}/Util/CMakeLists.txt; fi"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp"
    COMMENT "Modifying CMakeLists.txt from SIESTA/Util"
    VERBATIM
    )

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp
    COMMAND cmake -S "${SIESTA_HOME}" -B "${SIESTA_BUILD_DIR}"
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp
    COMMENT "Configuring SIESTA"
    )

  add_custom_target(siesta2ph.x ALL cmake --build "${SIESTA_BUILD_DIR}" -t siesta2ph_exe
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp
    BYPRODUCTS ${SIESTA_BUILD_DIR}/Util/Contrib/siesta2ph/siesta2ph.x
    )

  add_custom_target(siesta2fc.x ALL cmake --build "${SIESTA_BUILD_DIR}" -t siesta2fc_exe
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp
    BYPRODUCTS ${SIESTA_BUILD_DIR}/Util/Contrib/siesta2ph/siesta2fc.x
    )

  add_custom_target(siesta2dv.x ALL cmake --build "${SIESTA_BUILD_DIR}" -t siesta2dv_exe
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp
    BYPRODUCTS ${SIESTA_BUILD_DIR}/Util/Contrib/siesta2ph/siesta2dv.x
    )

  add_custom_target(restore-siesta2
    COMMAND rm -rf "${SIESTA_BUILD_DIR}/Util/Contrib/siesta2ph"
    COMMAND rm -rf "${SIESTA_HOME}/Util/Contrib/siesta2ph"
    COMMAND mv "${SIESTA_HOME}/Util/CMakeLists.txt" "${SIESTA_HOME}/Util/CMakeLists.txt_backup"
    COMMAND grep -v siesta2ph "${SIESTA_HOME}/Util/CMakeLists.txt_backup" > "${SIESTA_HOME}/Util/CMakeLists.txt"
    COMMAND rm -rf "${SIESTA_HOME}/Util/CMakeLists.txt_backup"
    COMMAND cmake -S "${SIESTA_HOME}" -B "${SIESTA_BUILD_DIR}"
    COMMAND rm -rf "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp"
    COMMENT "SIESTA restored"
    VERBATIM
  )

  add_dependencies(restore-siesta restore-siesta2)


endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/utils/run-disps.bash.in"
                "${CMAKE_CURRENT_BINARY_DIR}/${SIESTA_VERSION}/run-disps.x" @ONLY)
