

set(SIESTA2INTW_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION} CACHE PATH "Path to the siesta2intw src directory." FORCE)

set(SIESTA2INTW_SRC_FILES ${SIESTA2INTW_SRC_DIR}/cfftnd.f90
                          ${SIESTA2INTW_SRC_DIR}/siesta2intw.F
                          ${SIESTA2INTW_SRC_DIR}/siesta2intw_io.f90
                          ${SIESTA2INTW_SRC_DIR}/siesta2intw_ph.f90
                          ${SIESTA2INTW_SRC_DIR}/siesta2intw_symmetry.f90
                          ${SIESTA2INTW_SRC_DIR}/spglib_f08.f90
                          ${SIESTA2INTW_SRC_DIR}/siesta2intw_write.F90
                          ${SIESTA2INTW_SRC_DIR}/siesta2intw_fft.f90
                          ${SIESTA2INTW_SRC_DIR}/siesta2intw_pp.f90
                          ${SIESTA2INTW_SRC_DIR}/siesta2intw_utils.F90
                          )


# Add some debugging flags
string(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
if(BT STREQUAL "DEBUG")
  set(SIESTA2INTW_FFLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
  set(SIESTA2INTW_FPPFLAGS -DDEBUG)
endif()


if("${SIESTA_CMAKE}" STREQUAL "SIESTA_CMAKE-NOTFOUND")

  set(SIESTA2INTW_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${SIESTA_VERSION} CACHE PATH "Path to the siesta2intw build directory." FORCE)


  set(SIESTA2INTW_OBJ_FILES ${SIESTA2INTW_BUILD_DIR}/cfftnd.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_fft.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_io.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_ph.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_pp.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_symmetry.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_utils.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_write.o
                            ${SIESTA2INTW_BUILD_DIR}/spglib_f08.o
                            )
  set(SIESTA2INTW_MOD_FILES ${SIESTA2INTW_BUILD_DIR}/fftpack5.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_fft.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_io.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_ph.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_pp.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_symmetry.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_utils.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_write.mod
                            ${SIESTA2INTW_BUILD_DIR}/spglib_f08.mod
                            )



  list(JOIN SPGLIB_LIBRARIES " " SPGLIB_LIBRARIES_STR)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/Makefile.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/${SIESTA_VERSION}/Makefile" @ONLY)

  add_custom_command(OUTPUT ${SIESTA2INTW_EXE}
    COMMAND make -C ${SIESTA2INTW_BUILD_DIR} siesta2intw.x
    COMMENT "Building siesta2intw"
    DEPENDS ${SIESTA2INTW_SRC_FILES}
    BYPRODUCTS ${SIESTA2INTW_OBJ_FILES} ${SIESTA2INTW_MOD_FILES}
    )

  add_custom_target(siesta2intw.x ALL
    DEPENDS ${SIESTA2INTW_EXE}
    )

else()

  set(SIESTA2INTW_BUILD_DIR ${SIESTA_BUILD_DIR}/Util/Contrib/siesta2intw CACHE PATH "Path to the siesta2intw build directory." FORCE)

  set(SIESTA2INTW_OBJ_FILES ${SIESTA2INTW_BUILD_DIR}/CMakeFiles/siesta2intw_exe.dir/cfftnd.f90.o
                            ${SIESTA2INTW_BUILD_DIR}/CMakeFiles/siesta2intw_exe.dir/siesta2intw.F.o
                            ${SIESTA2INTW_BUILD_DIR}/CMakeFiles/siesta2intw_exe.dir/siesta2intw_fft.f90.o
                            ${SIESTA2INTW_BUILD_DIR}/CMakeFiles/siesta2intw_exe.dir/siesta2intw_io.f90.o
                            ${SIESTA2INTW_BUILD_DIR}/CMakeFiles/siesta2intw_exe.dir/siesta2intw_ph.f90.o
                            ${SIESTA2INTW_BUILD_DIR}/CMakeFiles/siesta2intw_exe.dir/siesta2intw_pp.f90.o
                            ${SIESTA2INTW_BUILD_DIR}/CMakeFiles/siesta2intw_exe.dir/siesta2intw_symmetry.f90.o
                            ${SIESTA2INTW_BUILD_DIR}/CMakeFiles/siesta2intw_exe.dir/siesta2intw_utils.F90.o
                            ${SIESTA2INTW_BUILD_DIR}/CMakeFiles/siesta2intw_exe.dir/siesta2intw_write.F90.o
                            ${SIESTA2INTW_BUILD_DIR}/CMakeFiles/siesta2intw_exe.dir/spglib_f08.f90.o
                            )
  set(SIESTA2INTW_MOD_FILES ${SIESTA2INTW_BUILD_DIR}/fftpack5.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_fft.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_io.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_ph.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_pp.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_symmetry.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_utils.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_write.mod
                            ${SIESTA2INTW_BUILD_DIR}/spglib_f08.mod
                            )


  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/CMakeLists.txt.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/CMakeLists.txt" @ONLY)


  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp
    COMMAND mkdir -p "${SIESTA_HOME}/Util/Contrib/siesta2intw"
    COMMAND cp ${SIESTA2INTW_SRC_FILES} "${SIESTA_HOME}/Util/Contrib/siesta2intw/"
    COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/CMakeLists.txt "${SIESTA_HOME}/Util/Contrib/siesta2intw/"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp"
    COMMENT "Copying siesta2intw to Util/Contrib"
    DEPENDS ${SIESTA2INTW_SRC_FILES}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/CMakeLists.txt
    )

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp
    COMMAND bash -c "if ! grep -q siesta2intw ${SIESTA_HOME}/Util/CMakeLists.txt; then echo \"add_subdirectory(\\\"Contrib/siesta2intw\\\")\" >> ${SIESTA_HOME}/Util/CMakeLists.txt; fi"
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

  add_custom_command(OUTPUT ${SIESTA2INTW_BUILD_DIR}/siesta2intw.x
    COMMAND cmake --build "${SIESTA_BUILD_DIR}" -t siesta2intw_exe
    COMMENT "Building siesta2intw.x"
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp
    BYPRODUCTS ${SIESTA2INTW_OBJ_FILES} ${SIESTA2INTW_MOD_FILES}
    )

  add_custom_target(siesta2intw.x ALL
    DEPENDS ${SIESTA2INTW_BUILD_DIR}/siesta2intw.x
    )

  add_custom_target(remove-siesta2intw
    COMMENT "Removing siesta2intw from SIESTA"
    COMMAND rm -rf "${SIESTA_BUILD_DIR}/Util/Contrib/siesta2intw"
    COMMAND rm -rf "${SIESTA_HOME}/Util/Contrib/siesta2intw"
    COMMAND mv "${SIESTA_HOME}/Util/CMakeLists.txt" "${SIESTA_HOME}/Util/CMakeLists.txt_backup"
    COMMAND grep -v siesta2intw "${SIESTA_HOME}/Util/CMakeLists.txt_backup" > "${SIESTA_HOME}/Util/CMakeLists.txt"
    COMMAND rm -rf "${SIESTA_HOME}/Util/CMakeLists.txt_backup"
    COMMAND rm -rf "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp"
    VERBATIM
    )

  add_custom_target(restore-siesta
    COMMAND cmake -S "${SIESTA_HOME}" -B "${SIESTA_BUILD_DIR}"
    COMMENT "Restoring SIESTA"
    )

  add_dependencies(restore-siesta remove-siesta2intw)

endif()

set(SIESTA2INTW_EXE ${SIESTA2INTW_BUILD_DIR}/siesta2intw.x CACHE PATH "Path to the siesta2intw.x executable." FORCE)
