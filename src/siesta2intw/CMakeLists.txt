
# Define siesta2intw source directories and files
set(SIESTA2INTW_MODULES_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/modules CACHE PATH "Path to the siesta2intw modules src directory." FORCE)
set(SIESTA2INTW_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION} CACHE PATH "Path to the siesta2intw src directory." FORCE)
set(SIESTA2INTW_SRC_FILES ${SIESTA2INTW_MODULES_SRC_DIR}/cfftnd.f90
                          ${SIESTA2INTW_SRC_DIR}/siesta2intw.F
                          ${SIESTA2INTW_MODULES_SRC_DIR}/siesta2intw_io.f90
                          ${SIESTA2INTW_MODULES_SRC_DIR}/siesta2intw_ph.f90
                          ${SIESTA2INTW_MODULES_SRC_DIR}/siesta2intw_symmetry.f90
                          ${SIESTA2INTW_MODULES_SRC_DIR}/siesta2intw_write.F90
                          ${SIESTA2INTW_MODULES_SRC_DIR}/siesta2intw_fft.f90
                          ${SIESTA2INTW_MODULES_SRC_DIR}/siesta2intw_pp.f90
                          ${SIESTA2INTW_MODULES_SRC_DIR}/siesta2intw_utils.F90
                          )

# Define siesta2intw build directory and object and module files
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
                            )
endif()
set(SIESTA2INTW_MOD_FILES ${SIESTA2INTW_BUILD_DIR}/fftpack5.mod
                          ${SIESTA2INTW_BUILD_DIR}/siesta2intw_fft.mod
                          ${SIESTA2INTW_BUILD_DIR}/siesta2intw_io.mod
                          ${SIESTA2INTW_BUILD_DIR}/siesta2intw_ph.mod
                          ${SIESTA2INTW_BUILD_DIR}/siesta2intw_pp.mod
                          ${SIESTA2INTW_BUILD_DIR}/siesta2intw_symmetry.mod
                          ${SIESTA2INTW_BUILD_DIR}/siesta2intw_utils.mod
                          ${SIESTA2INTW_BUILD_DIR}/siesta2intw_write.mod
                          )
# Define siesta2intw executable (only needed for testing)
set(SIESTA2INTW_EXE ${SIESTA2INTW_BUILD_DIR}/siesta2intw.x CACHE PATH "Path to the siesta2intw.x executable." FORCE)


# Add some flags
set(SIESTA2INTW_FFLAGS "${CMAKE_Fortran_FLAGS}")
#
string(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
if(BT STREQUAL "DEBUG")
  set(SIESTA2INTW_FFLAGS "${SIESTA2INTW_FFLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
  set(SIESTA2INTW_FPPFLAGS -DDEBUG)
endif()


set(OUTPUT_PREFIX "...... ")

if("${SIESTA_CMAKE}" STREQUAL "SIESTA_CMAKE-NOTFOUND")

  # Create Makefile for siesta2intw
  set(SPGLIB_LIBRARIES_LIST "-Wl,--start-group" ${SPGLIB_LIBRARIES} "-Wl,--end-group")
  list(JOIN SPGLIB_LIBRARIES_LIST " " SPGLIB_LIBRARIES_STR)

  set(SPGLIB_INCLUDE_DIRS_LIST "" ${SPGLIB_INCLUDE_DIRS})
  list(JOIN SPGLIB_INCLUDE_DIRS_LIST " -I" SPGLIB_INCLUDE_FLAGS_STR)

  configure_file("${SIESTA2INTW_SRC_DIR}/Makefile.in"
                 "${SIESTA2INTW_BUILD_DIR}/Makefile" @ONLY)

  # Create the targets for the siesta2intw executable
  add_custom_command(OUTPUT ${SIESTA2INTW_EXE}
    COMMAND make -C ${SIESTA2INTW_BUILD_DIR} siesta2intw.x | sed 's/^/${OUTPUT_PREFIX}/'
    COMMENT "Building siesta2intw.x"
    DEPENDS ${SIESTA2INTW_SRC_FILES}
            Spglib::symspg
            Spglib::fortran
    BYPRODUCTS ${SIESTA2INTW_OBJ_FILES} ${SIESTA2INTW_MOD_FILES}
    )

else()

  # Create the cmake configuration file for siesta2intw
  configure_file("${SIESTA2INTW_SRC_DIR}/CMakeLists.txt.in"
                 "${SIESTA2INTW_SRC_DIR}/CMakeLists.txt" @ONLY)

  # Create the targets to patch siesta2intw into ${SIESTA_HOME}/Util/Contrib/siesta2intw
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp
    COMMAND mkdir -p "${SIESTA_HOME}/Util/Contrib/siesta2intw"
    COMMAND cp ${SIESTA2INTW_SRC_FILES} "${SIESTA_HOME}/Util/Contrib/siesta2intw/"
    COMMAND cp ${SIESTA2INTW_SRC_DIR}/CMakeLists.txt "${SIESTA_HOME}/Util/Contrib/siesta2intw/"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp"
    COMMENT "Copying siesta2intw to Util/Contrib"
    DEPENDS ${SIESTA2INTW_SRC_FILES}
    DEPENDS ${SIESTA2INTW_SRC_DIR}/CMakeLists.txt
    )

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp
    COMMAND bash -c "if ! grep -q siesta2intw ${SIESTA_HOME}/Util/CMakeLists.txt; then echo \"add_subdirectory(\\\"Contrib/siesta2intw\\\")\" >> ${SIESTA_HOME}/Util/CMakeLists.txt; fi"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp"
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp
    COMMENT "Modifying CMakeLists.txt from SIESTA/Util"
    VERBATIM
    )

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/patch.stamp
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/patch.stamp
    COMMAND touch ${PROJECT_BINARY_DIR}/CMakeFiles/patch-siesta.stamp # This is needed to reconfigure SIESTA
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp
    COMMENT ""
    BYPRODUCTS ${PROJECT_BINARY_DIR}/CMakeFiles/patch-siesta.stamp
    )

  add_custom_target(patch-siesta-siesta2intw
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/patch.stamp
    )

  add_custom_target(patch-siesta)

  add_dependencies(patch-siesta patch-siesta-siesta2intw)

  # Create targets to configure SIESTA
  add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/CMakeFiles/configure-siesta.stamp
    COMMAND cmake -S "${SIESTA_HOME}" -B "${SIESTA_BUILD_DIR}" | sed 's/^/${OUTPUT_PREFIX}/'
    COMMAND touch ${PROJECT_BINARY_DIR}/CMakeFiles/configure-siesta.stamp
    DEPENDS ${PROJECT_BINARY_DIR}/CMakeFiles/patch-siesta.stamp
            patch-siesta
    COMMENT "Configuring SIESTA"
    )

  add_custom_target(configure-siesta
    DEPENDS ${PROJECT_BINARY_DIR}/CMakeFiles/configure-siesta.stamp
    )

  # Create the targets for the siesta2intw executables
  add_custom_command(OUTPUT ${SIESTA2INTW_EXE}
    COMMAND cmake --build "${SIESTA_BUILD_DIR}" -t siesta2intw_exe | sed 's/^/${OUTPUT_PREFIX}/'
    COMMENT "Building siesta2intw.x"
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp
            patch-siesta
            configure-siesta
            Spglib::symspg
            Spglib::fortran
    BYPRODUCTS ${SIESTA2INTW_OBJ_FILES} ${SIESTA2INTW_MOD_FILES}
    )

  # Create a target to remove all siesta2intw files and restore siesta
  add_custom_target(remove-siesta2intw
    COMMENT "Removing siesta2intw from SIESTA"
    COMMAND rm -rf "${SIESTA_BUILD_DIR}/Util/Contrib/siesta2intw"
    COMMAND rm -rf "${SIESTA_HOME}/Util/Contrib/siesta2intw"
    COMMAND mv "${SIESTA_HOME}/Util/CMakeLists.txt" "${SIESTA_HOME}/Util/CMakeLists.txt_backup"
    COMMAND grep -v siesta2intw "${SIESTA_HOME}/Util/CMakeLists.txt_backup" > "${SIESTA_HOME}/Util/CMakeLists.txt"
    COMMAND rm -rf "${SIESTA_HOME}/Util/CMakeLists.txt_backup"
    COMMAND rm -rf "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/patch.stamp" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp"
    VERBATIM
    )

  add_custom_target(restore-siesta
    COMMAND cmake -S "${SIESTA_HOME}" -B "${SIESTA_BUILD_DIR}"
    COMMAND rm -rf "${PROJECT_BINARY_DIR}/CMakeFiles/configure-siesta.stamp" "${PROJECT_BINARY_DIR}/CMakeFiles/patch-siesta.stamp"
    COMMENT "Restoring SIESTA"
    )

  add_dependencies(restore-siesta remove-siesta2intw)

endif()

# Create the targets to build siesta2intw
add_custom_target(siesta2intw.x DEPENDS ${SIESTA2INTW_EXE})

add_custom_target(siesta2intw ALL DEPENDS siesta2intw.x)
