
# Define siesta2ph source directories and files
set(SIESTA2PH_MODULES_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/modules CACHE PATH "Path to the siesta2ph modules src directory." FORCE)
set(SIESTA2PH_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION} CACHE PATH "Path to the siesta2ph src directory." FORCE)
set(SIESTA2PH_MODULES_SRC_FILES ${SIESTA2PH_MODULES_SRC_DIR}/siesta2ph_io.f90
                                ${SIESTA2PH_MODULES_SRC_DIR}/siesta2ph_symmetry.f90
                                ${SIESTA2PH_MODULES_SRC_DIR}/siesta2ph_utils.f90
                                ${SIESTA2PH_MODULES_SRC_DIR}/siesta2ph_linalg.F90
                                ${SIESTA2PH_MODULES_SRC_DIR}/siesta2ph_system.f90
                                ${SIESTA2PH_MODULES_SRC_DIR}/spglib_f08.f90
                                )
set(SIESTA2PH_SRC_FILES ${SIESTA2PH_MODULES_SRC_FILES}
                        ${SIESTA2PH_SRC_DIR}/siesta2ph.F90
                        ${SIESTA2PH_SRC_DIR}/siesta2fc.F90
                        ${SIESTA2PH_SRC_DIR}/siesta2dv.F90
                        )

# Define siesta2intw build directory and object and module files
if("${SIESTA_CMAKE}" STREQUAL "SIESTA_CMAKE-NOTFOUND")
  set(SIESTA2PH_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${SIESTA_VERSION} CACHE PATH "Path to the siesta2ph build directory." FORCE)
  set(SIESTA2PH_MODULES_OBJ_FILES ${SIESTA2PH_BUILD_DIR}/siesta2ph_io.o
                                  ${SIESTA2PH_BUILD_DIR}/siesta2ph_symmetry.o
                                  ${SIESTA2PH_BUILD_DIR}/siesta2ph_utils.o
                                  ${SIESTA2PH_BUILD_DIR}/siesta2ph_linalg.o
                                  ${SIESTA2PH_BUILD_DIR}/siesta2ph_system.o
                                  ${SIESTA2PH_BUILD_DIR}/spglib_f08.o
                                  )
else()
  set(SIESTA2PH_BUILD_DIR ${SIESTA_BUILD_DIR}/Util/Contrib/siesta2ph CACHE PATH "Path to the siesta2ph build directory." FORCE)
  set(SIESTA2PH_MODULES_OBJ_FILES ${SIESTA2PH_BUILD_DIR}/CMakeFiles/siesta2ph_modules.dir/siesta2ph_io.f90.o
                                  ${SIESTA2PH_BUILD_DIR}/CMakeFiles/siesta2ph_modules.dir/siesta2ph_symmetry.f90.o
                                  ${SIESTA2PH_BUILD_DIR}/CMakeFiles/siesta2ph_modules.dir/siesta2ph_utils.f90.o
                                  ${SIESTA2PH_BUILD_DIR}/CMakeFiles/siesta2ph_modules.dir/siesta2ph_linalg.F90.o
                                  ${SIESTA2PH_BUILD_DIR}/CMakeFiles/siesta2ph_modules.dir/siesta2ph_system.f90.o
                                  ${SIESTA2PH_BUILD_DIR}/CMakeFiles/siesta2ph_modules.dir/spglib_f08.f90.o
                                  )
endif()
set(SIESTA2PH_MODULES_MOD_FILES ${SIESTA2PH_BUILD_DIR}/siesta2ph_io.mod
                                ${SIESTA2PH_BUILD_DIR}/siesta2ph_linalg.mod
                                ${SIESTA2PH_BUILD_DIR}/spglib_f08.mod
                                ${SIESTA2PH_BUILD_DIR}/siesta2ph_symmetry.mod
                                ${SIESTA2PH_BUILD_DIR}/siesta2ph_system.mod
                                ${SIESTA2PH_BUILD_DIR}/siesta2ph_utils.mod
                                )


# Add some flags
set(SIESTA2PH_FFLAGS "${CMAKE_Fortran_FLAGS}")
#
string(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
if(BT STREQUAL "DEBUG")
  set(SIESTA2PH_FFLAGS "${SIESTA2PH_FFLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
  set(SIESTA2PH_FPPFLAGS -DDEBUG)
endif()


if("${SIESTA_CMAKE}" STREQUAL "SIESTA_CMAKE-NOTFOUND")

  # Create Makefile for siesta2ph
  list(JOIN SPGLIB_LIBRARIES " " SPGLIB_LIBRARIES_STR)
  configure_file("${SIESTA2PH_SRC_DIR}/Makefile.in"
                 "${SIESTA2PH_BUILD_DIR}/Makefile" @ONLY)

  # Create the target siesta2ph_modules to build the siesta2ph modules
  add_custom_command(OUTPUT ${SIESTA2PH_BUILD_DIR}/libsiesta2ph_modules.a
    COMMAND make -C ${SIESTA2PH_BUILD_DIR} siesta2ph_modules
    COMMENT "Building siesta2ph_modules"
    DEPENDS ${SIESTA2PH_MODULES_SRC_FILES}
    BYPRODUCTS ${SIESTA2PH_MODULES_OBJ_FILES}
              ${SIESTA2PH_MODULES_MOD_FILES}
    )

  add_custom_target(siesta2ph_modules
    DEPENDS ${SIESTA2PH_BUILD_DIR}/libsiesta2ph_modules.a
    )

  # Create the targets for the siesta2ph executables
  add_custom_command(OUTPUT ${SIESTA2PH_BUILD_DIR}/siesta2ph.x
    COMMAND make -C ${SIESTA2PH_BUILD_DIR} siesta2ph.x
    COMMENT "Building siesta2ph.x"
    DEPENDS ${SIESTA2PH_BUILD_DIR}/libsiesta2ph_modules.a
            siesta2ph_modules
            ${SIESTA2PH_SRC_DIR}/siesta2ph.F90
    BYPRODUCTS ${SIESTA2PH_BUILD_DIR}/siesta2ph.o
    )

  add_custom_command(OUTPUT ${SIESTA2PH_BUILD_DIR}/siesta2fc.x
    COMMAND make -C ${SIESTA2PH_BUILD_DIR} siesta2fc.x
    COMMENT "Building siesta2fc.x"
    DEPENDS ${SIESTA2PH_BUILD_DIR}/libsiesta2ph_modules.a
            siesta2ph_modules
            ${SIESTA2PH_SRC_DIR}/siesta2fc.F90
    BYPRODUCTS ${SIESTA2PH_BUILD_DIR}/siesta2fc.o
    )

  add_custom_command(OUTPUT ${SIESTA2PH_BUILD_DIR}/siesta2dv.x
    COMMAND make -C ${SIESTA2PH_BUILD_DIR} siesta2dv.x
    COMMENT "Building siesta2dv.x"
    DEPENDS ${SIESTA2PH_BUILD_DIR}/libsiesta2ph_modules.a
            siesta2ph_modules
            ${SIESTA2PH_SRC_DIR}/siesta2dv.F90
    BYPRODUCTS ${SIESTA2PH_BUILD_DIR}/siesta2dv.o
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

  add_custom_target(siesta2ph DEPENDS siesta2ph.x siesta2fc.x siesta2dv.x)

else()

  # Create the cmake configuration file for siesta2ph
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/CMakeLists.txt.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/CMakeLists.txt" @ONLY)

  # Create the targets to copy and configure all the siesta2ph files to ${SIESTA_HOME}/Util/Contrib/siesta2ph
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp
    COMMAND mkdir -p "${SIESTA_HOME}/Util/Contrib/siesta2ph"
    COMMAND cp ${SIESTA2PH_SRC_FILES} "${SIESTA_HOME}/Util/Contrib/siesta2ph/"
    COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/CMakeLists.txt "${SIESTA_HOME}/Util/Contrib/siesta2ph/"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp"
    COMMENT "Copying siesta2ph to Util/Contrib"
    DEPENDS ${SIESTA2PH_SRC_FILES}
            ${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/CMakeLists.txt
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
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp
    COMMENT "Configuring SIESTA"
    )

  add_custom_target(patch-siesta-siesta2ph
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp
    )

  # Create the target siesta2ph_modules to build the siesta2ph modules
  add_custom_command(OUTPUT ${SIESTA2PH_BUILD_DIR}/libsiesta2ph_modules.a
    COMMAND cmake --build "${SIESTA_BUILD_DIR}" -t siesta2ph_modules
    COMMENT "Building siesta2ph_modules"
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp
            patch-siesta-siesta2ph
    BYPRODUCTS ${SIESTA2PH_MODULES_OBJ_FILES}
               ${SIESTA2PH_MODULES_MOD_FILES}
    )

  add_custom_target(siesta2ph_modules DEPENDS ${SIESTA2PH_BUILD_DIR}/libsiesta2ph_modules.a)

  # Create the targets for the siesta2ph executables
  add_custom_command(OUTPUT ${SIESTA2PH_BUILD_DIR}/siesta2ph.x
    COMMAND cmake --build "${SIESTA_BUILD_DIR}" -t siesta2ph_exe
    COMMENT "Building siesta2ph.x"
    DEPENDS ${SIESTA2PH_BUILD_DIR}/libsiesta2ph_modules.a
            siesta2ph_modules
    BYPRODUCTS ${SIESTA2PH_BUILD_DIR}//CMakeFiles/siesta2ph_exe.dir/siesta2ph.F90.o
    )

  add_custom_command(OUTPUT ${SIESTA2PH_BUILD_DIR}/siesta2fc.x
    COMMAND cmake --build "${SIESTA_BUILD_DIR}" -t siesta2fc_exe
    COMMENT "Building siesta2fc.x"
    DEPENDS ${SIESTA2PH_BUILD_DIR}/libsiesta2ph_modules.a
            siesta2ph_modules
    BYPRODUCTS ${SIESTA2PH_BUILD_DIR}//CMakeFiles/siesta2fc_exe.dir/siesta2fc.F90.o
    )

  add_custom_command(OUTPUT ${SIESTA2PH_BUILD_DIR}/siesta2dv.x
    COMMAND cmake --build "${SIESTA_BUILD_DIR}" -t siesta2dv_exe
    COMMENT "Building siesta2dv.x"
    DEPENDS ${SIESTA2PH_BUILD_DIR}/libsiesta2ph_modules.a
            siesta2ph_modules
    BYPRODUCTS ${SIESTA2PH_BUILD_DIR}//CMakeFiles/siesta2dv_exe.dir/siesta2dv.F90.o
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

  add_custom_target(siesta2ph DEPENDS siesta2ph.x siesta2fc.x siesta2dv.x)

  # Create a target to remove all siesta2ph files and restore siesta
  add_custom_target(remove-siesta2ph
    COMMENT "Removing siesta2ph from SIESTA"
    COMMAND rm -rf "${SIESTA_BUILD_DIR}/Util/Contrib/siesta2ph"
    COMMAND rm -rf "${SIESTA_HOME}/Util/Contrib/siesta2ph"
    COMMAND mv "${SIESTA_HOME}/Util/CMakeLists.txt" "${SIESTA_HOME}/Util/CMakeLists.txt_backup"
    COMMAND grep -v siesta2ph "${SIESTA_HOME}/Util/CMakeLists.txt_backup" > "${SIESTA_HOME}/Util/CMakeLists.txt"
    COMMAND rm -rf "${SIESTA_HOME}/Util/CMakeLists.txt_backup"
    COMMAND rm -rf "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp"
    VERBATIM
  )

  add_dependencies(restore-siesta remove-siesta2ph)

endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/utils/run-disps.bash.in"
               "${CMAKE_CURRENT_BINARY_DIR}/${SIESTA_VERSION}/run-disps.x" @ONLY)
