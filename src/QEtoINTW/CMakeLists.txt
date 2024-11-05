
set(QE2INTW_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${QE_VERSION} CACHE PATH "Path to the pw2intw src directory." FORCE)


if("${QE_CMAKE}" STREQUAL "QE_CMAKE-NOTFOUND")

  set(QE2INTW_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${QE_VERSION} CACHE PATH "Path to the pw2intw build directory." FORCE)

  # Add some debugging flags
  string(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
  if(BT STREQUAL "DEBUG")
      set(QE2INTW_FFLAGS ${CMAKE_Fortran_FLAGS_DEBUG})
      set(QE2INTW_DFLAGS -DDEBUG)
  endif()

  configure_file("${QE2INTW_SRC_DIR}/Makefile.in"
                "${QE2INTW_BUILD_DIR}/Makefile" @ONLY)


  add_custom_command(OUTPUT ${QE2INTW_BUILD_DIR}/pw2intw.x
    COMMAND make -C ${QE2INTW_BUILD_DIR} pw2intw.x
    COMMENT "Building pw2intw"
    DEPENDS ${QE2INTW_SRC_DIR}/pw2intw.f90
    BYPRODUCTS ${QE2INTW_BUILD_DIR}/pw2intw.o
    )

  add_custom_target(pw2intw.x ALL
    DEPENDS ${QE2INTW_BUILD_DIR}/pw2intw.x
    )

else()

  set(QE2INTW_BUILD_DIR ${QE_BUILD_DIR}/bin CACHE PATH "Path to the pw2intw build directory." FORCE)

  # If QE is built with CMake
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/pw2intw.f90.stamp
    COMMAND cp "${QE2INTW_SRC_DIR}/pw2intw.f90" "${QE_HOME}/PP/src/pw2intw.f90"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/pw2intw.f90.stamp"
    COMMENT "Copying pw2intw.f90 to QE/PP"
    DEPENDS ${QE2INTW_SRC_DIR}/pw2intw.f90
    )

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/QE_CMakeLists.txt_backup.stamp
    COMMAND bash -c "if [ ! -f ${QE_HOME}/PP/CMakeLists.txt_qe ]; then mv ${QE_HOME}/PP/CMakeLists.txt ${QE_HOME}/PP/CMakeLists.txt_qe; fi"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/QE_CMakeLists.txt_backup.stamp"
    COMMENT "Creating a backup of the original QE/PP CMakeLists.txt"
    VERBATIM
    )
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/QE_CMakeLists.txt.stamp
    COMMAND cp "${QE2INTW_SRC_DIR}/CMakeLists.txt_intw" "${QE_HOME}/PP/CMakeLists.txt"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/QE_CMakeLists.txt.stamp"
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/QE_CMakeLists.txt_backup.stamp
    COMMENT "Copying CMakeLists.txt to QE/PP"
    )

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/QE_configure.stamp
    COMMAND cmake -S "${QE_HOME}" -B "${QE_BUILD_DIR}"
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/QE_configure.stamp
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/pw2intw.f90.stamp ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/QE_CMakeLists.txt.stamp
    COMMENT "Configuring QE"
    )

  add_custom_target(pw2intw.x ALL cmake --build "${QE_BUILD_DIR}" -t qe_pp_pw2intw_exe
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/QE_configure.stamp
    BYPRODUCTS ${QE_BUILD_DIR}/PP/mod/qe_pp_pw2intw_exe ${QE_BUILD_DIR}/PP/CMakeFiles/qe_pp_pw2intw_exe.dir ${QE_BUILD_DIR}/bin/pw2intw.x
    )

  add_custom_target(restore-qe
    COMMAND rm -rf "${QE_BUILD_DIR}/bin/pw2intw.x" "${QE_BUILD_DIR}/PP/mod/qe_pp_pw2intw_exe" "${QE_BUILD_DIR}/PP/CMakeFiles/qe_pp_pw2intw_exe.dir"
    COMMAND rm -rf "${QE_HOME}/PP/src/pw2intw.f90"
    COMMAND bash -c "if [ -f ${QE_HOME}/PP/CMakeLists.txt_qe ]; then cp ${QE_HOME}/PP/CMakeLists.txt_qe ${QE_HOME}/PP/CMakeLists.txt; fi"
    COMMAND bash -c "if [ -f ${QE_HOME}/PP/CMakeLists.txt_qe ]; then cmake -S ${QE_HOME} -B ${QE_BUILD_DIR}; fi"
    COMMAND rm -rf "${QE_HOME}/PP/CMakeLists.txt_qe"
    COMMAND rm -rf "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/QE_configure.stamp" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/QE_CMakeLists.txt_backup.stamp" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/QE_CMakeLists.txt.stamp" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/pw2intw.f90.stamp"
    COMMENT "QE restored"
    VERBATIM
    )

endif()
