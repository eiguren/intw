

if("${QE_CMAKE}" STREQUAL "QE_CMAKE-NOTFOUND")
  # If QE is built with autotools
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/pw2intw.f90.stamp
    COMMAND cp "${CMAKE_CURRENT_SOURCE_DIR}/${QE_VERSION}/pw2intw.f90" "${QE_HOME}/PP/src/pw2intw.f90"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/pw2intw.f90.stamp"
    COMMENT "Copying pw2intw.f90 to QE/PP"
    )

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Makefile_intw.stamp
    COMMAND cp "${CMAKE_CURRENT_SOURCE_DIR}/${QE_VERSION}/Makefile_intw" "${QE_HOME}/PP/src/Makefile_intw"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Makefile_intw.stamp"
    COMMENT "Copying Makefile_intw to QE/PP"
    )

  add_custom_target(pw2intw.x ALL make -C ${QE_HOME}/PP/src -f Makefile_intw pw2intw.x
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/pw2intw.f90.stamp ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Makefile_intw.stamp
    BYPRODUCTS ${QE_HOME}/PP/src/pw2intw.o ${QE_HOME}/PP/src/pw2intw.x
    )

  add_custom_target(restore-qe
    COMMAND rm -rf "${QE_HOME}/PP/src/pw2intw.o" "${QE_HOME}/PP/src/pw2intw.x"
    COMMAND rm -rf "${QE_HOME}/PP/src/pw2intw.f90" "${QE_HOME}/PP/src/Makefile_intw"
    COMMAND rm -rf "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/pw2intw.f90.stamp" "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/Makefile_intw.stamp"
    COMMENT "QE restored"
    )
else()
  # If QE is built with CMake
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/pw2intw.f90.stamp
    COMMAND cp "${CMAKE_CURRENT_SOURCE_DIR}/${QE_VERSION}/pw2intw.f90" "${QE_HOME}/PP/src/pw2intw.f90"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/pw2intw.f90.stamp"
    COMMENT "Copying pw2intw.f90 to QE/PP"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${QE_VERSION}/pw2intw.f90
    )

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/QE_CMakeLists.txt_backup.stamp
    COMMAND bash -c "if [ ! -f ${QE_HOME}/PP/CMakeLists.txt_qe ]; then mv ${QE_HOME}/PP/CMakeLists.txt ${QE_HOME}/PP/CMakeLists.txt_qe; fi"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/QE_CMakeLists.txt_backup.stamp"
    COMMENT "Creating a backup of the original QE/PP CMakeLists.txt"
    VERBATIM
    )
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/QE_CMakeLists.txt.stamp
    COMMAND cp "${CMAKE_CURRENT_SOURCE_DIR}/${QE_VERSION}/CMakeLists.txt_intw" "${QE_HOME}/PP/CMakeLists.txt"
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

