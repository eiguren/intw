
# Define qe2intw source directories and files
set(QE2INTW_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${QE_VERSION} CACHE PATH "Path to the pw2intw src directory." FORCE)
set(QE2INTW_SRC_FILES ${QE2INTW_SRC_DIR}/pw2intw.f90
                      ${QE2INTW_SRC_DIR}/create_intw_q_dirs.f90
                      )

# Define qe2intw build directory and object and module files
if("${QE_CMAKE}" STREQUAL "QE_CMAKE-NOTFOUND")
  set(QE2INTW_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${QE_VERSION} CACHE PATH "Path to the pw2intw build directory." FORCE)
else()
  set(QE2INTW_BUILD_DIR ${QE_BUILD_DIR}/bin CACHE PATH "Path to the pw2intw build directory." FORCE)
endif()


# Add some flags
set(QE2INTW_FFLAGS "${CMAKE_Fortran_FLAGS}")
#
string(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
if(BT STREQUAL "DEBUG")
  set(QE2INTW_FFLAGS "${SIESTA2INTW_FFLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
  set(QE2INTW_FPPFLAGS -DDEBUG)
endif()

# Add some debugging flags
string(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
if(BT STREQUAL "DEBUG")
    set(QE2INTW_FFLAGS ${CMAKE_Fortran_FLAGS_DEBUG})
    set(QE2INTW_FPPFLAGS -DDEBUG)
endif()


set(OUTPUT_PREFIX "...... ")

if("${QE_CMAKE}" STREQUAL "QE_CMAKE-NOTFOUND")

  # Create Makefile for qe2intw
  configure_file("${QE2INTW_SRC_DIR}/Makefile.in"
                 "${QE2INTW_BUILD_DIR}/Makefile" @ONLY)

  # Create the targets for the qe2intw executables
  add_custom_command(OUTPUT ${QE2INTW_BUILD_DIR}/pw2intw.x
    COMMAND make -C ${QE2INTW_BUILD_DIR} pw2intw.x | sed 's/^/${OUTPUT_PREFIX}/'
    COMMENT "Building pw2intw.x"
    DEPENDS ${QE2INTW_SRC_DIR}/pw2intw.f90
    BYPRODUCTS ${QE2INTW_BUILD_DIR}/pw2intw.o
    )

    add_custom_command(OUTPUT ${QE2INTW_BUILD_DIR}/create_intw_q_dirs.x
    COMMAND make -C ${QE2INTW_BUILD_DIR} create_intw_q_dirs.x | sed 's/^/${OUTPUT_PREFIX}/'
    COMMENT "Building create_intw_q_dirs.x"
    DEPENDS ${QE2INTW_SRC_DIR}/create_intw_q_dirs.f90
    BYPRODUCTS ${QE2INTW_BUILD_DIR}/create_intw_q_dirs.o
    )

else()

  # Create the cmake configuration file for qe2intw
  configure_file("${QE2INTW_SRC_DIR}/CMakeLists.txt_intw.in"
	         "${QE2INTW_SRC_DIR}/CMakeLists.txt_intw" @ONLY)

  # Create the targets to copy and configure all the qe2intw files to ${QE_HOME}/PP
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp
    COMMAND cp ${QE2INTW_SRC_FILES} "${QE_HOME}/PP/src/"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp"
    COMMENT "Copying qe2intw to QE/PP"
    DEPENDS ${QE2INTW_SRC_FILES}
    )

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake_backup.stamp
    COMMAND bash -c "if [ ! -f ${QE_HOME}/PP/CMakeLists.txt_qe ]; then mv ${QE_HOME}/PP/CMakeLists.txt ${QE_HOME}/PP/CMakeLists.txt_qe; fi"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake_backup.stamp"
    COMMENT "Creating a backup of the original QE/PP CMakeLists.txt"
    VERBATIM
    )
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp
    COMMAND cp "${QE2INTW_SRC_DIR}/CMakeLists.txt_intw" "${QE_HOME}/PP/CMakeLists.txt"
    COMMAND touch "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp"
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake_backup.stamp
            ${QE2INTW_SRC_DIR}/CMakeLists.txt_intw
    COMMENT "Copying CMakeLists.txt to QE/PP"
    )

  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp
    COMMAND cmake -S "${QE_HOME}" -B "${QE_BUILD_DIR}" | sed 's/^/${OUTPUT_PREFIX}/'
    COMMAND touch ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp
	    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp
    COMMENT "Configuring QE"
    )

  add_custom_target(patch-qe
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake_backup.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp
    )

  # Create the targets for the qe2intw executables
  add_custom_command(OUTPUT ${QE2INTW_BUILD_DIR}/pw2intw.x
    COMMAND cmake --build "${QE_BUILD_DIR}" -t qe_pp_pw2intw_exe | sed 's/^/${OUTPUT_PREFIX}/'
    COMMENT "Building pw2intw.x"
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake_backup.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp
	    patch-qe
    BYPRODUCTS ${QE_BUILD_DIR}/PP/mod/qe_pp_pw2intw_exe/*
               ${QE_BUILD_DIR}/PP/CMakeFiles/qe_pp_pw2intw_exe.dir/src/pw2intw.f90.o
    )

  add_custom_command(OUTPUT ${QE2INTW_BUILD_DIR}/create_intw_q_dirs.x
    COMMAND  cmake --build "${QE_BUILD_DIR}" -t qe_pp_create_intw_q_dirs_exe | sed 's/^/${OUTPUT_PREFIX}/'
    COMMENT "Building create_intw_q_dirs.x"
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake_backup.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp
            ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp
            patch-qe
    BYPRODUCTS ${QE_BUILD_DIR}/PP/mod/qe_pp_create_intw_q_dirs_exe/*
               ${QE_BUILD_DIR}/PP/CMakeFiles/qe_pp_create_intw_q_dirs_exe.dir/src/create_intw_q_dirs.f90.o
    )

  # Create a target to remove all qe2intw files and restore QE
  add_custom_target(restore-qe
    COMMAND rm -rf "${QE2INTW_BUILD_DIR}/pw2intw.x" "${QE_BUILD_DIR}/PP/mod/qe_pp_pw2intw_exe" "${QE_BUILD_DIR}/PP/CMakeFiles/qe_pp_pw2intw_exe.dir"
    COMMAND rm -rf "${QE_HOME}/PP/src/pw2intw.f90"
    COMMAND rm -rf "${QE2INTW_BUILD_DIR}/create_intw_q_dirs.x" "${QE_BUILD_DIR}/PP/mod/qe_pp_create_intw_q_dirs_exe" "${QE_BUILD_DIR}/PP/CMakeFiles/qe_pp_create_intw_q_dirs_exe.dir"
    COMMAND rm -rf "${QE_HOME}/PP/src/create_intw_q_dirs.f90"
    COMMAND bash -c "if [ -f ${QE_HOME}/PP/CMakeLists.txt_qe ]; then cp ${QE_HOME}/PP/CMakeLists.txt_qe ${QE_HOME}/PP/CMakeLists.txt; fi"
    COMMAND bash -c "if [ -f ${QE_HOME}/PP/CMakeLists.txt_qe ]; then cmake -S ${QE_HOME} -B ${QE_BUILD_DIR}; fi"
    COMMAND rm -rf "${QE_HOME}/PP/CMakeLists.txt_qe"
    COMMAND rm -rf "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.stamp"
    COMMAND rm -rf "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/configure.stamp"
    COMMAND rm -rf "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake_backup.stamp"
    COMMAND rm -rf "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cmake.stamp"
    COMMENT "Restoring QE"
    VERBATIM
    )

endif()

# Create the target to build qe2intw
add_custom_target(pw2intw.x DEPENDS ${QE2INTW_BUILD_DIR}/pw2intw.x)
add_custom_target(create_intw_q_dirs.x DEPENDS ${QE2INTW_BUILD_DIR}/create_intw_q_dirs.x)

add_custom_target(qe2intw ALL DEPENDS pw2intw.x create_intw_q_dirs.x)
