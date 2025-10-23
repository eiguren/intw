
cmake_minimum_required(VERSION 3.12)

# Add our local modlues to the module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

include(PreventInSourceBuilds)

project(intw
  VERSION 0.0.0
  LANGUAGES Fortran C
)


include(GNUInstallDirs)

include(SetCompilers)


# Define user options
option(USE_OPENMP "Use OpenMP for parallelization" OFF)
option(USE_MPI_ON_TESTS "Use the MPI parallelization to execute the tests " OFF)


# Find version of INTW with git describe
execute_process(COMMAND git describe --tags --match "v[0-9]*.[0-9]*.[0-9]*" --dirty OUTPUT_VARIABLE GIT_VERSION RESULTS_VARIABLE GIT_RESULT OUTPUT_STRIP_TRAILING_WHITESPACE)
if(GIT_RESULT EQUAL 0)
  # git describe command succeeded, use the retrieved version number
  set(VERSION_STRING "${GIT_VERSION}")
else()
  # git describe command failed, use version number specified in project()
  set(VERSION_STRING "v${PROJECT_VERSION}")
endif()


set(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> cr <TARGET> <LINK_FLAGS> <OBJECTS>") # This are the flags used in OpenBLAS and netlib lapack
# cr should be correct, c for crate, and r for replace (cmake uses q for quick-append instead of r).


# Find packages

# Find OpenMP flag
if(USE_OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# Find BLAS and LAPACK
find_package(LAPACK REQUIRED)

# Find triangle and tetgen
include(SetTriFS)


# Quantum Espresso
option(USE_QE "Use Quantum Espresso" OFF)

# Try to define the QE home directory from the environmental variable
if(DEFINED ENV{QE_HOME})
  set(QE_HOME $ENV{QE_HOME} CACHE PATH "Path to the Quantum Espresso home directory." FORCE)
endif()

# Turn on USE_QE if the QE home directory is defined
if(DEFINED QE_HOME)
  set(USE_QE ON CACHE BOOL "Use Quantum Espresso" FORCE)
endif()

if (USE_QE)
  #
  # Check if the QE home directory is defined
  if(NOT DEFINED QE_HOME)
    message(FATAL_ERROR "QE_HOME variable not set.")
  endif()
  #
  # Setup QE variables
  include(SetQuantumEspresso)
  #
endif()


# SIESTA
option(USE_SIESTA "Use SIESTA" OFF)

# Try to define the SIESTA home directory from the environmental variable
if(DEFINED ENV{SIESTA_HOME})
  set(SIESTA_HOME $ENV{SIESTA_HOME} CACHE PATH "Path to the SIESTA home directory." FORCE)
endif()

# Turn on USE_SIESTA if the SIESTA home directory is defined
if(DEFINED SIESTA_HOME)
  set(USE_SIESTA ON CACHE BOOL "Use SIESTA" FORCE)
endif()

if (USE_SIESTA)
  #
  # Check if the SIESTA home directory is defined
  if(NOT DEFINED SIESTA_HOME)
    message(FATAL_ERROR "SIESTA_HOME variable not set.")
  endif()
  #
  # Add Spglib library
  include(FetchSpglib)
  #
  # Setup SIESTA variables
  include(SetSiesta)
  #
endif()


# Wannier90
option(USE_WANNIER "Use Wannier90" OFF)

# Try to define the Wannier home directory from the environmental variable
if(DEFINED ENV{W_HOME})
  set(W_HOME $ENV{W_HOME} CACHE PATH "Path to the Wannier90 home directory." FORCE)
endif()

# Turn on USE_WANNIER if the Wannier home directory is defined
if(DEFINED W_HOME)
  set(USE_WANNIER ON CACHE BOOL "Use Wannier90" FORCE)
endif()

if (USE_WANNIER)
  #
  # Check if the Wannier home directory is defined
  if(NOT DEFINED W_HOME)
    message(FATAL_ERROR "W_HOME variable not set.")
  endif()
  #
  # Setup Wannier variables
  include(SetWannier)
  #
endif()


# Add source
add_subdirectory(src/modules)
add_subdirectory(src/utilities)
if (USE_QE)
  add_subdirectory(src/QEtoINTW)
endif()
if (USE_SIESTA)
  add_subdirectory(src/siesta2intw)
  add_subdirectory(src/siesta2ph)
endif()


# Add tests
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  enable_testing()
  add_subdirectory(test)
else()
  set(DISABLED_TESTS "all")
endif()


# Add uninstall target
if(NOT TARGET uninstall)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()


# Print a summary of all the variables

message("")
message("SUMMARY:")

message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_Fortran_COMPILER: ${CMAKE_Fortran_COMPILER}")
message(STATUS "Fortran compiler vendor: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "Fortran compiler version: ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "CMAKE_AR: ${CMAKE_AR}")
message(STATUS "CMAKE_RANLIB: ${CMAKE_RANLIB}")
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
if (USE_OPENMP)
  if(OpenMP_C_FOUND)
    if(OpenMP_C_FLAGS)
      message(STATUS "OpenMP_C_FLAGS: ${OpenMP_C_FLAGS}")
    else()
      message(STATUS "OpenMP_C_INCLUDE_DIRS: ${OpenMP_C_INCLUDE_DIRS}")
      message(STATUS "OpenMP_C_LIBRARIES: ${OpenMP_C_LIBRARIES}")
    endif()
  endif()
  if(OpenMP_Fortran_FOUND)
    if(OpenMP_Fortran_FLAGS)
    message(STATUS "OpenMP_Fortran_FLAGS: ${OpenMP_Fortran_FLAGS}")
  else()
    message(STATUS "OpenMP_Fortran_INCLUDE_DIRS: ${OpenMP_Fortran_INCLUDE_DIRS}")
    message(STATUS "OpenMP_Fortran_LIBRARIES: ${OpenMP_Fortran_LIBRARIES}")
  endif()
  endif()
endif()
message(STATUS "BLAS_LINKER_FLAGS: ${BLAS_LINKER_FLAGS}")
message(STATUS "BLAS_LIBRARIES: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK_LINKER_FLAGS: ${LAPACK_LINKER_FLAGS}")
message(STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_Fortran_FLAGS: ${CMAKE_Fortran_FLAGS}")
message(STATUS "CMAKE_Fortran_FLAGS_RELEASE: ${CMAKE_Fortran_FLAGS_RELEASE}")
message(STATUS "CMAKE_Fortran_FLAGS_DEBUG: ${CMAKE_Fortran_FLAGS_DEBUG}")
message(STATUS "CMAKE_Fortran_FLAGS_RELWITHDEBINFO: ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
message(STATUS "CMAKE_Fortran_FLAGS_MINSIZEREL: ${CMAKE_Fortran_FLAGS_MINSIZEREL}")
if(Python3_Interpreter_FOUND)
  message(STATUS "Python3_EXECUTABLE: ${Python3_EXECUTABLE}")
else()
  message(STATUS "Python3: NOT-FOUND (test suite disabled)")
endif()
if (USE_QE)
  message(STATUS "QE_HOME: ${QE_HOME}")
  message(STATUS "QE_VERSION: ${QE_VERSION}")
endif()
if (USE_SIESTA)
  message(STATUS "SIESTA_HOME: ${SIESTA_HOME}")
  message(STATUS "SIESTA_VERSION: ${SIESTA_VERSION}")
  message(STATUS "SPGLIB_LIBRARIES: ${SPGLIB_LIBRARIES}")
endif()
if (USE_WANNIER)
  message(STATUS "W_HOME: ${W_HOME}")
  message(STATUS "W_VERSION: ${W_VERSION}")
endif()
message(STATUS "TETGEN_EXE: ${TETGEN_EXE}")
message(STATUS "TRIANGLE_EXE: ${TRIANGLE_EXE}")
if(NOT TRIFS)
  message(STATUS "Disabled utilities: triFS.x")
endif()
if(DISABLED_TESTS)
  message(STATUS "Disabled tests: ${DISABLED_TESTS}")
endif()
message("")
