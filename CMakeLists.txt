
cmake_minimum_required(VERSION 3.12)

# Add our local modlues to the module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

#set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
include(PreventInSourceBuilds)
#include(PreventInBuildInstalls)
#if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
#  message(FATAL_ERROR "In-source builds are not allowed.")
#endif()


#set(CMAKE_RULE_MESSAGES OFF)
#set(CMAKE_VERBOSE_MAKEFILE ON)

project(intw
  VERSION 0.1
  LANGUAGES Fortran C
)

include(GNUInstallDirs)

if(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
  message(FATAL_ERROR "Fortran compiler does not support F90")
endif(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  if(DEFINED ENV{CPATH})
    string(REPLACE ":" ";" CPATH $ENV{CPATH})
  endif()
  include_directories(${CPATH})
endif()

option(USE_MPI "Use the MPI library for parallelization" OFF)
option(USE_OPENMP "Use OpenMP for parallelization" OFF)
option(USE_MPI_ON_TESTS "Use the MPI parallelization to execute the tests " OFF)

#message("CMAKE_Fortran_CREATE_SHARED_LIBRARY: ${CMAKE_Fortran_CREATE_SHARED_LIBRARY}")
#message("CMAKE_Fortran_CREATE_SHARED_MODULE: ${CMAKE_Fortran_CREATE_SHARED_MODULE}")
#message("CMAKE_Fortran_CREATE_STATIC_LIBRARY: ${CMAKE_Fortran_CREATE_STATIC_LIBRARY}")
#message("CMAKE_Fortran_COMPILE_OBJECT: ${CMAKE_Fortran_COMPILE_OBJECT}")
#message("CMAKE_Fortran_LINK_EXECUTABLE: ${CMAKE_Fortran_COMPILE_OBJECT}")
#message("CMAKE_Fortran_ARCHIVE_CREATE: ${CMAKE_Fortran_ARCHIVE_CREATE}")
#message("CMAKE_Fortran_ARCHIVE_APPEND: ${CMAKE_Fortran_ARCHIVE_APPEND}")
#message("CMAKE_Fortran_ARCHIVE_FINISH: ${CMAKE_Fortran_ARCHIVE_FINISH}")


#set(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> ruv <TARGET> <LINK_FLAGS> <OBJECTS>") # This are the flags used in QE and wannier
#set(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> cru <TARGET> <LINK_FLAGS> <OBJECTS>") # This are the flags used in fftw3 and siesta
set(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> cr <TARGET> <LINK_FLAGS> <OBJECTS>") # This are the flags used in OpenBLAS and netlib lapack
# cr should be correct, c for crate, and r for replace (cmake uses q for quick-append instead of r).
# v adds verbose output.
# And I don't know how works u, and I'm getting a warning also...

include(SetFortranFlags)

# Find MPI if needed
if(USE_MPI)
  find_package(MPI REQUIRED)
  add_definitions(-DMPI)
  # Instead of ussing target_include_directories (which it is prefered), I am ussing include_directories, because all the targets will depend on MPI,
  # and include_directories will add the include directories to all new targets.
  include_directories(${MPI_Fortran_INCLUDE_PATH})
  add_compile_options(${MPI_Fortran_COMPILE_FLAGS})
endif()

# Find OpenMP flag
if(USE_OPENMP)
  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
                   Fortran REQUIRED "-qopenmp"        # Intel
                                    "-fopenmp"        # GNU
                  )
endif()




# Find BLAS and LAPACK
find_package(LAPACK REQUIRED)




# Quantum Espresso
option(USE_QE "Use Quantum Espresso" OFF)
if(DEFINED ENV{QE_HOME})
  set(QE_HOME $ENV{QE_HOME} CACHE PATH "Path to the Quantum Espresso home directory." FORCE)
  set(USE_QE ON CACHE BOOL "Use Quantum Espresso" FORCE)
endif()
if (USE_QE)
  include(SetQuantumEspresso)
endif()

# Wannier90
option(USE_WANNIER "Use Wannier90" OFF)
if(DEFINED ENV{W_HOME})
  set(W_HOME $ENV{W_HOME} CACHE PATH "Path to the Wannier90 home directory." FORCE)
  set(USE_WANNIER ON CACHE BOOL "Use Wannier90" FORCE)
endif()
if (USE_WANNIER)
  include(SetWannier)
endif()


add_subdirectory(src/modules)
add_subdirectory(src/utilities)
if (USE_QE)
  add_subdirectory(src/QEtoINTW)
endif()


# Find python3 for test scripts
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

message("-- C compiler: ${CMAKE_C_COMPILER}")
message("-- Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message("-- Fortran compiler vendor: ${CMAKE_Fortran_COMPILER_ID}")
message("-- Archive: ${CMAKE_AR}")
message("-- Ranlib: ${CMAKE_RANLIB}")
message("-- Install prefix: ${CMAKE_INSTALL_PREFIX} (Use -DCMAKE_INSTALL_PREFIX=/install/prefix to change the default value)")
message("-- MPI C compiler: ${MPI_C_COMPILER}")
message("-- MPI Fortran compiler: ${MPI_Fortran_COMPILER}")
message("-- MPI include dirs: ${MPI_Fortran_INCLUDE_DIRS}")
message("-- MPI link flags: ${MPI_Fortran_LINK_FLAGS}")
message("-- MPI library: ${MPI_Fortran_LIBRARIES}")
message("-- BLAS flags: ${BLAS_LINKER_FLAGS}")
message("-- BLAS library: ${BLAS_LIBRARIES}")
message("-- LAPACK flags: ${LAPACK_LINKER_FLAGS}")
message("-- LAPACK library: ${LAPACK_LIBRARIES}")
message("-- CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message("-- FFLAGS: ${CMAKE_Fortran_FLAGS}")
message("-- FFLAGS_RELEASE: ${CMAKE_Fortran_FLAGS_RELEASE}")
message("-- FFLAGS_DEBUG: ${CMAKE_Fortran_FLAGS_DEBUG}")
message("-- FFLAGS_RELWITHDEBINFO: ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
# message("-- FFLAGS_MINSIZEREL: ${CMAKE_Fortran_FLAGS_MINSIZEREL}")
if(Python3_Interpreter_FOUND)
  message("-- Python3: ${Python3_EXECUTABLE}")
else()
  message("-- Python3: NOT-FOUND (test suite disabled)")
endif()
if (USE_QE)
  message("-- QE_HOME: ${QE_HOME}")
  message("-- QE_VERSION: ${QE_VERSION}")
endif()
if (USE_WANNIER)
  message("-- W_HOME: ${W_HOME}")
  message("-- W_VERSION: ${W_VERSION}")
endif()
if(DISABLED_TESTS)
  message("-- Disabled tests: ${DISABLED_TESTS}")
endif()
message("")
