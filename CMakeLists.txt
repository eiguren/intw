
cmake_minimum_required(VERSION 3.10)

# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

#set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
include(PreventInSourceBuilds)
#include(PreventInBuildInstalls)
#if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
#	message(FATAL_ERROR "In-source builds are not allowed.")
#endif()


#set(CMAKE_RULE_MESSAGES OFF)
#set(CMAKE_VERBOSE_MAKEFILE ON)

project(intw
	VERSION 0.1
	LANGUAGES Fortran C
	)

include(GNUInstallDirs)

IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

OPTION(USE_MPI "Use the MPI library for parallelization" ON)
OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)

#message("CMAKE_Fortran_CREATE_SHARED_LIBRARY: ${CMAKE_Fortran_CREATE_SHARED_LIBRARY}")
#message("CMAKE_Fortran_CREATE_SHARED_MODULE: ${CMAKE_Fortran_CREATE_SHARED_MODULE}")
#message("CMAKE_Fortran_CREATE_STATIC_LIBRARY: ${CMAKE_Fortran_CREATE_STATIC_LIBRARY}")
#message("CMAKE_Fortran_COMPILE_OBJECT: ${CMAKE_Fortran_COMPILE_OBJECT}")
#message("CMAKE_Fortran_LINK_EXECUTABLE: ${CMAKE_Fortran_COMPILE_OBJECT}")
#message("CMAKE_Fortran_ARCHIVE_CREATE: ${CMAKE_Fortran_ARCHIVE_CREATE}")
#message("CMAKE_Fortran_ARCHIVE_APPEND: ${CMAKE_Fortran_ARCHIVE_APPEND}")
#message("CMAKE_Fortran_ARCHIVE_FINISH: ${CMAKE_Fortran_ARCHIVE_FINISH}")


#set(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> ruv <TARGET> <LINK_FLAGS> <OBJECTS>") # This are the flags used in QE and wannier
set(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> cru <TARGET> <LINK_FLAGS> <OBJECTS>") # This are the flags used in fftw3 and siesta
#set(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> cr <TARGET> <LINK_FLAGS> <OBJECTS>") # This are the flags used in OpenBLAS and netlib lapack
# cr should be correct, c for crate, and r for replace (cmake uses q for quick-append instead of r).
# v adds verbose output.
# And I don't know how works u, and I'm getting a warning also...

INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)

# Find MPI if needed
if(USE_MPI)
	find_package(MPI REQUIRED)
	add_definitions(-DMPI)
	# Instead of ussing target_include_directories (which it is prefered), I am ussing include_directories, because all the targets will depend on MPI,
	# and include_directories will add the include directories to all new targets.
	include_directories(${MPI_Fortran_INCLUDE_PATH})
	add_compile_options(${MPI_Fortran_COMPILE_FLAGS})
endif()





# Find BLAS and LAPACK
find_package(LAPACK REQUIRED)

# FFTW3 (Asier)
#find_package(fftw3 REQUIRED)



# Quantum Espresso
#set(QE_HOME "/home/haritz/Codes/quantum-espresso/espresso-4.3.2-intel")
#set(QE_HOME "/home/asier/KODEAK/espresso-4.3.2")
set(QE-LAST_HOME "/home/asier/KODEAK/q-e/")
#set(QE-LAST_HOME "/home/asier/KODEAK/espresso-4.3.2")

# Wannier90
#set(W_HOME "/home/haritz/Codes/wannier90/wannier90-1.2-intel")
set(W_HOME "/home/asier/KODEAK/wannier90-3.0.0/")


add_subdirectory(src/modules)
add_subdirectory(src/utilities)
add_subdirectory(src/QEtoINTW)

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
message("-- FFLAGS_TESTING: ${CMAKE_Fortran_FLAGS_TESTING}")
message("-- FFLAGS_RELWITHDEBINFO: ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
message("-- FFLAGS_MINSIZEREL: ${CMAKE_Fortran_FLAGS_MINSIZEREL}")

message("")
