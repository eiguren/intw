#############################################################################
#
#############################################################################

# List of supported SIESTA versions
set(SUPPORTED_SIESTA
  "4.1.5"
  "5.0.0"
  )

if (NOT SIESTA_VERSION)

  message("-- Checking SIESTA.")


  # Find SIESTA build dir
  file(GLOB files_and_dirs RELATIVE ${SIESTA_HOME} ${SIESTA_HOME}/*)
  set(dirs "")
  foreach(file_or_dir ${files_and_dirs})
    if(IS_DIRECTORY ${SIESTA_HOME}/${file_or_dir})
      list(APPEND dirs ${file_or_dir})
    endif()
  endforeach()

  find_file(SIESTA_CMAKE CMakeCache.txt PATHS ${SIESTA_HOME} PATH_SUFFIXES ${dirs} NO_DEFAULT_PATH)
  set(SIESTA_CMAKE "${SIESTA_CMAKE}" CACHE INTERNAL "SIESTA cmake cache.")
  if("${SIESTA_CMAKE}" STREQUAL "SIESTA_CMAKE-NOTFOUND")
    message("-- SIESTA built with autotools.")
    find_program(_siesta_exe siesta PATHS ${SIESTA_HOME} PATH_SUFFIXES ${dirs} NO_DEFAULT_PATH)
    get_filename_component(SIESTA_BUILD_DIR ${_siesta_exe} DIRECTORY)
  else()
    message("-- SIESTA built with CMake.")
    get_filename_component(SIESTA_BUILD_DIR ${SIESTA_CMAKE} DIRECTORY)
  endif()
  set(SIESTA_BUILD_DIR "${SIESTA_BUILD_DIR}" CACHE PATH "SIESTA build directory.")
  message("-- Found SIESTA_BUILD_DIR: ${SIESTA_BUILD_DIR}")


  # Some checks
  message("-- Checking SIESTA build.")
  find_program(SIESTA_EXE siesta PATHS ${SIESTA_BUILD_DIR} PATH_SUFFIXES ${dirs} NO_DEFAULT_PATH)
  set(SIESTA_EXE "${SIESTA_EXE}" CACHE PATH "SIESTA executable.")
  if("${SIESTA_EXE}" STREQUAL "SIESTA_EXE-NOTFOUND")
    message(FATAL_ERROR "SIESTA_EXE NOT FOUND! SIESTA is not compiled or SIESTA_HOME is incorrect.")
  endif()




  # Find SIESTA version number
  find_file(SIESTA_VERSION_FILE NAMES version.info SIESTA.release SIESTA.version PATHS ${SIESTA_HOME} NO_DEFAULT_PATH)
  set(SIESTA_VERSION_FILE "${SIESTA_VERSION_FILE}" CACHE INTERNAL "SIESTA version file.")
  if("${SIESTA_VERSION_FILE}" STREQUAL "SIESTA_VERSION_FILE-NOTFOUND")
    message(FATAL_ERROR "SIESTA_VERSION_FILE NOT FOUND!")
  endif()

  file(STRINGS ${SIESTA_VERSION_FILE} SIESTA_VERSION)
  set(SIESTA_VERSION "${SIESTA_VERSION}" CACHE INTERNAL "SIESTA version.")
  message("-- Found SIESTA_VERSION: ${SIESTA_VERSION}")

endif()


# Check if the given SIESTA version is supported by intw
list(FIND SUPPORTED_SIESTA ${SIESTA_VERSION} SIESTA_SUPPORTED)
if("${SIESTA_SUPPORTED}" EQUAL "-1")
  message("SIESTA especified by SIESTA_HOME is not suported!")
  message("Supported SIESTA versions: ${SUPPORTED_SIESTA}")
  message(FATAL_ERROR)
endif()
