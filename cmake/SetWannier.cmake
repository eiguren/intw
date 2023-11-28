#############################################################################
#
#############################################################################

# List of supported Wannier90 versions
set(SUPPORTED_W
  "3.1.0"
  )

if (NOT W_VERSION)

  message("-- Checking Wannier90.")


  # Find Wannier90 build dir
  set(W_BUILD_DIR "${W_HOME}" CACHE PATH "Wannier90 build directory.")


  # Some checks
  message("-- Checking Wannier90 build.")
  find_program(W_EXE wannier90.x PATHS ${W_BUILD_DIR} NO_DEFAULT_PATH)
  set(W_EXE "${W_EXE}" CACHE PATH "Wannier90 executable.")
  if("${W_EXE}" STREQUAL "W_EXE-NOTFOUND")
  message(FATAL_ERROR "W_EXE NOT FOUND! Wannier90 is not compiled or WANNIER_HOME is incorrect.")
  endif()


  # Find Wannier90 version number
  find_file(W_VERSION_FILE io.F90 PATHS ${W_HOME}/src NO_DEFAULT_PATH)
  set(W_VERSION_FILE "${W_VERSION_FILE}" CACHE INTERNAL "Wannier90 version file.")
  if("${W_VERSION_FILE}" STREQUAL "W_VERSION_FILE-NOTFOUND")
     message(FATAL_ERROR "W_VERSION_FILE NOT FOUND!")
  endif()

  file(STRINGS ${W_VERSION_FILE} W_VERSION_LINE REGEX "w90_version *=")
  string(REGEX MATCH "' *([0-9]+)(\\.[0-9]+)?(\\.[0-9]+)? *'" W_VERSION ${W_VERSION_LINE})
  string(REPLACE "'" "" W_VERSION "${W_VERSION}")
  string(REPLACE " " "" W_VERSION "${W_VERSION}")
  set(W_VERSION "${W_VERSION}" CACHE INTERNAL "Wannier90 version.")
  message("-- Found W_VERSION: ${W_VERSION}")

endif()