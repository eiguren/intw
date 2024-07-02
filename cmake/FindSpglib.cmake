#[=======================================================================[.rst
FindSpglib
--------

Find Spglib library

This module finds an installed `Spglib library`_.
.

``C`` language must be enabled.

.. _`Spglib library`: https://spglib.readthedocs.io/

Input Variables
^^^^^^^^^^^^^^^

The following variables may be set to influence this module's behavior:

``SPGLIB_PREFER_PKGCONFIG``
  if set ``pkg-config`` will be used to search for a Spglib library first
  and if one is found that is preferred

Imported targets
^^^^^^^^^^^^^^^^

This module defines the following :prop_tgt:`IMPORTED` targets:

``Spglib::Spglib``
  .. versionadded:: 3.18

  The libraries to use for Spglib, if found.

Result Variables
^^^^^^^^^^^^^^^^

This module defines the following variables:

``Spglib_FOUND``
  library Spglib is found

``SPGLIB_LIBRARIES``
  uncached list of libraries (using full path name) to link against
  to use Spglib

#]=======================================================================]


# # Check the language being used
if(NOT CMAKE_C_COMPILER_LOADED)
  if(Spglib_FIND_REQUIRED)
    message(FATAL_ERROR "FindSpglib requires C to be enabled.")
  else()
    message(STATUS "Looking for Spglib... - NOT found (Unsupported languages)")
    return()
  endif()
endif()


function(_add_spglib_target)
  if(Spglib_FOUND AND NOT TARGET Spglib::Spglib)
    add_library(Spglib::Spglib INTERFACE IMPORTED)
    if(SPGLIB_LIBRARIES)
      set_target_properties(Spglib::Spglib PROPERTIES
        INTERFACE_LINK_LIBRARIES "${SPGLIB_LIBRARIES}"
      )
    endif()
  endif()
endfunction()


# If Spglib library is given by the user
if(SPGLIB_LIBRARIES)
  find_package_handle_standard_args(Spglib REQUIRED_VARS SPGLIB_LIBRARIES)
  mark_as_advanced(SPGLIB_LIBRARIES)
  _add_spglib_target()
  return()
endif()


# If pkg-config will be used to search Spglib
if(SPGLIB_PREFER_PKGCONFIG)

  find_package(PkgConfig QUIET)
  pkg_check_modules(PKGC_SPGLIB spglib)

  if(PKGC_SPGLIB_FOUND)

    set(Spglib_FOUND ${PKGC_SPGLIB_FOUND})
    set(SPGLIB_LIBRARIES "${PKGC_SPGLIB_LINK_LIBRARIES}")

    _add_spglib_target()

    return()

  endif()

endif()



# If Spglib library is not given by the user and pkg-config won't be used to search Spglib

# Include needed modules
include(CheckFunctionExists)
include(FindPackageHandleStandardArgs)

# Initialize SPGLIB_LIBRARIES
set(SPGLIB_LIBRARIES)

# Specify additional directories to search from CMAKE_C_IMPLICIT_LINK_DIRECTORIES and LD_LIBRARY_PATH environmental variable.
set(_extaddlibdir "")
list(APPEND _extaddlibdir "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
if(WIN32)
  list(APPEND _extaddlibdir ENV LIB)
elseif(APPLE)
  list(APPEND _extaddlibdir ENV DYLD_LIBRARY_PATH)
else()
  list(APPEND _extaddlibdir ENV LD_LIBRARY_PATH)
endif()

# Find libsymspg(.so|.a)
find_library(SPGLIB_LIBRARY NAMES symspg PATHS ${_extaddlibdir})
set(SPGLIB_LIBRARIES "${SPGLIB_LIBRARY}")

# Check if library can be linked
set(CMAKE_REQUIRED_LIBRARIES ${SPGLIB_LIBRARIES})
set(CMAKE_REQUIRED_QUIET ${SPGLIB_FIND_QUIETLY})
check_function_exists("spg_get_symmetry" SPGLIB_WORKS)

# If linking failed, add dependencies and try again
if(NOT SPGLIB_WORKS)
  unset(SPGLIB_WORKS CACHE)
  list(APPEND SPGLIB_LIBRARIES "-lm")
  set(CMAKE_REQUIRED_LIBRARIES ${SPGLIB_LIBRARIES})
  check_function_exists("spg_get_symmetry" SPGLIB_WORKS)
endif()
set(CMAKE_REQUIRED_LIBRARIES)

# Set SPGLIB_LIBRARIES to FALSE if linking failed
if(NOT SPGLIB_WORKS)
  set(SPGLIB_LIBRARIES FALSE)
endif()


find_package_handle_standard_args(Spglib REQUIRED_VARS SPGLIB_LIBRARIES)

_add_spglib_target()
