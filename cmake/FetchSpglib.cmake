include(FetchContent)

# Default spglib git repository and tag
set(REPOSITORY https://github.com/spglib/spglib)
set(TAG v2.6.0)

# Allow users to pass spglib git repository and tag via flags
if(DEFINED "SPGLIB_GIT_REPOSITORY")
  set(REPOSITORY ${SPGLIB_GIT_REPOSITORY})
endif()

if(DEFINED "SPGLIB_GIT_TAG")
  # Validate the tag format (simple regex for vX.Y.Z)
  if(NOT SPGLIB_GIT_TAG MATCHES "^v[0-9]+\\.[0-9]+\\.[0-9]+$")
    message(FATAL_ERROR "Invalid SPGLIB_GIT_TAG format. Please specify a valid release Git tag (e.g., v2.6.0).")
  endif()
  set(TAG ${SPGLIB_GIT_TAG})
endif()


# Function to compare version strings
function(version_is_newer_or_equal version1 version2)

    # Remove the 'v' prefix if it exists
    string(REPLACE "v" "" version1_no_v ${version1})
    string(REPLACE "v" "" version2_no_v ${version2})

    # Split the version strings into components
    string(REPLACE "." ";" version1_list ${version1_no_v})
    string(REPLACE "." ";" version2_list ${version2_no_v})

    # Compare each component
    foreach(i RANGE 0 2)
        list(GET version1_list ${i} v1_component)
        list(GET version2_list ${i} v2_component)

        # If version1 is greater, return true
        if(v1_component GREATER v2_component)
            return()
        endif()

        # If version1 is less, return false
        if(v1_component LESS v2_component)
            message(FATAL_ERROR "Version ${version1} is older than ${version2}.")
            return()
        endif()
    endforeach()

    # If we reach here, the versions are equal
endfunction()

# Check if the specified version is 2.1.0 or newer
version_is_newer_or_equal(${TAG} "v2.1.0")


FetchContent_Declare(
  Spglib

  GIT_REPOSITORY ${REPOSITORY}
  GIT_TAG ${TAG}


  # PREFIX "${PROJECT_BINARY_DIR}/_deps"
  SOURCE_DIR "${PROJECT_BINARY_DIR}/_deps/spglib-src"
  BINARY_DIR "${PROJECT_BINARY_DIR}/_deps/spglib-build"
  SUBBUILD_DIR "${PROJECT_BINARY_DIR}/_deps/spglib-subbuild"

  OVERRIDE_FIND_PACKAGE
  )

set(SPGLIB_USE_OMP OFF)
set(SPGLIB_WITH_TESTS OFF)
set(SPGLIB_WITH_Fortran ON)
set(SPGLIB_SHARED_LIBS OFF)

FetchContent_MakeAvailable(Spglib)

# Set SPGLIB_LIBRARIES and SPGLIB_INCLUDE_DIRS
get_target_property(SPGLIB_SYMSPG_BINARY_DIR Spglib_symspg BINARY_DIR)
get_target_property(SPGLIB_SYMSPG_OUTPUT_NAME Spglib_symspg OUTPUT_NAME)

get_target_property(SPGLIB_FORTRAN_BINARY_DIR Spglib_fortran BINARY_DIR)
get_target_property(SPGLIB_FORTRAN_OUTPUT_NAME Spglib_fortran OUTPUT_NAME)

set(SPGLIB_LIBRARIES "${SPGLIB_SYMSPG_BINARY_DIR}/lib${SPGLIB_SYMSPG_OUTPUT_NAME}.a"
                     "${SPGLIB_FORTRAN_BINARY_DIR}/lib${SPGLIB_FORTRAN_OUTPUT_NAME}.a"
                     )

get_target_property(SPGLIB_INCLUDE_DIRS Spglib_fortran Fortran_MODULE_DIRECTORY)
