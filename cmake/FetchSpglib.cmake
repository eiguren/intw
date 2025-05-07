include(FetchContent)

# Default spglib git repository and tag
set(REPOSITORY https://github.com/spglib/spglib)
set(TAG v2.6.0)

# Allow users to pass spglib git repository and tag via flags
if(DEFINED "SPGLIB_GIT_REPOSITORY")
  set(REPOSITORY ${SPGLIB_GIT_REPOSITORY})
endif()

if(DEFINED "SPGLIB_GIT_TAG")
  set(TAG ${SPGLIB_GIT_TAG})
endif()

FetchContent_Declare(
  Spglib

  GIT_REPOSITORY ${REPOSITORY}
  GIT_TAG ${TAG}


  SOURCE_DIR "${PROJECT_BINARY_DIR}/_deps/spglib_${TAG}-src"
  BINARY_DIR "${PROJECT_BINARY_DIR}/_deps/spglib_${TAG}-build"
  SUBBUILD_DIR "${PROJECT_BINARY_DIR}/_deps/spglib_${TAG}-subbuild"

  OVERRIDE_FIND_PACKAGE
  )

  set(SPGLIB_USE_OMP OFF)
  set(SPGLIB_WITH_TESTS OFF)
  set(SPGLIB_WITH_Fortran ON)
  set(SPGLIB_SHARED_LIBS OFF)

  FetchContent_MakeAvailable(Spglib)

set(SPGLIB_LIBRARIES "${Spglib_BINARY_DIR}/libsymspg.a")
