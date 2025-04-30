include(FetchContent)


FetchContent_Declare(
  Spglib

  GIT_REPOSITORY https://github.com/spglib/spglib
  GIT_TAG v2.6.0


  SOURCE_DIR "${PROJECT_BINARY_DIR}/External/spglib-src"
  BINARY_DIR "${PROJECT_BINARY_DIR}/External/spglib-build"
  SUBBUILD_DIR "${PROJECT_BINARY_DIR}/External/spglib-subbuild"

  OVERRIDE_FIND_PACKAGE
  )

  set(SPGLIB_USE_OMP OFF)
  set(SPGLIB_WITH_TESTS OFF)
  set(SPGLIB_WITH_Fortran ON)
  set(SPGLIB_SHARED_LIBS OFF)

  FetchContent_MakeAvailable(Spglib)

set(SPGLIB_LIBRARIES "${Spglib_BINARY_DIR}/libsymspg.a")
