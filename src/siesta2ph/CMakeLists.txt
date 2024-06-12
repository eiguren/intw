

if(USE_SIESTA)


  # If SPGLIB is not specified by the user try to find it
  if(NOT SPGLIB)

    # Try to use PkgConfig first
    find_package(PkgConfig)

    if(PKG_CONFIG_FOUND)
      pkg_check_modules(SPGLIB spglib)
    else()
      set(SPGLIB_FOUND 0)
    endif()

    if(SPGLIB_FOUND)
      set(SPGLIB ${SPGLIB_LINK_LIBRARIES})
    else()
      message(STATUS "Trying to find Spglib")
      set(_extaddlibdir "")
      list(APPEND _extaddlibdir ENV LD_LIBRARY_PATH)
      find_library(SPGLIB NAMES symspg PATHS ${_extaddlibdir})
      if(SPGLIB)
        set(SPGLIB_FOUND 1)
      endif()
    endif()

    if(SPGLIB_FOUND)
      set(SPGLIB ${SPGLIB} CACHE PATH "Spglib library.")
    else()
      message(FATAL_ERROR "Spglib not found!")
    endif()

  endif()


  set(SIESTA2PH_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION} CACHE PATH "Path to the siesta2ph src directory." FORCE)
  set(SIESTA2PH_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${SIESTA_VERSION} CACHE PATH "Path to the siesta2ph build directory." FORCE)


  set(SIESTA2PH_SRC_FILES ${SIESTA2PH_SRC_DIR}/siesta2ph_io.f90
                          ${SIESTA2PH_SRC_DIR}/siesta2ph_linalg.F90
                          ${SIESTA2PH_SRC_DIR}/spglib_f08.f90
                          ${SIESTA2PH_SRC_DIR}/siesta2ph_symmetry.f90
                          ${SIESTA2PH_SRC_DIR}/siesta2ph_system.f90
                          ${SIESTA2PH_SRC_DIR}/siesta2ph_utils.f90
                          )
  set(SIESTA2PH_OBJ_FILES ${SIESTA2PH_BUILD_DIR}/siesta2ph_io.o
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_linalg.o
                          ${SIESTA2PH_BUILD_DIR}/spglib_f08.o
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_symmetry.o
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_system.o
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_utils.o
                          )
  set(SIESTA2PH_MOD_FILES ${SIESTA2PH_BUILD_DIR}/siesta2ph_io.mod
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_linalg.mod
                          ${SIESTA2PH_BUILD_DIR}/spglib_f08.mod
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_symmetry.mod
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_system.mod
                          ${SIESTA2PH_BUILD_DIR}/siesta2ph_utils.mod
                          )


  # Add some debugging flags
  string(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
  if(BT STREQUAL "DEBUG")
      set(SIESTA2PH_FFLAGS ${CMAKE_Fortran_FLAGS_DEBUG})
      set(SIESTA2PH_FPPFLAGS -DDEBUG)
  endif()

  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/Makefile.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/${SIESTA_VERSION}/Makefile" @ONLY)

  add_custom_command(OUTPUT ${SIESTA2PH_BUILD_DIR}/siesta2ph.x
    COMMAND make -C ${SIESTA2PH_BUILD_DIR} siesta2ph.x
    COMMENT "Building siesta2ph.x"
    DEPENDS ${SIESTA2PH_SRC_FILES} ${SIESTA2PH_SRC_DIR}/siesta2ph.F90
    BYPRODUCTS ${SIESTA2PH_OBJ_FILES} ${SIESTA2PH_MOD_FILES} ${SIESTA2PH_BUILD_DIR}/siesta2ph.o
    )

  add_custom_command(OUTPUT ${SIESTA2PH_BUILD_DIR}/siesta2fc.x
    COMMAND make -C ${SIESTA2PH_BUILD_DIR} siesta2fc.x
    COMMENT "Building siesta2fc.x"
    DEPENDS ${SIESTA2PH_SRC_FILES} ${SIESTA2PH_SRC_DIR}/siesta2fc.F90
    BYPRODUCTS ${SIESTA2PH_OBJ_FILES} ${SIESTA2PH_MOD_FILES} ${SIESTA2PH_BUILD_DIR}/siesta2fc.o
    )

  add_custom_command(OUTPUT ${SIESTA2PH_BUILD_DIR}/siesta2dv.x
    COMMAND make -C ${SIESTA2PH_BUILD_DIR} siesta2dv.x
    COMMENT "Building siesta2dv.x"
    DEPENDS ${SIESTA2PH_SRC_FILES} ${SIESTA2PH_SRC_DIR}/siesta2dv.F90
    BYPRODUCTS ${SIESTA2PH_OBJ_FILES} ${SIESTA2PH_MOD_FILES} ${SIESTA2PH_BUILD_DIR}/siesta2dv.o
    )

  add_custom_target(siesta2ph.x ALL
    DEPENDS ${SIESTA2PH_BUILD_DIR}/siesta2ph.x
    )

  add_custom_target(siesta2fc.x ALL
    DEPENDS ${SIESTA2PH_BUILD_DIR}/siesta2fc.x
    )

  add_custom_target(siesta2dv.x ALL
    DEPENDS ${SIESTA2PH_BUILD_DIR}/siesta2dv.x
    )

endif()
