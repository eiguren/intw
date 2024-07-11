

if(USE_SIESTA)

  set(SIESTA2INTW_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION} CACHE PATH "Path to the siesta2intw src directory." FORCE)
  set(SIESTA2INTW_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${SIESTA_VERSION} CACHE PATH "Path to the siesta2intw build directory." FORCE)


  set(SIESTA2INTW_SRC_FILES ${SIESTA2INTW_SRC_DIR}/cfftnd.f90
                            ${SIESTA2INTW_SRC_DIR}/siesta2intw.F
                            ${SIESTA2INTW_SRC_DIR}/siesta2intw_io.f90
                            ${SIESTA2INTW_SRC_DIR}/siesta2intw_ph.f90
                            ${SIESTA2INTW_SRC_DIR}/siesta2intw_symmetry.f90
                            ${SIESTA2INTW_SRC_DIR}/spglib_f08.f90
                            ${SIESTA2INTW_SRC_DIR}/siesta2intw_write.F90
                            ${SIESTA2INTW_SRC_DIR}/siesta2intw_fft.f90
                            ${SIESTA2INTW_SRC_DIR}/siesta2intw_pp.f90
                            ${SIESTA2INTW_SRC_DIR}/siesta2intw_utils.F90
                            )
  set(SIESTA2INTW_OBJ_FILES ${SIESTA2INTW_BUILD_DIR}/cfftnd.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_io.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_ph.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_symmetry.o
                            ${SIESTA2INTW_BUILD_DIR}/spglib_f08.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_write.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_fft.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_pp.o
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_utils.o
                            )
  set(SIESTA2INTW_MOD_FILES ${SIESTA2INTW_BUILD_DIR}/fftpack5.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_fft.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_io.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_ph.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_pp.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_symmetry.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_utils.mod
                            ${SIESTA2INTW_BUILD_DIR}/spglib_f08.mod
                            ${SIESTA2INTW_BUILD_DIR}/siesta2intw_write.mod
                            )


  # Add some debugging flags
  string(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
  if(BT STREQUAL "DEBUG")
      set(SIESTA2INTW_FFLAGS ${CMAKE_Fortran_FLAGS_DEBUG})
      set(SIESTA2INTW_FPPFLAGS -DDEBUG)
  endif()

  list(JOIN SPGLIB_LIBRARIES " " SPGLIB_LIBRARIES_STR)
  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${SIESTA_VERSION}/Makefile.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/${SIESTA_VERSION}/Makefile" @ONLY)

  add_custom_command(OUTPUT ${SIESTA2INTW_BUILD_DIR}/siesta2intw.x
    COMMAND make -C ${SIESTA2INTW_BUILD_DIR} siesta2intw.x
    COMMENT "Building siesta2intw"
    DEPENDS ${SIESTA2INTW_SRC_FILES}
    BYPRODUCTS ${SIESTA2INTW_OBJ_FILES} ${SIESTA2INTW_MOD_FILES}
    )

  add_custom_target(siesta2intw.x ALL
    DEPENDS ${SIESTA2INTW_BUILD_DIR}/siesta2intw.x
    )

endif()
