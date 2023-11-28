#############################################################################
#
#############################################################################

function (mangle_fortran_name CNAME FNAME)
    set (TMP)
    if (WIN32)
        string (TOUPPER "${FNAME}" TMP)
    else ()
        string (TOLOWER "${FNAME}_" TMP)
    endif ()
    set (${CNAME} ${TMP} PARENT_SCOPE)
endfunction ()

function (mangle_fortran_filename_list MANGLED)
    set (TMP)
    foreach (TFILE ${ARGN})
        string (REGEX REPLACE ".f90$" "" TESTNAME ${TFILE})
        mangle_fortran_name (C_TESTNAME ${TESTNAME}_test)
        list (APPEND TMP ${C_TESTNAME})
    endforeach ()
    set (${MANGLED} ${TMP} PARENT_SCOPE)
endfunction()

function (add_fortran_test_executable TARGET)
    set (TEST_FILES ${ARGN})
    mangle_fortran_filename_list (TEST_FILES_MANGLED ${TEST_FILES})

    create_test_sourcelist (_ Testing/${TARGET}.c ${TEST_FILES_MANGLED})

    add_library (${TARGET}_lib ${TEST_FILES})
    target_link_libraries(${TARGET}_lib PRIVATE intw_modules)
  	target_link_libraries(${TARGET}_lib PRIVATE test_modules)
  	target_link_libraries(${TARGET}_lib PRIVATE ${LAPACK_LIBRARIES})
    if(USE_MPI)
        target_link_libraries(${TARGET}_lib PRIVATE MPI::MPI_Fortran)
    endif()

    add_executable (${TARGET}.x Testing/${TARGET}.c)
    target_link_libraries (${TARGET}.x PRIVATE ${TARGET}_lib)

    set (INDEX 0)
    list (LENGTH TEST_FILES LEN)
    while (${LEN} GREATER ${INDEX})
        list (GET TEST_FILES ${INDEX} TEST_FILE)
        list (GET TEST_FILES_MANGLED ${INDEX} TEST_FILE_MANGLED)
        add_test (
            NAME ${TARGET}/${TEST_FILE}
            COMMAND $<TARGET_FILE:${TARGET}.x> ${TEST_FILE_MANGLED}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
        math (EXPR INDEX "${INDEX} + 1")
    endwhile ()
endfunction ()
