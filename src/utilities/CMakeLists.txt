


set(PROGS
        ep_melements
        read_pp_test
	intw2W90_2
	)

##list(APPEND PROGS ${PROGS_HARITZ})
list(REMOVE_DUPLICATES PROGS)


foreach(prog IN ITEMS ${PROGS})
	add_executable(${prog}.x ${prog}.f90)
	install(TARGETS ${prog}.x DESTINATION bin)
	target_include_directories( ${prog}.x PRIVATE ${W_HOME}/src)
        target_include_directories( ${prog}.x PRIVATE ${W_HOME}/src/obj)
        target_include_directories( ${prog}.x PRIVATE ${W_HOME}/src/objp)

	target_link_libraries(${prog}.x PRIVATE intw_modules)
	#target_link_libraries(${prog}.x PRIVATE ${QE_HOME}/Modules/libqemod.a ${QE_HOME}/PW/libpw.a ${QE_HOME}/flib/flib.a ${QE_HOME}/flib/ptools.a ${QE_HOME}/clib/clib.a ${QE_HOME}/iotk/src/libiotk.a)
        # target_link_libraries(${prog}.x PRIVATE ${QE_HOME}/flib/flib.a ${QE_HOME}/flib/ptools.a ${QE_HOME}/Modules/libqemod.a  ${QE_HOME}/iotk/src/libiotk.a )

        
	target_link_libraries(${prog}.x PRIVATE ${W_HOME}/libwannier.a)
	target_link_libraries(${prog}.x PRIVATE ${LAPACK_LIBRARIES})
	# this only links the libraries specified by MPI_Fortran_LIBRARIES (adding their paths to rpath)
	#link_libraries(${MPI_Fortran_LIBRARIES})
	# this two forms are equivalent, and they add also the LINK_FLAGS (because all the other variables are not set)
	# if MPI_Fortran_COMPILER_FLAGS where also set, I think that the first way will also link the targets with that variable
	target_link_libraries(${prog}.x PRIVATE MPI::MPI_Fortran)
        target_compile_definitions( ${prog}.x PRIVATE __MPI)
	#link_libraries(${MPI_Fortran_LIBRARIES} ${MPI_Fortran_LINK_FLAGS})
endforeach()

