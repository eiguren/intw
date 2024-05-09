


set(PROGS
  nscf
  ep_melements
  read_pp_test
  intw2W90
  w902intw
  interpolate
  interpolate_phonons
  ep_interp_on_trFS_wannier
  ep_interp_on_trFS_nufft
  a2F_on_trFS
)

list(REMOVE_DUPLICATES PROGS)



foreach(prog IN ITEMS ${PROGS})
  add_executable(${prog}.x ${prog}.f90)
  install(TARGETS ${prog}.x DESTINATION bin)

  target_link_libraries(${prog}.x PRIVATE intw_modules)
  target_link_libraries(${prog}.x PRIVATE ${LAPACK_LIBRARIES})
  if(USE_MPI)
    target_link_libraries(${prog}.x PRIVATE MPI::MPI_Fortran)
  endif()
endforeach()
