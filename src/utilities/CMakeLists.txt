


set(PROGS
  ep_melements
  ep_melements_noW
  read_pp_test
  intw2W90_2
)

list(REMOVE_DUPLICATES PROGS)



foreach(prog IN ITEMS ${PROGS})
  add_executable(${prog}.x ${prog}.f90)
  install(TARGETS ${prog}.x DESTINATION bin)

  target_link_libraries( ${prog}.x PRIVATE W_lib)
  target_link_libraries(${prog}.x PRIVATE intw_modules)
  target_link_libraries(${prog}.x PRIVATE ${LAPACK_LIBRARIES})
  target_link_libraries(${prog}.x PRIVATE MPI::MPI_Fortran)
endforeach()
