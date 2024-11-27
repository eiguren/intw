


set(PROGS
  nscf
  ep_melements
  read_pp_test
  intw2W90
  w902intw
  interpolate
  interpolate_phonons
  ep_interp_on_trFS_wannier
  ep_interp_on_trFS_dV
  a2F_on_trFS
)

list(REMOVE_DUPLICATES PROGS)



foreach(prog IN ITEMS ${PROGS})
  add_executable(${prog}.x ${prog}.f90)
  install(TARGETS ${prog}.x DESTINATION bin)

  target_link_libraries(${prog}.x PRIVATE intw_modules)
  target_link_libraries(${prog}.x PRIVATE ${LAPACK_LIBRARIES})
  if(USE_OPENMP)
    target_link_libraries(${prog}.x PRIVATE OpenMP::OpenMP_Fortran)
  endif()
endforeach()


if(TRIFS)

  add_executable(triFS.x triFS.f90)
  install(TARGETS triFS.x DESTINATION bin)

  target_link_libraries(triFS.x PRIVATE intw_modules)
  target_link_libraries(triFS.x PRIVATE triFS_modules)
  target_link_libraries(triFS.x PRIVATE ${LAPACK_LIBRARIES})
  if(USE_OPENMP)
    target_link_libraries(triFS.x PRIVATE OpenMP::OpenMP_Fortran)
  endif()

endif()