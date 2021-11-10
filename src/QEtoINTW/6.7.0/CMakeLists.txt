include(ExternalProject)

add_executable(pw2intw.x pw2intw.f90)


target_link_libraries(pw2intw.x PRIVATE W_lib)
target_link_libraries(pw2intw.x PRIVATE QE_lib)
target_link_libraries(pw2intw.x PRIVATE ${LAPACK_LIBRARIES})
target_link_libraries(pw2intw.x PRIVATE MPI::MPI_Fortran)
