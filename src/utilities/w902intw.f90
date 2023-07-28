program w902intw

    use kinds, only: dp
    use intw_w90_setup, only: read_w90_chk, read_eig, allocate_and_build_u_mesh, write_formatted_u_mesh, &
                              allocate_and_build_ws_irvec, allocate_and_build_ham_k, &
                              allocate_and_build_ham_r, write_ham_r, &
                              eigenval_intw
    use intw_input_parameters, only: nk1, nk2, nk3, mesh_dir, prefix, read_input
    use intw_reading, only: read_parameters_data_file_xml, set_num_bands, &
                            num_bands_intw
    use intw_intw2wannier, only: read_nnkp_file, &
                                 nnkp_num_kpoints

!================================================================================
!       Declare the variables
!================================================================================
    implicit none

    character(256) :: nnkp_file

    logical        :: read_status


!--------------------------------------------------------------------------------
!================================================================================
!       Talk to the user
!================================================================================
write(*,'(A)') '====================================================='
write(*,'(A)') '|                  program w902intw                 |'
write(*,'(A)') '|        ---------------------------------          |'
write(*,'(A)') '====================================================='
write(*,'(A)') '|    waiting for input file...                      |'

!================================================================================
!       read the input file
!       Read in the necessary information from standard input
!================================================================================
    call read_input(read_status)
 
    if (read_status ) then
           stop
    end if


!================================================================================
!       read the parameters from the SCF QE calculation
!================================================================================
    call read_parameters_data_file_xml()
   
    ! Set the number of bands for the calculation
    call set_num_bands()


!================================================================================
!       read .nnkp file from w90
!================================================================================
    nnkp_file = trim(mesh_dir)//trim(prefix)//".nnkp"
    call read_nnkp_file(nnkp_file)


!================================================================================
!       read .chk file, build u matrix and write to file
!================================================================================
    allocate(eigenval_intw(num_bands_intw,nnkp_num_kpoints))
    call read_eig(eigenval_intw)
    call allocate_and_build_u_mesh()

    call write_formatted_u_mesh(eigenval_intw)


!================================================================================
!       build Wigner-Seitz cell
!================================================================================
    call allocate_and_build_ws_irvec(nk1,nk2,nk3)


!================================================================================
!       build H(k) and H(R), write H(R) to file
!================================================================================
    call allocate_and_build_ham_k(eigenval_intw)
    call allocate_and_build_ham_r()

    call write_ham_r()


!================================================================================
!       clean up and finish
!================================================================================
    deallocate(eigenval_intw)

    write(*,'(A)') '====================================================='
    write(*,'(A)') '|               end program w902intw                |'
    write(*,'(A)') '|        ---------------------------------          |'
    write(*,'(A)') '====================================================='


end program w902intw
