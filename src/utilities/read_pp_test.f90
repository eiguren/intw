program ep_melements

  use intw_input_parameters
  use intw_reading
  use intw_pseudo
  use intw_useful_constants

!  use intw_uspp, only: nkb, vkb, vkqb, nhtol, nhtolm, indv
!
  use intw_matrix_elements


  !================================================================================
  !       Declare the variables
  !================================================================================
  implicit none

  logical                  :: read_status
  character(256)           :: method

  integer :: is, ir
  !================================================================================
  !       Talk to the user
  !================================================================================
  !
  write(*,20) '====================================================='
  write(*,20) '|                  program me                       |'
  write(*,20) '|        ---------------------------------          |'
  write(*,20) '====================================================='
  write(*,20) '|    waiting for input file...                      |'
  !
  !
  call read_input(read_status)
  !
  method=trim(intw2W_method)
  !
  if (read_status) then
     !
     stop
     !
  endif

  call read_parameters_data_file_xml()

  call read_all_pseudo ()

20 format(A)
end program  ep_melements

     !
