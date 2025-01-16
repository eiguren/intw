!
! Copyright (C) 2024 INTW group
!
! This file is part of INTW.
!
! INTW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! INTW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!
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
