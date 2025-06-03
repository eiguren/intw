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
program w902intw

  use kinds, only: dp

  use intw_utility, only: get_timing

  use intw_w90_setup, only: read_w90_chk, read_eig, allocate_and_build_u_mesh, write_formatted_u_mesh, &
                            allocate_and_build_ws_irvec, allocate_and_build_ham_k, &
                            allocate_and_build_ham_r, write_ham_r

  use intw_input_parameters, only: nk1, nk2, nk3, mesh_dir, prefix, read_input

  use intw_reading, only: read_parameters_data_file_xml, set_num_bands

  use intw_intw2wannier, only: read_nnkp_file

  implicit none

  character(256) :: nnkp_file
  logical :: have_nnkp

  logical :: read_status

  ! timing
  real(dp) :: time1, time2


  20 format(A)
  30 format(A,F8.2,6X,A)


  !================================================================================
  ! Beginning
  !================================================================================

  call get_timing(time1)

  write(*,20) '====================================================='
  write(*,20) '|                  program w902intw                 |'
  write(*,20) '|         ---------------------------------         |'
  write(*,20) '====================================================='


  !================================================================================
  ! Read the input file
  !================================================================================

  call read_input(read_status)

  if (read_status) stop


  !================================================================================
  ! Check that $prefix.nnkp is present
  !================================================================================

  nnkp_file = trim(mesh_dir)//trim(prefix)//".nnkp"

  inquire(file=nnkp_file, exist=have_nnkp)

  if(.not. have_nnkp) then
    write(*,20) '**********************************************************'
    write(*,20) '* Could not find the file '//trim(nnkp_file)
    write(*,20) '* Did you run W90 -pp $seed to get the parameter file?   '
    write(*,20) '**********************************************************'
    stop
  end if
  write(*,20) '| - .nnkp file found                                |'
  write(*,20) '|         ---------------------------------         |'


  !================================================================================
  ! Read the parameters from the SCF calculation
  !================================================================================

  write(*,20) '| - Reading calculation parameters...               |'

  call read_parameters_data_file_xml()


  !================================================================================
  ! Set the number of wave functions
  !================================================================================

  call set_num_bands()


  !================================================================================
  ! Read the .nnkp file
  !================================================================================

  write(*,20) '| - Reading nnkp file...                            |'

  call read_nnkp_file(nnkp_file)

  write(*,20) '|         ---------------------------------         |'


  !================================================================================
  ! Read the .eig file
  !================================================================================

  write(*,20) '| - Reading eig file...                             |'

  call read_eig()


  !================================================================================
  ! Read the .chk file, build u matrix and write to file
  !================================================================================

  write(*,20) '| - Reading and building Wannier U matrix...        |'

  call allocate_and_build_u_mesh()

  write(*,20) '| - Saving Wannier U matrix...                      |'

  call write_formatted_u_mesh()


  !================================================================================
  ! Build Wigner-Seitz cell
  !================================================================================

  write(*,20) '| - Building WS mesh...                             |'

  call allocate_and_build_ws_irvec(nk1, nk2, nk3)


  !================================================================================
  ! Build H(k) and H(R), write H(R) to file
  !================================================================================

  write(*,20) '| - Building Wannier H(k) and H(R)...               |'

  call allocate_and_build_ham_k()
  call allocate_and_build_ham_r()

  write(*,20) '| - Saving Wannier H(R)...                          |'

  call write_ham_r()

  write(*,20) '====================================================='


  !================================================================================
  ! Finish
  !================================================================================

  call get_timing(time2)

  write(*,20) '|                      ALL DONE                     |'
  write(*,30) '|     total time: ',time2-time1,' seconds            |'
  write(*,20) '====================================================='

end program w902intw
