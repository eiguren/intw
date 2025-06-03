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
program interpolate

  ! MBR 2024
  ! Uses nnkp and ham_r files to write an interpolated band structure and DOS.
  ! Projections on Wannier functions are also provided, so that "fatband"
  ! and PDOS plots can be made.

  use kinds, only: dp

  use intw_utility, only: get_timing, find_free_unit, generate_and_allocate_kpath, fermi_dirac

  use intw_w90_setup, only: interpolate_1k, interpolated_DOS, &
                            allocate_and_read_ham_r

  use intw_input_parameters, only: nk1, nk2, nk3, mesh_dir, prefix, read_input, &
                                   read_cards, exist_kpath, nkpath, nkspecial, kspecial, &
                                   nk1_dos, nk2_dos, nk3_dos, ne_dos, eini_dos, efin_dos, esmear_dos, ktsmear, &
                                   chemical_potential

  use intw_reading, only: read_parameters_data_file_xml, set_num_bands, &
                          num_wann_intw, &
                          at, bg, tpiba

  implicit none

  character(256) :: bpath_file, dos_file
  logical :: read_status
  integer :: io_unit_b, io_unit_d, ik, ie, iw
  integer, allocatable :: kspecial_indices(:)
  real(dp) :: kpoint(3)
  real(dp) :: estep, ener, num_elec
  real(dp), allocatable :: eig_int(:), DOS(:), PDOS(:,:), weights(:)
  real(dp), allocatable :: kpath(:,:), dkpath(:)
  complex(dp), allocatable :: u_int(:,:)

  ! timing
  real(dp) :: time1, time2


  20 format(A)
  30 format(A,F8.2,6X,A)
  !
  !
  !================================================================================
  ! Beginning
  !============================================================================20
  !
  call get_timing(time1)
  !
  write(*,20) '====================================================='
  write(*,20) '|                program interpolate                |'
  write(*,20) '|         ---------------------------------         |'
  write(*,20) '====================================================='
  !
  !
  !================================================================================
  ! Read the input file
  !================================================================================
  !
  call read_input(read_status)
  !
  if (read_status) stop
  !
  ! K_PATH
  call read_cards()
  if (.not.exist_kpath) then
    write(*,*) 'K_PATH not found. Bands/DOS cannot be interpolated. Stopping.'
    stop
  end if
  !
  !
  !================================================================================
  ! Read the parameters from the SCF calculation
  !================================================================================
  !
  write(*,20) '| - Reading calculation parameters...               |'
  !
  ! From this, I only need seednames, actually,
  ! until intw.in is modified to include postprocessing options
  call read_parameters_data_file_xml()
  !
  !
  !================================================================================
  ! Set the number of bands for the calculation
  !================================================================================
  !
  ! this will read num_wann_intw and num_bands_intw dimensions from .nnkp if available
  call set_num_bands()
  !
  !
  !================================================================================
  ! Read u_mesh, ham_r file from w902intw
  !================================================================================
  !
  write(*,20) '| - Reading Wannier U matrix and H(R)...            |'
  !
  call allocate_and_read_ham_r()
  !
  write(*,20) '|         ---------------------------------         |'
  !
  !
  !================================================================================
  ! Build kpoint path to plot bands.
  ! The nkpath number of points from the input might fluctuate.
  ! Use kspecial_indices option to print out the special k-points
  ! along the path (useful for plotting).
  !================================================================================
  !
  write(*,20) '| - Building k-path...                              |'
  !
  call generate_and_allocate_kpath(at, bg, tpiba, nkpath, nkspecial, kspecial, &
                                   kpath, dkpath, kspecial_indices)
  !
  ! write(*,*) nkpath
  ! do ik=1,nkpath
  !   write(*,'(3f12.6)') kpath(:,ik)
  ! end do
  !
  !
  !================================================================================
  ! Interpolate bands over path
  !================================================================================
  !
  write(*,20) '| - Interpolating bands...                          |'
  !
  allocate(eig_int(num_wann_intw))
  allocate(weights(num_wann_intw))
  allocate(u_int(num_wann_intw,num_wann_intw))
  !
  io_unit_b = find_free_unit()
  bpath_file = trim(mesh_dir)//trim(prefix)//trim('.bnd_int')
  open(unit=io_unit_b, file=bpath_file, status='unknown')
  write(io_unit_b,'(A)') '# k-in-path    energy(eV)      weight(1:num_wann_intw)'
  !
  do iw = 1, num_wann_intw
    do ik = 1, nkpath
      kpoint = kpath(:,ik)
      call interpolate_1k(kpoint, eig_int, u_int)
      weights = (abs(u_int(iw,:)))**2
      write(io_unit_b,'(40e14.6)') dkpath(ik), eig_int(iw), weights
    end do
    write(io_unit_b,*)
  end do
  !
  ! Write special point indices
  write(io_unit_b,*) '#'
  write(io_unit_b,*) '#Special k-points in the .bnd_int file are:'
  do ik=1,nkspecial
    write(io_unit_b,'(a,3f10.4,a,i4,e14.6)') '#', kspecial(:,ik), ' --> ', kspecial_indices(ik), dkpath(kspecial_indices(ik))
  end do
  write(io_unit_b,*) '#'
  !
  close(io_unit_b)
  !
  write(*,20) '|   Band interpolation finished and written to:     |'
  write(*,20) "|   "//bpath_file(1:max(47,len(trim(bpath_file))))//" |"
  !
  write(*,20) '|         ---------------------------------         |'
  !
  !
  !================================================================================
  ! DOS over a fine grid
  ! Parameters of DOS plot from namelist /DOS/
  !================================================================================
  !
  write(*,20) '| - Computing DOS...                                |'
  !
  allocate(DOS(ne_dos), PDOS(ne_dos,num_wann_intw))
  !
  io_unit_d = find_free_unit()
  dos_file = trim(mesh_dir)//trim(prefix)//trim('.dos_int')
  open(unit=io_unit_d, file=dos_file, status='unknown')
  write(io_unit_d,'(A)') '# E-E_Fermi(eV)      DOS          PDOS(1:num_wann_intw)'
  !
  call interpolated_DOS(nk1_dos, nk2_dos, nk3_dos, eini_dos, efin_dos, esmear_dos, ne_dos, DOS, PDOS)
  !
  estep = (efin_dos-eini_dos)/real(ne_dos-1,dp)
  num_elec = 0.0_dp
  do ie = 1, ne_dos
    ener = eini_dos + (ie-1)*estep
    write(io_unit_d,'(40e14.6)') ener-chemical_potential, DOS(ie), PDOS(ie,1:num_wann_intw)
    num_elec = num_elec + DOS(ie)*fermi_dirac(ener-chemical_potential, ktsmear)
  end do
  !
  close(io_unit_d)
  !
  write(*,20) '| - DOS sum test:                                   |'
  write(*,'(A37,F10.6,5X,A1)') '|   DOS integral up to Fermi level = ', num_elec, '|'
  !
  write(*,20) '|   DOS computed and written to:                    |'
  write(*,20) '|   '//dos_file(1:max(47,len(trim(dos_file))))//' |'
  !
  write(*,20) '====================================================='
  !
  !
  !================================================================================
  ! Finish
  !================================================================================
  !
  call get_timing(time2)
  !
  write(*,20) '|                      ALL DONE                     |'
  write(*,30) '|     total time: ',time2-time1,' seconds            |'
  write(*,20) '====================================================='

end program interpolate
