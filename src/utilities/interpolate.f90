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
    use intw_utility, only: find_free_unit, generate_and_allocate_kpath, fermi_dirac
    use intw_w90_setup, only: interpolate_1k, interpolated_DOS, &
            allocate_and_read_ham_r, allocate_and_read_u_mesh
    use intw_input_parameters, only: nk1, nk2, nk3, mesh_dir, prefix, read_input, &
            read_cards, exist_kpath, nkpath, nkspecial, kspecial, &
            nk1_dos, nk2_dos, nk3_dos, ne_dos, eini_dos, efin_dos, esmear_dos, ktsmear, &
            chemical_potential
    use intw_reading, only: read_parameters_data_file_xml, set_num_bands, &
                            num_bands_intw, num_wann_intw, &
                            at, bg, tpiba

!================================================================================
!       Declare the variables
!================================================================================
    implicit none

    character(256) :: nnkp_file, bpath_file, dos_file
    logical        :: read_status
    integer :: io_unit_k, io_unit_b, io_unit_d, ik, ie, iw
    integer , allocatable :: kspecial_indices(:)
    real(dp) :: kpoint(3), kpoint_dist
    real(dp) :: estep, ener, num_elec
    real(dp) , allocatable :: eig_int(:), DOS(:), PDOS(:,:), weights(:)
    real(dp) , allocatable :: kpath(:,:), dkpath(:)
    complex(dp) , allocatable :: u_int(:,:)


!--------------------------------------------------------------------------------
!================================================================================
!       Talk to the user
!================================================================================
write(*,'(A)') '====================================================='
write(*,'(A)') '|              program interpolate                  |'
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

    ! K_PATH
    call read_cards
    if ( .not. exist_kpath) then
       write(*,*)' K_PATH not found. Bands/DOS cannot be interpolated. Stopping.'
       stop
    end if


!================================================================================
!       read the parameters from the SCF QE calculation
!
! From this, I only need seednames, actually,
! until intw.in is modified to include postprocessing options
!================================================================================
    call read_parameters_data_file_xml()

!================================================================================
!   Build kpoint path to plot bands.
!   The nkpath number of points from the input might fluctuate.
!   Use kspecial_indices option to print out the special k-points
!   along the path (useful for plotting).
!================================================================================
    call generate_and_allocate_kpath (at, bg, tpiba, nkpath, nkspecial, kspecial, kpath, dkpath, &
            kspecial_indices)
    !write(*,*) nkpath
    !do ik=1,nkpath
    !   write(*,'(3f12.6)') kpath(:,ik)
    !end do


!================================================================================
! Set the number of bands for the calculation
!================================================================================

    ! (this will read num_wann_intw and num_bands_intw dimensions from .nnkp if available)
    call set_num_bands()

!================================================================================
!       read u_mesh, ham_r file from w902intw
!================================================================================

    call allocate_and_read_ham_r()

!================================================================================
!   Test: interpolate bands over path
!================================================================================

    allocate (eig_int(num_wann_intw))
    allocate (weights(num_wann_intw))
    allocate (u_int(num_wann_intw,num_wann_intw))

    io_unit_b = find_free_unit()
    bpath_file = trim(mesh_dir)//trim(prefix)//trim('.bnd_int')
    open(unit=io_unit_b,file=bpath_file,status='unknown')
    write(io_unit_b,'(A)') '# k-in-path    energy(eV)      weight(1:num_wann_intw)'

    do ik = 1,nkpath
       kpoint=kpath(:,ik)
       call interpolate_1k (kpoint, eig_int, u_int)
       do iw = 1,num_wann_intw
          weights = (abs(u_int(iw,:)))**2
          !write(io_unit_b,'(2i4,20e14.6)') ik, iw, eig_int(iw),weights
          write(io_unit_b,'(40e14.6)') dkpath(ik), eig_int(iw), weights
       end do
    end do

    !
    write(io_unit_b,*) '#'
    write(io_unit_b,*) '#Special k-points in the .bnd_int file are:'
    do ik=1,nkspecial
         write(io_unit_b,'(a,3f10.4,a,i4,e14.6)') '#', kspecial(:,ik), ' --> ', kspecial_indices(ik), dkpath(kspecial_indices(ik))
    end do
    write(io_unit_b,*) '#'
    !

    close(io_unit_b)

    write(*,*)' Bands interpolation finished and written to file ', bpath_file

!================================================================================
!   DOS over a fine grid
!   Parameters of DOS plot from namelist /DOS/
!================================================================================

    estep = (efin_dos-eini_dos)/real(ne_dos-1,dp)
    allocate (DOS(ne_dos), PDOS(ne_dos,num_wann_intw))
    call interpolated_DOS (nk1_dos, nk2_dos, nk3_dos, eini_dos, efin_dos, esmear_dos, ne_dos, DOS, PDOS)
    io_unit_d = find_free_unit()
    dos_file = trim(mesh_dir)//trim(prefix)//trim('.dos_int')
    open(unit=io_unit_d,file=dos_file,status='unknown')
    write(io_unit_d,'(A)') '# E-E_Fermi(eV)      DOS          PDOS(1:num_wann_intw)'

    num_elec = 0.0_dp
    do ie = 1,ne_dos
       ener = eini_dos+(ie-1)*estep
       write(io_unit_d,'(40e14.6)') ener-chemical_potential, DOS(ie), PDOS(ie,1:num_wann_intw)
       num_elec = num_elec +  DOS(ie)*fermi_dirac(ener-chemical_potential, ktsmear)
    end do

    close(io_unit_d)

    write(*,*)' DOS sum test:'
    write(*,*)'       DOS integral up to Fermi (trapeze) = ', num_elec

    write(*,*)' DOS interpolation finished and written to file ', dos_file

!==========================
! Finally make sure allocate_and_read_u_mesh works (I do not need this for interpolation)
!==========================
!!!     call allocate_and_read_u_mesh()


!================================================================================
!       clean up and finish
!================================================================================

    deallocate (eig_int,DOS)

    write(*,'(A)') '====================================================='
    write(*,'(A)') '|               end program interpolate             |'
    write(*,'(A)') '|        ---------------------------------          |'
    write(*,'(A)') '====================================================='


end program interpolate
