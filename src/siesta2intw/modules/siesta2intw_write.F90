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
module siesta2intw_write

  use precision, only: dp

  implicit none

  ! k points
  integer :: nkirr
  real(kind=dp), allocatable, dimension(:,:) :: kirr_cryst
  real(kind=dp), allocatable, dimension(:,:) :: kirr_cart

  ! Unit cell
  real(kind=dp) :: alat ! Lattice parameter
  real(kind=dp) :: tpiba ! 2pi/alat
  real(kind=dp), dimension(3,3) :: at ! Lattice vectors
  real(kind=dp), dimension(3,3) :: bg ! Reciprocal lattice vectors

  ! Wave functions
  integer :: nbnd ! Number of wfc's (for colinear calculations we doubled the number of bands instead of dealing with an additional spin index)
  integer :: nspin ! nspin = 1 for non-polarized calculations, otherwise nspin = 2
  integer :: spinor_comps ! spinor_comps = 2 for non-colinear calculations, otherwise spinor_comps = 1
  complex(kind=dp), allocatable, dimension(:,:,:,:) :: wf ! Wave functions in the localized orbitals basis readed from *.WFSX wf(no_u, spinor_comps, nbnd, nk)
  real(kind=dp), allocatable, dimension(:,:) :: et ! Eigenvalues readed from *.WFSX et(nbnd, nk)

  ! Radial Fourier transform interpolation tables of the localized orbital basis set
  real(dp), allocatable, dimension(:) :: gg_list
  real(dp), allocatable, dimension(:,:,:) :: phi_g_table
  real(dp), allocatable, dimension(:,:,:) :: phi_g_y2


contains

  subroutine setup_write_intw()
    !
    !
    !

    ! variables
    use precision, only: dp
    use units, only: pi
    use m_spin, only: spin
    use atomlist, only: no_u
    use siesta_geom, only: ucell
    use writewave, only: nwk, wfk, gamma_wavefunctions
    use parallel, only: ionode
    use siesta2intw_io, only: stdout, nbnd_initial, nbnd_final, use_sym, cutoff, &
                              nk1, nk2, nk3, kmesh, nkpoints, kpoints
    use siesta2intw_fft, only: gamma_only
    ! functions and subroutines
    use alloc, only: re_alloc
    use fdf, only: fdf_get, fdf_block, block_fdf
    use writewave, only: setup_wfs_list
    use siesta2intw_io, only: find_free_unit
    use siesta2intw_symmetry, only: compute_symmetry, print_symmetry_data
    use siesta2intw_fft, only: compute_fft_info
    !
    implicit none
    !
    type(block_fdf) :: bfdf
    integer :: max_nwf
    logical :: WaveFuncPresent
    integer :: ik

    if (ionode) then
      write(stdout, *) ""
      write(stdout, *) "***************************** BEGIN : siesta2intw *****************************"
      write(stdout, *) ""
    endif

    ! Get unit cell data
    alat = fdf_get( 'LatticeConstant', 0.0_dp, 'Bohr' )
    tpiba = 2.0_dp*pi/alat
    at = ucell/alat
    call reclat( at, bg, -1 )

    ! Set maximum number of eigenvalues
    if ( spin%Grid == 4 ) then
      max_nwf = 2 * no_u
    else
      max_nwf = no_u
    end if

    ! Set default band limits
    if (nbnd_initial .lt. 1) then
      nbnd_initial = 1
    endif
    if (nbnd_final .lt. 1) then
      nbnd_final = max_nwf
    endif
    ! Check limits
    if (nbnd_final .gt. max_nwf) then
      if (ionode) write(stdout, *) "nbnd_final is grether than the number of eigenvalues. Set it to the number of eigenvalues."
      nbnd_final = max_nwf
    endif
    if (nbnd_initial .gt. nbnd_final) then
      if (ionode) write(stdout, *) "nbnd_initial is grether than nbnd_final. Set it to nbnd_final."
      nbnd_initial = nbnd_final
    endif

    if (spin%Col) then
      nbnd = 2*(nbnd_final - nbnd_initial + 1) ! We doubled the bands instead of dealing with an exadditionaltra index
    else
      nbnd = nbnd_final - nbnd_initial + 1
    endif


    ! Get symmetry operations
    call compute_symmetry(use_sym)
    if (ionode) call print_symmetry_data()

    if (kmesh) then
      ! Compute irreducible k-points
      allocate(kirr_cryst(3, nk1*nk2*nk3))
      allocate(kirr_cart(3, nk1*nk2*nk3))
      kirr_cryst = 0.0_dp
      !
      call irreduciblek(nk1, nk2, nk3, nkirr, kirr_cryst)
      !
    else
      ! Use KPOINTS list
      allocate(kirr_cryst(3, nkpoints))
      allocate(kirr_cart(3, nkpoints))
      nkirr = nkpoints
      kirr_cryst = kpoints
      !
    endif
    !
    do ik = 1, nkirr
      kirr_cart(:, ik) = tpiba*matmul(bg, kirr_cryst(:, ik))
    enddo


    ! Configure siesta to compute and save the wave functions in the irreducible k-points

    ! Check if WaveFuncKPoints is specified in the fdf
    WaveFuncPresent = fdf_block('WaveFuncKPoints', bfdf)
    if (ionode .and. WaveFuncPresent) write(stdout, *) "WARNING: The WaveFuncKPoints block will be overwritten by siesta2intw"

    ! Allocate and set k points structure
    nwk = nkirr
    call re_alloc( wfk, 1, 3, 1, nwk, 'kpoint', 'setup_write_intw' )
    do ik = 1, nkirr
      wfk(:, ik) = kirr_cart(:, ik)
    enddo

    ! gamma_wavefunctions indicates to siesta to write the wave functions as real
    if ( kmesh ) then
      !
      if ( nk1 == 1 .and. nk2 == 1 .and. nk3 == 1 ) then
        gamma_wavefunctions = .true.
      else
        gamma_wavefunctions = .false.
      endif
      !
    else
      !
      if ( nkpoints == 1 .and. all(kpoints(:,1) < 0.000001_dp) ) then
        gamma_wavefunctions = .true.
      else
        gamma_wavefunctions = .false.
      endif
      !
    endif

    ! gamma_only indicates to siesta2intw to write the wave functions only for the positive G vectors
    gamma_only = gamma_wavefunctions .and. (spin%Grid /= 4)
    gamma_only = .false. ! For compatibility with the first version of intw (it does not support gamma_only wave functions)

    ! Allocate and set nwflist and iwf structures
    call setup_wfs_list(nwk, max_nwf, nbnd_initial, nbnd_final, .false., .false.)

    ! Allocate and set all FFT stuff (G vectors, list_iG...)
    if (ionode) call compute_fft_info(tpiba, bg, cutoff) ! there is no need to call it for non ionodes because it sets variables used only by write_intw_file

  end subroutine setup_write_intw


  subroutine write_intw_file()
    !
    !
    !

    ! variables
    use siesta_geom, only: na_u, xa, isa
    use siesta_options, only: g2cut
    use m_spin, only: TrSym, spin
    use atm_types, only: nspecies, species
    use writewave, only: nwk
    use m_ntm, only: ntm
    use siesta2intw_io, only: stdout, intwdir, outdir, prefix, phonons
    use siesta2intw_fft, only: nGk_max, list_iG, ngk, gvec, gvec_cart
    use siesta2intw_symmetry, only: nsym, s, ftau, t_rev
    ! functions and subroutines
    use siesta2intw_io, only: find_free_unit
    use siesta2intw_fft, only: write_fft_info
    use siesta2intw_pp, only: write_PP_files

    implicit none

    !
    ! Make the save.intw directory
    intwdir = trim(outdir)//trim(prefix)//".save.intw/"
    call execute_command_line("mkdir -p "//trim(intwdir))

    ! Clear save.intw directory
    call execute_command_line("rm -rf "//trim(intwdir)//"*")

    ! Write crystal.dat file
    call write_crystal_dat()

    ! Write all FFT stuff (G vectors, list_iG...)
    call write_fft_info()

    ! Write PPs
    call write_PP_files()

    ! Read wf and et from slabel.selected.WFSX
    call read_siesta_wfcs()

    ! Compute WFS in G space and write them
    call write_wfc()

    ! Deallocate wf and et
    deallocate(wf)
    deallocate(et)

    ! Write all phonon stuff (dynamical matrices, induced potentials...)
    if (phonons) call write_phonon_info()

    ! Deallocate all FFT stuff
    deallocate(list_iG, ngk)
    deallocate(gvec)
    deallocate(gvec_cart)

    write(stdout, *) ""
    write(stdout, *) "****************************** END : siesta2intw ******************************"
    write(stdout, *) ""

  end subroutine write_intw_file


  subroutine write_crystal_dat()
    !
    ! Write the crystal.dat file that contains all information about the system.
    !

    ! variables
    use siesta_geom, only: na_u, xa, isa
    use siesta_options, only: g2cut
    use m_spin, only: TrSym, spin
    use atm_types, only: nspecies, species
    use writewave, only: nwk
    use m_ntm, only: ntm
    use siesta2intw_io, only: stdout, intwdir, cutoff
    use siesta2intw_fft, only: nG, nGk_max, gamma_only
    use siesta2intw_symmetry, only: nsym, s, ftau, t_rev
    ! functions and subroutines
    use siesta2intw_io, only: find_free_unit

    implicit none

    character(len=256) :: datafile
    integer :: io_unit, i, isym


    ! open file
    io_unit = find_free_unit()
    datafile = trim(intwdir)//"crystal.dat"
    open(unit=io_unit, file=trim(datafile), status="unknown", action="write", form="formatted")

    ! Lattice constant
    write(unit=io_unit,fmt=*) "ALAT"
    write(unit=io_unit,fmt="(f16.10)") alat

    ! Lattice vectors
    write(unit=io_unit,fmt=*) "AT"
    do i=1,3
      write(unit=io_unit, fmt="(3f12.6)") at(i, 1:3)
    enddo

    ! Reciprocal lattice vectors
    write(unit=io_unit,fmt=*) "BG"
    do i=1,3
      write(unit=io_unit, fmt="(3f12.6)") bg(i, 1:3)
    enddo

    ! Number of species
    write(unit=io_unit,fmt=*) "NTYP"
    write(unit=io_unit,fmt=*) nspecies

    ! Atom labels, mass and PP file
    write(unit=io_unit,fmt=*) "ATOM_LABELS, MASS AND PP_FILE (1:NTYP)"
    do i=1,nspecies
      write(unit=io_unit,fmt="(a,f16.5,x,a)") trim(species(i)%symbol), species(i)%mass, trim(species(i)%label)
    enddo

    ! Number of atoms
    write(unit=io_unit,fmt=*) "NAT"
    write(unit=io_unit,fmt=*) na_u

    ! Atomic positions
    write(unit=io_unit,fmt=*) "POSITIONS (1:NAT)"
    do i=1,na_u
      write(unit=io_unit,fmt="(a4,i4,3f16.8)") trim(species(isa(i))%symbol), isa(i), xa(1:3, i)/alat
    end do

    ! Number of symmetries
    write(unit=io_unit,fmt=*) "NSYM"
    write(unit=io_unit,fmt=*) nsym

    ! Symmetry operations
    do isym=1,nsym
      write(unit=io_unit,fmt=*) "SYM"
      write(unit=io_unit, fmt="(i8)") isym
      do i = 1, 3
        write(unit=io_unit, fmt="(3i8)") s(i, 1:3, isym)
      enddo
      write(unit=io_unit, fmt="(3f16.10)")  ftau(1:3, isym)
      write(unit=io_unit, fmt="(i3)") t_rev(isym)
    enddo

    ! Spin-polarized calculation
    write(unit=io_unit,fmt=*) "LSPIN"
    write(unit=io_unit,fmt=*)  spin%Col .or. spin%NCol .or. spin%SO ! We transform the colinear calculations into non-colinear format

    ! Non-colinear spin polarized calculation with spin-orbit coupling
    write(unit=io_unit,fmt=*) "LSPINORB"
    write(unit=io_unit,fmt=*)  spin%SO

    ! This variable means if the induced potential is a 2x2 matirx or not.
    ! With SIESTA, except if the calculation is non-polarized, we always will write the full 2x2 matrix.
    write(unit=io_unit,fmt=*) "LMAG"
    write(unit=io_unit,fmt=*) spin%Col .or. spin%NCol .or. spin%SO
    !
    ! TODO: What happens if we have a non-magnetic colinear calculation?
    !       This type of calculation has no sense, as it would be equivalent to a non-polarized calculation.
    !       Maybe we could check if this is the case, and print an error or warning message.
    !
    ! TODO: What happens if we have a non-magnetic non-colinear calculation?
    !       If the calculation is not magnetic, the system has time reversal symmetry. By default, SIESTA
    !       does not use time reversal symmetry to reduce the k-points for non-colinear calculations.
    !       However, QE does use time reversal symmetry in non-magnetic non-colinear calculations, and
    !       therefore, only writes one component of the induced potential (rho), as the 2x2 matrix will
    !       be rho.I2 in such a case. This is what the domag variable indicates in QE.
    !       We should check and write in the documentation what does imply time reversal symmetry on the
    !       wave functions, on the Hamiltonian and on the induced potential.
    !

    ! Number of k-points
    write(unit=io_unit,fmt=*) "NKS"
    write(unit=io_unit,fmt=*) nwk

    ! Number of bands
    write(unit=io_unit,fmt=*) "NBAND"
    write(unit=io_unit,fmt=*) nbnd

    ! FFT grid
    write(unit=io_unit,fmt=*) "FFT GRID"
    write(unit=io_unit,fmt=*) ntm(1:3)

    ! Wave function cutoff
    write(unit=io_unit,fmt=*) "ECUTWFC"
    write(unit=io_unit,fmt="(f12.6)") cutoff

    ! Charge density cutoff
    if (cutoff .gt. g2cut) then
      write(stdout, *) "WARNING: The wave function cutoff is higher than the MeshCutoff used in the siesta calculation:"
      write(stdout, *) "WARNING: The G-vectors could be limited by the size of the FFT grid imposed by MeshCutoff."
    endif
    write(unit=io_unit,fmt=*) "ECUTRHO"
    write(unit=io_unit,fmt="(f12.6)") g2cut

    ! Number of G vectors
    write(unit=io_unit,fmt=*) "NG"
    write(unit=io_unit,fmt=*) nG

    ! Max number of G vectors for the wave function
    write(unit=io_unit,fmt=*) "NGK_MAX"
    write(unit=io_unit,fmt=*) nGk_max

    close(unit=io_unit)

  end subroutine write_crystal_dat


  subroutine write_phonon_info()
    !
    !
    !

    ! variables
    use siesta_geom, only: na_u, isa, xa
    use m_ntm, only: ntm
    use atm_types, only: nspecies, species
    use siesta2intw_io, only: stdout, intwdir, outdir, prefix, phdir, dvscfdir, &
                              nqirr
    use siesta2intw_utils, only: cmplx_0
    ! functions and subroutines
    use siesta2intw_io, only: find_free_unit
    use siesta2intw_utils, only: ainv
    use siesta2intw_ph, only: allocate_vloc, deallocate_vloc, phq_init, calculate_local_part_dv

    implicit none

    integer :: iq
    character(len=6) :: iq_str

    integer :: io_unit, io_unit_read, io_unit_write, ios, rl
    integer :: imode, jmode, ispin, i, j, ia, ja
    complex(kind=dp) :: u_irr(3*na_u)
    character(len=256) :: filename_read, filename_write

    character(len=2) :: dv_precision
    integer :: nspin_dv, nq_star
    complex(kind=dp), allocatable :: dv_local(:), dvq(:, :)
    real(kind=dp), allocatable :: q_irr_cryst(:, :), q_irr_cart(:, :)
    real(kind=dp) :: qpoint(3)
    real(kind=dp) :: dynq_re(3,3), dynq_im(3,3)

    ! Read the dynamical matrix and q-points
    allocate(q_irr_cryst(3, nqirr))
    allocate(q_irr_cart(3, nqirr))
    q_irr_cryst = 0.0_dp
    q_irr_cart = 0.0_dp
    do iq = 1, nqirr
      !
      write(filename_read,"(a3,i4.4,a4)") "dyn", iq, ".dat"
      io_unit_read = find_free_unit()
      open( unit=io_unit_read, file=trim(outdir)//trim(phdir)//trim(filename_read), iostat=ios, &
            action="read", status="old" )
      if (ios .ne. 0) stop "ERROR: write_phonon_info: Error opening dynamical matrix file."

      if (                iq <   10) write(iq_str,"(a2,i1)") "_q", iq
      if ( 10 <= iq .and. iq <  100) write(iq_str,"(a2,i2)") "_q", iq
      if (100 <= iq .and. iq < 1000) write(iq_str,"(a2,i3)") "_q", iq

      filename_write = trim(prefix)//".dyn"//trim(iq_str)
      io_unit_write = find_free_unit()
      open( unit=io_unit_write, file=trim(intwdir)//trim(filename_write), iostat=ios, &
            action="write", status="replace" )
      if (ios .ne. 0) stop "ERROR: write_phonon_info: Error opening intw's dynamical matrix file."


      ! Number of q-points in the star
      read(io_unit_read, "(10x,i3)") nq_star

      ! Irreducible q-point
      read(io_unit_read, "(10x,3f10.6)") qpoint
      write(io_unit_write,"(a11,3f14.9,a2)") "q_cryst = ( ", qpoint, " )"
      q_irr_cryst(:, iq) = qpoint
      q_irr_cart(:, iq) = matmul(bg, qpoint)

      ! Dynamical matrix
      do ia=1,na_u
        do ja=1,na_u
          !
          read(io_unit_read, *) !ia, ja
          !
          do i=1,3
            read(io_unit_read, "(3(2f12.8),2x)") (dynq_re(i,j), dynq_im(i,j), j=1,3)
          enddo
          !
          do i=1,3
            write(io_unit_write, "(3(a,f16.10,a,f16.10,a))") ("(", dynq_re(i,j), ",", dynq_im(i,j), ") ", j=1,3)
          enddo
          !
        enddo
      enddo
      !
    enddo ! iq
    !
    close(unit=io_unit_read)
    close(unit=io_unit_write)

    ! Write irreducible q-points list
    io_unit = find_free_unit()
    filename_write = "qlist.txt"
    open( unit=io_unit, file=trim(outdir)//trim(filename_write), iostat=ios, form='formatted', &
          status='replace', action='write' )
    if (ios .ne. 0) stop "ERROR: write_phonon_info: Error opening qlist.txt file."
    !
    do iq = 1, nqirr
      write(io_unit, "(i4,3f18.12)") iq, q_irr_cart(:, iq) ! Remember that qlist_file is in Cartesian coordinates in units of 2pi/alat!!
    enddo
    !
    close(unit=io_unit)


    ! Write irreducible patterns
    io_unit = find_free_unit()
    filename_write = "irrq_patterns.dat"
    open(unit=io_unit, file=trim(intwdir)//trim(filename_write), status="replace", form="formatted", iostat=ios)
    if (ios /= 0) stop "ERROR: write_phonon_info: Error opening irrq_patterns.dat file."
    !
    do iq = 1, nqirr
      !
      write(unit=io_unit, fmt="(a,i4)") "q", iq
      do imode=1,3*na_u
        u_irr = cmplx_0
        u_irr(imode) = 1.0_dp
        write(unit=io_unit,fmt="(10000(a,f16.10,a,f16.10,a))") &
          ( "(", real(u_irr(jmode),dp), ",", aimag(u_irr(jmode)), ") ",jmode=1,3*na_u)
      enddo
      !
    enddo
    close(unit=io_unit)

    !
    ! Open the dvscf.info file
    io_unit = find_free_unit()
    filename_read = trim(dvscfdir)//"dvscf.info"
    open( unit=io_unit, file=trim(outdir)//trim(phdir)//trim(filename_read), iostat=ios, form='formatted', &
          status='old', action='read' )
    if (ios .ne. 0) stop "ERROR: write_phonon_info: Error opening dvscf.info file."
    !
    ! Read the data
    read(io_unit,*) !"DV_PRECISION"
    read(io_unit,*) dv_precision ! sp or dp for the potential
    if (dv_precision == "sp") stop "ERROR: write_phonon_info: dv_precision not implemented yet."
    !
    read(io_unit,*) !"NAT"
    read(io_unit,*)
    !
    read(io_unit,*) !"NIRRED_ATM"
    read(io_unit,*)
    !
    read(io_unit,*) !"NSPIN"
    read(io_unit,*) nspin_dv
    !
    close(unit=io_unit)

    !
    ! Write induced potential
    !
    ! TODO: Even if intw is not prepared to work with colinear calculations at this moment, the best option for a
    !       future implementation will be to write the induced potential in the same way as QE:
    !         - Non-polarized: (rho)
    !         - Colinear Magnetic (has to be magnetic): (rho_11, rho_22)
    !         - Non-Colinear Non-Magnetic (rho)
    !         - Non-Colinear Magnetic: (rho, m_x, m_y, m_z)
    !
    !       Meanwhile, we will transform the induced potential of colinear calculations to an equivalent non-colinear form:
    !        - Diagonal 2x2 matrix: (rho, 0, 0, m_z).
    !
    allocate(dvq(ntm(1)*ntm(2)*ntm(3), nspin_dv))
    dvq = cmplx_0
    allocate(dv_local(ntm(1)*ntm(2)*ntm(3)))
    dv_local = cmplx_0
    !
    inquire(iolength=rl) dvq(:, :)
    !
    call allocate_vloc()
    !
    do iq = 1, nqirr
      !
      qpoint = q_irr_cryst(:, iq)
      write(stdout, "(a5,i0,a6,3(f15.8))") "iq = ", iq, ', q = ', qpoint
      !
      call phq_init(tpiba, bg, qpoint)
      !
      write(filename_read,"(a5,i4.4,a4)") "dvscf", iq, ".dat"
      io_unit_read = find_free_unit()
      open( unit=io_unit_read, file=trim(outdir)//trim(phdir)//trim(dvscfdir)//trim(filename_read), iostat=ios, &
            form='unformatted', status='old', action='read', access='direct', recl=rl )
      if (ios .ne. 0) stop "ERROR: write_phonon_info: Error opening dv_file_read file."
      !
      if (                iq <   10) write(iq_str,"(a2,i1)") "_q", iq
      if ( 10 <= iq .and. iq <  100) write(iq_str,"(a2,i2)") "_q", iq
      if (100 <= iq .and. iq < 1000) write(iq_str,"(a2,i3)") "_q", iq
      filename_write = trim(prefix)//".dvscf"//trim(iq_str)
      io_unit_write = find_free_unit()
      open( unit=io_unit_write, file=trim(intwdir)//trim(filename_write), iostat=ios, form='unformatted', &
            status='replace', action='write', access='direct', recl=rl )
      if (ios .ne. 0) stop "ERROR: write_phonon_info: Error opening dv_file_write file."

      do imode = 1, 3*na_u
        !
        read(unit=io_unit_read, rec=imode, iostat=ios) dvq
        !
        call calculate_local_part_dv(tpiba, bg, qpoint, imode, dv_local)
        !
        dvq(:, 1) = dvq(:, 1) - dv_local ! the local part only affects to rho
        !
        if (nspin_dv == 2) then
          ! Transform the colinear potential into a diagonal 2x2 non-colinear potential for compatibility with intw
          write(unit=io_unit_write, rec=2*(imode-1)+1, iostat=ios) dvq(:, 1), 0.0_dp*dvq(:, 1)
          write(unit=io_unit_write, rec=2*(imode-1)+2, iostat=ios) 0.0_dp*dvq(:, 2), dvq(:, 2)
        else
          write(unit=io_unit_write, rec=imode, iostat=ios) dvq
        endif
        !
      enddo ! imode
      !
      close(io_unit_read)
      close(io_unit_write)
      !
    enddo ! iq
    !
    call deallocate_vloc()

  end subroutine write_phonon_info


  subroutine read_siesta_wfcs()
    !
    ! Read wave functions (wf) and eigenvalues (et) for all k points from the output SIESTA
    ! file (slabel.selected.WFSX)
    !

    ! variables
    use precision, only: sp
    use files, only: slabel
    use atomlist, only: no_u
    use writewave, only: gamma_wavefunctions
    use m_spin, only: spin
    use siesta2intw_io, only: stdout
    use siesta2intw_utils, only: cmplx_0
    ! functions and subroutines
    use siesta2intw_io, only: find_free_unit

    implicit none

    real(kind=sp), allocatable :: wf_single(:,:)
    integer :: number_of_wfns ! Number of wfc's for each k-point
    integer :: nspin_blocks ! nspin_blocks = 2 for colinear calculations, otherwise nspin_blocks = 1
    logical :: non_coll
    real(kind=dp), dimension(3) :: kpoint
    integer :: ik, ispin, iwf, iwf_orig, ibnd
    character (len=256) :: wfs_filename
    real(kind=dp) :: ener
    integer :: io, ios, idummy
    integer :: wf_unit
    integer :: no_u_wfsx, nspin_wfsx
    logical :: gamma_wfsx
    integer :: nk



    write(stdout, *) "- Reading wfcs..."
    !
    wfs_filename =trim(slabel)//".selected.WFSX"
    wf_unit=find_free_unit()
    open(wf_unit, file=trim(wfs_filename), form="unformatted", status='old', iostat=ios, position='rewind')
    if (ios /= 0) stop "ERROR: read_siesta_wfcs: error opening slabel.selected.WFSX file."
    !
    read(wf_unit) nk, gamma_wfsx
    if (nkirr /= nk) stop "ERROR: read_siesta_wfcs: nkirr /= nk"
    if (gamma_wfsx .neqv. gamma_wavefunctions) stop "ERROR: read_siesta_wfcs: gamma_wfsx /= gamma_wavefunctions"
    !
    read(wf_unit) nspin_wfsx
    if (nspin_wfsx /= spin%H) stop "ERROR: read_siesta_wfcs: nspin_wfsx /= spin%H"
    !
    if (nspin_wfsx >= 2) then
      nspin = 2
    else
      nspin = 1
    endif
    !
    ! Non-collinear or SOC files have a single "spin" block as opposed to collinear-spin
    ! files, which contain two spin blocks per k section.
    non_coll = (nspin_wfsx >= 4)
    if (non_coll) then
      nspin_blocks = 1
      spinor_comps = 2
    else
      nspin_blocks = nspin_wfsx
      spinor_comps = 1
    endif
    !
    read(wf_unit) no_u_wfsx
    if (no_u /= no_u_wfsx) then
      call die("Mismatch in no_u in WFSX and DIM/PLD files")
    endif
    read(wf_unit)   ! orbital labels
    ! write(iu) (iaorb(j),labelfis(isa(iaorb(j))),
    !  .            iphorb(j), cnfigfio(isa(iaorb(j)),iphorb(j)),
    !  .            symfio(isa(iaorb(j)),iphorb(j)), j=1,nuotot)

    if (allocated(et)) deallocate(et)
    allocate(et(nbnd, nk))
    et = 0.0_dp
    !
    if (allocated(wf)) deallocate (wf)
    allocate(wf(no_u, spinor_comps, nbnd, nk))
    wf = cmplx_0
    !
    if (allocated(wf_single)) deallocate (wf_single)
    if (non_coll) then
      allocate(wf_single(4, no_u))
    else
      if (gamma_wfsx) then
        allocate(wf_single(1, no_u))
      else
        allocate(wf_single(2, no_u))
      endif
    endif
    !
    !
    do ik=1,nk
      !
      do ispin=1,nspin_blocks
        !
        read(wf_unit) idummy, kpoint(1:3)
        if (idummy /= ik) then
          write(stdout, *) "ik index mismatch in WFS file"
          write(stdout, *) "ik in file, ik: ", idummy, ik
          stop
        endif
        if ( sum(abs(kirr_cart(:, ik) - kpoint)) > 0.0001_dp ) then
          write(stdout, *) "k point mismatch in WFS file"
          write(stdout, "(a,3f15.10)") "k in file: ", kpoint
          write(stdout, "(a,3f15.10)") "k: ", kirr_cart
          stop
        endif

        !
        read(wf_unit) idummy
        if (idummy /= ispin) then
          write(stdout, *) "ispin index mismatch in WFS file"
          write(stdout, *) "ispin in file, ispin: ", idummy, ispin
          stop
        endif
        !
        read(wf_unit) number_of_wfns
        if (number_of_wfns /= nbnd/nspin_blocks) then
          write(stdout, *) "number_of_wfns mismatch in WFS file"
          write(stdout, *) "nbnd in file, nbnd: ", number_of_wfns, nbnd/nspin_blocks
          stop
        endif
        !
        write(stdout, "(a,i0,a,i0)") "Processing kpoint ", ik, " / ", nk
        write(stdout, "(a,3f15.10)") "k = ", kpoint
        write(stdout, *) "--------------------------------"
        !
        do iwf=1,number_of_wfns
          !
          ibnd = nspin_blocks*(iwf-1) + ispin ! We doubled nbnd for colinear calculations to dont deal with an additional index
          !
          read(wf_unit) iwf_orig
          ! Note that we mean the *original* index
          !
          if (iwf_orig /= iwf) then
            ! The file holds a subset of wfs, with the original indexes...
            write(stdout, "(a6,i0,a3,i0,a13,i0,a1)") "iwf = ", iwf, " / ", number_of_wfns, &
                                                     " (iwf_orig = ", iwf_orig, ")"
          else
            write(stdout, "(a6,i0,a3,i0)") "iwf = ", iwf, " / ", number_of_wfns
          endif
          !
          read(wf_unit) ener
          et(ibnd,ik) = ener
          !
          read(wf_unit) (wf_single(:,io), io=1,no_u)
          ! Use a double precision complex form in what follows
          if (non_coll) then
            wf(:,1,ibnd,ik) = cmplx(wf_single(1,:), wf_single(2,:), kind=dp)
            wf(:,2,ibnd,ik) = cmplx(wf_single(3,:), wf_single(4,:), kind=dp)
          else
            if (gamma_wfsx) then
              wf(:,1,ibnd,ik) = cmplx(wf_single(1,:), 0.0_sp, kind=dp)
            else
              wf(:,1,ibnd,ik) = cmplx(wf_single(1,:), wf_single(2,:), kind=dp)
            endif
          endif
          !
        enddo ! wfn
        !
      enddo ! spin_block
      !
    enddo ! k-point

    deallocate(wf_single)

  end subroutine read_siesta_wfcs


  subroutine write_wfc()
    !
    !
    !

    ! variables
    use precision, only: dp
    use units, only: pi
    use siesta_geom, only: xa, isa, volume_of_some_cell
    use writewave, only: nwk
    use atomlist, only: indxuo, iaorb, no_u, iphorb
    use atmparams, only: lmaxd
    use m_spin, only: spin
    use atm_types, only: nspecies, species
    use siesta2intw_io, only: stdout, intwdir
    use siesta2intw_utils, only: cmplx_0, cmplx_i
    use siesta2intw_fft, only: nGk_max, list_iG, ngk, gvec_cart
#ifdef DEBUG
    use siesta2intw_fft, only: gamma_only
#endif
    ! functions and subroutines
    use spher_harm, only: rlylm
    use siesta2intw_io, only: find_free_unit

    implicit none

    integer :: ik, ibnd, jbnd, is, ig, ispin, ispinor
    character(256) :: wfc_file, datafile
    integer :: io_unit
    integer :: io, ia
    real(dp) :: kpoint(3), kr, kpg, kg_unit(1:3)
    integer :: ipol
    complex(dp) :: phase
    integer :: l,m, ind
    integer :: io_l, io_u
    real(kind=dp), allocatable :: rylm(:), grylm(:,:)
    type orbitals
      complex(kind=dp), allocatable :: o(:,:,:)
    end type orbitals
    type(orbitals), allocatable :: orbital_gk(:)

    complex(kind=dp), allocatable :: evc(:)

#ifdef DEBUG
    complex(kind=dp) :: ss
    ! Variables to print the orbitals or wave functions in real space
    integer, dimension(3) :: nri
    integer :: ir1, ir2, ir3
    real(kind=dp) :: r1, r2, r3, gr
    complex(kind=dp) :: or
#endif


    write(stdout, *) "- Writing k-points..."

    !
    ! Write k points
    datafile=trim(intwdir)//"kpoints.dat"
    io_unit=find_free_unit()
    open(unit=io_unit,file=datafile,status="unknown", action="write",form="formatted")
    do ik=1,nwk
     write(io_unit,"(3f18.10)") kirr_cart(1:3, ik) / tpiba
    end do
    close(unit=io_unit)


    write(stdout, *) "- Computing basis orbital's FFT..."

    call build_basis_g_table()

    ! First, compute the spherical Fourier transform for the basis set orbitals.
    ! Then, for each atomic possition xa(:,ia) we just have to add a phase.

    allocate(rylm( 0:(lmaxd+1)**2))
    allocate(grylm(1:3, 0:(lmaxd+1)**2-1))
    allocate(orbital_gk(nspecies))
    do is=1,nspecies
      allocate(orbital_gk(is)%o(nGk_max, species(is)%norbs, nwk))
      orbital_gk(is)%o = cmplx_0
    enddo
    rylm = 0.0_dp
    grylm = 0.0_dp
    !
    do ik = 1, nwk
      !
      write(stdout, "(a5,i0,a8,i0)") "ik = ", ik, ", ngk = ", ngk(ik)
      kpoint(1:3) = kirr_cart(1:3, ik)
      !
      do iG = 1, ngk(ik)
        !
        kpg = sqrt(sum( ( kpoint(:) + gvec_cart(:, list_iG(ig, ik)) )**2 )) ! |k+G| modulus.
        kg_unit = - (kpoint(:) + gvec_cart(:, list_iG(ig, ik))) ! (k+G)/|k+G| unit vector
        if (kpg>1.0E-9) kg_unit = kg_unit/kpg
        call rlylm( lmaxd, kg_unit, rylm, grylm ) ! compute all Ylm up to lamx
        !
        do is = 1, nspecies
          !
          do io = 1, species(is)%norbs
            !
            l = species(is)%orb_l(io) ! l of the orbital
            m = species(is)%orb_m(io) ! m of the orbital
            !
            ind = l**2 + l + m !+ 1 ! unified lm index
            phase = cmplx_i**l ! i^l factor exp(kr) = 4pi \sum_{lm} i^l j_l(kr) Y_lm(^k) Y_lm(^r)
            !
            orbital_gk(is)%o(iG, io, ik) = 4*pi * phase * rylm(ind) * spline_g_of_phi(is, io, kpg)
            !
          end do ! iG
          !
        end do ! io
        !
      end do ! is
      !
    end do ! ik
    !
    call deallocate_basis_g_table()


#ifdef DEBUG
    !
    ! Print the basis orbitals in real space
    !
    ik=1
    !
    nri(1) = 500
    !
    do is=1,nspecies
      do io=1,species(is)%norbs
        !
        l = species(is)%orb_l(io)    ! l of the orbital
        m = species(is)%orb_m(io)    ! m of the orbital
        ind = l**2 + l + m !+ 1     ! unified lm index
        write(stdout,*) "print orb:", is, io, l, m
        write(is*10000+io,"(a)") "#             r            mod           real           imag"
        !
        do ir1=1,nri(1)
          r1 = 0.0001_dp+1.0_dp*(ir1-1)/(nri(1)-1)*alat
          r2 = 0.0001_dp+1.0_dp*(ir1-1)/(nri(1)-1)*alat
          r3 = 0.0001_dp+1.0_dp*(ir1-1)/(nri(1)-1)*alat
          !
          or = cmplx_0
          do ig=1,ngk(ik)
            gr = (kirr_cart(1,ik)+gvec_cart(1,list_iG(ig,ik)))*r1 + &
                 (kirr_cart(2,ik)+gvec_cart(2,list_iG(ig,ik)))*r2 +&
                 (kirr_cart(3,ik)+gvec_cart(3,list_iG(ig,ik)))*r3
            or = or + orbital_gk(is)%o(iG,io,ik)*exp(cmplx_i*gr)
          enddo
          if (gamma_only) then
            or = 2.0_dp*real(or, dp) - orbital_gk(is)%o(1,io,ik)
          endif
          !
          call rlylm( l, (/r1,r2,r3/), rylm, grylm ) ! compute all Ylm up to lamx
          !
          or = or/volume_of_some_cell
          if (rylm(ind).gt.0.000001_dp) or = or/rylm(ind)
          !
          write(is*10000+io,"(4f15.10)") sqrt(sum((/r1,r2,r3/)**2)), abs(or), real(or), aimag(or)
          !
        enddo ! ir1
        !
      enddo ! io
    enddo ! is
    !
    stop "print orb"
    !
#endif

    deallocate(rylm, grylm)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This is the old way, using vectorization on the ibnd index. The problem
    ! with this approach is that all the wave functions are calculated at the
    ! same time, thus, a variable to store all Fourier components of all the
    ! wave functions at the same time is needed, implying a high RAM consumption
    ! in big systems.
    ! In all different tests that I have done, this is the fastest way to do the
    ! calculation. And it works fine in all situation (non-polarized, collinear
    ! and non-collinear).
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! allocate(evc(nspin*nGk_max, nbnd, nwk))
    ! evc = cmplx_0
    !
    ! ! Use Fourier transform of the orbitals to compute wave function Fourier transform
    ! call cpu_time(t0)
    ! do ik=1,nwk
    !   kpoint(1:3) = kirr_cart(1:3, ik)
    !   write(stdout,*)"ik", ik, ngk(ik)
    !   do iG=1,ngk(ik)
    !     do ipol=1,nspin
    !       do io=1,no_u
    !         ia = iaorb(io) ! Index of the atom to wich this orbital belongs
    !         is = isa(ia) ! Index of the atomic species to wich this orbital belongs
    !         io_l = iphorb(io) ! Orbital index of each orbital in its atom
    !         io_u = indxuo(io) ! Index of equivalent orbital in "u" cell
    !         kr = sum( (- gvec_cart(:,list_iG(ig,ik))) * xa(:,ia)) ! xa-k is multiplied by alat, and kpoit by 2pi/alat
    !         phase = cmplx(cos(kr),sin(kr))
    !         evc((ipol-1)*ngk(ik)+ig, :,ik) = evc((ipol-1)*ngk(ik)+ig, :,ik) + &
    !             phase * orbital_gk(is)%o(iG,io_l,ik) * wf(io_u, ipol, 1:nbnd, ik) / sqrt(volume_of_some_cell)
    !       end do !io
    !     end do !ipol
    !   end do !iG
    ! end do !ik
    ! call cpu_time(t1)
    ! write(stdout,*) t1-t0
    !
    ! ! Save the Fourier transform of the wave function
    ! do ik=1,nwk
    !   do ispin=1,nspin_blocks
    !     if (nspin_blocks==1) then
    !       write(wfc_file,"(a4,i5.5,a4)") "wfc_",ik,".dat"
    !     elseif (nspin_blocks==2) then
    !       if (ispin==1) write(wfc_file,"(a5,i5.5,a4)") "wfc1_",ik,".dat"
    !       if (ispin==2) write(wfc_file,"(a5,i5.5,a4)") "wfc2_",ik,".dat"
    !     endif
    !     datafile=trim(trim(intwdir)//trim(wfc_file))
    !     io_unit=find_free_unit()
    !     open(unit=io_unit,file=datafile,status="unknown", action="write",form="unformatted")
    !     write(unit=io_unit) ngk(ik)
    !     write(unit=io_unit) list_iG(1:ngk(ik), ik)
    !     write(unit=io_unit) ( et(ibnd,ispin,ik), ibnd=1,nbnd )
    !     do ibnd=1,nbnd
    !       if (nspin_blocks==1) then
    !         write(unit=io_unit) ( ( evc( (ipol-1)*ngk(ik)+ig, ibnd,ik), ig=1,ngk(ik) ), ipol=1,nspin )
    !       elseif (nspin_blocks==2) then
    !         write(unit=io_unit) ( evc( (ispin-1)*ngk(ik)+ig, ibnd,ik), ig=1,ngk(ik) )
    !       endif
    !     enddo ! ibnd
    !     close(unit=io_unit)
    !   enddo ! ispin
    ! enddo ! ik
    !
    ! ! Orthogonality test
    ! ispin = 1
    ! do ik=1,nwk
    !   write(stdout,*) "ik", ik
    !   do ibnd=1,nbnd
    !     do jbnd=1,nbnd
    !       ss=(0.0,0.0)
    !       do ig=1, ngk(ik)
    !         if (nspin_blocks==1) then
    !           do ipol=1,nspin
    !             ss=ss+ conjg(evc( (ipol-1)*ngk(ik)+ig, ibnd,ik)) * evc( (ipol-1)*ngk(ik)+ig, jbnd,ik)
    !           enddo
    !         elseif (nspin_blocks==2) then
    !           ss=ss+ conjg(evc( (ispin-1)*ngk(ik)+ig, ibnd,ik)) * evc( (ispin-1)*ngk(ik)+ig, jbnd,ik)
    !         endif
    !       end do
    !       write(*,"(a,2i3,2f14.8)") "ibnd jbnd", ibnd, jbnd, ss
    !     enddo ! jbnd
    !   enddo ! ibnd
    ! enddo ! ik


    write(stdout, *) "- Writing wfc's FFT..."

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This is the new approach that I implemented, using vectorization over
    ! the iG index. In this way I can calculate all the Fourier transform
    ! components of a wave function at a time. Thus, I can save the Fourier
    ! transform of each eigenstate at a time, and I don't need a variable to
    ! store all Fourier components of all the wave functions at the same time.
    ! On the one hand, this allows to maintain a low RAM consumption even in big
    ! systems. But on the other hand, it is slower than the original approach.
    ! In the tests that I have done I found this approach ~6 times slower than
    ! the original one.
    ! I haven't done any test for the non-collinear case, so this should be
    ! tested.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate(evc(spinor_comps*nGk_max))

    ! Use Fourier transform of the orbitals to compute wave function Fourier transform
    do ik=1,nwk
      !
      kpoint(1:3) = kirr_cart(1:3, ik)
      write(stdout, "(a5,i0,a8,i0)") "ik = ", ik, ", ngk = ", ngk(ik)
      !
      write(wfc_file, "(a3,i5.5,a4)") "wfc", ik, ".dat"
      datafile = trim(intwdir)//trim(wfc_file)
      io_unit = find_free_unit()
      open(unit=io_unit, file=datafile, status="unknown", action="write", form="unformatted")

      write(unit=io_unit) ngk(ik)
      write(unit=io_unit) list_iG(1:ngk(ik), ik)
      write(unit=io_unit) ( et(ibnd,ik), ibnd=1,nbnd )
      !
      do ibnd=1,nbnd
        !
        write(stdout, "(a6,i0,a3,i0)") "iwf = ", ibnd, " / ", nbnd
        !
        evc = cmplx_0
        !
        do io=1,no_u
          !
          ia = iaorb(io) ! Index of the atom to wich this orbital belongs
          is = isa(ia) ! Index of the atomic species to wich this orbital belongs
          io_l = iphorb(io) ! Orbital index of each orbital in its atom
          io_u = indxuo(io) ! Index of equivalent orbital in "u" cell
          !
          do ispinor = 1, spinor_comps
            !
            evc(ngk(ik)*(ispinor-1)+1:ngk(ik)*ispinor) = evc(ngk(ik)*(ispinor-1)+1:ngk(ik)*ispinor) + &
              exp(-cmplx_i * (  gvec_cart( 1, list_iG(1:ngk(ik), ik) )*xa(1, ia) &
              + gvec_cart( 2, list_iG(1:ngk(ik), ik) )*xa(2, ia) &
              + gvec_cart( 3, list_iG(1:ngk(ik), ik) )*xa(3, ia) ) ) &
              * orbital_gk(is)%o(1:ngk(ik), io_l, ik) &
              * wf(io_u, ispinor, ibnd, ik) / sqrt(volume_of_some_cell)
            !
          enddo ! ispinor
          !
        enddo ! io
        !
        if (spin%Col) then
          if (mod(ibnd, 2) /= 0) then ! odd bands are for spin up
            write(unit=io_unit) evc(1:ngk(ik)), 0.0_dp*evc(1:ngk(ik))
          else ! even bands for spin down
            write(unit=io_unit) 0.0_dp*evc(1:ngk(ik)), evc(1:ngk(ik))
          endif
        else
          write(unit=io_unit) ( ( evc( (ispinor-1)*ngk(ik)+ig), ig=1,ngk(ik) ), ispinor=1,spinor_comps )
        endif
        !
#ifdef DEBUG
    !
    ! Print the wave functions in real space
    !
    nri(1) = 120
    nri(2) = 120
    !
    do ispinor=1,spinor_comps
      !
      write(stdout,*) "print wfc:", "ik=", ik, "ibnd=", ibnd, "ispin=", ispin
      ! write(ik*10000+ibnd,"(a)") "#          x           y            mod           real           imag"
      !
      do ir1=1,nri(1)
        r1 = 1.0_dp*(ir1-1)/(nri(1)-1)*alat
        do ir2=1,nri(2)
          r2 = 1.0_dp*(ir2-1)/(nri(2)-1)*alat
            !
            ! 2D plot for z=0
            !
            or = cmplx_0
            do ig=1,ngk(ik)
              gr = (kirr_cart(1,ik)+gvec_cart(1,list_iG(ig,ik)))*r1 + &
                   (kirr_cart(2,ik)+gvec_cart(2,list_iG(ig,ik)))*r2
              or = or + evc((ispinor-1)*ngk(ik)+ig)*exp(cmplx_i*gr)
            enddo
            if (gamma_only) then
              or = 2.0_dp*real(or, dp)-evc((ispinor-1)*ngk(ik)+1)
            endif
            or = or/sqrt(volume_of_some_cell)
            !
            write(ik*10000+ibnd,"(2f12.5,2x,f10.5)") r1, r2, abs(or)
            !
          enddo
          write(ik*10000+ibnd,*) ""
      enddo
      !
    enddo
    !
    ! stop "print wfc"
    !
#endif
        !
      enddo ! ibnd
      !
      close(io_unit)
      !
    enddo ! ik


    deallocate(orbital_gk)

  end subroutine write_wfc


  subroutine deallocate_basis_g_table()

    deallocate(phi_g_table)
    deallocate(phi_g_y2)
    deallocate(gg_list)

  end subroutine deallocate_basis_g_table


  subroutine build_basis_g_table()
    !
    !
    !

    ! variables
    use siesta_options, only: g2cut
    use atm_types, only: species, nspecies
    ! functions and subroutines
    use siesta2intw_utils, only: spline

    implicit none

    integer, parameter :: nq = 500
    integer :: is, io, iq

    ! allocate with the dimensions of the atomic specie with the maximum number of orbitals
    allocate(phi_g_table(nspecies, maxval(species(:)%norbs), nq))
    allocate(phi_g_y2   (nspecies, maxval(species(:)%norbs), nq))
    allocate(gg_list(nq))
    phi_g_table = 0.0_dp
    phi_g_y2 = 0.0_dp
    gg_list = 0.0_dp

    ! ASIER: Here the Fourier transform for different m are the same,
    !        but it is computed anyway for simplicity of the code

    do is=1,nspecies
      !
       do io=1, species(is)%norbs
        !
          ! write(stdout,*) " -- io --", io, species(is)%orb_l(io), species(is)%orb_m(io)
          ! write(stdout,*) "rc", rcut(is,io)
          do iq=1,nq
             gg_list(iq) = (iq-1)*sqrt(g2cut)/nq
             phi_g_table(is,io,iq) = g_of_phi(is,io,gg_list(iq))
          enddo ! iq
          !
          call spline(gg_list, phi_g_table(is,io,:), 0.0_dp, 0.0_dp, phi_g_y2(is,io,:))
          !
       end do ! io
       !
    end do ! is

  end subroutine build_basis_g_table


  function spline_g_of_phi(is, io, q)
    !
    !
    !

    ! functions and subroutines
    use siesta2intw_utils, only: splint

    implicit none

    integer, intent(in) :: is, io ! specie and orbital index
    real(kind=dp), intent(in) :: q
    real(kind=dp) :: spline_g_of_phi


    if (q > maxval(gg_list)) then
       spline_g_of_phi = 0.0_dp
       return
    end if
    !
    spline_g_of_phi = splint(gg_list, phi_g_table(is,io,:), phi_g_y2(is,io,:), q)

  end function spline_g_of_phi


  function g_of_phi(is,io,gg)
    !
    ! Runge-Kutta4 (=Simpson) integration \int^rc_0 r**2 * bessel(l,k*r) * phi(is,io,r) dr
    !

    use atm_types, only: species
    use atmfuncs, only: rcut, rphiatm

    implicit none

    integer, intent(in) :: is, io ! specie and orbital index
    real(kind=dp), intent(in) :: gg ! module of g+k vector
    real(kind=dp) :: g_of_phi
    !-lokal
    integer :: ir
    integer,parameter :: nr = 500 ! 500 points to integrate [0,rc]. I don't know if it's worth generalizing
    real(kind=dp) :: r, k1, k2, phi, dphidr, h, rh1, rh2
    integer :: l
    real(kind=dp), parameter :: sqrt3 = 0.577350269189626_dp

    l = species(is)%orb_l(io)

    h = rcut(is,io)/nr

    g_of_phi = 0.0_dp
    r = 0.0_dp
    ! Gauss quadrature n=2 step by step
    do ir=1,nr
       r = (ir-1)*h ! (nr-1)*rc/nr=rc-h
       rh1 = r + h/2*( 1.0_dp - sqrt3 )
       call rphiatm(is,io,rh1,phi,dphidr)
       k1 = phi * rh1**2 * sphb(l,rh1*gg)
       rh2= r + h/2*( 1.0_dp + sqrt3 )
       call rphiatm(is,io,rh2,phi,dphidr)
       k2 = phi * rh2**2 * sphb(l,rh2*gg)
       g_of_phi = g_of_phi + h/2*( k1 + k2 )
    end do ! ir

  end function g_of_phi


  recursive function sphb(n,x) result(sphb_)
    !
    ! Spherical bessel functions
    !

    implicit none

    real(kind=dp), intent(in) :: x
    integer, intent(in) :: n
    real(kind=dp) :: sphb_

    if (n==0) then
       if (abs(x)<epsilon(x)) then
          sphb_ = 1.0_dp
       else
          sphb_ = sin(x)/x
       end if
    else if (n==1) then
       if (abs(x)<epsilon(x)) then
          sphb_ = 0.0_dp
       else
          sphb_ = sin(x)/x**2 - cos(x)/x
       end if
    else if (n>1) then
       if (abs(x)<epsilon(x)) then
          sphb_ = 0.0_dp
       else
          sphb_ = - sphb(n-2,x) + (2.0_dp* n - 1.0_dp ) * sphb(n-1,x)/x
       end if

    end if

  end function sphb


  subroutine irreduciblek(nk1, nk2, nk3, nirr, kirr)
    !
    ! Find irreducible k-points
    !

    ! variables
    use precision, only: dp
    use siesta2intw_symmetry, only: nsym, s

    implicit none

    integer, intent(in) :: nk1, nk2, nk3 ! k-mesh
    integer, intent(out) :: nirr ! number of irreducible k-points
    real(kind=dp), dimension(3, nk1*nk2*nk3), intent(out) :: kirr ! irreducible k-points incrystal coordinates
    logical, dimension(nk1, nk2, nk3) :: done

    integer, dimension(3) :: kv, skv
    integer :: i, j, k, mi, mj, mk
    integer :: isym

    done = .false.
    nirr = 0
    do i=1,nk1
      do j=1,nk2
        do k=1,nk3
          !
          kv = (/i-1, j-1, k-1/)
          !
          if (.not.done(i,j,k)) then
            !
            nirr = nirr + 1
            kirr(:, nirr) = real(kv, kind=dp)/(/nk1, nk2, nk3/)
            do isym=1,nsym
              skv = matmul(s(:,:,isym), kv)
              mi = modulo(skv(1), nk1) + 1
              mj = modulo(skv(2), nk2) + 1
              mk = modulo(skv(3), nk3) + 1
              done(mi,mj,mk) = .true.
            end do!isym
            !
          end if
          !
        end do!i
      end do!j
    end do!k

  end subroutine irreduciblek

end module siesta2intw_write
