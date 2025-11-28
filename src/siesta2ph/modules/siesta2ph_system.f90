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
module siesta2ph_system

  use precision, only: dp

  implicit none

  public :: slabel, nat, ntyp, alat, at, bg, tau, tau_cryst, ityp, atom_labels, amass, magnetic

  public :: read_unit_cell_data, print_unit_cell_data

  private

  ! system variables
  character(len=256) :: slabel
  integer :: nat, ntyp
  real(kind=dp) :: alat
  real(kind=dp), dimension(3,3) :: at, bg
  real(kind=dp), allocatable, dimension(:,:) :: tau
  real(kind=dp), allocatable, dimension(:,:) :: tau_cryst
  integer, allocatable, dimension(:) :: ityp
  character(len=3), allocatable, dimension(:) :: atom_labels
  real(kind=dp), allocatable, dimension(:) :: amass
  logical :: magnetic

contains

  subroutine read_unit_cell_data(fdffilename)
    !
    ! Read unit cell data from a fdf file.
    !
    use m_struct_init, only: struct_init
    use chemical, only: read_chemical_types, number_of_species, atomic_number, species_label
    use periodic_table, only: atmass
    use siesta_geom, only: na_u, ucell, xa, isa
    use fdf, only: block_fdf, parsed_line, fdf_init, fdf_string, fdf_physical, fdf_block, &
                   fdf_bline, fdf_bmatch, fdf_bintegers, fdf_bvalues, leqi, fdf_get, &
                   fdf_bnnames, fdf_bnintegers, fdf_bnreals, fdf_breals, fdf_bnames, fdf_defined
    !
    use siesta2ph_io, only: stdout, prefix, use_sym
    use siesta2ph_utils, only: pi
    !
    use siesta2ph_linalg, only: ainv
    !
    implicit none
    !
    character(len=*), optional, intent(in) :: fdffilename
    character(len=256) :: fdfin, fdfout
    type(block_fdf) :: bfdf
    type(parsed_line), pointer :: pline
    real(kind=dp), dimension(3) :: origin, originc
    logical :: lOrigin
    character(len=256) :: acf_default, acf
    integer :: ia, id
    character(len=32) :: opt, opt_old
    logical :: spin_none, spin_Col, spin_NCol, spin_SO
    logical :: noncol, peratm, badsyntax
    integer :: ni, nn, nr
    character(len=1) :: updo
    real(kind=dp), allocatable, dimension(:) :: spin, theta, phi
    integer, dimension(:), allocatable :: atom
    logical :: can_be_magnetic


    write(stdout,*) "- Reading fdf data..."
    !
    if (present(fdffilename)) then
      fdfin = trim(fdffilename)
    else
      fdfin = trim(prefix)
    endif
    fdfout = "fdf.log"
    call fdf_init(fdfin, fdfout)
    !
    slabel = fdf_string("SystemLabel","siesta")
    !
    alat = fdf_physical("LatticeConstant",0.0_dp,"Bohr")
    !
    call coor(na_u, ucell)  ! Sets na_u, xa, and isa
    !
    nat = na_u
    at = ucell/alat ! Lattice vectors by column
    bg = transpose(ainv(at)) ! Reciprocal lattice vectors by column
    !
    call read_chemical_types()
    ntyp = number_of_species()
    !
    allocate(ityp(nat))
    ityp(1:nat) = isa(1:nat)
    !
    allocate(atom_labels(ntyp))
    do ia=1,ntyp
      atom_labels(ia) = species_label(ia)
    enddo
    !
    allocate(amass(ntyp))
    do ia=1,ntyp
      amass(ia) = atmass(abs(int(atomic_number(ia))))
    enddo
    if (fdf_block("AtomicMass",bfdf)) then
      do while(fdf_bline(bfdf,pline))
        if (.not. fdf_bmatch(pline,"iv")) call die("Wrong format in AtomicMass")
        ia = fdf_bintegers(pline,1)
        if (ia .gt. ntyp .or. ia .lt. 1) call die("Wrong specnum in AtomicMass")
        amass(ia) = fdf_bvalues(pline,2)
        write(stdout,"(a, i4, a, f12.5)") "remass: Read atomic mass for species ", ia, " as ", amass
      enddo
    endif
    !
    ! In the symmetries calculated by phonopy the origin is not taken into
    ! account, thus it should be readed and removed.
    lOrigin = fdf_block("AtomicCoordinatesOrigin",bfdf)
    ! Read Origin if present
    if (lOrigin) then
      if (.not. fdf_bline(bfdf,pline)) call die("coor: ERROR in AtomicCoordinatesOrigin block")
      origin(1) = fdf_bvalues(pline,1)
      origin(2) = fdf_bvalues(pline,2)
      origin(3) = fdf_bvalues(pline,3)
    else
      origin = 0.0_dp
    endif
    !
    ! Transform the origin according to the AtomicCoordinatesFormat
    acf_default = "NotScaledCartesianBohr"
    acf = fdf_string("AtomicCoordinatesFormat",acf_default)
    if (leqi(acf,"NotScaledCartesianBohr") .or. leqi(acf,"Bohr") ) then
      ! Do nothing
    else if (leqi(acf,"NotScaledCartesianAng") .or. leqi(acf,"Ang") ) then
      origin = origin / 0.529177_dp
    else if (leqi(acf,"ScaledCartesian")) then
      origin = origin * alat
    else if (leqi(acf,"ScaledByLatticeVectors") .or.  leqi(acf,"Fractional") ) then
      do ia=1,nat
        do id=1,3
          originc(id) = origin(id)
        enddo
        do id=1,3
          origin(id) = ucell(id,1) * originc(1) + &
                       ucell(id,2) * originc(2) + &
                       ucell(id,3) * originc(3)
        enddo
      enddo
    endif
    !
    ! Remove the origin
    allocate(tau(3,nat))
    allocate(tau_cryst(3,nat))
    do ia=1,nat
      tau(:,ia) = (xa(:,ia)-origin(:))/alat
      tau_cryst(:,ia) = matmul(ainv(at), tau(:,ia))
    enddo
    !
    ! Read type of spin polarization
    spin_none = .false.
    spin_Col = .false.
    spin_NCol = .false.
    spin_SO = .false.

    spin_Col = fdf_get("SpinPolarized", spin_Col)
    spin_NCol = fdf_get("NonCollinearSpin", spin_NCol)
    spin_SO = fdf_get("SpinOrbit", spin_SO)

    if ( spin_SO ) then
      opt_old = "spin-orbit"
    else if ( spin_NCol ) then
      opt_old = "non-collinear"
    else if ( spin_Col ) then
      opt_old = "polarized"
    else
      opt_old = "none"
    end if
    !
    opt = fdf_get("Spin", opt_old)
    !
    if ( leqi(opt, "none") .or. &
         leqi(opt, "non-polarized") .or. &
         leqi(opt, "non-polarised") .or. &
         leqi(opt, "NP") .or. leqi(opt,"N-P") ) then
      !
      spin_none = .true.
      call warn_and_set_to_false(spin_Col, opt_old, opt)
      call warn_and_set_to_false(spin_NCol, opt_old, opt)
      call warn_and_set_to_false(spin_SO, opt_old, opt)
      !
    else if ( leqi(opt, "polarized") .or. &
              leqi(opt, "collinear") .or. leqi(opt, "colinear") .or. &
              leqi(opt, "polarised") .or. leqi(opt, "P") ) then
      !
      spin_Col = .true.
      call warn_and_set_to_false(spin_none, opt_old, opt)
      call warn_and_set_to_false(spin_NCol, opt_old, opt)
      call warn_and_set_to_false(spin_SO, opt_old, opt)
      !
    else if ( leqi(opt, "non-collinear") .or. leqi(opt, "non-colinear") .or. &
              leqi(opt, "NC") .or. leqi(opt, "N-C") ) then
      !
      spin_NCol = .true.
      call warn_and_set_to_false(spin_none, opt_old, opt)
      call warn_and_set_to_false(spin_Col, opt_old, opt)
      call warn_and_set_to_false(spin_SO, opt_old, opt)
      !
    else if ( leqi(opt, "spin-orbit") .or. leqi(opt, "S-O") .or. &
              leqi(opt, "SOC") .or. leqi(opt, "SO") ) then
      !
      spin_SO = .true.
      call warn_and_set_to_false(spin_none, opt_old, opt)
      call warn_and_set_to_false(spin_Col, opt_old, opt)
      call warn_and_set_to_false(spin_NCol, opt_old, opt)
      !
    else
      stop "ERROR: read_unit_cell_data: unknown 'Spin' flag"
    end if
    !
    if (spin_none) then
      can_be_magnetic = .false.
    else if (spin_Col .or. spin_NCol .or. spin_SO) then
      can_be_magnetic = .true.
    else
      stop "ERROR: read_unit_cell_data: Wrong 'Spin' flag"
    endif
    !
    ! Read initial spin configuration
    allocate(spin(nat), theta(nat), phi(nat))
    spin = 0.0_dp
    theta = 0.0_dp
    phi = 0.0_dp
    allocate(atom(nat))
    atom = 0
    !
    if (can_be_magnetic .and. fdf_defined("DM.InitSpin")) then
      noncol = .false.
      peratm = fdf_block("DM.InitSpin",bfdf)
      ia = 0
      badsyntax = .false.
      do while(fdf_bline(bfdf,pline) .and. (ia .lt. nat) .and. (.not. badsyntax))
        !
        nn = fdf_bnnames(pline)
        ni = fdf_bnintegers(pline)
        nr = fdf_bnreals(pline)
        !
        if (ni .eq. 1) then
          !
          ia = ia + 1
          atom(ia) = fdf_bintegers(pline,1)
          !
          if (nn .eq. 0) then
            !
            ! Read value of spin
            if (nr .eq. 3) then
              ! Read spin value and direction
              spin(ia) = fdf_breals(pline,1)
              theta(ia) = fdf_breals(pline,2) * pi/180.0d0
              phi(ia) = fdf_breals(pline,3) * pi/180.0d0
            elseif (nr .eq. 1) then
              ! Read spin value. Default direction.
              spin(ia) = fdf_breals(pline,1)
              theta(ia) = 0.d0
              phi(ia) = 0.d0
            else
              ! Print bad-syntax error and stop
              badsyntax = .true.
            endif
            !
          elseif (nn .eq. 1) then
            !
            ! Read spin as + or - (maximun value)
            updo = fdf_bnames(pline,1)
            if (updo .eq. "+") then
              spin(ia) = 100.d0
            elseif (updo .eq. "-") then
              spin(ia) = -100.d0
            else
              ! Print bad-syntax error and stop
              badsyntax = .true.
            endif
            !
            if (nr .eq. 2) then
              theta(ia) = fdf_breals(pline,1) * pi/180.0d0
              phi(ia) = fdf_breals(pline,2) * pi/180.0d0
            elseif (nr .eq. 0) then
              theta(ia) = 0.d0
              phi(ia) = 0.d0
            else
              ! Print bad-syntax error and stop
              badsyntax = .true.
            endif
            !
          else
            !
            ! Print bad-syntax error and stop
            badsyntax = .true.
          endif
          !
          if ((atom(ia) < 1) .or. (atom(ia) > nat)) stop "ERROR: read_unit_cell_data: Bad atom index in 'DM.InitSpin'"
          !
          if (abs(theta(ia)) > 1.d-12) noncol = .true.
          !
        else
          ! Print bad-syntax error and stop
          badsyntax = .true.
        endif
        !
      enddo
      !
      if (badsyntax) stop "ERROR: read_unit_cell_data: Bad syntax in 'DM.InitSpin'"
      !
      if (noncol .and. spin_Col) write(stdout,"(a)") "WARNING: read_unit_cell_data: noncollinear spins in 'DM.InitSpin' not used because 'Spin' flag is set to colinear"
      !
      if (any(abs(spin)>0.00000001_dp)) then
        magnetic = .true.
      else
        magnetic = .false.
      endif
      !
    else if (can_be_magnetic .and. (.not.fdf_defined("DM.InitSpin"))) then
      !
      ! From init_DM_atomic in m_new_dm.F90:
      ! All atoms with maximum polarization compatible with atomic configuration. In Ferromagnetic ordering (up).
      magnetic = .true.
      !
    else
      !
      magnetic = .false.
      !
    endif
    !
    ! Check if the spin configuration makes sense or not
    if ( (spin_Col .or. spin_NCol) .and. (.not. magnetic) ) then
      if (spin_Col) write(stdout,"(a)") "WARNING: read_unit_cell_data: A non-magnetic colinear calculation makes no sense!"
      if (spin_NCol) write(stdout,"(a)") "WARNING: read_unit_cell_data: A non-magnetic non-colinear calculation makes no sense!"
      stop "STOPPING..."
    endif
    !
    ! Disable symmetries for spin-orbit calculations for precaution
    ! until rotation of 2x2 induced potential is tested
    if (spin_SO .and. use_sym) then
      write(stdout,"(a)") "WARNING: Spin-orbit calculation with symmetries not implemented yet: Symmetries will not be used!"
      use_sym = .false.
    endif
    !
    ! ! INFO
    ! ! a1, a2 and a3 vectors are stored by column in at
    ! ! The same for b1, b2 and b3 vectors on bg
    ! ! ai.bj = delta_ij, thus, transpose(at).bg = I
    ! ! or transpose(bg).at = I
    ! print*, ""
    ! print"(3f)", transpose(at)
    ! print"(a,3f)", "a1=", at(:,1)
    ! print"(a,3f)", "a2=", at(:,2)
    ! print"(a,3f)", "a3=", at(:,3)
    ! print*, ""
    ! print"(3f)", transpose(bg)
    ! print"(a,3f)", "b1=", bg(:,1)
    ! print"(a,3f)", "b2=", bg(:,2)
    ! print"(a,3f)", "b3=", bg(:,3)
    ! print*, ""
    ! print"(a,3f)", "a1.b1=", dot_product(at(:,1),bg(:,1))
    ! print"(a,3f)", "a1.b2=", dot_product(at(:,1),bg(:,2))
    ! print"(a,3f)", "a1.b3=", dot_product(at(:,1),bg(:,3))
    ! print"(a,3f)", "a2.b1=", dot_product(at(:,2),bg(:,1))
    ! print"(a,3f)", "a2.b2=", dot_product(at(:,2),bg(:,2))
    ! print"(a,3f)", "a2.b3=", dot_product(at(:,2),bg(:,3))
    ! print"(a,3f)", "a3.b1=", dot_product(at(:,3),bg(:,1))
    ! print"(a,3f)", "a3.b2=", dot_product(at(:,3),bg(:,2))
    ! print"(a,3f)", "a3.b3=", dot_product(at(:,3),bg(:,3))
    ! print*, ""
    ! print*, "aT.b="
    ! print"(3f)", transpose(matmul(transpose(at),(bg)))
    ! !
    ! ! Since a1, a2 and a3 are stored by column in at,
    ! ! X_cart = at. X_cryst, and, X_cryst = ainv(at).X_cart
    ! write(stdout,"(3f)") tau(:,1)*alat ! tau is in cartesian coordinates, in units of alat
    ! write(stdout,"(3f)") matmul(ainv(at),tau(:,1)*alat) ! This gives x_cryst
    ! !INFO

    contains

    subroutine warn_and_set_to_false(var, opt_deprecated, opt_spin)

      implicit none

      logical, intent(inout) :: var
      character(len=32), intent(in) :: opt_spin, opt_deprecated


      if (var) then
          write(stdout,*) "WARNING","Deprecated spin keyword overridden by new-style 'Spin' input"
          write(stdout,*) "WARNING","Option from deprecated keyword: "//trim(opt_deprecated)
          write(stdout,*) "WARNING","Option from 'Spin' keyword input: "//trim(opt_spin)
          var = .false.
      endif

    end subroutine warn_and_set_to_false

  end subroutine read_unit_cell_data


  subroutine print_unit_cell_data()
    !
    ! Print system's unit cell data in a nice format.
    !
    use siesta2ph_io, only: stdout
    !
    implicit none
    !
    integer :: ia


    write(stdout,*) "- Unit cell:"
    !
    write(stdout,*) " alat"
    write(stdout,*) alat
    ! !
    write(stdout,"(a)") "         a1        a2        a3"
    write(stdout,"(a1,3f10.5,a1)") "(", at(1,:), ")"
    write(stdout,"(a1,3f10.5,a1)") "(", at(2,:), ")"
    write(stdout,"(a1,3f10.5,a1)") "(", at(3,:), ")"
    !
    write(stdout,"(a)") "         b1        b2        b3"
    write(stdout,"(a1,3f10.5,a1)") "(", bg(1,:), ")"
    write(stdout,"(a1,3f10.5,a1)") "(", bg(2,:), ")"
    write(stdout,"(a1,3f10.5,a1)") "(", bg(3,:), ")"
    !
    write(stdout,*) " ntyp"
    write(stdout,*) ntyp
    !
    write(stdout,*) " amass"
    do ia=1,ntyp
      write(stdout,"(i4,x,a2,3f10.5)") ia, atom_labels(ia), amass(ia)
    enddo
    !
    write(stdout,*) " nat"
    write(stdout,*) nat
    !
    write(stdout,*) " tau"
    do ia=1,nat
      write(stdout,"(i4,i2,x,a2,3f10.5)") ia, ityp(ia), atom_labels(ityp(ia)), tau(:,ia)
    enddo

  end subroutine print_unit_cell_data

end module siesta2ph_system
