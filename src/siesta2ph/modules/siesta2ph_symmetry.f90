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
module siesta2ph_symmetry

  use precision, only: dp
  use spglib_f08, only: spg_get_dataset, SpglibDataset

  implicit none

  public :: nsym, s, sinv, s_cart, ftau, nsite_sym, site_atm, nsite_rot, site_s, &
            site_sinv, site_s_cart, nirred, irred_disp, irred_atm, atom_mapping, &
            rtau_index, rtau_cryst

  public :: get_q_mesh, compute_symmetry, print_symmetry_data, &
            find_irreducible_atoms, find_site_symmetry, &
            set_displacements_direction, find_irreducible_displacements

  private

  ! symmetry variables
  integer :: nsym
  integer, allocatable, dimension(:,:,:) :: s, sinv
  real(kind=dp), allocatable, dimension(:,:,:) :: s_cart
  real(kind=dp), allocatable, dimension(:,:) :: ftau
  integer :: nsite_sym
  integer, allocatable, dimension(:) :: site_atm
  integer, allocatable, dimension(:) :: nsite_rot
  integer, allocatable, dimension(:,:,:,:) :: site_s, site_sinv
  real(kind=dp), allocatable, dimension(:,:,:,:) :: site_s_cart
  integer :: nirred
  logical, allocatable, dimension(:,:) :: irred_disp
  logical, allocatable, dimension(:) :: irred_atm
  integer, allocatable, dimension(:) :: atom_mapping
  integer, allocatable :: rtau_index(:,:)
  real(kind=dp), allocatable :: rtau_cryst(:,:,:)

contains

  subroutine get_q_mesh(nq1, nq2, nq3, qmesh_cryst, qirr_map, irreducible, nqirr, qirr_cryst)
    !
    ! Find q-mesh, qirr_map and optionally irreducible q-points
    !
    use siesta2ph_io, only: stdout
    use siesta2ph_system, only: nat, tau_cryst, at, ityp
    !
    use spglib_f08, only: spg_get_ir_reciprocal_mesh

    implicit none

    integer, intent(in) :: nq1, nq2, nq3 ! k-mesh
    real(kind=dp), dimension(3, nq1*nq2*nq3), intent(out) :: qmesh_cryst ! k-points mesh in crystal coordinates
    integer, intent(out) :: qirr_map(nq1*nq2*nq3)
    logical, intent(in), optional :: irreducible
    integer, intent(out), optional :: nqirr
    real(kind=dp), dimension(3, nq1*nq2*nq3), intent(out), optional :: qirr_cryst
    !
    integer :: mesh(3), is_shift(3), is_time_reversal
    integer :: spglib_qmesh(3, nq1*nq2*nq3), spglib_qirr_map(nq1*nq2*nq3)
    integer :: spglib_nqirr, iq, iqirr
    logical :: qirr_done(nq1*nq2*nq3)


    write(stdout,*) "- Computing q-mesh..."
    write(stdout,"(a16,3i6)") "nq1, nq2, nq3 = ", nq1, nq2, nq3

    ! Compute irreducible q-mesh
    mesh = (/nq1, nq2, nq3/)
    is_shift = (/0, 0, 0/)
    is_time_reversal = 0
    spglib_nqirr = spg_get_ir_reciprocal_mesh(spglib_qmesh, spglib_qirr_map, mesh, is_shift, is_time_reversal, transpose(at), tau_cryst, ityp, nat, 1.d-5)
    !
    qirr_map = spglib_qirr_map + 1
    !
    do iq = 1, nq1*nq2*nq3
      qmesh_cryst(:, iq) = real(spglib_qmesh(:, iq), kind=dp)/(/nq1, nq2, nq3/)
    enddo

    if (present(irreducible) .or. present(nqirr) .or. present(qirr_cryst)) then
      !
      if (present(irreducible) .and. present(nqirr) .and. present(qirr_cryst)) then
        !
        if (.not. irreducible) then
          nqirr = nq1*nq2*nq3
          qirr_cryst = qmesh_cryst
          return
        endif
        !
        write(stdout,*) "- Computing irreducible q-points..."
        !
        qirr_done = .false.
        iqirr = 0
        do iq = 1, nq1*nq2*nq3
          !
          if (qirr_done(qirr_map(iq))) cycle
          !
          iqirr = iqirr + 1
          qirr_cryst(:, iqirr) = qmesh_cryst(:, iq)
          qirr_done(qirr_map(iq)) = .true.
          !
        enddo
        !
        nqirr = iqirr
        !
      else
        stop "ERROR: get_q_mesh: irreducible, nqirr and qirr_cryst must be present"
      endif
      !
    endif

  end subroutine get_q_mesh


  subroutine compute_symmetry()
    !
    ! Compute the symmetry of the system using spglib
    !
    use siesta2ph_io, only: use_sym, stdout
    use siesta2ph_system, only: nat, tau_cryst, at, ityp, magnetic
    !
    use siesta2ph_linalg, only: ainv
    !
    implicit none
    !
    type(SpglibDataset) :: spg_sym
    !
    integer :: isym


    if (magnetic .and. use_sym) then
      use_sym = .false.
      write(stdout,"(a)") "WARNING: Initial magnetic moment detected: Symmetries will not be used!"
    endif
    !
    spg_sym = spg_get_dataset(transpose(at), tau_cryst, ityp, nat, 1.d-5)
    !
    if (use_sym) then
      nsym = spg_sym%n_operations
    else
      ! Keep only pure translations
      nsym = 0
      do isym=1,spg_sym%n_operations
        if (check_identity(spg_sym%rotations(:,:,isym))) nsym = nsym + 1
      enddo
    endif
    !
    allocate(s(3,3,nsym))
    allocate(sinv(3,3,nsym))
    allocate(ftau(3,nsym))
    allocate(s_cart(3,3,nsym))
    !
    nsym = 0
    do isym=1,spg_sym%n_operations
      if ((.not.use_sym) .and. (.not.check_identity(spg_sym%rotations(:,:,isym)))) cycle
      nsym = nsym + 1
      s(:,:,nsym) = transpose(spg_sym%rotations(:,:,isym))
      ftau(:,nsym) = -spg_sym%translations(:,isym) ! The minus sign is to agree with the QE convention
      sinv(:,:,nsym) = ainv(real(s(:,:,nsym),dp))
      s_cart(:,:,nsym) = matmul(at, matmul(s(:,:,nsym), ainv(at)))
    enddo
    !
    call rot_atoms()

    contains

    function check_identity(matrix)

      implicit none

      integer, dimension(3,3), intent(in) :: matrix
      logical :: check_identity


      check_identity = all(abs(matrix - reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/)))==0)

    end function check_identity

  end subroutine compute_symmetry


  subroutine print_symmetry_data()
    !
    ! Print the symmetry data of the system in a nice format.
    !
    use siesta2ph_io, only: stdout
    use siesta2ph_system, only: nat
    !
    implicit none
    !
    integer :: isym, ia, isite_sym


    write(stdout,*) "- Symmetry operations:"

    write(stdout,*) " nsym"
    write(stdout,*) nsym

    23 format("  - [",i2,", ",i2,",",i2,"]")
    24 format("  translation: [",f8.5,", ",f8.5,", ",f8.5," ]")
    do isym=1,nsym
      write(stdout,*) " s", isym
      write(stdout,23) s(1,:,isym)
      write(stdout,23) s(2,:,isym)
      write(stdout,23) s(3,:,isym)
      write(stdout,24) ftau(:,isym)
    enddo

    write(stdout,*) " atom_mapping"
    do ia=1,nat
      write(stdout,"(i4,a1,i4)") ia, ":", atom_mapping(ia)
    enddo

    write(stdout,*) " nsite_sym"
    write(stdout,*) nsite_sym
    do isite_sym=1,nsite_sym
      write(stdout,*) " site symmetry", isite_sym
      write(stdout,*) " atom", site_atm(isite_sym)
      write(stdout,*) " rotations", nsite_rot(isite_sym)
      do isym=1,nsite_rot(isite_sym)
        write(stdout,*) " s", isym
        write(stdout,23) site_s(1,:,isym,isite_sym)
        write(stdout,23) site_s(2,:,isym,isite_sym)
        write(stdout,23) site_s(3,:,isym,isite_sym)
      enddo
    enddo

  end subroutine print_symmetry_data


  subroutine rot_atoms()
    !
    ! Compute the rotated atomic positions (rtau_cryst) and equivalences (rtau_index).
    !
    use siesta2ph_system, only: nat, tau_cryst
    !
    use siesta2ph_utils, only: eqvect
    !
    implicit none
    !
    integer :: ia, ja
    integer :: isym
    !
    real(kind=dp), parameter :: epsat = 1E-3


    allocate(rtau_index(nat,nsym))
    allocate(rtau_cryst(3,nsym,nat))
    !
    rtau_index = 0
    !
    do ia=1,nat
      do isym=1,nsym
        rtau_cryst(:,isym,ia) = matmul(s(:,:,isym), tau_cryst(:,ia)) - dble(ftau(:,isym))
      enddo
    enddo
    !
    do ia=1,nat
      do isym=1,nsym
        do ja=1,nat
          if ( eqvect(rtau_cryst(:,isym,ia), tau_cryst(:,ja)) ) rtau_index(ia,isym) = ja
        enddo
      enddo
    enddo
    !
    do ia=1,nat
      do isym=1,nsym
        if (rtau_index(ia,isym).eq.0) then
          write(*,*) "ERROR in rot_at: At least one atom does not map properly under sym. op.", isym, "atom:", ia
        endif
      enddo
    enddo

  end subroutine rot_atoms


  subroutine find_irreducible_atoms()
    !
    ! TODO: Add description
    !
    use siesta2ph_system, only: nat
    !
    implicit none
    !
    integer, allocatable, dimension(:) :: tmp_irred_atm_list
    integer :: ia, isym


    allocate(irred_atm(nat))
    allocate(tmp_irred_atm_list(nat))
    irred_atm = .false.
    tmp_irred_atm_list = 0
    !
    do ia=1,nat
      !
      if ( any(ia==tmp_irred_atm_list) ) then
        cycle
      else
        irred_atm(ia) = .true.
      endif
      !
      do isym=1,nsym
        tmp_irred_atm_list(rtau_index(ia, isym)) = rtau_index(ia, isym)
      enddo
      !
    enddo ! ia
    !
    deallocate(tmp_irred_atm_list)
    !
    ! Calculate atom mapping
    allocate(atom_mapping(nat))
    atom_mapping = 0
    !
    do ia=1,nat
      !
      if ( .not. irred_atm(ia) ) cycle
      !
      atom_mapping(ia) = ia
      !
      do isym=1,nsym
        atom_mapping(rtau_index(ia, isym)) = ia
      enddo
      !
    enddo

  end subroutine find_irreducible_atoms


  subroutine find_site_symmetry()
    !
    ! TODO: Add description
    !
    use siesta2ph_system, only: nat, tau_cryst, at
    !
    implicit none
    !
    real(kind=dp), dimension(3) :: R, SR, dR
    integer :: ia, isym, isite_sym


    nsite_sym = count(irred_atm)
    allocate(nsite_rot(nsite_sym))
    allocate(site_atm(nsite_sym))
    allocate(site_s(3,3,nsym,nsite_sym))
    allocate(site_sinv(3,3,nsym,nsite_sym))
    allocate(site_s_cart(3,3,nsym,nsite_sym))
    nsite_rot = 0
    site_atm = 0
    site_s = 0
    site_sinv = 0
    site_s_cart = 0
    !
    isite_sym = 0
    do ia=1,nat
      !
      if ( .not. irred_atm(ia) ) cycle
      !
      isite_sym = isite_sym + 1
      site_atm(isite_sym) = ia
      R = tau_cryst(:,ia)
      !
      do isym=1,nsym
        !
        SR = matmul(s(:,:,isym),R) + ftau(:,isym)
        dR = R - SR
        dR = dR - nint(dR)
        dR = matmul(dR,at)
        !
        if ( sqrt(sum(dR**2)) < 1.d-5 ) then
          nsite_rot(isite_sym) = nsite_rot(isite_sym) + 1
          site_s(:,:,nsite_rot(isite_sym),isite_sym) = s(:,:,isym)
          site_sinv(:,:,nsite_rot(isite_sym),isite_sym) = sinv(:,:,isym)
          site_s_cart(:,:,nsite_rot(isite_sym),isite_sym) = s_cart(:,:,isym)
        endif
        !
      enddo ! isym
      !
    enddo ! ia

  end subroutine find_site_symmetry


  subroutine set_displacements_direction(disp_cart)
    !
    ! disp_cart returns in the output the displacement vectors in Cartesian coordinates.
    ! If disp_along_cart == .true., the displacement vectors will be along the Cartesian directions,
    ! if disp_along_cart == .false., the displacement vectors will be along the Crystal directions.
    !
    use siesta2ph_io, only: stdout, verbose, disp_along_cart, dx
    use siesta2ph_system, only: at
  !
    implicit none
    !
    real(kind=dp), dimension(3,3), intent(out) :: disp_cart


    if (disp_along_cart) then
      !
      disp_cart(:,1) = (/ 1.0_dp, 0.0_dp, 0.0_dp /) ! x direction
      disp_cart(:,2) = (/ 0.0_dp, 1.0_dp, 0.0_dp /) ! y direction
      disp_cart(:,3) = (/ 0.0_dp, 0.0_dp, 1.0_dp /) ! z direction
      !
    else
      !
      disp_cart(:,1) = at(:,1)/norm2(at(:,1)) ! a1 direction
      disp_cart(:,2) = at(:,2)/norm2(at(:,2)) ! a2 direction
      disp_cart(:,3) = at(:,3)/norm2(at(:,3)) ! a3 direction
      !
    endif
    !
    disp_cart = disp_cart*dx
    !
    if (verbose) then
      if (disp_along_cart) then
        write(stdout,"(a)") "Displacements along Cartesian directions will be used:"
      else
        write(stdout,"(a)") "Displacements along Crystal directions will be used:"
      endif
      write(stdout,"(3f)") disp_cart(:,1)
      write(stdout,"(3f)") disp_cart(:,2)
      write(stdout,"(3f)") disp_cart(:,3)
    endif

  end subroutine set_displacements_direction


  subroutine find_irreducible_displacements()
    !
    ! TODO: Add description
    !
    use siesta2ph_io, only: stdout
    use siesta2ph_system, only: nat, at
    !
    use siesta2ph_linalg, only: rank, ainv
    !
    implicit none
    !
    integer :: ia, id, isym, isite_sym
    real(kind=dp), dimension(3,3) :: disp
    real(kind=dp), dimension(3,3) :: basis
    integer :: nbasis, basis_rank


    write(stdout,*) "- Finding irreducible displacements..."
    !
    call set_displacements_direction(disp)
    !
    allocate(irred_disp(3,nat))
    irred_disp = .false.
    !
    do ia=1,nat
      !
      if ( .not. irred_atm(ia) ) cycle
      !
      ! Find site symmetry of the atom
      do isite_sym=1,nsite_sym
        if ( atom_mapping(ia) == site_atm(isite_sym)) exit
      enddo
      ! Check if site symmetry has been found correctly
      if ((isite_sym == nsite_sym) .and. ( atom_mapping(ia) /= site_atm(isite_sym))) stop "find_irreducible_displacements: ERROR: site symmetry not found"
      !
      !
      basis = 0.0_dp
      nbasis = 0
      basis_rank = 0
      do id=1,3
        !
        basis(:,nbasis+1) = disp(:,id)
        call rank(basis, basis_rank)
        !
        ! If the rank has increased, id is an irreducible direction for atom ia
        if ( basis_rank > nbasis ) then
          nirred = nirred + 1
          irred_disp(id,ia) = .true.
          nbasis = basis_rank
          write(stdout,"(a3,i3,x,a3,i3,x,a)") "ia=", ia, "id=", id, "is an irreducible displacement."
        elseif ( basis_rank == nbasis ) then
          irred_disp(id,ia) = .false.
          write(stdout,"(a3,i3,x,a3,i3,x,a)") "ia=", ia, "id=", id, "is not an irreducible displacement."
          cycle
        else
          stop "ERROR: find_irreducible_displacements: the rank has decreased"
        endif
        !
        do isym=1,nsite_rot(isite_sym)
          !
          if (nbasis == 3) exit
          !
          basis(:,nbasis+1) = matmul(transpose(site_s_cart(:,:,isym,isite_sym)), disp(:,id))
          call rank(basis,nbasis)
          !
        enddo ! isym
        !
      enddo ! id
      !
    enddo ! ia

  end subroutine find_irreducible_displacements

end module siesta2ph_symmetry
