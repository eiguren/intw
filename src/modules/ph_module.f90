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
module intw_ph

  !! display: none
  !!
  !! This module contains the main variables related to phonon modes.
  !!

  use kinds, only: dp

  implicit none
  !
  save
  ! variables
  public :: q_irr_cryst, u_irr, dvscf_cart, dvscf_irr, &
            nqmesh, qmesh, &
            QE_folder_nosym_q, QE_folder_sym_q, symlink_q
  !
  ! subroutines
  public :: rot_gep, read_ph_information, readfc, mat_inv_four_t, read_allq_dvr, &
            get_dv, rot_dvq, func_by_gr, wsinit, &
            read_dynq, set_asr_frc, set_asr_dynq, rotate_dyn
  !
  ! functions
  public :: wsweight, rot_k_index
  !
  private
  !
  !
  real(dp), allocatable :: q_irr_cryst(:,:) ! coo. of irr. q points
  complex(dp), allocatable :: u_irr(:,:,:) ! displacement patterns for the irr.  q.
  complex(dp), allocatable :: dvscf_cart(:,:,:,:,:), dvscf_irr(:,:,:,:)

  integer :: nqmesh
  real(dp), allocatable :: qmesh(:,:)

  integer, allocatable :: QE_folder_nosym_q(:)
  integer, allocatable :: QE_folder_sym_q(:)
  integer, allocatable :: symlink_q(:,:)


contains

  subroutine rot_gep(s_index, imq, qpoint_irr, nbnd, nspin, nat, g_matin, g_matout)

    use intw_symmetries, only: rtau_index, inverse_indices, rtau
    use intw_reading, only: bg, s, at
    use intw_useful_constants, only: tpi, cmplx_0, cmplx_i

    implicit none

    ! input
    real(dp), intent(in) :: qpoint_irr(3)
    integer, intent(in) :: imq, nbnd, nspin, nat

    complex(dp), intent(in) :: g_matin(nbnd,nbnd,nspin,nspin,3*nat)
    complex(dp), intent(out) :: g_matout(nbnd,nbnd,nspin,nspin,3*nat)

    ! local
    integer :: s_index, s_inv_index, na
    integer :: rna
    complex(dp) :: phase(nat)
    integer :: kpol, jpol, lpol, ipol, ibnd, jbnd
    real(dp) :: s_cart(3,3), qpoint_irr_cart(3)


    qpoint_irr_cart = matmul(bg, qpoint_irr)

    g_matout = cmplx_0

    s_inv_index = inverse_indices(s_index)

    do na = 1, nat
      phase(na) = &
          exp( - cmplx_i * tpi * dot_product(qpoint_irr_cart(:), rtau(:,s_index,na)) )
    enddo


    do ipol = 1, 3
      do jpol = 1, 3
        !
        s_cart(ipol,jpol) = 0.d0
        do kpol = 1, 3
          do lpol = 1, 3
            s_cart(ipol,jpol) = s_cart(ipol,jpol) + at(ipol,kpol) * &
                  real(s(lpol,kpol,s_index), dp) * bg(jpol,lpol)
          enddo
        enddo
        !
      enddo
    enddo


    do na = 1, nat

      rna = rtau_index(na,s_index)

      do ibnd = 1, nbnd
        do jbnd = 1, nbnd
          !
          do ipol = 1, 3
            do jpol = 1, 3
              g_matout(ibnd, jbnd, 1:nspin, 1:nspin, (rna-1)*3+ipol) = &
                  g_matout(ibnd, jbnd, 1:nspin, 1:nspin, (rna-1)*3+ipol) &
                  + phase(rna) * s_cart(ipol,jpol) * g_matin(ibnd, jbnd, 1:nspin, 1:nspin, (na-1)*3+jpol)
            enddo ! jpol
          enddo ! ipol
          !
        enddo ! jbnd
      enddo ! ibnd

    enddo ! na

    if (imq == 1) g_matout = conjg(g_matout)

  end subroutine rot_gep


  function rot_k_index(s_index, k_index, nk1, nk2, nk3, kmesh)

    use intw_reading, only: s
    use intw_utility, only: find_k_1BZ_and_G, triple_to_joint_index_g

    implicit none

    integer, intent(in) :: s_index, k_index, nk1, nk2, nk3
    real(dp), intent(in) :: kmesh(3,nk1*nk2*nk3)

    integer :: rot_k_index

    real(dp) :: kpoint_1bz(3), kpoint_rot(3)
    integer :: GKQ_bz(3), i, j, k


    kpoint_rot = matmul(s(:,:,s_index), kmesh(:,k_index))

    call find_k_1BZ_and_G(kpoint_rot, nk1, nk2, nk3, i, j, k, kpoint_1bz, GKQ_bz)

    call triple_to_joint_index_g(nk1, nk2, nk3, rot_k_index, i, j, k)

  end function rot_k_index


  subroutine read_ph_information()
    !
    !This subroutine reads the information related to phonons from files.
    !
    use intw_input_parameters, only: outdir, prefix, qlist, nqirr, nq1, nq2, nq3
    use intw_reading, only: nat, bg
    use intw_matrix_vector, only: ainv
    use intw_utility, only: find_free_unit

    implicit none

    integer :: io_unit, ios
    character(256) :: datafile
    real(dp) :: qpoint(3)
    integer :: imode, jmode, iq
    character(8) :: dummy


    write(*,"(A)") "|   PH BZ division is:                              |"
    write(*,"(3(A,I2),38X,A1)") "|   ", nq1, " x", nq2, " x", nq3, "|"
    write(*,"(A)") "|   Input n. of irreducible points:                 |"
    write(*,"(A1,3X,I3,45X,A1)") "|", nqirr, "|"


    write(*,"(A)") "| - Reading irreducible q-point list...             |"

    ! Read irreducible q list
    allocate(q_irr_cryst(3,nqirr))

    io_unit = find_free_unit()
    open(unit=io_unit, file=trim(outdir)//trim(qlist), status="old", iostat=ios)
    if (ios /= 0) stop "ERROR: read_ph_information: error opening qlist."

    do iq = 1, nqirr
      read(io_unit,*) dummy, qpoint(1:3)
      q_irr_cryst(:,iq) = matmul(ainv(bg), qpoint)
    enddo

    close(unit=io_unit)


    write(*,"(A)") "| - Reading displacement patterns...                |"

    io_unit = find_free_unit()
    datafile = trim(outdir)//trim(prefix)//".save.intw/"//"irrq_patterns.dat"
    open(unit=io_unit, file=datafile, status="old", iostat=ios)
    if (ios /= 0) stop "ERROR: read_ph_information: error opening irrq_patterns."

    ! Read displacement patters for each q.
    allocate(u_irr(3*nat,3*nat,nqirr))

    do iq = 1, nqirr
      read(unit=io_unit,fmt=*) dummy
      do imode = 1, 3*nat
        read(unit=io_unit,fmt=*) ( u_irr(jmode,imode,iq), jmode = 1, 3*nat )
      end do
    end do ! iq

    close(unit=io_unit)

  end subroutine read_ph_information


  subroutine readfc(flfrc, frc)
    !
    ! Read QE force constants file
    !
    ! This subroutine is originally distributed as part of the Quantum Espresso code and has
    ! been adapted to INTW:
    !   Copyright (C) 2001-2009 Quantum ESPRESSO group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    !  Modifications by INTW group, 2024:
    !   - Intent(in) variables specified explicitly.
    !   - Remove MPI parts.
    !   - Add some checks.
    !
    use intw_reading, only: nat, ntyp, at, alat, amass, ityp, tau
    use intw_utility, only: find_free_unit
    use intw_useful_constants, only: cmplx_1, eps_6, Ry_to_Ha, pmass
    use intw_input_parameters, only: nq1, nq2, nq3, apply_asr

    implicit none

    character(256), intent(in) :: flfrc
    complex(dp), intent(out) :: frc(nq1,nq2,nq3,3,3,nat,nat) ! force constants

    ! local variables
    integer :: ntyp_fc, nat_fc
    integer :: nr1_fc, nr2_fc, nr3_fc, ibrav_fc, ityp_fc(nat)
    real(dp) :: alat_fc, celldm_fc(6), at_fc(3,3), amass_fc(ntyp), tau_fc(3,nat), zeu_fc(3,3,nat)
    logical :: has_zstar
    real(dp) :: epsil_fc(3,3)
    character(3) :: atm
    character(9) :: symm_type
    integer :: nt
    integer :: i, j, na, nb, m1, m2, m3
    integer :: ibid, jbid, nabid, nbbid, m1bid,m2bid,m3bid
    integer :: io_unit, ios
    real(dp) :: frc_R(nq1,nq2,nq3,3,3,nat,nat)


    io_unit = find_free_unit()
    open(unit=io_unit, file=flfrc, status='old', form='formatted', iostat=ios)
    if ( ios /= 0 ) stop 'ERROR: readfc: error opening flfrc'
    !
    !  read cell data
    !
    read(io_unit,*) ntyp_fc, nat_fc, ibrav_fc, (celldm_fc(i), i = 1, 6)
    if (ibrav_fc == 0) then
      read(io_unit,'(a)') symm_type
      read(io_unit,*) ((at_fc(i,j), i = 1, 3), j = 1, 3)
    end if
    !
    alat_fc = celldm_fc(1)
    !
    ! Some checks
    if (nat_fc /= nat) stop "ERROR: readfc: Wrong number of atoms"
    if (ntyp_fc /= ntyp) stop "ERROR: readfc: Wrong number of species"
    if (ibrav_fc == 0 .and. any(abs(at_fc-at) > eps_6)) stop "ERROR: readfc: Wrong lattice vectors"
    if (abs(alat_fc-alat) > eps_6) stop "ERROR: readfc: Wrong lattice parameter"
    !
    !
    !  read atomic types, positions and masses
    !
    do nt = 1, ntyp
      read(io_unit,*) i, atm, amass_fc(nt)
    end do
    !
    do na = 1, nat
      read(io_unit,*) i, ityp_fc(na), (tau_fc(j,na), j = 1, 3)
    end do
    !
    ! More checks
    if (any(abs(amass_fc/pmass*2-amass) > eps_6)) stop "ERROR: readfc: Wrong atomic weights"
    if (any(abs(ityp_fc-ityp) /= 0)) stop "ERROR: readfc: Wrong atomic types"
    if (any(abs(tau_fc-tau) > eps_6)) stop "ERROR: readfc: Wrong atomic positions"
    !
    !  read macroscopic variable
    !
    read(io_unit,*) has_zstar
    if (has_zstar) then
      read(io_unit,*) ((epsil_fc(i,j), j = 1, 3), i = 1, 3)
      do na = 1, nat
        read(io_unit,*)
        read(io_unit,*) ((zeu_fc(i,j,na), j = 1, 3), i = 1, 3)
      end do
    endif
    !
    read(io_unit,*) nr1_fc, nr2_fc, nr3_fc
    !
    ! Even more checks
    if (nr1_fc /= nq1 .or. nr2_fc /= nq2 .or. nr3_fc /= nq3) stop "ERROR: readfc: Wrong supercell"
    !
    !  read real-space interatomic force constants
    !
    frc(:,:,:,:,:,:,:) = 0.d0
    do i = 1, 3
      do j = 1, 3
        do na = 1, nat
          do nb = 1, nat
            read(io_unit,*) ibid, jbid, nabid, nbbid
            read(io_unit,*) (((m1bid, m2bid, m3bid, frc_R(m1,m2,m3,i,j,na,nb), &
                                m1 = 1, nq1), m2 = 1, nq2), m3 = 1, nq3) ! in Ry/Bohr^2
          end do
        end do
      end do
    end do
    !
    close(unit=io_unit)
    !
    frc_R = frc_R * Ry_to_Ha ! Transform to a.u.
    !
    ! MBR 21/05/2024
    ! Instead of that set_asr module, use only the "simple" method
    ! (worked in the routine as real)
    if (apply_asr) then
      write(*,*) ' Applying ASR (simple) to force constant matrix'
      call set_asr_frc(nq1, nq2, nq3, nat, frc_R)
    end if
    !
    frc = frc_R*cmplx_1 ! make it complex

  end subroutine readfc

  !====================================================================
  ! MBR 21/05/24
  ! Apply ASR "simple" method:
  !   set_asr_frc --> to the real-space force constant matrix
  !                   (real elements)
  !   set_asr_dynr --> to the dynamical matrices after doing dynq_to_dynr
  !====================================================================

  subroutine set_asr_frc(n1, n2, n3, nat, frc_R)

    implicit none

    integer, intent(in) :: n1, n2, n3, nat
    real(dp), intent(inout) :: frc_R(n1,n2,n3,3,3,nat,nat)

    ! Internal
    integer :: i, j, iat, jat, i1, i2, i3
    real(dp) :: suma


    do i = 1, 3 ! coords.
      do j = 1, 3
        !
        do iat = 1, nat ! atom
          !
          suma = 0.0_dp
          do jat = 1, nat ! atom
            do i1 = 1, n1 ! lattice vectors
              do i2 = 1, n2
                do i3 = 1, n3
                  suma = suma + frc_R(i1,i2,i3,i,j,iat,jat)
                end do
              end do
            end do ! lattice vectors
            !
          end do ! jat atom
          !
          frc_R(1,1,1,i,j,iat,iat) = frc_R(1,1,1,i,j,iat,iat) - suma
          !
        end do ! iat atom
        !
      end do
    end do ! coords.

  end subroutine set_asr_frc


  subroutine set_asr_dynq(nq, iq, nat, dynq)
    !
    ! iq is the q=0 index in the mesh
    !
    use intw_useful_constants, only: cmplx_1

    implicit none

    integer, intent(in) :: nq, iq, nat
    complex(dp), intent(inout) :: dynq(3*nat,3*nat,nq)

    ! Internal
    integer :: i, j, iat, jat
    complex(dp) :: suma
    complex(dp) :: dynq_aux(3*nat,3*nat,nq)


    dynq_aux = dynq
    do i = 1, 3 ! coords.
      do j = 1, 3
        !
        do iat = 1, nat ! atom
          !
          suma = 0.0_dp
          do jat = 1, nat ! atom ! sum over cells
            suma = suma + dynq((iat-1)*3+i, (jat-1)*3+j, iq)
          end do ! jat atom
          !
          ! option 1 : apply only to q=0:
          ! dynq_aux((iat-1)*3+i, (iat-1)*3+j, iq) = dynq((iat-1)*3+i, (iat-1)*3+j, iq) - suma
          ! option 2 : generalize to all q:
          dynq_aux((iat-1)*3+i, (iat-1)*3+j, :) = dynq((iat-1)*3+i, (iat-1)*3+j, :) - suma
          !
        end do ! iat atom
        !
      end do
    end do ! coords.

    dynq = dynq_aux

  end subroutine set_asr_dynq


  subroutine read_dynq(dynq)
    !
    ! Read dynamical matrices of irreducible q-points from
    ! prefix.dyn_q(iq)_sym files in prefix.save.intw and store
    ! into dynq for all q-points by using symmetries to get
    ! dynamical matrices of the non-irreducible q-points
    !
    use intw_reading, only: nat, lmag
    use intw_useful_constants, only: cmplx_i, tpi, Ry_to_Ha
    use intw_utility, only: find_free_unit, &
                            triple_to_joint_index_g, find_k_1BZ_and_G
    use intw_input_parameters, only: outdir, prefix, nqirr, nq1, nq2, nq3, &
                                     apply_asr

    implicit none

    ! I/O

    complex(dp), intent(out) :: dynq(3*nat,3*nat,nqmesh)

    ! Local

    ! I/O variables
    character(4) :: iq_str
    character(256) :: dynq_file, intwdir
    integer :: dynq_unit, ierr
    ! q-point variables
    real(dp) :: qirr_cryst(3), qirr_1BZ(3)
    integer :: Gq(3), i, j, k, iqirr2iq
    logical :: iq_done(nqmesh)
    ! Symmetry variables
    integer :: isym, tr
    ! Dynamical matrix variables
    complex(dp) :: dynq_irr(3,3,nat,nat), dynRq_irr(3,3,nat,nat)
    ! Loop variables
    integer :: iqirr, iq, ia, ja


    iq_done = .false.

    ! INTW directory
    intwdir = trim(outdir)//trim(prefix)//".save.intw/"

    do iqirr = 1, nqirr
      !
      ! open file
      if (                   iqirr <   10) write(iq_str,"(i1)") iqirr
      if ( 10 <= iqirr .and. iqirr <  100) write(iq_str,"(i2)") iqirr
      if (100 <= iqirr .and. iqirr < 1000) write(iq_str,"(i3)") iqirr
      dynq_unit = find_free_unit()
      dynq_file = trim(intwdir)//trim(prefix)//".dyn_q"//trim(adjustl(iq_str))
      open(unit=dynq_unit, iostat=ierr, file=dynq_file, form='formatted', status='old')
      if (ierr /= 0 ) then
        write(*,*) 'Error opening .dyn_q file in ', trim(intwdir), ' . Stopping.'
        stop
      end if
      !
      !
      ! Read irreducible q-point
      read(dynq_unit,'(11x,3(f14.9))') qirr_cryst
      !
      ! Find index of irreducible q point in q-mesh
      call find_k_1BZ_and_G(qirr_cryst, nq1, nq2, nq3, i, j, k, qirr_1BZ, Gq)
      call triple_to_joint_index_g(nq1, nq2, nq3, iqirr2iq, i, j, k)
      !
      !
      ! And read dynamical matrix for the irreducible q-point
      do ia = 1, nat
        do ja = 1, nat
          !
          ! Cartesian 3x3 block of this atom pair in dynq matrix
          do i = 1, 3
            read(dynq_unit,*) (dynq_irr(i,j,ia,ja), j = 1, 3)
          enddo
          !
        end do ! ja
      end do ! ia
      !
      !
      ! Loop over symmetry equivalent q-points
      do iq = 1, nqmesh
        !
        if (QE_folder_sym_q(iq) /= QE_folder_sym_q(iqirr2iq)) cycle
        !
        ! Check for consistency
        if (QE_folder_sym_q(iq) /= iqirr) stop "ERROR: read_dynq"
        !
        ! Find symmetry operation
        isym = symlink_q(iq,1)
        tr = symlink_q(iq,2)
        !
        ! Rotate dynamical matrix
        call rotate_dyn(isym, qirr_cryst, dynq_irr(:,:,:,:), dynRq_irr)
        !
        ! Apply TR
        if ( (.not. lmag) .and. tr==1) dynRq_irr = conjg(dynRq_irr) ! For magnetic symmetries TR is applied by rotate_dyn
        !
        do ia = 1, nat
          do ja = 1, nat
            !
            ! Cartesian 3x3 block of this atom pair
            dynq((ia-1)*3+1:ia*3, (ja-1)*3+1:ja*3, iq) = dynRq_irr(:,:,ia,ja)
            !
          enddo ! ja
        enddo ! ia
        !
        iq_done(iq) = .true.
        !
      end do ! iq
      !
      close(dynq_unit)
      !
    end do ! iqirr
    !
    !
    ! Check that all the qpoints in the full mesh have been read
    do iq = 1, nqmesh
      if ( .not. iq_done(iq) ) then
        write(*,*) "ERROR: read_dynq: Failed to read dynq for iq=", iq
        stop
      end if
    end do
    !
    !
    ! Apply ASR to q=0 matrix:
    ! \sum_{ja} dynq( ia, i, ja, j, q=0) = 0
    if (apply_asr) then
      write(*,'(A)') '|   Applying ASR to all q vector indices            |'
      iq = 1 ! q=0 index in mesh
      call set_asr_dynq(nqmesh, iq, nat, dynq)
    end if

  end subroutine read_dynq


  subroutine rotate_dyn(isym, q_cryst, dynq, dynRq)
    !
    ! Rotate dynamical matrix to the symmetry equivalent q point
    !
    use intw_reading, only: nat, at, bg, lmag, s, TR
    use intw_useful_constants, only: eps_6, cmplx_0, cmplx_i, cmplx_1, tpi
    use intw_utility, only: triple_to_joint_index_g, find_k_1BZ_and_G
    use intw_symmetries, only: rtau_index, rtau
    use intw_matrix_vector, only: ainv, det

    implicit none

    ! I/O
    integer, intent(in) :: isym ! Index of the symmetry operation used
    real(dp), intent(in) :: q_cryst(3) ! q point of the input dynamical matrix
    complex(dp), intent(in) :: dynq(3,3,nat,nat) ! input dynamical matrix
    complex(dp), intent(out) :: dynRq(3,3,nat,nat) ! dynamical matrix of the symmetry equivalent q point

    ! Local
    real(dp) :: q_cart(3), Rq_cart(3), s_cryst(3,3), s_cart(3,3)
    integer :: ia, ja, Sia, Sja
    complex(dp) :: phase


    ! Get the rotation matrix of the symmetry operation
    s_cryst = dble(s(:,:,isym))
    s_cart = transpose(matmul(at, matmul(transpose(s_cryst), ainv(at))))

    ! Rotate the Q point
    q_cart = matmul(bg, q_cryst)
    Rq_cart = matmul(s_cart, q_cart)

    ! Rotate the dynamical matrix
    dynRq = cmplx_0
    do ia = 1, nat
      do ja = 1, nat
        !
        Sia = rtau_index(ia,isym)
        Sja = rtau_index(ja,isym)
        !
        phase = exp(-tpi*cmplx_i*dot_product(Rq_cart, rtau(:,isym,ia) - rtau(:,isym,ja)))
        !
        dynRq(:,:,Sia,Sja) = phase * matmul(s_cart, matmul(dynq(:,:,ia,ja), ainv(s_cart)))
        !
      enddo
    enddo

    ! Apply TR if the magnetic symmetry requires it
    if (lmag .and. TR(isym)) dynRq = conjg(dynRq)

  end subroutine rotate_dyn


  subroutine mat_inv_four_t(q_cryst, nkk1, nkk2, nkk3, in_mat, out_mat)
    !
    ! This subroutine is based on the matdyn code distributed as part of
    ! the Quantum Espresso project and has been adapted to INTW:
    !
    !   Copyright (C) 2001-2004 PWSCF group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    use intw_reading, only: tau, at, bg, nat
    use intw_useful_constants, only: tpi, cmplx_0, cmplx_i, Ry_to_eV

    implicit none

    real(dp), intent(in) :: q_cryst(3)
    integer, intent(in) :: nkk1, nkk2, nkk3
    complex(dp), intent(in) :: in_mat(nkk1,nkk2,nkk3,3,3,nat,nat) ! FC matrix
    complex(dp), intent(out) :: out_mat(3*nat,3*nat) ! Dynamical matrix (without the mass factor)

    real(dp) :: q_cart(3)
    integer :: i, j
    integer :: na, nb
    real(dp) :: r(3), r_ws(3), weight, atw(3,3)
    real(dp) :: rws_aux(0:3,200)
    real(dp), allocatable :: rws(:,:)
    integer :: n1, n2, n3, m1, m2, m3, nrws
    real(dp) :: arg, total_weight
    integer, parameter :: nrwsx = 200


    q_cart = matmul(bg, q_cryst)

    out_mat(:,:) = cmplx_0

    atw(:,1) = at(:,1)*dble(nkk1)
    atw(:,2) = at(:,2)*dble(nkk2)
    atw(:,3) = at(:,3)*dble(nkk3)


    call wsinit(rws_aux, nrwsx, nrws, atw)

    allocate( rws(0:3,nrws) )

    do i = 0, 3
      do j = 1, nrws
        rws(i,j) = rws_aux(i,j)
      enddo
    enddo

    do na = 1, nat
      do nb = 1, nat
        !
        total_weight = 0.0d0
        do n1 = -16, 16 ! TODO we will need a wider range if the q-grid is very fine
          do n2 = -16, 16
            do n3 = -16, 16
              !
              do i = 1, 3
                r(i) = n1*at(i,1) + n2*at(i,2) + n3*at(i,3)
                r_ws(i) = r(i) + tau(i,na) - tau(i,nb)
              end do
              weight = wsweight(r_ws, rws, nrws)
              if (weight > 0.0) then
                !
                m1 = mod(n1 + 1, nkk1)
                if(m1 <= 0) m1 = m1 + nkk1
                m2 = mod(n2 + 1, nkk2)
                if(m2 <= 0) m2 = m2 + nkk2
                m3 = mod(n3 + 1, nkk3)
                if(m3 <= 0) m3 = m3 + nkk3
                !
                arg = tpi * (q_cart(1)*r(1) + q_cart(2)*r(2) + q_cart(3)*r(3))
                do i = 1, 3
                    do j = 1, 3
                      out_mat((na-1)*3+i, (nb-1)*3+j) = out_mat((na-1)*3+i, (nb-1)*3+j) &
                          + in_mat(m1,m2,m3, i,j, na,nb) * exp(-cmplx_i*arg) * weight
                    end do
                end do
              end if
              !
              total_weight = total_weight + weight
              !
            end do
          end do
        end do
        !
        if (abs(total_weight-nkk1*nkk2*nkk3) > 1.0d-8) then
          write(*,*) 'ERROR in mat_inv_fourier'
          write(*,*) 'total weight for R sumation:', total_weight
          stop
        end if
        !
      end do
    end do

    out_mat = 0.5d0 * (out_mat + conjg(transpose(out_mat)))

    deallocate( rws )

  end subroutine mat_inv_four_t


  subroutine read_allq_dvr()

    use intw_utility, only: find_free_unit
    use intw_reading, only: nr1, nr2, nr3, lspin, nspin, lmag, nat
    use intw_useful_constants, only: I2, sig_x, sig_y, sig_z, cmplx_0
    use intw_input_parameters, only: outdir, prefix, nqirr

    implicit none

    ! local variables

    integer :: iq, record_length, mode, ios, jspin
    character(4) :: num
    integer :: nr(3), io_unit, ir, ispin, imode, jmode
    character(256) :: dv_name


    nr(1) = nr1
    nr(2) = nr2
    nr(3) = nr3
    !
    if (lmag) then
      allocate(dvscf_irr(nr1*nr2*nr3,nqirr,3*nat,nspin**2))
    else
      allocate(dvscf_irr(nr1*nr2*nr3,nqirr,3*nat,1))
    endif
    allocate(dvscf_cart(nr1*nr2*nr3,nqirr,3*nat,nspin,nspin))
    !
    dvscf_irr = cmplx_0
    dvscf_cart = cmplx_0
    !
    inquire(iolength=record_length) dvscf_irr(1:nr1*nr2*nr3, 1, 1, 1)
    !
    do iq = 1, nqirr
      !
      if (iq <= 9) then
        write(num,'(I1.1)') iq
      elseif (iq >= 10 .and. iq <= 99) then
        write(num,'(I2.2)') iq
      elseif ( iq >= 100 .and. iq <= 999 ) then
        write(num,'(I3.3)') iq
      else
        write(num,'(I4.4)') iq
      endif
      !
      io_unit = find_free_unit()
      !
      dv_name = trim(prefix)//".dvscf_q"//trim(num)
      !
      open(unit=io_unit, file=trim(outdir)//trim(prefix)//".save.intw/"//trim(dv_name), iostat = ios, &
           form='unformatted', status='old', access='direct', recl=record_length)
      if (ios /= 0) stop "ERROR: read_allq_dvr: error opening dv_name."
      !
      write(unit=*,fmt="(a,i2,a3)") "|   Reading file "// &
            dv_name(1:max(25,len(trim(dv_name))))//" (ios ", ios, ") |"
      !
      do mode = 1, 3*nat
        if (lmag) then
          read(io_unit, rec = mode, iostat = ios) (dvscf_irr(1:nr1*nr2*nr3,iq,mode,ispin), ispin = 1, nspin**2)
        else
          read(io_unit, rec = mode, iostat = ios) dvscf_irr(1:nr1*nr2*nr3,iq,mode,1)
        endif
      enddo ! mode
      !
      ! Below not "transpose(conjg(u_irr(:,:,iq))", instead only "conjg(u_irr(:,:,iq))" because u_irr is already transpose!
      !
      if (lspin) then
        !
        if (lmag) then
          !
          do ir = 1, nr1*nr2*nr3
            do imode = 1, 3*nat
              do jmode = 1, 3*nat
                do ispin = 1, nspin
                  do jspin = 1, nspin
                    !
                    dvscf_cart(ir,iq,imode,ispin,jspin) = dvscf_cart(ir,iq,imode,ispin,jspin) + &
                        conjg(u_irr(imode,jmode,iq)) * dvscf_irr(ir,iq,jmode,1) * I2(ispin,jspin) + &
                        conjg(u_irr(imode,jmode,iq)) * dvscf_irr(ir,iq,jmode,2) * sig_x(ispin,jspin) + &
                        conjg(u_irr(imode,jmode,iq)) * dvscf_irr(ir,iq,jmode,3) * sig_y(ispin,jspin) + &
                        conjg(u_irr(imode,jmode,iq)) * dvscf_irr(ir,iq,jmode,4) * sig_z(ispin,jspin)
                  enddo ! jspin
                enddo ! ispin
              enddo ! jmode
            enddo ! imode
          enddo ! ir
          !
        else
          !
          do ir = 1, nr1*nr2*nr3
            do imode = 1, 3*nat
              do jmode = 1, 3*nat
                !
                dvscf_cart(ir,iq,imode,:,:) = dvscf_cart(ir,iq,imode,:,:) + &
                    conjg(u_irr(imode,jmode,iq)) * dvscf_irr(ir,iq,jmode,1) * I2(:,:)

                !
              enddo ! jmode
            enddo ! imode
          enddo ! ir
          !
        endif ! lmag
        !
      else
        !
        do ir = 1, nr1*nr2*nr3
          do imode = 1, 3*nat
            do jmode = 1, 3*nat
              !
              dvscf_cart(ir,iq,imode,1,1) = dvscf_cart(ir,iq,imode,1,1) + &
                                        conjg(u_irr(imode,jmode,iq)) * dvscf_irr(ir,iq,jmode,1)
              !
            enddo ! jmode
          enddo ! imode
        enddo ! ir
        !
      end if ! lspin
      !
      close(io_unit)
      !
    enddo ! iq

  end subroutine read_allq_dvr


  subroutine get_dv(qpoint_cryst, dv)

    use intw_utility, only: triple_to_joint_index_g, find_k_1BZ_and_G
    use intw_reading, only: s, nr1, nr2, nr3, nat, lspin, nspin
    use intw_useful_constants, only: I2, sig_x, sig_y, sig_z, cmplx_0
    use intw_input_parameters, only: nq1, nq2, nq3
    use intw_matrix_vector, only: cmplx_trace

    implicit none

    ! I/O variables
    real(dp), intent(in) :: qpoint_cryst(3)
    complex(dp), intent(out) :: dv(nr1*nr2*nr3,3*nat,nspin,nspin)

    ! local variables
    complex(dp) :: dvr(nr1*nr2*nr3,3*nat,nspin,nspin)
    complex(dp) :: vr(nspin,nspin)
    complex(dp) :: mrx(nspin,nspin), mry(nspin,nspin), mrz(nspin,nspin)
    integer :: i, j, k
    real(dp) :: qpoint_1bz(3), q_irr_cryst_rot(3)
    integer :: q_index, q_index_irr, s_index, imq, ir, imode
    integer :: GKQ_bz(3), GKQ(3)


    dv = cmplx_0
    !
    ! We find the associated of qpoint in 1BZ (it usually is itself!)
    !
    call find_k_1BZ_and_G(qpoint_cryst, nq1, nq2, nq3, i, j, k, qpoint_1bz, GKQ_bz)
    call triple_to_joint_index_g(nq1, nq2, nq3, q_index, i, j, k)
    !
    ! We search the irr q related with qpoint and the sym.op. which gives
    ! qpoint from irr q
    !
    q_index_irr = QE_folder_sym_q(q_index)
    s_index = symlink_q(q_index,1)
    imq = symlink_q(q_index,2)
    !
    ! We rotate irr q using the symm.op. (s_index + imq)
    !
    q_irr_cryst_rot = matmul(s(:,:,s_index), q_irr_cryst(:,q_index_irr))
    if (imq==1) q_irr_cryst_rot = -q_irr_cryst_rot
    !
    ! We define the lattice vector relating q1BZ and qrot
    !
    GKQ = nint(qpoint_cryst-q_irr_cryst_rot)
    !
    if (sum(abs(qpoint_cryst-q_irr_cryst_rot-dble(GKQ))) > 1E-4) then
      write(*,*) "ERROR: get_dv: qpoint not recovered by symmetry"
      stop
    end if
    !
    ! We rotate dV_cart of q_irr for finding dV_cart of q_rot
    !
    call rot_dvq(q_irr_cryst(1:3,q_index_irr), nr1, nr2, nr3, s_index, &
                 dvscf_cart(1:nr1*nr2*nr3,q_index_irr,1:3*nat,1:nspin,1:nspin), dvr)
    !
    dv = dvr
    !
    if (imq == 1) then ! TR YES
      !
      if (lspin) then
        !
        do ir = 1, nr1*nr2*nr3
          do imode = 1, 3*nat
            !
            vr = 0.5d0*cmplx_trace(matmul(I2, dvr(ir,imode,:,:)))
            mrx = 0.5d0*cmplx_trace(matmul(sig_x, dvr(ir,imode,:,:)))
            mry = 0.5d0*cmplx_trace(matmul(sig_y, dvr(ir,imode,:,:)))
            mrz = 0.5d0*cmplx_trace(matmul(sig_z, dvr(ir,imode,:,:)))
            !
            dv(ir,imode,:,:) = conjg(vr(:,:)) * I2 - ( conjg(mrx(:,:))*sig_x &
                                                     + conjg(mry(:,:))*sig_y &
                                                     + conjg(mrz(:,:))*sig_z)
            !
          enddo ! imode
        enddo ! ir
        !
      else
        !
        dv = conjg(dvr)
        !
      endif ! lspin
      !
    endif ! TR
    !
    ! with this we add to dV the phase for being associated to q1BZ from qrot
    !
    call func_by_gr(nr1, nr2, nr3, 3*nat, -dble(GKQ), dv)

  end subroutine get_dv


  subroutine rot_dvq(qpoint_irr_cryst, nr1, nr2, nr3, s_index, dv_in, dv_out)

    use intw_symmetries, only: rtau_index, spin_symmetry_matrices
    use intw_utility, only: triple_to_joint_index_r
    use intw_reading, only: s, ftau, nspin, lmag, at, bg, tau, nat
    use intw_useful_constants, only: cmplx_i, cmplx_0, tpi
    use intw_matrix_vector, only: cmplx_ainv, ainv

    implicit none

    ! I/O variables

    real(dp), intent(in) :: qpoint_irr_cryst(3)
    integer, intent(in) :: nr1, nr2, nr3, s_index
    complex(dp), intent(in) :: dv_in(nr1*nr2*nr3,3*nat,nspin,nspin)
    complex(dp), intent(out) :: dv_out(nr1*nr2*nr3,3*nat,nspin,nspin)

    ! local variables

    integer :: ir, ispin, jspin, i, j, k, ri, rj, rk, rri, rrj, rrk
    integer :: ipol, jpol, lpol, kpol, na, rna
    integer :: r_index, rr_index, imode
    real(dp) :: s_cart(3,3)
    real(dp) :: qpoint_irr_cart(3), qpoint_irr_cart_rot(3)
    complex(dp) :: dv_aux(nr1*nr2*nr3,3*nat,nspin,nspin), phase(nat)


    !
    ! q irr cryst -> cart
    !
    qpoint_irr_cart = matmul(bg, qpoint_irr_cryst)
    !
    ! We define our sym op.
    !
    do ipol = 1, 3
      do jpol = 1, 3
        !
        s_cart(ipol,jpol) = 0.d0
        !
        do kpol = 1, 3
          do lpol = 1, 3
            !
            s_cart(ipol,jpol) = s_cart(ipol,jpol) + &
                                at(ipol,kpol) * real(s(lpol,kpol,s_index), dp) * bg(jpol,lpol)
            !
          enddo ! lpol
        enddo ! kpol
        !
      enddo ! jpol
    enddo ! ipol
    !
    ! We rotate q irr. q rot cryst -> cart
    !
    qpoint_irr_cart_rot = matmul(ainv(s_cart), qpoint_irr_cart)
    !
    ! phase of q irr related to the atoms. We must remove this phase and add the one related
    ! to q rot
    !
    do na=1,nat
      phase(na) = exp(cmplx_i*tpi*dot_product(qpoint_irr_cart(:), tau(:,na)))
    enddo ! na
    !
    ! dV_Sq(r,na,alpha) = dV_q(Sr,Sna,Salpha)
    !
    dv_out = cmplx_0
    do ispin = 1, nspin
      do jspin = 1, nspin
        !
        do k = 1, nr3
          do j = 1, nr2
            do i = 1, nr1
              !
              ri = s(1,1,s_index)*(i-1) + s(2,1,s_index)*(j-1) + s(3,1,s_index)*(k-1) - nint(ftau(1,s_index)*nr1)
              !
              rj = s(1,2,s_index)*(i-1) + s(2,2,s_index)*(j-1) + s(3,2,s_index)*(k-1) - nint(ftau(2,s_index)*nr2)
              !
              rk = s(1,3,s_index)*(i-1) + s(2,3,s_index)*(j-1) + s(3,3,s_index)*(k-1) - nint(ftau(3,s_index)*nr3)
              !
              rri = mod(ri, nr1) + 1
              rrj = mod(rj, nr2) + 1
              rrk = mod(rk, nr3) + 1
              !
              if (rrk < 1) rrk = rrk + nr3
              if (rrj < 1) rrj = rrj + nr2
              if (rri < 1) rri = rri + nr1
              !
              call triple_to_joint_index_r(nr1, nr2, nr3, r_index, i, j, k)
              call triple_to_joint_index_r(nr1, nr2, nr3, rr_index, rri, rrj, rrk)
              !
              do na = 1, nat
                !
                rna = rtau_index(na,s_index)
                !
                do ipol = 1, 3
                  do jpol = 1, 3
                    !
                    dv_out(r_index, (na-1)*3+ipol, ispin, jspin) = dv_out(r_index, (na-1)*3+ipol, ispin, jspin) &
                                                                 + phase(rna) * s_cart(jpol,ipol) * dv_in(rr_index, (rna-1)*3+jpol, ispin, jspin)
                    !
                  enddo ! jpol
                enddo ! ipol
                !
              enddo ! na
              !
            enddo ! i
          enddo ! j
        enddo ! k
        !
      enddo ! jspin
    enddo ! ispin
    !
    do na = 1, nat
      phase(na) = exp(-cmplx_i * tpi * dot_product(qpoint_irr_cart_rot(:), tau(:,na)))
    enddo ! na
    !
    do ir = 1, nr1*nr2*nr3
      !
      do ispin = 1, nspin
        do jspin = 1, nspin
          !
          do na = 1, nat
            do ipol = 1, 3
              !
              dv_out(ir, (na-1)*3+ipol, ispin, jspin) = dv_out(ir, (na-1)*3+ipol, ispin, jspin) * phase(na)
              !
            enddo ! ipol
          enddo ! na
          !
        enddo ! jspin
      enddo ! ispin
      !
    enddo ! ir
    !
    ! Rotate the spin using spin symmetry matrices (could be done also using quaternion properties)
    !
    if (lmag) then
      !
      dv_aux = dv_out
      dv_out = cmplx_0
      !
      do ir = 1, nr1*nr2*nr3
        do imode = 1, 3*nat
          !
          dv_out(ir,imode,:,:) = matmul(cmplx_ainv(spin_symmetry_matrices(:,:,s_index)), &
                                        matmul(dv_aux(ir,imode,:,:), spin_symmetry_matrices(:,:,s_index)))
          !
        enddo ! imode
      enddo ! ir
      !
    endif ! lmag

  end subroutine rot_dvq


  subroutine func_by_gr(nr1, nr2, nr3, nmode, q, dvr)

    use intw_reading, only: nspin
    use intw_useful_constants, only: cmplx_i, cmplx_0, tpi
    use intw_utility, only: joint_to_triple_index_r

    implicit none

    ! I/O variables

    integer, intent(in) :: nr1, nr2, nr3, nmode
    complex(dp), intent(inout) :: dvr(nr1*nr2*nr3,nmode,nspin,nspin)
    real(dp), intent(in) :: q(3)

    ! local variables

    integer :: ir, i, j, k, ispin, jspin, imode
    real(dp) :: phase


    do ir = 1, nr1*nr2*nr3
      !
      call joint_to_triple_index_r(nr1, nr2, nr3, ir, i, j, k)
      !
      phase = tpi * (q(1)*real(i-1, dp)/nr1 + q(2)*real(j-1, dp)/nr2 + q(3)*real(k-1, dp)/nr3)
      !
      do imode = 1, nmode
        do ispin = 1, nspin
          do jspin = 1, nspin
            !
            dvr(ir,imode,ispin,jspin) = dvr(ir,imode,ispin,jspin) * exp(cmplx_i*phase)
            !
          enddo ! jspin
        enddo ! ispin
      enddo ! imode
      !
    enddo ! ir

  end subroutine func_by_gr


  subroutine wsinit(rws, nrwsx, nrws, atw)
    !
    ! This subroutine is originally distributed as part of the Quantum Espresso code and has
    ! been adapted to INTW:
    !   Copyright (C) 2001 PWSCF group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    implicit none

    integer :: i, ii, ir, jr, kr, nrws, nrwsx
    real(dp) :: rws(0:3,nrwsx), atw(3,3)
    integer, parameter :: nx = 2
    real(dp), parameter :: eps = 1.0d-6

    ii = 1
    do ir = -nx, nx
      do jr = -nx, nx
        do kr = -nx, nx
          !
          do i = 1, 3
            rws(i,ii) = atw(i,1)*ir + atw(i,2)*jr + atw(i,3)*kr
          end do
          !
          rws(0,ii) = rws(1,ii)*rws(1,ii) + rws(2,ii)*rws(2,ii) + rws(3,ii)*rws(3,ii)
          rws(0,ii) = 0.5d0*rws(0,ii)
          !
          if (rws(0,ii) > eps) ii = ii + 1
          if (ii > nrwsx) write(*,*) 'ERROR in wsinit'
          !
        end do
      end do
    end do
    !
    nrws = ii - 1

  end subroutine wsinit


  function wsweight(r, rws, nrws)
    !
    ! This subroutine is originally distributed as part of the Quantum Espresso code and has
    ! been adapted to INTW:
    !   Copyright (C) 2001 PWSCF group
    !   Distributed under the terms of the GNU General Public License.
    !   See the LICENSE file in the original Quantum Espresso source for license details.
    !   For the original source visit: https://www.quantum-espresso.org/
    !
    implicit none

    integer :: ir, nreq, nrws
    real(dp) :: r(3), rrt, ck, rws(0:3,nrws), wsweight
    real(dp), parameter :: eps = 1.0d-6


    wsweight = 0.d0
    nreq = 1
    do ir = 1, nrws
      rrt = r(1)*rws(1,ir) + r(2)*rws(2,ir) + r(3)*rws(3,ir)
      ck = rrt-rws(0,ir)
      if ( ck > eps ) return
      if ( abs(ck) < eps ) nreq = nreq + 1
    end do
    wsweight = 1.d0/dble(nreq)

  end function wsweight

end module intw_ph
