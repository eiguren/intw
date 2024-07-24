!----------------------------------------------------------------------------!
!	intw project.
!
!----------------------------------------------------------------------------!
!
module intw_ph
  !
  use kinds, only : dp

  !----------------------------------------------------------------------------!
  !
  !  This module contains the main variables related to phonon modes.
  !
  !----------------------------------------------------------------------------!

  !
  implicit none
  !
  save
  ! variables
  public :: nmodes, q_irr, u_irr, dvscf_cart, dvscf_irr, frc, nqmesh, &
            q_irr_cryst, qmesh, QE_folder_nosym_q, QE_folder_sym_q, nosym_G_q, sym_G_q, &
            symlink_q, dvq_local, dvpsi, zeu, frc_R
  !
  ! subroutines
  public :: rot_gep, read_ph_information_xml, readfc, mat_inv_four_t, read_allq_dvr, &
            get_dv, rot_dvq, func_by_gr, wsinit, deallocate_ph, &
            read_dynq, set_asr_frc,  set_asr_dynq
  !
  ! functions
  public :: wsweight, rot_k_index
  !
  private
  !
  !
  integer                    :: nmodes               ! 3 x nat
  real(dp),allocatable       :: q_irr(:,:)           !coo. of irr. q points
  complex(dp),allocatable    :: u_irr(:,:,:)         !displacement patterns for the irr.  q.
  complex(dp),allocatable    :: fcmat(:,:,:)         !dymanical matrices for irr. q.
  complex(dp),allocatable    :: dvscf_cart(:,:,:,:,:), dvscf_irr(:,:,:,:)
  complex(dp), allocatable :: frc(:,:,:,:,:,:,:) ! force constants
  integer :: nqmesh

  real(dp),allocatable::  q_irr_cryst(:,:)
  real(dp),allocatable::  qmesh(:,:)

  integer, allocatable ::  QE_folder_nosym_q(:)
  integer, allocatable ::  QE_folder_sym_q(:)
  integer, allocatable ::  nosym_G_q      (:,:)
  integer, allocatable ::  sym_G_q      (:,:)
  integer, allocatable ::  symlink_q      (:,:)

  complex(dp), allocatable    :: dvq_local    (:,:,:,:)
  complex(dp), allocatable    :: dvpsi (:,:,:,:,:)

  REAL(DP), ALLOCATABLE ::  zeu(:,:,:), frc_R(:,:,:,:,:,:,:)

  !
contains

  subroutine rot_gep(  s_index, imq, qpoint_irr, nbnd, nspin, nat, g_matin, g_matout  )

    use intw_symmetries, only: rtau_index, inverse_indices, rtau
    use intw_reading, only: bg, s, at
    use intw_useful_constants, only: tpi, cmplx_0, cmplx_i

    implicit none
    !input
    real(kind=dp), intent(in) :: qpoint_irr(3)
    integer      , intent(in) :: imq, nbnd, nspin, nat

    complex(kind=dp), intent(in)  :: g_matin ( nbnd, nbnd, nspin, nspin, 3*nat)
    complex(kind=dp), intent(out) :: g_matout( nbnd, nbnd, nspin, nspin, 3*nat)
    !loca
    integer :: s_index, s_inv_index, na
    integer :: rna
    complex (dp) :: phase(nat)
    integer :: kpol, jpol, lpol, ipol, ibnd, jbnd
    real(dp) :: s_cart(3,3), s_crys(3,3), qpoint_irr_cart (3)

    qpoint_irr_cart = matmul(bg, qpoint_irr)


    g_matout=cmplx_0

    s_inv_index=inverse_indices(s_index)  !

    do na=1,nat
       phase(na) = &
            exp(  - cmplx_i * ( qpoint_irr_cart (1) * rtau (1, s_index, na)            &
                              + qpoint_irr_cart (2) * rtau (2, s_index, na)            &
                              + qpoint_irr_cart (3) * rtau (3, s_index, na) ) * tpi)
    enddo

    do ipol=1,3
       do jpol=1,3
          s_crys(ipol,jpol) = real(s (ipol,jpol,s_index),dp)
       enddo
    enddo


    do ipol = 1, 3
       do jpol = 1, 3
          s_cart (ipol, jpol) = 0.d0
          do kpol = 1, 3
             do lpol = 1, 3
                s_cart (ipol, jpol) = s_cart (ipol, jpol) + at (ipol, kpol) * &
                     s_crys (lpol, kpol)  * bg (jpol, lpol)
             enddo
          enddo
       enddo
    enddo


       do na = 1, nat

          rna= rtau_index(na,s_index)

          do ibnd=1,nbnd
             do jbnd=1, nbnd
                do ipol=1,3
                   do jpol=1,3
                     g_matout ( ibnd, jbnd, 1:nspin, 1:nspin, (rna-1)*3 + ipol) =  &
                     g_matout ( ibnd, jbnd, 1:nspin, 1:nspin, (rna-1)*3 + ipol)    &
                           + s_cart(ipol,jpol)* g_matin ( ibnd, jbnd, 1:nspin, 1:nspin, (na-1)*3 + jpol) * phase(rna)
                   enddo !jpol
                enddo !ipol
             enddo
          enddo

       enddo !na

    if  (imq==1) g_matout =conjg(g_matout)

    return
  end subroutine rot_gep

  function rot_k_index(s_index, k_index, nk1, nk2, nk3, kmesh )

    use intw_reading, only: s
    use intw_utility, only: find_k_1BZ_and_G, triple_to_joint_index_g

    implicit none

    integer, intent(in) :: s_index, k_index, nk1, nk2, nk3
    real(dp), intent(in) :: kmesh(3,nk1*nk2*nk3)

    integer :: rot_k_index


    real(dp) :: kpoint_1bz(3), kpoint_rot(3)
    integer :: GKQ_bz(3), i, j, k

    kpoint_rot = matmul(s(:,:,s_index), kmesh(:,k_index))

    call find_k_1BZ_and_G(kpoint_rot,nk1,nk2,nk3,i ,j, k, kpoint_1bz, GKQ_bz)

    call triple_to_joint_index_g(nk1, nk2, nk3, rot_k_index, i, j, k)

  end function rot_k_index

  subroutine read_ph_information_xml()
    !------------------------------------------------------------------
    !
    !This subroutine reads the information related to phonons from xml files.
    !
    !------------------------------------------------------------------
    use intw_utility, only: find_free_unit
    use intw_reading, only: nat
    use intw_input_parameters, only: mesh_dir, prefix, ph_dir, qlist, nqirr, nq1, nq2, nq3

    implicit none

    integer ::  io_unit, ios

    character(256) :: datafile ! full path of the data-file.xml file
    ! in the .xml file

    integer :: imode, jmode, iq
    !
    character(len=8) :: dummy


    nmodes = 3 * nat

    write(*,"(4(a,i4))     ")"PH BZ division is             : ", nq1," x",nq2," x",nq3
    write(*,"(a,i4)"     ) "Input n. of irreducible points: ", nqirr


    ! Read irreducible q list
    allocate(q_irr(3,nqirr))

    io_unit = find_free_unit()
    open(unit=io_unit, file=trim(mesh_dir)//trim(ph_dir)//trim(qlist), status="old", iostat=ios)
    if (ios /= 0) stop "ERROR: read_ph_information_xml: error opening qlist."

    do iq=1,nqirr
      read(io_unit,*) dummy, q_irr(1:3,iq)
    enddo

    close(unit=io_unit)


    io_unit = find_free_unit()
    datafile = trim(mesh_dir)//trim(prefix)//".save.intw/"//"irrq_patterns.dat"
    open(unit=io_unit, file=datafile, status="old", iostat=ios)
    if (ios /= 0) stop "ERROR: read_ph_information_xml: error opening irrq_patterns."

    !Read displacement patters for each q.
    allocate (u_irr(nmodes,nmodes,nqirr))

    write(*,"(a)") "Reading displacement patterns: "
    write(*,"(a)") datafile
    do iq=1, nqirr
       read(unit=io_unit,fmt=*)dummy
       do imode=1,3*nat
        read(unit=io_unit,fmt=*) ( u_irr(jmode,imode,iq), jmode=1,3*nat )
       end do
    end do !iq

    close(unit=io_unit)

  end subroutine read_ph_information_xml

  SUBROUTINE readfc ( flfrc, nr1, nr2, nr3, nat, alat, at, ntyp, amass )

    !use intw_set_asr, only: set_asr
    use intw_utility, only: find_free_unit
    use intw_useful_constants, only: cmplx_1
    use intw_input_parameters, only: nq1, nq2, nq3, apply_asr

    IMPLICIT NONE

    CHARACTER(LEN=256) :: flfrc

    INTEGER :: ibrav, nr1, nr2, nr3, nat, ntyp
    REAL(DP) :: alat, at(3,3), epsil(3,3)
    LOGICAL :: has_zstar
    !  ! local variables
    INTEGER :: i, j, na, nb, m1, m2, m3
    INTEGER :: ibid, jbid, nabid, nbbid, m1bid,m2bid,m3bid
    REAL(DP) :: amass(ntyp), amass_from_file, celldm(6)
    INTEGER :: nt
    CHARACTER(LEN=3) :: atm
    CHARACTER(LEN=9) :: symm_type
    real(dp) :: tau_fc(3, nat), zeu_fc(3,3,nat)
    integer  :: ityp_fc(nat)
    integer :: io_unit, ios

    !  !
    !  !
    io_unit = find_free_unit()
    OPEN(unit=io_unit, file=flfrc, status='old', form='formatted', iostat=ios)
    IF ( ios /= 0 ) STOP 'ERROR: readfc: error opening flfrc'
    !  !
    !  !  read cell data
    !  !
    READ(io_unit,*) ntyp,nat,ibrav,(celldm(i),i=1,6)
    if (ibrav==0) then
       read(io_unit,'(a)') symm_type
       read(io_unit,*) ((at(i,j),i=1,3),j=1,3)
    end if
    !

    alat = celldm(1)
    !  !
    !  !  read atomic types, positions and masses
    !  !
    DO nt = 1,ntyp
       READ(io_unit,*) i,atm,amass_from_file
    END DO
    !  !
    !  !
    DO na=1,nat
       READ(io_unit,*) i,ityp_fc(na),(tau_fc(j,na),j=1,3)
    END DO
    !  !
    !  !  read macroscopic variable
    !  !
    READ (io_unit,*) has_zstar
    IF (has_zstar) THEN
       READ(io_unit,*) ((epsil(i,j),j=1,3),i=1,3)
       DO na=1,nat
          READ(io_unit,*)
          READ(io_unit,*) ((zeu_fc(i,j,na),j=1,3),i=1,3)
       END DO
    ENDIF
    !  !
    READ (io_unit,*) nr1,nr2,nr3
    !  !
    !  !  read real-space interatomic force constants
    !  !
    ALLOCATE ( frc(nr1,nr2,nr3,3,3,nat,nat) )
    frc(:,:,:,:,:,:,:) = 0.d0
    allocate(frc_R(nq1,nq2,nq3,3,3,nat,nat))
    frc_R(:,:,:,:,:,:,:) = 0.d0
    DO i=1,3
       DO j=1,3
          DO na=1,nat
             DO nb=1,nat
                READ (io_unit,*) ibid, jbid, nabid, nbbid
                READ (io_unit,*) (((m1bid, m2bid, m3bid, frc_R(m1,m2,m3,i,j,na,nb), &
                                    m1=1,nr1),m2=1,nr2),m3=1,nr3)
             END DO
          END DO
       END DO
    END DO
    !  !
    CLOSE(unit=io_unit)
    frc=frc_R*cmplx_1
       ALLOCATE (zeu(3,3,nat))
       zeu =0.d0
       !MBR 21/05/2024
       ! Instead of that set_asr module, use only the "simple" method
       ! (worked in the routine as real)
       if (apply_asr) then
          write(*,*)' Applying ASR (simple) to force constant matrix'
          call set_asr_frc(nq1, nq2, nq3, nat, frc_R)
       end if   

       frc=frc_R*cmplx_1                   ! make it complex again
       deallocate (zeu)
       deallocate(frc_R)
       ! End IGG

  END SUBROUTINE readfc


  !====================================================================
  ! MBR 21/05/24
  ! Apply ASR "simple" method:
  !   set_asr_frc --> to the real-space force constant matrix 
  !                   (real elements)
  !   set_asr_dynr --> to the dynamical matrices after doing dynq_to_dynr
  !====================================================================

   subroutine set_asr_frc(n1, n2, n3, nat, frc_R)

   implicit none

   integer , intent(in) :: n1, n2, n3, nat 
   real(dp) , intent(inout) :: frc_R(n1,n2,n3,3,3,nat,nat)

   ! Internal
   integer :: i,j, iat, jat, i1,i2,i3
   real(dp) :: suma

   do i=1,3  !coords.
   do j=1,3  
        do iat=1,nat  !atom
           suma=0.0_dp
           do jat=1,nat  !atom
              do i1=1,n1  !lattice vectors
              do i2=1,n2
              do i3=1,n3
                   suma=suma+frc_R(i1,i2,i3,i,j,iat,jat)
              end do
              end do
              end do !lattice vectors
           end do !jat atom
           frc_R(1,1,1,i,j,iat,iat) = frc_R(1,1,1,i,j,iat,iat) - suma
        end do !iat atom
   end do 
   end do !coords.

   return 
   end subroutine set_asr_frc         



   subroutine set_asr_dynq(nq, iq, nat, dynq)
           ! iq is the q=0 index in the mesh

   use intw_useful_constants, only: cmplx_1
   implicit none

   integer , intent(in) :: nq, iq, nat
   complex(dp) , intent(inout) :: dynq(3*nat,3*nat, nq)

   ! Internal
   integer :: i,j, iat, jat,ir
   complex(dp) :: suma
   complex(dp) :: dynq_aux(3*nat,3*nat, nq)

   dynq_aux = dynq
      do i=1,3  !coords.
      do j=1,3
          do iat=1,nat  !atom
             suma=0.0_dp
             do jat=1,nat  !atom !sum over cells
                 suma=suma + dynq( (iat-1)*3+i, (jat-1)*3+j, iq )  
             end do !jat atom
             ! option 1 : apply only to q=0:
             ! dynq_aux( (iat-1)*3+i, (iat-1)*3+j, iq ) = dynq( (iat-1)*3+i, (iat-1)*3+j, iq ) - suma
             ! option 2 : generalize to all q:
             dynq_aux( (iat-1)*3+i, (iat-1)*3+j, : ) = dynq( (iat-1)*3+i, (iat-1)*3+j, : ) - suma
          end do !iat atom
      end do
      end do !coords.
   dynq = dynq_aux

   return
   end subroutine set_asr_dynq


  !====================================================================
      ! MBR 14/02/24
      ! read dynamical matrices info from prefix.dyn_q(iq) files in
      ! prefix.save.intw and store into dynq for all q-points
      ! (each files contains the info of an irreducible q-point)
  !====================================================================

      subroutine read_dynq (dynq)

      use intw_reading, only: ntyp, nat, at, amass, ityp
      use intw_useful_constants, only: eps_6, cmplx_0, cmplx_i, cmplx_1
      use intw_utility, only: cryst_to_cart, find_free_unit, &
              triple_to_joint_index_g, find_k_1BZ_and_G
      use intw_input_parameters, only: mesh_dir, prefix, nqirr, nq1, nq2, nq3, &
              apply_asr

      implicit none

      ! I/O
      complex(dp) , intent(out) :: dynq(3*nat,3*nat,nqmesh)

      ! Internal
      character(len=4) :: iq_loc
      character(len=100) :: comentario
      character(256) :: dynq_file, intwdir
      logical :: all_done
      integer :: ntyp_, nat_, ibrav_
      integer :: dynq_unit, ierr
      integer :: iat, iat1, iat2, iq1, iq2, iq, iqirr, Gq(3), i,j,k, iq_done(nqmesh)
      real(dp) :: celldm_(6), at_(3,3), dist, massfac
      real(dp) ::  fracpos(3), qpoint(3), qpoint1(3), qpoint_cart(3,nqmesh)
      real(dp) ::  dynq_re(3,3), dynq_im(3,3)
      ! From mat_inv_four_t, see IGG's comments
      ! TODO add to useful_constants
      real(dp), parameter :: pmass=1822.88848426_dp, aumev=  27211.396d0

      iq_done=0

      ! INTW directory
      intwdir = trim(mesh_dir)//trim(prefix)//".save.intw/"

      do iqirr = 1,nqirr

      ! open file
      if (                   iqirr <   10) write(iq_loc,"(i1)") iqirr
      if ( 10 <= iqirr .and. iqirr <  100) write(iq_loc,"(i2)") iqirr
      if (100 <= iqirr .and. iqirr < 1000) write(iq_loc,"(i3)") iqirr
      dynq_unit = find_free_unit()
      dynq_file = trim(intwdir)//trim(prefix)//".dyn_q"//adjustl(iq_loc)
      open(unit=dynq_unit, iostat=ierr, file=dynq_file, form='formatted', status='old')
      if (ierr /= 0 ) then
         write(*,*) 'Error opening .dyn_q file in ',  intwdir,' . Stopping.'
         stop
      end if

      read(dynq_unit,'(a)') comentario
      read(dynq_unit,'(a)') comentario

      ! Information about lattice and atoms (not used but checked for
      ! consistency in the parameters)
      read(dynq_unit,*) ntyp_, nat_, ibrav_, celldm_
      if (ntyp_ /= ntyp .or.  nat_ /= nat ) then
         write(*,*)' Error reading parameters in .dyn_q file ', dynq_file, '. Stopping.'
         stop
      end if
      if ( ibrav_ == 0 ) then
        read(dynq_unit,'(a)') comentario  ! symm_type  not used
        read(dynq_unit,*) ((at_(i,j),i=1,3),j=1,3)  ! not used
      end if
      ! atoms
      do i=1,ntyp
         read(dynq_unit,'(a)') comentario !i, atom, amass, not used, not stored
      end do
      ! positions
      do iat = 1,nat
        read(dynq_unit,*) i,j, fracpos(:)  !not used, not stored
      end do

      ! loop over symmetry equivalent qpoints:
      ! a priori, we do not know how many there are.
      ! The signal to stop is that the last qpoint line is the
      ! one for diagonalization information, where the initial
      ! q-point is stated again.

      do iq1 = 1,nqmesh

         read(dynq_unit,'(/,a,/)') comentario
         read(dynq_unit,'(11x, 3(f14.9), / )') qpoint_cart(:,iq1)

         ! compare to first read qpoint to see if we are at the
         ! end of the part of the file that contains the dynq matrices
         dist = ( qpoint_cart(1,iq1) - qpoint_cart(1,1) )**2 + &
                ( qpoint_cart(2,iq1) - qpoint_cart(2,1) )**2 + &
                ( qpoint_cart(3,iq1) - qpoint_cart(3,1) )**2
         if (iq1 > 1 .and. dist < eps_6) exit

         !cartesian to crystalline
         qpoint=qpoint_cart(:,iq1)
         call cryst_to_cart (1, qpoint, at, -1)
         ! find index in full list: iq
         call find_k_1BZ_and_G(qpoint,nq1,nq2,nq3,i,j,k,qpoint1,Gq)
         call triple_to_joint_index_g(nq1,nq2,nq3,iq,i,j,k)

         ! block: dynq matrix
         ! loop over atoms
         do iat1 = 1,nat
           do iat2 = 1,nat

              ! mass factor sqrt(m1*m2) to be used below
              ! (QE routine elph.f90 adds it when reading dynq and then
              ! again after diagonalization)
              massfac = sqrt( amass(ityp(iat1)) * amass(ityp(iat2))  )

              read(dynq_unit,*) i,j ! not used

              ! cartesian 3x3 block of this atom pair in dynq matrix
              read(dynq_unit,*) ( (dynq_re(i,j), dynq_im(i,j), j=1,3),i=1,3)

              ! Add mass factor and unit conversion
              ! (this will give the same as mat_inv_four_t)
              dynq( (iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3, iq ) = &
                     ( dynq_re*cmplx_1 + dynq_im*cmplx_i ) &
          !           / massfac &  !sqrt(m1*m2)
                     / (pmass/2.d0) & !Up to this in Ry.
                     * (aumev/2.d0)**2 !Now in meV.


           end do
         end do !atoms

         iq_done(iq)=1

      end do !iq1

      ! read next block in .dyn file: diagonalization of dynq
      ! (skipped, as we usually do this step elsewhere in INTW)

      close(dynq_unit)

      end do ! iqirr
            ! Check that all the qpoints in the full mesh have been read
      all_done=.true.
      do iq=1,nqmesh
         if (iq_done(iq) == 0) then
            all_done=.false.
            write(*,*)' Failed to read dynq for iq= ',iq
         end if
      end do
      if ( .not. all_done ) stop

      ! Apply ASR to q=0 matrix:
      ! \sum_{jat} dynq( iat, i, jat, j, q=0) = 0
      if (apply_asr) then 
           write(*,*)' Applying ASR to all q vector indices'
           !iq = 1 !q=0 index in mesh
           ! just in case, find q=0 index in mesh 
           qpoint = (/ 0._dp, 0._dp, 0._dp /)
           call find_k_1BZ_and_G(qpoint,nq1,nq2,nq3,i,j,k,qpoint1,Gq)
           call triple_to_joint_index_g(nq1,nq2,nq3,iq,i,j,k)
           call set_asr_dynq(nqmesh, iq, nat, dynq)
      end if        

      !add mass factor
      do iat1 = 1,nat
           do iat2 = 1,nat
              massfac = sqrt( amass(ityp(iat1)) * amass(ityp(iat2))  )
              dynq( (iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3, : ) = &
                      dynq( (iat1-1)*3+1:iat1*3, (iat2-1)*3+1:iat2*3, : ) / massfac  !sqrt(m1*m2)
           end do   
      end do   

      return
      end subroutine read_dynq
  !====================================================================
      


  subroutine mat_inv_four_t(q_point, nkk1, nkk2, nkk3, nnmode, in_mat, out_mat)

    USE intw_reading, only : tau, ityp, at, bg, amass, nat
    use intw_useful_constants, only: tpi, cmplx_0, cmplx_i, Ry_to_eV

    implicit none

    complex(dp), intent(in) :: in_mat(nkk1,nkk2,nkk3,3,3,nat,nat)
    complex(dp), intent(out) :: out_mat(nnmode, nnmode)

    real(dp) :: q_point(3),q(3)
    integer :: nkk1, nkk2, nkk3, nnmode
    integer :: i, j
    integer :: na,nb
    real(dp) :: r(3), r_ws(3), weight, atw(3,3)
    real(dp) :: rws_aux(0:3,200)
    real(dp), allocatable :: rws(:,:)
    integer :: n1, n2, n3, m1, m2, m3, nrws
    real(dp) :: arg, total_weight
    integer, parameter :: nrwsx=200
    real(dp), parameter :: pmass = 1822.88848426_dp

    out_mat(:,:)=cmplx_0

    atw(:,1) = at(:,1)*dble(nkk1)
    atw(:,2) = at(:,2)*dble(nkk2)
    atw(:,3) = at(:,3)*dble(nkk3)


    call wsinit(rws_aux,nrwsx,nrws,atw)

    allocate( rws(0:3, nrws) )

    do i=0,3
       do j=1,nrws
          rws(i,j)= rws_aux(i,j)
       enddo
    enddo

    q=matmul(bg,q_point) ! q_point kristaletan etorri behar da.

    do na=1, nat
       do nb=1, nat
          total_weight=0.0d0
          do n1=-16,16   ! TODO we will need a wider range if the q-grid is very fine
             do n2=-16,16
                do n3=-16,16
                   !
                   do i=1, 3
                      r(i) =    n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
                      r_ws(i) = r(i) + tau (i,na)-tau (i,nb)
                   end do
                   weight = wsweight(r_ws,rws,nrws)
                   if (weight .gt. 0.0) then
                      !
                      m1 = mod(n1+1,nkk1)
                      if(m1.le.0) m1=m1+nkk1
                      m2 = mod(n2+1,nkk2)
                      if(m2.le.0) m2=m2+nkk2
                      m3 = mod(n3+1,nkk3)
                      if(m3.le.0) m3=m3+nkk3
                      !
                      arg = tpi*(q(1)*r(1) + q(2)*r(2) + q(3)*r(3))
                      do i=1, 3
                         do j=1, 3
                            out_mat((na-1)*3+i, (nb-1)*3+j)= &
                                 out_mat((na-1)*3+i, (nb-1)*3+j) + &
                                 in_mat(m1,m2,m3, i,j, na,nb) * exp(-cmplx_i*arg)*weight / &
                                 sqrt(amass(ityp(na))*amass(ityp(nb))) / &
                                 (pmass/2.d0) * & !Up to this in Ry.
                                 (Ry_to_eV*1000)**2 !Now in meV.
                         end do
                      end do
                   end if
                   total_weight=total_weight + weight
                end do
             end do
          end do
          if (abs(total_weight-nkk1*nkk2*nkk3).gt.1.0d-8) then
             write(*,*)'ERROR in mat_inv_fourier'
             write(*,*) 'total weight for R sumation:',total_weight
             stop
          end if
       end do
    end do

    out_mat = 0.5d0 * (out_mat + conjg(transpose(out_mat)))

    deallocate( rws )


  end subroutine mat_inv_four_t

!***********************************************************************
!-----------------------------------------------------------------------
  subroutine read_allq_dvr()
!-----------------------------------------------------------------------

    use intw_utility, only: find_free_unit
    use intw_reading, only: nr1, nr2, nr3, nspin, spinorb_mag, nat
    use intw_useful_constants, only: I2, sig_x, sig_y, sig_z, cmplx_0
    use intw_input_parameters, only: mesh_dir, prefix, dvscf_name, nqirr

    implicit none

    !local variables

    integer :: iq, record_length, mode, ios, jspin
    character(len=4) ::  num
    integer :: nr(3), io_unit, ir, ispin, imode, jmode
    character(len=256) :: dv_name


    !
    nr(1)=nr1
    nr(2)=nr2
    nr(3)=nr3
    !
    allocate(dvscf_irr(nr1*nr2*nr3,nqirr,3*nat,1:nspin**2),dvscf_cart(nr1*nr2*nr3,nqirr,3*nat,1:nspin,1:nspin))
    !
    dvscf_irr=cmplx_0
    dvscf_cart=cmplx_0
    !
    if (spinorb_mag) then
      inquire(iolength=record_length) dvscf_irr(1:nr1*nr2*nr3, 1, 1,1:nspin**2)
    else
      inquire(iolength=record_length) dvscf_irr(1:nr1*nr2*nr3, 1, 1,1)
    endif
    !
    do iq=1,nqirr
       !
       if (iq.le.9) then
          write(num,'(I1.1)')iq!+4
       elseif (iq.ge.10.and.iq.le.99) then
          write(num,'(I2.2)')iq!+4
       elseif ( iq.ge.100.and.iq.le.999 ) then
          write(num,'(I3.3)')iq!+4
       else
          write(num,'(I4.4)')iq!+4
       endif
       !
       io_unit=find_free_unit()
       !
       dv_name=trim(mesh_dir)//trim(prefix)//".save.intw/"//trim(prefix)//"."//trim(dvscf_name)//trim(num)
       !
       open(unit=io_unit, file=dv_name, iostat = ios, &
            form='unformatted', status='old', access='direct', recl=record_length)
       if (ios /= 0) stop "ERROR: read_allq_dvr: error opening dv_name."
!
       write(unit=*,fmt="(2a,x,a,i4)")"Reading irreducible displacement dvscr file ", &
            trim(dv_name),",ios: ",ios
       do mode=1, 3*nat
          if (spinorb_mag) then
             read (io_unit, rec = mode, iostat = ios) (dvscf_irr(1:nr1*nr2*nr3, iq, mode,ispin),ispin=1,nspin**2)
          else
             read (io_unit, rec = mode, iostat = ios) dvscf_irr(1:nr1*nr2*nr3, iq, mode,1) ! beste guztiak 0
          endif
       enddo !mode
       !
       ! Below not "transpose(conjg(u_irr(:,:,iq))", instead only "conjg(u_irr(:,:,iq))" because u_irr is already transpose !!
       !

       if (nspin==2) then
          !
          if (spinorb_mag) then
             !
             do ir=1,nr1*nr2*nr3
                do imode=1,3*nat
                   do jmode=1,3*nat
                      do ispin=1,nspin
                         do jspin=1,nspin
                           !
                            dvscf_cart(ir,iq,imode,ispin,jspin)=dvscf_cart(ir,iq,imode,ispin,jspin) + &
                               conjg(u_irr(imode,jmode,iq))*dvscf_irr(ir,iq,jmode,1)*I2(ispin,jspin) + &
                               conjg(u_irr(imode,jmode,iq))*dvscf_irr(ir,iq,jmode,2)*sig_x(ispin,jspin) + &
                               conjg(u_irr(imode,jmode,iq))*dvscf_irr(ir,iq,jmode,3)*sig_y(ispin,jspin) + &
                               conjg(u_irr(imode,jmode,iq))*dvscf_irr(ir,iq,jmode,4)*sig_z(ispin,jspin)
                         enddo !jspin
                      enddo !ispin
                   enddo !jmode
                enddo !imode
             enddo !ir
             !
          else
             !
             do ir=1,nr1*nr2*nr3
                do imode=1,3*nat
                   do jmode=1,3*nat
                      !
                      dvscf_cart(ir,iq,imode,:,:)=dvscf_cart(ir,iq,imode,:,:) + &
                         conjg(u_irr(imode,jmode,iq))*dvscf_irr(ir,iq,jmode,1)*I2(:,:)

                      !
                   enddo !jmode
                enddo !imode
             enddo !ir
             !
          endif !spinorb_mag
          !
       else
          !
          do ir=1,nr1*nr2*nr3
             do imode=1,3*nat
                do jmode=1,3*nat
                   !
                   dvscf_cart(ir,iq,imode,1,1)=dvscf_cart(ir,iq,imode,1,1) + &
                                              conjg(u_irr(imode,jmode,iq))*dvscf_irr(ir,iq,jmode,1)
                   !
                enddo !jmode
             enddo !imode
          enddo !ir
          !
       end if !nspin==2
       !
       close(io_unit)
       !
    enddo !iq
!stop
    !
    return
    !
  end subroutine read_allq_dvr
!*******************************************************************
!-------------------------------------------------------------------
  subroutine get_dv(qpoint, nmode, nspin, dv)
!-------------------------------------------------------------------

    use intw_utility, only: cmplx_trace, triple_to_joint_index_g, find_k_1BZ_and_G
    use intw_reading, only: s, nr1, nr2,nr3
    use intw_useful_constants, only: I2, sig_x, sig_y, sig_z, cmplx_0
    use intw_input_parameters, only: nq1, nq2, nq3

    implicit none

    !I/O variables
    integer, intent(in) :: nmode, nspin
    real(dp), intent(in) :: qpoint(1:3) ! in crystal
    complex(dp), intent(out) :: dv(nr1*nr2*nr3,nmode,nspin,nspin)

    !local variables
    complex(dp) :: dvr(nr1*nr2*nr3,nmode,nspin,nspin)
    complex(dp) :: vr(nspin,nspin)
    complex(dp) :: mrx(nspin,nspin),mry(nspin,nspin),mrz(nspin,nspin)
    integer :: i, j, k
    real(dp) :: qpoint_1bz(1:3), qpoint_rot(1:3)
    integer :: q_index, q_index_irr, s_index, imq, ir, imode
    integer :: GKQ_bz(1:3), GKQ (1:3)


    dv = cmplx_0
    !
    ! We find the associated of qpoint in 1BZ (it usually is itself!)
    !
    call find_k_1BZ_and_G(qpoint,nq1,nq2,nq3,i,j,k,qpoint_1bz,GKQ_bz)
    call triple_to_joint_index_g(nq1,nq2,nq3,q_index,i,j,k)
    !
    ! We search the irr q related with qpoint and the sym.op. which gives
    ! qpoint from irr q
    !
    q_index_irr=QE_folder_sym_q(q_index)
    s_index=symlink_q(q_index,1)
    imq=symlink_q(q_index,2)
    !
    ! We rotate irr q using the symm.op. (s_index + imq)
    !
    qpoint_rot=matmul(s(:,:,s_index),q_irr_cryst(:,q_index_irr))
    if (imq==1) qpoint_rot=-qpoint_rot
    !
    ! We define the lattice vector relating q1BZ and qrot
    !
    GKQ=nint(qpoint-qpoint_rot)
    !
    if (sum(abs(qpoint-qpoint_rot-dble(GKQ)))>1E-4) then
       write(*,*)"ERROREA get_dv: qpoint not recovered by symmetry"
       stop
    end if
    !
    ! We rotate dV_cart of q_irr for finding dV_cart of q_rot
    !
    call rot_dvq(qpoint,q_irr_cryst(1:3,q_index_irr),nr1,nr2,nr3,nmode,s_index, &
         (/0,0,0/),dvscf_cart(1:nr1*nr2*nr3,q_index_irr,1:nmode,1:nspin,1:nspin),dvr)
    !
    dv = dvr
    !
    if (imq.eq.1) then !TR YES
       !
       if (nspin==1) then !SOC NO
          !
          dv = conjg(dvr)
          !
       else if (nspin==2) then !SOC YES
          !
          do ir=1,nr1*nr2*nr3
             do imode=1,nmode
                !
                vr=0.5d0*cmplx_trace(matmul(I2,dvr(ir,imode,:,:)))
                mrx=0.5d0*cmplx_trace(matmul(sig_x,dvr(ir,imode,:,:)))
                mry=0.5d0*cmplx_trace(matmul(sig_y,dvr(ir,imode,:,:)))
                mrz=0.5d0*cmplx_trace(matmul(sig_z,dvr(ir,imode,:,:)))
                !
                dv(ir,imode,:,:) = conjg(vr(:,:))*I2-(conjg(mrx(:,:))*sig_x+conjg(mry(:,:))*sig_y+conjg(mrz(:,:))*sig_z)
                !
             enddo !imode
          enddo !ir
          !
       endif !SOC
       !
    endif !TR
    !
    ! with this we add to dV the phase for being associated to q1BZ from qrot
    !
    call func_by_gr(nr1, nr2, nr3, nmode, -dble(GKQ), dv)
    !
    return

  end subroutine get_dv
!--------------------------------------------------------------------------------------------------------
  subroutine rot_dvq(q_point_crys,q_point_crys_irr,nr1,nr2,nr3,nmode,s_index,GKQ,dv_in,dv_out)
!--------------------------------------------------------------------------------------------------------
    use intw_symmetries, only: rtau_index, spin_symmetry_matrices
    use intw_utility, only: cmplx_ainv, triple_to_joint_index_r
    use intw_reading, only: s, ftau, nspin, spinorb_mag, at, bg, tau, nat
    use intw_useful_constants, only: cmplx_i, cmplx_0, tpi

    implicit none

    !I/O variables

    real(dp), intent(in) :: q_point_crys(3), q_point_crys_irr(3) ! crystal
    integer, intent(in) :: nr1, nr2, nr3, nmode, s_index ,GKQ(3)
    complex(dp), intent(in) :: dv_in(nr1*nr2*nr3,nmode,nspin,nspin)
    complex(dp), intent(out) :: dv_out(nr1*nr2*nr3,nmode,nspin,nspin)

    !local variables

    integer :: ir, ispin, jspin, i, j, k, ri, rj, rk, rri, rrj, rrk
    integer :: ipol, jpol, lpol, kpol, na, rna
    integer :: r_index, rr_index, imode
    real(dp) :: s_cart(3,3), s_crys(3,3)
    real(dp) :: q_point(3), q_point_r(3) ! cart
    complex(dp) :: dv_aux(nr1*nr2*nr3,nmode,nspin,nspin), phase(nat)

    ! q irr cryst -> cart
    !
    q_point = matmul(bg, q_point_crys_irr)
    !
    ! We define our sym op. label
    !
    do ipol=1,3
       do jpol=1,3
          s_crys(ipol,jpol) = real(s(ipol,jpol,s_index),dp)
       enddo
    enddo
    !
    ! We rotate q irr. q rot cryst -> cart
    !
    q_point_r = matmul(s(:,:,s_index),q_point_crys_irr (:))
    q_point_r = matmul(bg,q_point_r)
    !
    ! phase of q irr related to the atoms. We must remove this phase and add the one related
    ! to q rot
    !
    do na=1,nat
       phase(na)=exp(cmplx_i*(q_point(1)*tau(1,na)+q_point(2)*tau(2,na)+q_point(3)*tau(3,na))*tpi)
    enddo !na
    !
    ! Sym op. cryst -> cart (as QE do)
    !
    do ipol=1,3
       do jpol=1,3
          !
          s_cart(ipol,jpol) = 0.d0
          !
          do kpol=1,3
             do lpol=1,3
                !
                s_cart(ipol,jpol) = s_cart(ipol,jpol) + &
                                   at(ipol,kpol) * s_crys(lpol,kpol) * bg(jpol,lpol)
                !
             enddo !lpol
          enddo !kpol
       enddo !jpol
    enddo !ipol
    !
    ! dV_Sq(r,na,alpha)=dV_q(Sr,Sna,Salpha)
    !
    dv_out = cmplx_0
    do ispin=1,nspin
       do jspin=1,nspin
          !
          do k = 1, nr3
             do j = 1, nr2
                do i = 1, nr1
                   !
                   ri=s(1,1,s_index)*(i-1)+s(2,1,s_index)*(j-1) &
                     +s(3,1,s_index)*(k-1)-nint(ftau(1,s_index)*nr1)
                   !
                   rj=s(1,2,s_index)*(i-1)+s(2,2,s_index)*(j-1) &
                     +s(3,2,s_index)*(k-1)-nint(ftau(2,s_index)*nr2)
                   !
                   rk=s(1,3,s_index)*(i-1)+s(2,3,s_index)*(j-1) &
                     +s(3,3,s_index)*(k-1)-nint(ftau(3,s_index)*nr3)
                   !
                   rri=mod(ri,nr1)+1
                   rrj=mod(rj,nr2)+1
                   rrk=mod(rk,nr3)+1
                   !
                   if (rrk < 1) rrk=rrk+nr3
                   if (rrj < 1) rrj=rrj+nr2
                   if (rri < 1) rri=rri+nr1
                   !
                   call triple_to_joint_index_r(nr1,nr2,nr3,r_index,i,j,k)
                   call triple_to_joint_index_r(nr1,nr2,nr3,rr_index,rri,rrj,rrk)
                   !
                   do na = 1, nat
                      !
                      rna=rtau_index(na,s_index)
                      !
                      do ipol=1,3
                         do jpol=1,3
                            !
                            dv_out(r_index,(na-1)*3+ipol,ispin,jspin)=dv_out(r_index,(na-1)*3+ipol,ispin,jspin)    &
                                                               +s_cart(jpol,ipol)*dv_in(rr_index,(rna-1)*3+jpol,ispin,jspin)*phase(rna)
                            !
                         enddo !jpol
                      enddo !ipol
                   enddo !na
                enddo !i
             enddo !j
          enddo !k
          !
       enddo !jspin
    enddo !ispin
    !
    do na=1,nat
       phase(na)=exp(-cmplx_i*(q_point_r(1)*tau(1,na)+q_point_r(2)*tau(2,na)+q_point_r(3)*tau(3,na))*tpi)
    enddo !na
    !
    do ir=1,nr1*nr2*nr3
       do ispin=1,nspin
          do jspin=1,nspin
             !
             do na=1,nat
                do ipol=1,3
                   !
                   dv_out(ir,(na-1)*3+ipol,ispin,jspin)=dv_out(ir,(na-1)*3+ipol,ispin,jspin)*phase(na)
                   !
                enddo !ipol
             enddo !na
             !
          enddo !jspin
       enddo !ispin
    enddo !ir
    !
    ! Spinarekin errotatu behar da. Hemen spin matrizeen bidez egina, egin daiteke quaternioi propietateak erabiliz ere.
    !
    if (nspin==2.and.spinorb_mag) then
       !
       dv_aux=dv_out
       dv_out=cmplx_0
       !
       do ir=1,nr1*nr2*nr3
          do imode=1,3*nat
             !
             dv_out(ir,imode,:,:)=matmul(cmplx_ainv(spin_symmetry_matrices(:,:,s_index)), &
                  matmul(dv_aux(ir,imode,:,:),spin_symmetry_matrices(:,:,s_index)))
             !
          enddo !imode
       enddo !ir
       !
    endif !nspin=2
    !
    return

  end subroutine rot_dvq
!***********************************************************************
!-----------------------------------------------------------------------
  subroutine func_by_gr(nr1,nr2,nr3,nmode,q,dvr)
!-----------------------------------------------------------------------

    use intw_reading, only: nspin
    use intw_useful_constants, only: cmplx_i, cmplx_0, tpi
    use intw_utility, only: joint_to_triple_index_r

    implicit none

    !I/O variables

    integer,intent(in) :: nr1,nr2,nr3,nmode
    complex(dp),intent(inout) :: dvr(nr1*nr2*nr3,nmode,nspin,nspin)
    real(dp),intent(in) :: q(3)

    !local variables

    integer :: ir, i, j, k, ispin, jspin, imode
    real(dp) :: phase

    do ir=1,nr1*nr2*nr3
       !
       call joint_to_triple_index_r(nr1,nr2,nr3,ir,i,j,k)
       !
       phase=tpi*(q(1)*real(i-1,dp)/nr1 + q(2)*real(j-1,dp)/nr2 + q(3)*real(k-1,dp)/nr3)
       !
       do imode=1,nmode
          do ispin=1,nspin
             do jspin=1,nspin
                !
                dvr(ir,imode,ispin,jspin)=dvr(ir,imode,ispin,jspin)*exp(cmplx_i*phase)
                !
             enddo !jspin
          enddo !ispin
       enddo !imode
       !
    enddo !ir
    !
    return

  end subroutine func_by_gr
!************************************************************************
  !-----------------------------------------------------------------------
  subroutine wsinit(rws,nrwsx,nrws,atw)
    !-----------------------------------------------------------------------
    !
    implicit none
    integer :: i, ii, ir, jr, kr, nrws, nrwsx, nx
    real(dp) :: eps, rws(0:3,nrwsx), atw(3,3)
    parameter (eps=1.0d-6,nx=2)

    ii = 1
    do ir=-nx,nx
       do jr=-nx,nx
          do kr=-nx,nx
             do i=1,3
                rws(i,ii) = atw(i,1)*ir + atw(i,2)*jr + atw(i,3)*kr
             end do
             rws(0,ii)=rws(1,ii)*rws(1,ii)+rws(2,ii)*rws(2,ii)+            &
                  rws(3,ii)*rws(3,ii)
             rws(0,ii)=0.5d0*rws(0,ii)
             if (rws(0,ii).gt.eps) ii = ii + 1
             if (ii.gt.nrwsx) write(*,*) 'ERROR in wsinit'!call errore('wsinit', 'ii.gt.nrwsx',1)
          end do
       end do
    end do
    nrws = ii - 1
    return
  end subroutine wsinit
  !
  !-----------------------------------------------------------------------
  function wsweight(r,rws,nrws)
    !-----------------------------------------------------------------------
    !
    implicit none
    integer :: ir, nreq, nrws
    real(dp) :: r(3), rrt, ck, eps, rws(0:3,nrws), wsweight
    parameter (eps=1.0d-6)
    !
    wsweight = 0.d0
    nreq = 1
    do ir =1,nrws
       rrt = r(1)*rws(1,ir) + r(2)*rws(2,ir) + r(3)*rws(3,ir)
       ck = rrt-rws(0,ir)
       if ( ck .gt. eps ) return
       if ( abs(ck) .lt. eps ) nreq = nreq + 1
    end do
    wsweight = 1.d0/dble(nreq)
    return
  end function wsweight

  subroutine deallocate_ph()
    deallocate(q_irr)
    deallocate(u_irr)
  end subroutine deallocate_ph

end module intw_ph
