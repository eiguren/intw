!----------------------------------------------------------------------------!
!	intw project.
!
!----------------------------------------------------------------------------!
!
module intw_ph
  !
  use kinds, only : dp
  use intw_utility, only: find_free_unit, cmplx_trace, cmplx_ainv_2, find_k_1BZ_and_G, &
                          switch_indices, switch_indices_zyx
  use intw_reading, only: at, bg, s, ftau, nr1, nr2, nr3, nat, tau, amass, ityp, &
                          npol, spinorb_mag, tpiba, write_tag
  use intw_useful_constants, only: i2, sig_x, sig_y, sig_z
  use intw_useful_constants, only: cmplx_0, cmplx_1, cmplx_i, tpi
  use intw_input_parameters, only: mesh_dir, prefix, ph_dir, dvscf_name, qlist, fc_mat
  use intw_input_parameters, only: nq1, nq2, nq3, nqirr
  use intw_symmetries, only: rtau_index, inverse_indices, rtau, spin_symmetry_matrices
  use intw_fft, only: mill, gvec_cart


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
  public :: nmodes, q_irr, u_irr, fcmat, dvscf_cart, dvscf_irr, amass_ph, frc, nqmesh, &
            q_irr_cryst, qmesh, QE_folder_nosym_q, QE_folder_sym_q, nosym_G_q, sym_G_q, &
            symlink_q, eigqts, dvloc, dvq_local, dvpsi, asr, zeu, frc_R, dvscf_cart_comps
  !
  ! subroutines
  public :: rot_gep, read_ph_information_xml, readfc, mat_inv_four_t, read_allq_dvr, &
            calculate_local_part_dv, get_dv, rot_dvq, func_by_gr, wsinit, &
            deallocate_ph, read_fc_from_XML
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
  real(dp),allocatable       :: amass_ph(:)
  complex(dp), allocatable :: frc(:,:,:,:,:,:,:) ! force constants
  integer :: nqmesh

  real(dp),allocatable::  q_irr_cryst(:,:)
  real(dp),allocatable::  qmesh(:,:)

  integer, allocatable ::  QE_folder_nosym_q(:)
  integer, allocatable ::  QE_folder_sym_q(:)
  integer, allocatable ::  nosym_G_q      (:,:)
  integer, allocatable ::  sym_G_q      (:,:)
  integer, allocatable ::  symlink_q      (:,:)

  complex(dp), allocatable :: eigqts(:)

  complex(dp), allocatable :: dvloc (:)

  complex(dp), allocatable    :: dvq_local    (:,:,:,:)
  complex(dp), allocatable    :: dvpsi (:,:,:,:,:)

!  integer :: fc_ibrav
  CHARACTER(LEN=10)  :: asr
  REAL(DP), ALLOCATABLE ::  zeu(:,:,:), frc_R(:,:,:,:,:,:,:)
  complex(dp),allocatable    :: dvscf_cart_comps(:,:,:,:)

  !
contains

  subroutine rot_gep(  s_index, imq, qpoint_irr, nbnd, npol, nat, g_matin, g_matout  )
    implicit none
    !input
    real(kind=dp), intent(in) :: qpoint_irr(3)
    integer      , intent(in) :: imq, nbnd, npol, nat

    complex(kind=dp), intent(in)  :: g_matin ( nbnd, nbnd, npol, npol, 3*nat)
    complex(kind=dp), intent(out) :: g_matout( nbnd, nbnd, npol, npol, 3*nat)
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
                     g_matout ( ibnd, jbnd, 1:npol, 1:npol, (rna-1)*3 + ipol) =  &
                     g_matout ( ibnd, jbnd, 1:npol, 1:npol, (rna-1)*3 + ipol)    &
                           + s_cart(ipol,jpol)* g_matin ( ibnd, jbnd, 1:npol, 1:npol, (na-1)*3 + jpol) * phase(rna)
                   enddo !jpol
                enddo !ipol
             enddo
          enddo

       enddo !na

    if  (imq==1) g_matout =conjg(g_matout)

    return
  end subroutine rot_gep

  function rot_k_index(s_index, k_index, nk1, nk2, nk3, kmesh )

    implicit none

    integer, intent(in) :: s_index, k_index, nk1, nk2, nk3
    real(dp), intent(in) :: kmesh(3,nk1*nk2*nk3)

    integer :: rot_k_index


    real(dp) :: kpoint_1bz(3), kpoint_rot(3)
    integer :: GKQ_bz(3), i, j, k

    kpoint_rot = matmul(s(:,:,s_index), kmesh(:,k_index))

    call find_k_1BZ_and_G(kpoint_rot,nk1,nk2,nk3,i ,j, k, kpoint_1bz, GKQ_bz)

    call switch_indices (nk1, nk2, nk3, rot_k_index, i, j, k, +1)

  end function rot_k_index

  subroutine read_ph_information_xml()
    !------------------------------------------------------------------
    !
    !This subroutine reads the information related to phonons from xml files.
    !
    !------------------------------------------------------------------
    implicit none

    integer        ::  io_unit  ! input/output file unit

    character(256) :: datafile ! full path of the data-file.xml file
    ! in the .xml file

    integer        :: imode, jmode, iq
    !
    character(len=8)  :: dummy


    nmodes =3 * nat

    write(*,"(4(a,i4))     ")"PH BZ division is             : ", nq1," x",nq2," x",nq3
    write(*,"(a,i4)"     ) "Input n. of irreducible points: ", nqirr


    ! Read irreducible q list
    allocate(q_irr(3,nqirr))

    io_unit = find_free_unit()
    open(unit=io_unit, file=trim(trim(mesh_dir)//trim(ph_dir)//trim(qlist)),status="old")

    do iq=1,nqirr
      read(io_unit,*) dummy, q_irr(1:3,iq)
    enddo

    close(unit=io_unit)


    io_unit = find_free_unit()
    datafile=trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//"irrq_patterns.dat")
    open(unit=io_unit, file=datafile,status="old")

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

    use intw_set_asr, only: set_asr

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
    integer :: io_unit

    !  !
    !  !
    io_unit = find_free_unit()
    OPEN (unit=io_unit,file=flfrc,status='old',form='formatted')
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
                READ (io_unit,*) (((m1bid, m2bid, m3bid,             &
                     frc_R(m1,m2,m3,i,j,na,nb),                  &
                     m1=1,nr1),m2=1,nr2),m3=1,nr3)
!                DO m3=1,nr3
!                   DO m2=1,nr2
!                      DO m1=1,nr3
!                         READ (1,*) m1bid, m2bid, m3bid, frc(m1,m2,m3,i,j,na,nb)
!                      ENDDO
!                   ENDDO
!                ENDDO
                !
             END DO
          END DO
       END DO
    END DO
    !  !
    CLOSE(unit=io_unit)
    frc=frc_R*cmplx_1
    !  !
       !IGG
       ! =========================================
       ! Set Acoustic Sum rule -- **  crystal*** --
       ! ==========================================
       ! Apply ASR as in QE
       ! Beware. In this module (ph_module) frc is defined as complex.
       ! But in subroutine set_asr it is real(dp).
       ! And then in mat_inv_four_t it is complex again.
       asr='crystal'
       ALLOCATE (zeu(3,3,nat))
       zeu =0.d0
       write(*,*)
       write(*,*) 'Applying ASR : ' ,asr

       call set_asr(asr,nq1,nq2,nq3,frc_R,zeu,nat,ibrav,tau_fc)

       frc=frc_R*cmplx_1                   ! make it complex again
       deallocate (zeu)
       deallocate(frc_R)
       ! End IGG

!! IGG-k komentatuta behekoa
    !ASR aplikatu
!    do i=1,3
!       do j=1,3
!          do na=1,nat
!             batura=0.0d0
!             do nb=1,nat
!                do m1=1,nq1
!                   do m2=1,nq2
!                      do m3=1,nq3
!                         batura=batura+ frc(m1,m2,m3, i ,j, na, nb)
!                      end do
!                   end do
!                end do
!             end do
!             frc(1,1,1, i ,j, na, na) = frc(1,1,1,i ,j, na, na) - batura
!          end do
!       end do
!    end do
!! IGG-k honaino komentatuta

    RETURN
  END SUBROUTINE readfc


  subroutine mat_inv_four_t ( q_point, nkk1, nkk2, nkk3,  nnmode , frc , out_mat )

    USE intw_reading, only : tau, ityp
    use intw_useful_constants, only: tpi

    implicit none

    complex(dp), intent(in) :: frc(nkk1, nkk2, nkk3, 3,3, nat,nat)
    complex(dp), intent(out) :: out_mat(nnmode, nnmode)

    real(dp) :: q_point(3),q(3)
    integer :: nkk1, nkk2, nkk3, nnmode
    integer :: i, j
    complex(dp), parameter :: II=(0.d0,1.d0)
    integer :: na,nb
    real(dp) :: r(3), r_ws(3), weight, atw(3,3)
    real(dp) :: rws_aux(0:3,200)
    real(dp), allocatable :: rws(:,:)
    integer :: n1, n2, n3, m1, m2, m3, nrws
    real(dp) :: arg, total_weight
    integer, parameter :: nrwsx=200
! IGG : aumev zuzenduta from 27211.6d0  to  27211.396d0
    real(dp), parameter :: pmass=1822.88848426_dp, aumev=  27211.396d0

    out_mat(:,:)=(0.d0,0.d0)

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
    !call cryst_to_cart (1, q , bg, 1)

    do na=1, nat
       do nb=1, nat
          total_weight=0.0d0
          do n1=-16,16
             do n2=-16,16
                do n3=-16,16
                   !
                   do i=1, 3
                      r(i) =    n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
                      r_ws(i) = r(i) + tau (i,na)-tau (i,nb)
                   end do
                   weight = wsweight(r_ws,rws,nrws)
                   !if (weight.gt.0.001) write(*,'(a,3i4,f12.6)')'weight',n1,n2,n3,weight
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
                                 frc(m1,m2,m3, i,j, na, nb) * exp(-II*arg)*weight / &
                                 sqrt(amass(ityp(na))*amass(ityp(nb))) / &
                                 (pmass/2.d0) * & !Up to this in Ry.
                                 (aumev/2.d0)**2 !Now in meV.
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
  subroutine   read_allq_dvr(nqirr,nmodes)
!-----------------------------------------------------------------------

    implicit none

    !I/O varibales

    integer, intent (in) :: nqirr, nmodes

    !local variables

    integer :: iq, record_length, mode, ios, is, js
    character(len=20) ::  num
    integer :: nr(3), io_unit, ir, ipol, imode, jmode
    character(len=256) :: dv_name


    if (nmodes/=3*nat) then
       write(*,"(a)")"ERROR(read_allq_dvr): nmodes = ", nmodes, " and 3*nat =", 3*nat
       stop
    endif
    !
    nr(1)=nr1
    nr(2)=nr2
    nr(3)=nr3
    !
    allocate(dvscf_irr(nr1*nr2*nr3,nqirr,3*nat,1:npol**2),dvscf_cart(nr1*nr2*nr3,nqirr,3*nat,1:npol,1:npol))
    !
    dvscf_irr=cmplx_0
    dvscf_cart=cmplx_0
    !
    if (spinorb_mag) then
      inquire(iolength=record_length) dvscf_irr(1:nr1*nr2*nr3, 1, 1,1:npol**2)
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
       dv_name=trim(trim(mesh_dir)//trim(prefix)//".save.intw/"//trim(prefix)//"."//trim(dvscf_name)//trim(num))
       !
       open (io_unit, file = dv_name, iostat = ios, &
            form = 'unformatted', status = 'old', access = 'direct', recl= record_length )
!
       write(unit=*,fmt="(2a,x,a,i4)")"Reading irreducible displacement dvscr file ", &
            trim(dv_name),",ios: ",ios
       do mode=1, 3*nat
          if (spinorb_mag) then
             read (io_unit, rec = mode, iostat = ios) (dvscf_irr(1:nr1*nr2*nr3, iq, mode,ipol),ipol=1,npol**2)
          else
             read (io_unit, rec = mode, iostat = ios) dvscf_irr(1:nr1*nr2*nr3, iq, mode,1) ! beste guztiak 0
          endif
       enddo !mode
       !
       ! Below not "transpose(conjg(u_irr(:,:,iq))", instead only "conjg(u_irr(:,:,iq))" because u_irr is already transpose !!
       !

       if (npol==2) then
          !
          if (spinorb_mag) then
             !
             do ir=1,nr1*nr2*nr3
                do imode=1,3*nat
                   do jmode=1,3*nat
                      do is=1,npol
                         do js=1,npol
                           !
                            dvscf_cart(ir,iq,imode,is,js)=dvscf_cart(ir,iq,imode,is,js) + &
                               conjg(u_irr(imode,jmode,iq))*dvscf_irr(ir,iq,jmode,1)*I2(is,js) + &
                               conjg(u_irr(imode,jmode,iq))*dvscf_irr(ir,iq,jmode,2)*sig_x(is,js) + &
                               conjg(u_irr(imode,jmode,iq))*dvscf_irr(ir,iq,jmode,3)*sig_y(is,js) + &
                               conjg(u_irr(imode,jmode,iq))*dvscf_irr(ir,iq,jmode,4)*sig_z(is,js)
                         enddo !js
                      enddo !is
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
       end if !npol==2
       !
       close(io_unit)
       !
    enddo !iq
!stop
    !
    return
    !
  end subroutine read_allq_dvr
!*********************************************************************************
!---------------------------------------------------------------------------------
  subroutine calculate_local_part_dv(qpoint, nat, npol, dvq_local )
!---------------------------------------------------------------------------------
!
!======================================================================
! We have dV_scf as input and we add to it the derivative of the PP   !
!======================================================================

    USE intw_reading, ONLY : nr1, nr2, nr3, ngm
    USE intw_fft, ONLY : eigts1, eigts2, eigts3, nl
    USE intw_pseudo, ONLY : vlocq

    implicit none

    external :: cfftnd

    !I/O variables

    real(dp),intent(in) :: qpoint(1:3) ! crystal coord.
    integer,intent(in) :: nat,npol
    complex(dp),intent(inout) :: dvq_local(nr1*nr2*nr3,3*nat,npol,npol) ! spin idependentea da baina koherentzia mantenduko dugu.

    !local variables

    integer :: imode,na,ipol,ig,nt,is,ir
    complex(dp) :: aux (nr1*nr2*nr3), fact, gtau
    real(dp) :: qcart(3) ! qpoint in cart.

    qcart=matmul(bg,qpoint)
    !
    do imode=1,3*nat
       !
       na=(imode-1)/3+1
       ipol=modulo(imode-1,3)+1
       nt=ityp (na)
       !
       aux(:)=cmplx_0
       !
       fact=-tpiba*cmplx_i*eigqts(na)
       !
       do ig=1,ngm
          !
          gtau=eigts1(mill(1,ig),na)*eigts2(mill(2,ig),na)*eigts3(mill(3,ig),na)
          aux(nl(ig))=aux(nl(ig))+vlocq(ig,nt)*(qcart(ipol)+gvec_cart(ipol,ig))*gtau*fact
          !
       enddo !ig
       !
       call cfftnd(3,(/nr1,nr2,nr3/),1,aux)
       !
       do ir=1,nr1*nr2*nr3
          do is=1,npol
             !
             dvq_local(ir,imode,is,is)=dvq_local(ir,imode,is,is)+aux(ir)
             !
          enddo !is
       enddo !ir
       !
    enddo !imode osagai kanonikoetan zehar goaz hemen ..!
    !
    return
    !
  end subroutine calculate_local_part_dv
!*******************************************************************
!-------------------------------------------------------------------
  subroutine get_dv(iq,qpoint,nmode,npol,dvq_local)
!-------------------------------------------------------------------

    implicit none

    !I/O variables
    integer, intent(in) :: iq, nmode, npol
    real(dp), intent(in) :: qpoint(1:3) ! in crystal
    complex(dp), intent(out) :: dvq_local(nr1*nr2*nr3,nmode,npol,npol)

    !local variables
    complex(dp) :: dvr(nr1*nr2*nr3,nmode,npol,npol)
    complex(dp) :: vr(npol,npol)
    complex(dp) :: mrx(npol,npol),mry(npol,npol),mrz(npol,npol)
    integer :: i, j, k
    real(dp) :: qpoint_1bz(1:3), qpoint_rot(1:3)
    integer :: q_index, q_index_irr, s_index, imq, ir, imode
    integer :: GKQ_bz(1:3), GKQ (1:3)


    dvq_local=cmplx_0
    !
    ! We find the associated of qpoint in 1BZ (it usually is itself!)
    !
    call find_k_1BZ_and_G(qpoint,nq1,nq2,nq3,i,j,k,qpoint_1bz,GKQ_bz)
    call switch_indices(nq1,nq2,nq3,q_index,i,j,k,+1)
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
         (/0,0,0/),dvscf_cart(1:nr1*nr2*nr3,q_index_irr,1:nmode,1:npol,1:npol),dvr)
    !
    dvq_local = dvr
    !
    if (imq.eq.1) then !TR YES
       !
       if (npol==1) then !SOC NO
          !
          dvq_local=conjg(dvr)
          !
       else if (npol==2) then !SOC YES
          !
          do ir=1,nr1*nr2*nr3
             do imode=1,nmode
                !
                vr=0.5d0*cmplx_trace(matmul(I2,dvr(ir,imode,:,:)))
                mrx=0.5d0*cmplx_trace(matmul(sig_x,dvr(ir,imode,:,:)))
                mry=0.5d0*cmplx_trace(matmul(sig_y,dvr(ir,imode,:,:)))
                mrz=0.5d0*cmplx_trace(matmul(sig_z,dvr(ir,imode,:,:)))
                !
                dvq_local(ir,imode,:,:)=conjg(vr(:,:))*I2-(conjg(mrx(:,:))*sig_x+conjg(mry(:,:))*sig_y+conjg(mrz(:,:))*sig_z)
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
    call func_by_gr(nr1,nr2,nr3,nmode,-dble(GKQ),dvq_local)
    !
    return

  end subroutine get_dv
!--------------------------------------------------------------------------------------------------------
  subroutine rot_dvq(q_point_crys,q_point_crys_irr,nr1,nr2,nr3,nmode,s_index,GKQ,dv_in,dv_out)
!--------------------------------------------------------------------------------------------------------

    implicit none

    !I/O variables

    real(dp), intent(in) :: q_point_crys(3), q_point_crys_irr(3) ! crystal
    integer, intent(in) :: nr1, nr2, nr3, nmode, s_index ,GKQ(3)
    complex(dp), intent(in) :: dv_in(nr1*nr2*nr3,nmode,npol,npol)
    complex(dp), intent(out) :: dv_out(nr1*nr2*nr3,nmode,npol,npol)

    !local variables

    integer :: ir, is, js, i, j, k, ri, rj, rk, rri, rrj, rrk
    integer :: ipol, jpol, lpol, kpol, na, rna
    integer :: r_index, rr_index, imode
    real(dp) :: s_cart(3,3), s_crys(3,3)
    real(dp) :: q_point(3), q_point_r(3) ! cart
    complex(dp) :: dv_aux(nr1*nr2*nr3,nmode,npol,npol), phase(nat)

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
    dv_out = (0.d0, 0.d0)
    do is=1,npol
       do js=1,npol
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
                   call switch_indices_zyx(nr1,nr2,nr3,r_index,i,j,k,+1)
                   call switch_indices_zyx(nr1,nr2,nr3,rr_index,rri,rrj,rrk,+1)
                   !
                   do na = 1, nat
                      !
                      rna=rtau_index(na,s_index)
                      !
                      do ipol=1,3
                         do jpol=1,3
                            !
                            dv_out(r_index,(na-1)*3+ipol,is,js)=dv_out(r_index,(na-1)*3+ipol,is,js)    &
                                                               +s_cart(jpol,ipol)*dv_in(rr_index,(rna-1)*3+jpol,is,js)*phase(rna)
                            !
                         enddo !jpol
                      enddo !ipol
                   enddo !na
                enddo !i
             enddo !j
          enddo !k
          !
       enddo !js
    enddo !is
    !
    do na=1,nat
       phase(na)=exp(-cmplx_i*(q_point_r(1)*tau(1,na)+q_point_r(2)*tau(2,na)+q_point_r(3)*tau(3,na))*tpi)
    enddo !na
    !
    do ir=1,nr1*nr2*nr3
       do is=1,npol
          do js=1,npol
             !
             do na=1,nat
                do ipol=1,3
                   !
                   dv_out(ir,(na-1)*3+ipol,is,js)=dv_out(ir,(na-1)*3+ipol,is,js)*phase(na)
                   !
                enddo !ipol
             enddo !na
             !
          enddo !js
       enddo !is
    enddo !ir
    !
    ! Spinarekin errotatu behar da. Hemen spin matrizeen bidez egina, egin daiteke quaternioi propietateak erabiliz ere.
    !
    if (npol==2.and.spinorb_mag) then
       !
       dv_aux=dv_out
       dv_out=cmplx_0
       !
       do ir=1,nr1*nr2*nr3
          do imode=1,3*nat
             !
             dv_out(ir,imode,:,:)=matmul(cmplx_ainv_2(spin_symmetry_matrices(:,:,s_index)), &
                  matmul(dv_aux(ir,imode,:,:),spin_symmetry_matrices(:,:,s_index)))
             !
          enddo !imode
       enddo !ir
       !
    endif !npol=2
    !
    return

  end subroutine rot_dvq
!***********************************************************************
!-----------------------------------------------------------------------
  subroutine func_by_gr(nr1,nr2,nr3,nmode,q,dvr)
!-----------------------------------------------------------------------

    implicit none

    !I/O variables

    integer,intent(in) :: nr1,nr2,nr3,nmode
    complex(dp),intent(inout) :: dvr(nr1*nr2*nr3,nmode,npol,npol)
    real(dp),intent(in) :: q(3)

    !local variables

    integer :: ir, i, j, k, is, js, imode
    real(dp) :: phase

    do ir=1,nr1*nr2*nr3
       !
       call switch_indices_zyx(nr1,nr2,nr3,ir,i,j,k,-1)
       !
       phase=tpi*(q(1)*real(i-1,dp)/nr1 + q(2)*real(j-1,dp)/nr2 + q(3)*real(k-1,dp)/nr3)
       !
       do imode=1,nmode
          do is=1,npol
             do js=1,npol
                !
                dvr(ir,imode,is,js)=dvr(ir,imode,is,js)*exp(cmplx_i*phase)
                !
             enddo !js
          enddo !is
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
!
! ----------------------------------------------------
!IGG

! ----------------------------------------------------
  subroutine read_fc_from_XML ()

    use intw_set_asr, only: set_asr

    implicit none

    ! local variables
    integer num_atom
    real(dp), allocatable :: mass_atom(:), pos_atom(:), masa(:)
    integer :: io_unit
    integer :: i, n
    integer :: ntyp, meshq(3), fc_ibrav
    character(256) :: fc_file
    character(256) :: type_string
    INTEGER :: na, nb, nn, m1, m2, m3
    REAL(dp) :: aux(3,3)
    CHARACTER(LEN=10) :: asr
    REAL(DP), ALLOCATABLE :: zeu(:,:,:), frc_R(:,:,:,:,:,:,:)

    real(dp), parameter :: pmass = 1822.88848426_dp, aumev = 27211.6d0

    fc_file = trim(mesh_dir)//trim(ph_dir)//trim(fc_mat)
    io_unit = find_free_unit()
    print*, 'Reading from fc file:', trim(fc_file)
    write(*,*)
    print*, 'checking consistency with the electronic calculation'
    ! Check consistency with the electronic calculation
    !call iotk_scan_begin (io_unit,"GEOMETRY_INFO")
    !  call iotk_scan_dat (io_unit,"BRAVAIS_LATTICE_INDEX",fc_ibrav)
    !IGG : I am not sure ibrav is read anywhere else (?)
    !if (fc_ibrav.ne.ibrav) then
    !print*, 'WARNING!, ibrav mismatch', ibrav,fc_ibrav
    !end if
    print*, 'Bravais lattice type (from FC file)', fc_ibrav
    if (num_atom.ne.nat) then
      print*, 'WARNING!, number of atoms do not match', num_atom, nat
      print*,'STOPPING'
      STOP
    end if

    allocate(masa(ntyp))
    do i = 1,ntyp
      call write_tag("MASS.",i, type_string)
    end do

    allocate(mass_atom(nat), pos_atom(nat))
    do i=1,nat
      write(*,*) 'Atom: ' ,i
      call write_tag("ATOM.",i, type_string)
      mass_atom(i)=masa(n)
      !write(*,'(2(a7,i3,1x),a10,2f16.4)')  &
      !    'Atom :',i, '; Type:',n,'; mass : ', mass_atom(i), &
      !                     mass_atom(i)/(pmass/2)
      if ((abs(amass(ityp(i))-mass_atom(i)/(pmass/2))).gt.(1.d-01)) then
        print*, '! ! ! W A R N I N G ! ! ! ! '
        print*, 'Masses do not match for atom', i
        write(*,'(a20, 3f16.4)') 'Mass data.xml file', amass(ityp(i))
        write(*,'(a20, 3(f16.4,2x))') 'Mass fc.xml file', mass_atom(i)/(pmass/2)
      else
        print*, 'Masses consistent in data.xml and fc.xml files'
      end if
      !write(*,'(2(a7,i3,1x),a10,3f16.4)')  &
      !    'Atom :',i, '; Type:',n,'position', pos_atom(:)
      if ( ( abs(tau(1,i) - pos_atom(1)).gt.1d-04) .or. &
           ( abs(tau(2,i) - pos_atom(2)).gt.1d-04) .or. &
           ( abs(tau(3,i) - pos_atom(3)).gt.1d-04) ) then
        print*, '! ! ! W A R N I N G ! ! ! ! '
        print*, 'Atomic positions do not match for atom', i
      else
        print*, 'Atomic positions consistent in data.xml and fc.xml files'
      end if
    end do !nat
    deallocate(masa)
    deallocate(mass_atom, pos_atom)

    if ( (meshq(1).ne.nq1).or. &
         (meshq(2).ne.nq2).or. &
         (meshq(3).ne.nq3) ) then
      print*, '! ! ! W A R N I N G ! ! ! ! '
      print*, 'q-mesh not consistent. FCfile:', meshq
      print*, '             Input file:', nq1,nq2,nq3
    else
      print*, 'q-mesh consistent in input and FC file:',  meshq
    end if

    print*, 'Reading FC matrices'
    !Zuzenean QEko read_ifc subrutinatik
    !allocate( phid(nq1*nq2*nq3,3,3,nat,nat) )
    allocate( frc(nq1,nq2,nq3,3,3,nat,nat) )
    !phid(:,:,:,:,:)=0d0
    frc(:,:,:,:,:,:,:)=0d0
    DO na=1,nat
      DO nb=1,nat
        nn=0
        DO m3=1,nq3
          DO m2=1,nq2
            DO m1=1,nq1
              nn=nn+1
              frc(m1,m2,m3,:,:,na,nb) = aux(:,:)
            ENDDO
          ENDDO
        ENDDO
        print*,na,nb,nn
      ENDDO
    ENDDO

    ! =========================================
    ! Set Acoustic Sum rule -- **  crystal*** --
    ! ==========================================
    ! Apply ASR as in QE
    ! Beware. In this module (ph_module) frc is defined as complex.
    ! But in subroutine set_asr it is real(dp).
    ! And then in mat_inv_four_t it is complex again.

    asr='crystal'
    ALLOCATE (zeu(3,3,nat))
    allocate(frc_R(nq1,nq2,nq3,3,3,nat,nat))
    zeu =0.d0
    frc_R =real(frc,dp)         ! make it real
    write(*,*)
    write(*,*) 'Applying ASR : ' ,asr

    call set_asr(asr,nq1,nq2,nq3,frc_R,zeu,nat,fc_ibrav,tau)

    frc=frc_R*cmplx_1           !make it complex again
    deallocate (zeu)
    deallocate(frc_R)

    return
  end subroutine read_fc_from_XML

!End IGG

end module intw_ph
