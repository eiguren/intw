module intw_pseudo
  !
  use kinds, only: dp
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! variables
  public :: INTWPSEUDO, upf, vlocq
  public :: nqxq, nqx, dq, qrad, tab, spline_ps, tab_d2y, npsx, nh, nhm, &
            nbetam, lmaxkb, lmaxx, nkb, indv, nhtol, nhtolm, ijtoh, vkb, vkqb, &
            nhtoj, DKB, beta
  !
  ! subroutines
  public :: read_all_pseudo, init_KB_projectors, init_pp, phq_init
  !
  private
  !
  TYPE INTWPSEUDO
    CHARACTER(LEN=2) :: psd=' '      ! Element label
    CHARACTER(LEN=20) :: typ=' '     ! Pseudo type ( NC or US or PAW)
    CHARACTER(len=6) :: rel=' '      ! relativistic: {no|scalar|full}
    LOGICAL :: nlcc                  ! Non linear core corrections
    REAL(DP) :: zp                   ! z valence
    REAL(DP) :: ecutwfc              ! suggested cut-off for wfc
    !
    INTEGER :: lmax                  ! maximum l component in beta
    INTEGER :: nbeta                 ! number of projectors
    INTEGER, DIMENSION(:), ALLOCATABLE :: kbeta   ! kbeta(nbeta): number of grid points used for betas
                                    ! this defines the cutoff radius for each of them.
    !
    INTEGER :: kkbeta
    INTEGER, DIMENSION(:), ALLOCATABLE :: lll     ! lll(nbeta) l of each projector
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: beta   ! beta(mesh,nbeta) projectors
    !
    INTEGER :: mesh                  ! number of points in the radial mesh
    REAL(DP), DIMENSION(:), ALLOCATABLE :: r      ! r(mesh)  radial grid
    REAL(DP), DIMENSION(:), ALLOCATABLE :: rab    ! rab(mesh) dr(x)/dx (x=linear grid)
    INTEGER :: lloc                 ! L of channel used to generate local potential
    REAL(DP) :: rcloc               ! vloc = v_ae for r > rcloc
    REAL(DP), DIMENSION(:), ALLOCATABLE :: vloc    ! vloc(mesh) local atomic potential
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dion  ! dion(nbeta,nbeta) atomic D_{mu,nu}

    LOGICAL :: has_so             ! if .true. includes spin-orbit
    REAL(DP), DIMENSION(:), ALLOCATABLE :: jjj   ! jjj(nbeta) j=l+1/2 or l-1/2 of beta
  END TYPE INTWPSEUDO

  TYPE (INTWPSEUDO), DIMENSION(:), ALLOCATABLE :: UPF

  REAL (DP), ALLOCATABLE :: vlocq(:,:)

  LOGICAL :: spline_ps=.true.

  !Former US in QE
  INTEGER :: &
       nqxq,             &! size of interpolation table
       nqx                ! number of interpolation points
  REAL(DP), PARAMETER:: &
       dq = 0.01D0           ! space between points in the pseudopotential tab.
  REAL(DP), ALLOCATABLE :: &
       qrad(:,:,:,:),         &! radial FT of Q functions
       tab(:,:,:)           ! interpolation table for PPs

  REAL(DP), ALLOCATABLE :: &
       tab_d2y(:,:,:)            ! for cubic splines
  !Former USPP_PARAM in QE (INTW)
  !TYPE (INTWPSEUDO),  ALLOCATABLE, TARGET :: upf(:)

  INTEGER, PARAMETER :: npsx=10

  INTEGER :: &
       nh(npsx),             &! number of beta functions per atomic type
       nhm,                  &! max number of different beta functions per atom
       nbetam                 ! max number of beta functions
  INTEGER :: &
       lmaxkb                ! max angular momentum
  !
  !Former USPP
  INTEGER, PARAMETER :: &
       lmaxx  = 3      ! max non local angular momentum (l=0 to lmaxx)
  !
  INTEGER :: nkb        ! total number of beta functions, with struct.fact.
  !
  INTEGER, ALLOCATABLE ::&
       indv(:,:),        &! indes linking  atomic beta's to beta's in the solid
       nhtol(:,:),       &! correspondence n <-> angular momentum l
       nhtolm(:,:),      &! correspondence n <-> combined lm index for (l,m)
       ijtoh(:,:,:)       ! correspondence beta indexes ih,jh -> composite index ijh
  !
  !
  COMPLEX(DP), ALLOCATABLE, TARGET :: &
       vkb(:,:), vkqb(:,:)             ! all beta functions in reciprocal space
  REAL(DP), ALLOCATABLE :: &
       nhtoj(:,:)              ! correspondence n <-> total angular momentum
  !
  COMPLEX(DP), ALLOCATABLE :: & ! variables for spin-orbit/noncolinear case:
       DKB(:,:,:,:,:)            !
  !
  REAL(DP), ALLOCATABLE :: &
       beta(:,:,:)           ! beta functions for CP (without struct.factor)

contains

  !---------------------------------------------------------------------
  subroutine read_all_pseudo ()
    !
    use intw_utility, only: find_free_unit
    USE intw_reading, only: ntyp
    use intw_input_parameters, only: mesh_dir, prefix
    !
    IMPLICIT NONE
    !
    !     Local variables
    !
    INTEGER :: ios, nr, is, nb, ir, nb1
    INTEGER :: io_unit, ierr
    CHARACTER(256) :: file_pseudo
    CHARACTER(256) :: dum
    CHARACTER(1)   :: tag1
    CHARACTER(2)   :: tag2


    ALLOCATE( UPF(ntyp) )
    !
    ierr = 0
    do is=1,ntyp

      io_unit = find_free_unit()

      if ((is>0).and.(is<9)) then
        write(tag1,"(i1)")is
        write(*,20)"|       - Reading:   "//tag1//"-KBPP.txt"//" ..                  |"
        file_pseudo=trim(trim(adjustl(mesh_dir)))//trim(prefix)//".save.intw/"//tag1//"-KBPP.txt"
      else if ((is>9).and.(is<19) ) then
        write(tag2,"(i2)")is
        file_pseudo=trim(trim(adjustl(mesh_dir)))//trim(prefix)//".save.intw/"//tag2//"-KBPP.txt"
      else
        print*, "ERROR: The num. of species is bigger than 19 (or <0)"
      end if

      OPEN(UNIT=io_unit,FILE=file_pseudo,STATUS='old',FORM='formatted', IOSTAT=ios)

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%psd
      write(*,20)"|                 .. for the specie "//trim(UPF(is)%psd)//"              |"

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%rel

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%has_so

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%nlcc

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%zp

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%ecutwfc

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%lloc

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%lmax

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%nbeta

      allocate( UPF(is)%kbeta( UPF(is)%nbeta ) )
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%kbeta
      UPF(is)%kkbeta = MAXVAL( UPF(is)%kbeta(:) )

      allocate( UPF(is)%lll( UPF(is)%nbeta ) )
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%lll

      allocate( UPF(is)%jjj( UPF(is)%nbeta ) )
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      if (UPF(is)%has_so) then
        read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%jjj
      else
        read(unit=io_unit,fmt=*,iostat=ierr) dum
      end if

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum

      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%mesh

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) UPF(is)%rcloc

      allocate( UPF(is)%dion( UPF(is)%nbeta,UPF(is)%nbeta ) )
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum

      do nb=1,UPF(is)%nbeta
        read(unit=io_unit,fmt=*,iostat=ierr) (UPF(is)%dion(nb,nb1),nb1=1,UPF(is)%nbeta)
      end do

      nr=UPF(is)%mesh
      nb=UPF(is)%nbeta
      allocate( UPF(is)%r(nr), UPF(is)%rab(nr), UPF(is)%vloc(nr) )
      allocate( UPF(is)%beta(nr,nb) )

      read(unit=io_unit,fmt="(a)",iostat=ierr)dum
      do ir=1,nr
        read(unit=io_unit,fmt=*,iostat=ierr) &
        UPF(is)%r(ir), UPF(is)%rab(ir), UPF(is)%vloc(ir), (UPF(is)%beta(ir,nb), nb=1,UPF(is)%nbeta)
      end do

      IF (ierr /= 0) then
         write(unit=*,fmt=*)"ERROR reading PP, ", file_pseudo
         stop
      END IF
      !
      CLOSE(io_unit)

    END DO !is
    !

20 format(A)
30 format(A,F8.2,6X,A)

  END subroutine read_all_pseudo


  subroutine init_KB_projectors(npw, npwx, igk, qpoint_cryst, vkb_)
    !----------------------------------------------------------------------
    !
    !   Calculates beta functions (Kleinman-Bylander projectors), with
    !   structure factor, for all atoms, in reciprocal space
    !
    USE kinds, ONLY: dp
    USE intw_reading, ONLY: nat, ntyp, ityp, tau, bg, tpiba
    USE intw_useful_constants, ONLY : tpi, cmplx_0, cmplx_i
    USE intw_fft, ONLY: gvec_cart
    USE mcf_spline, only: splint_mcf

    implicit none

    !I/O variables

    integer, intent(in) :: npw, npwx !number of PW's
    integer, intent(in) :: igk(npwx) !G list of q vector
    real(dp), intent(in) :: qpoint_cryst(3) !q vector
    complex(dp), intent(out) :: vkb_(npwx,nkb) !beta functions

    !local variables

    integer :: ig, l, lm, na, nt, nb, ih, jkb, iq
    real(DP) :: arg, qpoint_cart(3), gk(3,npw), vkb1(npw,nhm), xdata(nqx), ylm(npw, (lmaxkb+1)**2)
    real(DP), dimension(npw) :: qg(npw), vq(npw)
    complex(DP) :: phase, pref, sk(npw)


    vkb_ = cmplx_0
    !
    qpoint_cart = matmul(bg, qpoint_cryst)
    !
    if (lmaxkb.lt.0) return
    !
    do ig=1,npw
      !
      if (igk(ig)==0) exit
      !
      gk(1,ig) = qpoint_cart(1) + gvec_cart(1, igk(ig))
      gk(2,ig) = qpoint_cart(2) + gvec_cart(2, igk(ig))
      gk(3,ig) = qpoint_cart(3) + gvec_cart(3, igk(ig))
      !
      qg(ig) = gk(1,ig)**2 + gk(2,ig)**2 + gk(3,ig)**2
      !
    enddo
    !
    call intw_real_ylmr2((lmaxkb+1)**2, npw, gk, qg, ylm)
    !
    ! set now qg=|q+G| in atomic units
    !
    do ig=1,npw
      !
      if (igk(ig)==0) exit
      !
      qg(ig) = tpiba * sqrt(qg(ig))
      !
    enddo !ig
    !
    do iq=1,nqx
          !
      xdata(iq) = (iq - 1) * dq
          !
    enddo !ig
    !
    ! |beta_lm(q)> = (4pi/omega).Y_lm(q).f_l(q).(i^l).S(q)
    !
    jkb = 0
    vq = 0.d0
    !
    do nt=1,ntyp
      !
      ! calculate beta in G-space using an interpolation table f_l(q)=\int _0 ^\infty dr r^2 f_l(r) j_l(q.r)
      !
      do nb=1,upf(nt)%nbeta
          !
          do ig=1,npw
            !
            if (igk(ig)==0) exit
            !
            call splint_mcf(xdata, tab(:,nb,nt), tab_d2y(:,nb,nt), nqx, qg(ig), vq(ig))
          enddo !ig
          !
          ! add spherical harmonic part  (Y_lm(q)*f_l(q))
          !
          do ih=1,nh(nt)
            !
            if (nb.eq.indv(ih, nt)) then
                !
                l  = nhtol(ih,nt)
                lm = nhtolm(ih,nt)
                !
                do ig=1,npw
                  !
                  if (igk(ig)==0) exit
                  !
                  vkb1(ig,ih) = ylm(ig,lm) * vq(ig)
                  !
                enddo !ig
            endif !nb
          enddo !ih
      enddo !nt
      !
      ! vkb1 contains all betas including angular part for type nt
      ! now add the structure factor and factor (-i)^l
      !
      do na=1,nat
          !
          ! ordering: first all betas for atoms of type 1
          !           then  all betas for atoms of type 2  and so on
          !
          if (ityp(na).eq.nt) then
            !
            ! qpoint_cart is cart. coordinates
            !
            arg = (  qpoint_cart(1) * tau(1,na) &
                   + qpoint_cart(2) * tau(2,na) &
                   + qpoint_cart(3) * tau(3,na) ) * tpi
            !
            phase = CMPLX(cos(arg), - sin(arg) ,kind=DP)
            !
            do ig=1,npw
                !
                if (igk(ig)==0) exit
                !
                sk (ig) = exp( -tpi*cmplx_i*(  gvec_cart(1,igk(ig))*tau(1,na) &
                                             + gvec_cart(2,igk(ig))*tau(2,na) &
                                             + gvec_cart(3,igk(ig))*tau(3,na) ) )
                !
            enddo !ig
            !
            do ih=1,nh(nt)
                !
                jkb = jkb + 1
                pref = phase * (-cmplx_i)**nhtol(ih,nt)
                !
                do ig=1,npw
                  !
                  if (igk(ig)==0) exit
                  !
                  vkb_(ig, jkb) = vkb1(ig,ih) * sk(ig) * pref
                  !
                enddo
            enddo
          endif
      enddo !n bet
    enddo !ntyp

  end subroutine init_KB_projectors


  subroutine init_pp()
    !----------------------------------------------------------------------
    !
    USE kinds, only: dp
    USE intw_useful_constants, only: fpi, sqrt2, cmplx_0, cmplx_1, cmplx_i
    USE intw_reading, only: ntyp, volume0, lspinorb
    USE intw_spin_orb, only: rot_ylm, fcoef
    USE mcf_spline, only: spline_mcf
    USE intw_utility, only: intgr_spline_gaussq !, simpson
    !
    implicit none
    !
    !     here a few local variables
    !
    integer :: nt, ih, jh, nb, ijv, l, m, ir, iq, is, ndm
    ! various counters
    real(dp), allocatable :: aux(:), aux1(:), besr(:), qtot(:,:)
    ! various work space
    real(dp) :: prefr, pref, qi
    ! the prefactor of the q functions
    ! the prefactor of the beta functions
    ! the modulus of g for each shell
    ! q-point grid for interpolation
    ! the spherical harmonics
    real(dp) ::  vqint, j
    ! interpolated value
    ! J=L+S (noninteger!)
    integer :: n1, m0, m1, n, li, mi, vi, vj, ijs, ispin, jspin, lk, mk, vk, kh
    complex(dp) :: coeff
    real(dp) :: xdata(nqx)
    real(dp) :: d1, ji, jk
    !
    real(dp), external :: spinor
    integer, external :: sph_ind
    !
    !
    !    Initialization of the variables
    !
    ndm = MAXVAL( upf(:)%kkbeta )
    allocate(aux(ndm))
    allocate(aux1(ndm))
    allocate(besr(ndm))
    allocate(qtot(ndm, nbetam*(nbetam+1)/2))
    !
    ! the following prevents an out-of-bound error: upf(nt)%nqlc=2*lmax+1
    ! but in some versions of the PP files lmax is not set to the maximum
    ! l of the beta functions but includes the l of the local potential
    !
    prefr = fpi / volume0

    DKB = cmplx_0

    if (lspinorb) then
      !
      !  In the spin-orbit case we need the unitary matrix u which rotates the
      !  real spherical harmonics and yields the complex ones.
      !
      rot_ylm = cmplx_0
      l = lmaxx
      rot_ylm(l+1, 1) = cmplx_1
      do n1 = 2, 2*l+1, 2
        m = n1/2
        n = l + 1 - m
        rot_ylm(n, n1) = cmplx_1/sqrt2 * (-1)**m
        rot_ylm(n, n1+1) = -cmplx_i/sqrt2 * (-1)**m
        n = l + 1 + m
        rot_ylm(n, n1) = cmplx_1/sqrt2
        rot_ylm(n, n1+1) = cmplx_i/sqrt2
      enddo
      fcoef = cmplx_0
    endif
    !
    !   For each pseudopotential we initialize the indices nhtol, nhtolm,
    !   nhtoj, indv, and if the pseudopotential is of KB type we initialize the
    !   atomic D terms
    !
    do nt = 1, ntyp
      ih = 1
      do nb = 1, upf(nt)%nbeta
        l = upf(nt)%lll(nb)
        do m = 1, 2 * l + 1
          nhtol(ih,nt) = l
          nhtolm(ih,nt) = l*l+m
          indv(ih,nt) = nb
          ih = ih + 1
        enddo
      enddo
      if ( upf(nt)%has_so ) then
        ih = 1
        do nb = 1, upf(nt)%nbeta
          l = upf(nt)%lll(nb)
          j = upf(nt)%jjj(nb)
          do m = 1, 2 * l + 1
            nhtoj(ih, nt) = j
            ih = ih + 1
          enddo
        enddo
      endif
      ! ijtoh map augmentation channel indexes ih and jh to composite
      ! "triangular" index ijh
      ijtoh(:,:,nt) = -1
      ijv = 0
      do ih = 1, nh(nt)
        do jh = ih, nh(nt)
          ijv = ijv + 1
          ijtoh(jh,ih,nt) = ijv
        enddo
      enddo
      !
      !    From now on the only difference between KB and US pseudopotentials
      !    is in the presence of the q and Q functions.
      !
      !    Here we initialize the D of the solid
      !
      if (upf(nt)%has_so .and. lspinorb) then
        !
        !  first calculate the fcoef coefficients
        !
        do ih = 1, nh(nt)
          li = nhtol(ih,nt)
          ji = nhtoj(ih,nt)
          mi = nhtolm(ih,nt) - li*li
          vi = indv(ih,nt)
          do kh=1,nh(nt)
            lk = nhtol(kh,nt)
            jk = nhtoj(kh,nt)
            mk = nhtolm(kh,nt) - lk*lk
            vk = indv(kh,nt)
            if (li == lk .and. abs(ji-jk) < 1.d-7) then
              do ispin=1,2
                do jspin=1,2
                  coeff = cmplx_0
                  do m = -li-1, li
                    m0 = sph_ind(li,ji,m,ispin) + lmaxx + 1
                    m1 = sph_ind(lk,jk,m,jspin) + lmaxx + 1
                    coeff = coeff + rot_ylm(m0,mi)*spinor(li,ji,m,ispin) &
                                    * CONJG(rot_ylm(m1,mk))*spinor(lk,jk,m,jspin)
                  enddo
                  fcoef(ih,kh,ispin,jspin,nt) = coeff
                enddo
              enddo
            endif
          enddo
        enddo
        !
        !   and calculate the bare coefficients
        !
        do ih = 1, nh(nt)
          vi = indv(ih,nt)
          do jh = 1, nh(nt)
            vj = indv(jh,nt)
            ijs = 0
            do ispin=1,2
              do jspin=1,2
                ijs = ijs + 1
                DKB(ih,jh,ispin,jspin,nt) = upf(nt)%dion(vi,vj) * fcoef(ih,jh,ispin,jspin,nt)
                if (vi.ne.vj) fcoef(ih,jh,ispin,jspin,nt) = cmplx_0
              enddo
            enddo
          enddo
        enddo
      else
        do ih = 1, nh(nt)
          do jh = 1, nh(nt)
            if ( nhtol(ih,nt) == nhtol(jh,nt) .and. nhtolm(ih,nt) == nhtolm(jh,nt) ) then
              ir = indv(ih,nt)
              is = indv(jh,nt)
              if (lspinorb) then
                DKB(ih,jh,1,1,nt) = upf(nt)%dion(ir,is)
                DKB(ih,jh,2,2,nt) = upf(nt)%dion(ir,is)
              else
                DKB(ih,jh,1,1,nt) = upf(nt)%dion(ir,is)
              endif
            endif
          enddo
        enddo
      endif
    enddo
    !
    !
    !     fill the interpolation table tab
    !
    pref = fpi / sqrt(volume0)
    tab(:,:,:) = 0.d0
    do nt = 1, ntyp
      do nb = 1, upf(nt)%nbeta
        l = upf(nt)%lll(nb)
        do iq = 1,nqx
          qi = (iq - 1) * dq
          call sph_bes(upf(nt)%kkbeta, upf(nt)%r, qi, l, besr)
          do ir = 1, upf(nt)%kkbeta
            aux(ir) = upf(nt)%beta(ir,nb) * besr(ir) * upf(nt)%r(ir)
          enddo
          !ASIER 29/07/2021
          !call simpson (upf(nt)%kkbeta, aux, upf(nt)%rab, vqint)

          !Integrating by spline + gauss 2. order
          vqint = intgr_spline_gaussq( upf(nt)%r(1:upf(nt)%kkbeta), aux )

          tab(iq,nb,nt) = vqint * pref

        enddo
      enddo
    enddo

    ! initialize spline interpolation
    do iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    enddo
    do nt = 1, ntyp
      do nb = 1, upf(nt)%nbeta
        d1 = (tab(2,nb,nt) - tab(1,nb,nt)) / dq
        call spline_mcf(xdata, tab(:,nb,nt), nqx, tab_d2y(:,nb,nt))
      enddo
    enddo

    deallocate(qtot)
    deallocate(besr)
    deallocate(aux1)
    deallocate(aux)

  end subroutine init_pp


  SUBROUTINE phq_init(q_cart)
    !----------------------------------------------------------------------------
    !
    use kinds, only: dp
    use intw_ph, only: eigqts
    use intw_reading, only: nat, tau, ntyp, tpiba2, ngm, volume0
    use intw_useful_constants, only: tpi
    use intw_fft, only: gvec_cart
    !
    implicit none
    !
    real(dp), intent(in) :: q_cart(3)
    !
    ! local variables
    real(dp) :: arg ! the argument of the phase
    integer :: nt, na


    DO na = 1, nat
      !
      arg = (  q_cart(1) * tau(1,na) &
             + q_cart(2) * tau(2,na) &
             + q_cart(3) * tau(3,na) ) * tpi
      !
      eigqts(na) = CMPLX( COS(arg), -SIN(arg), kind=DP)
      !
    END DO
    !
    vlocq(:,:) = 0.D0
    !
    DO nt = 1, ntyp
      CALL setlocq( q_cart, upf(nt)%mesh, upf(nt)%mesh, upf(nt)%rab, upf(nt)%r,&
                    upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngm, gvec_cart, volume0, &
                    vlocq(1,nt) )
    END DO

  END SUBROUTINE phq_init

end module intw_pseudo
