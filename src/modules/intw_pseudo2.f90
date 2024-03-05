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
  public :: read_all_pseudo, init_KB_projectors, init_pp, phq_init, setlocq, &
            setlocq_coul, allocate_nlpot, set_nqxq, allocate_phq, deallocate_phq
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
    COMPLEX (DP), ALLOCATABLE :: fcoef(:,:,:,:,:) ! function needed to account for spinors.
    COMPLEX (DP) :: rot_ylm(2*lmaxx+1,2*lmaxx+1)  ! transform real spherical harmonics into complex ones
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
    if (lspinorb) allocate(fcoef(nhm,nhm,2,2,ntyp))

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
    use intw_reading, only: nat, tau, ntyp, tpiba2, ngm, volume0
    use intw_useful_constants, only: tpi
    use intw_fft, only: gvec_cart
    !
    implicit none
    !
    real(dp), intent(in) :: q_cart(3)
    !
    ! local variables
    integer :: nt, na


    !
    vlocq(:,:) = 0.D0
    !
    DO nt = 1, ntyp
      CALL setlocq( q_cart, upf(nt)%mesh, upf(nt)%mesh, upf(nt)%rab, upf(nt)%r,&
                    upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngm, gvec_cart, volume0, &
                    vlocq(1,nt) )
    END DO

  END SUBROUTINE phq_init


  subroutine setlocq(q_cart, mesh, msh, rab, r, vloc_at, zp, tpiba2, ngm, g, omega, vloc)
    !----------------------------------------------------------------------
    !
    !    This routine computes the Fourier transform of the local
    !    part of the pseudopotential in the q+G vectors.
    !
    !    The local pseudopotential of the US case is always in
    !    numerical form, expressed in Ry units.
    !
    USE kinds, only : DP
    USE intw_useful_constants, ONLY : e2, fpi, pi
    USE intw_utility, ONLY : simpson, qe_erf, qe_erfc
    !
    implicit none
    !
    !    first the dummy variables
    !
    integer :: ngm, mesh, msh
    ! input: the number of G vectors
    ! input: the dimensions of the mesh
    ! input: mesh points for radial integration

    real(DP), intent(in) :: q_cart(3), zp, rab(mesh), r(mesh), vloc_at(mesh), tpiba2, omega, g(3,ngm)
    ! input: the q point
    ! input: valence pseudocharge
    ! input: the derivative of mesh points
    ! input: the mesh points
    ! input: the pseudo on the radial
    ! input: 2 pi / alat
    ! input: the volume of the unit cell
    ! input: the g vectors coordinates
    real(DP), intent(out) :: vloc(ngm)
    ! output: the fourier transform of the potential
    !
    !    and the local variables
    !
    real(DP), parameter :: eps = 1.d-8
    real(DP) :: vlcp, vloc0, fac, g2a, aux(mesh), aux1(mesh), gx
    ! auxiliary variables
    ! gx = modulus of g vectors
    !real(DP), external :: qe_erf
    ! the erf function
    integer :: ig, ir
    ! counters
    !
    ! Pseudopotentials in numerical form (Vnl(lloc) contain the local part)
    ! in order to perform the Fourier transform, a term erf(r)/r is
    ! subtracted in real space and added again in G space
    !
    ! first the G=0 term
    !
    do ir = 1, msh
      aux(ir) = r(ir) * (r(ir) * vloc_at(ir) + zp * e2)
    enddo
    call simpson(msh, aux, rab, vloc0)
    !
    !   here the G<>0 terms, we first compute the part of the integrand func
    !   indipendent of |G| in real space
    !
    do ir = 1, msh
      aux1(ir) = r(ir) * vloc_at(ir) + zp * e2 * qe_erf(r(ir))
    enddo
    fac = zp * e2 / tpiba2
    !
    !    and here we perform the integral, after multiplying for the |G|
    !    dependent  part
    !
    do ig = 1, ngm
      g2a = (q_cart(1) + g(1,ig))**2 + (q_cart(2) + g(2,ig) )**2 + (q_cart(3) + g(3,ig))**2
      if (g2a < eps) then
        vloc (ig) = vloc0
      else
        gx = sqrt (g2a * tpiba2)
        do ir = 1, msh
          aux(ir) = aux1(ir) * sin(gx * r(ir)) / gx
        enddo
        call simpson(msh, aux, rab, vlcp)
        !
        !     here we add the analytic fourier transform of the erf function
        !
        vlcp = vlcp - fac * exp( - g2a * tpiba2 * 0.25d0) / g2a
        vloc(ig) = vlcp
      endif
    enddo

    vloc(:) = vloc(:) * fpi / omega

  end subroutine setlocq


  subroutine setlocq_coul(q_cart, zp, tpiba2, ngm, g, omega, vloc)
    !----------------------------------------------------------------------
    !
    !    Fourier transform of the Coulomb potential - For all-electron
    !    calculations, in specific cases only, for testing purposes
    !
    USE kinds, ONLY: DP
    USE intw_useful_constants, ONLY : fpi, e2, eps_8
    implicit none
    !
    integer, intent(in) :: ngm
    real(DP), intent(in) :: q_cart(3), zp, tpiba2, omega, g(3,ngm)
    real(DP), intent (out) :: vloc(ngm)
    !
    real(DP) :: g2a
    integer :: ig

    do ig = 1, ngm
      g2a = (q_cart(1) + g(1,ig))**2 + (q_cart(2) + g(2,ig))**2 + (q_cart(3) + g(3,ig))**2
      if (g2a < eps_8) then
          vloc(ig) = 0.d0
      else
          vloc(ig) = - fpi * zp *e2 / omega / tpiba2 / g2a
      endif
    enddo

  end subroutine setlocq_coul


  subroutine allocate_nlpot()
    !-----------------------------------------------------------------------
    ! Adapted from QEspresso4.3.2 Asier&&Idoia
    !
    !     ngk           !  number of plane waves (for each k point)
    !     npwx          !  maximum number of plane waves
    !     nqx           !  number of points of the interpolation table
    !     nqxq          !  as above, for q-function interpolation table
    !
    !
    use kinds, only: dp
    use intw_reading, only: nat, ntyp, ityp, ecutwfc, nspin, lspinorb, ng_max

    implicit none

    !local variables

    integer :: nt, nb, na


    lmaxkb = - 1
    DO nt = 1, ntyp
      !
      nh(nt) = 0
      !
      DO nb = 1, upf(nt)%nbeta
        nh(nt) = nh(nt) + 2 * upf(nt)%lll(nb) + 1
        lmaxkb = MAX(lmaxkb, upf(nt)%lll(nb))
      ENDDO
      !
    ENDDO
    !
    ! calculate the maximum number of beta functions
    !
    nhm = MAXVAL(nh(1:ntyp))
    nbetam = MAXVAL(upf(:)%nbeta)
    !
    ! Number of beta functions
    !
    nkb = 0
    do na = 1, nat
      nt = ityp(na)
      nkb = nkb + nh(nt)
    enddo

    allocate(indv( nhm,ntyp))
    allocate(nhtol(nhm,ntyp))
    allocate(nhtolm(nhm,ntyp))
    allocate(nhtoj(nhm,ntyp))
    allocate(ijtoh(nhm,nhm,ntyp))
    if (lspinorb) then
      allocate (DKB(nhm,nhm,nspin,nspin,ntyp))
    else
      allocate (DKB(nhm,nhm,1,1,ntyp))
    endif
    !
    !
    ! Calculate dimensions for array tab (including a possible factor
    ! coming from cell contraction during variable cell relaxation/MD)
    !
    nqx = INT( (sqrt(2*ecutwfc) / dq + 4) ) ! x2 zeren Ry -> Hartree egin behar

    allocate(tab(2*nqx,nbetam,ntyp)) ! ASIER originala nqx, errepasatu

    ! d2y is for the cubic splines
    allocate(tab_d2y(nqx,nbetam,ntyp))

    if (allocated(vkb)) deallocate(vkb)
    if (allocated(vkqb)) deallocate(vkqb)

    if (nkb > 0) then
      allocate(vkb(nG_max,nkb))
      allocate(vkqb(nG_max,nkb))
    end if
    !
    return

  end subroutine allocate_nlpot


  subroutine set_nqxq(qpoint)
    !
    USE kinds, only: dp
    USE intw_reading, ONLY: ecutrho, tpiba2, bg

    implicit none

    !I/O variables

    real(dp), intent(in) :: qpoint(3)
    real(dp) :: xq(3), qnorm
    !
    xq = matmul(bg, qpoint)
    !
    qnorm =sqrt(xq(1)**2+xq(2)**2+xq(3)**2)
    !
    ! This routine is called also by the phonon code, in which case it should
    ! allocate an array that includes q+G vectors up to |q+G|_max <= |Gmax|+|q|
    !
    nqxq = INT( (( (sqrt(2*ecutrho/tpiba2) + qnorm ) / dq + 4.0_dp) )) ! x2 zeren Ry -> Hartree egin behar
    !
    !
    return

  end subroutine set_nqxq


  subroutine allocate_phq
    !-----------------------------------------------------------------------
    !JLB: Commented out all bands references, not used.
    !     Many variables are declared and not used, I don't understand why.
    !     Discuss with Haritz/Asier, and then remove and clean up.
    !
    USE intw_reading, only : nat, ntyp, nr1, nr2, nr3, ngm

    implicit none

    allocate(vlocq (ngm, ntyp))

  end subroutine allocate_phq


  subroutine deallocate_phq
    !-----------------------------------------------------------------------
    !
    implicit none
    !
    deallocate(vlocq)

    return
  end subroutine deallocate_phq

end module intw_pseudo
