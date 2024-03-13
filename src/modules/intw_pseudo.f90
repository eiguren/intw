module intw_pseudo
  !
  use kinds, only: dp
  !
  implicit none
  !
  save
  !
  ! variables
  public :: intwpseudo, upf, &
            nt_max, l_max, &
            dq, nqx, tab, tab_d2y, &
            nh, nhm, nbetam, lmaxkb, &
            nkb, indv, nhtol, nhtolm, nhtoj, &
            DKB, vlocq, vkb, vkqb
  !
  ! subroutines
  public :: read_all_pseudo, init_KB_projectors, init_pp, phq_init, setlocq, &
            setlocq_coul, allocate_nlpot, allocate_phq, deallocate_phq
  !
  private
  !
  type intwpseudo
    character(len=2)                           :: psd = ' ' ! Element label
    character(len=6)                           :: rel = ' ' ! relativistic: {no|scalar|full}
    logical                                    :: nlcc      ! Non linear core corrections
    real(kind=dp)                              :: zp        ! z valence
    real(kind=dp)                              :: ecutwfc   ! Suggested cut-off for wfc
    logical                                    :: has_so    ! If .true. includes spin-orbit
    !
    integer                                    :: mesh      ! Number of points in the radial mesh
    real(kind=dp), dimension(:),   allocatable :: r         ! r(mesh)  radial grid
    real(kind=dp), dimension(:),   allocatable :: rab       ! rab(mesh) dr(i)/di
    !
    real(kind=dp)                              :: rcloc     ! vloc = v_ae for r > rcloc
    integer                                    :: lloc      ! l of channel used to generate local potential
    real(kind=dp), dimension(:),   allocatable :: vloc      ! vloc(mesh) local atomic potential
    !
    integer                                    :: nbeta     ! Number of beta projectors
    integer                                    :: lmax      ! Max l component in beta
    integer,       dimension(:),   allocatable :: kbeta     ! kbeta(nbeta): number of grid points used for betas
                                                            ! this defines the cutoff radius for each of them.
    integer                                    :: kkbeta    ! Max number of grid points used for betas
    integer,       dimension(:),   allocatable :: lll       ! lll(nbeta) l of each projector
    real(kind=dp), dimension(:),   allocatable :: jjj       ! jjj(nbeta) j=l+1/2 or l-1/2 of beta
    real(kind=dp), dimension(:,:), allocatable :: beta      ! beta(mesh,nbeta) projectors
    real(kind=dp), dimension(:,:), allocatable :: dion      ! dion(nbeta,nbeta) atomic D_{mu,nu}
  end type intwpseudo
  !
  !
  type(intwpseudo), dimension(:), allocatable :: upf
  !
  !
  integer, parameter :: nt_max = 10 ! Max number of different atomic types
  integer, parameter :: l_max = 3   ! Max non local angular momentum (l=0 to l_max)
  !
  !
  real(kind=dp), parameter   :: dq = 0.01d0    ! Space between points in the pseudopotential tab
  integer                    :: nqx            ! Number of interpolation points
  real(kind=dp), allocatable :: tab(:,:,:)     ! Interpolation table for PPs
  real(kind=dp), allocatable :: tab_d2y(:,:,:) ! For cubic splines
  !
  !
  integer :: nh(nt_max)  ! Number of beta(lm) functions per atomic type
  integer :: nhm         ! Max number of beta(lm) functions per atomic type
  integer :: nbetam      ! Max number of beta functions per atomic type
  integer :: lmaxkb      ! Max angular momentum of beta functions
  !
  !
  integer :: nkb ! Total number of beta(lm) functions in the solid
  !
  !
  integer,       allocatable :: indv(:,:)   ! Link between index of beta(lm) function in the solid -> index of beta function in the atomic type
  integer,       allocatable :: nhtol(:,:)  ! Link between index of beta(lm) function in the atomic type -> angular momentum l
  integer,       allocatable :: nhtolm(:,:) ! Link between index of beta(lm) function in the atomic type -> combined lm angular momentum index l*l+m
  real(kind=dp), allocatable :: nhtoj(:,:)  ! Link between index of beta(lm) function in the atomic type -> total angular momentum j
  !
  !
  complex(kind=dp), allocatable :: DKB(:,:,:,:,:) ! D_{mu,nu} matrix for beta(lm) functions for each atomic type
  !
  !
  real(kind=dp), allocatable :: vlocq(:,:) ! Local potential in reciprocal space
  complex(kind=dp), allocatable, target :: vkb(:,:), vkqb(:,:) ! All beta functions in reciprocal space

contains

  !---------------------------------------------------------------------
  subroutine read_all_pseudo()
    !
    use intw_utility, only: find_free_unit
    use intw_reading, only: ntyp
    use intw_input_parameters, only: mesh_dir, prefix
    !
    implicit none
    !
    !     Local variables
    !
    integer :: ios, nr, nt, nb, ir, nb1
    integer :: io_unit, ierr
    character(256) :: file_pseudo
    character(256) :: dum
    character(1)   :: tag1
    character(2)   :: tag2


    allocate(upf(ntyp))
    !
    ierr = 0
    do nt = 1, ntyp

      io_unit = find_free_unit()

      if ( nt>0 .and. nt<9 ) then
        write(tag1,"(i1)")nt
        write(*,20)"|       - Reading:   "//tag1//"-KBPP.txt"//" ..                  |"
        file_pseudo=trim(trim(adjustl(mesh_dir)))//trim(prefix)//".save.intw/"//tag1//"-KBPP.txt"
      else if ( nt>9 .and. nt<19 ) then
        write(tag2,"(i2)")nt
        file_pseudo=trim(trim(adjustl(mesh_dir)))//trim(prefix)//".save.intw/"//tag2//"-KBPP.txt"
      else
        print*, "ERROR: The num. of species is bigger than 19 (or <0)"
      end if

      open(unit=io_unit ,file=file_pseudo, status='old', form='formatted', iostat=ios)

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%psd
      write(*,20)"|                 .. for the specie "//trim(upf(nt)%psd)//"              |"

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%rel

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%has_so

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%nlcc

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%zp

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%ecutwfc

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%lloc

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%lmax

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%nbeta

      allocate(upf(nt)%kbeta(upf(nt)%nbeta))
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%kbeta
      upf(nt)%kkbeta = maxval( upf(nt)%kbeta(:) )

      allocate(upf(nt)%lll(upf(nt)%nbeta))
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%lll

      allocate(upf(nt)%jjj(upf(nt)%nbeta))
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      if (upf(nt)%has_so) then
        read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%jjj
      else
        read(unit=io_unit,fmt=*,iostat=ierr) dum
      end if

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum

      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%mesh

      read(unit=io_unit,fmt="(a)",iostat=ierr) dum
      read(unit=io_unit,fmt=*,iostat=ierr) upf(nt)%rcloc

      allocate(upf(nt)%dion( upf(nt)%nbeta,upf(nt)%nbeta))
      read(unit=io_unit,fmt="(a)",iostat=ierr) dum

      do nb = 1, upf(nt)%nbeta
        read(unit=io_unit,fmt=*,iostat=ierr) (upf(nt)%dion(nb,nb1),nb1=1,upf(nt)%nbeta)
      end do

      nr=upf(nt)%mesh
      nb=upf(nt)%nbeta
      allocate(upf(nt)%r(nr), upf(nt)%rab(nr), upf(nt)%vloc(nr))
      allocate(upf(nt)%beta(nr,nb))

      read(unit=io_unit,fmt="(a)",iostat=ierr)dum
      do ir = 1, nr
        read(unit=io_unit,fmt=*,iostat=ierr) &
        upf(nt)%r(ir), upf(nt)%rab(ir), upf(nt)%vloc(ir), (upf(nt)%beta(ir,nb), nb=1,upf(nt)%nbeta)
      end do

      if (ierr /= 0) then
         write(unit=*,fmt=*)"ERROR reading PP, ", file_pseudo
         stop
      end if
      !
      close(io_unit)

    end do !nt
    !

20 format(A)
30 format(A,F8.2,6X,A)

  end subroutine read_all_pseudo


  subroutine init_KB_projectors(npw, npwx, igk, qpoint_cryst, vkb_)
    !----------------------------------------------------------------------
    !
    !   Calculates beta functions (Kleinman-Bylander projectors), with
    !   structure factor, for all atoms, in reciprocal space
    !
    use kinds, only: dp
    use intw_reading, only: nat, ntyp, ityp, tau, bg, tpiba
    use intw_useful_constants, only: tpi, cmplx_0, cmplx_i
    use intw_fft, only: gvec_cart
    use mcf_spline, only: splint_mcf

    implicit none

    !I/O variables

    integer, intent(in) :: npw, npwx !number of PW's
    integer, intent(in) :: igk(npwx) !G list of q vector
    real(kind=dp), intent(in) :: qpoint_cryst(3) !q vector
    complex(kind=dp), intent(out) :: vkb_(npwx,nkb) !beta functions

    !local variables

    integer :: ig, l, lm, na, nt, nb, ih, jkb, iq
    real(kind=dp) :: arg, qpoint_cart(3), gk(3,npw), vkb1(npw,nhm), xdata(nqx), ylm(npw, (lmaxkb+1)**2)
    real(kind=dp), dimension(npw) :: qg(npw), vq(npw)
    complex(kind=dp) :: phase, pref, sk(npw)


    vkb_ = cmplx_0
    !
    qpoint_cart = matmul(bg, qpoint_cryst)
    !
    if (lmaxkb.lt.0) return
    !
    do ig = 1, npw
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
    do ig = 1, npw
      !
      if (igk(ig)==0) exit
      !
      qg(ig) = tpiba * sqrt(qg(ig))
      !
    enddo !ig
    !
    do iq = 1, nqx
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
    do nt = 1, ntyp
      !
      ! calculate beta in G-space using an interpolation table f_l(q)=\int _0 ^\infty dr r^2 f_l(r) j_l(q.r)
      !
      do nb = 1, upf(nt)%nbeta
          !
          do ig = 1, npw
            !
            if (igk(ig)==0) exit
            !
            call splint_mcf(xdata, tab(:,nb,nt), tab_d2y(:,nb,nt), nqx, qg(ig), vq(ig))
          enddo !ig
          !
          ! add spherical harmonic part  (Y_lm(q)*f_l(q))
          !
          do ih = 1, nh(nt)
            !
            if (nb.eq.indv(ih, nt)) then
                !
                l  = nhtol(ih,nt)
                lm = nhtolm(ih,nt)
                !
                do ig = 1, npw
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
      do na = 1, nat
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
            phase = cmplx(cos(arg), - sin(arg) ,kind=dp)
            !
            do ig = 1, npw
                !
                if (igk(ig)==0) exit
                !
                sk (ig) = exp( -tpi*cmplx_i*(  gvec_cart(1,igk(ig))*tau(1,na) &
                                             + gvec_cart(2,igk(ig))*tau(2,na) &
                                             + gvec_cart(3,igk(ig))*tau(3,na) ) )
                !
            enddo !ig
            !
            do ih = 1, nh(nt)
                !
                jkb = jkb + 1
                pref = phase * (-cmplx_i)**nhtol(ih,nt)
                !
                do ig = 1, npw
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
    use kinds, only: dp
    use intw_useful_constants, only: fpi, sqrt2, cmplx_0, cmplx_1, cmplx_i
    use intw_reading, only: ntyp, volume0, lspinorb
    use mcf_spline, only: spline_mcf
    use intw_utility, only: intgr_spline_gaussq !, simpson
    !
    implicit none
    !
    !     here a few local variables
    !
    integer :: nt, ih, jh, nb, l, m, ir, iq, is, ndm
    ! various counters
    real(kind=dp), allocatable :: aux(:), aux1(:), besr(:), qtot(:,:)
    ! various work space
    real(kind=dp) :: prefr, pref, qi
    ! the prefactor of the q functions
    ! the prefactor of the beta functions
    ! the modulus of g for each shell
    ! q-point grid for interpolation
    ! the spherical harmonics
    real(kind=dp) ::  vqint, j
    ! interpolated value
    ! J=L+S (noninteger!)
    integer :: n1, m0, m1, n, li, mi, vi, vj, ijs, ispin, jspin, lk, mk, vk, kh
    complex(kind=dp) :: coeff
    real(kind=dp) :: xdata(nqx)
    real(kind=dp) :: d1, ji, jk
    !
    complex(kind=dp), allocatable :: fcoef(:,:,:,:,:) ! function needed to account for spinors.
    complex(kind=dp) :: rot_ylm(2*l_max+1,2*l_max+1)  ! transform real spherical harmonics into complex ones
    !
    real(kind=dp), external :: spinor
    integer, external :: sph_ind
    !
    !
    !    Initialization of the variables
    !
    ndm = maxval(upf(:)%kkbeta)
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
      l = l_max
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
          do kh = 1, nh(nt)
            lk = nhtol(kh,nt)
            jk = nhtoj(kh,nt)
            mk = nhtolm(kh,nt) - lk*lk
            vk = indv(kh,nt)
            if (li == lk .and. abs(ji-jk) < 1.d-7) then
              do ispin = 1, 2
                do jspin = 1, 2
                  coeff = cmplx_0
                  do m = -li-1, li
                    m0 = sph_ind(li,ji,m,ispin) + l_max + 1
                    m1 = sph_ind(lk,jk,m,jspin) + l_max + 1
                    coeff = coeff +         rot_ylm(m0,mi)  * spinor(li,ji,m,ispin) &
                                    * conjg(rot_ylm(m1,mk)) * spinor(lk,jk,m,jspin)
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
            do ispin = 1, 2
              do jspin = 1, 2
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
          !call simpson(upf(nt)%kkbeta, aux, upf(nt)%rab, vqint)

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


  subroutine phq_init(q_cart)
    !----------------------------------------------------------------------------
    !
    use kinds, only: dp
    use intw_reading, only: ntyp, tpiba2, ngm, volume0
    use intw_useful_constants, only: tpi
    use intw_fft, only: gvec_cart
    !
    implicit none
    !
    real(kind=dp), intent(in) :: q_cart(3)
    !
    ! local variables
    integer :: nt


    !
    vlocq(:,:) = 0.d0
    !
    do nt = 1, ntyp
      call setlocq( q_cart, upf(nt)%mesh, upf(nt)%mesh, upf(nt)%rab, upf(nt)%r,&
                    upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngm, gvec_cart, volume0, &
                    vlocq(1,nt) )
    end do

  end subroutine phq_init


  subroutine setlocq(q_cart, mesh, msh, rab, r, vloc_at, zp, tpiba2, ngm, g, omega, vloc)
    !----------------------------------------------------------------------
    !
    !    This routine computes the Fourier transform of the local
    !    part of the pseudopotential in the q+G vectors.
    !
    !    The local pseudopotential of the US case is always in
    !    numerical form, expressed in Ry units.
    !
    use kinds, only: dp
    use intw_useful_constants, only: e2, fpi, pi
    use intw_utility, only: simpson, qe_erf, qe_erfc
    !
    implicit none
    !
    !    first the dummy variables
    !
    integer :: ngm, mesh, msh
    ! input: the number of G vectors
    ! input: the dimensions of the mesh
    ! input: mesh points for radial integration

    real(kind=dp), intent(in) :: q_cart(3), zp, rab(mesh), r(mesh), vloc_at(mesh), tpiba2, omega, g(3,ngm)
    ! input: the q point
    ! input: valence pseudocharge
    ! input: the derivative of mesh points
    ! input: the mesh points
    ! input: the pseudo on the radial
    ! input: 2 pi / alat
    ! input: the volume of the unit cell
    ! input: the g vectors coordinates
    real(kind=dp), intent(out) :: vloc(ngm)
    ! output: the fourier transform of the potential
    !
    !    and the local variables
    !
    real(kind=dp), parameter :: eps = 1.d-8
    real(kind=dp) :: vlcp, vloc0, fac, g2a, aux(mesh), aux1(mesh), gx
    ! auxiliary variables
    ! gx = modulus of g vectors
    !real(kind=dp), external :: qe_erf
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
      g2a = (q_cart(1) + g(1,ig))**2 + (q_cart(2) + g(2,ig))**2 + (q_cart(3) + g(3,ig))**2
      if (g2a < eps) then
        vloc(ig) = vloc0
      else
        gx = sqrt(g2a * tpiba2)
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
    use kinds, only: dp
    use intw_useful_constants, only: fpi, e2, eps_8
    implicit none
    !
    integer, intent(in) :: ngm
    real(kind=dp), intent(in) :: q_cart(3), zp, tpiba2, omega, g(3,ngm)
    real(kind=dp), intent(out) :: vloc(ngm)
    !
    real(kind=dp) :: g2a
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
    !
    !
    use kinds, only: dp
    use intw_reading, only: nat, ntyp, ityp, ecutwfc, nspin, lspinorb, ng_max

    implicit none

    !local variables

    integer :: nt, nb, na


    lmaxkb = - 1
    do nt = 1, ntyp
      !
      nh(nt) = 0
      !
      do nb = 1, upf(nt)%nbeta
        nh(nt) = nh(nt) + 2 * upf(nt)%lll(nb) + 1
        lmaxkb = max(lmaxkb, upf(nt)%lll(nb))
      end do
      !
    end do
    !
    ! calculate the maximum number of beta functions
    !
    nhm = maxval(nh(1:ntyp))
    nbetam = maxval(upf(:)%nbeta)
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
    if (lspinorb) then
      allocate(DKB(nhm,nhm,nspin,nspin,ntyp))
    else
      allocate(DKB(nhm,nhm,1,1,ntyp))
    endif
    !
    !
    ! Calculate dimensions for array tab (including a possible factor
    ! coming from cell contraction during variable cell relaxation/MD)
    !
    nqx = int( sqrt(2*ecutwfc) / dq + 4 ) ! x2 zeren Ry -> Hartree egin behar

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


  subroutine allocate_phq
    !-----------------------------------------------------------------------
    !
    use intw_reading, only: ngm, ntyp

    implicit none

    allocate(vlocq(ngm, ntyp))

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
