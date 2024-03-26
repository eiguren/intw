module intw_pseudo_non_local

  use kinds, only: dp

  implicit none

  public :: l_kb_max, nh, nhm, nbetam, lmaxkb, &
            nkb, indv, nhtol, nhtolm, nhtoj, &
            DKB, vkb, vkqb

  public :: init_KB_projectors, init_pp, allocate_nlpot, &
            multiply_psi_by_vKB, multiply_psi_by_dvKB

  private


  integer, parameter :: l_kb_max = 3 ! Max non local angular momentum (l=0 to l_kb_max)

  real(kind=dp), parameter   :: dq = 0.01d0    ! Space between points in the pseudopotential tab
  integer                    :: nqx            ! Number of interpolation points
  real(kind=dp), allocatable :: tab(:,:,:)     ! Interpolation table for PPs
  real(kind=dp), allocatable :: tab_d2y(:,:,:) ! For cubic splines

  integer, allocatable :: nh(:)  ! Number of beta(lm) functions per atomic type
  integer              :: nhm         ! Max number of beta(lm) functions per atomic type
  integer              :: nbetam      ! Max number of beta functions per atomic type
  integer              :: lmaxkb      ! Max angular momentum of beta functions

  integer :: nkb ! Total number of beta(lm) functions in the solid
  integer,       allocatable :: indv(:,:)   ! Link between index of beta(lm) function in the solid -> index of beta function in the atomic type
  integer,       allocatable :: nhtol(:,:)  ! Link between index of beta(lm) function in the atomic type -> angular momentum l
  integer,       allocatable :: nhtolm(:,:) ! Link between index of beta(lm) function in the atomic type -> combined lm angular momentum index l*l+m
  real(kind=dp), allocatable :: nhtoj(:,:)  ! Link between index of beta(lm) function in the atomic type -> total angular momentum j

  complex(kind=dp), allocatable :: DKB(:,:,:,:,:) ! D_{mu,nu} matrix for beta(lm) functions for each atomic type

  complex(kind=dp), allocatable, target :: vkb(:,:), vkqb(:,:) ! All beta functions in reciprocal space


contains

  subroutine init_KB_projectors(npw, npwx, igk, qpoint_cryst, vkb_)
    !----------------------------------------------------------------------
    !
    !   Calculates beta functions (Kleinman-Bylander projectors), with
    !   structure factor, for all atoms, in reciprocal space
    !
    use intw_reading, only: nat, ntyp, ityp, tau, bg, tpiba
    use intw_useful_constants, only: tpi, cmplx_0, cmplx_i
    use intw_fft, only: gvec_cart
    use mcf_spline, only: splint_mcf
    use intw_pseudo, only: upf
    use intw_utility, only: real_ylmr2

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
    end do
    !
    call real_ylmr2(lmaxkb, npw, gk, qg, ylm)
    !
    ! set now qg=|q+G| in atomic units
    !
    do ig = 1, npw
      !
      if (igk(ig)==0) exit
      !
      qg(ig) = tpiba * sqrt(qg(ig))
      !
    end do !ig
    !
    do iq = 1, nqx
          !
      xdata(iq) = (iq - 1) * dq
          !
    end do !ig
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
          end do !ig
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
                end do !ig
            end if !nb
          end do !ih
      end do !nt
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
            end do !ig
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
                end do
            end do
          end if
      end do !n bet
    end do !ntyp

  end subroutine init_KB_projectors


  subroutine init_pp()
    !----------------------------------------------------------------------
    !
    use intw_useful_constants, only: fpi, sqrt2, cmplx_0, cmplx_1, cmplx_i
    use intw_reading, only: ntyp, volume0, lspinorb
    use mcf_spline, only: spline_mcf
    use intw_utility, only: sphb, intgr_spline_gaussq !, simpson
    use intw_pseudo, only: upf
    !
    implicit none
    !
    !     here a few local variables
    !
    integer :: nt, ih, jh, nb, l, m, ir, iq, is, ndm
    ! various counters
    real(kind=dp), allocatable :: aux(:), aux1(:), qtot(:,:)
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
    complex(kind=dp) :: rot_ylm(2*l_kb_max+1,2*l_kb_max+1)  ! transform real spherical harmonics into complex ones
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
      l = l_kb_max
      rot_ylm(l+1, 1) = cmplx_1
      do n1 = 2, 2*l+1, 2
        m = n1/2
        n = l + 1 - m
        rot_ylm(n, n1) = cmplx_1/sqrt2 * (-1)**m
        rot_ylm(n, n1+1) = -cmplx_i/sqrt2 * (-1)**m
        n = l + 1 + m
        rot_ylm(n, n1) = cmplx_1/sqrt2
        rot_ylm(n, n1+1) = cmplx_i/sqrt2
      end do
      fcoef = cmplx_0
    end if
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
        end do
      end do
      if ( upf(nt)%has_so ) then
        ih = 1
        do nb = 1, upf(nt)%nbeta
          l = upf(nt)%lll(nb)
          j = upf(nt)%jjj(nb)
          do m = 1, 2 * l + 1
            nhtoj(ih, nt) = j
            ih = ih + 1
          end do
        end do
      end if
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
                    m0 = sph_ind(li,ji,m,ispin) + l_kb_max + 1
                    m1 = sph_ind(lk,jk,m,jspin) + l_kb_max + 1
                    coeff = coeff +         rot_ylm(m0,mi)  * spinor(li,ji,m,ispin) &
                                    * conjg(rot_ylm(m1,mk)) * spinor(lk,jk,m,jspin)
                  end do
                  fcoef(ih,kh,ispin,jspin,nt) = coeff
                end do
              end do
            end if
          end do
        end do
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
              end do
            end do
          end do
        end do
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
              end if
            end if
          end do
        end do
      end if
    end do
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
          aux = upf(nt)%beta(:,nb) * sphb(l, qi*upf(nt)%r) * upf(nt)%r
          !ASIER 29/07/2021
          !call simpson(upf(nt)%kkbeta, aux, upf(nt)%rab, vqint)

          !Integrating by spline + gauss 2. order
          vqint = intgr_spline_gaussq( upf(nt)%r(1:upf(nt)%kkbeta), aux )

          tab(iq,nb,nt) = vqint * pref

        end do
      end do
    end do

    ! initialize spline interpolation
    do iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    end do
    do nt = 1, ntyp
      do nb = 1, upf(nt)%nbeta
        d1 = (tab(2,nb,nt) - tab(1,nb,nt)) / dq
        call spline_mcf(xdata, tab(:,nb,nt), nqx, tab_d2y(:,nb,nt))
      end do
    end do

    deallocate(qtot)
    deallocate(aux1)
    deallocate(aux)

  end subroutine init_pp


  subroutine allocate_nlpot()
    !-----------------------------------------------------------------------
    ! Adapted from QEspresso4.3.2 Asier&&Idoia
    !
    !     ngk           !  number of plane waves (for each k point)
    !     npwx          !  maximum number of plane waves
    !     nqx           !  number of points of the interpolation table
    !
    !
    use intw_reading, only: nat, ntyp, ityp, ecutwfc, nspin, lspinorb, ng_max
    use intw_pseudo, only: upf

    implicit none

    !local variables

    integer :: nt, nb, na

    allocate(nh(ntyp))

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
    nhm = maxval(nh)
    nbetam = maxval(upf(:)%nbeta)
    !
    ! Number of beta functions
    !
    nkb = 0
    do na = 1, nat
      nt = ityp(na)
      nkb = nkb + nh(nt)
    end do

    allocate(indv(nhm,ntyp))
    allocate(nhtol(nhm,ntyp))
    allocate(nhtolm(nhm,ntyp))
    allocate(nhtoj(nhm,ntyp))
    if (lspinorb) then
      allocate(DKB(nhm,nhm,nspin,nspin,ntyp))
    else
      allocate(DKB(nhm,nhm,1,1,ntyp))
    end if
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


  subroutine multiply_psi_by_vKB(nbands, list_iGk, psi, dvnl_psi)
    !INTW project: KB projection by wave functions.
    !

    use intw_useful_constants, only: cmplx_0, cmplx_i
    use intw_reading, only: nat, nspin, nG_max, lspinorb, ityp
    use intw_pseudo, only: upf

    implicit none

    integer, intent(in)             :: nbands, list_iGk(nG_max)
    complex(kind=dp), intent(inout) :: psi(nG_max,nbands,nspin), dvnl_psi(nG_max,nbands,nspin)

    complex(kind=dp)                :: projec_d(nat,nhm,nspin)
    integer                         :: ispin, jspin, na, ih, ig
    integer                         :: nt, ikb, iGk, iband


    do iband = 1, nbands
      !
      projec_d = cmplx_0
      !
      ikb = 0
      do na = 1, nat
        do ih = 1, nh(ityp(na))
          !
          ikb = ikb + 1
          !
          do ig = 1, nG_max
            !
            iGk = list_iGk(iG)
            if (iGk==0) exit
            !
            do ispin = 1, nspin
              ! DKB(nhm,nhm,nspin,nspin,ntyp)
              projec_d(na,ih,ispin) = projec_d(na,ih,ispin) + conjg(vkb(iG,ikb)) * psi(iG,iband,ispin)
            end do
            !
          end do !ig
          !
        end do !ih
      end do !na

      ikb = 0
      do na = 1, nat
        do ih = 1, nh(ityp(na))
          !
          ikb = ikb + 1
          !
          do ig = 1, nG_max
            !
            iGk = list_iGk(iG)
            if (iGk==0) exit
            !
            if ( upf(nt)%has_so .and. lspinorb ) then
              !
              do ispin = 1, nspin
                do jspin = 1, nspin
                  dvnl_psi(iG,iband,ispin) = dvnl_psi(iG,iband,ispin) &
                      + DKB(ih,ih,ispin,jspin,ityp(na)) * projec_d(na,ih,jspin) * vkb(iG,ikb)
                end do
              end do
              !
            else ! no SO
              !
              do ispin = 1, nspin
                dvnl_psi(iG,iband,ispin) = dvnl_psi(iG,iband,ispin) &
                    + DKB(ih,ih,ispin,ispin,ityp(na)) * projec_d(na,ih,ispin) * vkb(iG,ikb)
              end do
            end if
            !
          end do !ig
          !
        end do!ih
      end do! na
      !
    end do !iband

  end subroutine multiply_psi_by_vKB


  subroutine multiply_psi_by_dvKB(k_cryst, q_cryst, nbands, list_iGk, list_iGkq, psi, dvnl_psi)

    use intw_useful_constants, only: cmplx_0, cmplx_i
    use intw_reading, only: nat, nspin, nG_max, tpiba, bg, lspinorb, ntyp, ityp
    use intw_fft, only: gvec_cart
    use intw_pseudo, only: upf

    implicit none

    real(kind=dp), intent(in)       :: k_cryst(3), q_cryst(3)
    integer, intent(in)             :: nbands, list_iGk(nG_max), list_iGkq(nG_max)
    complex(kind=dp), intent(inout) :: psi(nG_max,nbands,nspin), dvnl_psi(nG_max,nbands,nspin,nspin,3*nat)

    complex(kind=dp)                :: Dij(nkb,nkb,nspin,nspin)
    complex(kind=dp)                :: projec_1(nkb,3,nspin), projec_2(nkb,3,nspin)
    complex(kind=dp)                :: Dij_projec_1(nkb,3,nspin,nspin), Dij_projec_2(nkb,3,nspin,nspin)

    real(kind=dp)                   :: k_cart(3), q_cart(3)
    integer                         :: nt, ntj, na, naj, ih, jh, ikb, jkb, ispin, jspin, ipol, ig
    integer                         :: imode, iband, iGk, iGkq


    k_cart = matmul(bg, k_cryst)
    q_cart = matmul(bg, q_cryst)

    ! build D matrix with ikb index
    Dij = cmplx_0
    ikb = 0
    do nt = 1, ntyp
      do na = 1, nat
        !
        if (ityp(na)==nt) then
          !
          do ih = 1, nh(ityp(na))
            ikb = ikb + 1
            jkb = 0
            do ntj = 1, ntyp
              do naj = 1, nat
                if (ityp(naj)==ntj) then
                  do jh = 1, nh(ityp(naj))
                    jkb = jkb + 1
                    if (na==naj) then
                      if (lspinorb) then
                        Dij(ikb,jkb,:,:) = DKB(ih,jh,:,:,ityp(na))
                      else
                        do ispin=1,nspin
                          Dij(ikb,jkb,ispin,ispin) = DKB(ih,jh,1,1,ityp(na))
                        end do
                      end if
                    end if
                  end do !jh
                end if
              end do !naj
            end do !ntj
          end do !ih
          !
        end if
        !
      end do !na
    end do !nt




    do iband = 1, nbands

      projec_1 = cmplx_0
      projec_2 = cmplx_0

      ! Asier: KB potentziala hurrengo eran emanik dago: sum_l |b(l)> <b(l)|
      !               beraz deribatuak bi gai dauzka:
      !               sum_l d|b(l,r)> <b(l,r)| + |b(l,r)> d<b(l,r)| ~ Fourier ~
      !               sum_l i(k+G)|b(l,G)> <b(l,G)| + |b(l,G)> <b(l,G)|
      !

      do ipol = 1, 3 ! Cart. coord.
        !
        ikb = 0
        do nt = 1, ntyp
          do na = 1, nat
            !
            if (ityp(na)==nt) then
              !
              do ih = 1, nh(ityp(na))
                !
                ikb=ikb+1
                !
                do ig = 1, nG_max
                !
                  iGk = list_iGk(iG)
                  if (iGk==0) exit
                  !
                  do ispin = 1, nspin

                    projec_1(ikb,ipol,ispin) = projec_1(ikb,ipol,ispin) &
                        + conjg(vkb(iG,ikb)) * psi(iG,iband,ispin)

                    projec_2(ikb,ipol,ispin) = projec_2(ikb,ipol,ispin) &
                        + conjg(vkb(iG,ikb)) * psi(iG,iband,ispin) * &
                          tpiba * cmplx_i * ( k_cart(ipol) + gvec_cart(ipol,iGk) )

                  end do !ispin
                  !
                end do !ig
                !
              end do !ih
              !
            end if
            !
          end do !na
        end do ! nt
        !
      end do !ipol


      ! multiplay the projections <\beta_j|\psi_n> by the matrix Dij
      Dij_projec_1 = cmplx_0
      Dij_projec_2 = cmplx_0
      do ipol = 1, 3 ! Cart. coord.
        !
        if (lspinorb) then
          do ispin = 1, nspin
            do jspin = 1, nspin
              Dij_projec_1(:,ipol,ispin,jspin) = matmul( Dij(:,:,ispin,jspin), projec_1(:,ipol,jspin) )
              Dij_projec_2(:,ipol,ispin,jspin) = matmul( Dij(:,:,ispin,jspin), projec_2(:,ipol,jspin) )
            end do !jspin
          end do !ispin
        else
          do ispin = 1, nspin
            Dij_projec_1(:,ipol,ispin,ispin) = matmul( Dij(:,:,ispin,ispin), projec_1(:,ipol,ispin) )
            Dij_projec_2(:,ipol,ispin,ispin) = matmul( Dij(:,:,ispin,ispin), projec_2(:,ipol,ispin) )
          end do !ispin
        end if
        !
      end do !ipol


      do ipol = 1, 3 ! Cart. coord.
        !
        ikb = 0
        do nt = 1, ntyp
          do na = 1, nat
            !
            if (ityp(na)==nt) then
              !
              imode = (na-1)*3 + ipol
              !
              do ih = 1, nh(ityp(na))
                !
                ikb = ikb + 1
                !
                do ig = 1, nG_max
                  !
                  iGkq = list_iGkq(iG)
                  if (iGkq==0) exit
                  !
                  if ( upf(nt)%has_so .and. lspinorb ) then
                    !
                    do ispin = 1, nspin
                      do jspin = 1, nspin

                        dvnl_psi(iG,iband,ispin,jspin,imode) = dvnl_psi(iG,iband,ispin,jspin,imode) &
                            + Dij_projec_2(ikb,ipol,ispin,jspin) * vkqb(iG,ikb)
                        dvnl_psi(iG,iband,ispin,jspin,imode) = dvnl_psi(iG,iband,ispin,jspin,imode) &
                            - Dij_projec_1(ikb,ipol,ispin,jspin) * vkqb(iG,ikb) * &
                              tpiba * cmplx_i * ( k_cart(ipol) + q_cart(ipol) + gvec_cart(ipol,iGkq) )

                      end do !jspin
                    end do !ispin
                    !
                  else ! no SO
                    !
                    do ispin = 1, nspin

                      dvnl_psi(iG,iband,ispin,ispin,imode) = dvnl_psi(iG,iband,ispin,ispin,imode) &
                          + Dij_projec_2(ikb,ipol,ispin,ispin) * vkqb(iG,ikb)
                      dvnl_psi(iG,iband,ispin,ispin,imode) = dvnl_psi(iG,iband,ispin,ispin,imode) &
                          - Dij_projec_1(ikb,ipol,ispin,ispin) * vkqb(iG,ikb) * &
                            tpiba * cmplx_i * ( k_cart(ipol) + q_cart(ipol) + gvec_cart(ipol,iGkq) )

                    end do !ispin
                    !
                  end if

                end do !ig
                !
              end do !ih
              !
            end if
            !
          end do !na
        end do !nt
        !
      end do !ipol

    end do !iband

  end subroutine multiply_psi_by_dvKB

end module intw_pseudo_non_local