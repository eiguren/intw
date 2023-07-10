!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine init_KB_projectors (npw_, npwx_, igk_, qpoint_, vkb_)
  !----------------------------------------------------------------------
  !
  !   Calculates beta functions (Kleinman-Bylander projectors), with
  !   structure factor, for all atoms, in reciprocal space
  !
  USE kinds, ONLY: dp                                       !-QE :USE kinds,      ONLY : DP
  USE intw_reading, ONLY: nat, ntyp, ityp, tau, bg, ngm, tpiba
  USE intw_useful_constants, ONLY : tpi, cmplx_0, cmplx_i            !-QE :USE constants,  ONLY : tpi
  USE intw_pseudo, ONLY: nqx,nqxq, dq, tab, tab_d2y
  USE splinelib
  USE intw_pseudo, ONLY: nkb, vkb, nhtol, nhtolm, indv
  USE intw_pseudo, ONLY : upf, lmaxkb, nhm, nh
  USE intw_fft, ONLY: gvec_cart

  USE mcf_spline

  implicit none

  !I/O variables

  integer,intent(in) :: npw_,npwx_           !number of PW's
  integer,intent(in) :: igk_(npwx_)          !G list of q vector
  real(dp),intent(in) :: qpoint_(3)          !q vector
  complex(dp),intent(out) :: vkb_(npwx_,nkb) !beta functions

  !local variables

  integer :: i0,i1,i2,i3, ig, l, lm, na, nt, nb, ih, jkb
  real(DP) :: px, ux, vx, wx, arg, bat, q_(3)
  real(DP),allocatable :: gk (:,:), qg (:), vq (:), ylm (:,:), vkb1(:,:)
  complex(DP) :: phase, pref
  complex(DP),allocatable :: sk(:)
  real(DP),allocatable :: xdata(:)
  integer :: iq

  q_=matmul(bg,qpoint_)
  !
  vkb_=cmplx_0
  allocate (vkb1( npw_,nhm))
  vkb1=0.d0
  !
  if (lmaxkb.lt.0) return
  !
  allocate (  sk( npw_))
  allocate (  qg( npw_))
  allocate (  vq( npw_))
  allocate ( ylm( npw_, (lmaxkb + 1) **2))
  allocate (  gk( 3, npw_))
  !
  qg=0.d0
  !
  do ig=1,npw_
     !
     if (igk_(ig)==0) exit
     !
     gk (1,ig) = q_(1) + gvec_cart(1,igk_(ig))
     gk (2,ig) = q_(2) + gvec_cart(2,igk_(ig))
     gk (3,ig) = q_(3) + gvec_cart(3,igk_(ig))
     !
     qg (ig) = gk(1,ig)**2 +  gk(2,ig)**2 + gk(3,ig)**2
     !
  enddo
  !
  call intw_real_ylmr2 ((lmaxkb+1)**2, npw_, gk, qg, ylm)
  !
  ! set now qg=|q+G| in atomic units
  !
  do ig=1,npw_
     !
     if (igk_(ig)==0) exit
     !
     qg(ig)=sqrt(qg(ig))*tpiba
     !
  enddo !ig
  !
  allocate(xdata(nqx))
     !
  do iq = 1, nqx
        !
     xdata(iq) = (iq - 1) * dq
        !
  enddo !ig
  !
  ! |beta_lm(q)> = (4pi/omega).Y_lm(q).f_l(q).(i^l).S(q)
  !
  jkb=0
  vq=0.d0
  !
  do nt=1,ntyp
     !
     ! calculate beta in G-space using an interpolation table f_l(q)=\int _0 ^\infty dr r^2 f_l(r) j_l(q.r)
     !
     do nb=1,upf(nt)%nbeta
        !
        do ig=1,npw_
           !
           if (igk_(ig)==0) exit
           !
           call splint_mcf (xdata, tab(:,nb,nt), tab_d2y(:,nb,nt), nqx, qg(ig), vq(ig))
        enddo !ig
        !
        ! add spherical harmonic part  (Y_lm(q)*f_l(q))
        !
        do ih=1,nh(nt)
           !
           if (nb.eq.indv(ih,nt)) then
              !
              l  = nhtol (ih,nt)
              lm = nhtolm(ih,nt)
              !
              do ig=1,npw_
                 !
                 if (igk_(ig)==0) exit
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
           ! q_ is cart. coordinates
           !
           arg = (q_(1) * tau (1, na) + &
                  q_(2) * tau (2, na) + &
                  q_(3) * tau (3, na) ) * tpi
           !
           phase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
           !
           do ig=1,npw_
              !
              if (igk_(ig)==0) exit
              !
              sk (ig) = exp ( -tpi*cmplx_i* (  gvec_cart(1,igk_(ig))*tau(1,na) &
                                             + gvec_cart(2,igk_(ig))*tau(2,na) &
                                             + gvec_cart(3,igk_(ig))*tau(3,na) ) )
              !
           enddo !ig
           !
           do ih=1,nh(nt)
              !
              jkb = jkb + 1
              pref = (0.d0, -1.d0) **nhtol (ih, nt) * phase
              !
              do ig = 1, npw_
                 !
                 if (igk_(ig)==0) exit
                 !
                 vkb_(ig, jkb) = vkb1 (ig,ih) * sk (ig) * pref
                 !
              enddo
           enddo
        endif
     enddo !n bet
  enddo !ntyp
  !
  deallocate (xdata)
  deallocate (gk)
  deallocate (ylm)
  deallocate (vq)
  deallocate (qg)
  deallocate (sk)
  deallocate (vkb1)
  !
  return

end subroutine init_KB_projectors

