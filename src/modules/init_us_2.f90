!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine init_us_2 (npw_, npwx_, igk_, qpoint_, vkb_)
  !----------------------------------------------------------------------
  !
  !   Calculates beta functions (Kleinman-Bylander projectors), with
  !   structure factor, for all atoms, in reciprocal space
  !
  USE kinds, ONLY: dp                                       !-QE :USE kinds,      ONLY : DP
!haritz
!  USE intw_reading, ONLY: nat, ntyp, ityp, tau, bg, ngm     !-QE :USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
!  USE intw_pseudo, ONLY: tpiba                              !-QE :USE cell_base,  ONLY : tpiba
  USE intw_reading, ONLY: nat, ntyp, ityp, tau, bg, ngm, tpiba
!haritz
  USE intw_useful_constants, ONLY : tpi, cmplx_0            !-QE :USE constants,  ONLY : tpi
  USE intw_pseudo, ONLY: nqx,nqxq, dq, tab, tab_d2y, spline_ps
  USE splinelib
  USE intw_pseudo, ONLY: nkb, vkb, nhtol, nhtolm, indv
  USE intw_pseudo, ONLY : upf, lmaxkb, nhm, nh
  USE intw_fft, ONLY: gvec_cart, eigts1, eigts2, eigts3, mill

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
  vkb1=0.d0
  !
  if (lmaxkb.lt.0) return
  !
  allocate (vkb1( npw_,nhm))    
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
  call intw_ylmr2 ((lmaxkb+1)**2, npw_, gk, qg, ylm)
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
  if (spline_ps) then
     !
     allocate(xdata(nqx))
     !
     do iq = 1, nqx
        !
        xdata(iq) = (iq - 1) * dq
        !
     enddo !ig
  endif !spline_ps
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
           if (spline_ps) then
              !
              vq(ig) = splint(xdata, tab(:,nb,nt), tab_d2y(:,nb,nt), qg(ig))
              !
           else
              !             
              px = qg (ig) / dq - int (qg (ig) / dq)
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = INT( qg (ig) / dq ) + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              !
              vq (ig) = tab (i0, nb, nt) * ux * vx * wx / 6.d0 + &
                        tab (i1, nb, nt) * px * vx * wx / 2.d0 - &
                        tab (i2, nb, nt) * px * ux * wx / 2.d0 + &
                        tab (i3, nb, nt) * px * ux * vx / 6.d0
              !
           endif !spline_ps
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
              sk (ig) = eigts1 (mill(1,igk_(ig)), na) * &
                        eigts2 (mill(2,igk_(ig)), na) * &
                        eigts3 (mill(3,igk_(ig)), na)
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
  deallocate (gk)
  deallocate (ylm)
  deallocate (vq)
  deallocate (qg)
  deallocate (sk)
  deallocate (vkb1)
  !
  return

end subroutine init_us_2

