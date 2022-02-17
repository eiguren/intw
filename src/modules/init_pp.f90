
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine init_pp
  !----------------------------------------------------------------------
  !
  USE kinds,        ONLY : DP
  USE intw_pseudo,   ONLY : lmaxx
  USE intw_useful_constants,    ONLY : fpi, sqrt2
  !  ASIER
  !  USE intw_atom,         ONLY : rgrid
  USE intw_reading,    ONLY : ntyp, volume0, lspinorb, tpiba

  !haritz

  USE intw_pseudo,           ONLY : nqxq, dq, nqx, tab, tab_d2y, qrad
  USE splinelib
  USE intw_pseudo,         ONLY : nhtol, nhtoj, nhtolm, ijtoh, dvan, indv,&
       dvan_so
  USE intw_pseudo,   ONLY : upf, nbetam, nh, nhm, lmaxkb

  USE intw_fft,        ONLY : gvec_cart, gg
  USE intw_spin_orb,     ONLY : rot_ylm, fcoef
  USE mcf_spline

  !  ASIER
  USE intw_utility, ONLY: intgr_spline_gaussq, simpson
  !
  implicit none
  !
  !     here a few local variables
  !
  integer :: nt, ih, jh, nb, mb, ijv, l, m, ir, iq, is, startq, &
       lastq, ilast, ndm
  ! various counters
  real(DP), allocatable :: aux (:), aux1 (:), besr (:), qtot (:,:)
  ! various work space
  real(DP) :: prefr, pref, q, qi
  ! the prefactor of the q functions
  ! the prefactor of the beta functions
  ! the modulus of g for each shell
  ! q-point grid for interpolation
  ! the spherical harmonics
  real(DP) ::  vqint, j
  ! interpolated value
  ! J=L+S (noninteger!)
  integer :: n1, m0, m1, n, li, mi, vi, vj, ijs, is1, is2, &
       lk, mk, vk, kh, lh
  integer, external :: sph_ind
  complex(DP) :: coeff, qgm(1)
  real(DP) :: spinor, ji, jk
  !
  real(DP), allocatable :: xdata(:)
  real(DP) :: d1
  !
  !
  !    Initialization of the variables
  !
  ndm = MAXVAL ( upf(:)%kkbeta )
  allocate (aux ( ndm))
  allocate (aux1( ndm))
  allocate (besr( ndm))
  allocate (qtot( ndm , nbetam*(nbetam+1)/2 ))
  !
  ! the following prevents an out-of-bound error: upf(nt)%nqlc=2*lmax+1
  ! but in some versions of the PP files lmax is not set to the maximum
  ! l of the beta functions but includes the l of the local potential
  !
  prefr = fpi / volume0
  if (lspinorb) then
     !
     !  In the spin-orbit case we need the unitary matrix u which rotates the
     !  real spherical harmonics and yields the complex ones.
     !
     rot_ylm=(0.d0,0.d0)
     l=lmaxx
     rot_ylm(l+1,1)=(1.d0,0.d0)
     do n1=2,2*l+1,2
        m=n1/2
        n=l+1-m
        rot_ylm(n,n1)=CMPLX((-1.d0)**m/sqrt2,0.0_dp,kind=DP)
        rot_ylm(n,n1+1)=CMPLX(0.d0,-(-1.d0)**m/sqrt2,kind=DP)
        n=l+1+m
        rot_ylm(n,n1)=CMPLX(1.0_dp/sqrt2,0.d0,kind=DP)
        rot_ylm(n,n1+1)=CMPLX(0.d0, 1.0_dp/sqrt2,kind=DP)
     enddo
     fcoef=(0.d0,0.d0)
     dvan_so = (0.d0,0.d0)
  else
     dvan = 0.d0
  endif
  !
  !   For each pseudopotential we initialize the indices nhtol, nhtolm,
  !   nhtoj, indv, and if the pseudopotential is of KB type we initialize the
  !   atomic D terms
  !
  do nt = 1, ntyp
     ih = 1
     do nb = 1, upf(nt)%nbeta
        l = upf(nt)%lll (nb)
        do m = 1, 2 * l + 1
           nhtol (ih, nt) = l
           nhtolm(ih, nt) = l*l+m
           indv  (ih, nt) = nb
           ih = ih + 1
        enddo
     enddo
     if ( upf(nt)%has_so ) then
        ih = 1
        do nb = 1, upf(nt)%nbeta
           l = upf(nt)%lll (nb)
           j = upf(nt)%jjj (nb)
           do m = 1, 2 * l + 1
              nhtoj (ih, nt) = j
              ih = ih + 1
           enddo
        enddo
     endif
     ! ijtoh map augmentation channel indexes ih and jh to composite
     ! "triangular" index ijh
     ijtoh(:,:,nt) = -1
     ijv = 0
     do ih = 1,nh(nt)
        do jh = ih,nh(nt)
           ijv = ijv+1
           ijtoh(jh,ih,nt) = ijv
        enddo
     enddo
     !
     !    From now on the only difference between KB and US pseudopotentials
     !    is in the presence of the q and Q functions.
     !
     !    Here we initialize the D of the solid
     !
     if (upf(nt)%has_so.and.lspinorb) then
        !
        !  first calculate the fcoef coefficients
        !
        do ih = 1, nh (nt)
           li = nhtol(ih, nt)
           ji = nhtoj(ih, nt)
           mi = nhtolm(ih, nt)-li*li
           vi = indv (ih, nt)
           do kh=1,nh(nt)
              lk = nhtol(kh, nt)
              jk = nhtoj(kh, nt)
              mk = nhtolm(kh, nt)-lk*lk
              vk = indv (kh, nt)
              if (li == lk .and. abs(ji-jk) < 1.d-7) then
                 do is1=1,2
                    do is2=1,2
                       coeff = (0.d0, 0.d0)
                       do m=-li-1, li
                          m0= sph_ind(li,ji,m,is1) + lmaxx + 1
                          m1= sph_ind(lk,jk,m,is2) + lmaxx + 1
                          coeff=coeff + rot_ylm(m0,mi)*spinor(li,ji,m,is1)* &
                               CONJG(rot_ylm(m1,mk))*spinor(lk,jk,m,is2)
                       enddo
                       fcoef(ih,kh,is1,is2,nt)=coeff
                    enddo
                 enddo
              endif
           enddo
        enddo
        !
        !   and calculate the bare coefficients
        !
        do ih = 1, nh (nt)
           vi = indv (ih, nt)
           do jh = 1, nh (nt)
              vj = indv (jh, nt)
              ijs=0
              do is1=1,2
                 do is2=1,2
                    ijs=ijs+1
                    dvan_so(ih,jh,ijs,nt) = upf(nt)%dion(vi,vj) * &
                         fcoef(ih,jh,is1,is2,nt)
                    if (vi.ne.vj) fcoef(ih,jh,is1,is2,nt)=(0.d0,0.d0)
                 enddo
              enddo
           enddo
        enddo
     else
        do ih = 1, nh (nt)
           do jh = 1, nh (nt)
              if (nhtol (ih, nt) == nhtol (jh, nt) .and. &
                   nhtolm(ih, nt) == nhtolm(jh, nt) ) then
                 ir = indv (ih, nt)
                 is = indv (jh, nt)
                 if (lspinorb) then
                    dvan_so (ih, jh, 1, nt) = upf(nt)%dion (ir, is)
                    dvan_so (ih, jh, 4, nt) = upf(nt)%dion (ir, is)
                 else
                    dvan (ih, jh, nt) = upf(nt)%dion (ir, is)
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
  !         vqint = intgr_spline_gaussq( upf(nt)%r, aux )
  pref = fpi / sqrt (volume0)
  tab (:,:,:) = 0.d0
  do nt = 1, ntyp
     do nb = 1, upf(nt)%nbeta
        l = upf(nt)%lll (nb)
        !ASIER 29/07/2021
        do iq = 1,nqx!startq, lastq
           qi = (iq - 1) * dq
           !ASIER 29/07/2021
            call sph_bes (upf(nt)%kkbeta, upf(nt)%r, qi, l, besr)
           !call sph_bes (upf(nt)%kkbeta, rgrid(nt)%r, qi, l, besr)
           do ir = 1, upf(nt)%kkbeta
              !ASIER 29/07/2021
              aux (ir) = upf(nt)%beta (ir, nb) * besr (ir) * upf(nt)%r(ir)
              !aux (ir) = upf(nt)%beta (ir, nb) * besr (ir) * rgrid(nt)%r(ir)
           enddo
           !ASIER 29/07/2021
           !call simpson (upf(nt)%kkbeta, aux, upf(nt)%rab, vqint)

           !Integrating by spline + gauss 2. order
           vqint  =  intgr_spline_gaussq( upf(nt)%r(1:upf(nt)%kkbeta), aux )

           tab (iq, nb, nt) = vqint * pref

        enddo
     enddo
  enddo

  ! initialize spline interpolation
     allocate( xdata(nqx) )
     do iq = 1, nqx
        xdata(iq) = (iq - 1) * dq
     enddo
     do nt = 1, ntyp
        do nb = 1, upf(nt)%nbeta
           d1 = (tab(2,nb,nt) - tab(1,nb,nt)) / dq
           !call spline(xdata, tab(:,nb,nt), 0.d0, d1, tab_d2y(:,nb,nt))
            call spline_mcf(xdata,tab(:,nb,nt), nqx, tab_d2y(:,nb,nt))
        enddo
     enddo
     deallocate(xdata)

  if (allocated(xdata)) deallocate(xdata)

  deallocate (qtot)
  deallocate (besr)
  deallocate (aux1)
  deallocate (aux)

  return
end subroutine init_pp

!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
function sph_ind(l,j,m,spin)
  ! This function calculates the m index of the spherical harmonic
  ! in a spinor with orbital angular momentum l, total angular
  ! momentum j, projection along z of the total angular momentum m+-1/2.
  ! Spin selects the up (spin=1) or down (spin=2) coefficient.
  !
  use kinds
  use intw_utility, ONLY : errore
  implicit none

  integer :: sph_ind
  integer :: l, &            ! orbital angular momentum
       m, &            ! projection of the total angular momentum+-1/2
       spin            ! 1 or 2 select the component

  real(DP) :: j         ! total angular momentum

  if (spin.ne.1.and.spin.ne.2) call errore('sph_ind','spin direction unknown',1)
  if (m.lt.-l-1.or.m.gt.l) call errore('sph_ind','m not allowed',1)

  if (abs(j-l-0.5d0).lt.1.d-8) then
     if (spin.eq.1) sph_ind= m
     if (spin.eq.2) sph_ind= m+1
  elseif (abs(j-l+0.5d0).lt.1.d-8) then
     if (m.lt.-l+1) then
        sph_ind=0
     else
        if (spin.eq.1) sph_ind= m-1
        if (spin.eq.2) sph_ind= m
     endif
  else
     write(6,*) l, j
     call errore('sph_ind','l and j not compatible',1)
  endif
  if (sph_ind.lt.-l.or.sph_ind.gt.l) sph_ind=0

  return
end function sph_ind

function spinor(l,j,m,spin)
  ! This function calculates the numerical coefficient of a spinor
  ! with orbital angular momentum l, total angular momentum j,
  ! projection along z of the total angular momentum m+-1/2. Spin selects
  ! the up (spin=1) or down (spin=2) coefficient.

  use kinds
  use intw_utility, ONLY : errore

  implicit none

  real(DP) :: spinor
  integer :: l, &            ! orbital angular momentum
       m, &            ! projection of the total angular momentum+-1/2
       spin            ! 1 or 2 select the component

  real(DP) :: j         ! total angular momentum
  real(DP) :: denom     ! denominator

  if (spin.ne.1.and.spin.ne.2) call errore('spinor','spin direction unknown',1)
  if (m.lt.-l-1.or.m.gt.l) call errore('spinor','m not allowed',1)

  denom=1.d0/(2.d0*l+1.d0)
  if (abs(j-l-0.5d0).lt.1.d-8) then
     if (spin.eq.1) spinor= sqrt((l+m+1.d0)*denom)
     if (spin.eq.2) spinor= sqrt((l-m)*denom)
  elseif (abs(j-l+0.5d0).lt.1.d-8) then
     if (m.lt.-l+1) then
        spinor=0.d0
     else
        if (spin.eq.1) spinor= sqrt((l-m+1.d0)*denom)
        if (spin.eq.2) spinor= -sqrt((l+m)*denom)
     endif
  else
     call errore('spinor','j and l not compatible',1)
  endif

  return
end function spinor
