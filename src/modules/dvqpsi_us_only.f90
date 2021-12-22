!------------------------------------------------------------------------------
subroutine dvqpsi_us_only (xk,xq,nat,nG_max,nbands,npol,igk,igkq,wfc_k,dvpsi)
!------------------------------------------------------------------------------
!
!========================================================================
! This routine calculates dV_bare/dtau * psi for one perturbation       !
! with a given q. The displacements are described by a vector uact.     !
! The result is stored in dvpsi. The routine is called for each k point !
! and for each pattern u. It computes simultaneously all the bands.     !
! This routine implements Eq. B29 of PRB 64, 235118 (2001).             !
! Only the contribution of the nonlocal potential is calculated here.   !
!========================================================================

  use kinds, only: DP
  use intw_reading, only: tpiba
  use intw_fft, only: gvec_cart
  use intw_useful_constants, only: cmplx_0,cmplx_i
  use intw_reading, only:  ityp, ntyp, noncolin 
  use intw_pseudo, only: nkb, vkb, vkqb, nhtol, nhtoj, nhtolm, ijtoh, dvan, indv, dvan_so
  use intw_pseudo, only: nh, nhm
!
  use intw_phus, only: alphap, becp1
  use intw_becmod, only: calbec

  implicit none

  !I/O variables

  integer,intent(in) :: nat,nG_max,nbands,npol,igk(ng_max),igkq(ng_max)
  real(dp),intent(in) :: xq(3),xk(3)
  complex(dp),intent(in) :: wfc_k(nG_max,nbands,npol)
  complex(dp),intent(inout) :: dvpsi(nG_max,nbands,npol,npol,3*nat)  

  !local variables
  
  complex(dp) :: aux1(nG_max*npol,nbands),uact(3*nat) !uact: pattern of displ.
  integer :: nsp, na, nb, mu, nu, ikk, ikq, ig, igg, nt, ibnd
  integer :: ijkb0, ikb, jkb, ih, jh, ipol, is, js, ijs, imode
  ! counter on atoms: na, nb
  ! counter on modes: mu, nu
  ! the point k: ikk
  ! the point k+q: ikq
  ! counter on G vectors: ig
  ! auxiliary counter on G vectors: igg
  ! counter on atomic types: nt
  ! counter on bands: ibnd
  ! auxiliary variable for counting: ijkb0
  ! counter on becp functions: ikb
  ! counter on becp functions: jkb
  ! counter on n index: ih
  ! counter on m index: jh
  ! counter on polarizations

  real(dp), parameter :: eps = 1.d-12
  complex(DP), allocatable :: ps1(:,:),ps2(:,:,:),aux(:),deff_nc(:,:,:,:)
  complex(DP), allocatable :: ps1_nc (:,:,:), ps2_nc (:,:,:,:)
  real(DP), allocatable :: deff(:,:,:)
  ! work space

  logical :: ok
!Peio
  integer :: nG,npw
!Peio

  nsp=ntyp
  !
  if (noncolin) then
     allocate (ps1_nc(nkb,npol,nbands))
     allocate (ps2_nc(nkb,npol,nbands,3))
  else
     allocate (ps1(nkb,nbands))
     allocate (ps2(nkb,nbands,3))
  end if
  allocate (aux(nG_max))
  !
  aux1=(0.d0,0.d0)
  !
  nG=0
  !
  ! How much G vectors for wfc_k
  !
  do iG=1,nG_max
     if (igk(iG)==0) exit
     nG=nG+1
  enddo
  npw=nG
  !
  ! We need wfc_k in QE fashion -> aux1
  !
  do ibnd=1,nbands
     do ig=1,nG_max
        do ipol=1,npol
           aux1(ig+(ipol-1)*nG_max,ibnd)=wfc_k(ig,ibnd,ipol)
        enddo !ipol
     enddo !ig
  enddo !ibnd
  !
  ! Calculation of becp1, necessary for the next
  !
  CALL calbec (npw,vkb,aux1,becp1(1) )
  !
  ! We prepare dwfc_k/dxi as QE (taking into account spinor too)
  ! Calculate alphap, necessary for the next
  !
  do ipol=1,3
     !
     aux1=(0.d0,0.d0)
     !
     do ibnd=1,nbands
        do ig=1,nG_max
           !
           if ((igk(ig)>0).and.(igk(ig)<nG_max)) then
              !
              aux1(ig,ibnd)=wfc_k(ig,ibnd,1)*tpiba*cmplx_i* &
                           (xk(ipol)+gvec_cart(ipol,igk(ig)))
              !
           else
              !
              aux1(ig,ibnd)=cmplx_0
              !
           endif
           !
        enddo
        !
        if (npol==2) then
           !
           do ig=1,nG_max
              !
              if ((igk(ig)>0).and.(igk(ig)<nG_max)) then
                 !
                 aux1(ig+ng_max,ibnd)=wfc_k(ig,ibnd,2)*tpiba*cmplx_i* &
                                     (xk(ipol)+gvec_cart(ipol,igk(ig)))
                 !
              else
                 !
                 aux1(ig+ng_max,ibnd)=cmplx_0
                 !
              endif
              !
           enddo !ig
        endif !npol
     enddo !ibnd
     !
     CALL calbec (nG_max,vkb,aux1,alphap(ipol,1))
     !
  enddo !ipol
  !
  ! We calculate dV_non_local|psi_k> for each mode
  !
  do imode=1,3*nat
     !
     uact(:)=(0.d0,0.d0)
     uact(imode)=(1.d0,0.d0) 
     !
     if (noncolin) then
        ps1_nc(:,:,:)=(0.d0,0.d0)
        ps2_nc(:,:,:,:)=(0.d0,0.d0)
     else
        ps1(:,:)=(0.d0,0.d0)
        ps2(:,:,:)=(0.d0,0.d0)
     end if
     !
     do ibnd=1,nbands
        !
        ijkb0=0
        !
        do nt=1,ntyp
           do na=1,nat
              !
              if (ityp(na).eq.nt) then
                 !
                 mu=3*(na - 1)
                 !
                 do ih=1,nh(nt)
                    !
                    ikb=ijkb0+ih
                    !
                    do jh=1,nh(nt)
                       !
                       jkb=ijkb0+jh
                       !
                       do ipol = 1, 3
                          !
                          if (abs(uact(mu+1))+abs(uact(mu+2))+abs(uact(mu+3))>eps) then
                             !
                             if (noncolin) then
                                !
                                ijs=0
                                !
                                do is=1,npol
                                   do js=1,npol
                                      !
                                      !if (is.ne.js) cycle
                                      !
                                      ijs=ijs+1
                                      !
                                      ps1_nc(ikb,is,ibnd)=ps1_nc(ikb,is,ibnd)+ &
                                             dvan_so(ih,jh,ijs,nt)*alphap(ipol,1)%nc(jkb,js,ibnd)* &
                                             uact(mu+ipol)
                                      !
!haritz
!                                      ps2_nc(ikb,is,ibnd,ipol)=ps2_nc(ikb,is,ibnd,ipol)+ &
!                                             dvan_so(ih,jh,ijs,nt)*becp1(1)%nc(jkb,js,ibnd)* &
!                                             -cmplx_i*uact(mu+ipol)*tpiba
                                      ps2_nc(ikb,is,ibnd,ipol)=ps2_nc(ikb,is,ibnd,ipol)+ &
                                             dvan_so(ih,jh,ijs,nt)*becp1(1)%nc(jkb,js,ibnd)* &
                                             (-cmplx_i*uact(mu+ipol)*tpiba)
!haritz
                                      !
                                   enddo !js
                                enddo !is
                                !
                             else !noncolin
                             !
                             ps1(ikb,ibnd)=ps1(ikb,ibnd)+dvan(ih,jh,nt)*alphap(ipol,1)%k(jkb,ibnd)* &
                                           uact(mu+ipol)
                             !
!haritz
!                             ps2(ikb,ibnd,ipol)=ps2(ikb,ibnd,ipol)+dvan(ih,jh,nt)*becp1(1)%k(jkb,ibnd)* &
!                                           -cmplx_i*uact(mu+ipol)*tpiba
                             ps2(ikb,ibnd,ipol)=ps2(ikb,ibnd,ipol)+dvan(ih,jh,nt)*becp1(1)%k(jkb,ibnd)* &
                                           (-cmplx_i*uact(mu+ipol)*tpiba)
!haritz
                             endif !noncolin
                             !
                          endif !uact>0
                          !
                       enddo !ipol
                    enddo !jh
                 enddo !ih
                 !
                 ijkb0=ijkb0+nh(nt)
                 !
              endif !ityp=type(na)
              !
           enddo !na
        enddo !nt
     enddo !nbands
     !
     ! This term is proportional to beta(k+q+G)
     !
     if (nkb.gt.0) then
        !
        if (noncolin) then
           !
           do is=1,npol
              !
              call zgemm ('N','N',ng_max,nbands,nkb,(1.d0,0.d0),vkqb,ng_max, &
                          ps1_nc(:,is,:),nkb,(1.d0,0.d0),dvpsi(:,:,is,is,imode),ng_max)
              !
           enddo !is
           !
        else
        !
        call zgemm ('N','N',ng_max,nbands,nkb,(1.d0, 0.d0),vkqb,ng_max, &
                    ps1,nkb,(1.d0,0.d0),dvpsi(:,:,1,1,imode),ng_max)
        !
        endif !noncolin
        !
     endif !nkb
     !
     ! This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
     !
     do ikb=1,nkb
        do ipol=1,3
           !
           ok=.false.
           !
           if (noncolin) then
              !
              do ibnd=1,nbands
                 !
                 ok=ok.or.(abs(ps2_nc(ikb,1,ibnd,ipol)).gt.eps).or. &
                          (abs(ps2_nc(ikb,2,ibnd,ipol)).gt.eps)
                 !
              enddo
              !
           else
              !
              do ibnd=1,nbands
                 !
                 ok=ok.or.(abs(ps2(ikb,ibnd,ipol)).gt.eps)
                 !
              enddo
              !
           endif !noncolin
           !
           if (ok) then
              !
              do ig=1,nG_max
                 !
                 igg=igkq(ig)
                 !
                 if ((igg>0).and.(igg<nG_max)) then 
                    !
                    aux(ig)=vkqb(ig,ikb)*(xk(ipol)+xq(ipol)+gvec_cart(ipol,igg))
                    !
                 else
                    !
                    aux(ig)=cmplx_0
                    !
                 endif !igg>0.and.<nG_max
                 !
              enddo
              !
              do ibnd=1,nbands
                 !
                 if (noncolin) then
                    !
                    do is=1,npol
                       !
                       call zaxpy(ng_max,ps2_nc(ikb,is,ibnd,ipol),aux,1,dvpsi(1,ibnd,is,is,imode),1)
                       !
                    enddo !is
                    !
                 else
                    !
                    call zaxpy(ng_max,ps2(ikb,ibnd,ipol),aux,1,dvpsi(1,ibnd,1,1,imode),1)
                    !
                 endif !noncolin
                 !
              enddo !ibnd
              !
           endif !ok
           !
        enddo !ipol
     enddo !ikb
  enddo !imode
  !
  deallocate (aux)
  if (noncolin) then
     deallocate (ps2_nc)
     deallocate (ps1_nc)
  else
     deallocate (ps2)
     deallocate (ps1)
  endif
  !
  return

end subroutine dvqpsi_us_only
