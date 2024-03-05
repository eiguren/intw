!---------------------------------------------------------------------------------------------------------------------
subroutine force_linearity_G(nqf1,nqf2,nqf3,qfmesh,omega)
!---------------------------------------------------------------------------------------------------------------------

  use kinds, only: dp
  use intw_utility, only: find_k_1bz_and_g
  use intw_reading, only: nat, bg
  use intw_useful_constants, only: eps_6

  implicit none

  !I/O variables

  integer,intent(in) :: nqf1,nqf2,nqf3
  real(dp),intent(in) :: qfmesh(3,nqf1*nqf2*nqf3)
  real(dp),intent(inout) :: omega(3*nat,nqf1*nqf2*nqf3)

  !local variables

  real(dp) :: w_min,wqv,w_min_K,w_min_M,w_trshld_K,w_trshld_M
  real(dp) :: diff_min_K,diff_min_M,w_trshld,mod_trshld
  real(dp) :: qpoint(3),qpoint_cart(3),qcart_0(3),qcart_nqfmesh(3),qcart_gora(3),qcart_behera(3)
  real(dp) :: mod_0,mod_gora,mod_behera,mod_nqfmesh
  real(dp) :: mod_min,modul,mod_min_K,mod_min_M,mod_qcart,mod_trshld_K,mod_trshld_M
  integer :: iq_w_min_K,iq_w_min_M,iq_trshld_K,iq_trshld_M
  integer :: iq,nqfmesh,i,j,k
  real(dp) :: kpoint_in_1bz(3)
  integer :: GKQ_bz(3)


  nqfmesh=nqf1*nqf2*nqf3
  !
  ! We look for the minimum value of phonon modes
  !
  w_min=0.0d0
  !
  do iq=1,nqfmesh
     !
     wqv=omega(1,iq)
     !
     if (iq==1) w_min=wqv
     !
     if (wqv.lt.w_min) then
        !
        w_min=wqv
        qpoint(:)=qfmesh(:,iq)
        qpoint_cart=matmul(bg,qpoint)
        mod_min=sqrt(qpoint_cart(1)**2+qpoint_cart(2)**2+qpoint_cart(3)**2)
        !
     endif
     !
  enddo !iq
  !
  ! We look for the q-point phonon mode is again equal to 0 in M and K directions
  ! from Gamma
  !
  do iq=1,nqfmesh
     !
     wqv=omega(1,iq)
     qpoint(:)=qfmesh(:,iq)
     call find_k_1BZ_and_G(qpoint,nqf1,nqf2,nqf3,i,j,k,kpoint_in_1bz,GKQ_bz)
     qpoint_cart=matmul(bg,qpoint)
     modul=sqrt(qpoint_cart(1)**2+qpoint_cart(2)**2+qpoint_cart(3)**2)
     !
     if (modul.gt.8*mod_min) cycle
     !
     if (i==j) then
        !
        if (modul.gt.mod_min.and.wqv.gt.w_min.and.wqv.lt.eps_6) then
           !
           iq_w_min_K=iq
           w_min_K=wqv
           mod_min_K=sqrt(qpoint_cart(1)**2+qpoint_cart(2)**2+qpoint_cart(3)**2)
           !
        endif
        !
     endif
        !
     if (i==1) then
        !
        if (modul.gt.mod_min.and.wqv.gt.w_min.and.wqv.lt.eps_6) then
           !
           iq_w_min_M=iq
           w_min_M=wqv
           mod_min_M=sqrt(qpoint_cart(1)**2+qpoint_cart(2)**2+qpoint_cart(3)**2)
           !
        endif
        !
     endif
     !
  enddo !iq
  !
  ! We multiply this distance by 5, and we get the qpoint and energy from where
  ! we will apply our linearity of the lowest phonon mode
  !
  do iq=1,nqfmesh
     !
     qpoint(:)=qfmesh(:,iq)
     call find_k_1BZ_and_G(qpoint,nqf1,nqf2,nqf3,i,j,k,kpoint_in_1bz,GKQ_bz)
     qpoint_cart=matmul(bg,qpoint)
     mod_qcart=sqrt(qpoint_cart(1)**2+qpoint_cart(2)**2+qpoint_cart(3)**2)
     wqv=omega(1,iq)
     !
     if (i==j) then
        !
        if (iq==1) then
           !
           diff_min_K=abs(mod_qcart-3*mod_min_K)
           !
        elseif (iq.gt.1.and.abs(mod_qcart-3*mod_min_K).lt.diff_min_K) then
           !
           iq_trshld_K=iq
           diff_min_K=abs(mod_qcart-3*mod_min_K)
           w_trshld_K=wqv
           mod_trshld_K=mod_qcart
           !
        endif
        !
     endif
        !
     if (i==1) then
        !
        if (iq==1) then
           !
           diff_min_M=abs(mod_qcart-3*mod_min_M)
           !
        elseif (iq.gt.1.and.abs(mod_qcart-3*mod_min_M).lt.diff_min_M) then
           !
           iq_trshld_M=iq
           diff_min_M=abs(mod_qcart-3*mod_min_M)
           w_trshld_M=wqv
           mod_trshld_M=mod_qcart
           !
        endif
        !
     endif
     !
  enddo !iq
  !
  ! We get the lowest treshold value of them
  !
  if (w_trshld_M.lt.w_trshld_K) then
     !
     w_trshld=w_trshld_M
     mod_trshld=mod_trshld_M
     !
  else
     !
     w_trshld=w_trshld_K
     mod_trshld=mod_trshld_K
     !
  endif
  !
  ! We apply our forcing linearity. Be carefull, linearity must be applied in the
  ! 4 edges of the SBZ
  !
  qpoint=qfmesh(:,1)
  qcart_0=matmul(bg,qpoint)
  !
  qpoint=qfmesh(:,1)
  qpoint(1)=qpoint(1)+1.d0
  qcart_behera=matmul(bg,qpoint)
  !
  qpoint=qfmesh(:,1)
  qpoint(2)=qpoint(2)+1.d0
  qcart_gora=matmul(bg,qpoint)
  !
  qpoint=qfmesh(:,nqfmesh)
  qcart_nqfmesh=matmul(bg,qpoint)
  !
  do iq=1,nqfmesh
     !
     qpoint(:)=qfmesh(:,iq)
     wqv=omega(1,iq)
     qpoint_cart=matmul(bg,qpoint)
     !
     mod_0 = (qpoint_cart(1)-qcart_0(1))**2 + &
             (qpoint_cart(2)-qcart_0(2))**2 + &
             (qpoint_cart(3)-qcart_0(3))**2
     mod_0 = sqrt( mod_0 )

     mod_behera = (qpoint_cart(1)-qcart_behera(1))**2 + &
                  (qpoint_cart(2)-qcart_behera(2))**2 + &
                  (qpoint_cart(3)-qcart_behera(3))**2
     mod_behera = sqrt( mod_behera )

     mod_gora = (qpoint_cart(1)-qcart_gora(1))**2 + &
                (qpoint_cart(2)-qcart_gora(2))**2 + &
                (qpoint_cart(3)-qcart_gora(3))**2
     mod_gora = sqrt( mod_gora )

     mod_nqfmesh = (qpoint_cart(1)-qcart_nqfmesh(1))**2 + &
                   (qpoint_cart(2)-qcart_nqfmesh(2))**2 + &
                   (qpoint_cart(3)-qcart_nqfmesh(3))**2
     mod_nqfmesh = sqrt( mod_nqfmesh )
     !
     if (mod_0.gt.10*mod_min.and.mod_gora.gt.10*mod_min.and.mod_behera.gt.10*mod_min.and.mod_nqfmesh.gt.10*mod_min) cycle
     !
     if (wqv.lt.w_trshld) then
        !
        if (mod_0.lt.mod_nqfmesh.and.mod_0.lt.mod_gora.and.mod_0.lt.mod_behera) then
           !
           omega(1,iq)=w_trshld*(mod_0/mod_trshld)
           !
        elseif (mod_behera.lt.mod_nqfmesh.and.mod_behera.lt.mod_gora.and.mod_behera.lt.mod_0) then
           !
           omega(1,iq)=w_trshld*(mod_behera/mod_trshld)
           !
        elseif (mod_gora.lt.mod_nqfmesh.and.mod_gora.lt.mod_behera.and.mod_gora.lt.mod_0) then
           !
           omega(1,iq)=w_trshld*(mod_gora/mod_trshld)
           !
        elseif (mod_nqfmesh.lt.mod_gora.and.mod_nqfmesh.lt.mod_behera.and.mod_nqfmesh.lt.mod_0) then
           !
           omega(1,iq)=w_trshld*(mod_nqfmesh/mod_trshld)
           !
        endif
        !
     endif
     !
  enddo !iq
  !
  return

end subroutine force_linearity_G
