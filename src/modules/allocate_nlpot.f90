!-----------------------------------------------------------------------
subroutine allocate_nlpot1 ()
  !-----------------------------------------------------------------------
  ! Adapted from QEspresso4.3.2 Asier&&Idoia
  !
  !     ngk           !  number of plane waves (for each k point)
  !     npwx          !  maximum number of plane waves
  !     nqx           !  number of points of the interpolation table
  !     nqxq          !  as above, for q-function interpolation table
  !
  !
  use kinds
  use intw_reading, only: nat,ntyp,ityp,ecutwfc,ecutrho,dual,gcutm,tpiba2,bg,ngm,nspin
  use intw_reading, only: noncolin,lspinorb,ng_max
  USE intw_pseudo, only: tab,tab_d2y,dq,nqx,nqxq,spline_ps
  use intw_pseudo, only: indv,nhtol,nhtolm,ijtoh,dvan,vkb,vkqb
  use intw_pseudo, only: nkb,nhtoj,becsum,dvan_so
  use intw_pseudo, only: upf,lmaxkb,nh,nhm,nbetam
  use intw_spin_orb, only: fcoef

  implicit none

  !local variables

  integer :: nwfcm,nsp,npwx 
  real(dp),parameter :: cell_factor=1.0_dp
  real(dp) :: qnorm

  nsp=ntyp
  !
  npwx=ng_max
  !
  call pre_init()
  !
  allocate (indv( nhm, nsp))    
  allocate (nhtol(nhm, nsp))    
  allocate (nhtolm(nhm, nsp))    
  allocate (nhtoj(nhm, nsp))    
  allocate (ijtoh(nhm, nhm, nsp))
  if (lspinorb) then
    allocate (dvan_so(nhm,nhm,nspin*nspin,nsp)) ! dvan_so needs 4 components, 1 for all spin couple possible
    allocate (fcoef(nhm,nhm,2,2,nsp))
  else
    allocate (dvan( nhm, nhm, nsp))
  endif
  !
  allocate (becsum( nhm * (nhm + 1)/2, nat, nspin))    
  !
  ! Calculate dimensions for array tab (including a possible factor
  ! coming from cell contraction during variable cell relaxation/MD)
  !
  nqx = INT( (sqrt (2*ecutwfc) / dq + 4) * cell_factor ) ! x2 zeren Ry -> Hartree egin behar

  allocate (tab( 2*nqx , nbetam , nsp)) ! ASIER originala nqx, errepasatu 

  ! d2y is for the cubic splines
  if (spline_ps) allocate (tab_d2y( nqx , nbetam , nsp))

  if (allocated(vkb)) deallocate (vkb)
  if (allocated(vkqb)) deallocate (vkqb)

  if (nkb > 0) then 
      allocate (vkb ( nG_max,  nkb))
      allocate (vkqb( nG_max,  nkb))  
  end if
  !
  return

end subroutine allocate_nlpot1
!*************************************************************************
!-------------------------------------------------------------------------
subroutine allocate_nlpot2 (qpoint)
  !-----------------------------------------------------------------------
  ! Adapted from QEspresso4.3.2 Asier&&Idoia
  !
  !     ngk           !  number of plane waves (for each k point)
  !     npwx          !  maximum number of plane waves
  !     nqx           !  number of points of the interpolation table
  !     nqxq          !  as above, for q-function interpolation table
  !
  !
  USE kinds
  USE intw_reading, ONLY: nat,ntyp,ityp,ecutwfc,ecutrho,dual,gcutm,tpiba2,bg
  USE intw_reading, ONLY: ngm,nspin,noncolin,ng_max
  USE intw_pseudo, ONLY: tab,tab_d2y,dq,nqx,nqxq,spline_ps
  USE intw_pseudo, ONLY: indv,nhtol,nhtolm,ijtoh,dvan, &
                       nkb,nhtoj,becsum,dvan_so
  USE intw_pseudo, ONLY: upf,lmaxkb,nh,nhm,nbetam
  USE intw_spin_orb, ONLY: fcoef
  !USE intw_pseudo

  implicit none

  !I/O variables 

  real(dp),intent(in) :: qpoint(3)

  !local variables
  
  integer :: nwfcm,nsp,npwx 
  real(dp), parameter :: cell_factor=1.0_dp
  real(dp) :: xq(3),qnorm

  nsp=ntyp
  !
  npwx=ng_max
  !
  xq = matmul(bg, qpoint) 
  !
  qnorm =sqrt(xq(1)**2+xq(2)**2+xq(3)**2)
  !
  ! This routine is called also by the phonon code, in which case it should
  ! allocate an array that includes q+G vectors up to |q+G|_max <= |Gmax|+|q|
  !
  nqxq = INT( (( (sqrt(2*gcutm) + qnorm ) / dq + 4) * cell_factor )) ! x2 zeren Ry -> Hartree egin behar
  !
  !
  return

end subroutine allocate_nlpot2
!*********************************************************************
!---------------------------------------------------------------------
SUBROUTINE pre_init()
  !-------------------------------------------------------------------
  !
  USE intw_reading,     ONLY : nat, ntyp, ityp
  USE intw_pseudo,  ONLY : upf, lmaxkb, nh, nhm, nbetam
  USE intw_pseudo,        ONLY : nkb 
  IMPLICIT NONE
  INTEGER :: na, nt, nb, nsp
  !
  !     calculate the number of beta functions for each atomic type
  !
  nsp=ntyp
  !
  lmaxkb = - 1
  DO nt = 1, nsp
     !
     nh (nt) = 0
     !
     ! do not add any beta projector if pseudo in 1/r fmt (AF)
     !
     DO nb = 1, upf(nt)%nbeta
        nh (nt) = nh (nt) + 2 * upf(nt)%lll(nb) + 1
        lmaxkb = MAX (lmaxkb, upf(nt)%lll(nb) )
     ENDDO
     !
  ENDDO
  !
  ! calculate the maximum number of beta functions
  !
  nhm = MAXVAL (nh (1:nsp))
  nbetam = MAXVAL (upf(:)%nbeta)
  !
  ! calculate the number of beta functions of the solid
  !
  nkb = 0
  do na = 1, nat
     nt = ityp(na)
     nkb = nkb + nh (nt)
  enddo
 
END SUBROUTINE pre_init
