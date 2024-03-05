!-----------------------------------------------------------------------
subroutine allocate_nlpot ()
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
  USE intw_pseudo, only: tab, tab_d2y, dq, nqx, indv, nhtol, nhtolm, &
                         ijtoh, vkb, vkqb, nkb, nhtoj, DKB, upf, &
                         lmaxkb, nh, nhm, nbetam

  implicit none

  !local variables

  integer :: nt, nb, na

  real(dp),parameter :: cell_factor=1.0_dp


  lmaxkb = - 1
  DO nt = 1, ntyp
     !
     nh (nt) = 0
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
  nhm = MAXVAL (nh (1:ntyp))
  nbetam = MAXVAL (upf(:)%nbeta)
  !
  ! Number of beta functions
  !
  nkb = 0
  do na = 1, nat
     nt = ityp(na)
     nkb = nkb + nh (nt)
  enddo

  allocate (indv( nhm, ntyp))
  allocate (nhtol(nhm, ntyp))
  allocate (nhtolm(nhm, ntyp))
  allocate (nhtoj(nhm, ntyp))
  allocate (ijtoh(nhm, nhm, ntyp))
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
  nqx = INT( (sqrt (2*ecutwfc) / dq + 4) * cell_factor ) ! x2 zeren Ry -> Hartree egin behar

  allocate (tab( 2*nqx , nbetam , ntyp)) ! ASIER originala nqx, errepasatu

  ! d2y is for the cubic splines
  allocate (tab_d2y( nqx , nbetam , ntyp))

  if (allocated(vkb)) deallocate (vkb)
  if (allocated(vkqb)) deallocate (vkqb)

  if (nkb > 0) then
      allocate (vkb ( nG_max,  nkb))
      allocate (vkqb( nG_max,  nkb))
  end if
  !
  return

end subroutine allocate_nlpot
!*************************************************************************
!-------------------------------------------------------------------------
subroutine set_nqxq (qpoint)
  !
  USE kinds, only: dp
  USE intw_pseudo,  ONLY: dq, nqxq
  USE intw_reading, ONLY: ecutrho, tpiba2, bg

  implicit none

  !I/O variables

  real(dp),intent(in) :: qpoint(3)
  real(dp) :: xq(3),qnorm
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
