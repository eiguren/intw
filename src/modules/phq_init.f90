!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE phq_init( xq)
  !----------------------------------------------------------------------------
  !
  !ASIER july 21
  !USE intw_atom,             ONLY :  msh, rgrid
  USE intw_ph,               ONLY :  eigqts
  USE kinds, ONLY : dp
  USE intw_reading, ONLY : nat, tau, ntyp, tpiba2, ngm, volume0 
  USE intw_useful_constants,      ONLY: tpi 
  USE intw_pseudo, ONLY :  vlocq
  USE intw_pseudo,       ONLY : upf
  USE intw_fft, ONLY : gvec_cart

  IMPLICIT NONE
   REAL(DP), INTENT(IN)  :: xq(3) !kartestarrak.
  !
  ! ... local variables
  !
  INTEGER :: nt, ik, ikq, ipol, ibnd, ikk, na, ig, irr, imode0
    ! counter on atom types
    ! counter on k points
    ! counter on k+q points
    ! counter on polarizations
    ! counter on bands
    ! index for wavefunctions at k
    ! counter on atoms
    ! counter on G vectors
  REAL(DP) :: arg
    ! the argument of the phase
  COMPLEX(DP), EXTERNAL :: zdotc
  !
   
  DO na = 1, nat
     !
     arg = ( xq(1) * tau(1,na) + &
             xq(2) * tau(2,na) + &
             xq(3) * tau(3,na) ) * tpi
     !
     eigqts(na) = CMPLX( COS( arg ), - SIN( arg ) ,kind=DP)
     !
  END DO
  !
  !
  vlocq(:,:) = 0.D0
  !
 
  DO nt = 1, ntyp
     !ASIER july 21
     !CALL setlocq( xq, rgrid(nt)%mesh, msh(nt), rgrid(nt)%rab, rgrid(nt)%r,&
     !              upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngm, gvec_cart, volume0, &
     !              vlocq(1,nt) )

     CALL setlocq( xq, upf(nt)%mesh, upf(nt)%mesh, upf(nt)%rab, upf(nt)%r,&
                   upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngm, gvec_cart, volume0, &
                   vlocq(1,nt) )

     !
  END DO
  !
  !
  RETURN
  !
END SUBROUTINE phq_init
