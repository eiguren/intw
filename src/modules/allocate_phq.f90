!-----------------------------------------------------------------------
subroutine allocate_phq
  !-----------------------------------------------------------------------
  !
  USE kinds, only : DP
  USE intw_reading, only : nat, ntyp, nbands, nr1,nr2,nr3
  USE intw_reading, ONLY : ngm
  USE intw_pseudo, ONLY: nkb
  USE intw_pseudo, ONLY: nhm

  USE intw_pseudo, ONLY :  vlocq

  USE intw_ph, ONLY : eigqts,dvloc
!haritz
  use intw_input_parameters, only: nk1, nk2, nk3
!haritz

!Peio
!Number of Bloch Original states for the wannierization
  use w90_parameters, only: num_bands 
!Peio

  implicit none
  INTEGER :: ipol, ik

  integer :: nbands_loc

  !
  nbands_loc=num_bands
  allocate (eigqts(nat))
  allocate (vlocq ( ngm , ntyp))
  allocate (dvloc(nr1*nr2*nr3))

  return
end subroutine allocate_phq

!-----------------------------------------------------------------------
subroutine deallocate_phq
  !-----------------------------------------------------------------------
  !
  USE kinds, only : DP
  USE intw_reading, only : nat, ntyp, nbands, ngm

  USE intw_pseudo, ONLY: nkb
  USE intw_pseudo, ONLY: nhm

  USE intw_pseudo, ONLY :  vlocq

  USE intw_ph, ONLY : eigqts,dvloc  

  implicit none
  INTEGER :: ipol
  !
  deallocate (eigqts)
  deallocate (vlocq )
  deallocate (dvloc) 

  return
end subroutine deallocate_phq
