!-----------------------------------------------------------------------
subroutine allocate_phq
  !-----------------------------------------------------------------------
  !JLB: Commented out all bands references, not used.
  !     Many variables are declared and not used, I don't understand why.
  !     Discuss with Haritz/Asier, and then remove and clean up.
  !
  USE intw_reading, only : nat, ntyp, nr1, nr2, nr3
  USE intw_reading, ONLY : ngm
  USE intw_pseudo, ONLY : vlocq
  USE intw_ph, ONLY : dvloc

  implicit none

  allocate(vlocq (ngm, ntyp))
  allocate(dvloc(nr1*nr2*nr3))

end subroutine allocate_phq

!-----------------------------------------------------------------------
subroutine deallocate_phq
  !-----------------------------------------------------------------------
  !
  USE intw_pseudo, ONLY : vlocq
  USE intw_ph, ONLY : dvloc

  implicit none
  !
  deallocate(vlocq)
  deallocate(dvloc)

  return
end subroutine deallocate_phq
