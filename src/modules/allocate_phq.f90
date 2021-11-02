!-----------------------------------------------------------------------
subroutine allocate_phq
  !-----------------------------------------------------------------------
  !
  USE kinds, only : DP
  USE intw_reading, only : nat, ntyp, nbands, nr1,nr2,nr3
!haritz
!  USE intw_reading, ONLY : ngm, nk1, nk2, nk3
  USE intw_reading, ONLY : ngm
!haritz
  USE intw_becmod, ONLY:  becp, allocate_bec_type, deallocate_bec_type
  USE intw_pseudo, ONLY: nkb
  USE intw_pseudo, ONLY: nhm

  USE intw_phus, ONLY : becp1, alphap
  USE intw_eqv, ONLY :  vlocq

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

!Peio
!The variable we use instead of num_bands
  integer :: nbands_loc
!Peio

  !
  nbands_loc=num_bands
  allocate (eigqts(nat))
  allocate (vlocq ( ngm , ntyp))
  allocate (becp1(nk1*nk2*nk3))
  allocate (alphap(3,nk1*nk2*nk3))
  allocate (dvloc(nr1*nr2*nr3))


  do ik=1,nk1*nk2*nk3

     call allocate_bec_type ( nkb, nbands_loc, becp1(ik) )
!     call allocate_bec_type ( nkb, nbands, becp1(ik) )
     DO ipol=1,3
        call allocate_bec_type ( nkb, nbands_loc, alphap(ipol,ik) )
!        call allocate_bec_type ( nkb, nbands, alphap(ipol,ik) )
     ENDDO
  end do !ik

  CALL allocate_bec_type ( nkb, nbands_loc, becp )
!  CALL allocate_bec_type ( nkb, nbands, becp )
  
  return
end subroutine allocate_phq

!-----------------------------------------------------------------------
subroutine deallocate_phq
  !-----------------------------------------------------------------------
  !
  USE kinds, only : DP
  USE intw_reading, only : nat, ntyp, nbands, ngm

  USE intw_becmod, ONLY:  becp, deallocate_bec_type
  USE intw_pseudo, ONLY: nkb
  USE intw_pseudo, ONLY: nhm

  USE intw_phus, ONLY : becp1, alphap
  USE intw_eqv, ONLY :  vlocq

  USE intw_ph, ONLY : eigqts,dvloc  

  implicit none
  INTEGER :: ipol
  !
  deallocate (eigqts)
  deallocate (vlocq )
  deallocate (becp1)
  deallocate (alphap)
  deallocate (dvloc) 

     call deallocate_bec_type ( becp1(1) )
     DO ipol=1,3
        call deallocate_bec_type (alphap(ipol,1) )
     ENDDO

  CALL deallocate_bec_type (  becp )

  return
end subroutine deallocate_phq
