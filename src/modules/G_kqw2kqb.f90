!----------------------------------------------------------------------------------
subroutine G_kqw2kqb(pel1,pel2,pph,g_wannier,g_elph)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 08/03/2016
!
!============================================================================================!
! 1) We transform matrix elements from Recipocal Wannier representation to Reciprocal        !
! Bloch representation, so we Wannier unrotate them:                                         !
!                                                                                            ! 
! g_int(q,k,m,n,si,sj,l)=Sum_m'_n'[Uel+(m',m,k+q)*g_wannier(q,k,m',n',si,sj,l)*Uel(n',n,k)]  !
!                                                                                            !
! 2) We transform matrix elements from Cartesian (Canonical) Phonon representation to        !
! real physical meaning Phonon representation:                                               !
!                                                                                            !
! g_elph(q,k,m,n,si,sj,nu)=Sum_l[Uph(l,nu)*g_wannier(q,k,m,n,si,sj,l)]                       !
!============================================================================================!

!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants
  use w90_parameters, only: num_wann

  implicit none

  !I/O variables

  complex(dp),intent(in) :: pel1(num_wann,num_wann),pel2(num_wann,num_wann)
  complex(dp),intent(in) :: pph(3*nat,3*nat)
  complex(dp),intent(in) :: g_wannier(num_wann,num_wann,npol,npol,3*nat)
  complex(dp),intent(inout) :: g_elph(num_wann,num_wann,npol,npol,3*nat)
  
  !local variables

  integer :: imode,jmode,ipol,jpol,ibnd,jbnd,iibnd,jjbnd
  complex(dp) :: g_int(num_wann,num_wann,npol,npol,3*nat)

  ! 1) Reciprocal Wannier -> Reciprocal Bloch (canonical phonons)
  !
  g_int=cmplx_0
  !
  do imode=1,3*nat
     do ipol=1,npol
        do jpol=1,npol
           !
           do ibnd=1,num_wann
              do jbnd=1,num_wann
                 do iibnd=1,num_wann
                    do jjbnd=1,num_wann
                       !
                       g_int(ibnd,jbnd,ipol,jpol,imode)=g_int(ibnd,jbnd,ipol,jpol,imode)+ &
                       conjg(pel1(iibnd,ibnd))*g_wannier(iibnd,jjbnd,ipol,jpol,imode)*pel2(jjbnd,jbnd)
                       !
                    enddo !jjbnd
                 enddo !!iibnd
              enddo !! jbnd
           enddo !ibnd
           !
        enddo !jpol
     enddo !ipol
  enddo !imode
  !
  ! 2) Canonical phonon -> Real phonon (Bloch gauge)
  !
  do ibnd=1,num_wann
     do jbnd=1,num_wann
        do ipol=1,npol
           do jpol=1,npol
              !
              do imode=1,3*nat
                 do jmode=1,3*nat
                    !
                    g_elph(ibnd,jbnd,ipol,jpol,imode)=g_elph(ibnd,jbnd,ipol,jpol,imode)+ &
                                       pph(jmode,imode)*g_int(ibnd,jbnd,ipol,jpol,jmode)
                    !
                 enddo !jmode
              enddo !imode
              !
           enddo !jpol
        enddo !ipol
     enddo !jbnd
  enddo !ibnd
  !
  return

end subroutine G_kqw2kqb
!------------------------------------------------------------------------
!************************************************************************
!------------------------------------------------------------------------
subroutine G_kqw2kqb_no_spin(pel1,pel2,pph,g_wannier,g_elph)
!----------------------------------------------------------------------------------
!
! Created by Peio G. Goiricelaya 08/03/2016
!
!============================================================================================!
! 1) We transform matrix elements from Recipocal Wannier representation to Reciprocal        !
! Bloch representation, so we Wannier unrotate them:                                         !
!                                                                                            ! 
! g_int(q,k,m,n,si,sj,l)=Sum_m'_n'[Uel+(m',m,k+q)*g_wannier(q,k,m',n',si,sj,l)*Uel(n',n,k)]  !
!                                                                                            !
! 2) We transform matrix elements from Cartesian (Canonical) Phonon representation to        !
! real physical meaning Phonon representation:                                               !
!                                                                                            !
! g_elph(q,k,m,n,si,sj,nu)=Sum_l[Uph(l,nu)*g_wannier(q,k,m,n,si,sj,l)]                       !
!============================================================================================!

!haritz
  use kinds, only: dp
!haritz
  use intw_input_parameters
  use intw_reading
  use intw_useful_constants
  use w90_parameters, only: num_wann

  implicit none

  !I/O variables

  complex(dp),intent(in) :: pel1(num_wann,num_wann),pel2(num_wann,num_wann)
  complex(dp),intent(in) :: pph(3*nat,3*nat)
  complex(dp),intent(in) :: g_wannier(num_wann,num_wann,3*nat)
  complex(dp),intent(inout) :: g_elph(num_wann,num_wann,3*nat)
  
  !local variables

  integer :: imode,jmode,ibnd,jbnd,iibnd,jjbnd
  complex(dp) :: g_int(num_wann,num_wann,3*nat)

  ! 1) Reciprocal Wannier -> Reciprocal Bloch (canonical phonons)
  !
  g_int=cmplx_0
  !
  do imode=1,3*nat
     !
     do ibnd=1,num_wann
        do jbnd=1,num_wann
           do iibnd=1,num_wann
              do jjbnd=1,num_wann
                 !
                 g_int(ibnd,jbnd,imode)=g_int(ibnd,jbnd,imode)+ &
                 conjg(pel1(iibnd,ibnd))*g_wannier(iibnd,jjbnd,imode)*pel2(jjbnd,jbnd)
                 !
              enddo !jjbnd
           enddo !iibnd
        enddo !jbnd
     enddo !ibnd
     !
  enddo !imode
  !
  ! 2) Canonical phonon -> Real phonon (Bloch gauge)
  !
  do ibnd=1,num_wann
     do jbnd=1,num_wann
        !
        do imode=1,3*nat
           do jmode=1,3*nat
              !
              g_elph(ibnd,jbnd,imode)=g_elph(ibnd,jbnd,imode)+ &
                                 pph(jmode,imode)*g_int(ibnd,jbnd,jmode)
              !
           enddo !jmode
        enddo !imode
        !
     enddo !jbnd
  enddo !ibnd
  !
  return

end subroutine G_kqw2kqb_no_spin



