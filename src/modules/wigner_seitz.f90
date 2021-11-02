!--------------------------------------------------------------------------------------------------
subroutine count_wigner_seitz_points(nk1,nk2,nk3,nrws)
!--------------------------------------------------------------------------------------------------
!
! Adapted from Wannier90 by Peio G. Goiricelaya 07/03/2016
!
!================================================================================!
! Calculates a grid of points that fall inside of (and eventually on the         !
! surface of) the Wigner-Seitz supercell centered on the origin of the B         !
! lattice with primitive translations nmonkh(1)*a_1+nmonkh(2)*a_2+nmonkh(3)*a_3  !
!================================================================================!

  use kinds, only: dp
  use w90_parameters, only: real_metric
  use w90_constants, only: eps7,eps8

  implicit none

  !I/O variables

  integer,intent(in) :: nk1,nk2,nk3
  integer,intent(inout) :: nrws

  !local variables
  
  integer       :: ndiff (3)
  real(kind=dp) :: dist(125),dist_min
  integer       :: n1,n2,n3,i1,i2,i3,icnt,i,j

  ! The Wannier functions live in a supercell of the real space unit cell
  ! this supercell is nk unit cells long in each direction
  !
  ! We loop over grid points r on a unit cell that is 8 times larger than this
  ! primitive supercell. 
  !
  ! One of these points is in the W-S cell if it is closer to R=0 than any of the
  ! other points, R (where R are the translation vectors of the supercell)

  ! In the end nrpts contains the total number of grid
  ! points that have been found in the Wigner-Seitz cell

  nrws=0  
  do n1=-nk1,nk1
     do n2=-nk2,nk2
        do n3=-nk3,nk3
           ! Loop over the 125 points R. R=0 corresponds to i1=i2=i3=1, or icnt=14
           icnt=0  
           do i1=-2,2  
              do i2=-2,2  
                 do i3=-2,2  
                    icnt=icnt+1  
                    ! Calculate distance squared |r-R|^2
                    ndiff(1)=n1-i1*nk1
                    ndiff(2)=n2-i2*nk2
                    ndiff(3)=n3-i3*nk3
                    dist(icnt)=0.0d0
                    do i=1,3  
                       do j=1,3  
                          dist(icnt)=dist(icnt)+real(ndiff(i),dp)*real_metric(i,j) &
                               *real(ndiff(j),dp)
                       enddo !j
                    enddo !i
                 enddo !i3
              enddo !i2
           enddo !i1
           !
           dist_min=minval(dist)
           if (abs(dist(63)-dist_min).lt.eps7) then
              nrws=nrws+1  
           end if
        enddo !n3
     enddo !n2
  enddo !n1
  !
  return  

end subroutine
!--------------------------------------------------------------------------------------------------
!**************************************************************************************************
!--------------------------------------------------------------------------------------------------
subroutine wigner_seitz_points(nk1,nk2,nk3,nrws,irws,degws)
!--------------------------------------------------------------------------------------------------
!
! Adapted from Wannier90 by Peio G. Goiricelaya 07/03/2016
!
!================================================================================!
! Calculates a grid of points that fall inside of (and eventually on the         !
! surface of) the Wigner-Seitz supercell centered on the origin of the B         !
! lattice with primitive translations nmonkh(1)*a_1+nmonkh(2)*a_2+nmonkh(3)*a_3  !
!================================================================================!

  use kinds, only: dp
  use w90_parameters, only: real_metric
  use w90_constants, only: eps7,eps8

  implicit none

  !I/O variables

  integer,intent(in) :: nk1,nk2,nk3
  integer,intent(in) :: nrws
  integer,intent(inout) :: irws(3,nrws)
  integer,intent(inout) :: degws(nrws)

  !local variables

  integer       :: ndiff (3),nrws2
  real(kind=dp) :: dist(125),tot,dist_min
  integer       :: n1,n2,n3,i1,i2,i3,icnt,i,j

  ! The Wannier functions live in a supercell of the real space unit cell
  ! this supercell is nk unit cells long in each direction
  !
  ! We loop over grid points r on a unit cell that is 8 times larger than this
  ! primitive supercell. 
  !
  ! One of these points is in the W-S cell if it is closer to R=0 than any of the
  ! other points, R (where R are the translation vectors of the supercell)

  ! In the end nrpts contains the total number of grid
  ! points that have been found in the Wigner-Seitz cell

  nrws2=0  
  do n1=-nk1,nk1
     do n2=-nk2,nk2
        do n3=-nk3,nk3
           ! Loop over the 125 points R. R=0 corresponds to i1=i2=i3=1, or icnt=14
           icnt=0  
           do i1=-2,2  
              do i2=-2,2  
                 do i3=-2,2  
                    icnt=icnt+1  
                    ! Calculate distance squared |r-R|^2
                    ndiff(1)=n1-i1*nk1
                    ndiff(2)=n2-i2*nk2
                    ndiff(3)=n3-i3*nk3
                    dist(icnt)=0.0d0  
                    do i=1,3  
                       do j=1,3  
                          dist(icnt)=dist(icnt)+real(ndiff(i),dp)*real_metric(i,j) &
                               *real(ndiff(j),dp)
                       enddo !j
                    enddo !i
                 enddo !i3
              enddo !i2
           enddo !i1
           !
           dist_min=minval(dist)
           if (abs(dist(63)-dist_min).lt.eps7 ) then
              nrws2=nrws2+1  
              degws(nrws2)=0
              do i=1,125
                 if (abs(dist(i)-dist_min).lt.eps7) degws(nrws2)=degws(nrws2)+1
              end do
              irws(1,nrws2)=n1  
              irws(2,nrws2)=n2   
              irws(3,nrws2)=n3   
           end if
        enddo !nr3
     enddo !nr2
  enddo !nr1
  !
  if (nrws2.ne.nrws) then
     write(*,*)'ERROR: nrws is different to nrws2 when being equal expected.'
  endif
  !
  ! Check the "sum rule"
  tot=0.0d0  
  do i=1,nrws2
     tot=tot+1.0_dp/real(degws(i),dp)  
  enddo
  if (abs(tot-real(nk1*nk2*nk3,dp))>eps8) then
     write(*,*)'ERROR in wigner_seitz: error in finding Wigner-Seitz points'
  endif
  !
  return  

end subroutine
!--------------------------------------------------------------------------------------------------
