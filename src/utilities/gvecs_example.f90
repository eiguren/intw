program proba

  use kinds, only: dp
  use intw_useful_constants, only: pi
  use intw_reading, only: bg, tpiba2

  implicit none

  ! integer, PARAMETER :: nr1=32, nr2=32, nr3=480
  ! real(kind=dp), PARAMETER :: ecut=80.0_dp, alat=5.51_dp
  ! logical :: gamma_only=.false.
  integer, PARAMETER :: nr1=180, nr2=180, nr3=648
  real(kind=dp), PARAMETER :: ecut=150.0_dp, alat=22.04_dp
  logical :: gamma_only=.true.

  real(kind=dp) :: ecutrho

  real(kind=dp) :: G(3)
  integer :: nri1, nrj2, nrk3
  integer :: ngm

  integer :: i, j, k

  tpiba2=(2.0_dp*pi/alat)**2
  ecutrho=4.0_dp*ecut


  ! bg(:,1)=(/ 1.0 , 0.0 , 0.0                      /)
  ! bg(:,2)=(/ 0.0 , 1.0 , 0.0                      /)
  ! bg(:,3)=(/ 0.0 , 0.0 , 0.0666666666666666666668 /)
  bg(:,1)=(/ 1.0 , 0.0 , 0.0                /)
  bg(:,2)=(/ 0.0 , 1.0 , 0.0                /)
  bg(:,3)=(/ 0.0 , 0.0 , 0.2666666666666667 /)

  nri1=int(nr1/2.0)
  nrj2=int(nr2/2.0)
  nrk3=int(nr3/2.0)


  ngm=0
  iloop: do  i=-nri1,nri1! contar cuantos gs caen en Ecut
    !
    IF ( gamma_only .and. i < 0) CYCLE iloop
    !
    jloop: do  j=-nrj2,nrj2
      !
      IF ( gamma_only .and. i == 0 .and. j < 0) CYCLE jloop
      !
      kloop: do  k=-nrk3,nrk3
        !
        IF ( gamma_only .and. i == 0 .and. j == 0 .and. k < 0) CYCLE kloop
        !
        G = i*bg(:,1) + j*bg(:,2) + k*bg(:,3)
        if ( G(1)**2 + G(2)**2 + G(3)**2 <= ecutrho/tpiba2 ) then
          ngm=ngm+1
        end if
        !
      enddo kloop
    enddo jloop
  enddo iloop

  print*, ngm
  if (gamma_only) print*, ngm*2-1


end program proba
