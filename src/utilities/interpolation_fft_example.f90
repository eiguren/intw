program proba
  !
  use kinds, only: dp
  use intw_reading, only: bg, tpiba2
  use intw_fft_interp
  use intw_utility, only: switch_indices_zyx
  use intw_useful_constants, only: pi

  !
  implicit none
  !
  integer, parameter :: nr1=10, nr2=10, nr3=10, &
                        nr1s=15, nr2s=15, nr3s=15
  !
  real(kind=dp) :: fr(nr1*nr2*nr3), frs(nr1s*nr2s*nr3s)
  real(kind=dp) :: x, y, z
  integer :: i, j, k, ir



  bg(:,1)=(/1.0,0.0,0.0/)
  bg(:,2)=(/0.0,1.0,0.0/)
  bg(:,3)=(/0.0,0.0,1.0/)

  tpiba2 = (2*pi/1.0)**2


  ! open(unit=123,file="fr.dat",status="replace")
  ! write(123,*) '#   x         y         z              f(x)'
  do k=1,nr3
    do j=1,nr2
      do i=1,nr1
        !
        x= real((i-1),dp)/nr1
        y= real((j-1),dp)/nr2
        z= real((k-1),dp)/nr3
        call switch_indices_zyx(nr1,nr2,nr3,ir,i,j,k,+1)
        !
        ! print"(3f10.8,i5)", x,y,z,ir
        fr(ir)=f(x,y,z)
        !
        ! write(123,"(i10,x,f15.10)") ir, fr(ir) ! funtzioa bera
        !
      enddo
    enddo
  enddo

  do k=1,nr3
    do j=1,nr2
      do i=1,nr1
        !
        x= real((i-1),dp)/nr1
        y= real((j-1),dp)/nr2
        z= real((k-1),dp)/nr3
        call switch_indices_zyx(nr1,nr2,nr3,ir,i,j,k,+1)
        !
        fr(ir)=f(x,y,z)
        !
      enddo
    enddo
  enddo


  open(unit=123,file="fr.dat",status="replace")
  write(123,*) '#   ind              f(x)'
  do k=1,nr3
    do j=1,nr2
      do i=1,nr1
        !
        x= real((i-1),dp)/nr1
        y= real((j-1),dp)/nr2
        z= real((k-1),dp)/nr3
        call switch_indices_zyx(nr1,nr2,nr3,ir,i,j,k,+1)
        !
        write(123,"(i10,x,f15.10)") ir, f(x,y,z) ! funtzioa bera
        !
      enddo
    enddo
  enddo
  close(123)
  !
  open(unit=123,file="frs.dat",status="replace")
  write(123,*) '#   ind              f(x)'
  do k=1,nr3s
    do j=1,nr2s
      do i=1,nr1s
        !
        x= real((i-1),dp)/nr1s
        y= real((j-1),dp)/nr2s
        z= real((k-1),dp)/nr3s
        call switch_indices_zyx(nr1s,nr2s,nr3s,ir,i,j,k,+1)
        !
        write(123,"(i10,x,f15.10)") ir, f(x,y,z) ! funtzioa bera sare finean
        !
      enddo
    enddo
  enddo
  close(123)

  !
  call fft_interp_3d_real(nr1, nr2, nr3, nr1s, nr2s, nr3s, huge(1.0_dp), fr, frs)
  !

  open(unit=123,file="FFTfgs.dat",status="replace")
  write(123,*) '#   x         y         z              f(x)'
  do k=1,nr3s
    do j=1,nr2s
      do i=1,nr1s
        !
        call switch_indices_zyx(nr1s,nr2s,nr3s,ir,i,j,k,+1)
        !
        write(123,"(i10,x,2f15.10)") ir, frs(ir) ! funtzioa interpolatuta sare finean
        !
      enddo
    enddo
  enddo
  close(123)


contains

  function f (xx,yy,zz)

    real(dp), intent(in) :: xx,yy,zz
    real(dp)             :: f

    f = -sin(2*pi*xx)*exp(-xx**2)*sin(2*pi*yy)*exp(-2*yy**2)*sin(4*pi*zz)*exp(-zz**2)

  end function f


end program proba
