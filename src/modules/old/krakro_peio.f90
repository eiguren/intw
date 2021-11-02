!----------------------------------------------------------------------------------
subroutine krakro_peio(e_min,e_max,npts,sigma,listin,listout)
!----------------------------------------------------------------------------------

  use intw_reading
  use intw_useful_constants
  use intw_W90
  use w90_parameters, only: num_wann,num_bands

  implicit none

  real(dp) :: e_min,e_max,batuketa,eps,epsp,sigma
  integer :: npts,i,j
  real(dp) :: listin(npts),listout(npts)

  listout(:)=0.0d0
  !
  do i=1,npts
     !
     eps=e_min+abs(e_max-e_min)*(i-1)/(npts-1)
     batuketa=0.0d0
     !
     do j=1,npts
        !
        epsp=e_min+abs(e_max-e_min)*(j-1)/(npts-1)
        !
        batuketa=batuketa+listin(j)*(eps-epsp)/((eps-epsp)**2.d0+sigma**2.d0)
        !
     enddo
     !
     listout(i)=batuketa/(npts*pi)
     !
  enddo

end subroutine krakro_peio
