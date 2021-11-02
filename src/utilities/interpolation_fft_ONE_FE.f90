program proba
  !
  use kinds, only: dp
  use intw_reading, only: bg!, nr1, nr2, nr3, &
!                          read_parameters_data_file_xml
  use intw_fft_interp
  use intw_utility, only: switch_indices_zyx
  use intw_useful_constants, only: pi
  !
  implicit none
  !
  integer :: nr1=10, nr2=10, nr3=10, &
             nr1s=15, nr2s=15, nr3s=15
  !
  complex(kind=dp), allocatable :: dv_sc_r(:)
  real(kind=dp), allocatable :: dv_sc_rs(:)
  real(kind=dp) :: x, y, z
  integer :: rl_dv_sc_r, rl_dv_sc_rs
  integer :: ia, id, imode, nat, ia_hand, imode_hand
  integer :: ios


  bg(:,1)=(/1.0,0.0,0.0/)
  bg(:,2)=(/0.0,1.0,0.0/)
  bg(:,3)=(/0.0,0.0,1.0/)


  nr1=128
  nr2=128
  nr3=480
  nr1s=180
  nr2s=180
  nr3s=648
  nat=273
  allocate( dv_sc_r ( nr1* nr2* nr3) )
  allocate( dv_sc_rs(nr1s*nr2s*nr3s) )
  !
  inquire(iolength=rl_dv_sc_r) dv_sc_r(:)
  inquire(iolength=rl_dv_sc_rs) dv_sc_rs(:)
  print*, rl_dv_sc_r, rl_dv_sc_rs
  !
  ! Open the files to read the dv
  open ( unit = 1234, file = '/scratch/hgarai/Kalkuluak/Fe-MgO-Ag100/1ML/M_ELEMENTS/mgo-ag100_R441.dvscf_q1', iostat = ios, &
         form = 'unformatted', status = 'old', access = 'direct', &
         action = 'read', recl = rl_dv_sc_r )
  open ( unit = 1235, file = '/scratch/hgarai/Kalkuluak/Fe-MgO-Ag100/1ML/M_ELEMENTS/fe-mgo-ag100.dvscf_hand_only', iostat = ios, &
         form = 'unformatted', status = 'old', access = 'direct', &
         action = 'read', recl = rl_dv_sc_rs )

  !
  ! Open the file where the dv will be written
  open ( unit = 4321, file = '/scratch/hgarai/Kalkuluak/Fe-MgO-Ag100/1ML/ONE_FE/SPIN/PH/fe-mgo-ag100.dvscf_q1', iostat = ios, &
         form = 'unformatted', status = 'new', access = 'direct', &
         action = 'write', recl = rl_dv_sc_rs )


  ! The dv calculated by hand are ordered in this way:
  ! 1.: atm = 273 ! Fe_1
  ! 3.: atm = 172 ! O_1
  ! 4.: atm = 186 ! O_2
  do ia=1,nat
    do id=1,3
      !
      imode = id + 3*(ia-1)
      !
      print*, imode, "/", 3 + 3*(nat-1)
      !
      if (ia==172) then
        ! O_1
        ia_hand=3
        imode_hand = id + 3*(ia_hand-1)
        !
        ! read the dv (this is already on the fine grid)
        read(unit=1235, rec=imode_hand, iostat=ios) dv_sc_rs
        !
      elseif (ia==186) then
        ! O_2
        ia_hand=4
        imode_hand = id + 3*(ia_hand-1)
        !
        ! read the dv (this is already on the fine grid)
        read(unit=1235, rec=imode_hand, iostat=ios) dv_sc_rs
        !
      elseif (ia==273) then
        ! Fe_1
        ia_hand=1
        imode_hand = id + 3*(ia_hand-1)
        !
        ! read the dv (this is already on the fine grid)
        read(unit=1235, rec=imode_hand, iostat=ios) dv_sc_rs
        !
      else
        !
        ! read dv on the coarse grid
        read(unit=1234, rec=imode, iostat=ios) dv_sc_r
        !
        ! interpolate to the fine grid using Fourier interpolation
        call fft_interp_3d_real(nr1, nr2, nr3, nr1s, nr2s, nr3s, huge(1.0_dp), real(dv_sc_r), dv_sc_rs)
        !
      endif
      !
      write(unit=4321, rec=imode, iostat=ios) dv_sc_rs(:)
      !
    enddo
    !
  enddo
  !
  deallocate(dv_sc_r, dv_sc_rs)
  !
  close(1234)
  close(1235)
  close(4321)
  !
  !



contains

  function f (xx,yy,zz)

    real(dp), intent(in) :: xx,yy,zz
    real(dp)             :: f

    f = -sin(2*pi*xx)*exp(-xx**2)*sin(2*pi*yy)*exp(-2*yy**2)*sin(4*pi*zz)*exp(-zz**2)

  end function f


end program proba
