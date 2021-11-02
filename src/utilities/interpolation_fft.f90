program fourier_interpolation
  !
  ! This program reads a change of potential of a gamma_only supercell (real dv)
  ! in a coarse (nr1,nr2,nr3) FFT grid, and interpolates it using Fourier
  ! interpolation to a (nr1s,nr2s,nr3s) finer grid.
  !
  use kinds, only: dp
  !
  use intw_reading, only: bg
  !
  use intw_fft_interp
  !
  implicit none
  !
  !
  real(kind=dp), allocatable :: dv_r(:)
  real(kind=dp), allocatable :: dv_rs(:)
  complex(kind=dp), allocatable :: dv_aux(:)
  integer :: rl_dv_r, rl_dv_rs
  !
  ! Loop variables
  integer :: ia, id, imode
  !
  ! I/O variables
  integer :: ios
  !
  ! Input variables
  character(len=256) :: coarse_dv_file ! name of the coarse dv file (to read)
  character(len=256) :: fine_dv_file   ! name of the fine dv file (to write)
  logical :: coarse_dv_complex ! specify if the coarse dv is complex
  integer :: nr1=10, nr2=10, nr3=10    ! original coarse FFT grid dimensions
  integer :: nr1s=15, nr2s=15, nr3s=15 ! dimensions of the fine FFT grid to interpolate
  integer :: nat ! number of atoms

  !
  ! Specifying bg is necessary only if you use a cutoff
  bg(:,1)=(/1.0,0.0,0.0/)
  bg(:,2)=(/0.0,1.0,0.0/)
  bg(:,3)=(/0.0,0.0,1.0/)
  !
  !
  call read_input_file()
  !
  !
  allocate( dv_r( nr1* nr2* nr3) )
  inquire(iolength=rl_dv_r) dv_r(:)
  if ( coarse_dv_complex ) then
    allocate( dv_aux( nr1* nr2* nr3) )
    inquire(iolength=rl_dv_r) dv_aux(:)
  endif
  !
  allocate( dv_rs(nr1s*nr2s*nr3s) )
  inquire(iolength=rl_dv_rs) dv_rs(:)
  !
  print*, rl_dv_r, rl_dv_rs
  !
  ! Open the file to read the dv
  open ( unit = 1234, file = trim(coarse_dv_file), iostat = ios, &
         form = 'unformatted', status = 'old', access = 'direct', &
         action = 'read', recl = rl_dv_r )
  !
  ! Open the file where the dv will be written
  open ( unit = 4321, file = trim(fine_dv_file), iostat = ios, &
         form = 'unformatted', status = 'new', access = 'direct', &
         action = 'write', recl = rl_dv_rs )
  !
  !
  do ia=1,nat
    do id=1,3
      !
      imode = id + 3*(ia-1)
      !
      print*, imode, "/", 3*nat
      !
      ! read dv on the coarse grid
      if (coarse_dv_complex) then
        read(unit=1234, rec=imode, iostat=ios) dv_aux
        dv_r = real(dv_aux)
      else
        read(unit=1234, rec=imode, iostat=ios) dv_r
      endif
      if (ios/=0) print*, "read ios:", ios
      !
      ! interpolate to the fine grid using Fourier interpolation
      call fft_interp_3d_real(nr1, nr2, nr3, nr1s, nr2s, nr3s, huge(1.0_dp), dv_r, dv_rs)
      !
      !
      write(unit=4321, rec=imode, iostat=ios) dv_rs(:)
      if (ios/=0) print*, "write ios:", ios
      !
    enddo ! id
  enddo ! ia
  !
  deallocate(dv_r, dv_rs)
  if (coarse_dv_complex) deallocate(dv_aux)
  !
  close(1234)
  close(4321)
  !
  !



contains

  subroutine read_input_file()
    !
    ! This subroutine sets default values for variables of input namelist
    ! and then reads them from input file if there are specified
    !
    use intw_haritz, only: stdin
    !
    implicit none
    !
    character(len=256) :: input_file_name
    ! Namelists
    !
    ! input variables namelist
    NAMELIST / input_interpolation_fft / coarse_dv_file, fine_dv_file, coarse_dv_complex, &
                                         nr1, nr2, nr3, nr1s, nr2s, nr3s, nat
    !
    ! Set default values
    coarse_dv_file = "/scratch/hgarai/Kalkuluak/Fe-MgO-Ag100/1ML/M_ELEMENTS/mgo-ag100_R441.dvscf_q1"
    fine_dv_file = "/scratch/hgarai/Kalkuluak/Fe-MgO-Ag100/1ML/RELAX/Fe1/ph/mgo-ag100_R441_fine.dvscf_q1"
    coarse_dv_complex = .true.
    nr1=128
    nr2=128
    nr3=480
    nr1s=180
    nr2s=180
    nr3s=648
    nat=272
    !
    INQUIRE(stdin, NAME=input_file_name)
    if (input_file_name(1:4)=="/dev") then
      ! there is no input file: use default values
      return
    else
      read(stdin,input_interpolation_fft)
    endif
    !
    ! Reopen terminal as input file
    ! close(unit=stdin)
    ! open(unit=stdin,file='/dev/tty')
    !
    !
  end subroutine read_input_file



end program fourier_interpolation
