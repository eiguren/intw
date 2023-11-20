module intw_haritz
  !
  ! This module contains general purpose variables and subroutines
  ! that I use in many programs.
  !
  use kinds, only: dp
  !
  implicit none
  !
  save
  !
  ! Declare public variables
  ! Conversion parameters
  public :: eV_to_Ha, eV_to_Ry, Ha_to_eV, Ha_to_Ry, Ry_to_eV, Ry_to_Ha
  ! I/O variables
  public :: stdout, stdin, CR, LF
  ! FFT variables
  !ASIER jULY 2021
   public ::  nlm, all_ngm

  !public :: gamma_only, nlm, all_ngm
  !
  ! Declare public subroutines
  public :: wfc_from_g_to_r_haritz, &
            wfc_from_r_to_g_haritz, &
            generate_nl_haritz, &
            get_gamma_only_haritz, &
            get_spin_haritz
  !
  private
  !
  ! Energy unit conversion parameters
  real(kind=dp), parameter :: eV_to_Ha = 1.0_dp/27.211383860484776_dp
  real(kind=dp), parameter :: eV_to_Ry = 2.0_dp/27.211383860484776_dp
  real(kind=dp), parameter :: Ha_to_eV = 27.211383860484776_dp
  real(kind=dp), parameter :: Ha_to_Ry = 2.0_dp
  real(kind=dp), parameter :: Ry_to_eV = 27.211383860484776_dp/2.0_dp
  real(kind=dp), parameter :: Ry_to_Ha = 0.5_dp
  !
  ! I/O variables
  integer         , parameter :: stdout=6    ! standard output unit
  integer         , parameter :: stdin=5     ! standard input unit
  character(len=1), parameter :: CR=CHAR(13) ! carriage return character
  character(len=1), parameter :: LF=CHAR(10) ! carriage return character
  !
  ! FFT variables
  !ASIER July 2021
  !logical :: gamma_only ! specifies if the QE calculation uses only Gamma point
  integer, allocatable :: nlm(:)
  integer :: all_ngm



contains


  subroutine wfc_from_g_to_r_haritz( list_iG, wfc_g, wfc_r )
    !
    !  This subroutine is a driver which uses the 3D-FFT
    !  code to transform a wavefunction in G space to
    !  a wavefunction in r space.
    !
    use intw_reading, only: nG_max, nr1, nr2, nr3
    use intw_useful_constants, only: cmplx_0
    use intw_fft, only: nl

    !ASier
    use intw_reading, only:  gamma_only

    !
    implicit none
    !
    ! I/O variables
    integer         , intent(in)  :: list_iG(nG_max)    ! the indices of the G vectors used in the wave function
    complex(kind=dp), intent(in)  :: wfc_g(nG_max)      ! the u_{nk}(G) coefficients
    complex(kind=dp), intent(out) :: wfc_r(nr1*nr2*nr3) ! The periodic part of the wave functions in real space
    !
    ! local variables
    integer :: i, iG
    !
    !
    ! initialize work array
    wfc_r(:)  =  cmplx_0
    !
    ! put wfc_g in wfc_r
    do i=1,nG_max
        ! identify the G vector by its index, as stored in list_iG
        iG = list_iG(i)
        if (iG == 0) exit
        ! use nl to identify which G_fft vector G corresponds to,
        ! and assign the value of the wave function in the aux array
        wfc_r(nl(iG)) = wfc_g(iG)
    enddo

    if (gamma_only) then
      do i = 1, nG_max
        iG = list_iG(i)
        if (iG == 0) exit
        wfc_r(nlm(ig)) = conjg(wfc_r(nl(ig)))
      enddo
    end if

    ! perform fourier transform in place wfc_g(G) -> wfc_r(r)
    ! CONVENTION BY ASIER
     call cfftnd(3,(/nr1,nr2,nr3/),1,wfc_r) !
                               ! this convention reproduces
                               ! the results of pw2wannier EXACTLY


  end subroutine wfc_from_g_to_r_haritz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wfc_from_r_to_g_haritz( list_iG, wfc_r, wfc_g )
    !
    !  This subroutine is a driver which uses the 3D-FFT
    !  code to transform a wavefunction in r space to
    !  a wavefunction in G space.
    !
    use intw_reading, only: nG_max, nr1, nr2, nr3
    use intw_useful_constants, only: cmplx_0
    use intw_fft, only: nl
    !
    implicit none
    !
    ! I/O variables
    integer         , intent(in)  :: list_iG(nG_max)    ! the indices of the G vectors used in the wave function
    complex(kind=dp), intent(in)  :: wfc_r(nr1*nr2*nr3) ! The periodic part of the wave functions in real space
    complex(kind=dp), intent(out) :: wfc_g(nG_max)      ! the u_{nk}(G) coefficients
    !
    ! local variables
    complex(kind=dp) :: wfc_r_tmp(nr1*nr2*nr3)
    integer :: i, iG
    !
    !
    ! initialize work arrays
    wfc_g(:) = cmplx_0
    wfc_r_tmp(:) = wfc_r(:)
    !
    ! perform fourier transform in place wfc_g(G) -> wfc_r(r)
    ! CONVENTION BY ASIER
     call cfftnd(3,(/nr1,nr2,nr3/),-1,wfc_r_tmp) !
                               ! this convention reproduces
                               ! the results of pw2wannier EXACTLY
    !
    ! put wfc_r_tmp in wfc_g
    do ig=1,ng_max
         wfc_g(ig)=wfc_r_tmp(nl(ig))
    enddo
    !
  end subroutine wfc_from_r_to_g_haritz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine generate_nl_haritz()
    !
    ! nl is an array that gives the indices into the FFT grid
    ! for a particular g vector. This subroutine calculates the nl array
    !
    use intw_useful_constants, only: tpi
    use intw_reading, only: ngm, gvec, nat, bg, nr1, nr2, nr3, tau
    use intw_utility, only: switch_indices_zyx
    use intw_fft, only: nl, gvec_cart, mill, eigts1, eigts2, eigts3

    use intw_reading, only: gamma_only
    !
    implicit none
    !
    ! local variables
    integer :: n1, n2, n3, ng, na, ipol, igm_, switch
    logical :: assigned(nr1,nr2,nr3)
    real(kind=dp) :: arg
    real(kind=dp), dimension(3) :: bgtau
    !
    assigned = .false.
    !
    nl = 0
    if (gamma_only) nlm = 0
    !
    switch    =   1   ! triplet - to - singlet index
    !
    ! loop on all G vectors in the global array gvec
    do igm_=1,ngm
      !
      n1 = modulo( gvec(1,igm_), nr1 ) + 1
      n2 = modulo( gvec(2,igm_), nr2 ) + 1
      n3 = modulo( gvec(3,igm_), nr3 ) + 1
      !
      if (.not. assigned(n1,n2,n3) ) then
        !
        assigned(n1,n2,n3)  =  .true.
        !
        ! compute the index corresponding to n1,n2,n3 and
        ! assign it to nl(ng)
        !
        call switch_indices_zyx(nr1,nr2,nr3,nl(igm_),n1,n2,n3,switch)
        !
      else
        !
        write(stdout,*) 'ERROR in generate_nl. FFT mesh too small?'
        write(stdout,*) '    More than one G-vector in the gvec array are being'
        write(stdout,*) '    assigned to the same FFT triplet (n1,n2,n3);      '
        write(stdout,*) '    this suggests that the FFT mesh (nr1,nr2,nr3) is  '
        write(stdout,*) '    too small.                                        '
        stop
        !
      endif
      !
      !
      if (gamma_only) then
        !
        if ( gvec(1,igm_)==-gvec(1,igm_) .and. &
             gvec(2,igm_)==-gvec(2,igm_) .and. &
             gvec(3,igm_)==-gvec(3,igm_) ) assigned(n1,n2,n3)=.false.
        !
        n1 = modulo( -gvec(1,igm_), nr1 ) + 1
        n2 = modulo( -gvec(2,igm_), nr2 ) + 1
        n3 = modulo( -gvec(3,igm_), nr3 ) + 1
        !
        if (.not. assigned(n1,n2,n3) ) then
          !
          assigned(n1,n2,n3)  =  .true.
          !
          ! compute the index corresponding to n1,n2,n3 and
          ! assign it to nl(ng)
          !
          call switch_indices_zyx(nr1,nr2,nr3,nlm(igm_),n1,n2,n3,switch)
          !
        else
          !
          write(stdout,*) 'ERROR in generate_nl. FFT mesh too small?'
          write(stdout,*) '    More than one G-vector in the gvec array are being'
          write(stdout,*) '    assigned to the same FFT triplet (n1,n2,n3);      '
          write(stdout,*) '    this suggests that the FFT mesh (nr1,nr2,nr3) is  '
          write(stdout,*) '    too small.                                        '
          stop
          !
        endif
        !
      endif
      !
    enddo
    !
    ! Obtain the G vectors in cartesian coordinates.
    allocate(gvec_cart(3,ngm))
    gvec_cart(1:3,1:ngm) = gvec(1:3,1:ngm)
    call cryst_to_cart (ngm, gvec_cart , bg, 1)
    !
    ! Generate the mill coefficients
    allocate( mill(3,ngm) )
    do ng = 1, ngm
      !
      mill (1, ng) = gvec(1,ng)
      mill (2, ng) = gvec(2,ng)
      mill (3, ng) = gvec(3,ng)
      !
    enddo
    !
    ! Calculate the eigts phases
    allocate (eigts1(-nr1:nr1, nat))
    allocate (eigts2(-nr2:nr2, nat))
    allocate (eigts3(-nr3:nr3, nat))
    do na = 1, nat
       do ipol = 1, 3
          bgtau (ipol) = bg (1, ipol) * tau (1, na) + &
                         bg (2, ipol) * tau (2, na) + &
                         bg (3, ipol) * tau (3, na)
       enddo
       do n1 = - nr1, nr1
          arg = tpi * n1 * bgtau (1)
          eigts1 (n1, na) = CMPLX(cos (arg), - sin (arg) ,kind=DP)
       enddo
       do n2 = - nr2, nr2
          arg = tpi * n2 * bgtau (2)
          eigts2 (n2, na) = CMPLX(cos (arg), - sin (arg) ,kind=DP)
       enddo
       do n3 = - nr3, nr3
          arg = tpi * n3 * bgtau (3)
          eigts3 (n3, na) = CMPLX(cos (arg), - sin (arg) ,kind=DP)
       enddo
    enddo
    !
  end subroutine generate_nl_haritz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_gamma_only_haritz()
    !
    ! This subroutine reads from the data-file if the calculation is gamma_only
    !
    use intw_input_parameters, only: mesh_dir, prefix
    use intw_utility, only: find_free_unit
    !use iotk_module
    !
    implicit none
    !
    ! Local variables
    integer :: io_unit
    character(len=256) :: datafile
    !
    !io_unit = find_free_unit()
    !
    !datafile  = trim(mesh_dir)//trim(prefix)//".save/"//trim("data-file.xml")
    !
    !call iotk_open_read( io_unit, datafile )
      !
    !  call iotk_scan_begin( io_unit, "PLANE_WAVES" )
        !
    !    call iotk_scan_dat( io_unit, "GAMMA_ONLY", gamma_only )
        !
    !  call iotk_scan_end( io_unit, "PLANE_WAVES" )
      !
    !call iotk_close_read( io_unit )
    !
  end subroutine get_gamma_only_haritz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_spin_haritz()
    !
    ! intw reads only NON-COLINEAR_CALCULATION from data-file.
    ! If the calculation is not noncolinear but is a colinear one
    ! this subroutine checks it to set the correct value to nspin
    !
    use intw_input_parameters, only: mesh_dir, prefix
    use intw_reading, only: nspin, noncolin
    use intw_utility, only: find_free_unit
    !use iotk_module

    !ASIER JUly 2021
    use intw_reading, only: lsda
    !
    implicit none
    !
    integer        :: io_unit  ! input/output file unit
    character(256) :: datafile ! full path of the data-file.xml file
    !logical        :: lsda
    !
    !
    if (nspin==2 .and. noncolin) then
      !
      ! the value of nspin is correct
      return
      !
    elseif (nspin==1 .and. (.not.noncolin) ) then
      !
      ! the value of nspin could be incorrect
      !
      !io_unit = find_free_unit()
      !datafile  = trim(mesh_dir)//trim(prefix)//".save/"//trim("data-file.xml")
      !call iotk_open_read(io_unit,datafile)
      !  call iotk_scan_begin (io_unit,"SPIN")
      !    call iotk_scan_dat (io_unit,"LSDA",lsda)
      !  call iotk_scan_end (io_unit,"SPIN")
      !call iotk_close_read(io_unit)
      !
      ! check the values
      if (lsda) then
        !
        ! nspin should be 2
        nspin=2
        return
        !
      else
        !
        ! the value of nspin is correct
        return
        !
      endif
      !
    else
      !
      ! this situation should not happen
      stop "ERROR: get_spin"
      !
    endif
    !
  end subroutine get_spin_haritz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




end module intw_haritz

module intw_mpi_haritz
  !
  ! This module contains MPI variables that I use in many MPI programs.
  !
  implicit none
  !
  save
  !
  ! Declare public variables
  public :: master, nNodes, myNode
  !
  private

  integer, parameter :: master = 0      ! ID of the master node
  integer            :: nNodes          ! Number of MPI nodes
  integer            :: myNode          ! ID of each MPI node

end module intw_mpi_haritz
