program ep_melements

  ! Module variables
  use kinds, only: dp

  ! use intw_useful_constants, only: cmplx_i, cmplx_0, cmplx_1, tpi
  use intw_useful_constants, only: cmplx_0, Ha_in_eV, zero

  use intw_input_parameters, only: mesh_dir, prefix, nk1, nk2, nk3, &
                                   ph_dir, dvscf_dir, dvscf_name, ep_mat_file, &
                                   calc_epmat, nq1, nq2, nq3, nqirr

  ! use intw_reading, only:  at, bg, alat, npol s, ftau, nsym

  use intw_reading, only: nspin, noncolin, lspinorb, ntyp, ngm, gvec, nG_max, &
                          nr1, nr2, nr3, nat, nbands, ityp, atom_labels, tau, at, alat, volume0

  use intw_pseudo, only: read_all_pseudo, average_pp

  use intw_fft, only: nl

  ! Module subroutines and functions
  use intw_reading, only: read_parameters_data_file_xml, get_gvec

  use intw_utility, only: find_free_unit, switch_indices_zyx
  !================================================================================
  !       Declare the variables
  !================================================================================

  implicit none

  ! Parameters
  integer, parameter :: stdout=6
  character(len=1), parameter :: CR=CHAR(13), LF=CHAR(10)


  ! Local variables

  ! FFT variables
  integer :: all_ngm
  integer, allocatable :: nlm(:)


  ! wfc related variables
  ! real(kind=dp) :: eig

  complex(kind=dp), allocatable :: wfc_kq_g(:) ! wave function k+q in G space (nG_max)
  complex(kind=dp), allocatable :: wfc_kq_r(:) ! wave function k+q in real space (nr1*nr2*nr3)
  complex(kind=dp), allocatable :: wfc_k_g(:)  ! wave function k in G space (nG_max)
  complex(kind=dp), allocatable :: wfc_k_r(:)  ! wave function k in real space (nr1*nr2*nr3)

  integer, allocatable :: list_iG(:) ! ZER DA HAU?

  ! dv related variables
  real(kind=dp), allocatable :: dv_local(:) ! local part of the variation of the potential in real space (nr1*nr2*nr3)

  ! ep_mat_el related variables
  real(kind=dp) :: ep_mat_el ! conjg(wfc_kq)*dv_q*wfc_k
  integer :: ep_unit

  ! I/O variables
  integer :: ierr, record_length, ios, ounit
  real(kind=dp), allocatable, dimension(:,:,:) :: rho


  ! Loop variables
  integer :: is    ! Loop on spin components (is=1,nspin)
  integer :: iband ! Loop on band index (iband=1,nbands)
  integer :: jband ! Loop on band index (jband=1,nbands)
  integer :: imode ! Loop on vibrational modes (imode=1,3*nat)

  ! Timing variables
  real(kind=dp) :: t0, t1, t2

  integer :: i, x, y, z
  real(kind=dp) :: e_fermi
  real(kind=dp), allocatable :: eig(:,:), eig_up(:,:), eig_down(:,:), buf(:)
  integer, allocatable :: energy_list(:,:)

  20 format(A)
  30 format(A,F8.2,6X,A)
  !
  !
  !================================================================================
  !       Talk to the user
  !================================================================================
  !
  write(stdout,20) '================================================================='
  write(stdout,20) '|                         program me                            |'
  write(stdout,20) '|---------------------------------------------------------------|'
  write(stdout,20) '================================================================='
  !
  !
  !================================================================================
  !       Specify input parameters
  !================================================================================
  !
  ! This is the nspin=1 calculation
  ! Then we will be interested with the nspin=2 calculation
  mesh_dir = '/scratch/hgarai/Kalkuluak/Fe-MgO-Ag100/1ML/PH/smearing_0.01/V0/'
  mesh_dir = '/home/haritz/Kalkuluak/Ag/Fe-MgO-Ag100/1ML/ph/Fe-MgO-Ag100/spin/smearing_0.01/'
  prefix = 'Fe-MgO-Ag100'
  !mesh_dir = '/home/haritz/Kalkuluak/Ag/Fe-MgO-Ag100/1ML/ph/Fe-MgO-Ag100/smearing_0.01/V0/'
  nk1 = 1
  nk2 = 1
  nk3 = 1
  calc_epmat=.true.
  ph_dir='./'
  dvscf_dir='./'
  dvscf_name='fe-mgo-ag100.dvscf_q1'
  nq1=1
  nq2=1
  nq3=1
  nqirr=1
  !
  !
  !================================================================================
  !       read the parameters from the SCF QE calculation
  !================================================================================
  !
  ! Read from data-file:
  ! -alat (Bohr)
  ! -at (in units of alat), a1=at(:,1)
  ! -volume0 (Bohr)
  ! -bg (2pi/alat)
  ! -ecutwfc (Hartree)
  ! -ecutrho (Hartree)
  ! -nr1, nr2, nr3 (FFT_GRID)
  ! -noncolin
  ! -lspinorb
  ! -spinorb_mag
  !   if (noncolin) then
  !     nspin = 2
  !   else
  !     nspin = 1
  !   end if
  ! -nat
  ! -ntyp
  ! allocate(atom_labels(ntyp))
  ! allocate(atom_pfile(ntyp))
  ! allocate(ityp(nat))
  ! allocate(tau(3,nat)) (alat)
  ! allocate(amass(ntyp)) (amu)
  ! -nsym
  ! allocate(s(3,3,nsym))
  ! allocate(ftau(3,nsym))
  ! allocate(can_use_TR(nsym))
  ! -nkpoints_QE
  ! -nG_max (MAX_NUMBER_OF_GK-VECTORS)
  ! -nbands
  !
  ! Remember that this subroutine reads nspin=2 only if noncolin=.true. !!!!!!!
  !
  call read_parameters_data_file_xml()
  !
  nspin=2
  !
  !================================================================================
  !       Some checks
  !================================================================================
  !
  if (nspin==1) then
     !
     write(stdout,20) '|      - The calculation is paramagnetic nspin=1                |'
     write(stdout,20) '|---------------------------------------------------------------|'
     !
  elseif (nspin==2) then
     !
     write(stdout,20) '|      - Spin calculation nspin = 2                             |'
     !
     if (noncolin) then
        !
        write(stdout,20) '|      - Non-collinear Spin calculation                         |'
        !
     endif
     !
     write(stdout,20) '|---------------------------------------------------------------|'
     !
  else
     !
     write(stdout,20) '*****************************************************************'
     write(stdout,20) '*          ERROR: Allowed values for nspin are 1 or 2           *'
     write(stdout,20) '*                 program stops.                                *'
     write(stdout,20) '*****************************************************************'
     !
     stop
     !
  endif
  !
  !
  !
  !================================================================================
  !       Read PP files
  !================================================================================
  !
  write(stdout,20) '|      - Reading pseudopotentials from UPF files                |'
  !
  ! call read_all_pseudo ()
  ! if (.not.lspinorb) call average_pp(ntyp)
  !
  write(stdout,20) '|                          PPs are OK                           |'
  write(stdout,20) '|---------------------------------------------------------------|'
  !
  !
  !================================================================================
  !       Allocate and read G vectors
  !================================================================================
  !
  write(stdout,20) '|      - Reading G vectors                                      |'
  !
  allocate(gvec(3,ngm)) ! it is deallocated in the bottom of this main code.
  !
  all_ngm = 2*ngm-1
  !
  call get_gvec() ! read them from gvectors.dat file
  !
  write(stdout,20) '|---------------------------------------------------------------|'
  !
  !
  !================================================================================
  !       Allocate and calculate nl
  !================================================================================
  !
  write(stdout,20) '|      - Calculating nl                                         |'
  !
  allocate (nl(ngm), nlm(ngm))
  call generate_nl_haritz()
  !
  write(stdout,20) '|---------------------------------------------------------------|'
  !
  !
  !================================================================================
  !       Allocate wfc related variables
  !================================================================================
  !
  write(stdout,20) '|      - Allocate wfc related variables                         |'
  !
  allocate( list_iG(nG_max) )
  allocate(eig(nbands,nspin), eig_up(nbands,2), eig_down(nbands,2), buf(1))
  allocate(energy_list(nbands,2))
  allocate( wfc_k_g(nG_max), wfc_kq_g(nG_max) )
  allocate( wfc_k_r(nr1*nr2*nr3), wfc_kq_r(nr1*nr2*nr3) )
  list_iG  = 0
  energy_list = 0
  eig      = ZERO
  wfc_k_g  = cmplx_0
  wfc_kq_g = cmplx_0
  wfc_k_r  = cmplx_0
  wfc_kq_r = cmplx_0
  !
  write(stdout,20) '|---------------------------------------------------------------|'
  !
  !
  !================================================================================
  !       Allocate dv related variables
  !================================================================================
  !
  write(stdout,20) '|      - Allocate dv related variables                          |'
  !
  allocate( dv_local(nr1*nr2*nr3) )
  dv_local=cmplx_0
  !
  write(stdout,20) '|---------------------------------------------------------------|'
  !


  call get_e_fermi()
  call get_GAMMA_eigenvalues()
  ! do i=1,nbands
  !   print*, i, eig(i,1)-e_fermi, eig(i,2)-e_fermi
  ! enddo
  !
  do iband=1,nbands
    eig(iband,1) = eig(iband,1) - e_fermi
    energy_list(iband,1) = iband
    eig(iband,2) = eig(iband,2) - e_fermi
    energy_list(iband,2) = iband
  enddo
  !
  ! Sort

  do iband = 1, nbands
    !
    ! Spin up
    jband = minloc( abs(eig( iband:nbands, 1 )), dim=1 ) + iband - 1
    !
    buf( 1 )     = energy_list( iband, 1 )
    energy_list( iband, 1 ) = energy_list( jband, 1 )
    energy_list( jband, 1 ) = buf( 1 )
    buf( 1 )     = eig( iband, 1 )
    eig( iband, 1 ) = eig( jband, 1 )
    eig( jband, 1 ) = buf( 1 )
    !
    ! Spin down
    jband = minloc( abs(eig( iband:nbands, 1 )), dim=1 ) + iband - 1
    !
    buf( 1 )     = energy_list( iband, 2 )
    energy_list( iband, 2 ) = energy_list( jband, 2 )
    energy_list( jband, 2 ) = buf( 1 )
    buf( 1 )     = eig( iband, 2 )
    eig( iband, 2 ) = eig( jband, 2 )
    eig( jband, 2 ) = buf( 1 )
    !
  enddo
  !
  ! do i=1,nbands
  !   print*, i, eig(i,1), energy_list(i,1)
  ! enddo
  !
  ! stop
  ! iband=47
  ! is=1
  ! call read_wfc( iband, is, wfc_kq_g )
  !
  ! call wfc_from_g_to_r( list_iG, wfc_kq_g, wfc_kq_r )
  !
  ! allocate(rho(nr1,nr2,nr3))
  !
  ! do i=1,nr1*nr2*nr3
  !   call switch_indices_zyx(nr1,nr2,nr3,i,x,y,z,-1)
  !   rho(x,y,z)=abs(wfc_kq_r(i))**2.0_dp/volume0*sign(-1.0_dp,real(wfc_kq_r(i)))
  ! enddo
  !
  ! ounit = find_free_unit()
  ! open(unit=ounit, file='fe.xsf', iostat=ios, status="replace", action="write")
  ! if ( ios /= 0 ) stop "Error opening file fe.xsf"
  ! !
  ! call xsf_struct(alat,at,nat,tau,atom_labels,ityp,ounit)
  ! !
  ! call xsf_fast_datagrid_3d(rho,nr1,nr2,nr3,nr1,nr2,nr3,at,alat,ounit)
  ! !
  ! close(unit=ounit, iostat=ios)
  ! if ( ios /= 0 ) stop "Error closing file unit ounit"
  !
  !
  ! stop

  !
  if (calc_epmat) then
    !
    ep_unit=find_free_unit()
    inquire(iolength=record_length) ep_mat_el
    open( unit=ep_unit, file=trim(trim(mesh_dir)//trim(ep_mat_file)), iostat=ierr, &
          form='unformatted', status='new', access='direct', recl=record_length )
    !
    !
    ! For each mode
    do imode=1,3*nat
      !
      ! Read from the dvscf file the induced part of the variation of the potential
      call read_dv(imode, dv_local)
      call RANDOM_NUMBER(dv_local)
      !
      ! For each state
      do iband=1,nbands
        do jband=1,nbands
          !
          ! For each spin component
          do is=1,nspin
            !
            ! Talk to the user
            write(unit=stdout,fmt="(a7,i3,a1,i3,x,a5,i7,a1,i7)", advance="no") &
              CR//"Mode: ", imode,"/",3*nat, "wfc: ", is*iband*jband,"/",(nspin*nbands**2)
            flush(stdout)
            write(unit=stdout,fmt="(a1)", advance="no") LF
            !
            call CPU_TIME(t0)
            !
            ! Read the wave functions in G space
            call read_wfc( iband, is, wfc_kq_g )
            call read_wfc( jband, is,  wfc_k_g )
            !
            ! Fourier transform the wfc to r space
            call wfc_from_g_to_r( list_iG, wfc_kq_g, wfc_kq_r )
            call wfc_from_g_to_r( list_iG, wfc_k_g, wfc_k_r )
            !
            ! Calculate the matrix element
            call CPU_TIME(t1)
            call matrix_element(real(wfc_kq_r),dv_local,real(wfc_k_r))
            !
            call CPU_TIME(t2)
            print*, t1-t0, t2-t1, t2-t0
            stop
            !
          enddo ! is
          !
        enddo ! jband
      enddo ! iband
      !
    enddo ! imode
    !
    close(ep_unit)
    !
  else
    ! Read them
  endif

contains

  subroutine matrix_element(bra_wfc, dv, ket_wfc)
    !
    !
    !
    implicit none
    !
    ! I/O variables
    real(kind=dp), intent(in) :: bra_wfc(:), ket_wfc(:)
    real(kind=dp), intent(in)    :: dv(:)
    !
    ! Local variables
    integer :: ir
    !
    20 format(A)
    !
    ! Check if the dimensions are equal
    if ( size(bra_wfc) == size(ket_wfc) ) then
      if ( size(bra_wfc) == size(dv) ) then
      else
        write(stdout,20) "matrix_element: The dimensions of dv and wfc's are not the same"
        stop
      endif
    else
      write(stdout,20) "matrix_element: The dimensions of the wfc's are not the same"
      stop
    endif
    !
    ! Do the integral
    ep_mat_el = ZERO
    do ir=1,size(bra_wfc)
      ep_mat_el = ep_mat_el + (bra_wfc(ir))*dv(ir)*(ket_wfc(ir))
    enddo
    !
    ep_mat_el = ep_mat_el/size(bra_wfc)
    !
  end subroutine matrix_element

  subroutine get_e_fermi()
  !
  ! This subroutine reads the Fermi energy from the .save
  !
  use iotk_module, only: iotk_open_read, iotk_close_read, iotk_scan_dat, iotk_scan_begin, iotk_scan_end
  !
  implicit none
  !
  character(256) :: datafile ! full path of the data-file.xml file
  integer :: io_unit
  !
  datafile  = trim(mesh_dir)//trim(prefix)//".save/"//trim("data-file.xml")
  !
  io_unit = find_free_unit()
  call iotk_open_read (io_unit,trim(datafile))

  call iotk_scan_begin (io_unit,"BAND_STRUCTURE_INFO")
   call iotk_scan_dat (io_unit,"FERMI_ENERGY",e_fermi)
  call iotk_scan_end (io_unit,"BAND_STRUCTURE_INFO")

  call iotk_close_read (io_unit)

  e_fermi = e_fermi*Ha_in_eV

end subroutine get_e_fermi


  subroutine generate_nl_haritz()
    !
    ! nl is an array that gives the indices into the FFT grid
    ! for a particular g vector. This subroutine calculates the nl array
    !
    implicit none
    !
    ! local variables
    integer  :: n1, n2, n3
    integer  :: ng
    integer  :: switch
    logical  :: assigned(nr1,nr2,nr3)
    !
    nl(:)            = 0
    assigned         = .false.
    !
    switch    =   1   ! triplet - to - singlet index
    !
    ! loop on all G vectors in the global array gvec
    do ng = 1, ngm
      !
      n1       = modulo(gvec(1,ng), nr1)+1 !modulo(gvec(1,ng), nr1)+1
      n2       = modulo(gvec(2,ng), nr2)+1
      n3       = modulo(gvec(3,ng), nr3)+1
      !
      if (.not. assigned(n1,n2,n3) ) then
        !
        assigned(n1,n2,n3)  =  .true.
        !
        ! compute the index corresponding to n1,n2,n3 and
        ! assign it to nl(ng)
        !
        call switch_indices_zyx(nr1,nr2,nr3,nl(ng),n1,n2,n3,switch)
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
      if ( gvec(1,ng)==-gvec(1,ng) .and. gvec(2,ng)==-gvec(2,ng) .and. gvec(3,ng)==-gvec(3,ng) ) cycle
      !
      n1       = modulo(-gvec(1,ng), nr1)+1 !modulo(gvec(1,ng), nr1)+1
      n2       = modulo(-gvec(2,ng), nr2)+1
      n3       = modulo(-gvec(3,ng), nr3)+1
      !
      if (.not. assigned(n1,n2,n3) ) then
        !
        assigned(n1,n2,n3)  =  .true.
        !
        ! compute the index corresponding to n1,n2,n3 and
        ! assign it to nl(ng)
        !
        call switch_indices_zyx(nr1,nr2,nr3,nlm(ng),n1,n2,n3,switch)
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
    end do
    !
  end subroutine generate_nl_haritz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_dv(mode, dv)
    !
    ! This subroutine reads the variation of the potential induced by
    ! the cartesian displacement # mode
    !
    implicit none
    !
    ! I/O variables
    integer, intent(in) :: mode
    real(kind=dp), intent(out) :: dv(:)
    !
    ! local variables
    integer :: read_unit, rl, ios
    character(len=256) :: dv_file
    !
    !
    ! Initialize the output variable
    dv = ZERO
    !
    ! Complete path to the file
    dv_file= trim(mesh_dir) // trim(ph_dir) // trim(dvscf_dir) // trim(dvscf_name)
    !
    ! Record length of the variable to read
    inquire(iolength=rl) dv
    !
    ! Open the file
    read_unit=find_free_unit()
    open( unit=read_unit, file=trim(dv_file), iostat=ios, form='unformatted', &
          status='old', action='read', access='direct', recl=rl )
    !
    ! Read the variation of the potential
    read(unit=read_unit,rec=mode,iostat=ios) dv
    !
    ! Close the file
    close(unit=read_unit)
    !
  end subroutine read_dv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_wfc(band, spin, wfc)
    !
    ! Reads the wfc's from QE .save
    !
    use iotk_module, only: iotk_open_read, iotk_close_read, iotk_scan_dat
    use intw_reading, only: write_tag
    !
    implicit none

    !I/O variables

    integer,intent(in) :: band, spin
    complex(kind=dp), intent(out) :: wfc(:)

    !logical variables

    character(256) :: K_directory, Kdir
    character(256) :: tag, gvec_file
    character(256), dimension(nspin) :: wfc_file
    integer :: io_unit
    integer :: nG
    !
    !
    ! Initialize the output variable
    wfc = cmplx_0
    !
    !
    ! Locate all the files we need
    100 format('K'I5.5'/')
    write(K_directory,100) 1
    Kdir=trim(mesh_dir)//trim(prefix)//".save/"//trim(K_directory)
    gvec_file=trim(trim(Kdir)//'gkvectors.dat')
    !
    if (nspin==1) then
       !
       wfc_file(1)=trim(Kdir)//'evc.dat'
       !
    else
       !
       wfc_file(1)=trim(Kdir)//'evc1.dat'
       wfc_file(2)=trim(Kdir)//'evc2.dat'
       !
    endif
    !
    !
    ! Read in relevant parameters (nG,list_iG)
    io_unit=find_free_unit()
    call iotk_open_read(io_unit,gvec_file)
    !
    call iotk_scan_dat (io_unit,"NUMBER_OF_GK-VECTORS",nG)
    call iotk_scan_dat (io_unit,"INDEX",list_iG(1:nG))
    !
    call iotk_close_read(io_unit)
    !
    !
    ! Read the wave function
    io_unit = find_free_unit()
    call iotk_open_read(io_unit,wfc_file(spin))
    !
    call write_tag("evc.",band,tag)
    call iotk_scan_dat (io_unit,tag,wfc(1:nG))
    !
    call iotk_close_read(io_unit)
    !
    return
    !
  end subroutine read_wfc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine wfc_from_g_to_r (list_iG,wfc_g, wfc_r)
    !--------------------------------------------------------
    !  This subroutine is a driver which uses the 3D-FFT
    !  code to transform a wavefunction in G space to
    !  a wavefunction in r space.
    !
    !     in ::        wfc_g(nG_max)         : the u_{nk}(G) coefficients
    !                  list_iG               : the indices of the G vectors used
    !                                          in the wave function
    !
    !     out::  wfc_r(nr1*nr2*nr3)          : The periodic part of the wave functions
    !                                          in real space, with the space index
    !                                          represented by a scalar.
    !--------------------------------------------------------
    implicit none

    integer       :: i,  iG
    integer       :: list_iG(nG_max)

    complex(dp)   :: wfc_g(nG_max)
    complex(dp)   :: wfc_r(nr1*nr2*nr3)



    ! initialize work array
    wfc_r(:)  =  cmplx_0

    ! put wfc_g in wfc_r
    do i=1,nG_max
        ! identify the G vector by its index, as stored in list_iG
        iG = list_iG(i)

        if (iG == 0) exit
        ! use nl to identify which G_fft vector G corresponds to,
        ! and assign the value of the wave function in the aux array
        wfc_r(nl(iG)) = wfc_g(iG)
    enddo
    do i=2,nG_max ! do not use iG=1, which is Gamma
        ! identify the G vector by its index, as stored in list_iG
        iG = list_iG(i)

        if (iG == 0) exit
        ! use nl to identify which G_fft vector G corresponds to,
        ! and assign the value of the wave function in the aux array
        wfc_r(nlm(iG)) = conjg(wfc_g(iG))
    enddo

    ! perform fourier transform in place wfc_g(G) -> wfc_r(r)
    ! CONVENTION BY ASIER
     call cfftnd(3,(/nr1,nr2,nr3/),1,wfc_r) !
                               ! this convention reproduces
                               ! the results of pw2wannier EXACTLY


  end subroutine wfc_from_g_to_r



  subroutine get_GAMMA_eigenvalues()
    !
    ! Reads the eigenvalues from QE .save
    !

    use iotk_module

    implicit none

    !local variables
    character(256) :: K_directory, Kdir
    character(256) :: tag
    character(256), dimension(nspin) :: eig_file
    integer :: io_unit,is,i,n_yes,ibnd
    integer :: nG
    !
    ! initialize the arrays to zero (zero will be broadcasted)
    !
    write(K_directory,100) 1
    !
    ! locate all the files we need for the next
    !
    Kdir=trim(mesh_dir)//trim(prefix)//".save/"//trim(K_directory)
    !
    !
    if (nspin==1) then
       !
       eig_file(1)=trim(Kdir)//'eigenval.xml'
       !
    else
       !
       eig_file(1)=trim(Kdir)//'eigenval1.xml'
       eig_file(2)=trim(Kdir)//'eigenval2.xml'
       !
    endif
    !
    ! read the eigenvalues for every band
    !
    write(stdout,*) "Reading eigenvalues"
    do is=1,nspin
       !
       io_unit = find_free_unit()
       !
       call iotk_open_read(io_unit,eig_file(is))
       !
       call iotk_scan_dat (io_unit,"EIGENVALUES",eig(:,is))
       !
       call iotk_close_read(io_unit)
       !
    enddo
    !
    ! convert eigenvalues from a.u. to eV
    !
    eig=eig*Ha_in_eV
    !
    100 format('K'I5.5'/')
    !
    return

  end subroutine get_GAMMA_eigenvalues

end program ep_melements
