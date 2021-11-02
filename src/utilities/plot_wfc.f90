program wfc
  !
  ! This program reads the wfc (in G space) and transforms it to r space,
  ! writing it on xsf format, suitable for xcrysden or vesta.
  !
  ! Import modules
  !
  ! Module variables
  use kinds, only: dp
  !
  use intw_useful_constants, only: cmplx_0, Ha_in_eV
  !
  use intw_input_parameters, only: mesh_dir, prefix, nk1, nk2, nk3
  !
  use intw_reading, only: nspin, ngm, gvec, nG_max, nr1, nr2, nr3, &
                          nat, nbands, ityp, atom_labels, tau, at, &
                          alat, volume0, nkpoints_QE
  !
  use intw_fft, only: nl
  !
  use intw_haritz, only: stdout, stdin, CR, gamma_only, nlm, all_ngm
  !
  ! Module subroutines and functions
  use intw_reading, only: read_parameters_data_file_xml, get_gvec
  use intw_reading, only: write_tag
  !
  use intw_utility, only: find_free_unit, switch_indices_zyx
  !
  use iotk_module
  !
  use intw_haritz, only: get_spin_haritz, get_e_fermi_haritz, &
                         generate_nl_haritz, wfc_from_g_to_r_haritz, &
                         get_gamma_only_haritz
  !
  ! Declare variables
  !
  implicit none
  !
  ! wfc variables
  complex(kind=dp), allocatable :: wfc_g(:)   ! wave function in G space (nG_max)
  complex(kind=dp), allocatable :: wfc_r(:)   ! wave function in real space (nr1*nr2*nr3)
  integer         , allocatable :: list_iG(:) !
  !
  ! eigenvalue variables
  real(kind=dp), allocatable :: eig(:,:,:)
  real(kind=dp) :: e_fermi

  ! FFT variables
  integer :: nG_K

  ! input variables
  integer :: ikpoint
  integer :: ibands
  integer :: ispin
  integer :: nbands_initial
  integer :: nbands_final
  integer :: nspin_initial
  integer :: nspin_final

  ! loop variables
  integer :: ir

  ! variables to save the wavefunction in xsf format
  real(kind=dp), allocatable :: rho(:,:,:)
  integer                    :: x, y, z
  integer                    :: xsf_unit    !
  character(len=256)         :: xsf_file

  ! I/O variables
  integer                     :: wfc_unit    !
  integer                     :: ios         !
  !
  !
  20 format(A)
  !
  !================================================================================
  !       Talk to the user
  !================================================================================
  !
  write(stdout,20) '================================================================='
  write(stdout,20) '|                      program  plot_wfc.x                      |'
  write(stdout,20) '|              -----------------------------------              |'
  write(stdout,20) '================================================================='
  !
  !================================================================================
  !
  ! Read stdin file
  call read_input_file()
  !
  ! Print input parameters
  write(stdout,"(2a)") "mesh_dir: ", trim(mesh_dir)
  write(stdout,"(2a)") "prefix: ", trim(prefix)
  !
  !
  !================================================================================
  !       read the parameters from the SCF QE calculation
  !================================================================================
  !
  write(stdout,20) '|      - Reading QE parameters                                  |'
  !
  call read_parameters_data_file_xml()
  !
  if (nbands_initial==-1) nbands_initial=1
  if (nbands_final==-1) nbands_final=nbands
  !
  ! Check the limits
  if (nbands_initial<1) stop "ERROR: nbands_initial"
  if (nbands_final>nbands) stop "ERROR: nbands_final"
  if (nbands_final<nbands_initial) stop "ERROR: nbands_initial and nbands_final"
  !
  call get_gamma_only_haritz()
  !
  call get_spin_haritz()
  !
  if (nspin==1) then
    nspin_initial=1
    nspin_final=1
  elseif (nspin==2) then
    if (nspin_initial==-1) nspin_initial=1
    if (nspin_final==-1) nspin_final=2
  endif
  !
  ! Check the limits
  if (nspin>2) stop "ERROR: nspin"
  if (nspin_initial<1) stop "ERROR: nspin_initial"
  if (nspin_final>nspin) stop "ERROR: nspin_final"
  if (nspin_final<nspin_initial) stop "ERROR: nspin_initial and nspin_final"
  !
  !================================================================================
  !       Allocate and read G vectors
  !================================================================================
  !
  write(stdout,20) '|      - Reading G vectors                                      |'
  !
  allocate( gvec(3,ngm) )
  !
  call get_gvec()
  !
  !
  !================================================================================
  !       Allocate and calculate nl
  !================================================================================
  !
  write(stdout,20) '|      - Calculating nl                                         |'
  !
  allocate( nl(ngm) )
  if (gamma_only) allocate( nlm(ngm) )
  if (gamma_only) all_ngm = 2*ngm-1
  !
  call generate_nl_haritz()
  !
  !
  !================================================================================
  !       Allocate and read eigenvalues
  !================================================================================
  !
  write(stdout,20) '|      - Reading eigenvalues                                    |'
  !
  allocate( eig(nbands,nspin,nkpoints_QE) )
  !
  call get_eigenvalues()
  !
  call get_e_fermi_haritz(e_fermi)
  !
  !
  !================================================================================
  !       Allocate and read wavefunctions
  !================================================================================
  !
  write(stdout,20) '|      - Reading wavefunctions                                  |'
  !
  allocate( wfc_g(nG_max) )
  allocate( list_iG(nG_max) )
  allocate( wfc_r(nr1*nr2*nr3) )
  allocate(rho(nr1,nr2,nr3))
  !
  xsf_unit = find_free_unit()
  !
  write(stdout,"(4a14)") 'kpoint', 'spin', 'band', 'E-Ef'
  !
  do ikpoint=1,nkpoints_QE
    !
    do ispin=nspin_initial,nspin_final
      !
      do ibands=nbands_initial,nbands_final
        !
        write(stdout,"(3i14,f14.6)") ikpoint, ispin, ibands, eig(ibands,ispin,ikpoint)-e_fermi
        !
        wfc_unit=find_free_unit()
        call get_K_wfc(wfc_unit,ibands)
        !
        call wfc_from_g_to_r_haritz(list_iG,wfc_g,wfc_r)
        !
        do ir=1,nr1*nr2*nr3
          !
          call switch_indices_zyx(nr1,nr2,nr3,ir,x,y,z,-1)
          rho(x,y,z)=abs(wfc_r(ir))**2.0_dp/volume0*sign(-1.0_dp,real(wfc_r(ir)))
          !
        enddo
        !
        write(xsf_file,200) ikpoint, ispin, ibands
        200 format("wfc_k"i2.2"_s"i1"_b"i5.5".xsf")
        !
        open(unit=xsf_unit, file=xsf_file, iostat=ios, status="replace", action="write")
        if ( ios /= 0 ) stop "Error opening xsf_unit file"
        !
        call xsf_struct(alat,at,nat,tau,atom_labels,ityp,xsf_unit)
        call xsf_fast_datagrid_3d(rho,nr1,nr2,nr3,nr1,nr2,nr3,at,alat,xsf_unit)
        !
        close(unit=xsf_unit, iostat=ios)
        if ( ios /= 0 ) stop "Error closing xsf_unit file"
        !
      enddo ! ibands
      !
    enddo ! ispin
    !
  enddo ! ikpoint





contains

  subroutine get_eigenvalues()
    !
    ! Reads the eigenvalues from QE .save
    !
    !
    use iotk_module
    !
    implicit none
    !
    ! local variables
    character(len=256) :: K_directory, Kdir
    character(len=256) :: eig_file(nspin)
    integer :: io_unit, ik, is
    !
    !
    do ik=1,nkpoints_QE
      !
      ! directory corresponding to the K point
      write(K_directory,100) ik
      Kdir=trim(mesh_dir)//trim(prefix)//".save/"//trim(K_directory)
      !
      ! files we need to read
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
      ! read the eigenvalues for each spin component
      do is=1,nspin
         !
         io_unit = find_free_unit()
         !
         call iotk_open_read(io_unit,eig_file(is))
         !
         call iotk_scan_dat (io_unit,"EIGENVALUES",eig(:,is,ik))
         !
         call iotk_close_read(io_unit)
         !
      enddo
      !
    enddo
    !
    ! convert eigenvalues from a.u. to eV
    eig=eig*Ha_in_eV
    !
    100 format('K'I5.5'/')
    !
    return

  end subroutine get_eigenvalues


  subroutine get_K_wfc(unit,band)
    !
    ! Reads the wfc's from QE .save
    !
    use intw_reading, only: write_tag
    use iotk_module
    !
    implicit none
    !
    ! I/O variables
    integer,intent(in) :: unit
    integer,intent(in) :: band
    !
    ! local variables
    character(256) :: tag
    !
    !
    call open_K_wfc_file(unit,ikpoint,ispin)
    ! read the wavefunction
    call write_tag("evc.",band,tag)
    call iotk_scan_dat (unit,tag,wfc_g(1:nG_K))
    wfc_g(nG_K+1:nG_max) = cmplx_0
    call close_K_wfc_file(unit)
    !
  end subroutine get_K_wfc

  subroutine open_K_wfc_file(unit,i_kpoint,spin)
    !
    ! Reads the wfc's from QE .save
    !
    use intw_reading, only: write_tag
    use iotk_module
    !
    implicit none
    !
    ! I/O variables
    integer,intent(in) :: unit
    integer,intent(in) :: i_kpoint
    integer,intent(in) :: spin
    !
    ! local variables
    character(256) :: K_directory, Kdir
    character(256) :: wfc_file, gvec_file
    integer :: io_unit
    !
    !
    ! initialize the arrays to zero (zero will be broadcasted)
    !
    write(K_directory,100) i_kpoint
    !
    ! locate all the files we need for the next
    !
    Kdir=trim(mesh_dir)//trim(prefix)//".save/"//trim(K_directory)
    gvec_file=trim(trim(Kdir)//'gkvectors.dat')
    !
    if (nspin==1) then
       !
       wfc_file=trim(Kdir)//'evc.dat'
       !
    else
       !
       if (spin==1) wfc_file=trim(Kdir)//'evc1.dat'
       if (spin==2) wfc_file=trim(Kdir)//'evc2.dat'
       !
    endif
    !
    !
    ! read in relevant parameters (nG,list_iG) for this K point
    io_unit=find_free_unit()
    call iotk_open_read(io_unit,gvec_file)
      call iotk_scan_dat (io_unit,"NUMBER_OF_GK-VECTORS",nG_K)
      call iotk_scan_dat (io_unit,"INDEX",list_iG(1:nG_K))
      list_iG(nG_K+1:nG_max) = 0
    call iotk_close_read(io_unit)
    !
    ! open the wfc file
    call iotk_open_read(unit,wfc_file)
    !
    100 format('K'I5.5'/')
    !
  end subroutine open_K_wfc_file

  subroutine close_K_wfc_file(unit)
    !
    ! Reads the wfc's from QE .save
    !
    use iotk_module
    !
    implicit none
    !
    ! I/O variables
    integer,intent(in) :: unit
    !
    !
    call iotk_close_read(unit)
    !
  end subroutine close_K_wfc_file

  subroutine read_input_file()
    !
    ! This subroutine sets default values for variables of input namelist
    ! and then reads them from input file if there are specified
    !
    implicit none
    !
    character(len=256) :: input_file_name
    ! Namelists
    !
    ! input variables namelist
    NAMELIST / input_plot_wfc / mesh_dir, prefix, nbands_initial, nbands_final, nspin_initial, nspin_final
    !
    ! Set default values
    mesh_dir='./'
    prefix='Fe-MgO-Ag100'
    nbands_initial=-1
    nbands_final=-1
    nspin_initial=-1
    nspin_final=-1
    !
    INQUIRE(stdin, NAME=input_file_name)
    if (input_file_name(1:4)=="/dev") then
      ! there is no input file: use default values
      return
    else
      read(stdin,input_plot_wfc)
    endif
    !
    ! Reopen terminal as input file
    close(unit=stdin)
    open(unit=stdin,file='/dev/tty')
    !
  end subroutine read_input_file


end program wfc
