program wfc
  !
  ! This program reads the wfc (in G space) and transforms it to r space
  ! to calculate it's weight in the different parts of the unit cell
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
                          nbands, tau, at, alat, nkpoints_QE
  !
  use intw_fft, only: nl
  !
  use intw_haritz, only: stdout, stdin, CR, gamma_only, nlm, all_ngm
  !
  ! Module subroutines and functions
  use intw_reading, only: read_parameters_data_file_xml, get_gvec, &
                          write_tag
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
  !
  ! FFT variables
  integer :: nG_K
  !
  ! loop variables
  integer :: ikpoint
  integer :: ibands
  integer :: ispin
  integer :: ir
  !
  ! weights
  real(kind=dp) :: weight_iron_1
  real(kind=dp) :: weight_iron_2
  real(kind=dp) :: weight_mgo_1
  real(kind=dp) :: weight_mgo_2
  real(kind=dp) :: weight_silver
  real(kind=dp), allocatable :: weight_z(:)
  ! weight ranges
  ! real(kind=dp), parameter :: Fe1 = 32.098 ! Position of the superior Fe atom in bohr
  ! real(kind=dp), parameter :: Fe2 = 50.552 ! Position of the inferior Fe atom in bohr
  ! real(kind=dp), parameter :: MgO1= 28.271 ! Position of the superior MgO layer in bohr
  ! real(kind=dp), parameter :: MgO2= 54.378 ! Position of the inferior MgO layer in bohr
  real(kind=dp), parameter :: Fe1 = 32.098 ! Position of the superior Fe atom in bohr
  real(kind=dp), parameter :: Fe2 = 56.062 ! Position of the inferior Fe atom in bohr
  real(kind=dp), parameter :: MgO1= 28.271 ! Position of the superior MgO layer in bohr
  real(kind=dp), parameter :: MgO2= 59.888 ! Position of the inferior MgO layer in bohr
  !
  ! other variables
  integer :: i, j, k
  real(kind=dp) :: x, y, z
  !
  ! I/O variables

  !
  !
  ! Begining of the program
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
  !================================================================================
  !       read the parameters from the SCF QE calculation
  !================================================================================
  !
  write(stdout,20) '|      - Reading QE parameters                                  |'
  !
  call read_parameters_data_file_xml()
  !
  call get_gamma_only_haritz()
  !
  call get_spin_haritz()
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
  allocate(weight_z(nr3))
  !
  open(unit=91233, file='weights.dat', status="replace", action="write")
  write(91233,'(a15,x,a15,x,a15,x,a15,x,a15,x,a15,x,a15,x,a15)') 'ispin', 'iband', 'E-Ef', 'weight_iron_1', 'weight_iron_2', 'weight_mgo_1', 'weight_mgo_2', 'weight_silver'
  !
  do ikpoint=1,nkpoints_QE
    !
    do ispin=1,nspin
      !
      do ibands=1,nbands
        !
        call get_wfc()
        !
        call wfc_from_g_to_r_haritz(list_iG,wfc_g,wfc_r)
        !
        weight_iron_1 = 0.0_dp
        weight_iron_2 = 0.0_dp
        weight_mgo_1  = 0.0_dp
        weight_mgo_2  = 0.0_dp
        weight_silver = 0.0_dp
        weight_z      = 0.0_dp
        !
        do ir=1,nr1*nr2*nr3
          !
          call switch_indices_zyx(nr1,nr2,nr3,ir,i,j,k,-1)
          !
          x = ((i-1.0_dp)/nr1*at(1,1) + (j-1.0_dp)/nr2*at(1,2) + (k-1.0_dp)/nr3*at(1,3))*alat
          y = ((i-1.0_dp)/nr1*at(2,1) + (j-1.0_dp)/nr2*at(2,2) + (k-1.0_dp)/nr3*at(2,3))*alat
          z = ((i-1.0_dp)/nr1*at(3,1) + (j-1.0_dp)/nr2*at(3,2) + (k-1.0_dp)/nr3*at(3,3))*alat
          !
          weight_z(k) = weight_z(k) + abs(wfc_r(ir))**2/(nr1*nr2*nr3)
          !
          if ( abs(z-Fe1) < 2.0_dp ) then
            ! 2.0 Bhor up and down from the superior iron atom
            weight_iron_1 = weight_iron_1 + abs(wfc_r(ir))**2/(nr1*nr2*nr3)
            !
          elseif ( abs(z-Fe2) < 2.0_dp ) then
            ! 2.0 Bhor up and down from the inferior iron atom
            weight_iron_2 = weight_iron_2 + abs(wfc_r(ir))**2/(nr1*nr2*nr3)
            !
          elseif ( abs(z-MgO1) < 2.5_dp ) then
            ! 2.5 Bhor up and down from the superior MgO layer
            weight_mgo_1 = weight_mgo_1 + abs(wfc_r(ir))**2/(nr1*nr2*nr3)
            !
          elseif ( abs(z-MgO2) < 2.5_dp ) then
            ! 2.5 Bhor up and down from the inferior MgO layer
            weight_mgo_2 = weight_mgo_2 + abs(wfc_r(ir))**2/(nr1*nr2*nr3)
            !
          else
            !
            weight_silver  = weight_silver  + abs(wfc_r(ir))**2/(nr1*nr2*nr3)
            !
          endif
          !
        enddo ! ir
        !
        ! Write to a file the average weight on XY plane for each z
        do k=1,nr3
          if (ispin==1) write(10000+ibands,*) k, weight_z(k)
          if (ispin==1) write(20000+ibands,*) k, weight_z(k)
        enddo
        !
        write(91233,'(i15,x,i15,x,f15.8,x,f15.8,x,f15.8,x,f15.8,x,f15.8,x,f15.8)') ispin, ibands, eig(ibands,ispin,ikpoint)-e_fermi, weight_iron_1, weight_iron_2, weight_mgo_1, weight_mgo_2, weight_silver
        !
        write(stdout,   '(a,i2)') 'ispin: '  , ispin
        write(stdout,   '(a,i4)') 'iband: '  , ibands
        write(stdout,'(a,f12.8)') 'E-Ef: ', eig(ibands,ispin,ikpoint)-e_fermi
        write(stdout,'(a,f12.8)') 'iron 1: ', weight_iron_1
        write(stdout,'(a,f12.8)') 'iron 2: ', weight_iron_2
        write(stdout,'(a,f12.8)') 'mgo 1: ', weight_mgo_1
        write(stdout,'(a,f12.8)') 'mgo 2: ', weight_mgo_2
        write(stdout,'(a,f12.8)') 'silver: ' ,weight_silver
        write(stdout,*) ''
        !
        ! if (ibands==4) stop
        !
      enddo ! ibands
      !
      ! call close_K_wfc_file(wfc_unit)
      !
    enddo ! ispin
    !
  enddo ! ikpoint
  !
  close(91233)

  do k=1,nr3
    write(214978,*) k, weight_z(k)
  enddo




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


  subroutine get_wfc()
    !
    ! Reads the wfc's from QE .save
    !
    use intw_reading, only: write_tag
    use iotk_module
    !
    implicit none
    !
    ! local variables
    integer :: unit
    character(256) :: tag
    !
    !
    unit = find_free_unit()
    !
    call open_wfc_file(unit)
    ! read the wavefunction
    call write_tag("evc.",ibands,tag)
    call iotk_scan_dat (unit,tag,wfc_g(1:nG_K))
    wfc_g(nG_K+1:nG_max) = cmplx_0
    call close_wfc_file(unit)
    !
  end subroutine get_wfc


  subroutine open_wfc_file(unit)
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
    !
    ! local variables
    character(256) :: K_directory, Kdir
    character(256) :: wfc_file, gvec_file
    integer :: io_unit
    !
    !
    ! initialize the arrays to zero (zero will be broadcasted)
    !
    write(K_directory,100) ikpoint
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
       if (ispin==1) wfc_file=trim(Kdir)//'evc1.dat'
       if (ispin==2) wfc_file=trim(Kdir)//'evc2.dat'
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
  end subroutine open_wfc_file


  subroutine close_wfc_file(unit)
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
  end subroutine close_wfc_file


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
    NAMELIST / input_wfc_weight / mesh_dir, prefix
    !
    ! Set default values
    mesh_dir='./'
    prefix='Fe-MgO-Ag100'
    !
    INQUIRE(stdin, NAME=input_file_name)
    if (input_file_name(1:4)=="/dev") then
      ! there is no input file: use default values
      return
    else
      read(stdin,input_wfc_weight)
    endif
    !
  end subroutine read_input_file

end program wfc
