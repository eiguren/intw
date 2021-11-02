program ep_melements

  ! MPI library
  use mpi

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

  use intw_mpi_haritz, only: master, nNodes, myNode

  ! Module subroutines and functions
  use intw_reading, only: read_parameters_data_file_xml, get_gvec

  use intw_utility, only: find_free_unit, switch_indices_zyx
  !
  use intw_haritz, only: gamma_only, get_gamma_only_haritz, get_spin_haritz, &
                         get_e_fermi_haritz, nlm, generate_nl_haritz, wfc_from_g_to_r_haritz
  !================================================================================
  !       Declare the variables
  !================================================================================

  implicit none

  ! Parameters
  integer, parameter :: stdout=6
  character(len=1), parameter :: CR=CHAR(13), LF=CHAR(10)

  ! Local variables

  ! FFT variables
  integer, allocatable :: nl_all(:)


  ! wfc related variables
  ! real(kind=dp) :: eig

  complex(kind=dp), allocatable :: all_wfc(:,:) ! all wave functions in G space (nG_max)
  complex(kind=dp), allocatable :: wfc_kq_g(:)  ! wave function k+q in G space (nG_max)
  complex(kind=dp), allocatable :: wfc_kq_r(:)  ! wave function k+q in real space (nr1*nr2*nr3)
  complex(kind=dp), allocatable :: wfc_k_g(:)   ! wave function k in G space (nG_max)
  complex(kind=dp), allocatable :: wfc_k_r(:)   ! wave function k in real space (nr1*nr2*nr3)

  integer, allocatable :: list_iG(:) ! ZER DA HAU?

  ! dv related variables
  real(kind=dp), allocatable :: dv_local(:) ! local part of the variation of the potential in real space (nr1*nr2*nr3)

  ! ep_mat_el related variables
  real(kind=dp), allocatable :: ep_mat_el(:,:) ! conjg(wfc_kq)*dv_q*wfc_k
  integer :: ep_unit

  ! I/O variables
  integer :: ierr, record_length, ios, ounit
  ! real(kind=dp), allocatable, dimension(:,:,:) :: rho


  ! Loop variables
  integer :: is    ! Loop on spin components (is=1,nspin)
  integer :: iband ! Loop on band index (iband=1,nbands)
  integer :: jband ! Loop on band index (jband=1,nbands)
  integer :: imode ! Loop on vibrational modes (imode=1,3*nat)
  integer :: iNode ! Loop on Nodes (iNode=1,nNodes)

  ! Timing variables
  real(kind=dp) :: t0, t1, t2

  integer :: i, x, y, z
  real(kind=dp) :: e_fermi
  real(kind=dp), allocatable :: eig(:,:)

  ! MPI variables
  integer            :: MPIerror
  integer            :: nbands_initial
  integer            :: nbands_final
  integer            :: nSize




  20 format(A)
  30 format(A,F8.2,6X,A)
  !
  !
  !================================================================================
  !       Talk to the user
  !================================================================================
  !
  !
  !
  ! Launch MPI
  call MPI_INIT(MPIerror)
  !
  call MPI_COMM_RANK(MPI_COMM_WORLD, myNode, MPIerror)
  !
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nNodes, MPIerror)
  !
  ! Close std output and open mep.out file to this purpose
  if (myNode==master) then
    write(stdout,20) '================================================================='
    write(stdout,20) '|                         program me                            |'
    write(stdout,20) '|---------------------------------------------------------------|'
    write(stdout,20) '================================================================='
    write(stdout,'(a30,i2,a33)') '|      - MPI initialized with ', nNodes, ' nodes                          |'
    write(stdout,20) '|---------------------------------------------------------------|'
  endif
  !
  !================================================================================
  !       Specify input parameters
  !================================================================================
  !
  ! This is the nspin=1 calculation
  ! Then we will be interested with the nspin=2 calculation
  mesh_dir = '/scratch/hgarai/Kalkuluak/Fe-MgO-Ag100/1ML/ONE_FE/SPIN/'
  prefix = 'Fe-MgO-Ag100'
  nk1 = 1
  nk2 = 1
  nk3 = 1
  calc_epmat=.true.
  ph_dir='PH/'
  dvscf_dir='./'
  dvscf_name='fe-mgo-ag100.dvscf_q1_PP'
  ep_mat_file=trim("ep_mat_is_1_iband_1274_VL.dat")
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
  ! call MPI_Bcast(       nspin,      1,          MPI_integer, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast(    noncolin,      1,          MPI_logical, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast(    lspinorb,      1,          MPI_logical, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast(        ntyp,      1,          MPI_integer, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast(         ngm,      1,          MPI_integer, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast(      nG_max,      1,          MPI_integer, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast(         nr1,      1,          MPI_integer, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast(         nr2,      1,          MPI_integer, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast(         nr3,      1,          MPI_integer, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast(         nat,      1,          MPI_integer, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast(      nbands,      1,          MPI_integer, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast(        ityp,    nat,          MPI_integer, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast( atom_labels, 3*ntyp,             MPI_char, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast(         tau,  3*nat,          MPI_integer, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast(          at,    3*3, MPI_double_precision, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast(        alat,      1, MPI_double_precision, master, MPI_COMM_WORLD, MPIerror )
  ! call MPI_Bcast(     volume0,      1, MPI_double_precision, master, MPI_COMM_WORLD, MPIerror )
  !
  !
  !
  call get_gamma_only_haritz()
  !
  call get_spin_haritz()
  !
  !
  !
  !================================================================================
  !       Some checks
  !================================================================================
  !
  if (myNode==master) then
    if (gamma_only) then
      write(stdout,20) '|      - gamma_only: .TRUE. calculation                         |'
    else
      write(stdout,20) '|      - gamma_only: .FALSE. calculation                        |'
    endif
    write(stdout,20) '|---------------------------------------------------------------|'
  endif
  !
  if (myNode==master) then
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
       call MPI_FINALIZE(MPIerror)
       stop
       !
    endif
    !
  endif
  !
  !
  !
  !================================================================================
  !       Read PP files
  !================================================================================
  !
  if (myNode==master) then
    write(stdout,20) '|      - Reading pseudopotentials from UPF files                |'
  endif
  !
  ! call read_all_pseudo ()
  ! if (.not.lspinorb) call average_pp(ntyp)
  !
  if (myNode==master) then
    write(stdout,20) '|                          PPs are OK                           |'
    write(stdout,20) '|---------------------------------------------------------------|'
  endif
  !
  !
  !================================================================================
  !       Allocate and read G vectors
  !================================================================================
  !
  if (myNode==master) then
    write(stdout,20) '|      - Reading G vectors                                      |'
  endif
  !
  ! Read all g vectors
  allocate(gvec(3,ngm))
  if (myNode==master) call get_gvec() ! read them from gvectors.dat file
  nSize=size(gvec)
  call MPI_Bcast( gvec, nSize, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
  !
  ! Calculate nl
  allocate (nl(ngm), nlm(ngm))
  call generate_nl_haritz()
  !
  if (myNode==master) then
    write(stdout,20) '|---------------------------------------------------------------|'
  endif
  !
  !
  !
  !================================================================================
  !       Allocate wfc related variables
  !================================================================================
  !
  if (myNode==master) then
    write(stdout,20) '|      - Allocate wfc related variables                         |'
  endif
  !
  allocate( list_iG(nG_max) )
  allocate( eig(nbands,nspin) )
  allocate( all_wfc(nG_max,nbands), wfc_k_g(nG_max), wfc_kq_g(nG_max) )
  allocate( wfc_k_r(nr1*nr2*nr3), wfc_kq_r(nr1*nr2*nr3) )
  list_iG  = 0
  eig      = ZERO
  all_wfc  = cmplx_0
  wfc_k_g  = cmplx_0
  wfc_kq_g = cmplx_0
  wfc_k_r  = cmplx_0
  wfc_kq_r = cmplx_0
  !
  if (myNode==master) then
    write(stdout,20) '|---------------------------------------------------------------|'
  endif
  !
  !
  !================================================================================
  !       Allocate dv related variables
  !================================================================================
  !
  if (myNode==master) then
    write(stdout,20) '|      - Allocate dv related variables                          |'
  endif
  !
  allocate( dv_local(nr1*nr2*nr3) )
  allocate( ep_mat_el(nbands,3*nat) )
  dv_local=cmplx_0
  ep_mat_el=ZERO
  !
  if (myNode==master) then
    write(stdout,20) '|---------------------------------------------------------------|'
  endif
  !
  !
  !================================================================================
  !       Read eigenvalues and Fermi energy
  !================================================================================
  !
  if (myNode==master) then
    write(stdout,20) '|      - Read eigenvalues and Fermi energy                      |'
  endif

  if (myNode==master) call get_GAMMA_eigenvalues()
  if (myNode==master) call get_e_fermi_haritz(e_fermi)
  call MPI_Bcast( e_fermi,     1, MPI_double_precision, master, MPI_COMM_WORLD, MPIerror )
  nSize=size(eig)
  call MPI_Bcast(     eig, nSize, MPI_double_precision, master, MPI_COMM_WORLD, MPIerror )
  !
  if (myNode==master) then
    write(stdout,20) '|---------------------------------------------------------------|'
  endif
  !
  !
  ! Define the initial state
  iband = 1274
  is = 1
  !
  !
  !================================================================================
  !       Read the initial wave function
  !================================================================================
  !
  if (myNode==master) then
    write(stdout,20) '|      - Reading wave functions                                 |'
    call cpu_time(t0)
    call read_wfc()
    call cpu_time(t1)
    print*, t1-t0
    write(stdout,20) '|      - Bcast wave functions                                   |'
  endif
  !
  call cpu_time(t0)
  nSize=size(all_wfc)
  call MPI_Bcast( all_wfc, nSize, MPI_double_complex, master, MPI_COMM_WORLD, MPIerror )
  call cpu_time(t1)
  if (myNode==master) print*, t1-t0
  nSize=size(list_iG)
  call MPI_Bcast( list_iG, nSize,        MPI_integer, master, MPI_COMM_WORLD, MPIerror )
  wfc_kq_g = all_wfc(:,iband)
  !
  if (myNode==master) then
    write(stdout,20) '|---------------------------------------------------------------|'
  endif
  !
  !
  !================================================================================
  !       Fourier transform the initial wave function to r space
  !================================================================================
  !
  if (myNode==master) then
    write(stdout,20) '|      - Fourier transform the initial wave function to r space |'
  endif
  !
  call wfc_from_g_to_r_haritz( list_iG, wfc_kq_g, wfc_kq_r )
  !
  if (myNode==master) then
    write(stdout,20) '|---------------------------------------------------------------|'
  endif
  !
  !
  !================================================================================
  !       Calculate the matrix elements
  !================================================================================
  !
  if (myNode==master) then
    write(stdout,20) '|      - Calculate the matrix elements                          |'
  endif
  !
  ! Find the band ranges for each thread
  nbands_initial = nint(0.5+1.0*nbands/nNodes)*myNode + 1
  nbands_final   = nint(0.5+1.0*nbands/nNodes)*(myNode + 1)
  if (myNode+1==nNodes) nbands_final=nbands
  write(stdout,'(a,i2,a,i4,a,i4,a)') '|      - Node: ',myNode,' bands: ',nbands_initial,'-',nbands_final,'                 |'
  !
  if (calc_epmat) then
    !
    !
    ! For each mode
    do imode=412,3*nat
      !
      ! Read from the dvscf file the induced part of the variation of the potential
      if (myNode==master) call read_dv(imode, dv_local)
      nSize=size(dv_local)
      call MPI_Bcast( dv_local, nSize, MPI_double_precision, master, MPI_COMM_WORLD, MPIerror )
      if (myNode==master) write(stdout,20) '|      - Induced part of dV readed                              |'
      !
      ! For each state
      do jband=nbands_initial,nbands_final
        !
        ! Talk to the user
        ! write(unit=stdout,fmt="(a7,i3,a1,i3,x,a5,i7,a1,i7)", advance="no") &
        !       CR//"Mode: ", imode,"/",3*nat, "wfc: ", is*jband,"/",(nbands)
        ! flush(stdout)
        ! write(unit=stdout,fmt="(a1)", advance="no") LF
        !
        wfc_k_g = all_wfc(:,jband)
        !
        ! Fourier transform the wave function to r space
        call wfc_from_g_to_r_haritz( list_iG, wfc_k_g, wfc_k_r )
        !
        ! Calculate the matrix element
        call matrix_element(real(wfc_kq_r),dv_local,real(wfc_k_r))
        !
        write(stdout,"(i2,a8,i3,a1,i3,x,a5,i7,a1,i7)") &
              myNode, " Mode: ", imode,"/",3*nat, "wfc: ", jband,"/",(nbands_final)

        !
      enddo ! jband
      !
    enddo ! imode
    !
    ! Bcast matrix elements
    do iNode=0,nNodes-1
      nbands_initial = nint(0.5+1.0*nbands/nNodes)*iNode + 1
      nbands_final   = nint(0.5+1.0*nbands/nNodes)*(iNode + 1)
      if (iNode+1==nNodes) nbands_final=nbands
      nSize=size(ep_mat_el(nbands_initial:nbands_final,:))
      call MPI_Bcast( ep_mat_el(nbands_initial:nbands_final,:), nSize, MPI_double_precision, iNode, MPI_COMM_WORLD, MPIerror )
    enddo
    !
    if (myNode==master) then
      !
      ep_unit=find_free_unit()
      inquire(iolength=record_length) ep_mat_el(:,1)
      open( unit=ep_unit, file=trim(trim(mesh_dir)//trim(ph_dir)//trim(ep_mat_file)//trim("_2")), iostat=ierr, &
            form='unformatted', status='new', access='direct', recl=record_length )
      !
      do imode=1,3*nat
        write( unit=ep_unit, rec=imode ) ep_mat_el(:,imode)
      enddo
      close(ep_unit)
      !
    endif
    !
  else
    ! Read them
  endif


  call MPI_FINALIZE(MPIerror)

contains

  subroutine stop_mpi()
    !
    !
    !
    implicit none
    !
    stop

  end subroutine stop_mpi

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
    ep_mat_el(jband,imode) = ZERO
    do ir=1,size(bra_wfc)
      ep_mat_el(jband,imode) = ep_mat_el(jband,imode) + (bra_wfc(ir))*dv(ir)*(ket_wfc(ir))
    enddo
    !
    ep_mat_el(jband,imode) = ep_mat_el(jband,imode)/size(bra_wfc)
    !
  end subroutine matrix_element

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

  subroutine read_wfc()
    !
    ! Reads the wfc's from QE .save
    !
    use iotk_module, only: iotk_open_read, iotk_close_read, iotk_scan_dat
    use intw_reading, only: write_tag
    !
    implicit none

    ! Local variables

    character(256) :: K_directory, Kdir
    character(256) :: tag, gvec_file
    character(256), dimension(nspin) :: wfc_file
    integer :: io_unit
    integer :: nG, band
    !
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
    call iotk_open_read(io_unit,wfc_file(is))
    !
    do band=1,nbands
      !
      call write_tag("evc.",band,tag)
      call iotk_scan_dat (io_unit,tag,all_wfc(1:nG,band))
      !
    enddo
    !
    call iotk_close_read(io_unit)
    !
    return
    !
  end subroutine read_wfc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
