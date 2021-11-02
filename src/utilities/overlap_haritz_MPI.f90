program overlap_matrix

  ! MPI library
  use mpi

  ! Module variables
  use kinds, only: dp
  !
  use intw_useful_constants, only: cmplx_0, Ha_in_eV, zero
  !
  use intw_input_parameters, only: mesh_dir, prefix, nk1, nk2, nk3, &
                                   ph_dir, dvscf_dir, dvscf_name, &
                                   nq1, nq2, nq3, nqirr
  !
  use intw_reading, only: nspin, noncolin, ngm, gvec, nG_max, &
                          nbands, nkpoints_QE
  !
  use intw_pseudo, only: read_all_pseudo, average_pp
  !
  use intw_fft, only: nl
  !
  use intw_symmetries, only: full_mesh, IBZ
  !
  use intw_haritz, only: stdin, stdout, gamma_only, nlm

  ! Module subroutines and functions
  use intw_reading, only: read_parameters_data_file_xml, get_gvec, &
                          read_kpoints_data_file_xml
  !
  use intw_utility, only: find_free_unit, switch_indices_zyx, generate_kmesh
  !
  use w90_parameters, only: num_bands
  !
  use intw_symmetries, only: find_size_of_irreducible_k_set, find_the_irreducible_k_set
  !
  use intw_haritz, only: get_gamma_only_haritz, get_spin_haritz, &
                         generate_nl_haritz, wfc_from_g_to_r_haritz, &
                         wfc_from_g_to_r_haritz, get_e_fermi_haritz
  !
  use intw_mpi_haritz, only: master, nNodes, myNode

  !================================================================================
  !       Declare the variables
  !================================================================================

  implicit none

  ! Parameters
  character(len=1), parameter :: CR=CHAR(13), LF=CHAR(10)

  ! Local variables

  ! wfc related variables
  complex(kind=dp), allocatable :: all_wfc_k(:,:) ! all wave functions in G space (nG_max)
  complex(kind=dp), allocatable :: wfc_k(:,:,:) ! all selected wave functions in G space (nG_max)

  integer, allocatable :: list_iG(:) ! ZER DA HAU?

  ! ov_mat_el related variables
  complex(kind=dp), allocatable :: ov_mat_el(:,:,:,:) ! conjg(wfc_kq)*wfc_k
  integer :: ov_unit

  ! I/O variables
  integer :: ierr, record_length

  ! Loop variables
  integer :: jband ! Loop on band index (jband=1,nbands)
  integer :: ik ! Loop on k points
  integer :: jk ! Loop on k points
  integer :: iNode ! Loop on Nodes (iNode=1,nNodes)
  ! integer :: iG, nG

  ! Timing variables
  real(kind=dp) :: t0, t1

  !
  real(kind=dp) :: e_fermi
  real(kind=dp), allocatable :: eig(:,:,:)

  ! MPI variables
  integer            :: MPIerror
  integer            :: MPI_nbands_initial
  integer            :: MPI_nbands_final
  integer            :: nSize


  ! Input parameters
  integer :: nat_move ! The number of atoms moved in the calculation
  integer, dimension(100) :: ia_move_to_ia
  integer :: ispin ! Spin of the initial state
  integer :: iband ! Loop on band index (iband=1,nbands)
  integer :: nbands_initial, nbands_final, nbands_selected
  character(len=256) :: ov_mat_file

  real(kind=dp), allocatable, dimension(:,:) :: kpoints_irr
  real(kind=dp), allocatable, dimension(:,:) :: kmesh
  real(kind=dp), allocatable, dimension(:,:) :: kpoints_QE
  integer :: nk_irr, nkmesh


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
  call read_input_file()
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
  if (myNode==master) call get_gamma_only_haritz()
  call MPI_Bcast(gamma_only, 1, MPI_logical, master, MPI_COMM_WORLD, MPIerror )
  !
  if (myNode==master) call get_spin_haritz()
  call MPI_Bcast(nspin, 1, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
  !
  !================================================================================
  !       read in the kpoints from the QE folders
  !================================================================================
  !
  allocate(kpoints_QE(3,nkpoints_QE))
  !
  if (myNode==master) call read_kpoints_data_file_xml(kpoints_QE)
  call MPI_Bcast( kpoints_QE,    3*nkpoints_QE, MPI_double_precision, master, MPI_COMM_WORLD, MPIerror )
  !
  !
  !================================================================================
  ! Build the kmesh corresponding to the parameters in the input file
  !================================================================================
  !
  nkmesh = nk1*nk2*nk3
  allocate(kmesh(3,nkmesh))
  call generate_kmesh(kmesh,nk1,nk2,nk3)
  !
  !
  !================================================================================
  !    Find the size of the irreducible set of k-points (IBZ)
  !    and check that the number of kpoints corresponds to either
  !    a full mesh or the IBZ.
  !================================================================================
  !
  call find_size_of_irreducible_k_set(nk1,nk2,nk3,kmesh,nk_irr)
  !
  !This is only for testing: The result for nk_irr is different in both.
  !
  allocate(kpoints_irr(3,nkmesh))
  call find_the_irreducible_k_set (nk1,nk2,nk3,kmesh,kpoints_irr,nk_irr)
  !
  !
  if (nkpoints_QE/=nkmesh .and. nkpoints_QE/=nk_irr) then
     !
     ! the points in the folder are not consistent with the
     ! input file. Stop the program!
     !
     write(*,20) '*****************************************************'
     write(*,20) '*      The number of kpoints present in the QE      *'
     write(*,20) '*      folders are not consistent with a full       *'
     write(*,20) '*      Brillouin Zone or an irreducible Brillouin   *'
     write(*,20) '*      zone! Review your input...                   *'
     write(*,20) '*                   Program stops.                  *'
     write(*,20) '*****************************************************'
     write(*,20) '* debug information:                                *'
     write(*,*) '*        nkpoints_QE = ',nkpoints_QE
     write(*,*) '*        nkmesh      = ',nkmesh
     write(*,*) '*        nk_irr      = ',nk_irr
     !
     stop
     !
  elseif (nkpoints_QE==nk_irr) then
     !
     ! The points in the QE folders *could* be a valid choice for the IBZ;
     ! this must be checked!
     !
     full_mesh = .false.
     IBZ       = .true.
     !
  elseif (nkpoints_QE==nkmesh) then
     !
     ! The points in the QE folders *could* be consistent with a full mesh;
     ! this must be checked!
     !
     full_mesh = .true.
     IBZ       = .false.
     !
  endif
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
  !       Read eigenvalues and Fermi energy
  !================================================================================
  !
  if (myNode==master) then
    write(stdout,20) '|      - Read eigenvalues and Fermi energy                      |'
  endif
  !
  allocate( eig(nbands,nspin,nkpoints_QE) )
  !
  if (myNode==master) call get_eigenvalues()
  if (myNode==master) call get_e_fermi_haritz(e_fermi)
  call MPI_Bcast( e_fermi,     1, MPI_double_precision, master, MPI_COMM_WORLD, MPIerror )
  nSize=size(eig)
  call MPI_Bcast(     eig, nSize, MPI_double_precision, master, MPI_COMM_WORLD, MPIerror )
  !
  if (myNode==master) then
    write(stdout,20) '|---------------------------------------------------------------|'
  endif
  !
  ! do ik=1,nkpoints_QE
  !   print*, kpoints_QE(:,ik)
  !   print*, eig(:,1,ik)
  ! enddo
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
  allocate( all_wfc_k(nG_max,nbands) )
  allocate( wfc_k(nG_max,nbands_selected,nkpoints_QE) )
  list_iG  = 0
  eig      = ZERO
  all_wfc_k  = cmplx_0
  wfc_k = cmplx_0
  !
  allocate( ov_mat_el(nbands_selected,nbands_selected,nkpoints_QE,nkpoints_QE) )
  ov_mat_el=cmplx_0
  !
  if (myNode==master) then
    write(stdout,20) '|---------------------------------------------------------------|'
  endif
  !
  !
  !
  !================================================================================
  !       Read wave functions
  !================================================================================
  !
  if (myNode==master) write(stdout,20) '|      - Reading wave functions                                 |'
  !
  !
  do ik=1,nkpoints_QE
    !
    if (myNode==master) then
        write(stdout,"(a17,i4,a44)") '|      - k point ',ik,'                                           |'
      ! call cpu_time(t0)
      call read_wfc(ik)
      ! call cpu_time(t1)
      ! print*, t1-t0
      write(stdout,20) '|      - Bcast wave functions                                   |'
    endif
    !
    ! call cpu_time(t0)
    nSize=size(all_wfc_k)
    call MPI_Bcast( all_wfc_k, nSize, MPI_double_complex, master, MPI_COMM_WORLD, MPIerror )
    ! call cpu_time(t1)
    ! if (myNode==master) print*, t1-t0
    nSize=size(list_iG)
    call MPI_Bcast( list_iG, nSize,        MPI_integer, master, MPI_COMM_WORLD, MPIerror )
    !
    wfc_k(:,:,ik) = all_wfc_k(:,nbands_initial:nbands_final)
    !
  enddo
  !
  deallocate(all_wfc_k)
  !
  !
  if (myNode==master) then
    write(stdout,20) '|---------------------------------------------------------------|'
  endif
  !
  !
  !
  !
  !
  !================================================================================
  !       Calculate the overlap matrix
  !================================================================================
  !
  if (myNode==master) then
    write(stdout,20) '|      - Calculate the overlap matrix                           |'
  endif
  !
  !
  ! call cpu_time(t0)
  !
  !
  ! Calculate the matrix element < psi_kq | psi_k >
  do jk=1,nkpoints_QE
    do ik=1,nkpoints_QE
      do jband=MPI_nbands_initial,MPI_nbands_final
        do iband=1,nbands_selected
          call matrix_element(wfc_k(:,iband,ik),wfc_k(:,jband,jk),ov_mat_el(iband,jband,ik,jk))
        enddo ! iband
      enddo ! jband
    enddo ! ik
  enddo ! jk
  !
  !
  ! call cpu_time(t1)
  ! print*, "Time:", t1-t0
  !
  !
  call MPI_Barrier(MPI_COMM_WORLD, MPIerror)
  !
  !
  ! Bcast matrix elements
  do iNode=0,nNodes-1
    if ( mod(nbands_selected,nNodes)==0 ) then
      MPI_nbands_initial = int(0.5+1.0*nbands_selected/nNodes)*iNode + 1
      MPI_nbands_final   = int(0.5+1.0*nbands_selected/nNodes)*(iNode + 1)
    else
      MPI_nbands_initial = nint(0.5+1.0*nbands_selected/nNodes)*iNode + 1
      MPI_nbands_final   = nint(0.5+1.0*nbands_selected/nNodes)*(iNode + 1)
    endif
    if (iNode+1==nNodes) MPI_nbands_final=nbands_selected
    nSize=size(ov_mat_el(:,MPI_nbands_initial:MPI_nbands_final,:,:))
    call MPI_Bcast( ov_mat_el(:,MPI_nbands_initial:MPI_nbands_final,:,:), nSize, MPI_double_complex, iNode, MPI_COMM_WORLD, MPIerror )
  enddo
  call MPI_Barrier(MPI_COMM_WORLD, MPIerror)
  !
  if (myNode==master) then
    !
    do ik=1,nkpoints_QE
      do jk=1,nkpoints_QE
        print"(a,3f7.4,a,3f7.4)", "kpoint1:", kmesh(:,ik), " kpoint2:", kmesh(:,jk)
        do jband=1,nbands_selected
          do iband=1,nbands_selected
            print"(a6,2i5,4x,2f10.5,4x,f10.5)", "haritz", nbands_initial+iband-1, nbands_initial+jband-1, ov_mat_el(iband,jband,ik,jk), abs(ov_mat_el(iband,jband,ik,jk))
          enddo
        enddo
      enddo
    enddo
    !
    !
    ov_unit=find_free_unit()
    inquire(iolength=record_length) ov_mat_el(:,:,1,1)
    open( unit=ov_unit, file=trim(trim(mesh_dir)//trim(ph_dir)//trim(ov_mat_file)//trim("_1")), iostat=ierr, &
          form='unformatted', status='new', access='direct', recl=record_length )
    if (ierr/=0) then
      if (ierr==10) then
        ! Specified file already exists
        write(stdout,"(a)") "WARNING: ov_mat_file already exists, open with another name."
        ! open a new file with another name
        open( unit=ov_unit, file=trim(trim(mesh_dir)//trim(ph_dir)//trim(ov_mat_file)//trim("_1_BERRIA")), iostat=ierr, &
              form='unformatted', status='new', access='direct', recl=record_length )
        if (ierr/=0) then
          write(stdout,"(a)") "ERROR: ov_mat_file could not be opened."
          stop
        endif
      else
        write(stdout,"(a)") "ERROR: ov_mat_file could not be opened."
        stop
      endif
    endif
    !
    do ik=1,nkpoints_QE
      do jk=1,nkpoints_QE
        write( unit=ov_unit, rec=ik ) ov_mat_el(:,:,ik,jk)
      enddo
    enddo
    !
    close(ov_unit)
    !
  endif
  !

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine matrix_element(bra_wfc, ket_wfc, ep)
    !
    !
    !
    implicit none
    !
    ! I/O variables
    complex(kind=dp), intent(in)  :: bra_wfc(:), ket_wfc(:)
    complex(kind=dp), intent(out) :: ep
    !
    ! Local variables
    integer :: ig
    !
    20 format(A)
    !
    ! Check if the dimensions are equal
    if ( size(bra_wfc) == size(ket_wfc) ) then
      if ( size(bra_wfc) /= nG_max ) then
        write(stdout,20) "matrix_element: The dimension of the wfc's is not nG_max"
        stop
      endif
    else
      write(stdout,20) "matrix_element: The dimensions of the wfc's are not the same"
      stop
    endif
    !
    if (gamma_only) then
      do ig=1,nG_max
        ep = ep + 2.0_dp*real( conjg(bra_wfc(ig))*ket_wfc(ig) )
      enddo
      ep = ep - conjg(bra_wfc(1))*ket_wfc(1)
    else
      do ig=1,nG_max
        ep = ep + conjg(bra_wfc(ig))*ket_wfc(ig)
      enddo
    endif
    !
  end subroutine matrix_element

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_wfc(k_point)
    !
    ! Reads the wfc's from QE .save
    !
    use iotk_module, only: iotk_open_read, iotk_close_read, iotk_scan_dat
    use intw_reading, only: write_tag
    !
    implicit none

    ! Input variables

    integer, intent(in) :: k_point

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
    write(K_directory,100) k_point
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
    call iotk_open_read(io_unit,wfc_file(ispin))
    !
    do band=1,nbands
      !
      call write_tag("evc.",band,tag)
      call iotk_scan_dat (io_unit,tag,all_wfc_k(1:nG,band))
      !
    enddo
    !
    call iotk_close_read(io_unit)
    !
    return
    !
  end subroutine read_wfc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_eigenvalues()
    !
    ! Reads the eigenvalues from QE .save
    !

    use iotk_module

    implicit none

    !local variables
    character(256) :: K_directory, Kdir
    ! character(256) :: tag
    character(256), dimension(nspin) :: eig_file
    integer :: io_unit, is, ik
    ! integer :: nG
    !
    !
    do ik=1,nkpoints_QE
      !
      write(K_directory,100) ik
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
         call iotk_scan_dat (io_unit,"EIGENVALUES",eig(:,is,ik))
         !
         call iotk_close_read(io_unit)
         !
      enddo ! is
      !
    enddo ! ik
    !
    ! convert eigenvalues from a.u. to eV
    !
    eig=eig*Ha_in_eV
    !
    100 format('K'I5.5'/')
    !
    return

  end subroutine get_eigenvalues

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
    NAMELIST / input_ov_melements / mesh_dir, prefix, ph_dir, dvscf_dir, dvscf_name, &
                                    ov_mat_file, iband, ispin, nat_move, ia_move_to_ia, &
                                    nbands_initial, nbands_final, nk1, nk2, nk3
    !
    ! Set default values
    mesh_dir = '/home/haritz/Kalkuluak/Probak/Fe-O/K/'
    prefix = 'fe-o'
    nk1 = 1
    nk2 = 1
    nk3 = 1
    ph_dir='./'
    dvscf_dir='./dvscf_dir/'
    dvscf_name='fe-o.dvscf_q1'
    ov_mat_file=trim("ov_mat_is_1_iband_1274_VL.dat")
    nq1=1
    nq2=1
    nq3=1
    nqirr=1
    nat_move = 0
    ia_move_to_ia = 0
    ! Define the initial state
    iband = 1
    ispin = 1
    ! Select final states
    nbands_initial = 1
    nbands_final = 1
    ! way to calculate the matrix elements
    !
    if (myNode==master) then
      INQUIRE(stdin, NAME=input_file_name)
      if (input_file_name(1:4)=="/dev") then
        ! there is no input file: use default values
        return
      else
        read(stdin,input_ov_melements)
      endif
      !
      ! Reopen terminal as input file
      ! close(unit=stdin)
      ! open(unit=stdin,file='/dev/tty')
    endif
    !
    nSize=len(mesh_dir)
    call MPI_Bcast( mesh_dir, nSize, MPI_char, master, MPI_COMM_WORLD, MPIerror )
    nSize=len(prefix)
    call MPI_Bcast( prefix, nSize, MPI_char, master, MPI_COMM_WORLD, MPIerror )
    nSize=len(ph_dir)
    call MPI_Bcast( ph_dir, nSize, MPI_char, master, MPI_COMM_WORLD, MPIerror )
    nSize=len(dvscf_dir)
    call MPI_Bcast( dvscf_dir, nSize, MPI_char, master, MPI_COMM_WORLD, MPIerror )
    nSize=len(dvscf_name)
    call MPI_Bcast( dvscf_name, nSize, MPI_char, master, MPI_COMM_WORLD, MPIerror )
    nSize=len(ov_mat_file)
    call MPI_Bcast( ov_mat_file, nSize, MPI_char, master, MPI_COMM_WORLD, MPIerror )
    call MPI_Bcast( iband, 1, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
    call MPI_Bcast( ispin, 1, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
    call MPI_Bcast( nat_move, 1, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
    nSize=size(ia_move_to_ia)
    call MPI_Bcast( ia_move_to_ia, nSize, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
    call MPI_Bcast( nbands_initial, 1, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
    call MPI_Bcast( nbands_final, 1, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
    call MPI_Bcast( nk1, 1, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
    call MPI_Bcast( nk2, 1, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
    call MPI_Bcast( nk3, 1, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
    !
    nbands_selected = nbands_final - nbands_initial + 1
    num_bands = nbands_selected
    !
    ! Find the band ranges for each thread
    if ( mod(nbands_selected,nNodes)==0 ) then
      MPI_nbands_initial = int(0.5+1.0*nbands_selected/nNodes)*myNode + 1
      MPI_nbands_final   = int(0.5+1.0*nbands_selected/nNodes)*(myNode + 1)
    else
      MPI_nbands_initial = nint(0.5+1.0*nbands_selected/nNodes)*myNode + 1
      MPI_nbands_final   = nint(0.5+1.0*nbands_selected/nNodes)*(myNode + 1)
    endif
    if (myNode+1==nNodes) MPI_nbands_final=nbands_selected
    num_bands = MPI_nbands_final - MPI_nbands_initial + 1
    call MPI_Barrier(MPI_COMM_WORLD, MPIerror)
    write(stdout,'(a15,i2,a8,i4,a1,i4,a6,i4,a21)') '|      - Node: ',myNode,' bands: ',MPI_nbands_initial,'-',MPI_nbands_final,' tot: ',num_bands,'                    |'
    call MPI_Barrier(MPI_COMM_WORLD, MPIerror)
    !
  end subroutine read_input_file


end program overlap_matrix
