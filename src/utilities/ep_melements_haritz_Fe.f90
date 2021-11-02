program ep_melements

  ! MPI library
  use mpi

  ! Module variables
  use kinds, only: dp
  !
  use intw_useful_constants, only: cmplx_0, Ha_in_eV, zero
  !
  use intw_input_parameters, only: mesh_dir, prefix, nk1, nk2, nk3, &
                                   ph_dir, dvscf_dir, dvscf_name, ep_mat_file, &
                                   calc_epmat, nq1, nq2, nq3, nqirr
  !
  use intw_reading, only: nspin, noncolin, lspinorb, ntyp, ngm, gvec, nG_max, &
                          nat, nbands
  !
  use intw_pseudo, only: read_all_pseudo, average_pp
  !
  use intw_fft, only: nl
  !
  use intw_haritz, only: stdin, stdout, gamma_only, nlm
  !
  use intw_mpi_haritz, only: master, nNodes, myNode

  ! Module subroutines and functions
  use intw_reading, only: read_parameters_data_file_xml, get_gvec
  !
  use intw_utility, only: find_free_unit, switch_indices_zyx
  !
  use intw_haritz, only: get_gamma_only_haritz, get_spin_haritz, &
                         generate_nl_haritz, wfc_from_g_to_r_haritz, &
                         wfc_from_g_to_r_haritz, get_e_fermi_haritz
  !
  use intw_uspp, only: vkb, vkqb
  !
  use intw_ph, only: dvpsi
  !
  use w90_parameters, only: num_bands

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
  complex(kind=dp), allocatable :: wfc_kq(:,:)  ! wave function k+q in G space (nG_max)

  integer, allocatable :: list_iG(:) ! ZER DA HAU?

  ! ep_mat_el related variables
  complex(kind=dp), allocatable :: ep_mat_el(:,:,:) ! conjg(wfc_kq)*dv_q*wfc_k
  integer :: ep_unit

  ! I/O variables
  integer :: ierr, record_length

  ! Loop variables
  integer :: jband ! Loop on band index (jband=1,nbands)
  integer :: imode ! Loop on vibrational modes (imode=1,3*nat)
  integer :: iNode ! Loop on Nodes (iNode=1,nNodes)
  integer :: iG, nG

  ! Timing variables
  real(kind=dp) :: t0, t1

  !
  real(kind=dp) :: e_fermi
  real(kind=dp), allocatable :: eig(:,:)

  ! MPI variables
  integer            :: MPIerror
  integer            :: MPI_nbands_initial
  integer            :: MPI_nbands_final
  integer            :: nSize

  !
  integer :: npw

  ! Input parameters
  integer :: nat_move ! The number of atoms moved in the calculation
  integer, dimension(100) :: ia_move_to_ia
  integer :: ispin ! Spin of the initial state
  integer :: iband ! Loop on band index (iband=1,nbands)
  integer :: nbands_initial, nbands_final, nbands_selected

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
  call read_all_pseudo ()
  if (.not.lspinorb) call average_pp(ntyp)
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
  !       Read eigenvalues and Fermi energy
  !================================================================================
  !
  if (myNode==master) then
    write(stdout,20) '|      - Read eigenvalues and Fermi energy                      |'
  endif
  !
  allocate( eig(nbands,nspin) )
  !
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
  allocate( wfc_kq(nG_max,nbands_selected) )
  allocate( wfc_k(nG_max,num_bands,1) )
  list_iG  = 0
  eig      = ZERO
  all_wfc_k  = cmplx_0
  wfc_kq = cmplx_0
  wfc_k = cmplx_0
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
  nSize=size(all_wfc_k)
  call MPI_Bcast( all_wfc_k, nSize, MPI_double_complex, master, MPI_COMM_WORLD, MPIerror )
  call cpu_time(t1)
  if (myNode==master) print*, t1-t0
  nSize=size(list_iG)
  call MPI_Bcast( list_iG, nSize,        MPI_integer, master, MPI_COMM_WORLD, MPIerror )
  !
  wfc_k(:,:,1) = all_wfc_k(:,MPI_nbands_initial:MPI_nbands_final)
  wfc_kq(:,:) = all_wfc_k(:,nbands_initial:nbands_final)
  deallocate(all_wfc_k)
  !
  if (myNode==master) then
    write(stdout,20) '|---------------------------------------------------------------|'
  endif
  !
  !
  !
  !================================================================================
  !       Allocate dv related variables
  !================================================================================
  !
  if (myNode==master) then
    write(stdout,20) '|      - Allocate and calculate dv related variables            |'
  endif
  !
  !
  call allocate_nlpot1
  !
  call allocate_phq
  !
  call init_us_1
  !
  call allocate_nlpot2 ( (/ 0.0_dp, 0.0_dp, 0.0_dp /) )
  !
  call phq_init( (/ 0.0_dp, 0.0_dp, 0.0_dp /) )
  !
  nG=0
  do iG=1,nG_max
     if (list_iG(iG)==0) exit
     nG=nG+1
  enddo
  npw=nG
  !
  vkb=cmplx_0
  call init_us_2(npw, nG_max, list_iG, (/ 0.0_dp, 0.0_dp, 0.0_dp /), vkb)
  vkqb = vkb
  !
  allocate( ep_mat_el(nbands_selected,nbands_selected,3*nat) )
  ep_mat_el=cmplx_0
  !
  allocate(dvpsi(nG_max,num_bands,1,1,3*nat))
  dvpsi = cmplx_0
  !
  if (myNode==master) then
    write(stdout,20) '|---------------------------------------------------------------|'
  endif
  !
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
  if (calc_epmat) then
    !
    call cpu_time(t0)
    !
    !
    ! Calculate  dV_loc | psi_k >
    call calculate_dvloc_psi( nat, nG_max, num_bands, list_iG, wfc_k, dvpsi )
    !
    ! Calculate  dV_nl | psi_k >
    call dvqpsi_us_only_haritz( nat, nG_max, num_bands, 1, list_iG, wfc_k, dvpsi )
    !
    ! Calculate the matrix element < psi_kq | dV_nl | psi_k >
    do jband=1,num_bands
      !
      do imode=1,3*nat

        do iband=1,nbands_selected
          !
          call matrix_element(wfc_kq(:,iband),dvpsi(:,jband,1,1,imode),ep_mat_el(iband,MPI_nbands_initial+jband-1,imode))
          !
        enddo
        !
      enddo
      !
    enddo
    !
    !
    call cpu_time(t1)
    print*, "Time:", t1-t0
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
      do imode=1,3*nat
        nSize=size(ep_mat_el(:,MPI_nbands_initial:MPI_nbands_final,imode))
        call MPI_Bcast( ep_mat_el(:,MPI_nbands_initial:MPI_nbands_final,imode), nSize, MPI_double_complex, iNode, MPI_COMM_WORLD, MPIerror )
      enddo
    enddo
    call MPI_Barrier(MPI_COMM_WORLD, MPIerror)
    !
    if (myNode==master) then
      !
      do imode=1,3*nat
        do jband=1,nbands_selected
          do iband=1,nbands_selected
            print"(a6,3i5,4x,2f15.10,4x,f15.10)", "haritz", nbands_initial+iband-1, nbands_initial+jband-1, imode, ep_mat_el(iband,jband,imode), abs(ep_mat_el(iband,jband,imode))
          enddo
        enddo
      enddo
      !
      !
      ep_unit=find_free_unit()
      inquire(iolength=record_length) ep_mat_el(:,:,1)
      open( unit=ep_unit, file=trim(trim(mesh_dir)//trim(ph_dir)//trim(ep_mat_file)//trim("_1")), iostat=ierr, &
            form='unformatted', status='new', access='direct', recl=record_length )
      if (ierr/=0) then
        if (ierr==10) then
          ! Specified file already exists
          write(stdout,"(a)") "WARNING: ep_mat_file already exists, open with another name."
          ! open a new file with another name
          open( unit=ep_unit, file=trim(trim(mesh_dir)//trim(ph_dir)//trim(ep_mat_file)//trim("_1_BERRIA")), iostat=ierr, &
                form='unformatted', status='new', access='direct', recl=record_length )
          if (ierr/=0) then
            write(stdout,"(a)") "ERROR: ep_mat_file could not be opened."
            stop
          endif
        else
          write(stdout,"(a)") "ERROR: ep_mat_file could not be opened."
          stop
        endif
      endif
      !
      do imode=1,3*nat
        write( unit=ep_unit, rec=imode ) ep_mat_el(:,:,imode)
      enddo
      !
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calculate_dvloc_psi(nat,nG_max,nbands,igk,wfc_k_g,dvpsi)
    !
    !
    !
    use intw_reading, only: nr1, nr2, nr3
    use intw_haritz, only: wfc_from_g_to_r_haritz, wfc_from_r_to_g_haritz
    !
    implicit none
    !
    integer, intent(in) :: nat, nG_max, nbands
    integer, dimension(nG_max), intent(in) :: igk
    complex(kind=dp), dimension(nG_max,nbands,1), intent(in) :: wfc_k_g
    complex(kind=dp), dimension(nG_max,nbands,1,1,3*nat), intent(inout) :: dvpsi
    !
    complex(kind=dp), dimension(nr1*nr2*nr3,nbands,1) :: wfc_k_r
    complex(kind=dp), dimension(nr1*nr2*nr3) :: dv_local
    real(kind=dp), dimension(nr1*nr2*nr3) :: dv_scf
    complex(kind=dp), dimension(nr1*nr2*nr3) :: aux
    real(kind=dp), dimension(3) :: qpoint
    !
    integer :: iband, ir, imode, ia, id, ia_move, imode_move

    integer :: rl, rec, read_unit, read_unit_1, read_unit_2, ios
    character(len=256) :: dvfile
    !
    !
    qpoint=0.0_dp
    !
    ! Fourier transform the wave functions to r space
    do iband=1,nbands
      call wfc_from_g_to_r_haritz( igk, wfc_k_g(:,iband,1), wfc_k_r(:,iband,1) )
    enddo
    !
    dv_local = cmplx_0
    dv_scf = ZERO
    !
    ! Record length of the variable to read
    inquire(iolength=rl) dv_scf
    !
    ! Open the file
    read_unit_1=find_free_unit()
    dvfile = trim(mesh_dir) // trim(ph_dir) // trim(dvscf_dir) // "mgo-ag100_R441_fine.dvscf_q1"
    open( unit=read_unit_1, file=trim(dvfile), iostat=ios, form='unformatted', &
          status='old', action='read', access='direct', recl=rl )
    !
    read_unit_2=find_free_unit()
    dvfile = trim(mesh_dir) // trim(ph_dir) // trim(dvscf_dir) // "fe-mgo-ag100.dvscf_hand_only"
    open( unit=read_unit_2, file=trim(dvfile), iostat=ios, form='unformatted', &
          status='old', action='read', access='direct', recl=rl )
    !
    do imode=1,3*nat
      !
      write(stdout,"(a,i4,a)") "calculate_dvloc_psi: calculate local dv (mode:",imode,")"
      !
      ! Identify the atom and the displacement from the mode
      id = mod(imode-1,3) + 1
      ia = (imode-1)/3 + 1
      !
      ! Identify from which file must be readed the induced change of potential
      read_unit = read_unit_1
      do ia_move=1,nat_move
        if ( ia_move_to_ia(ia_move)==ia ) then
          read_unit = read_unit_2
          exit
        endif
      enddo
      !
      ! Identify the record to be readed
      if (read_unit==read_unit_1) then
        ! Read from mgo-ag100_R441_fine.dvscf_q1
        ! This is a non-magnetic potential
        rec = imode
      elseif (read_unit==read_unit_2) then
        ! Read from fe-mgo-ag100.dvscf_hand_only
        ! This is a magnetic potential (spin up, spin down)
        imode_move = id + 3*(ia_move-1)
        rec = nspin*(imode_move-1) + ispin
      else
        stop "This should not be possible. Stop the program."
      endif
      !
      ! Read the variation of the potential
      read(unit=read_unit,rec=rec,iostat=ios) dv_scf
      !
      call calculate_local_part_dv( qpoint, imode, dv_local )
      !
      dv_local = dv_local + dv_scf
      !
      do iband=1,nbands
        !
        ! dv | psi > in real space
        do ir=1,nr1*nr2*nr3
          aux(ir) = dv_local(ir)*wfc_k_r(ir,iband,1)
        enddo
        !
        ! Fourier transform to G space
        call wfc_from_r_to_g_haritz( igk, aux, dvpsi(:,iband,1,1,imode) )
        !
      enddo
      !
    enddo
    !
    ! Close the file
    close(unit=read_unit)
    !
  end subroutine calculate_dvloc_psi


  subroutine calculate_local_part_dv( qpoint, mode, dvq_local )
    !
    !======================================================================
    ! We have dV_scf as input and we add to it the derivative of the PP   !
    !======================================================================
    !
    use intw_reading , only:  nr1, nr2, nr3, ngm, tau, bg, ityp
    use intw_fft, only: nl
    use intw_reading, only: tpiba
    use intw_useful_constants, only: cmplx_i, tpi
    use intw_fft, only: gvec_cart
    use intw_eqv, only:  vlocq
    use intw_ph, only: eigqts
    !
    implicit none
    !
    !I/O variables
    !
    real(kind=dp), intent(in) :: qpoint(1:3) ! crystal coord.
    integer, intent(in) :: mode
    complex(kind=dp), intent(inout) :: dvq_local(nr1*nr2*nr3)
    !
    !local variables
    !
    integer :: na, ipol, ig, nt, ir
    complex(kind=dp) :: aux(nr1*nr2*nr3), fact, gtau
    real(kind=dp) :: qcart(3) ! qpoint in cart.
    real(kind=dp) :: arg
    !
    qcart = matmul(bg,qpoint)
    !
    !
    na = (mode-1)/3+1
    ipol = modulo(mode-1,3)+1
    nt = ityp (na)
    !
    aux(:) = cmplx_0
    !
    fact = -tpiba*cmplx_i*eigqts(na)
    !
    do ig=1,ngm
      !
      arg = tpi * ( gvec_cart(1,ig)*tau(1,na) + gvec_cart(2,ig)*tau(2,na) + gvec_cart(3,ig)*tau(3,na) ) ! G*tau
      gtau = exp(-cmplx_i*arg)
      !
      aux(nl(ig)) = aux(nl(ig)) + vlocq(ig,nt)*( qcart(ipol) + gvec_cart(ipol,ig) )*gtau*fact
      !
    enddo !ig
    !
    if (gamma_only) then
      do ig=1,ngm
        aux(nlm(ig)) = CONJG(aux(nl(ig)))
      enddo
    end if
    !
    call cfftnd(3,(/nr1,nr2,nr3/),1,aux)
    !
    do ir=1,nr1*nr2*nr3
      !
      dvq_local(ir) = dvq_local(ir) + aux(ir)
      !
    enddo !ir
    !
    return
    !
  end subroutine calculate_local_part_dv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine matrix_element(bra_wfc, dv_ket_wfc, ep)
    !
    !
    !
    implicit none
    !
    ! I/O variables
    complex(kind=dp), intent(in)  :: bra_wfc(:), dv_ket_wfc(:)
    complex(kind=dp), intent(out) :: ep
    !
    ! Local variables
    integer :: ig
    !
    20 format(A)
    !
    ! Check if the dimensions are equal
    if ( size(bra_wfc) == size(dv_ket_wfc) ) then
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
        ep = ep + 2.0_dp*real( conjg(bra_wfc(ig))*dv_ket_wfc(ig) )
      enddo
      ep = ep - conjg(bra_wfc(1))*dv_ket_wfc(1)
    else
      do ig=1,nG_max
        ep = ep + conjg(bra_wfc(ig))*dv_ket_wfc(ig)
      enddo
    endif
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

  subroutine get_GAMMA_eigenvalues()
    !
    ! Reads the eigenvalues from QE .save
    !

    use iotk_module

    implicit none

    !local variables
    character(256) :: K_directory, Kdir
    ! character(256) :: tag
    character(256), dimension(nspin) :: eig_file
    integer :: io_unit,is
    ! integer :: nG
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
    NAMELIST / input_ep_melements / mesh_dir, prefix, ph_dir, dvscf_dir, dvscf_name, &
                                    ep_mat_file, iband, ispin, nat_move, ia_move_to_ia, &
                                    nbands_initial, nbands_final
    !
    ! Set default values
    mesh_dir = '/home/haritz/Kalkuluak/Probak/Fe-O/K/'
    prefix = 'fe-o'
    nk1 = 1
    nk2 = 1
    nk3 = 1
    calc_epmat=.true.
    ph_dir='./'
    dvscf_dir='./dvscf_dir/'
    dvscf_name='fe-o.dvscf_q1'
    ep_mat_file=trim("ep_mat_is_1_iband_1274_VL.dat")
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
    nbands_final = 15
    ! way to calculate the matrix elements
    !
    if (myNode==master) then
      INQUIRE(stdin, NAME=input_file_name)
      if (input_file_name(1:4)=="/dev") then
        ! there is no input file: use default values
        return
      else
        read(stdin,input_ep_melements)
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
    nSize=len(ep_mat_file)
    call MPI_Bcast( ep_mat_file, nSize, MPI_char, master, MPI_COMM_WORLD, MPIerror )
    call MPI_Bcast( iband, 1, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
    call MPI_Bcast( ispin, 1, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
    call MPI_Bcast( nat_move, 1, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
    nSize=size(ia_move_to_ia)
    call MPI_Bcast( ia_move_to_ia, nSize, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
    call MPI_Bcast( nbands_initial, 1, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
    call MPI_Bcast( nbands_final, 1, MPI_integer, master, MPI_COMM_WORLD, MPIerror )
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
    write(stdout,'(a,i2,a,i4,a,i4,a,i4,a)') '|      - Node: ',myNode,' bands: ',MPI_nbands_initial,'-',MPI_nbands_final,' tot: ',num_bands,'            |'
    call MPI_Barrier(MPI_COMM_WORLD, MPIerror)
    !
  end subroutine read_input_file


end program ep_melements
