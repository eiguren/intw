program siesta2dv
  !
  ! This program reads the total potential (V_local + V_hartree + V_xc) of each
  ! irreducible displacement, and calculates the change of potential for the whole
  ! set of atomic displacements.
  !
  use precision, only: dp
  use m_timestamp, only: timestamp
  use fdf, only: fdf_parallel
  use mpi_siesta, only: MPI_Comm_World
  use parallel, only: Node, Nodes
  !
  use siesta2ph_io, only: stdout, outdir, prefix, v0dir, verbose
  !
  use siesta2ph_io, only: read_input
  use siesta2ph_system, only: read_unit_cell_data, print_unit_cell_data
  use siesta2ph_symmetry, only: compute_symmetry, print_symmetry_data, &
                                find_irreducible_atoms, find_irreducible_displacements, &
                                find_site_symmetry
  !
  implicit none
  !
  ! MPI variables
  logical :: initialized
  integer :: MPIerror
  !
  integer, dimension(3) :: mesh
  integer :: nsm, maxp, nspin
  !
  !
  ! Initialise environment
  !
#ifdef MPI
  call MPI_Initialized( initialized, MPIerror )
  if (.not.initialized) then
    call MPI_Init( MPIerror )
  endif ! (.not.initialized)
  call MPI_Comm_Rank( MPI_Comm_World, Node, MPIerror )
  call MPI_Comm_Size( MPI_Comm_World, Nodes, MPIerror )

  if (.not. fdf_parallel()) then
    call die("siesta_init: ERROR: FDF module doesn't have parallel support")
  endif

  if (Nodes>1) then
    call die("siesta2dv: ERROR: This code is not prepared to run in multiple nodes")
  endif
#endif
  !
  Node = 1 ! to turn off output from siesta subroutines
  !
  call timestamp("Start of run")
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Read input file
  call read_input()
  !
  ! Read the supercell data from the fdf created by siesta2ph
  call read_unit_cell_data(trim(outdir)//trim(v0dir)//"/supercell-"//trim(prefix))
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Compute symmetries
  call compute_symmetry()
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Find the irreducible atoms and atom mapping
  call find_irreducible_atoms()
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Find site symmetries for each irreducible atom
  call find_site_symmetry()
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Print all data
  if (verbose) then
    call print_unit_cell_data()
    call print_symmetry_data()
  endif
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Find the irreducible displacements
  call find_irreducible_displacements()
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Read the v0 files to get real space mesh grid
  call get_mesh()
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Write the dvscf.info file
  call write_dvscf_info()
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Read the output files for each irreducible displacements and compute dV
  call compute_dv()
  !
#ifdef MPI
  call MPI_Finalize( MPIerror )
#endif


contains

  subroutine compute_dv()
    !
    ! TODO: Add description
    !
    use precision, only: sp, grid_p
    use m_iorho, only: read_rho
    !
    use siesta2ph_io, only: outdir, v0dir, phdir, dx, nr1, nr2, nr3, dv_precision, lpm, phdir, dvscfdir, &
                            irreducible_q, xsf, full_xsf
    use siesta2ph_system, only: nat, slabel, tau, tau_cryst, at, alat, ityp, atom_labels
    use siesta2ph_symmetry, only: nsite_sym, nsite_rot, irred_atm, irred_disp, site_atm, atom_mapping, site_s_cart, site_sinv, &
                                  nsym, s, ftau, s_cart, rtau_index
    use siesta2ph_utils, only: cmplx_0
    !
    use siesta2ph_io, only: find_free_unit
    use siesta2ph_linalg, only: ainv, rank
    use siesta2ph_symmetry, only: get_q_mesh
    use siesta2ph_utils, only: xsf_datagrid_2d
    !
    implicit none
    !
    ! loop variables
    integer :: ia, ja
    integer :: id
    integer :: isym, isite_sym, iirred
    integer :: ir, i, j, k, Sir, Si, Sj, Sk
    integer :: ispin
    integer :: iq
    !
    ! irreducible displacement variables
    real(kind=dp), dimension(3,3) :: disp
    character(len=256) :: dispp_folder, dispn_folder, disp_filename
    real(kind=dp), dimension(3,3) :: cell
    real(kind=grid_p), allocatable, dimension(:,:) :: rhop
    real(kind=grid_p), allocatable, dimension(:,:) :: rhon
    !
    ! induced potential rotation variables
    real(kind=dp), dimension(3,3) :: basis
    logical, allocatable, dimension(:,:) :: nbasis_irred
    integer :: nbasis, nbasis_count
    real(kind=dp), dimension(3) :: coefx, coefy, coefz
    real(kind=dp), dimension(3) :: R, SR
    !
    ! Cartesian direction rotation variables
    real(kind=dp), dimension(3) :: coefx_cart, coefy_cart, coefz_cart
    !
    ! induced potential variables
    real(kind=dp), allocatable, dimension(:,:) :: dv ! Potential induced by the displacement of an irreducible atom along an irreducible Cartesian direction
    real(kind=dp), allocatable, dimension(:,:,:) :: dv_rot ! Potential induced by the displacement of an irreducible atom along all Cartesian directions
    real(kind=dp), allocatable, dimension(:,:,:) :: dvvq
    complex(kind=dp), allocatable, dimension(:,:,:) :: dvq
    complex(kind=dp), allocatable, dimension(:,:,:) :: dvq_write
    !
    ! unit cell variables
    integer :: nat_uc
    real(kind=dp), dimension(3,3) :: at_uc
    integer, dimension(3) :: mesh_uc
    !
    ! q point variables
    real(kind=dp), dimension(3, nr1*nr2*nr3) :: qmesh_cryst, qirr_cryst
    integer, dimension(nr1*nr2*nr3) :: qirr_map
    integer :: nqirr
    !
    ! dvscf file variables
    character(len=256) :: dvq_file
    integer, allocatable, dimension(:) :: dvq_unit
    character(len=256) :: filename
    integer :: rl
    integer :: ios
    !
    ! xsf file variables
    character(len=256) :: xsf_file
    integer :: xsfiounit
    real(kind=dp), dimension(3) :: x0, e1, e2
    real(kind=dp), allocatable, dimension(:,:) :: rho_xsf
    integer :: kind


    ! Get q-mesh
    qmesh_cryst = 0.0_dp
    qirr_map = 0
    qirr_cryst = 0.0_dp
    nqirr = 0
    !
    call get_q_mesh(nr1, nr2, nr3, qmesh_cryst, qirr_map, irreducible_q, nqirr, qirr_cryst)
    !
    write(stdout,*) "- Computing induced potential..."
    !
    disp(:,1) = (/dx, 0.0_dp, 0.0_dp/)!/alat
    disp(:,2) = (/0.0_dp, dx, 0.0_dp/)!/alat
    disp(:,3) = (/0.0_dp, 0.0_dp, dx/)!/alat
    !
    allocate(dv(maxp,nspin))
    allocate(dv_rot(maxp,nspin,3))
    allocate(dvvq(maxp,nspin,3))
    allocate(dvq(maxp/nr1/nr2/nr3,nspin,3))
    allocate(dvq_write(maxp/nr1/nr2/nr3,nspin,3))
    dv = 0.0_dp
    dv_rot = 0.0_dp
    dvvq = 0.0_dp
    dvq = cmplx_0
    dvq_write = cmplx_0
    !
    allocate(rhop(maxp,nspin))
    allocate(rhon(maxp,nspin))
    rhop = 0.0_grid_p
    rhon = 0.0_grid_p
    !
    ! Record length to write
    if (dv_precision == "sp") then
      inquire(iolength=rl) cmplx(dvq_write(:,:,:), kind=sp)
    else
      inquire(iolength=rl) dvq_write(:,:,:)
    endif
    !
    ! Open the dvscf files
    allocate(dvq_unit(nqirr))
    do iq = 1, nqirr
      !
      write(filename,"(a5,i4.4,a4)") "dvscf", iq, ".dat"
      dvq_unit(iq) = find_free_unit()
      dvq_file = trim(outdir)//trim(phdir)//trim(dvscfdir)//trim(filename)
      open( unit=dvq_unit(iq), file=trim(dvq_file), iostat=ios, form="unformatted", &
            status="replace", action="write", access="direct", recl=rl )
      if (ios .ne. 0) stop "Error opening dvq_file"
      !
    enddo ! iq
    !
    ! Compute the induced potential for each irreducible displacement
    iirred = 0
    do ia=1,nat
      !
      if ( .not. irred_atm(ia) ) cycle
      !
      ! Find site symmetry of the atom
      do isite_sym=1,nsite_sym
        if ( ia .eq. site_atm(isite_sym)) exit
      enddo
      !
      ! Find the rotated displacement vectors
      allocate(nbasis_irred(3,nsite_rot(isite_sym)))
      nbasis_irred = .false.
      basis = 0.0_dp
      nbasis_count = 0
      nbasis = 0
      do id=1,3
        !
        if ( .not. irred_disp(id,ia) ) cycle
        !
        do isym=1,nsite_rot(isite_sym)
          !
          ! Apply the symmetry operation to the displacement vector
          basis(:,nbasis+1) = matmul(transpose(site_s_cart(:,:,isym,isite_sym)), disp(:,id))
          !
          ! Calculate the rank of the basis vector's matrix
          call rank(basis,nbasis)
          !
          ! Save id and symmetry operation used
          if (nbasis > nbasis_count) then
            nbasis_irred(id,isym) = .true.
            nbasis_count = nbasis_count + 1
          endif
          !
          ! If the rank is 3, we have found a complete basis
          if (nbasis == 3) exit
          !
        enddo ! isym
        !
      enddo ! id
      !
      ! Find the linear combination of the basis vectors to obtain the cartesian displacements
      coefx = matmul((ainv(basis)), disp(:,1))
      coefy = matmul((ainv(basis)), disp(:,2))
      coefz = matmul((ainv(basis)), disp(:,3))
      !
      ! Calculate dV for atom ia
      nbasis_count = 0
      dv_rot = 0.0_dp
      do id=1,3
        !
        if ( .not. irred_disp(id,ia) ) cycle
        !
        iirred = iirred + 1
        write(stdout,"(a,i3,2(a5,i3))") "irreducible displacement", iirred, ": ia=", ia, ", id=", id
        !
        ! Read total potential for the irreducible displacement
        Node = 0 ! required by read_rho
        if (lpm) then
          !
          write(dispp_folder,"(a5,i4.4,a1)") "disp-", 2*iirred-1, "/"
          disp_filename = trim(outdir)//trim(phdir)//trim(dispp_folder)//trim(slabel)//".VT"
          call read_rho( disp_filename, cell, mesh, nsm, maxp, nspin, rhop)
          !
          write(dispn_folder,"(a5,i4.4,a1)") "disp-", 2*iirred, "/"
          disp_filename = trim(outdir)//trim(phdir)//trim(dispn_folder)//trim(slabel)//".VT"
          call read_rho( disp_filename, cell, mesh, nsm, maxp, nspin, rhon)
          !
        else
          !
          write(dispp_folder,"(a5,i4.4,a1)") "disp-", iirred, "/"
          disp_filename = trim(outdir)//trim(phdir)//trim(dispp_folder)//trim(slabel)//".VT"
          call read_rho( disp_filename, cell, mesh, nsm, maxp, nspin, rhop)
          !
          disp_filename = trim(outdir)//trim(v0dir)//trim(slabel)//".VT"
          call read_rho( disp_filename, cell, mesh, nsm, maxp, nspin, rhon)
          !
        endif
        Node = 1 ! turn off again all siesta output
        !
        ! Transform to Pauli matrix representation on the NC case (as in QE)
        if (nspin>=2) then
          write(stdout,"(a)") "Pauli matrices representation"
          !
          call rho_from_spin_components_2_pauli_matrices(rhop)
          call rho_from_spin_components_2_pauli_matrices(rhon)
          !
        endif
        !
        ! Compute dV
        do ispin=1,nspin
          if (lpm) then
            dv(:,ispin) = (rhop(:, ispin) - rhon(:, ispin)) / (2.0_dp*dx)
          else
            dv(:,ispin) = (rhop(:, ispin) - rhon(:, ispin)) / dx
          endif
        enddo
        !
        ! Calculate dV for the irreduble displacement and rotate to three cartesian directions
        do isym=1,nsite_rot(isite_sym)
          !
          if (.not.nbasis_irred(id,isym)) cycle
          !
          nbasis_count = nbasis_count + 1
          !
          ! Rotate to cartesian directions
          do k=1,mesh(3)
            do j=1,mesh(2)
              do i=1,mesh(1)
                !
                ! r on crystal coordinates
                R(1) = real(i-1,dp)/mesh(1)
                R(2) = real(j-1,dp)/mesh(2)
                R(3) = real(k-1,dp)/mesh(3)
                !
                ! Rotate and translate the r vector
                SR = matmul(site_sinv(:,:,isym,isite_sym), R-tau_cryst(:,ia)) + tau_cryst(:,ia)
                !
                ! Find the (i,j,k) indices of the rotated r
                Si = mod(nint(mod(SR(1)+1.0_dp,1.0_dp)*mesh(1)), mesh(1)) + 1
                Sj = mod(nint(mod(SR(2)+1.0_dp,1.0_dp)*mesh(2)), mesh(2)) + 1
                Sk = mod(nint(mod(SR(3)+1.0_dp,1.0_dp)*mesh(3)), mesh(3)) + 1
                ir = i + mesh(1)*(j-1) + mesh(1)*mesh(2)*(k-1)
                Sir = Si + mesh(1)*(Sj-1) + mesh(1)*mesh(2)*(Sk-1)
                !
                ! Linear combination for each cartesian direction
                do ispin=1,nspin
                  dv_rot(Sir,ispin,1) = dv_rot(Sir,ispin,1) + coefx(nbasis_count)*dv(ir,ispin)
                  dv_rot(Sir,ispin,2) = dv_rot(Sir,ispin,2) + coefy(nbasis_count)*dv(ir,ispin)
                  dv_rot(Sir,ispin,3) = dv_rot(Sir,ispin,3) + coefz(nbasis_count)*dv(ir,ispin)
                enddo
                !
              enddo
            enddo
          enddo
          !
        enddo ! isym
        !
      enddo ! id
      !
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Compute and save dvq to the dvscf file
      !
      write(stdout,*) "- Saving dvscf..."
      !
      nat_uc = nat/nr1/nr2/nr3
      at_uc(:,1) = at(:,1)/nr1
      at_uc(:,2) = at(:,2)/nr2
      at_uc(:,3) = at(:,3)/nr3
      mesh_uc = mesh / (/nr1, nr2, nr3/)
      !
      do iq = 1, nqirr
        !
        ! q point in the unit cell crytal basis
        write(stdout,"(a11,3f10.6)") "  q_cryst: ", qirr_cryst(:, iq)
        !
        call dVR_to_dVq(ia, dv_rot, qirr_cryst(:, iq), dvq)
        !
        ! Rotate to all equivalent unit cell atoms
        do ja=1,nat_uc
          !
          if (atom_mapping(ja) .ne. ia) cycle
          !
          ! Find symmetry operation
          do isym=1,nsym
            if (rtau_index(ia,isym) == ja) exit
            ! if (all(mod(abs(rtau_cryst(:,isym,ia) - tau_cryst(:,ja)), 1.0_dp) < 0.000001_dp)) exit
          enddo
          if (isym == nsym .and. rtau_index(ia,isym) /= ja) stop "ERROR: compute_dv: symmetry operation not found"
          !
          ! Find the coefficients to recover Cartesian directions after rotation
          coefx_cart = matmul((ainv(s_cart(:,:,isym))), (/1.0_dp, 0.0_dp, 0.0_dp/))
          coefy_cart = matmul((ainv(s_cart(:,:,isym))), (/0.0_dp, 1.0_dp, 0.0_dp/))
          coefz_cart = matmul((ainv(s_cart(:,:,isym))), (/0.0_dp, 0.0_dp, 1.0_dp/))
          !
          ! Rotate potential
          do k=1,mesh_uc(3)
            do j=1,mesh_uc(2)
              do i=1,mesh_uc(1)
                !
                ! r on crystal coordinates
                R(1) = real(i-1,dp)/mesh_uc(1)
                R(2) = real(j-1,dp)/mesh_uc(2)
                R(3) = real(k-1,dp)/mesh_uc(3)
                !
                ! Rotate and translate the r vector
                SR = matmul(s(:,:,isym), R) - ftau(:,isym)
                !
                ! Find the (i,j,k) indices of the rotated r
                Si = mod(nint(mod(SR(1)+1.0_dp,1.0_dp)*mesh(1)), mesh(1)) + 1
                Sj = mod(nint(mod(SR(2)+1.0_dp,1.0_dp)*mesh(2)), mesh(2)) + 1
                Sk = mod(nint(mod(SR(3)+1.0_dp,1.0_dp)*mesh(3)), mesh(3)) + 1
                ir = i + mesh(1)*(j-1) + mesh(1)*mesh(2)*(k-1)
                Sir = Si + mesh(1)*(Sj-1) + mesh(1)*mesh(2)*(Sk-1)
                !
                dvq_write(Sir,:,1) = coefx_cart(1)*dvq(ir,:,1) + coefx_cart(2)*dvq(ir,:,2) + coefx_cart(3)*dvq(ir,:,3)
                dvq_write(Sir,:,2) = coefy_cart(1)*dvq(ir,:,1) + coefy_cart(2)*dvq(ir,:,2) + coefy_cart(3)*dvq(ir,:,3)
                dvq_write(Sir,:,3) = coefz_cart(1)*dvq(ir,:,1) + coefz_cart(2)*dvq(ir,:,2) + coefz_cart(3)*dvq(ir,:,3)
                !
              enddo
            enddo
          enddo
          !
          ! Write to dvscf file
          if (dv_precision=="sp") write(unit=dvq_unit(iq),rec=ja,iostat=ios) cmplx(dvq_write, kind=sp)
          if (dv_precision=="dp") write(unit=dvq_unit(iq),rec=ja,iostat=ios) dvq_write
          if (ios .ne. 0) stop "Error writing dv"
          !
        enddo ! ja
        !
      enddo ! iq
      !
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Print the potential to a xsf file
      if (xsf) then
        !
        allocate(rho_xsf(mesh(1),mesh(2)))
        !
        write(stdout,*) "- Writing dVR xsf files..."
        !
        do ja=1,nat
          !
          if (full_xsf) then
            if (atom_mapping(ja) .ne. ia) cycle
          else
            if (ja.ne.ia) cycle
          endif
          !
          ! Find symmetry operation
          do isym=1,nsym
            if (rtau_index(ia,isym) == ja) exit
            ! if (all(mod(abs(rtau_cryst(:,isym,ia) - tau_cryst(:,ja)), 1.0_dp) < 0.000001_dp)) exit
          enddo
          if (isym == nsym .and. rtau_index(ia,isym) /= ja) stop "ERROR: compute_dv: symmetry operation not found"
          !
          ! Find the coefficients to recover Cartesian directions after rotation
          coefx_cart = matmul((ainv(s_cart(:,:,isym))), (/1.0_dp, 0.0_dp, 0.0_dp/))
          coefy_cart = matmul((ainv(s_cart(:,:,isym))), (/0.0_dp, 1.0_dp, 0.0_dp/))
          coefz_cart = matmul((ainv(s_cart(:,:,isym))), (/0.0_dp, 0.0_dp, 1.0_dp/))
          !
          do ispin=1,nspin
            !
            do id=1,3
              !
              if (nspin==1) then
                write(xsf_file,"(a6,i3.3,a3,i1,a4)") "dVR_ia", ja, "_id", id, ".xsf"
              else
                write(xsf_file,"(a6,i3.3,a3,i1,a3,i1,a4)") "dVR_ia", ja, "_id", id, "_is", ispin, ".xsf"
              endif
              !
              xsfiounit = find_free_unit()
              !
              x0 = modulo(tau_cryst(:,ja), 1.0_dp)
              x0(1:2) = 0.0_dp
              !
              kind = nint( x0(3)*mesh(3) ) + 1
              !
              x0 = matmul(at, x0) ! we need x0 in cartesian coordinates
              !
              e1 = at(:,1)
              e2 = at(:,2)
              !
              rho_xsf = 0.0_dp
              do k=1,mesh(3)
                do j=1,mesh(2)
                  do i=1,mesh(1)
                    !
                    ! r on crystal coordinates
                    R(1) = real(i-1,dp)/mesh(1)
                    R(2) = real(j-1,dp)/mesh(2)
                    R(3) = real(k-1,dp)/mesh(3)
                    !
                    ! Rotate and translate the r vector
                    SR = matmul(s(:,:,isym), R) - ftau(:,isym)
                    !
                    ! Find the (i,j,k) indices of the rotated r
                    Si = mod(nint(mod(SR(1)+1.0_dp,1.0_dp)*mesh(1)), mesh(1)) + 1
                    Sj = mod(nint(mod(SR(2)+1.0_dp,1.0_dp)*mesh(2)), mesh(2)) + 1
                    Sk = mod(nint(mod(SR(3)+1.0_dp,1.0_dp)*mesh(3)), mesh(3)) + 1
                    ir = i + mesh(1)*(j-1) + mesh(1)*mesh(2)*(k-1)
                    Sir = Si + mesh(1)*(Sj-1) + mesh(1)*mesh(2)*(Sk-1)
                    !
                    if (Sk.ne.kind) cycle
                    !
                    if (id == 1) rho_xsf(Si,Sj) = coefx_cart(1)*dv_rot(ir,ispin,1) + coefx_cart(2)*dv_rot(ir,ispin,2) + coefx_cart(3)*dv_rot(ir,ispin,3)
                    if (id == 2) rho_xsf(Si,Sj) = coefy_cart(1)*dv_rot(ir,ispin,1) + coefy_cart(2)*dv_rot(ir,ispin,2) + coefy_cart(3)*dv_rot(ir,ispin,3)
                    if (id == 3) rho_xsf(Si,Sj) = coefz_cart(1)*dv_rot(ir,ispin,1) + coefz_cart(2)*dv_rot(ir,ispin,2) + coefz_cart(3)*dv_rot(ir,ispin,3)
                    !
                  enddo
                enddo
              enddo
              !
              call xsf_datagrid_2d(rho_xsf, x0, e1, e2, mesh(1), mesh(2), at, alat, nat, tau, atom_labels, ityp, xsfiounit, trim(outdir)//trim(phdir)//trim(xsf_file))
              !
            enddo ! id
            !
          enddo ! ispin
          !
        enddo ! ja
        !
        write(stdout,*) "- Writing dVq xsf files..."
        !
        do iq = 1, nqirr
          !
          call dVR_to_dVq_sc(ia, dv_rot, qirr_cryst(:, iq), dvvq)
          !
          do ja=1,nat_uc
            !
            if (full_xsf) then
              if (atom_mapping(ja) .ne. ia) cycle
            else
              if (ja.ne.ia) cycle
            endif
            !
            ! Find symmetry operation
            do isym=1,nsym
              if (rtau_index(ia,isym) == ja) exit
              ! if (all(mod(abs(rtau_cryst(:,isym,ia) - tau_cryst(:,ja)), 1.0_dp) < 0.000001_dp)) exit
            enddo
            if (isym == nsym .and. rtau_index(ia,isym) /= ja) stop "ERROR: compute_dv: symmetry operation not found"
            !
            ! Find the coefficients to recover Cartesian directions after rotation
            coefx_cart = matmul((ainv(s_cart(:,:,isym))), (/1.0_dp, 0.0_dp, 0.0_dp/))
            coefy_cart = matmul((ainv(s_cart(:,:,isym))), (/0.0_dp, 1.0_dp, 0.0_dp/))
            coefz_cart = matmul((ainv(s_cart(:,:,isym))), (/0.0_dp, 0.0_dp, 1.0_dp/))
            !
            do ispin=1,nspin
              !
              do id=1,3
                !
                if (nspin==1) then
                  write(xsf_file,"(a6,i2.2,a3,i3.3,a3,i1,a4)") "dVq_iq", iq, "_ia", ja, "_id", id, ".xsf"
                else
                  write(xsf_file,"(a6,i2.2,a3,i3.3,a3,i1,a3,i1,a4)") "dVq_iq", iq, "_ia", ja, "_id", id, "_is", ispin, ".xsf"
                endif
                !
                xsfiounit = find_free_unit()
                !
                x0 = modulo(tau_cryst(:,ja), 1.0_dp)
                x0(1:2) = 0.0_dp
                !
                kind = nint( x0(3)*mesh(3) ) + 1
                !
                x0 = matmul(at, x0) ! we need x0 in cartesian coordinates
                !
                e1 = at(:,1)
                e2 = at(:,2)
                !
                rho_xsf = 0.0_dp
                do k=1,mesh(3)
                  do j=1,mesh(2)
                    do i=1,mesh(1)
                      !
                      ! r on crystal coordinates
                      R(1) = real(i-1,dp)/mesh(1)
                      R(2) = real(j-1,dp)/mesh(2)
                      R(3) = real(k-1,dp)/mesh(3)
                      !
                      ! Rotate and translate the r vector
                      SR = matmul(s(:,:,isym), R) - ftau(:,isym)
                      !
                      ! Find the (i,j,k) indices of the rotated r
                      Si = mod(nint(mod(SR(1)+1.0_dp,1.0_dp)*mesh(1)), mesh(1)) + 1
                      Sj = mod(nint(mod(SR(2)+1.0_dp,1.0_dp)*mesh(2)), mesh(2)) + 1
                      Sk = mod(nint(mod(SR(3)+1.0_dp,1.0_dp)*mesh(3)), mesh(3)) + 1
                      ir = i + mesh(1)*(j-1) + mesh(1)*mesh(2)*(k-1)
                      Sir = Si + mesh(1)*(Sj-1) + mesh(1)*mesh(2)*(Sk-1)
                      !
                      if (Sk.ne.kind) cycle
                      !
                      if (id == 1) rho_xsf(Si,Sj) = coefx_cart(1)*dvvq(ir,ispin,1) + coefx_cart(2)*dvvq(ir,ispin,2) + coefx_cart(3)*dvvq(ir,ispin,3)
                      if (id == 2) rho_xsf(Si,Sj) = coefy_cart(1)*dvvq(ir,ispin,1) + coefy_cart(2)*dvvq(ir,ispin,2) + coefy_cart(3)*dvvq(ir,ispin,3)
                      if (id == 3) rho_xsf(Si,Sj) = coefz_cart(1)*dvvq(ir,ispin,1) + coefz_cart(2)*dvvq(ir,ispin,2) + coefz_cart(3)*dvvq(ir,ispin,3)
                      !
                    enddo
                  enddo
                enddo
                !
                call xsf_datagrid_2d(rho_xsf, x0, e1, e2, mesh(1), mesh(2), at, alat, nat, tau, atom_labels, ityp, xsfiounit, trim(outdir)//trim(phdir)//trim(xsf_file))
                !
              enddo ! id
              !
            enddo ! ispin
            !
          enddo ! ja
          !
        enddo ! iq
        !
        write(stdout,*) "- Writing dvq xsf files..."
        !
        do iq = 1, nqirr
          !
          call dVR_to_dVq(ia, dv_rot, qirr_cryst(:, iq), dvq)
          !
          do ja=1,nat_uc
            !
            if (full_xsf) then
              if (atom_mapping(ja) .ne. ia) cycle
            else
              if (ja.ne.ia) cycle
            endif
            !
            ! Find symmetry operation
            do isym=1,nsym
              if (rtau_index(ia,isym) == ja) exit
              ! if (all(mod(abs(rtau_cryst(:,isym,ia) - tau_cryst(:,ja)), 1.0_dp) < 0.000001_dp)) exit
            enddo
            if (isym == nsym .and. rtau_index(ia,isym) /= ja) stop "ERROR: compute_dv: symmetry operation not found"
            !
            ! Find the coefficients to recover Cartesian directions after rotation
            coefx_cart = matmul((ainv(s_cart(:,:,isym))), (/1.0_dp, 0.0_dp, 0.0_dp/))
            coefy_cart = matmul((ainv(s_cart(:,:,isym))), (/0.0_dp, 1.0_dp, 0.0_dp/))
            coefz_cart = matmul((ainv(s_cart(:,:,isym))), (/0.0_dp, 0.0_dp, 1.0_dp/))
            !
            do ispin=1,nspin
              !
              do id=1,3
                !
                if (nspin==1) then
                  write(xsf_file,"(a6,i2.2,a3,i3.3,a3,i1,a4)") "dvq_iq", iq, "_ia", ja, "_id", id, ".xsf"
                else
                  write(xsf_file,"(a6,i2.2,a3,i3.3,a3,i1,a3,i1,a4)") "dvq_iq", iq, "_ia", ja, "_id", id, "_is", ispin, ".xsf"
                endif
                !
                xsfiounit = find_free_unit()
                !
                x0 = modulo(tau_cryst(:,ja), 1.0_dp)
                x0(1:2) = 0.0_dp
                !
                kind = nint( x0(3)*mesh(3) ) + 1
                !
                x0 = matmul(at, x0) ! we need x0 in cartesian coordinates
                !
                e1 = at(:,1)
                e2 = at(:,2)
                !
                rho_xsf = 0.0_dp
                do k=1,mesh_uc(3)
                  do j=1,mesh_uc(2)
                    do i=1,mesh_uc(1)
                      !
                      ! r on crystal coordinates
                      R(1) = real(i-1,dp)/mesh_uc(1)
                      R(2) = real(j-1,dp)/mesh_uc(2)
                      R(3) = real(k-1,dp)/mesh_uc(3)
                      !
                      ! Rotate and translate the r vector
                      SR = matmul(s(:,:,isym), R) - ftau(:,isym)
                      !
                      ! Find the (i,j,k) indices of the rotated r
                      Si = mod(nint(mod(SR(1)+1.0_dp,1.0_dp)*mesh(1)), mesh(1)) + 1
                      Sj = mod(nint(mod(SR(2)+1.0_dp,1.0_dp)*mesh(2)), mesh(2)) + 1
                      Sk = mod(nint(mod(SR(3)+1.0_dp,1.0_dp)*mesh(3)), mesh(3)) + 1
                      ir = i + mesh(1)*(j-1) + mesh(1)*mesh(2)*(k-1)
                      Sir = Si + mesh(1)*(Sj-1) + mesh(1)*mesh(2)*(Sk-1)
                      !
                      if (Sk.ne.kind) cycle
                      !
                      if (id == 1) rho_xsf(Si,Sj) = coefx_cart(1)*real(dvq(ir,ispin,1)) + coefx_cart(2)*real(dvq(ir,ispin,2)) + coefx_cart(3)*real(dvq(ir,ispin,3))
                      if (id == 2) rho_xsf(Si,Sj) = coefy_cart(1)*real(dvq(ir,ispin,1)) + coefy_cart(2)*real(dvq(ir,ispin,2)) + coefy_cart(3)*real(dvq(ir,ispin,3))
                      if (id == 3) rho_xsf(Si,Sj) = coefz_cart(1)*real(dvq(ir,ispin,1)) + coefz_cart(2)*real(dvq(ir,ispin,2)) + coefz_cart(3)*real(dvq(ir,ispin,3))
                      !
                    enddo
                  enddo
                enddo
                !
                call xsf_datagrid_2d(rho_xsf, x0, e1, e2, mesh(1), mesh(2), at, alat, nat, tau, atom_labels, ityp, xsfiounit, trim(outdir)//trim(phdir)//trim(xsf_file))
                !
              enddo ! id
              !
            enddo ! ispin
            !
          enddo ! ja
          !
        enddo ! iq
        !
        deallocate(rho_xsf)
        !
      endif
      !
      deallocate(nbasis_irred)
      !
    enddo ! ia
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Close the dvscf file
    close(unit=dvq_unit(1))

  end subroutine compute_dv


  subroutine write_dvscf_info()
    !
    ! TODO: Add description
    !
    use siesta2ph_io, only: outdir, phdir, dvscfdir, dv_precision
    use siesta2ph_system, only: nat
    use siesta2ph_symmetry, only: irred_atm, nsite_sym, atom_mapping
    !
    use siesta2ph_io, only: find_free_unit
    !
    implicit none
    !
    character(len=256) :: filename
    integer :: iounit, iostat
    integer :: ia


    write(stdout,*) "- Writing the dvscf.info file..."
    !
    ! Open the dvscf.info file
    iounit = find_free_unit()
    filename = trim(outdir)//trim(phdir)//trim(dvscfdir)//"dvscf.info"
    call system("mkdir -p "//trim(outdir)//trim(phdir)//trim(dvscfdir))
    open( unit=iounit, file=trim(filename), iostat=iostat, form="formatted", &
          status="replace", action="write" )
    if (iostat .ne. 0) stop "ERROR: write_dvscf_info: Error opening file."
    !
    ! Write the data
    write(iounit,*) "DV_PRECISION"
    write(iounit,*) dv_precision ! sp or dp for the potential
    !
    write(iounit,*) "NAT"
    write(iounit,*) nat ! Total number of atoms on the supercell TODO: should we save the number of atoms of the iounit cell?
    !
    write(iounit,*) "NIRRED_ATM"
    write(iounit,*) nsite_sym ! This is equal to the number of irreducible atoms
    !
    write(iounit,*) "NSPIN"
    write(iounit,*) nspin ! Number of spin components of the potential
    !
    write(iounit,*) "IRRED_ATM(1:NIRRED_ATM)"
    do ia=1,nat
      if (irred_atm(ia)) write(iounit,*) ia
    enddo
    !
    write(iounit,*) "ATM_MAPPING(1:NAT)"
    do ia=1,nat
      write(iounit,*) atom_mapping(ia)
    enddo
    !
    close(unit=iounit)

  end subroutine write_dvscf_info


  subroutine get_mesh()
    !
    ! TODO: Add description
    !
    use siesta2ph_io, only: outdir, v0dir, verbose
    use siesta2ph_system, only: at, alat, slabel
    !
    use siesta2ph_io, only: find_free_unit
    !
    implicit none
    !
    character(len=256) :: filename
    logical :: found
    integer :: iounit
    real(kind=dp), dimension(3,3) :: cell


    write(stdout,*) "- Reading grid dimensions..."
    !
    filename = trim(outdir)//trim(v0dir)//trim(slabel)//".VT"
    !
    ! Check if file exists
    inquire(file=filename, exist=found )
    if (.not. found) stop "ERROR: Cannot find v0dir/SystemLavel.VT file"
    !
    ! Open file
    iounit = find_free_unit()
    open(iounit, file=filename, form="unformatted", status="old")
    !
    ! Read cell vectors and mesh grid
    read(iounit) cell ! in Bohr
    read(iounit) mesh, nspin
    !
    ! Close file
    close(iounit)
    !
    ! Some checks
    if (sum(abs(cell-at*alat)) > 1.0d-8) stop "ERROR: the unit cell from v0dir/SystemLavel.VT doesn't correspond with the supercell"
    !
    maxp = mesh(1)*mesh(2)*mesh(3)
    nsm = 1
    !
    if (verbose) then
      write(stdout,"(a26,3i6)") "FFT grid: nr1, nr2, nr3 = ", mesh
      write(stdout,"(a24,i10)") "FFT grid: nr1*nr2*nr3 = ", maxp
    endif

  end subroutine get_mesh


  subroutine rho_from_spin_components_2_pauli_matrices(rho)
    !
    ! Transform the non-colinear 4x4 potentialn from the spin components
    ! to the pauli matrices representation
    !
    use precision, only: grid_p
    !
    implicit none
    !
    real(kind=grid_p), dimension(:,:), intent(inout) :: rho
    !
    real(kind=grid_p) :: average, difference
    integer :: nr, ns
    integer :: ir


    ! Check input/output array dimensions
    nr = size(rho, dim=1)
    ns = size(rho, dim=2)
    !
    if (ns/=2 .and. ns/=4) stop "ERROR: rho_from_spin_components_2_pauli_matrices: Invalid ns"
    !
    do ir=1,nr
      !
      ! Average of the diagonal components
      average = (rho(ir,1) + rho(ir,2)) / 2.0_grid_p
      ! Difference of the diagonal components
      difference = (rho(ir,1) - rho(ir,2)) / 2.0_grid_p
      !
      if (ns == 2) then
        rho(ir,1) = average ! I
        rho(ir,2) = difference ! sig_z
      else
        rho(ir,1) = average ! I
        rho(ir,2) = rho(ir,3) ! sig_x
        rho(ir,3) = rho(ir,4) ! sig_y
        rho(ir,4) = difference ! sig_z
      endif
      !
    enddo

  end subroutine rho_from_spin_components_2_pauli_matrices


  subroutine dVR_to_dVq(ia, dvv, q_cryst, dvq)
    !
    ! Compute the periodic part of the induced potential dv_q:
    ! dv_q(r) = e^{-iqr} dV_q(r) = e^{-iqr} \sum_{R} dV_R(r) e^{iqR}
    !
    use siesta2ph_io, only: nr1, nr2, nr3
    use siesta2ph_system, only: nat, at, tau_cryst
    use siesta2ph_utils, only: tpi, cmplx_0, cmplx_i
    !
    use siesta2ph_utils, only: eqvect
    !
    implicit none
    !
    integer, intent(in) :: ia
    real(kind=dp), dimension(:,:,:), intent(in) :: dvv ! dV_R(r): Potential induced by the displacement of atom ia (u_ia) in the supercell (r in sc)
    real(kind=dp), dimension(3), intent(in) :: q_cryst ! in unit cell crystal coordinates
    complex(kind=dp), dimension(:,:,:), intent(out) :: dvq ! dv_q(r): Periodic part of the potential induced by the periodic displacement of atom ia (u_ia*e^{iqR}) in the unit cell (r in uc)
    !
    real(kind=dp), dimension(3,27) :: vecs_cryst
    real(kind=dp), dimension(27) :: vecs_lenght
    real(kind=dp), dimension(3) :: R_cryst, vec_cart, r, Sr
    complex(kind=dp) :: phase
    real(kind=dp) :: min_lenght
    integer :: ja
    integer :: ir1, ir2, ir3, ir_uc
    integer :: Sir1, Sir2, Sir3, Sir
    integer :: ivec, n_vecs
    integer, dimension(3) :: mesh_uc


    ! Check input/output array dimensions
    if (size(dvv, dim=1) /= mesh(1)*mesh(2)*mesh(3)) stop "ERROR: dVR_to_dVq: wrong mesh dimensions in input"
    if (size(dvv, dim=2) /= nspin) stop "ERROR: dVR_to_dVq: wrong spin dimensions in input"
    if (size(dvv, dim=3) /= 3) stop "ERROR: dVR_to_dVq: wrong spatial dimensions in input"
    !
    mesh_uc = mesh / (/nr1, nr2, nr3/)
    if (size(dvq, dim=1) /= mesh_uc(1)*mesh_uc(2)*mesh_uc(3)) stop "ERROR: dVR_to_dVq: wrong mesh dimensions in output"
    if (size(dvq, dim=2) /= nspin) stop "ERROR: dVR_to_dVq: wrong spin dimensions in output"
    if (size(dvq, dim=3) /= 3) stop "ERROR: dVR_to_dVq: wrong spatial dimensions in output"
    !
    dvq = cmplx_0
    do ja=1,nat
      !
      if ( .not. eqvect(tau_cryst(:,ia)*(/nr1, nr2, nr3/), tau_cryst(:,ja)*(/nr1, nr2, nr3/)) ) cycle ! i.e. if the atoms are equivalent by translation
      !
      ivec = 0
      do ir1=-1,1
        do ir2=-1,1
          do ir3=-1,1
            !
            ivec = ivec + 1
            !
            R_cryst = (/ir1, ir2, ir3/) ! lattice vector of the supercell in supercell crystal coordinates
            !
            vecs_cryst(:, ivec) = ((tau_cryst(:,ja) + R_cryst) - tau_cryst(:,ia)) ! vector from tau(ia) to tau(ja)+R in supercell crystal coordinates
            vec_cart = matmul(at, vecs_cryst(:, ivec))
            vecs_lenght(ivec) = sqrt(sum(vec_cart**2))
            vecs_cryst(:, ivec) = vecs_cryst(:, ivec)*(/nr1, nr2, nr3/) ! transform to the unit cell crystal coordinates
            !
          enddo
        enddo
      enddo
      !
      min_lenght = minval(vecs_lenght)
      !
      n_vecs = 0
      phase = cmplx_0
      do ivec=1,27
        if ( abs( vecs_lenght(ivec) - min_lenght ) < 1.0d-5 ) then
          n_vecs = n_vecs + 1
          phase = phase + exp(tpi*cmplx_i*dot_product(q_cryst, vecs_cryst(:, ivec)))
        endif
      enddo
      !
      do ir3=1,mesh_uc(3)
        do ir2=1,mesh_uc(2)
          do ir1=1,mesh_uc(1)
            !
            ir_uc = ir1 + mesh_uc(1)*(ir2-1) + mesh_uc(1)*mesh_uc(2)*(ir3-1)
            !
            ! r on supercell crystal coordinates
            r(1) = real(ir1-1,dp)/mesh(1)
            r(2) = real(ir2-1,dp)/mesh(2)
            r(3) = real(ir3-1,dp)/mesh(3)
            !
            ! Translate the r vector
            Sr = r - (tau_cryst(:,ja) - tau_cryst(:,ia))
            !
            ! Find the supercell (i,j,k) indices of the translated r
            Sir1 = nint(mod(Sr(1)+1.0_dp,1.0_dp)*mesh(1)) + 1
            Sir2 = nint(mod(Sr(2)+1.0_dp,1.0_dp)*mesh(2)) + 1
            Sir3 = nint(mod(Sr(3)+1.0_dp,1.0_dp)*mesh(3)) + 1
            !
            ! Find the global index
            Sir = Sir1 + mesh(1)*(Sir2-1) + mesh(1)*mesh(2)*(Sir3-1)
            !
            ! r on unit cell crystal coordinates
            r(1) = real(ir1-1,dp)/mesh_uc(1)
            r(2) = real(ir2-1,dp)/mesh_uc(2)
            r(3) = real(ir3-1,dp)/mesh_uc(3)
            !
            dvq(ir_uc,:,:) = dvq(ir_uc,:,:) + dvv(Sir,:,:) * phase / n_vecs * exp(-cmplx_i*tpi*dot_product(q_cryst, r))
            !
          enddo
        enddo
      enddo
      !
    enddo ! ja

  end subroutine dVR_to_dVq


  subroutine dVR_to_dVq_sc(ia, dvv, q_cryst, dvvq)
    !
    ! Compute the super-cell periodic induced potential dV_q:
    ! dV_q(r) = e^{iqr} dv_q(r) = \sum_{R} dV_R(r) e^{iqR}
    !

    use siesta2ph_io, only: nr1, nr2, nr3
    use siesta2ph_system, only: nat, at, tau_cryst
    use siesta2ph_utils, only: tpi, cmplx_0, cmplx_i
    !
    use siesta2ph_utils, only: eqvect
    !
    implicit none
    !
    integer, intent(in) :: ia
    real(kind=dp), dimension(:,:,:), intent(in) :: dvv ! dV_R(r): Potential induced by the displacement of atom ia (u_ia) in the supercell (r in sc)
    real(kind=dp), dimension(3), intent(in) :: q_cryst ! in unit cell crystal coordinates
    real(kind=dp), dimension(:,:,:), intent(out) :: dvvq ! dV_q: Potential induced by the periodic displacement of atom ia (u_ia*e^{iqR}) in the supercell cell (r in sc)
    !
    real(kind=dp), dimension(3,27) :: vecs_cryst
    real(kind=dp), dimension(27) :: vecs_lenght
    real(kind=dp), dimension(3) :: R_cryst, vec_cart, r, Sr
    complex(kind=dp) :: phase
    real(kind=dp) :: min_lenght
    integer :: ja
    integer :: ir1, ir2, ir3, ir
    integer :: Sir1, Sir2, Sir3, Sir
    integer :: ivec, n_vecs


    ! Check input/output array dimensions
    if (size(dvv, dim=1) /= mesh(1)*mesh(2)*mesh(3)) stop "ERROR: dVR_to_dVq_sc: wrong mesh dimensions in input"
    if (size(dvv, dim=2) /= nspin) stop "ERROR: dVR_to_dVq_sc: wrong spin dimensions in input"
    if (size(dvv, dim=3) /= 3) stop "ERROR: dVR_to_dVq_sc: wrong spatial dimensions in input"
    !
    if (size(dvvq, dim=1) /= mesh(1)*mesh(2)*mesh(3)) stop "ERROR: dVR_to_dVq_sc: wrong mesh dimensions in output"
    if (size(dvvq, dim=2) /= nspin) stop "ERROR: dVR_to_dVq_sc: wrong spin dimensions in output"
    if (size(dvvq, dim=3) /= 3) stop "ERROR: dVR_to_dVq_sc: wrong spatial dimensions in output"
    !
    dvvq = 0.0_dp
    do ja=1,nat
      !
      if ( .not. eqvect(tau_cryst(:,ia)*(/nr1, nr2, nr3/), tau_cryst(:,ja)*(/nr1, nr2, nr3/)) ) cycle ! i.e. if the atoms are equivalent by translation
      !
      ivec = 0
      do ir1=-1,1
        do ir2=-1,1
          do ir3=-1,1
            !
            ivec = ivec + 1
            !
            R_cryst = (/ir1, ir2, ir3/) ! lattice vector of the supercell in supercell crystal coordinates
            !
            vecs_cryst(:, ivec) = ((tau_cryst(:,ja) + R_cryst) - tau_cryst(:,ia)) ! vector from tau(ia) to tau(ja)+R in supercell crystal coordinates
            vec_cart = matmul(at, vecs_cryst(:, ivec))
            vecs_lenght(ivec) = sqrt(sum(vec_cart**2))
            vecs_cryst(:, ivec) = vecs_cryst(:, ivec)*(/nr1, nr2, nr3/) ! transform to th unit cell crystal coordinates
            !
          enddo
        enddo
      enddo
      !
      min_lenght = minval(vecs_lenght)
      !
      n_vecs = 0
      phase = cmplx_0
      do ivec=1,27
        if ( abs( vecs_lenght(ivec) - min_lenght ) < 1.0d-5 ) then
          n_vecs = n_vecs + 1
          phase = phase + exp(tpi*cmplx_i*dot_product(q_cryst, vecs_cryst(:, ivec)))
        endif
      enddo
      !
      do ir3=1,mesh(3)
        do ir2=1,mesh(2)
          do ir1=1,mesh(1)
            !
            ! r on crystal coordinates
            r(1) = real(ir1-1,dp)/mesh(1)
            r(2) = real(ir2-1,dp)/mesh(2)
            r(3) = real(ir3-1,dp)/mesh(3)
            !
            ! Translate the r vector
            Sr = r - (tau_cryst(:,ja) - tau_cryst(:,ia))
            !
            ! Find the (i,j,k) indices of the translated r
            Sir1 = nint(mod(Sr(1)+1.0_dp,1.0_dp)*mesh(1)) + 1
            Sir2 = nint(mod(Sr(2)+1.0_dp,1.0_dp)*mesh(2)) + 1
            Sir3 = nint(mod(Sr(3)+1.0_dp,1.0_dp)*mesh(3)) + 1
            !
            ! Find the global indices
            ir = ir1 + mesh(1)*(ir2-1) + mesh(1)*mesh(2)*(ir3-1)
            Sir = Sir1 + mesh(1)*(Sir2-1) + mesh(1)*mesh(2)*(Sir3-1)
            !
            dvvq(ir,:,1) = dvvq(ir,:,1) + dvv(Sir,:,1) * phase / n_vecs
            dvvq(ir,:,2) = dvvq(ir,:,2) + dvv(Sir,:,2) * phase / n_vecs
            dvvq(ir,:,3) = dvvq(ir,:,3) + dvv(Sir,:,3) * phase / n_vecs
            !
          enddo
        enddo
      enddo
      !
    enddo ! ja

  end subroutine dVR_to_dVq_sc


  subroutine dVq_sc_2_dVq_uc(dvvq, q_cryst, dvq)
    !
    ! TODO: Add description
    !
    use siesta2ph_io, only: nr1, nr2, nr3
    use siesta2ph_utils, only: tpi, cmplx_0, cmplx_i
    !
    implicit none
    !
    real(kind=dp), dimension(:,:,:), intent(in) :: dvvq ! dV_q: Potential induced by the periodic displacement of an atom (u_ia*e^{iqR}) in the supercell cell (r in sc)
    real(kind=dp), dimension(3), intent(in) :: q_cryst ! in unit cell crystal coordinates
    complex(kind=dp), dimension(:,:,:), intent(inout) :: dvq ! dV_q: Periodic part of the potential induced by the periodic displacement of an atom (u_ia*e^{iqR}) in the unit cell (r in uc)
    !
    real(kind=dp), dimension(3) :: r_cryst
    complex(kind=dp) :: phase
    integer :: ir1, ir2, ir3, ir_uc, ir
    integer, dimension(3) :: mesh_uc


    ! Check input/output array dimensions
    if (size(dvvq, dim=1) /= mesh(1)*mesh(2)*mesh(3)) stop "ERROR: dVq_sc_2_dVq_uc: wrong mesh dimensions in input"
    if (size(dvvq, dim=2) /= nspin) stop "ERROR: dVq_sc_2_dVq_uc: wrong spin dimensions in input"
    if (size(dvvq, dim=3) /= 3) stop "ERROR: dVq_sc_2_dVq_uc: wrong spatial dimensions in input"
    !
    mesh_uc = mesh / (/nr1, nr2, nr3/)
    if (size(dvq, dim=1) /= mesh_uc(1)*mesh_uc(2)*mesh_uc(3)) stop "ERROR: dVq_sc_2_dVq_uc: wrong mesh dimensions in output"
    if (size(dvq, dim=2) /= nspin) stop "ERROR: dVq_sc_2_dVq_uc: wrong spin dimensions in output"
    if (size(dvq, dim=3) /= 3) stop "ERROR: dVq_sc_2_dVq_uc: wrong spatial dimensions in output"
    !
    dvq = cmplx_0
    do ir3=1,mesh_uc(3)
      do ir2=1,mesh_uc(2)
        do ir1=1,mesh_uc(1)
          !
          ! r in unit cell crystal coordinates
          r_cryst(1) = real(ir1-1,dp)/mesh_uc(1)
          r_cryst(2) = real(ir2-1,dp)/mesh_uc(2)
          r_cryst(3) = real(ir3-1,dp)/mesh_uc(3)
          !
          ! e^{-iqr}
          phase = exp(-cmplx_i*tpi*dot_product(q_cryst, r_cryst))
          !
          ! Find global indices
          ir = ir1 + mesh(1)*(ir2-1) + mesh(1)*mesh(2)*(ir3-1)
          ir_uc = ir1 + mesh_uc(1)*(ir2-1) + mesh_uc(1)*mesh_uc(2)*(ir3-1)
          !
          !
          dvq(ir_uc,:,1) = dvvq(ir,:,1) * phase
          dvq(ir_uc,:,2) = dvvq(ir,:,2) * phase
          dvq(ir_uc,:,3) = dvvq(ir,:,3) * phase
          !
        enddo
      enddo
    enddo

  end subroutine dVq_sc_2_dVq_uc

end program siesta2dv
