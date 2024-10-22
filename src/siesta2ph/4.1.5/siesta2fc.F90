program siesta2fc
  !
  ! This program reads the forces for each irreducible displacements and
  ! calculates the full force constant matrix.0Finally, it applies all
  ! symmetries to the force constants and writtes them to a file.
  !
  use precision, only: dp
  use m_timestamp, only: timestamp
  use fdf, only: fdf_parallel
  use mpi_siesta, only: MPI_Comm_World
  use parallel, only: Node, Nodes
  !
  use siesta2ph_io, only: stdout, outdir, prefix, v0dir, kpath_file, verbose
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
  !
  real(kind=dp), allocatable, dimension(:,:,:,:) :: fc

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
    call die("siesta2fc: ERROR: This code is not prepared to run in multiple nodes")
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
  ! Read the output files for each irreducible displacements and compute FC
  !
  call compute_force_constants()
  !
  ! Save FC to a file
  call write_fc_phonopy("fc.dat_1")
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Apply permutation symmetry to the FC
  call apply_permutation_sym()
  !
  ! Save FC to a file
  call write_fc_phonopy("fc.dat_2")
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Apply space group symmetry to FC
  call apply_sym()
  !
  ! Save FC to a file
  call write_fc_phonopy("fc.dat_3")
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Compute dynamical matrix
  write(stdout,*) "- Computing dynamical matrices..."
  !
  call dynamical_matrix()
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Compute phonon dispersion
  if (trim(kpath_file) /= "") then
    !
    write(stdout,*) "- Computing phonon dispersion for "//trim(kpath_file)//"..."
    !
    call phonon_bands()
    !
  endif
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  call apply_sym()
  !
  call write_fc_phonopy("fc.dat_4")
  !
  !
  call apply_asr()
  !
  call write_fc_phonopy("fc.dat_5")
  !
  !
  call apply_asr()
  !
  call write_fc_phonopy("fc.dat_6")
  !
  !
  call apply_permutation_sym()
  call apply_sym()
  call apply_permutation_sym()
  !
  call write_fc_phonopy("fc.dat_7")
  !
  !
#ifdef MPI
  call MPI_Finalize( MPIerror )
#endif
  stop

contains

  subroutine apply_asr()
    !
    ! TODO: Add description
    !
    use siesta2ph_system, only: nat
    !
    implicit none
    !
    real(kind=dp) :: asr_tmp
    integer :: ia_, ja_, id_, jd_


    do id_=1,3
      do jd_=1,3
        do ia_=1,nat
          asr_tmp = 0.0_dp
          do ja_=1,nat
            asr_tmp = asr_tmp + fc(id_,jd_,ia_,ja_)
          end do
          fc(id_,jd_,ia_,:) = fc(id_,jd_,ia_,:) - asr_tmp/nat
        end do
      end do
    end do

  end subroutine apply_asr


  subroutine write_fc_phonopy(filename)
    !
    ! Write the force constant matrix in the format used by phonpy
    !
    use siesta2ph_io, only: outdir, phdir
    use siesta2ph_system, only: nat
    use siesta2ph_utils, only: Hartree2eV, Bohr2Ang
    !
    use siesta2ph_io, only: find_free_unit
    !
    implicit none
    !
    character(*), intent(in) :: filename
    !
    integer :: iounit, ios
    character(len=10) :: ia_c, ja_c
    integer :: ia_, ja_


    iounit = find_free_unit()
    open(unit=iounit, name=trim(outdir)//trim(phdir)//trim(filename), &
         action="write", status="replace", iostat=ios)
    if (ios/=0) stop "ERROR: write_fc_phonopy: Error opening force constants file."
    !
    write(iounit,"(i4,x,i4)") nat, nat
    do ia_=1,nat
      write(ia_c,"(i10)") ia_
      do ja_=1,nat
        write(ja_c,"(i10)") ja_
        write(iounit,"(a,x,a)") trim(adjustl(ia_c)),trim(adjustl(ja_c))
        write(iounit,"(2(3f22.15))") fc(1,:,ia_,ja_)*Hartree2eV/Bohr2Ang
        write(iounit,"(2(3f22.15))") fc(2,:,ia_,ja_)*Hartree2eV/Bohr2Ang
        write(iounit,"(2(3f22.15))") fc(3,:,ia_,ja_)*Hartree2eV/Bohr2Ang
      enddo
    enddo
    !
    close(iounit)

  end subroutine write_fc_phonopy


  subroutine write_fc_QE(filename)
    !
    ! Write the force constant matrix in the format used by Quantum Espresso
    ! QE writtes the force constants in Ry/Boher^2
    !
    use siesta2ph_io, only: outdir, phdir
    use siesta2ph_system, only: nat
    use siesta2ph_utils, only: Hartree2Ry
    !
    use siesta2ph_io, only: find_free_unit
    !
    implicit none
    !
    character(*), intent(in) :: filename
    !
    integer :: iounit, ios
    character(len=10) :: ia_c, ja_c
    integer :: ia_, ja_


    iounit = find_free_unit()
    open(unit=iounit, name=trim(outdir)//trim(phdir)//trim(filename), &
         action="write", status="replace", iostat=ios)
    if (ios/=0) stop "ERROR: write_fc_QE: Error opening force constants file."
    !
    write(iounit,"(i4,x,i4)") nat, nat
    do ia_=1,nat
      write(ia_c,"(i10)") ia_
      do ja_=1,nat
        write(ja_c,"(i10)") ja_
        write(iounit,"(a,x,a)") trim(adjustl(ia_c)),trim(adjustl(ja_c))
        write(iounit,"(2(3f22.15))") fc(1,:,ia_,ja_)*Hartree2Ry
        write(iounit,"(2(3f22.15))") fc(2,:,ia_,ja_)*Hartree2Ry
        write(iounit,"(2(3f22.15))") fc(3,:,ia_,ja_)*Hartree2Ry
      enddo
    enddo
    !
    close(iounit)

  end subroutine write_fc_QE


  subroutine compute_force_constants()
    !
    ! Compute the force constants for the irreducible atoms using symmetries and the pseudo-inverse
    !
    use siesta2ph_io, only: outdir, v0dir, phdir, lpm, dx
    use siesta2ph_utils, only: eV2Hartree, Ang2Bohr
    use siesta2ph_system, only: nat, slabel, tau_cryst
    use siesta2ph_symmetry, only: irred_atm, nsite_sym, nsite_rot, irred_disp, site_atm, site_s_cart, site_sinv, nsym, s_cart, rtau_index
    !
    use siesta2ph_io, only: find_free_unit
    use siesta2ph_linalg, only: pseudoinverse, ainv
    !
    implicit none
    !
    ! loop variables
    integer :: ia, ja, ka, Sia, Sja
    integer :: id
    integer :: isym, isite_sym, isym_disp_p, isym_disp_n, iirred
    integer :: idisp
    !
    real(kind=dp), dimension(3) :: R, SR
    real(kind=dp), dimension(3,3) :: FC_ij
    !
    integer :: nat_fa
    real(kind=dp), dimension(3,3) :: disp
    character(len=256) :: dispp_filename
    character(len=256) :: dispn_filename
    integer :: dispp_unit, dispn_unit
    logical :: found
    logical, dimension(nat,nat) :: done
    real(kind=dp), allocatable, dimension(:,:) :: fp, fn
    real(kind=dp), dimension(3) :: drift
    real(kind=dp), allocatable, dimension(:,:,:) :: Sdisp, Sf
    integer :: nsym_disp
    !
    integer :: ios


    write(stdout,*) "- Computing force constants..."
    !
    disp(:,1) = (/dx, 0.0_dp, 0.0_dp/)
    disp(:,2) = (/0.0_dp, dx, 0.0_dp/)
    disp(:,3) = (/0.0_dp, 0.0_dp, dx/)
    !
    allocate(fp(3,nat))
    allocate(fn(3,nat))
    fp = 0.0_dp
    fn = 0.0_dp
    ! fc(id,jd,ia,ja) is the force that the atom ja feels in the direction jd
    ! due to a displacement of the atom ia in the direction id
    allocate(fc(3,3,nat,nat))
    fc = 0.0_dp
    !
    done = .false.
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
      nsym_disp = 2*count(irred_disp(:,ia))*nsite_rot(isite_sym)
      allocate(Sf(3,nsym_disp,nat))
      allocate(Sdisp(3,nsym_disp,nat))
      Sf = 0.0_dp
      Sdisp = 0.0_dp
      !
      idisp = 0
      !
      do id=1,3
        !
        if ( .not. irred_disp(id,ia) ) cycle
        !
        idisp = idisp + 1
        iirred = iirred + 1
        write(stdout,"(a,i3,2(a5,i3))") "irreducible displacement", iirred, ": ia=", ia, ", id=", id
        !
        !
        if (lpm) then
          write(dispp_filename,"(a5,i4.4,a1)") "disp-", 2*iirred-1, "/"
          !
          dispp_unit = find_free_unit()
          open(unit=dispp_unit, file=trim(outdir)//trim(phdir)//trim(dispp_filename)//trim(slabel)//".FA", &
          status="old", action="read", iostat=ios)
          if ( ios /= 0 ) stop "Error opening file dispp SystemLabel.FA"
          !
          write(dispn_filename,"(a5,i4.4,a1)") "disp-", 2*iirred, "/"
          dispn_unit = find_free_unit()
          open(unit=dispn_unit, file=trim(outdir)//trim(phdir)//trim(dispn_filename)//trim(slabel)//".FA", &
          status="old", action="read", iostat=ios)
          if ( ios /= 0 ) stop "Error opening file dispn SystemLabel.FA"
        else
          write(dispp_filename,"(a5,i4.4,a1)") "disp-", iirred, "/"
          !
          dispp_unit = find_free_unit()
          open(unit=dispp_unit, file=trim(outdir)//trim(phdir)//trim(dispp_filename)//trim(slabel)//".FA", &
          status="old", action="read", iostat=ios)
          if ( ios /= 0 ) stop "Error opening file dispp SystemLabel.FA"
          !
          dispn_unit = find_free_unit()
          open(unit=dispn_unit, file=trim(outdir)//trim(v0dir)//trim(slabel)//".FA", status="old", action="read", iostat=ios)
          if ( ios /= 0 ) stop "Error opening file dispn SystemLabel.FA"
        endif
        !
        ! Read number of atoms and check if it is correct
        read(dispp_unit,*) nat_fa
        if (nat_fa .ne. nat) stop "nat wrong in dispp SystemLabel.FA"
        read(dispn_unit,*) nat_fa
        if (nat_fa .ne. nat) stop "nat wrong in dispn SystemLabel.FA"
        !
        ! Read the forces
        do ja=1,nat
          read(dispp_unit,*) nat_fa, fp(:,ja)
          read(dispn_unit,*) nat_fa, fn(:,ja)
        enddo
        fp = fp * eV2Hartree / Ang2Bohr ! from eV/Ang to Hartree/Bohr
        fn = fn * eV2Hartree / Ang2Bohr ! from eV/Ang to Hartree/Bohr
        !
        ! Substract drift force
        drift = sum(fp, dim=2)/nat
        do ja=1,nat
          fp(:,ja) = fp(:,ja) - drift
        enddo
        drift = sum(fn, dim=2)/nat
        do ja=1,nat
          fn(:,ja) = fn(:,ja) - drift
        enddo
        !
        !
        ! Apply site symmetries
        do isym=1,nsite_rot(isite_sym)
          !
          isym_disp_p = 2*nsite_rot(isite_sym)*(idisp-1)+isym
          isym_disp_n = 2*nsite_rot(isite_sym)*(idisp-1)+nsite_rot(isite_sym)+isym
          !
          do ja=1,nat
            !
            ! Find interatomic distance
            R = tau_cryst(:,ja) - tau_cryst(:,ia)
            !
            SR = matmul(site_sinv(:,:,isym,isite_sym), R)
            !
            ! Find the original atom with SR interatomic distance
            found=.false.
            do ka=1,nat
              R = ( tau_cryst(:,ia) + SR(:) ) - tau_cryst(:,ka) - nint(( tau_cryst(:,ia) + SR(:) ) - tau_cryst(:,ka))
              if ( sum(R**2) .lt. 0.0000001_dp ) then
                found=.true.
                exit
              endif
            enddo
            if (.not.found) stop "STOP: Atom index not found"
            !
            ! Positive displacement
            Sf(:,isym_disp_p,ja) = matmul(site_s_cart(:,:,isym,isite_sym), fp(:,ka))
            Sdisp(:,isym_disp_p,ja) = matmul(site_s_cart(:,:,isym,isite_sym), disp(:,id))
            !
            ! Negative displacement
            Sf(:,isym_disp_n,ja) = matmul(site_s_cart(:,:,isym,isite_sym), fn(:,ka))
            if (lpm) then
              Sdisp(:,isym_disp_n,ja) = matmul(site_s_cart(:,:,isym,isite_sym), -disp(:,id))
            else
              Sdisp(:,isym_disp_n,ja) = (/0.0_dp, 0.0_dp, 0.0_dp/)
            endif
            !
          enddo ! ja
          !
        enddo ! isym
        !
        close(dispp_unit)
        close(dispn_unit)
        !
      enddo ! id
      !
      !
      do ja=1,nat
        !
        call pseudoinverse(Sdisp(:,:,ja))
        !
        ! in Hartree/(Bohr^2)
        FC_ij = -matmul(Sdisp(:,:,ja), transpose(Sf(:,:,ja)))
        !
        ! Asign FC to all symmetry equivalent atoms
        do isym=1,nsym
          !
          Sia = rtau_index(ia,isym)
          Sja = rtau_index(ja,isym)
          !
          if (done(Sia,Sja)) cycle
          !
          fc(:,:,Sia,Sja) = matmul( s_cart(:,:,isym), matmul( FC_ij(:,:), ainv(s_cart(:,:,isym)) ) )
          !
          done(Sia,Sja) = .true.
          !
        enddo ! isym
        !
      enddo ! ja
      !
      deallocate(Sf)
      deallocate(Sdisp)
      !
    enddo ! ia

  end subroutine compute_force_constants


  subroutine dynamical_matrix()
    !
    ! TODO: Add description
    !
    use siesta2ph_io, only: outdir, phdir, nr1, nr2, nr3, irreducible_q
    use siesta2ph_system, only: nat, at, tau_cryst, amass, ityp
    use siesta2ph_utils, only: tpi, cmplx_0, cmplx_i, pmass, Hartree2meV, au2THz
    !
    use siesta2ph_io, only: find_free_unit
    use siesta2ph_utils, only: eqvect
    use siesta2ph_linalg, only: diagonalize_cmat
    use siesta2ph_symmetry, only: get_q_mesh
    !
    implicit none
    !
    complex(kind=dp), allocatable, dimension(:,:,:) :: dm
    real(kind=dp), allocatable, dimension(:) :: w
    complex(kind=dp), dimension(3,3) :: dm_local
    real(kind=dp), dimension(3,27) :: R_vecs_cryst
    real(kind=dp), dimension(27) :: vecs_lenght
    real(kind=dp), dimension(3) :: R_cryst, vec_cryst, vec_cart
    complex(kind=dp) :: phase
    real(kind=dp) :: sqrt_mm, min_lenght
    integer :: ia, ja, ka
    integer :: id, jd
    integer :: ir1, ir2, ir3
    integer :: iq, jq, iqirr
    integer :: ivec, n_vecs
    integer :: nat_uc
    character(len=256) :: filename
    integer :: iounit
    real(kind=dp), dimension(3, nr1*nr2*nr3) :: qmesh_cryst
    integer, dimension(nr1*nr2*nr3) :: qirr_map
    integer :: nq_star
    logical, dimension(nr1*nr2*nr3) :: qirr_done


    ! Get q-mesh
    qmesh_cryst = 0.0_dp
    qirr_map = 0
    !
    call get_q_mesh(nr1, nr2, nr3, qmesh_cryst, qirr_map)
    !
    nat_uc = nat/nr1/nr2/nr3 ! The number of atoms in the unit cell
    !
    allocate(w(3*nat_uc))
    allocate(dm(3*nat_uc, 3*nat_uc, nr1*nr2*nr3))
    dm = cmplx_0
    !
    if (.not.irreducible_q) then
      do iq=1,nr1*nr2*nr3
        qirr_map(iq) = iq
      enddo
    endif
    !
    qirr_done = .false.
    iqirr = 0
    do iq = 1, nr1*nr2*nr3
      !
      if (qirr_done(qirr_map(iq))) cycle
      !
      iqirr = iqirr + 1
      !
      iounit = find_free_unit()
      write(filename,"(a3,i4.4,a4)") "dyn", iqirr, ".dat"
      open(unit=iounit, file=trim(outdir)//trim(phdir)//trim(filename), action="write", status="replace")
      !
      ! Count symmetry equivalent q points
      nq_star = 0
      do jq=1, nr1*nr2*nr3
        !
        if (qirr_map(jq) /= iq) cycle
        nq_star = nq_star + 1
        !
      enddo
      !
      write(iounit,"(a10,i3)") "nq_star = ", nq_star
      !
      do jq=1, nr1*nr2*nr3
        !
        if (qirr_map(jq) /= iq) cycle
        !
        do ia=1,nat_uc
          do ja=1,nat_uc
            !
            sqrt_mm = sqrt(amass(ityp(ia)) * amass(ityp(ja))) * pmass ! in a.u.
            !
            dm_local = cmplx_0
            do ka=1,nat
              !
              if ( .not. eqvect(tau_cryst(:,ja)*(/nr1, nr2, nr3/), tau_cryst(:,ka)*(/nr1, nr2, nr3/)) ) cycle ! i.e. if the atoms are equivalent by translation
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
                    ! lattice vector from tau(ja) to tau(ka)+R
                    R_vecs_cryst(:, ivec) = ((tau_cryst(:,ka) + R_cryst) - tau_cryst(:,ja)) ! in supercell crystal coordinates
                    R_vecs_cryst(:, ivec) = R_vecs_cryst(:, ivec)*(/nr1, nr2, nr3/) ! transform to unit cell crystal coordinates
                    !
                    ! lenght of the vector from tau(ia) to tau(ka)+R
                    vec_cryst = ((tau_cryst(:,ka) + R_cryst) - tau_cryst(:,ia)) ! in supercell crystal coordinates
                    vec_cart = matmul(at, vec_cryst)
                    vecs_lenght(ivec) = sqrt(sum(vec_cart**2))
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
                  phase = phase + exp(tpi*cmplx_i*dot_product(qmesh_cryst(:, jq), R_vecs_cryst(:, ivec)))
                endif
              enddo
              !
              dm_local = dm_local + fc(:,:,ia,ka) * phase / sqrt_mm / n_vecs ! in a.u. (Hartree/Bohr^2*/me = Hartree^2/hbar^2)
              !
            enddo ! ka
            !
            dm(((ia-1)*3+1):((ia-1)*3+3), ((ja-1)*3+1):((ja-1)*3+3), jq) = dm(((ia-1)*3+1):((ia-1)*3+3), ((ja-1)*3+1):((ja-1)*3+3), jq) + dm_local
            !
          enddo ! ja
        enddo ! ia
        !
        dm(:,:, jq) = 0.5_dp*(dm(:,:, jq) + conjg(transpose(dm(:,:, jq))))
        !
        write(iounit,"(a10,3f10.6)") "q_cryst = ", qmesh_cryst(:, jq)
        do ia=1,nat_uc
          do ja=1,nat_uc
            sqrt_mm = sqrt(amass(ityp(ia)) * amass(ityp(ja))) * pmass ! in a.u.
            write(iounit,"(2i5)") ia, ja
            do id=1,3
              write(iounit,"(3(2f12.8),2x)") ( dm((ia-1)*3+id, (ja-1)*3+jd, jq)*sqrt_mm, jd=1,3 )
            enddo
          enddo
        enddo
        !
        qirr_done(jq) = .true.
        !
      enddo ! jq
      !
      close(iounit)
      !
      call diagonalize_cmat(dm(:,:,iq), w)
      print"(a12,x,a2,x,i4,x,100f18.6)", "freq. (meV):", "iq", iq, sign(sqrt(abs(w)), w)*Hartree2meV
      ! print"(a12,x,a2,x,i4,x,100f18.6)", "freq. (THz):", "iq", jq, sign(sqrt(abs(w)), w)*au2THz ! QE
      !
    enddo ! iq

  end subroutine dynamical_matrix


  subroutine phonon_bands()
    !
    ! TODO: Add description
    !
    use siesta2ph_io, only: outdir, phdir, kpath_file, nr1, nr2, nr3
    use siesta2ph_system, only: nat, at, tau_cryst, amass, ityp, bg, alat
    use siesta2ph_utils, only: tpi, cmplx_0, cmplx_i, pmass, Hartree2meV, au2cm1, au2THz
    !
    use siesta2ph_io, only: find_free_unit
    use siesta2ph_utils, only: eqvect
    use siesta2ph_linalg, only: diagonalize_cmat
    !
    implicit none
    !
    real(kind=dp), allocatable, dimension(:,:) :: qpath_cryst, qpath_cart
    integer :: nqpath
    complex(kind=dp), allocatable, dimension(:,:) :: dm
    real(kind=dp), allocatable, dimension(:) :: w
    complex(kind=dp), dimension(3,3) :: dm_local
    real(kind=dp), dimension(3,27) :: R_vecs_cryst
    real(kind=dp), dimension(27) :: vecs_lenght
    real(kind=dp), dimension(3) :: q_cryst, R_cryst, vec_cryst, vec_cart
    complex(kind=dp) :: phase
    real(kind=dp) :: sqrt_mm, min_lenght, dist
    integer :: ia, ja, ka
    integer :: ir1, ir2, ir3
    integer :: iq
    integer :: ivec, n_vecs
    integer :: nat_uc
    integer :: iounit, iostat
    character(len=256) :: filename


    ! Read kpath_file
    iounit = find_free_unit()
    filename = trim(outdir)//trim(kpath_file)
    open( unit=iounit, file=trim(filename), iostat=iostat, form="formatted", &
          status="old", action="read" )
    if (iostat .ne. 0) stop "ERROR: phonon_bands: Error opening kpath_file."
    !
    read(iounit, *) nqpath
    !
    allocate(qpath_cryst(3, nqpath))
    allocate(qpath_cart(3, nqpath))
    !
    do iq = 1, nqpath
      read(iounit, *) qpath_cryst(:, iq)
      qpath_cart(:, iq) = tpi/alat*matmul(bg, qpath_cryst(:, iq))*(/nr1, nr2, nr3/) ! We have to multiply by nr1, nr2, nr3 because bg is for the supercell
    enddo
    !
    close(iounit)
    !
    ! Compute phonon dispersion
    nat_uc = nat/nr1/nr2/nr3 ! The number of atoms in the unit cell
    !
    allocate(w(3*nat_uc))
    allocate(dm(3*nat_uc, 3*nat_uc))
    dm = cmplx_0
    !
    iounit = find_free_unit()
    filename = trim(outdir)//trim(phdir)//"phonons.dat"
    open(unit=iounit, file=trim(filename), action="write", status="replace", iostat=iostat)
    if (iostat .ne. 0) stop "ERROR: phonon_bands: Error opening phonons.dat."
    !
    do iq=1,nqpath
      !
      if (iq==1) then
        dist = 0.0_dp
      else
        dist = dist + sqrt(sum((qpath_cart(:,iq)-qpath_cart(:,iq-1))**2))
      endif
      !
      q_cryst = qpath_cryst(:,iq)
      dm = cmplx_0
      !
      do ia=1,nat_uc
        do ja=1,nat_uc
          !
          sqrt_mm = sqrt(amass(ityp(ia)) * amass(ityp(ja))) * pmass ! in a.u.
          !
          dm_local = cmplx_0
          do ka=1,nat
            !
            if ( .not. eqvect(tau_cryst(:,ja)*(/nr1, nr2, nr3/), tau_cryst(:,ka)*(/nr1, nr2, nr3/)) ) cycle ! i.e. if the atoms are equivalent by translation
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
                  ! lattice vector from tau(ja) to tau(ka)+R
                  R_vecs_cryst(:, ivec) = ((tau_cryst(:,ka) + R_cryst) - tau_cryst(:,ja)) ! in supercell crystal coordinates
                  R_vecs_cryst(:, ivec) = R_vecs_cryst(:, ivec)*(/nr1, nr2, nr3/) ! transform to unit cell crystal coordinates
                  !
                  ! lenght of the vector from tau(ia) to tau(ka)+R
                  vec_cryst = ((tau_cryst(:,ka) + R_cryst) - tau_cryst(:,ia)) ! in supercell crystal coordinates
                  vec_cart = matmul(at, vec_cryst)
                  vecs_lenght(ivec) = sqrt(sum(vec_cart**2))
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
                phase = phase + exp(tpi*cmplx_i*dot_product(q_cryst, R_vecs_cryst(:, ivec)))
              endif
            enddo
            !
            dm_local = dm_local + fc(:,:,ia,ka) * phase / sqrt_mm / n_vecs ! in a.u. (Hartree/Bohr^2*hbar*/me = Hartree**2)
            !
          enddo ! ka
          !
          dm(((ia-1)*3+1):((ia-1)*3+3), ((ja-1)*3+1):((ja-1)*3+3)) = dm(((ia-1)*3+1):((ia-1)*3+3), ((ja-1)*3+1):((ja-1)*3+3)) + dm_local
          !
        enddo ! ia
      enddo ! ja
      !
      dm = 0.5_dp*(dm + conjg(transpose(dm)))
      !
      call diagonalize_cmat(dm, w)
      write(iounit,"(f12.8,x,100f12.6)") dist, sign(sqrt(abs(w)), w)*Hartree2meV
      ! write(iounit,"(f12.8,x,100f12.6)") dist, sign(sqrt(abs(w)), w)*au2cm1 ! matdyn
      ! write(iounit,"(f12.8,x,100f12.6)") dist/tpi, sign(sqrt(abs(w)), w)*au2THz ! phonopy
      !
    enddo ! iq
    close(iounit)

  end subroutine phonon_bands


  subroutine apply_permutation_sym()
    !
    ! TODO: Add description
    !
    use siesta2ph_system, only: nat
    !
    implicit none
    !
    integer :: ia, ja, i, j


    do ia=1,nat
      do ja=1,nat
        do i=1,3
          do j=1,3
            fc(i,j,ia,ja) = 0.5_dp * ( fc(i,j,ia,ja) + fc(j,i,ja,ia) )
            fc(j,i,ja,ia) = fc(i,j,ia,ja)
          enddo
        enddo
      enddo
    enddo
    !
  end subroutine apply_permutation_sym


  subroutine apply_sym()
    !
    ! Symetrizy the FC's
    !
    use siesta2ph_system, only: nat
    use siesta2ph_symmetry, only: nsym, rtau_index, s_cart
    !
    use siesta2ph_linalg, only: ainv
    !
    implicit none
    !
    integer :: isym, ia, ja
    integer :: Sia, Sja
    real(kind=dp), dimension(3,3) ::FC_ij
    logical, dimension(nat,nat) :: done


    done = .false.
    do ia=1,nat
      do ja=1,nat
        !
        if (done(ia,ja)) cycle
        !
        FC_ij = 0.0_dp
        do isym=1,nsym
          !
          Sia = rtau_index(ia,isym)
          Sja = rtau_index(ja,isym)
          !
          FC_ij(:,:) = FC_ij(:,:) + matmul( ainv(s_cart(:,:,isym)), matmul( fc(:,:,Sia,Sja), s_cart(:,:,isym) ) )/nsym
          !
        enddo ! isym
        !
        do isym=1,nsym
          !
          Sia = rtau_index(ia,isym)
          Sja = rtau_index(ja,isym)
          !
          fc(:,:,Sia,Sja) = matmul( s_cart(:,:,isym), matmul( FC_ij(:,:), ainv(s_cart(:,:,isym)) ) )
          !
          done(Sia,Sja) = .true.
          !
        enddo ! isym
        !
      enddo ! ja
    enddo ! ia
    !
  end subroutine apply_sym

end program siesta2fc
