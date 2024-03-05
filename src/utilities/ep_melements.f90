program ep_melements

  use kinds, only: dp
  use intw_input_parameters, only: nk1, nk2, nk3, &
                                   nq1, nq2, nq3, nqirr, mesh_dir, &
                                   ep_mat_file, TR_symmetry, &
                                   read_input
  use intw_reading, only: nkpoints_QE, kpoints_QE, nspin, noncolin, nsym, &
                          s, nbands, nG_max, nat, nspin, tau, bg, &
                          nr1, nr2, nr3, num_bands_intw, &
                          read_parameters_data_file_xml, &
                          get_gvec, &
                          read_kpoints_data_file_xml, &
                          set_num_bands
  use intw_pseudo, only: vkb, vkqb, &
                         read_all_pseudo, &
                         init_KB_projectors
  use intw_utility, only: get_timing, &
                          find_free_unit, &
                          switch_indices, &
                          generate_kmesh, &
                          conmesurate_and_coarser, &
                          ainv, &
                          find_k_1BZ_and_G
  use intw_useful_constants, only: cmplx_0, cmplx_1
  use intw_symmetries, only: full_mesh, IBZ, QE_folder_nosym, QE_folder_sym, nosym_G, sym_G, symlink, &
                             symtable, rtau_index, rtau, rtau_cryst, &
                             rot_atoms, &
                             find_size_of_irreducible_k_set, &
                             find_the_irreducible_k_set, &
                             allocate_symmetry_related_k, &
                             find_inverse_symmetry_matrices_indices, &
                             allocate_and_build_spin_symmetry_matrices, &
                             set_symmetry_relations, multable
  use intw_fft, only: nl, &
                      generate_nl, &
                      allocate_fft
  use intw_ph, only: nqmesh, qmesh, dvq_local, dvpsi, QE_folder_nosym_q, QE_folder_sym_q, &
                     nosym_G_q, sym_G_q, symlink_q, q_irr, q_irr_cryst, &
                     read_ph_information_xml, &
                     read_allq_dvr, &
                     get_dv, &
                     calculate_local_part_dv
  use intw_allwfcs, only: allocate_and_get_all_irreducible_wfc, &
                          get_psi_general_k_all_wfc

  !================================================================================
  ! Declare the variables
  !================================================================================
  implicit none

  !k point related variables
  logical                  :: k_points_consistent
  integer                  :: ik, nk_irr, nkmesh
  real(dp),allocatable     :: kmesh(:,:)
  real(dp)                 :: kpoint(3), kpoint_cart(3), kqpoint_cart(3)

  !q point related variables
  logical                  :: q_points_consistent
  real(dp)                 :: qpoint(3)
  integer                  :: iq, nq_irr
  logical                  :: full_mesh_q, IBZ_q

  !phonon related variables
  integer                  :: imode

  !time related variables
  real(dp)                 :: time1, time2

  !ep interaction related variables
  integer                  :: ep_unit
  complex(dp), allocatable :: ep_mat_el(:,:,:,:,:,:)

  !wave function realted variables
  integer, allocatable     :: list_igk(:)
  integer, allocatable     :: list_igkq(:)
  complex(dp), allocatable :: wfc_k (:,:,:)
  complex(dp), allocatable :: wfc_kq (:,:,:)

  !ep fitxa
  integer                  :: record_lengh, ierr

  !local/aux variables
  integer                  :: npw, npwq

  integer                  :: ig, ibnd, jbnd, ispin, jspin
  logical                  :: read_status
  character(len=4)         :: iq_loc

  complex(dp), external :: zdotc


  20 format(A)
  30 format(A,F8.2,6X,A)
  !
  !
  !================================================================================
  ! Begining
  !================================================================================
  !
  call get_timing(time1)
  !
  write(*,20) '====================================================='
  write(*,20) '|                  program me                       |'
  write(*,20) '|        ---------------------------------          |'
  write(*,20) '====================================================='
  !
  !
  !================================================================================
  ! Read the necessary information from standard input file
  !================================================================================
  !
  write(*,20) '|       - Waiting for input file...                 |'
  !
  call read_input(read_status)
  !
  if (read_status) stop
  !
  !
  !================================================================================
  ! Check if the k-mesh and q-mesh are conmesurate
  !================================================================================
  !
  if (.not.conmesurate_and_coarser(nk1,nk2,nk3,nq1,nq2,nq3)) then
    write(*,20) '**********************************************************'
    write(*,20) '*ERROR                                                  '
    write(*,20) '*   the electron k and phonon q are not                 '
    write(*,20) '*   conmesurate and the k grid does not contain         '
    write(*,20) '*   the phonon q grid                                   '
    write(*,20) '**********************************************************'
    stop
  endif
  !
  !
  !================================================================================
  ! Read the parameters from the SCF QE calculation
  !================================================================================
  !
  write(*,20) '|       - Reading calculation parameters...         |'
  !
  call read_parameters_data_file_xml()
  !
  !
  !================================================================================
  ! Set the number of wave functions
  !================================================================================
  !
  call set_num_bands()
  !
  !
  !================================================================================
  ! Print spin information
  !================================================================================
  !
  if (nspin==1) then
    write(*,20) '|       - Paramagnetic calculation nspin=1          |'
  elseif (nspin==2) then
    write(*,20) '|       - Spin-polarized calculation nspin = 2      |'
    if (noncolin) write(*,20) '|       - Non-collinear spin calculation            |'
  else
    write(*,20) '*****************************************************'
    write(*,20) '* ERROR: Allowed values for nspin are 1 or 2        *'
    write(*,20) '*            program stops.                         *'
    write(*,20) '*****************************************************'
    stop
  endif
  !
  write(*,20) '|           ---------------------------------       |'
  !
  !
  !================================================================================
  ! Set symmetry arrays
  !================================================================================
  !
  write(*,20) '|       - Setting symmetry arrays...                |'
  !
  ! Set the rotation table for each atom and symmetry
  allocate(rtau_index(nat,nsym))
  allocate(rtau(3,nsym,nat))
  allocate(rtau_cryst(3,nsym,nat))
  !
  call rot_atoms(nat,nsym,tau)
  !
  ! Compute the indices of the inverse rotation matrices
  call find_inverse_symmetry_matrices_indices()
  !
  ! Calculate the multiplication talble for symmetry operations
  call multable(nsym,s,symtable)
  !
  ! Set up spin_symmetry_matrices, needed to rotate wave functions and indueced potential for non-colinear calculations
  call allocate_and_build_spin_symmetry_matrices(nsym)
  !
  write(*,20) '|           ---------------------------------       |'
  !
  !
  !================================================================================
  ! Set up the gvec array and all FFT variables
  !================================================================================
  !
  write(*,20) '|       - Reading G vectors...                      |'
  !
  call get_gvec()
  !
  ! Allocate useful variables
  call allocate_fft()
  !
  ! Generate some important indices for FFT
  call generate_nl()
  !
  write(*,20) '|           ---------------------------------       |'
  !
  !
  !================================================================================
  ! Read PPs
  !================================================================================
  !
  write(*,20) '|       - Reading pseudopotentials...               |'
  !
  call read_all_pseudo ()
  !
  write(*,20) '|                    PPs are OK                     |'
  !
  ! Allocate and set PP variables
  call allocate_nlpot
  !
  call allocate_phq
  !
  call init_pp
  !
  write(*,20) '|           ---------------------------------       |'
  !
  !
  !================================================================================
  ! Read the kpoints from the calculation
  !================================================================================
  !
  write(*,20) '|       - Reading k-points...                       |'
  !
  allocate(kpoints_QE(3,nkpoints_QE))
  !
  call read_kpoints_data_file_xml(kpoints_QE)
  !
  !
  !================================================================================
  ! Build the wave function's k-mesh
  !================================================================================
  !
  write(*,20) '|       - Building k-mesh...                        |'
  !
  nkmesh = nk1*nk2*nk3
  allocate(kmesh(3,nkmesh))
  call generate_kmesh(kmesh,nk1,nk2,nk3)
  !
  ! Find the size of the irreducible set of k-points (IBZ)
  call find_size_of_irreducible_k_set(nk1,nk2,nk3,kmesh,nk_irr)
  !
  !
  !================================================================================
  ! Set symmetry relations between irreducible k-points and full k-mesh
  !================================================================================
  !
  ! Allocate arrays
  call allocate_symmetry_related_k(nk1,nk2,nk3,nsym)
  !
  ! Fill the symmetry arrays
  call set_symmetry_relations(nk1,nk2,nk3,nkpoints_QE,kpoints_QE,kmesh, k_points_consistent, &
                              QE_folder_nosym, QE_folder_sym, &
                              nosym_G, sym_G, symlink, full_mesh, IBZ)
  !
  !
  !================================================================================
  ! Check that the number of kpoints corresponds to either a full mesh or the IBZ
  !================================================================================
  !
  if (full_mesh .and. IBZ) then
    write(*,20) '|       - the kpoints present in the QE folders     |'
    write(*,20) '|         are consistent with a full 1BZ and a      |'
    write(*,20) '|         IBZ has also been found.                  |'
    write(*,20) '|           ---------------------------------       |'
  else if(IBZ) then
    write(*,20) '|       - the kpoints present in the QE folders     |'
    write(*,20) '|         are consistent with an IBZ.               |'
    write(*,20) '|           ---------------------------------       |'
  else
    write(*,20) '**********************************************************'
    write(*,20) '* The kpoints present in the QE folders are not consistent'
    write(*,20) '* with the parameters of the input file!                 '
    write(*,20) '**********************************************************'
    write(*,20) '* debug information:                                *'
    write(*,*) '*        nkpoints_QE = ',nkpoints_QE
    write(*,*) '*        nkmesh      = ',nkmesh
    write(*,*) '*        nk_irr      = ',nk_irr
    stop
  end if
  !
  !
  !================================================================================
  ! Read phonon information
  !================================================================================
  !
  write(*,20) '|       - Reading q-points...                      |'
  !
  ! Read q-points and irreducible patterns
  call read_ph_information_xml()
  !
  allocate(q_irr_cryst(3,nqirr))
  do iq=1,nqirr
    q_irr_cryst(:,iq) = matmul(ainv(bg),q_irr(:,iq))
  enddo
  !
  !
  !================================================================================
  ! Build the phonon q-mesh
  !================================================================================
  !
  write(*,20) '|       - Building q-mesh...                        |'
  !
  nqmesh = nq1*nq2*nq3
  allocate(qmesh(3,nqmesh))
  !
  call generate_kmesh(qmesh,nq1,nq2,nq3)
  !
  ! Find the size of the irreducible set of q-points (IBZ)
  call find_size_of_irreducible_k_set(nq1,nq2,nq3,qmesh,nq_irr)
  !
  !
  !================================================================================
  ! Symmetry relations between irreducible q-points and full q-mesh
  !================================================================================
  !
  allocate(QE_folder_nosym_q(nqmesh))
  allocate(QE_folder_sym_q(nqmesh))
  allocate(nosym_G_q(3,nqmesh))
  allocate(sym_G_q(3,nqmesh))
  allocate(symlink_q(nqmesh,2))
  !
  call set_symmetry_relations(nq1,nq2,nq3, nqirr, q_irr_cryst, qmesh, q_points_consistent, &
                              QE_folder_nosym_q, QE_folder_sym_q, &
                              nosym_G_q, sym_G_q, symlink_q, full_mesh_q, IBZ_q)
  !
  !
  !================================================================================
  ! Check that the number of kpoints corresponds to either a full mesh or the IBZ
  !================================================================================
  !
  if (full_mesh_q .and. IBZ_q) then
    write(*,20) '|       - the qpoints present in the QE folders     |'
    write(*,20) '|         are consistent with a full 1BZ and a      |'
    write(*,20) '|         IBZ has also been found.                  |'
    write(*,20) '|           ---------------------------------       |'
  else if(IBZ_q) then
    write(*,20) '|       - the qpoints present in the QE folders     |'
    write(*,20) '|         are consistent with an IBZ.               |'
    write(*,20) '|           ---------------------------------       |'
  else
    write(*,20) '**********************************************************'
    write(*,20) '* The qpoints present in the QE folders are not consistent'
    write(*,20) '* with the parameters of the input file!                 '
    write(*,20) '**********************************************************'
    write(*,20) '* debug information:                                *'
    write(*,*) '*        nqpoints_QE = ',nqirr
    write(*,*) '*        nqmesh      = ',nqmesh
    write(*,*) '*        nq_irr      = ',nq_irr
    stop
  end if
  !
  !
  !================================================================================
  ! Read all wave functions
  !================================================================================
  !
  write(*,20) '|       - Reading wave functions...                 |'
  !
  call allocate_and_get_all_irreducible_wfc()
  !
  write(*,20) '|           ---------------------------------       |'
  !
  !
  !================================================================================
  ! Read the induced potentials
  !================================================================================
  !
  write(*,20) '|       - Reading induced potentials...             |'
  !
  call read_allq_dvr(nqirr,3*nat)
  !
  write(*,20) '====================================================='
  !
  !
  !================================================================================
  ! Compute matrix elements
  !================================================================================
  !
  write(*,20) '|       - Computing matrix elements...              |'
  !
  ! Allocate wfc related variables
  allocate(list_igk(nG_max))
  allocate(list_igkq(nG_max))
  allocate(wfc_k(nG_max,num_bands_intw,nspin))
  allocate(wfc_kq(nG_max,num_bands_intw,nspin))
  !
  ! Allocate induced potential related variables (hauek behekoak ph_module.mod-n definituta daude eta aldagai globalak dira kontuz)
  allocate(dvq_local(nr1*nr2*nr3,3*nat,nspin,nspin))
  allocate(dvpsi(nG_max,num_bands_intw,nspin,nspin,3*nat))
  !
  ! Allocate matrix elements variable
  allocate(ep_mat_el(nk1*nk2*nk3,num_bands_intw,num_bands_intw,nspin,nspin,3*nat))
  !
  do iq=1,nqmesh
    !
    ep_mat_el = cmplx_0
    !
    qpoint = qmesh(:,iq)
    write(*,"(a,i4,100f12.6)") "qpoint", iq, qpoint
    !
    if (                iq <   10) write(iq_loc,"(i1)")iq
    if ( 10 <= iq .and. iq <  100) write(iq_loc,"(i2)")iq
    if (100 <= iq .and. iq < 1000) write(iq_loc,"(i3)")iq
    !
    ep_unit = find_free_unit()
    inquire(iolength=record_lengh) ep_mat_el
    open(unit=ep_unit, iostat=ierr, &
          file=trim(mesh_dir)//trim(ep_mat_file)//trim('_')//adjustl(iq_loc), &
          form='unformatted', status='unknown', access='direct', recl=record_lengh)
    if (ierr /= 0 ) stop 'Error opening ep_mat_file'
    !
    call set_nqxq(qpoint)
    !
    ! Potentzialaren alde induzitua kalkulatu simetria erabiliz (errotazioz beharrezkoa izanez).
    dvq_local = cmplx_0
    call get_dv(iq,qpoint,3*nat,nspin,dvq_local)
    !
    ! Alde induzituari (goian), KB pseudopotentzialaren(pp) deribatuaren ALDE LOKALA gehitu.
    call phq_init(matmul(bg,qpoint))
    call calculate_local_part_dv(qpoint,nat,nspin,dvq_local)
    !
    ! Bi subroutina hauek (goikoak), biak batera joan behar dira beti).
    do ik=1,nkmesh
      !
      kpoint=kmesh(:,ik)
      kpoint_cart=matmul(bg,kpoint)
      kqpoint_cart=matmul(bg,kpoint+qpoint)
      !
      write(*,'(a,i4,a,3(f15.8))') "ik= ", ik, ' k= ', kpoint
      !
      ! Uhina lortu RAM memorian dauden uhin irreduzibleak errotatuta
      call get_psi_general_k_all_wfc(kpoint       , list_iGk , wfc_k )
      call get_psi_general_k_all_wfc(kpoint+qpoint, list_iGkq, wfc_kq)
      !
      ! Ordenan jartzen ditugu G bektoreak, k=0-n nola dauden ordenatuta arabera
      npw=0
      do iG=1,nG_max
        if (list_iGk(iG)==0) exit
        npw=npw+1
      enddo
      !
      npwq=0
      do iG=1,nG_max
        if (list_iGkq(iG)==0) exit
        npwq=npwq+1
      enddo
      !
      ! Hemen KB potentzial ez lokaleko |beta> funtzioak kalkulatzen dira (k eta k+q puntuetarako hurrenez hurren)
      vkb =cmplx_0
      vkqb=cmplx_0
      !
      call init_KB_projectors(npw ,nG_max,list_iGk ,kpoint       ,vkb )
      call init_KB_projectors(npwq,nG_max,list_iGkq,kpoint+qpoint,vkqb)
      !
      ! psi_k uhinak, potentzial induzitua + KB pp-ren ALDE LOKALAREN
      ! batuketarekin biderkatu (output-a:dvpsi): dv_local x |psi_k> (G)
      call dvqpsi_local(3*nat, nG_max, num_bands_intw, nspin, list_iGk, list_iGkq, wfc_k, dvq_local, &
                                              dvpsi(1:nG_max,1:num_bands_intw,1:nspin,1:nspin,1:3*nat))
      !
      ! psi_k uhinak KB potentzialaren alde ez lokalarekin biderkatu eta emaitza dvpsi aldagaiari gehitu:
      !                    dvpsi^q_k --> dvpsi^q_k + D^q_mode [ KB ] |psi_k> (G)
      !                                  (lokala) + (ez lokala)
      call multiply_psi_by_dvKB( kpoint, qpoint, nspin, nbands, nG_max, list_iGk, list_iGkq, wfc_k, dvpsi)
      !
      do imode=1,3*nat ! Osagai kanonikoak, ez dira moduak, kontuz
        !
        ! Matrize elementuak kalkulatu
        do jspin=1,nspin
          do ispin=1,nspin
            do jbnd=1,num_bands_intw
              do ibnd=1,num_bands_intw
                !
                ep_mat_el(ik,ibnd,jbnd,ispin,jspin,imode) = zdotc( nG_max, wfc_kq(:,ibnd,ispin), 1, dvpsi(:,jbnd,ispin,jspin,imode), 1 )
                !
              enddo !ibnd
            enddo !jbnd
          enddo !is
        enddo !js
        !
      enddo !imode
      !
    enddo !ik
    !
    write(unit=ep_unit, rec = 1, iostat = ierr) ep_mat_el(:,:,:,:,:,:)
    !
    close(unit=ep_unit)
    !

#ifdef DEBUG
      ! Write the matrix elements to a formatted file to compare with QE matrix elements
      write(123400+iq,"(3f10.6)") qpoint
      write(123400+iq,"(i4)") nkmesh
      do ik=1,nkmesh
        kpoint = kmesh(:,ik)
        write(123400+iq,"(i4,3f10.6)") ik, kpoint
        write(123400+iq,"(i4,3f10.6)") 2*(ik-1)+2, kpoint+qpoint
      enddo
      write(123400+iq,"(4i4,2f20.15)")((((ik, ibnd, jbnd, imode, sum(ep_mat_el(ik,ibnd,jbnd,:,:,imode)), ibnd=1,num_bands_intw), jbnd=1,num_bands_intw), ik=1,nkmesh), imode=1,3*nat)
#endif

    !
  enddo !iq


  deallocate (ep_mat_el)
  !
  !
  !================================================================================
  ! Finish
  !================================================================================
  !
  call get_timing(time2)
  !
  write(*,20) '|                     ALL DONE                       |'
  write(*,30) '|     total time: ',time2-time1,' seconds            |'
  write(*,20) '====================================================='



contains

  subroutine dvqpsi_local(nmode,nG_max,nbands,nspin,list_iGk,list_iGkq,wfc_k,dvq_local,dvpsi_local)

    implicit none

    !I/O variables

    integer,intent(in) :: nmode,nG_max,nbands,nspin,list_iGk(nG_max),list_iGkq(nG_max)
    complex(kind=dp),intent(in) :: dvq_local(nr1*nr2*nr3,nmode,nspin,nspin),wfc_k(nG_max,nbands,nspin)
    complex(kind=dp),intent(inout) :: dvpsi_local(nG_max,nbands,nspin,nspin,nmode)

    !local variables

    integer :: ibnd,ispin,ig,imode,ir
    complex(kind=dp) :: wfc_r(nr1*nr2*nr3,nspin,nspin),wfc_r1(nr1*nr2*nr3,nspin)

    dvpsi_local(:,:,:,:,:)=cmplx_0
    !
    do imode =1,nmode !3*nat
      do ibnd=1,nbands
        !
        wfc_r1(:,:)=cmplx_0
        wfc_r(:,:,:)=cmplx_0
        !
        do ispin=1,nspin !elektroi uhin funtzioak espazio errelarela pasatu. SO kasuan ispin osagaiz osagai.
          do ig=1,nG_max
            !
            if (list_iGk(ig)==0) exit
            wfc_r1(nl(list_iGk(ig)),ispin)=wfc_k(ig,ibnd,ispin)
            !
          enddo !ig
          !
          call cfftnd(3,(/nr1,nr2,nr3/),1,wfc_r1(:,ispin))
          !
        enddo !ispin

        !
        if (nspin==2) then
          !
          do ir=1,nr1*nr2*nr3
            !
            do ispin=1,nspin
              do jspin=1,nspin
                !
                wfc_r(ir,ispin,jspin)=dvq_local(ir,imode,ispin,jspin)*wfc_r1(ir,jspin)
                !
              enddo !jspin
            enddo !ispin
            !
          enddo !ir
          !
        else !nspin
          !
          do ir=1,nr1*nr2*nr3
            !
            wfc_r(ir,1,1)=dvq_local(ir,imode,1,1)*wfc_r1(ir,1)
            !
          enddo !ir
          !
        endif !nspin
        !
        do ispin=1,nspin
          do jspin=1,nspin
            !
            call cfftnd(3,(/nr1,nr2,nr3/),-1,wfc_r(:,ispin,jspin))
            !
            do ig=1,nG_max
              !
              if (list_iGkq(ig)==0) exit
              !
              ! dvpsi_local(ig,ibnd,ispin,jspin,imode)=wfc_r(nl(ig),ispin,jspin)
              dvpsi_local(ig,ibnd,ispin,jspin,imode)=wfc_r(nl(list_iGkq(ig)),ispin,jspin)
              !
            enddo !ig
          enddo !jspin
        enddo !ispin
        !
      enddo !ibnd
      !
    enddo !imode

  end subroutine dvqpsi_local

end program ep_melements
