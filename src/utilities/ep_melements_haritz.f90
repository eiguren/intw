program ep_melements

  use kinds, only: dp

  use intw_input_parameters, only: mesh_dir, prefix, ph_dir, nk1, nk2, nk3, &
                                   nq1, nq2, nq3, nqirr, ep_mat_file, calc_epmat

  use intw_reading, only: at, bg, alat, nat, ntyp, nr1, nr2, nr3, nG_max, nspin, &
                          npol, nbands, lspinorb, noncolin, gvec, ngm, &
                          read_parameters_data_file_xml, get_gvec

  use intw_pseudo, only: deallocate_upfeak, read_all_pseudo, average_pp

  use intw_useful_constants, only: cmplx_i, cmplx_0, cmplx_1, tpi

  use intw_symmetries, only: rtau, rtau_index, inverse_indices, nsym, s, ftau, &
                             rotate_wfc_test

  use intw_fft, only: allocate_fft, nl

  use intw_ph, only: nqmesh, qmesh, dvq_local, dvpsi, q_irr_cryst, rot_k_index, &
                     symlink_q, QE_folder_sym_q, get_dv, calculate_local_part_dv

  use intw_w90, only: u_mesh

  use w90_parameters, only: num_wann, num_bands

  use intw_allwfcs, only: get_psi_general_k_all_wfc

  use intw_uspp, only: vkb, vkqb, nhtol

  use intw_matrix_elements, only: get_elec_phon_matrix_element_convolution, &
                                  get_plane_wave_matrix_element_convolution
  use intw_utility, only: find_free_unit, generate_kmesh, conmesurate_and_coarser, &
                          switch_indices, find_k_1bz_and_g, get_timing


  !================================================================================
  !       Declare the variables
  !================================================================================
  implicit none

  !k point related variables
  integer                  :: ikpt, iqpt, ikqpt, ikpt_1, ikpt_2, ikpt_k, ikpt_kq
  real(dp)                 :: kpoint(3)
  real(dp), allocatable    :: kmesh(:,:)
  real(dp)                 :: kpoint_in_1bz(3)

  !q point related variables
  real(dp)                 :: qpoint(3), qpoint_1bz(3)

  !phonon related variables
  integer                  :: mode, q_index_irr, q_index, imode

  !symmetry variables
  integer                  :: isym, jsym, s_index, imq, i_sym, j_sym
  integer, allocatable     :: nsym_sgk(:), sindex_sgk(:,:)
  complex(dp), allocatable :: unit_sym_sgk( :,:,:,:), umat(:,:)

  !time related variables
  real(dp)                 :: time, time1, time2

  !q point related variables
  integer                  :: iq, iqq, nq1_, nq2_, nq3_

  !ep interaction related variables
  integer                  :: ep_unit
  complex(dp), allocatable :: aep_mat_el(:,:,:,:,:,:,:), ep_mat_el(:,:,:,:,:,:)
  complex(dp), allocatable :: aW_rot_ep_mat_el(:,:,:,:,:,:,:), W_rot_ep_mat_el(:,:,:,:,:,:)
  complex(dp), allocatable :: fr(:,:,:), fg(:,:,:)

  !wave function realted variables information
  integer, allocatable     :: list_igk (:)
  integer, allocatable     :: list_igkq(:)
  integer, allocatable     :: list_igk_aux (:)
  integer, allocatable     :: list_igk_orig (:)

  complex(dp), allocatable :: wfc_k (:,:,:) ! nG_max is defined in reading
  complex(dp), allocatable :: wfc_kq (:,:,:)
  complex(dp), allocatable :: wfc_k_aux (:,:,:)
  complex(dp), allocatable :: wfc_k_orig (:,:,:)

  real(dp), allocatable    :: QE_eig_k(:)
  real(dp), allocatable    :: QE_eig_kq(:)

  !fft related
  integer                  :: ngk, ngkq, nr(3)
  integer                  :: G(1:3), GKQ_bz(3), G_plusk(3), G_pluskq(3)

  !ep fitxa
  integer                  :: record_lengh, ierr

  !local/aux variables
  integer                  :: na
  complex(dp)              :: suma_c, term
  integer                  :: LEN
  integer                  :: nbands_loc
  integer                  :: npw, npwq

  integer                  :: i, j, k
  integer                  :: ii, jj, kk
  integer                  :: jr, kr
  integer                  :: ig, ibnd, jbnd, ipol, jpol
  integer                  :: switch

  complex(dp),allocatable :: wfc_k_r(:)
  integer :: gg1(3), gg2(3), nG
  logical :: found
  !
  !
  complex(dp), allocatable    :: dv_local    (:,:,:,:)



  !
  !================================================================================
  !       Talk to the user
  !================================================================================
  !
  write(*,20) '====================================================='
  write(*,20) '|                  program me                       |'
  write(*,20) '|        ---------------------------------          |'
  write(*,20) '====================================================='
  write(*,20) '|    waiting for input file...                      |'
  !
  !================================================================================
  !       read the input file
  !       Read in the necessary information from standard input
  !================================================================================
  !
  !================================================================================
  ! READ INPUT
  !================================================================================
  ! This is the nspin=1 calculation
  ! Then we will be interested with the nspin=2 calculation
  mesh_dir = '/scratch/hgarai/Kalkuluak/Fe-MgO-Ag100/1ML/PH/smearing_0.01/V0/'
  prefix = 'Fe-MgO-Ag100'
  mesh_dir = '/home/haritz/Kalkuluak/Ag/Fe-MgO-Ag100/1ML/ph/Fe-MgO-Ag100/smearing_0.01/V0/'
  prefix = 'Fe-MgO-Ag100'
  nk1 = 1
  nk2 = 1
  nk3 = 1
  ph_dir= '/scratch/hgarai/Kalkuluak/Fe-MgO-Ag100/1ML/PH/smearing_0.01/V0/'
  nq1=1
  nq2=1
  nq3=1
  nqirr=1
  !
  !================================================================================
  !       read the parameters from the SCF QE calculation, in
  !       particular, read in the symmetry matrices!!!
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
  !
  !haritz: ngm irakurri da read_parameters_data_file_xml()-en, eztao get_ngm() erabili beharrik
  allocate(gvec(3,ngm)) ! it is deallocated in the bottom of this main code.
  !
  call get_gvec() ! read them from gvectors.dat file
  !
  !allocate useful variables
  !
!  call allocate_fft()
!  !
!  !generate some important indices for FFT
!  !
!  call generate_nl()
  !
  !
  !
  if (nspin==1) then
     !
     write(*,20) '|       - The calculation is paramagnetic nspin=1   |'
     write(*,20) '|                                                   |'
     write(*,20) '|           ---------------------------------       |'
     !
  elseif (nspin==2) then
     !
     write(*,20) '|       - Spin calculation nspin = 2                |'
     !
     if (noncolin) then
        !
        write(*,20) '|         Non-collinear Spin calculation            |'
        !
     endif
     !
     write(*,20) '|           ---------------------------------       |'
     !
  else
     !
     write(*,20) '*****************************************************'
     write(*,20) '* ERROR: Allowed values for nspin are 1 or 2        *'
     write(*,20) '*            program stops.                         *'
     write(*,20) '*****************************************************'
     !
     stop
     !
  endif
  !
  !
  !================================================================================
  !       Tell the user what is in the QE folders
  !================================================================================
  !
  !
  !
  write(*,20) '|       - reading pseudopotentials from UPF files   |'
  write(*,20) '|            (defined in .save data files)          |'
  write(*,20) '|                                                   |'
  !
  call read_all_pseudo ()
  if (.not.lspinorb) call average_pp(ntyp)
  !
  write(*,20) '|                    PPs are OK                     |'
  write(*,20) '|           ---------------------------------       |'
  !
  if (.not.conmesurate_and_coarser(nk1,nk2,nk3,nq1,nq2,nq3)) then
     !
     write(*,20)      '**********************************************************'
     write(*,20)       '*ERROR                                                  '
     write(*,20)       '*   the electron k and phonon q are not                 '
     write(*,20)       '*   conmesurate and the k grid does not contain         '
     write(*,20)       '*   the phonon q grid                                   '
     write(*,20)      '**********************************************************'
     !
     stop
     !
  endif
  stop
  !
  !
  !
  !
  !================================================================================
  !       Allocate wfc related files
  !================================================================================
  !
  ! haritz
  num_bands = nbands
  ! haritz
  nbands_loc=num_bands
  !
  allocate (list_igk(nG_max))
  allocate (list_igkq(nG_max))
  allocate (list_igk_aux(nG_max))
  allocate (list_igk_orig(nG_max))
  !
  allocate (wfc_k(nG_max,num_bands,nspin))
  allocate (wfc_kq(nG_max,num_bands,nspin))
  allocate (wfc_k_aux(nG_max,num_bands,nspin))
  allocate (wfc_k_orig(nG_max,num_bands,nspin))
  !
  allocate (wfc_k_r(nr1*nr2*nr3))
  !
  allocate (QE_eig_k(num_bands))
  allocate (QE_eig_kq(num_bands))
  !
  allocate (fr(nr1*nr2*nr3,num_bands,nspin))
  allocate (fg(nG_max,num_bands,nspin))
  !
  !================================================================================
  !hauek behekoak ph_module.mod-n definituta daude eta aldagai globalak dira kontuz.
  !
  allocate(dvq_local(nr1*nr2*nr3,3*nat,npol,npol))
  allocate(dvpsi(nG_max,nbands_loc,npol,npol,3*nat))
  !
  allocate( kmesh(3,1), qmesh(3,1) )
  call generate_kmesh(kmesh,nk1,nk2,nk3)
  call generate_kmesh(qmesh,nq1,nq2,nq3)
  !
  allocate(ep_mat_el(nk1*nk2*nk3,num_bands,num_bands,nspin,nspin,3*nat))
  allocate(W_rot_ep_mat_el(nk1*nk2*nk3,num_wann,num_wann,nspin,nspin,3*nat)) ! matrize elementuak Wannier espazioan (matrize unitarioarekin biderkatuta)
  allocate(aep_mat_el(nqmesh,nk1*nk2*nk3,num_bands,num_bands,nspin,nspin,3*nat))
  allocate(aW_rot_ep_mat_el(nqmesh,nk1*nk2*nk3,num_wann,num_wann,nspin,nspin,3*nat))
  !
  aep_mat_el(:,:,:,:,:,:,:)=cmplx_0
  aW_rot_ep_mat_el(:,:,:,:,:,:,:)=cmplx_0
  !
  !
  call allocate_nlpot1
  !
  call allocate_phq
  !
  call init_us_1
  !
  ep_mat_el = cmplx_0
  !
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
  !
  if (calc_epmat) then
    !
    ep_unit=find_free_unit()
    inquire(iolength=record_lengh) ep_mat_el
    open( unit=ep_unit, file=trim(trim(mesh_dir)//trim(ep_mat_file)), iostat=ierr, &
          form='unformatted', status='unknown', access='direct', recl=record_lengh )
    !
    !
    !
    !-q puntuaren menpeko aldagai bat allocatu (hau kenduko dugu laster).
    !
    call allocate_nlpot2 (qpoint)
    !
    !-potentzialaren alde induzitua kalkulatu simetria erabiliz (errotazioz beharrezkoa izanez).
    !
    dvq_local=cmplx_0
    call get_dv(iq,qpoint,3*nat,dvq_local)
    !
    !-alde induzituari (goian), KB pseudopotentzialaren(pp) deribatuaren ALDE LOKALA gehitu.
    !
    call phq_init(matmul(bg,qpoint))
    call calculate_local_part_dv(qpoint,nat,npol,dvq_local)
    !
    !-bi subroutina hauek (goikoak), biak batera joan behar dira beti).
    !
    !
    kpoint=kmesh(:,1)
    !
    !
    !-uhina lortu RAM memorian dauden uhin irreduzibleak errotatuta. (k+q)
    !
    call get_psi_general_k_all_wfc(.true.,kpoint       ,list_iGk ,wfc_k ,QE_eig_k ,G_plusk)
    call get_psi_general_k_all_wfc(.true.,kpoint+qpoint,list_iGkq,wfc_kq,QE_eig_kq,G_pluskq)
    !
    !-k eta k+q puntuen indizeak topatu. ikpt_1 ik-aldagaiaren berdina izan
    ! behar da, baina hurrengo eran segurtasuna irabazten dugu.
    ! k-ren indizea: ikpt_k
    !
    call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i,j,k,kpoint_in_1bz,GKQ_bz)
    call switch_indices(nk1,nk2,nk3,ikpt_k,i,j,k,+1)
    !
    !
    ! Ordenan jartzen ditugu G bektoreak, k=0-n nola dauden ordenatuta arabera
    !
    nG=0
    do iG=1,nG_max
      if (list_iGk(iG)==0) exit
      nG=nG+1
    enddo
    npw=nG
    nG=0
    do iG=1,nG_max
      if (list_iGkq(iG)==0) exit
      nG=nG+1
    enddo
    npwq=nG
    !
    !-hemen KB potentzial ez lokaleko |beta> funtzioak kalkulatzen dira (k eta k+q puntuetarako hurrenez hurren).
    !
    vkb=cmplx_0
    vkqb=cmplx_0
    call init_us_2(npw ,nG_max,list_iGk ,kpoint       ,vkb)
    call init_us_2(npwq,nG_max,list_iGkq,kpoint+qpoint,vkqb)
    !
    !-psi_k uhinak, potentzial induzitua + KB pp-ren ALDE LOKALAREN
    ! batuketarekin biderkatu (output-a:dvpsi): dv_local x |psi_k> (G)
    !
    call dvqpsi_local (3*nat, nG_max, nbands_loc, npol, list_iGk, list_iGkq, wfc_k, dvq_local, &
                                            dvpsi(1:nG_max,1:nbands_loc,1:npol,1:npol,1:3*nat))
    !
    !-psi_k uhinak KB potentzialaren alde ez lokalarekin biderkatu eta emaitza dvpsi aldagaiari gehitu:
    !                    dvpsi^q_k --> dvpsi^q_k + D^q_mode [ KB ] |psi_k> (G)
    !                                  (lokala) + (ez lokala)
    call dvqpsi_us_only (matmul(bg, kpoint), matmul(bg, qpoint), nat, nG_max, nbands_loc, &
                                                  npol, list_iGk, list_iGkq, wfc_k, dvpsi)
    !
    !-QE-ren subroutina goikoaren berdina egiteko.
    !
    do imode=1,3*nat ! Osagai kanonikoak, ez dira moduak, kontuz.
      !matrize elementuak kalkulatu
      call get_elec_phon_matrix_element_convolution((/0,0,0/),list_iGkq, &
           list_iGkq, wfc_kq, dvpsi(1:nG_max,1:nbands_loc,1:npol,1:npol,imode), &
           aep_mat_el(iq,ikpt_k,1:nbands_loc,1:nbands_loc,1:npol,1:npol,imode))
    enddo !imode
    !
    !
    write(unit=ep_unit,  rec = 1, iostat = ierr)  aep_mat_el      (iq,:,:,:,:,:,:)
    !
    close(unit=ep_unit)
    !
  else
    !
    ep_unit=find_free_unit()
    inquire(iolength=record_lengh) ep_mat_el
    open( unit=ep_unit, file=trim(trim(mesh_dir)//trim('./')//trim(ep_mat_file)), &
          iostat=ierr, form='unformatted', status='unknown', access='direct', recl=record_lengh )
    !
    write(*,*)'Reading Reciprocal Bloch matrix elements in the coarse mesh for iq = ',iq,'/',nqmesh,ierr
    read(unit=ep_unit,rec=1,iostat=ierr) aep_mat_el(iq,:,:,:,:,:,:)
    !
    close(unit=ep_unit)
    !
    write(222,*) aep_mat_el(iq,:,:,:,:,:,:)
    !
  end if ! calc_epmat

  deallocate (ep_mat_el)
  deallocate (W_rot_ep_mat_el)
  deallocate (aep_mat_el)
  deallocate (aW_rot_ep_mat_el)

  deallocate(nsym_sgk)
  deallocate(sindex_sgk)

  stop


  call deallocate_upfeak ()
  !================================================================================
  !       Finish
  !================================================================================
  call get_timing(time2)
  write(*,20) '|                     ALL DONE                       |'
  write(*,30) '|     total time: ',time2-time1,' seconds            |'
  write(*,20) '====================================================='


20 format(A)
30 format(A,F8.2,6X,A)

contains

  subroutine calculate_rotation_unitary_matrices()
    !
    complex(dp)  :: unit_mat(nbands_loc, nbands_loc)
    integer :: ik
    !
    unit_mat = cmplx_0
    !
    do ibnd=1,nbands_loc
      unit_mat(ibnd, ibnd) = cmplx_1
    enddo
    !
    ik=1
    !
    kpoint(:) = kmesh(:,ik)
    !
    call get_psi_general_k_all_wfc(.true., kpoint, list_iGk_orig , wfc_k_orig ,  QE_eig_k,  G_plusk)
    !
    do isym=1,nsym
      !
      !
      call find_k_1BZ_and_G(kpoint,nk1, nk2, nk3, i ,j, k, kpoint_in_1bz, GKQ_bz)
      call switch_indices (nk1, nk2, nk3, ikpt_k, i, j, k, +1)
      !
      call get_psi_general_k_all_wfc(.true., kpoint_in_1bz, list_iGk , wfc_k ,  QE_eig_k,  G_plusk)
      !
      list_iGk_aux = list_iGk_orig
      wfc_k_aux    = wfc_k_orig
      !
      G_pluskq=nint(kmesh(:,ikpt_k)-kpoint)
      !
      call rotate_wfc_test (wfc_k_aux,list_iGk_aux, wfc_kq, list_iGkq,         &
                 isym, s(:,:,isym),  ftau(:,isym) , (/0,0,0/))
      !
      call get_plane_wave_matrix_element_convolution  ((G_pluskq-G_plusk)*0, list_iGk, list_iGkq, wfc_k, wfc_kq , &
                 unit_sym_sgk( ik, isym, 1:nbands_loc,1:nbands_loc)  )
      !
      umat=  unit_sym_sgk( ikpt_k, isym, 1:nbands_loc,1:nbands_loc)
      !
      write(12334,*)ik,isym
      do i=1,nbands_loc
        write(12334,"(100(x(f8.3,a,f8.3,a)))") (real(umat(i,j)), " + ",aimag(umat(i,j)),"*I,", j=1,nbands_loc)
      enddo
      !
      !Test. Determinant of Det[u^{-1}.u] must be = 1.
      umat=matmul(umat,conjg(transpose(umat)))
      !
      write(10000,*)ik,isym
      do i=1,nbands_loc
        write(10000,"(100(x(f8.3,a,f8.3,a)))") (real(umat(i,j)), " + ",aimag(umat(i,j)),"*I,", j=1,nbands_loc)
      enddo
      !
      !call diagonalize_cmat (nbands_loc, umat, QE_eig_k )
      !suma_c=cmplx_1
      !do i=1,nbands_loc
      !   suma_c=suma_c*QE_eig_k(i)
      !enddo
      !if (abs(suma_c-cmplx_1)>1E-3) then
      !write(*,*)"ERROREA: The rotation is not described by a unitary matrix."
      !stop
      !end if
      !
      if (sum(abs(unit_mat - umat))>0.001) then
        write(*,*)"ERROREA: The rotation is not described by a unitary matrix."
        write(*,*) ik, isym
        write(*,*) sum(abs(unit_mat - umat))
        write(10000,*) '*****ERROREA*****'
        write(4325,*) ik, isym
  !      stop
      end if
      !
    enddo ! isym

  end subroutine calculate_rotation_unitary_matrices

  subroutine calculate_rotation_unitary_matrices_wannier()
    !
    complex(dp)  :: unit_mat(nbands_loc, nbands_loc)
    complex(dp)  :: wfc_k_W(nG_max,nbands_loc,nspin)
    complex(dp)  :: wfc_kq_W(nG_max,nbands_loc,nspin)
    complex(dp)  :: wfc_k_orig_W(nG_max,nbands_loc,nspin)
    complex(dp)  :: U_k(nbands_loc,num_wann)
    integer :: ik
    !
    unit_mat = cmplx_0
    !
    do ibnd=1,num_wann
      unit_mat(ibnd, ibnd) = cmplx_1
    enddo
    !
    ik=1
    !
    kpoint(:) = kmesh(:,ik)
    !
    call get_psi_general_k_all_wfc(.true., kpoint, list_iGk_orig , wfc_k_orig ,  QE_eig_k,  G_plusk)
    !
    call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3, i ,j, k, kpoint_in_1bz, GKQ_bz)
    call switch_indices (nk1, nk2, nk3, ikpt_k, i, j, k, +1)
    !
    U_k(:,:)=u_mesh(:,:,ikpt_k)
    call wfc_bands_to_wann(wfc_k_orig,U_k,wfc_k_orig_W)
    !
    do isym=1,nsym
      !
      !
      call find_k_1BZ_and_G(kpoint,nk1, nk2, nk3, i ,j, k, kpoint_in_1bz, GKQ_bz)
      call switch_indices (nk1, nk2, nk3, ikpt_k, i, j, k, +1)
      !
      call get_psi_general_k_all_wfc(.true., kpoint_in_1bz, list_iGk , wfc_k ,  QE_eig_k,  G_plusk)
      !
      U_k(:,:)=u_mesh(:,:,ikpt_k)
      call wfc_bands_to_wann(wfc_k,U_k,wfc_k_W)
      !
      list_iGk_aux = list_iGk_orig
      wfc_k_aux    = wfc_k_orig_W
      !
      G_pluskq=nint(kmesh(:,ikpt_k)-kpoint)
      !
      call rotate_wfc_test (wfc_k_aux,list_iGk_aux, wfc_kq_W, list_iGkq,         &
                 isym, s(:,:,isym),  ftau(:,isym) , (/0,0,0/))
      !
      do i=num_wann+1,nbands_loc
        wfc_kq_W(:,i,:)=(0.d0,0.d0)
        wfc_k_W(:,i,:)=(0.d0,0.d0)
      enddo
      !
      call get_plane_wave_matrix_element_convolution  ((G_pluskq-G_plusk)*0, list_iGk, list_iGkq, wfc_k_W, wfc_kq_W , &
                 unit_sym_sgk( ik, isym, 1:num_wann,1:num_wann)  )
      !
      ! umat=  unit_sym_sgk( ikpt_k, isym, 1:nbands_loc,1:nbands_loc)
      umat=  unit_sym_sgk( ik, isym, 1:num_wann,1:num_wann)
      !
      write(12334,*)ik,isym
      do i=1,num_wann
        write(12334,"(100(x(f8.3,a,f8.3,a)))") (real(umat(i,j)), " + ",aimag(umat(i,j)),"*I,", j=1,num_wann)
      enddo
      !
      !Test. Determinant of Det[u^{-1}.u] must be = 1.
      umat=matmul(umat,conjg(transpose(umat)))
      !
      write(10000,*)ik,isym
      do i=1,num_wann
        write(10000,"(100(x(f8.3,a,f8.3,a)))") (real(umat(i,j)), " + ",aimag(umat(i,j)),"*I,", j=1,num_wann)
      enddo
      !
      !call diagonalize_cmat (nbands_loc, umat, QE_eig_k )
      !suma_c=cmplx_1
      !do i=1,nbands_loc
      !   suma_c=suma_c*QE_eig_k(i)
      !enddo
      !if (abs(suma_c-cmplx_1)>1E-3) then
      !write(*,*)"ERROREA: The rotation is not described by a unitary matrix."
      !stop
      !end if
      !
      if (sum(abs(unit_mat - umat))>0.001) then
        write(*,*)"ERROREA: The rotation is not described by a unitary matrix."
        write(*,*) ik, isym
        write(*,*) sum(abs(unit_mat - umat))
        write(10000,*) '*****ERROREA*****'
        write(4326,*) ik, isym
        !stop
      end if
      !
    enddo ! isym

  end subroutine calculate_rotation_unitary_matrices_wannier

  subroutine get_ep_mat_el( nk1, nk2, nk3, kmesh, qpoint, num_wann, npol, ep_mat_el) !

    implicit none
    !input
    integer, intent(in)          :: num_wann, npol, nk1, nk2, nk3

    real   (kind=dp),intent(in ) :: qpoint(3), kmesh(3, nk1*nk2*nk3)
    !output

    complex(kind=dp), intent(out) :: ep_mat_el( nk1*nk2*nk3,num_wann, num_wann, npol, npol, 3*nat)
    !loca
    integer :: nkmesh, i,j,k, q_index, q_index_irr, s_index, imq, na
    real   (kind=dp) :: qpoint_1bz(3), qpoint_irr_cart(3), qpoint_rot(3), kqpoint_1bz(3)
    integer :: GKQ_bz(3), GKQ(3), rna, rotik, ik
    complex (dp) :: phase(nat)
    integer ::  s_inv_index, kpol, jpol, lpol, ipol, ibnd, jbnd, ikq
    real(dp) ::  s_cart(3,3),  s_crys(3,3)

    nkmesh=nk1*nk2*nk3

    ep_mat_el=cmplx_0

    call find_k_1BZ_and_G(qpoint,nq1,nq2,nq3,i ,j, k, qpoint_1bz, GKQ_bz)

    call switch_indices (nq1, nq2, nq3, q_index, i, j, k, +1)

    q_index_irr = QE_folder_sym_q(q_index)

    s_index = symlink_q(q_index,1) ; imq = symlink_q(q_index,2)

    qpoint_rot = matmul ( s(:,:,inverse_indices(symlink_q(q_index,1))), q_irr_cryst(:,q_index_irr))

    if (imq==1) qpoint_rot = - qpoint_rot

    qpoint_irr_cart=matmul(bg,q_irr_cryst(:,q_index_irr))

    GKQ = nint(qpoint - qpoint_rot)

    if (sum(abs(qpoint - qpoint_rot - dble(GKQ)))> 1E-4) then
       write(*,*)"ERROREA : qpoint not recovered by symmetry"; stop
    end if

    s_inv_index=inverse_indices(s_index)  ! Kontuz, hau behar bada alderantziz da.

    do na=1,nat
       phase(na) = &
            exp(  - cmplx_i * ( qpoint_irr_cart (1) * rtau (1, s_index, na)            &
            + qpoint_irr_cart (2) * rtau (2, s_index, na)            &
            + qpoint_irr_cart (3) * rtau (3, s_index, na) ) * tpi)
    enddo

    do ipol=1,3
       do jpol=1,3
          s_crys(ipol,jpol) = real(s (ipol,jpol,s_index),dp)
       enddo
    enddo

    do ipol = 1, 3
       do jpol = 1, 3
          s_cart (ipol, jpol) = 0.d0
          do kpol = 1, 3
             do lpol = 1, 3
                s_cart (ipol, jpol) = s_cart (ipol, jpol) + at (ipol, kpol) * &
                     s_crys (lpol, kpol)  * bg (jpol, lpol)
             enddo
          enddo
       enddo
    enddo


    !s_cart= 0.d0
    !do ipol = 1, 3
    !            s_cart (ipol, ipol) = 1.d0
    !enddo

    do ik=1,nkmesh

       rotik = rot_k_index( s_inv_index, ik, nk1, nk2, nk3, kmesh)

       call find_k_1BZ_and_G(kmesh(:,ik)+qpoint,nk1,nk2,nk3,i ,j, k, kqpoint_1bz, GKQ_bz)
       call switch_indices (nk1, nk2, nk3, ikq, i, j, k, +1)


       do na = 1, nat

          rna= rtau_index(na,s_index)

          !k_ep_mat_el( num_wann, num_wann, npol, npol, 3*nat)
          do ibnd=1,num_wann
             do jbnd=1, num_wann
                do ipol=1,3
                   do jpol=1,3



                      ep_mat_el (rotik, ibnd, jbnd, 1:npol, 1:npol, (rna-1)*3 + ipol) =  &
                           ep_mat_el (rotik, ibnd, jbnd, 1:npol, 1:npol, (rna-1)*3 + ipol)    &
                           + s_cart(ipol,jpol) * aep_mat_el(q_index_irr,  ik, i, j, 1:npol, 1:npol, (na-1)*3 + jpol) * phase(rna)
                   enddo !jpol
                enddo !ipol
             enddo
          enddo
       enddo !na
    enddo !ik

    if  (imq==1) ep_mat_el =conjg(ep_mat_el)


    return
  end subroutine get_ep_mat_el
!********************************************************************************************************
!-----------------------------------------------------------------------------------------------------
  subroutine dvqpsi_local(nmode,nG_max,nbands,npol,list_iGk,list_iGkq,wfc_k,dvq_local,dvpsi_local)
!-----------------------------------------------------------------------------------------------------

    implicit none

    !I/O variables

    integer,intent(in) :: nmode,nG_max,nbands,npol,list_iGk(nG_max),list_iGkq(nG_max)
    complex(kind=dp),intent(in) :: dvq_local(nr1*nr2*nr3,nmode,npol,npol),wfc_k(nG_max,nbands,npol)
    complex(kind=dp),intent(inout) :: dvpsi_local(nG_max,nbands,npol,npol,nmode)

    !local variables

    integer :: ibnd,ipol,ig,imode,ir
    complex(kind=dp) :: wfc_r(nr1*nr2*nr3,npol,npol),wfc_r1(nr1*nr2*nr3,npol)

    dvpsi_local(:,:,:,:,:)=cmplx_0
    !
    do imode =1,nmode !3*nat
       do ibnd=1,nbands
          !
          wfc_r1(:,:)=cmplx_0
          wfc_r(:,:,:)=cmplx_0
          !
          do ipol=1,npol !elektroi uhin funtzioak espazio errelarela pasatu. SO kasuan ipol osagaiz osagai.
             do ig=1,nG_max
                !
                if (list_iGk(ig)==0) exit
                wfc_r1(nl(list_iGk(ig)),ipol)=wfc_k(ig,ibnd,ipol)
                !
             enddo !ig
             !
             call cfftnd(3,(/nr1,nr2,nr3/),1,wfc_r1(:,ipol))
             !
          enddo !ipol
          !
!          if ((npol==2).and.spinorb_mag) then
          if (npol==2) then
             !
             do ir=1,nr1*nr2*nr3
                !
                do ipol=1,npol
                   do jpol=1,npol
                   !
                   wfc_r(ir,ipol,jpol)=dvq_local(ir,imode,ipol,jpol)*wfc_r1(ir,jpol)
                   !
                   enddo !jpol
                enddo !ipol
                !
             enddo !ir
             !
          else !npol
             !
             do ir=1,nr1*nr2*nr3
                !
                wfc_r(ir,1,1)=dvq_local(ir,imode,1,1)*wfc_r1(ir,1)
                !
             enddo !ir
             !
          endif !npol
          !
          do ipol=1,npol
             do jpol=1,npol
                !
                call cfftnd(3,(/nr1,nr2,nr3/),-1,wfc_r(:,ipol,jpol))
                !
                do ig=1,nG_max
                   !
                   if (list_iGkq(ig)==0) exit
                   !
!                   dvpsi_local(ig,ibnd,ipol,jpol,imode)=wfc_r(nl(ig),ipol,jpol)
                   dvpsi_local(ig,ibnd,ipol,jpol,imode)=wfc_r(nl(list_iGkq(ig)),ipol,jpol)
                   !
                enddo !ig
             enddo !jpol
          enddo !ipol
          !
       enddo !ibnd
    enddo !imode
    !
    return

  end subroutine dvqpsi_local
!********************************************************************************
  !----------------------------------------
  subroutine diagonalize_cmat (n,a,w)
    !----------------------------------------
    integer, intent(in)  :: n
    complex(dp),intent(inout) :: a(n,n)
    real(dp),intent(out) :: w(n)

    complex(dp) :: a_pack(n*(n+1)/2)

    integer :: i,j, nfound

    complex(dp) :: cwork(2*n)
    real   (dp) :: rwork(7*n)
    integer     :: iwork(5*n), ifail(n), info


    do j=1,n
       do i=1,j
          a_pack(i+((j-1)*j)/2)=a(i,j)
       enddo
    enddo

    ! Diagonalizing routine from lapack. Go see QE/.../lapack_atlas.f
    call ZHPEVX('V', 'A', 'U', n, a_pack, 0.0_dp, 0.0_dp, 0, 0, -1.0_dp,  &
         nfound, w, a, n,  cwork,  rwork,       iwork,          ifail, info)

  end subroutine diagonalize_cmat

  !--------------------------------------
  subroutine wfc_bands_to_wann(wfc_k,U_k,wfc_k_W)
   !----------------------------------------

   integer :: is, nG, i, j
   complex(dp),intent(in) :: wfc_k(nG_max,nbands_loc,nspin)
   complex(dp),intent(in) :: U_k(nbands_loc,num_wann)
   complex(dp),intent(out) :: wfc_k_W(nG_max,nbands_loc,nspin)

   wfc_k_W=(0.d0,0.d0)

   do nG=1,nG_max
   do is=1,nspin
   do i=1,num_wann
      do j=1,nbands_loc
         wfc_k_W(nG,i,is)=wfc_k_W(nG,i,is)+U_k(j,i)*wfc_k(nG,j,is)
      enddo
   enddo
   enddo
   enddo

  end subroutine wfc_bands_to_wann

end program ep_melements
