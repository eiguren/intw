program wanndenphon

  !
  !--------------------------------------------------------------------------------!
  !       mblanco 3/12/2020
  ! Read wannier functions and construct density using a Fermi-Dirac function.
  ! Then, following the previous "vacancy test" on real space, move atoms
  ! according to a q-phonon and displacement (we use the irreducible displacements
  ! found by QE) to create the distorted density.
  !
  !--------------------------------------------------------------------------------!
  !

  !================================================================================
  !       Load modules
  !================================================================================

  use intw_setup
  use intw_utility
  use intw_reading
  use intw_allwfcs
  use intw_setup
  use intw_ph
  use intw_matrix_elements
  use intw_fft
  use kinds, only: dp

  use intw_w90
  !use intw_input_parameters, only: chemical_potential

  use w90_parameters, only: num_wann,num_bands
  use w90_hamiltonian, only: irvec, nrpts, ndegen

  !================================================================================
  !       Declare the variables
  !================================================================================
  implicit none


  ! Local and aux variables
  integer                  :: io_ierr, wfc_unit_G, wfc_unit_r, i,j,k
  real(kind=dp)            :: time
  logical                  :: readstat
  logical, parameter       :: minus_k=.false.

  !k point related variables
  integer                  :: ik, iq, ikp, nkmesh, is, isp, nkpt, nkpt_sym, ik1,ik2,ik3
  integer, allocatable     :: ikqcorr(:), ikcorr(:), iqcorr(:)
  character(len=4)         :: ik_loc
  real(dp)                 :: kpoint(3), kpoint_1bz(3), qpoint_1bz(3)
  real(dp), allocatable    :: kpoints(:,:), kpoints_sym(:,:), kpoints_1bz(:,:), &
                              kqpoints_1bz(:,:), qpoints_1bz(:,:)


  !FOURIER TRANSFORMS related variables
  integer                  :: ig, jg, ii, jj, irr, ir, jr, kr
  integer                  :: Gq(3), Gk(3), Gkq(3), G(3)
  integer, allocatable     :: Gplus_ks(:,:), Gplus_kqs(:,:), Gplus_qs(:,:)
  integer, parameter       :: nrr1=50, nrr2=50, nrr3=50
  real(dp)                 :: xcrys(3), kvec(3), rcart(3)!, gvec_cart(3)
  complex(dp)              :: expiGr

  !wave function related variables information
  integer                  :: ibnd, jbnd
  integer, allocatable     :: list_igk(:)
  complex(dp), allocatable :: wfc_k_G(:,:,:)
  real(dp), allocatable    :: QE_eig_k(:)
  integer, allocatable     :: list_iGks(:,:)
  complex(dp), allocatable :: wfc_k_Gs(:,:,:,:)
  real(dp), allocatable    :: QE_eig_ks(:,:)

  real(kind=dp), parameter :: eV2Ry = 2.0_dp/27.211383860484776_dp

  !k-point symmetries
  integer, allocatable     :: equiv (:)

  ! wannier densities
  ! and phonon
  character(len=200)       :: filename
  integer                  :: record_length, ph_unit, io_err
  real(dp)                 :: ferwei, ktsmear, eiener, efermi, argum, u_amp
  character(len=1)         :: stencil
  integer                  :: nqpt, nmode, imode
  real(dp)                 :: qpoint_cart(3), qpoint(3)
  real(dp), allocatable    :: qpoints(:,:)
  real(dp), allocatable    :: uirr(:,:)
  complex(dp), allocatable :: u_mesh_1bz(:,:,:)
  complex(dp), allocatable :: wfcconv(:,:,:,:), drho(:,:,:,:),  drho_rs(:,:), &
                              drho_rs_iGr(:,:,:)

  !================================================================================
  ! Setup calculation and print xml info
  ! (includes wannier initialization)
  !================================================================================
  !

  call setup(time)

  !
  !================================================================================
  ! Read phonon data (here I only need the irreducible displacement patterns)
  !================================================================================
  !
  ! TODO Todavia no funciona en esta version, creo...
  ! como solo necesito el q y uirr (irreducible isplacement patterns), los tomo a mano
  ! de momento
!!  call  read_ph_information_xml()

  ! A single q point and three modes

  nqpt = 1
  allocate (qpoints(3,nqpt))
  !qpoint = (/ 0.5_dp, 0.5_dp, 0.5_dp /)   ! cart / alat
  qpoint = (/ 0.42857143, 0.42857143, 0.42857143 /)
  call cryst_to_cart (1, qpoint, at, -1)  ! transform to crystal coordinates DUDA: at,bg???
  qpoints(:,1) = qpoint


  nmode = 3
  allocate (uirr(3,nmode))  !cart / alat
  uirr(:,1) = (/ -0.5773502691896256_dp, -0.5773502691896256_dp, -0.5773502691896256_dp /)
  uirr(:,2) = (/ -0.4082482904638630_dp, -0.4082482904638630_dp, 0.8164965809277259_dp /)
  uirr(:,3) = (/ -0.7071067811865476_dp, 0.7071067811865476_dp, 0.0_dp /)
  uirr=uirr*0.1   ! instead of amplitude=1, use a smaller one, 0.1

  !
  !================================================================================
  !       Allocate wfc related files
  !================================================================================
  !
  allocate (list_igk(nG_max))
  allocate (wfc_k_G(nG_max,nbands,nspin))
  allocate (QE_eig_k(nbands))


  !
  !================================================================================
  ! Build and allocate the kmesh corresponding to the whole 1BZ
  ! instead of generate_kmesh I generate them by hand with the same switch as
  ! used in real space in next section.
  !
  ! For each kpoint, I use the 1bz contribution. 
  ! I read and store wfc(k) without add_G. I will add them "manually" later on
  ! when doing the Fourier transform
  !
  ! For each k I obtain k-q and look for the ik index of the vector in the 1bz.
  ! Thus I can just read once the wfc.
  ! I also store separately q point correspondence.
  !================================================================================
  !

  nkmesh = nk1*nk2*nk3
  call generate_kmesh(kmesh,nk1,nk2,nk3)
  nkpt = nkmesh

  allocate (kpoints(3,nkpt), kpoints_1bz(3,nkpt), kqpoints_1bz(3,nkpt))
  allocate (qpoints_1bz(3,nqpt))
  allocate (ikcorr(nkpt), ikqcorr(nkpt), iqcorr(nqpt))
  allocate (Gplus_ks(3,nkpt), Gplus_kqs(3,nkpt), Gplus_qs(3,nqpt))
  allocate (list_iGks(nkpt,nG_max))
  allocate (wfc_k_Gs(nkpt,nG_max,nbands,nspin))
  allocate (QE_eig_ks(nkpt,nbands))

  !   
  !================================================================================
  ! Set number of bands. intW is written to read this from W90 data
  ! Allocate array for rotation matrices in the 1bz list
  !================================================================================
  !
  num_bands = nbands
  allocate (u_mesh_1bz(nbands,num_wann,nkpt))



  !   
  !================================================================================
  ! TODO tengo duda sobre el nivel de Fermi en cuanto a unidades en que lee INTW... 
  !================================================================================
  !   
  print *, chemical_potential 

  ! eigenvalues are obtained in eV
  efermi =  14.2522_dp
  ktsmear = 0.2_dp 


  !
  !================================================================================
  ! Read all wavefunctions at irreducible k-mesh
  ! using root node only
  !================================================================================
  !

    call allocate_and_get_all_irreducible_wfc()



  !================================================================================
  ! Get wavefunctions and eigenvalues on full kmesh
  !================================================================================
  !

  do iq = 1,nqpt

    qpoint = qpoints(:,iq)

    ! Locate q-point in 1bz and store point information
    ! iqcorr is index in the 1bz
    call find_k_1BZ_and_G(qpoint,nk1,nk2,nk3,i,j,k,qpoints_1bz(:,iq),Gplus_qs(:,iq))
    call switch_indices_zyx(nk1,nk2,nk3,iqcorr(iq),i,j,k,+1)

    print *, qpoint
    print *, qpoints_1bz(:,iq)
    print *, iqcorr(iq), Gplus_qs(:,iq)

    do ik=1,nkpt  ! k-loop

     print*, "ik= ", ik, " of ", nkpt

     kpoint=kmesh(:,ik)

     ! Locate k-point in 1bz and store point information
     ! (not really necessary, since kmesh was generated in the 1bz)
     call find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i,j,k,kpoint_1bz,G)
     call switch_indices_zyx(nk1,nk2,nk3,ikp,i,j,k,+1)

     kpoints(:,ik) = kpoint
     ikcorr(ik)=ikp
     kpoints_1bz(:,ik)=kpoint_1bz
     Gplus_ks(:,ik)=G

     ! Get k-wavefunction and its list_iG at kpoint
     ! (this is like calling with kpoint and .false. 
     ! because I have moved it to 1bz)
     call get_psi_general_k_all_wfc(.true., kpoint_1bz, list_iGk, wfc_k_G, QE_eig_k, G)
     !
     print*, "nGmax=",nG_max
     print*, "size(list_iGk)=",size(list_iGk)
     print*, "eig_k=", QE_eig_k
     !
     ! Store Bloch state and Wannier rotation matrix of this kpoint.
     ! For this, use index of 1bz, since I want these information 
     ! to be used also by k-q later, and the 1bz-index is a common reference.
     QE_eig_ks(ikp,:) = QE_eig_k
     wfc_k_Gs(ikp,:,:,:) = wfc_k_G
     list_iGks(ikp,:) = list_iGk
     u_mesh_1bz(:,:,ikp)=u_mesh(:,:,ik)
     ! TODO: ojo! posible bug! que pasaría si no lleno estos arrays con 
     ! todos los indices de la 1bz en este loop ik???

     ! k-q point
     ! Locate in 1bz and store point information 
     ! ikqcorr is index in the 1bz
     call find_k_1BZ_and_G(kpoints_1bz(:,ik)-qpoints_1bz(:,iq),nk1,nk2,nk3,i,j,k,kpoint_1bz,G)
     call switch_indices_zyx(nk1,nk2,nk3,ikp,i,j,k,+1)
     ikqcorr(ik)=ikp
     kqpoints_1bz(:,ik)=kpoint_1bz
     Gplus_kqs(:,ik)=G

     print*, 'k-q=', kqpoints_1bz(:,ik)
     print*, Gplus_kqs(:,ik)

    end do ! kpoints
    print *, nkmesh,' wfcs read for qpoint ', iq
   
  end do ! qpoints


  ! not needed any more
  deallocate(wfc_k_G,QE_eig_k)

  !print *, qpoints(:,1)
  !print *, qpoints_1bz(:,1)
  !print *, Gplus_qs(:,1)


  print *, 'all wfc read in'

  !
  !================================================================================
  ! Get one-atom density in reciprocal space using the f-weighted wanniers 
  !================================================================================
  !

  allocate (drho(nqpt,nG_max,nspin,nmode))
  allocate (wfcconv(nbands,nbands,nspin,nspin))


  do iq = 1,nqpt
     write(*,*) ' Calculating drho(iq) ', iq

  do iG = 1,nG_max    ! list for q+G

      if (list_iGks(iqcorr(iq),iG) == 0 ) exit

      print *, iG, list_iGks(iqcorr(iq),iG)

      ! For the transform at q+G (G from the list of q in the 1bz,
      ! Gqplus will be added in the transform to real space)
      Gq = gvec(:, list_iGks(iqcorr(iq),iG) )

      do ik = 1,nkpt

         ! k' = k-q

         ! Each periodic part:
         ! u(r) =  sum_G e^{i G r} wfc(G)
         ! The G_plus vectors are added to the e^-iGr
         Gk=Gplus_ks(:,ik)
         Gkq=Gplus_kqs(:,ik)
         ! < wfc_1 | e^{-i G r} | wfc_2 > = sum_G1 wfc_1(G1)^* wfc_2(G1+G)
         ! Here, wfc(k') is taken for cc. Thus, in current notation:
         ! wfcconv(ibp,ib,isp,is) = < wfc(ikp) | e^{-i Gq r} | wfc(ik) > 

         call get_plane_wave_matrix_element_convolution (Gq+Gkq-Gk, &
                list_iGks(ikqcorr(ik),:),list_iGks(ikcorr(ik),:),&
                wfc_k_Gs(ikqcorr(ik),:,:,:), wfc_k_Gs(ikcorr(ik),:,:,:), wfcconv)

         ! (q+G) factor to be multiplied by mode vector uirr,
         ! so transform to cart coordinates DUDA: at,bg???
         xcrys = real(Gq,dp) + qpoints_1bz(:,iq)
         call cryst_to_cart (1, xcrys, bg, +1) 
 

         do jbnd = 1,nbands
         do ibnd = 1,nbands
          
           ferwei = (1.0_dp+exp((QE_eig_ks(ikcorr(ik),ibnd)-efermi)/ktsmear)) * &
                   (1.0_dp+exp((QE_eig_ks(ikqcorr(ik),jbnd)-efermi)/ktsmear))

           ! same spin selected in u*.u sum to get the total density
           ! TODO think of spinors... am I doing this right for the future?
           ! sum over num_wann index in the equation
           do is=1,nspin
           do imode=1,nmode
             drho(iq,iG,is,imode) = drho(iq,iG,is,imode) + &
                cmplx(0.0_dp,dot_product(uirr(:,imode),xcrys))* &
                sum ( conjg(u_mesh_1bz(jbnd,:,ikqcorr(ik)))*u_mesh_1bz(ibnd,:,ikcorr(ik)) ) &
                * wfcconv(jbnd,ibnd,is,is) /cmplx(ferwei*real(nkpt,dp),0.0_dp)
            end do
            end do

         end do
         end do

       end do   !ik

     print *, drho(iq,iG,:,:)
  end do   !iG
  end do   !iq


  !
  !================================================================================
  ! Transform to real space.
  ! Print out a drho each q vector. 
  !================================================================================
  !

  write(*,*)' Transform to real space'

  allocate (drho_rs(nr1*nr2*nr3,nspin), drho_rs_iGr(nr1*nr2*nr3,nspin,nmode))

  do iq = 1,nqpt

    ph_unit=find_free_unit()
    inquire(iolength=record_length) drho_rs(:,:)
    write(stencil,'(i1)') iq
    filename=trim("wanndrho_q")//trim(adjustl(stencil))//trim(".dat")
    print *, filename
    open( unit=ph_unit, file=filename, iostat=io_err, &
        form='unformatted', status='new', access='direct', recl=record_length )
    do imode=1,nmode

       do is = 1,nspin
        ! For current iq and mode, Fourier transform to real space G-->r and 
        ! store spins, as if it were a wfn
        call  wfc_from_g_to_r (list_iGks(iqcorr(iq),:), drho(iq,:,is,imode), drho_rs(:,is))
       end do

       call r_function_by_exp_igr (-Gplus_qs(:,iq),1, nr1,nr2,nr3, drho_rs, drho_rs_iGr(:,:,imode))

       write(*,*)' Write mode = ', imode

       ! test
       !  do i=1,nr1*nr2*nr3
       !    write(unit=1000*imode,fmt="(2f18.10)") drho_rs_iGr(i,1,imode)
       !  end do
       !end test


    end do   
    close(ph_unit)

         ! test (by the time being, I only use this output test for nspin=1,nq=1)
         do ik=1,nr1*nr2*nr3
           call switch_indices_zyx(nr1,nr2,nr3,ik,i,j,k,-1)
           write(unit=1004,fmt="(3i4,6f18.10)") i,j,k,drho_rs_iGr(ik,1,1:nmode)
         end do
       !end test

  end do ! qpots

  !
  !================================================================================
  ! Finish up
  !================================================================================
  !



  stop
end program wanndenphon

