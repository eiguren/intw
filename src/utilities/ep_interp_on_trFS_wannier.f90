      ! MBR 24/04/2024

      ! This utility reads the electron-phonon matrix elements calculated by the utility
      ! ep_melements.f90 of INTW and interpolates them on a triangulated Fermi surface,
      ! following the method of:
      ! F. Giustino et al, Phys. Rev. B 76, 165108 (2007)
      ! Finally, elements are printed to file.

program ep_on_trFS_wannier

        use kinds, only: dp
    use intw_useful_constants, only: cmplx_1, cmplx_0, cmplx_i, Ha_to_eV, tpi, eps_8
    use intw_utility, only: find_free_unit, find_k_1BZ_and_G, triple_to_joint_index_g, &
            generate_kmesh, cryst_to_cart
    use intw_matrix_vector, only: area_vec
    use intw_w90_setup, only: nrpts, irvec, ndegen, &
            allocate_and_read_ham_r, allocate_and_read_u_mesh, &
            allocate_and_build_ws_irvec, &
            wann_rotate_matrix, wann_IFT_1index, wann_FT_1index_1k, &
            interpolate_1k, &
            eigenval_intw, u_mesh
    use intw_input_parameters, only: nk1, nk2, nk3, mesh_dir, prefix, read_input, &
            nq1, nq2, nq3, ph_dir, chemical_potential, &
            ep_bands, ep_bands_initial, ep_bands_final, &
            use_exclude_bands, &
            ep_interp_method, ep_interp_bands, nfs_sheets_initial, nfs_sheets_final
    use intw_reading, only: read_parameters_data_file_xml, set_num_bands, &
                            nat, nspin, &
                            num_bands_intw, num_wann_intw, &
                            ntyp, ityp, amass, alat, bg, at
    use intw_intw2wannier, only: nnkp_num_kpoints, nnkp_kpoints

    use intw_ph, only: nqmesh, qmesh,  read_ph_information_xml
    use intw_ph_interpolate, only: allocate_and_build_ws_irvec_q, &
            !allocate_and_build_dyn_qmesh, allocate_and_build_dyn_qmesh2, dyn_q_to_dyn_r, &
            irvec_q, nrpts_q, ndegen_q

implicit none

        ! read triangulation
        logical :: read_status, have_nscf , have_wfcN, have_wfc1
        character(5) :: is_loc, ik_loc, comenta
        character(70) :: linea, lleft
        character(100) :: file_off
        character(200) :: comando
        integer :: unit_off
        integer :: nkpt_tr_tot, nkpt_tr_ibz_tot
        integer :: ik1, iface
        integer :: nfs_sheets_tot ! number of sheets considered
        integer :: nfs_sheets_initial_, nfs_sheets_final_ ! those selected sheets
        integer , allocatable :: nfs_sheet(:), & ! band indices of the sheets (num_bands_intw set)
                                 nkpt_tr(:), &  ! number of kpoints in each FS sheet
                                 nkpt_tr_ibz(:), &  ! number of kpoints in each FS sheet irreducible BZ wedge
                                 nface_tr(:)  ! number of faces in each FS sheet
        real(dp) :: k1(3), k2(3), k3(3), kwei
        real(dp), allocatable :: kpts_tr(:,:), kpts_tr_area(:), vk_tr(:,:), vabsk_tr(:)

        ! wannier
    character(256) :: nnkp_file, ep_mat_file, fc_file_name, qband_file_name
    character(4)   :: iq_loc
    integer        :: record_lengh, ierr, ep_unit, nomega, qpt_unit
    integer :: Gkq_1bz(3), is1,is2, iq, ik, ikq, i,j,k,ik_int, iq_int, ir, irq, nq
    integer :: nkmesh, is, js, ikp
    integer :: nq1_int,nq2_int,nq3_int, nk1_int,nk2_int,nk3_int, nkmesh_int, nqmesh_int, &
            iomega, iw,jw, imode, jmode, iat, nint_k, nint_q, nener
    integer , allocatable :: kqmap(:,:), kqmap_int(:,:)
    real(dp) :: kpoint(3), qpoint(3), kqpoint(3), kq_1bz(3), qpoint_cart(3), rcart(3)
    real(dp), allocatable :: kmesh(:,:), kqmesh(:,:,:), kmesh_int(:,:), qmesh_int(:,:), kqmesh_int(:,:,:)
    real(dp) , allocatable :: eig_kint(:), eig_kqint(:), eig_kint_all(:,:)
    real(dp) :: ener
    complex(dp) :: fac, facq, fack, gep_ij, cdum
    complex(dp) , allocatable ::  u_kint(:,:), u_kqint(:,:), u_kint_all(:,:,:)
    complex(dp), allocatable :: gmatkqk_wann(:,:,:,:), gmatL_wann(:,:,:,:,:,:,:), gmat_aux(:,:,:,:), &
            gmat_int(:,:), gmat_int_rot(:,:), gmat_aux1(:,:,:), gep_int(:,:)
    complex(dp), allocatable :: ep_mat_el_coarse(:,:,:,:,:,:,:)

    real(dp), external :: smeared_delta, fermi_dirac

         ! elements
         logical                  :: have_ep
        character(256) :: altprefix, file_ep
        integer :: ib,jb, unit_ep, ibp
        integer :: ipb, iksp, iks, ish, ishp, ir1,ir2,ir3
        real(dp) :: kpoint_p(3)
        complex(dp), allocatable :: aep_mat_el(:,:,:,:,:)


    ! From mat_inv_four_t, see IGG's comments
      ! TODO add to useful_constants
      real(dp), parameter :: pmass=1822.88848426_dp, aumev=  27211.396d0


      !--------------------------------------------------------------------------------
!================================================================================
!       Talk to the user
!================================================================================
write(*,'(A)') '====================================================='
write(*,'(A)') '|     Wannier interpolation of e-p elements         |'
write(*,'(A)') '|               on thetriangulated FS               |'
write(*,'(A)') '|        ---------------------------------          |'
write(*,'(A)') '====================================================='
write(*,'(A)') '|    waiting for input file...                      |'

!================================================================================
!       read the input files
!================================================================================
    call read_input(read_status)

    if (read_status ) then
       stop
    end if

    if ( use_exclude_bands .ne. "wannier") then
       write(*,*) ' use_exclude_bands /= wannier in input. Stopping.'
       stop
    end if

    if ( ep_interp_method .ne. 'wannier') then
       write(*,*) ' ep_interp_method /= wannier in input. Stopping.'
       stop
    end if

    !  Establish number of bands
    ! this needs to be invoked before set_num_bands
    call read_parameters_data_file_xml()
    ! this will read num_wann_intw and num_bands_intw dimensions from .nnkp if available
    call set_num_bands()
    write(*,*)' num_bands_intw = ', num_bands_intw

    ! choose Fermi surface sheets according to ep_interp_bands
    if ( ep_interp_bands == 'intw_bands' ) then
       nfs_sheets_tot = num_bands_intw
       nfs_sheets_initial_ = 1
       nfs_sheets_final_ = num_bands_intw
       allocate ( nfs_sheet(nfs_sheets_tot) )
       do ib = 1, num_bands_intw
           nfs_sheet(ib) = ib
       end do
    else if ( ep_interp_bands == 'ef_crossing' ) then
       nfs_sheets_tot = nfs_sheets_final - nfs_sheets_initial + 1
       nfs_sheets_initial_ = nfs_sheets_initial
       nfs_sheets_final_ = nfs_sheets_final
       allocate ( nfs_sheet(nfs_sheets_tot) )
       do ib = 1, nfs_sheets_tot
           nfs_sheet(ib) = nfs_sheets_initial + ib-1
       end do
    end if


    ! Control: for the Wannier interpolation we need all coarse ep elements in the
    ! ep_bands == 'intw' set (otherwise we cannot rotate with the U matrices).
    ! Therefore, check if we will have all the elements, even though custom_bands
    ! was set.

    if (ep_bands == 'custom' .and. &
       ( ep_bands_initial .ne. 1 .or. ep_bands_final .ne. num_bands_intw) ) then
       write(*,*)' Coarse ep elements where calculated in the', &
               ep_bands_initial, ' to', ep_bands_final,' band subset only.'
       write(*,*)' This means that we cannot Wannier interpolate. Stopping.'
         stop
    end if



!================================================================================
!       read u_mesh file from w902intw
!================================================================================

    call allocate_and_read_u_mesh()

!================================================================================
!       read ham_r file from w902intw, to be used in band interpolation
!       (theis also allocates and reads the irvec list)
!================================================================================

    call allocate_and_read_ham_r()

!================================================================================
!       Read all the information about phonons and pseudopotentials
!================================================================================

  call read_ph_information_xml()
!================================================================================
!   generate coarse meshes
!================================================================================

    nqmesh = nq1*nq2*nq3
    nkmesh = nk1*nk2*nk3

    allocate(qmesh(3,nqmesh))
    call generate_kmesh(qmesh,nq1,nq2,nq3)

    allocate(kmesh(3,nkmesh))
    call generate_kmesh(kmesh,nk1,nk2,nk3)

    ! DUDA kmesh is equal to nnkp_kpoints, but I do not if this is general
    ! si no, hay que usar kmesh explicitamente en las (I)FTs
    !do ik=1,nkmesh
    !   write(*,'(i3,3f7.3)') ik, kmesh(:,ik)-nnkp_kpoints(:,ik)
    !end do

    ! k+q mesh and index correspondence
    allocate(kqmesh(3,nqmesh,nkmesh))
    allocate(kqmap(nqmesh,nkmesh))
    do iq=1,nqmesh
            qpoint = qmesh(:,iq)
        do ik=1,nkmesh
            kpoint = kmesh(:,ik) ! nnkp_kpoints(:,ik) is equivalent
            kqpoint = kpoint+qpoint
            kqmesh(:,iq,ik) = kqpoint
            ! elements k+q,k are: ep_mat_el_coarse(iq,ik,:,:,is1,is2,imode)
            ! locate (k+q)-point index ikq in the kmesh
            call find_k_1BZ_and_G(kqpoint,nk1,nk2,nk3,i,j,k,kq_1bz,Gkq_1bz)
            call triple_to_joint_index_g(nk1,nk2,nk3,ikq,i,j,k)
            kqmap(iq,ik) = ikq
        end do
    end do

!================================================================================
! Read the force constant matrix from this file.
! Calculate dynamical matrices in qmesh and transform to real space
!================================================================================

  ! get the WS vectors for the coarse q-mesh Wannier interpolation
  call allocate_and_build_ws_irvec_q()

!================================================================================
!   open and read ep files
!================================================================================

    allocate(ep_mat_el_coarse(nqmesh,nkmesh,num_bands_intw,num_bands_intw,nspin,nspin,3*nat))
    inquire(iolength=record_lengh) ep_mat_el_coarse(1,:,:,:,:,:,:)

    write(*,*) 'reading ep_mat.dat files...'
    do iq = 1,nqmesh
      ep_unit = find_free_unit()
      if (                iq <   10) write(iq_loc,"(i1)")iq
      if ( 10 <= iq .and. iq <  100) write(iq_loc,"(i2)")iq
      if (100 <= iq .and. iq < 1000) write(iq_loc,"(i3)")iq
      open(unit=ep_unit,file=trim(ep_mat_file)//trim('_')//adjustl(iq_loc), &
          status='old',form='unformatted',access='direct',recl=record_lengh)
      read(unit=ep_unit, rec=1, iostat=ierr) ep_mat_el_coarse(iq,:,:,:,:,:,:)
      close(unit=ep_unit)
    end do
    write(*,*) '... done'

!================================================================================
!  interpolation meshes: given by the triangulated FS
!       Read prefix.off and velocity file(s)
!================================================================================


write(*,'(A)') '====================================================='
write(*,'(A)') '|    reading .off and v_k files...                  |'

allocate (nkpt_tr(nfs_sheets_tot), nface_tr(nfs_sheets_tot))
allocate (nkpt_tr_ibz(nfs_sheets_tot))

!open all sheet files just to see dimensions of kpoint lists
unit_off=find_free_unit()

do is=1,nfs_sheets_tot

   if (                is <   10) write(is_loc,"(i1)") nfs_sheet(is)
   if ( 10 <= is .and. is <  100) write(is_loc,"(i2)") nfs_sheet(is)

   file_off=trim(mesh_dir)//trim(prefix)//trim('.')//trim(adjustl(is_loc))//trim('_FS_tri.off')
   write(*,*) is, file_off

   open(unit_off, file=file_off, status='old')
   read(unit_off,*) comenta
   read(unit_off,*) nkpt_tr(is), nface_tr(is)  !number of vertices and faces (ignore edges)
   write(*,*) 'sheet, vertices, faces = ', nfs_sheet(is), nkpt_tr(is), nface_tr(is)
   close(unit_off)

   !open the IBZ off file and search for dimension nkpt_tr_ibz(is).
   !Its vertices coincide with the first nkpt_tr_ibz(is) vertices of the full off vertex list.

   file_off=trim(mesh_dir)//trim(prefix)//trim('.')//trim(adjustl(is_loc))//trim('_newton_IBZ_FS_tri.off')
   write(*,*) is, file_off

   open(unit_off, file=file_off, status='old')
   read(unit_off,*) comenta
   read(unit_off,*) nkpt_tr_ibz(is)  !number of vertices (ignore rest)
   write(*,*) 'sheet, vertices = ', nfs_sheet(is), nkpt_tr_ibz(is)
   close(unit_off)

end do

! total number of k-points to be calculated
nkpt_tr_tot=sum(nkpt_tr)
nkpt_tr_ibz_tot=sum(nkpt_tr_ibz)

write(*,'(A)') '|   Number of k-points (total vertices):            |'
write(*,*)  nkpt_tr_tot
write(*,'(A)') '|   Number of k-points in IBZ (total vertices):     |'
write(*,*)  nkpt_tr_ibz_tot
write(*,'(A)') '|   Number of total faces:                          |'
write(*,*) sum(nface_tr)


allocate(kpts_tr(3,nkpt_tr_tot), kpts_tr_area(nkpt_tr_tot))
allocate(vk_tr(3,nkpt_tr_tot), vabsk_tr(nkpt_tr_tot) )

!open .off's again, read k-points and read velocity files
ik1 = 0
do is=1,nfs_sheets_tot

   if (                is <   10) write(is_loc,"(i1)") nfs_sheet(is)
   if ( 10 <= is .and. is <  100) write(is_loc,"(i2)") nfs_sheet(is)

   ! .off file for this sheet

   file_off=trim(mesh_dir)//trim(prefix)//trim('.')//trim(adjustl(is_loc))//trim('_FS_tri.off')
   open(unit_off, file=file_off, status='old')

   read(unit_off,*) comenta
   read(unit_off,*) i,j,k  !number vertices, faces and edges (I will ignore edges)
   !read(unit_off,'(/)') !DUDA... esto dependdera de como está escrito el salto de línea en el fichero, creo...
   if  ( i .ne. nkpt_tr(is) .or. j .ne. nface_tr(is) ) then
           write(*,*)'Error reading ', file_off, '. Stopping.'
           stop
   end if

      !read vertices
   do ik=1,nkpt_tr(is)
      read(unit_off,*) kpts_tr(:,ik1+ik)  !units in the trFS.off file are cartesian 2pi/alat ("tpiba" for QE)
   end do

   !Read (triangular) faces on this sheet.
   !Each face contributes with 1/3 of its area to the effective area of each of its vertices.
   !Calculate the are on the go and add the contribution to each vertex, storing for global indices (i.e. ik1+ik).
   do iface = 1,nface_tr(is)
      read(unit_off,*) i, ir1,ir2,ir3 !indices ik of the vertices of the face, indexed from 0
      ir1=ir1+1
      ir2=ir2+1
      ir3=ir3+1  !now, ik of the vertices of the face, indexed from 1
      if  ( i .ne. 3 ) then
              write(*,*)'Error reading ', file_off, 'Only triangles allowed. Stopping.'
           stop
      end if
      !triangle vertex vectors (cartesian 2pi/alat)
      k1 = kpts_tr(:,ik1+ir1)
      k2 = kpts_tr(:,ik1+ir2)
      k3 = kpts_tr(:,ik1+ir3)
      !get spanned area and add contribution to each vertex
      !function copied from FSH/modules/geometry.f90
      kwei = area_vec(k2-k1,k3-k1)/3.0_dp
      kpts_tr_area(ik1+ir1) = kpts_tr_area(ik1+ir1) + kwei
      kpts_tr_area(ik1+ir2) = kpts_tr_area(ik1+ir2) + kwei
      kpts_tr_area(ik1+ir3) = kpts_tr_area(ik1+ir3) + kwei
   end do

   close(unit_off)

      ! velocity for this sheet (use same unit)

   file_off=trim(mesh_dir)//trim(prefix)//trim('.')//trim(adjustl(is_loc))//trim('_FS_v_k.dat')
   open(unit_off, file=file_off, status='old')

   read(unit_off,*) i
   if  ( i .ne. nkpt_tr(is) ) then
           write(*,*)'Error reading ', file_off, '. Stopping.'
           stop
   end if
   do ik=1,nkpt_tr(is)
      read(unit_off,*) i, vk_tr(:,ik1+ik), vabsk_tr(ik1+ik)  !velocity xyz and its modulus. DUDA 2pi/alat???
   end do

   close(unit_off)

   !accumulate ik global index for the reading of next sheet
   ik1=ik1+nkpt_tr(is)

end do

write(*,'(A)') '|    .... done reading                              |'
write(*,'(A)') '|================================================== |'

print *,  kpts_tr_area(1), vabsk_tr(1)
print *,  kpts_tr_area(nkpt_tr_tot), vabsk_tr(nkpt_tr_tot)
print *, 'total area FS = ', sum(kpts_tr_area) ,' (2 x pi / alat)^2 ' !OK



!================================================================================
!   ep elements in real space
!   For each nspin,nspin,3*nat:
!      1. rotate matrix g(k+q,k,m,n) to wannier gauge with u_mesh
!      2. For k+q and k indices, inverse Fourier index by index.
!         Note, k and q are fouriered on different grids.
!================================================================================

    write(*,*) 'nrpts = ', nrpts

    allocate (gmatkqk_wann(num_wann_intw,num_wann_intw,nqmesh,nkmesh))
    allocate (gmat_aux(num_wann_intw,num_wann_intw,nqmesh,nrpts) )
    allocate (gmatL_wann(3*nat,num_wann_intw,num_wann_intw,nspin,nspin,nrpts_q,nrpts) )

    gmatL_wann=cmplx_0

    do is1=1,nspin
    do is2=1,nspin

    do iat=1,3*nat

        k=((iat-1)/3)+1  !atom index

        ! 1.
        ! rotate U^dagger(k+q) * gmat * U(k)
        gmatkqk_wann=cmplx_0
        do iq=1,nqmesh
           do ik=1,nkmesh
              ikq = kqmap(iq,ik)
              call wann_rotate_matrix(ikq,ik, ep_mat_el_coarse(iq,ik,:,:,is1,is2,iat), &
                                      gmatkqk_wann(:,:,iq,ik))
           end do
        end do

        !2.
        gmat_aux=cmplx_0
        !
        ! Inverse Fourier over k-index, using kmesh and irvec grids
        do iq=1,nqmesh
           call wann_IFT_1index (nkmesh, kmesh, gmatkqk_wann(:,:,iq,:), gmat_aux(:,:,iq,:) )
        end do

        !
        ! Inverse Fourier over q-index, using qmesh and irvec_q grids
        ! (this is made "by hand")
        do iq=1,nqmesh
        do irq=1,nrpts_q
           facq = exp(-cmplx_i*tpi*dot_product(qmesh(:,iq),irvec_q(:,irq)))/real(nqmesh,dp)
           do ir=1,nrpts
              gmatL_wann(iat,:,:,is1,is2,irq,ir) = &
                  gmatL_wann(iat,:,:,is1,is2,irq,ir) + gmat_aux(:,:,iq,ir) * facq
           end do
        end do
        end do
        write(*,*) 'IFT on WS done for displacement', iat

    end do !iat

    end do !spin
    end do !spin


!================================================================================
!  Interpolation variables
!================================================================================

! ep arrays
allocate (gmat_aux1(num_wann_intw,num_wann_intw,nrpts) )
allocate (gmat_int(num_wann_intw,num_wann_intw) )
allocate (gmat_int_rot(num_wann_intw,num_wann_intw) )
allocate (gep_int(num_wann_intw,num_wann_intw) )
allocate(aep_mat_el(nkpt_tr_tot,nkpt_tr_ibz_tot,nspin,nspin,3*nat))

!band arrays
allocate (eig_kint(num_wann_intw), eig_kqint(num_wann_intw))
allocate (u_kint(num_wann_intw,num_wann_intw), u_kqint(num_wann_intw,num_wann_intw))


!================================================================================
! Interpolate elements
!================================================================================

 ! the interpolation grid is formed by the triangulated FS points
 ! In principle, in the trianglulated k-point list up to nkpt_tr_tot
 ! the index carries the information about the band implicitly.
 ! We will retrieve U rotation matrices and eigenvalues on the
 ! full num_wann_intw list for those k-points, as the w90_setup
 ! routines give that info, but then we will intepolate the ep element only
 ! for the bands index at Fermi given by nfs_sheet().


 ! interpolated bands in nkpt_tr_tot vertices
 ! u_kint_all, contains the rotation matrices (row is bloch and column is wannier,
 ! i.e. the "orbital weights" for each band is in the columns)
 ! Store frequencies in Ry.

 allocate (u_kint_all(num_wann_intw,num_wann_intw,nkpt_tr_tot))
 allocate ( eig_kint_all(num_wann_intw,nkpt_tr_tot))

 do ik=1,nkpt_tr_tot
     ! this is cartesians x 2pi/alat. Transform to crystal before interpolation
     kpoint = kpts_tr(:,ik)
     call cryst_to_cart (1, kpoint, at, -1)
     call interpolate_1k (kpoint, eig_kint, u_kint)
     eig_kint = eig_kint * 2.0_dp / Ha_to_eV
     eig_kint_all(:,ik) = eig_kint
     u_kint_all(:,:,ik) = u_kint
 end do

 !================================================================================
 ! Loop over k'=k+q, k electrons. Get interpolated ep elements and write out
 !================================================================================

  unit_ep=find_free_unit()
  file_ep=trim(mesh_dir)//trim(prefix)//trim('_ep_interp.dat')

  inquire(file=file_ep,exist=have_ep)
  print *, 'have_ep?', have_ep

  if (.not.have_ep) then ! calculate interpolated ep elements and write to file_ep

  write(*,*) '!      e-p elements to be interpolated and written to file:   !'
  write(*,*) file_ep

  open(unit_ep, file=file_ep, status='unknown')
  write(unit_ep,*)'# ik(irr)   jk(full)    is js   g(canonical modes)'

  ! ik, ikp indices implicitly contain the FS sheet index, i.e. the band indices ib, ib'
! to be selected, so instead of iterating over nkpt_tr_tot, I separate over sheets
! (aep_mat_el elements are stored only for the needed pair of sheets for a given kk')

  ik = 0  ! kpoint
  ik1 = 0
  do ish = 1, nfs_sheets_tot
      ib = nfs_sheet(ish) !band index for k
  do iks = 1, nkpt_tr_ibz(ish)
      ! ik is the k-index over nkpt_tr_tot in the Irreducible BZ
      ! ik1 is the corresponding k-index in the full BZ kpts_tr list (at the end of the loop I add the rest of nkpt_tr(ish))
      ik = ik + 1
      ik1 = ik1 + 1

      kpoint = kpts_tr(:,ik1)  ! this is cartesians x 2pi/alat. Transform to cryst.
      call cryst_to_cart (1, kpoint, at, -1)

      ! interpolated bands in k and k+q (in Ry) calculated above
      ! u_kint, u_kqint contain the rotation matrices (row is band and column is wannier,
      ! i.e. the "orbital weights" for each band is in the columns)
      u_kint = transpose(conjg(u_kint_all(:,:,ik1)))  ! we will need "dagger" below for U(k+q) g U^dagger(k)

  ikp = 0 ! kpoint+qpoint
  do ishp = 1, nfs_sheets_tot
      ibp = nfs_sheet(ishp) !band index for k'
  do iksp = 1, nkpt_tr(ishp)
      ikp = ikp + 1   !k'-index over nkpt_tr_tot

      kpoint_p = kpts_tr(:,ikp)  ! this is cartesians x 2pi/alat. Transform to cryst.
      call cryst_to_cart (1, kpoint_p, at, -1)

      u_kqint = u_kint_all(:,:,ikp)

      qpoint = kpoint_p-kpoint

      ! In the loop over atoms and directions below, I will
      ! interpolate matrix elements for each displacement.
      ! 1. interpolate irvec_q on this qpoint and store in gmat_aux1(:,:,nrpts)
      ! 2. interpolate irvec on kpoint and store in gmat_int.
      ! 3. Rotate Wannier as gmat_int_rot = U(k+q) * gmat_int * U(k)^dagger
      ! 4. The ibp,ib elements are the ep element needed (in canonical atom displacement coordinates)

      ! TODO some of these can be done out of the 3nat and spin loops to go faster
      do js=1,nspin
      do is=1,nspin

        do iat=1,3*nat
           gmat_aux1 = cmplx_0
           gmat_int = cmplx_0
           gmat_int_rot = cmplx_0
           do irq=1,nrpts_q
              facq = exp(cmplx_i*tpi*dot_product(qpoint(:),irvec_q(:,irq)))/real(ndegen_q(irq),dp)
              ! TODO DUDA the spin order in gmatL_wann is correct???
              gmat_aux1(:,:,:) = gmat_aux1(:,:,:) + facq*gmatL_wann(iat,:,:,js,is,irq,:)
           end do
           call wann_FT_1index_1k (kpoint, gmat_aux1(:,:,:), gmat_int(:,:))
           gmat_int_rot(:,:) = matmul(u_kqint, matmul(gmat_int, u_kint))
           aep_mat_el(ikp,ik,js,is,iat) = gmat_int_rot(ibp,ib)  !TODO DUDA exclude bands ?? (ver junto al comentario JLB en w90_setup)
        end do

!!        write(unit_ep,fmt="(2i6,3f10.4,100e16.6)") ikp,ik, kpts_tr(:,ikp),  &
          write(unit_ep,fmt="(6i6,100e16.6)") ibp, iksp, ikp, ib, iks, ik,  &
              (aep_mat_el(ikp,ik, js,is,iat), iat=1,3*nat)

      end do
      end do !spins

   end do  ! k'
  end do  ! sheet'

  print *, ik, ik1
  end do  ! k
  ik1 = ik1 + nkpt_tr(ish) - nkpt_tr_ibz(ish)
  end do  ! sheet

  close(unit_ep)

  else ! have_ep--> read interpolated matrix elements

  write(*,*) '!     interpolated e-p file already exists. Check!            !'
  write(*,*) file_ep

  end if

write(*,'(A)') '|================================================== |'
write(*,'(A)') '|               ALL DONE                            |'
write(*,'(A)') '|================================================== |'

stop

end program ep_on_trFS_wannier
