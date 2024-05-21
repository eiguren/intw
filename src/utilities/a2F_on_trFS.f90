      ! MBR 26/04/2024

      ! Calculate Eliashberg function from electron-phonon coupling
      ! matrix elements interpolated on a triangulated Fermi surface.

      ! a2F integral = double loop on kpoints of the FS triangulation.
      ! The ep element calculated for q=k'-k is read from a file.
      ! The needed phonons and  dynamical matrices calculated by QE are
      ! read in and then interpolated as in the method of:
      ! F. Giustino et al, Phys. Rev. B 76, 165108 (2007)


program a2F_on_trFS

        use kinds, only: dp    
        use intw_useful_constants, only: cmplx_1, cmplx_0, cmplx_i, Ha_to_eV, tpi, eps_8
        use intw_utility, only: find_free_unit, cryst_to_cart, area_vec, smeared_delta, &
                generate_kmesh

        use intw_input_parameters, only: mesh_dir, prefix, read_input,&
                nk1, nk2, nk3, nq1, nq2, nq3, nqirr, ph_dir, & !, fc_mat
                nomega, omega_ini, omega_fin, osmear_q, &
                ep_interp_bands, nfs_sheets_initial, nfs_sheets_final
        use intw_reading, only: read_parameters_data_file_xml, set_num_bands, &
                nspin, at, bg, volume0, alat, nat, ntyp, ityp, amass, &
                num_bands_intw
        use intw_ph, only: nqmesh, qmesh, frc, readfc, read_ph_information_xml, &
                q_irr, q_irr_cryst

        use intw_ph_interpolate, only: allocate_and_build_dyn_qmesh2,  dyn_q_to_dyn_r, &
                allocate_and_build_ws_irvec_q, dyn_interp_1q, dyn_diagonalize_1q, &
                w2_q, u_q

    implicit none

        logical :: read_status
        character(5) :: is_loc, ik_loc, comenta
        character(70) :: linea, lleft
        character(100) :: file_off, dir_scf, dir_nscf, file_a2f
        character(200) :: comando
        integer :: unit_off, unit_a2f
        integer :: nkpt_tr_tot, nkpt_tr_ibz_tot
        integer :: is, js, ikpt, ik, ik1, i,j,k, iface, ir, ir1, ir2, ir3
        integer  :: nfs_sheets_tot ! number of sheets considered
        integer , allocatable :: nfs_sheet(:), & ! band indices of the sheets (num_bands_intw set)
                                 nkpt_tr(:), &  ! number of kpoints in each FS sheet
                                 nkpt_tr_ibz(:), &  ! number of kpoints in each FS sheet irreducible BZ wedge
                                 nface_tr(:)  ! number of faces in each FS sheet
        real(dp) :: k1(3), k2(3), k3(3), kwei, vol1bz
        real(dp), allocatable :: kpts_tr(:,:), kpts_tr_area(:), vk_tr(:,:), vabsk_tr(:)
        real(dp) , allocatable :: area_ibz(:), area_fbz(:), factor_area_ibz(:)

        logical                  :: have_ep
        character(256) :: altprefix, file_ep
        integer :: ib,jb,iat,imode,ikp, iG, unit_ep

        real(dp) , allocatable :: w_q(:,:)

        integer :: nener, iomega, iks, ish,  iksp, ishp, ibp
        real(dp) :: omega, omega_step, rfacq, dosef, dsk_vk_2
        real(dp) :: qpoint(3), kpoint(3), kpoint_p(3)
        real(dp) , allocatable :: alpha2F(:,:,:,:), w2_qint(:), w_qint(:), lambda(:)
        complex(dp), allocatable :: gep_int(:)
        complex(dp), allocatable :: dyn_qint(:,:), u_qint(:,:)
        complex(dp), allocatable :: aep_mat_el(:,:,:,:,:)

        ! From mat_inv_four_t, see IGG's comments
        ! TODO add to useful_constants
        real(dp), parameter :: pmass=1822.88848426_dp, aumev=  27211.396d0


        20 format(A)
        30 format(A,F8.2,6X,A)

!--------------------------------------------------------------------------------
!================================================================================
!       Talk to the user
!================================================================================
write(*,'(A)') '====================================================='
write(*,'(A)') '|         Eliashberg function calculation           |'
write(*,'(A)') '|    from interpolated ep elements on the FS        |'
write(*,'(A)') '|        ---------------------------------          |'
write(*,'(A)') '====================================================='
write(*,'(A)') '|    waiting for input file...                      |'

!================================================================================
!       read the input file
!       Read in the necessary information from standard input 
!       (for this utility, only prefix needed, really)
!================================================================================

    call read_input(read_status)

    if (read_status ) then
       stop
    end if

    ! read the parameters file (for nat, ntyp)
    call read_parameters_data_file_xml()

    ! choose Fermi surface sheets according to ep_interp_bands
    if ( ep_interp_bands == 'intw_bands' ) then
       nfs_sheets_tot = num_bands_intw
       allocate ( nfs_sheet(nfs_sheets_tot) )
       do ib = 1, num_bands_intw
           nfs_sheet(ib) = ib
       end do
    else if ( ep_interp_bands == 'ef_crossing' ) then
       nfs_sheets_tot = nfs_sheets_final - nfs_sheets_initial + 1
       allocate ( nfs_sheet(nfs_sheets_tot) )
       do ib = 1, nfs_sheets_tot
           nfs_sheet(ib) = nfs_sheets_initial + ib-1
       end do
    end if

!================================================================================
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

   open(unit_off, file=file_off, status='old')
   read(unit_off,*) comenta
   read(unit_off,*) nkpt_tr(is), nface_tr(is)  !number of vertices and faces (ignore edges) 
   close(unit_off)

   !open the IBZ off file and search for dimension nkpt_tr_ibz(is). 
   !Its vertices coincide with the first nkpt_tr_ibz(is) vertices of the full off vertex list.

   file_off=trim(mesh_dir)//trim(prefix)//trim('.')//trim(adjustl(is_loc))//trim('_newton_IBZ_FS_tri.off')

   open(unit_off, file=file_off, status='old')
   read(unit_off,*) comenta
   read(unit_off,*) nkpt_tr_ibz(is)  !number of vertices (ignore rest)
   close(unit_off)

end do

! total number of k-points to be calculated
nkpt_tr_tot=sum(nkpt_tr)
nkpt_tr_ibz_tot=sum(nkpt_tr_ibz)

write(*,'(A)') '|   .off file                                       |'
write(*,*) is, file_off
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


! Calculate triangle areas for calculating the a2F integral later
  !N(EF) from sum over vertices, all sheets, using full triangulated mesh
  write(*,'(A)') '|   DOS at Fermi:                                   |'
  dosef = 0.0_dp
  do ik = 1,nkpt_tr_tot
     dosef = dosef + kpts_tr_area(ik) * (tpi/alat)**2 / vabsk_tr(ik)
  end do
  vol1bz = tpi**3 / volume0
  dosef = 2.0_dp * dosef / vol1bz
  write(*,*) 'Volume of BZ (bohr**-3) =', vol1bz
  write(*,*) 'DOS at FS = ', dosef, 'a.u.'

   !Ratio of areas of the FS in the full and irreducible BZ
  allocate(area_fbz(nfs_sheets_tot), area_ibz(nfs_sheets_tot), factor_area_ibz(nfs_sheets_tot))
  area_ibz=0.0_dp
  area_fbz=0.0_dp
  ik = 0
  do ish = 1, nfs_sheets_tot
  do iks = 1, nkpt_tr_ibz(ish)
      ik = ik + 1
      area_ibz(ish) = area_ibz(ish) + kpts_tr_area(ik)
  end do
  end do
  ik = 0
  do ish = 1, nfs_sheets_tot
  do iks = 1, nkpt_tr(ish)
      ik = ik + 1
      area_fbz(ish) = area_fbz(ish) + kpts_tr_area(ik)
  end do
  end do
  do ish=1,nfs_sheets_tot
     factor_area_ibz(ish) = area_fbz(ish)/area_ibz(ish)
  end do
  write(*,*) 'factor_area_ibz = ', factor_area_ibz(:)


!================================================================================
! Read ep elements file
!================================================================================

  unit_ep=find_free_unit()
  file_ep=trim(mesh_dir)//trim(prefix)//trim('_ep_interp.dat')
  inquire(file=file_ep,exist=have_ep)

  if (.not.have_ep) then

     write(*,*) '!      e-p elements interpolated and written to file:         !'
     write(*,*) file_ep
     write(*,*) '!      e-p elements file does not exist:                      !'
     write(*,*) '!      Stopping                                               !'
     stop

  else       

    ! Note, indices ikp,ik include all the calculated sheets
    allocate(aep_mat_el(nkpt_tr_tot,nkpt_tr_ibz_tot,nspin,nspin,3*nat))

    open(unit_ep, file=file_ep, status='old')
    read(unit_ep,*) comenta

    do i = 1,nkpt_tr_tot*nspin
    do j = 1,nkpt_tr_ibz_tot*nspin
        do js=1,nspin
        do is=1,nspin
           read(unit_ep,fmt="(2i6,3f10.4,100e16.6)") ikp,ik, kpoint, &
                 (aep_mat_el(ikp,ik, js,is,iat), iat=1,3*nat)
        end do
        end do 
    end do
    end do        

    close(unit_ep)

    write(*,*) '!     interpolated e-p elements read from file:               !'
    write(*,*) file_ep         

  end if  

!================================================================================
! Read all the information about phonons and prepare for interpolation of phonons:
! Calculate dynamical matrices in qmesh and transform to real space
!================================================================================

  ! Generate full zone coarse qmesh
  nqmesh = nq1*nq2*nq3
  allocate(qmesh(3,nqmesh))
  call generate_kmesh(qmesh,nq1,nq2,nq3)
  ! This Wigner-Seitz mesh will be needed to interpolate w_q below
  ! (q->R space will be done on that grid of nrpts_q lattice points).
  call allocate_and_build_ws_irvec_q()

  ! q_irr vectors info is obtained through this (it is contained in qlist.txt file)
  call read_ph_information_xml()

  ! INTW has two options to get dyn_q:
  ! 1. read force constants
  !  fc_file_name = trim(mesh_dir)//trim(ph_dir)//trim(fc_mat)
  !  call allocate_and_build_dyn_qmesh(fc_file_name)
  ! or:
  ! 2. read dyn files directly (new by MBR, not fully tested)
  call allocate_and_build_dyn_qmesh2()

  ! transform to R space
  call dyn_q_to_dyn_r()

  ! Now w2_q contains phonon frequencies and u_q contains the polarization vectors.
  ! Phonon frequencies (negative means unstable mode)
  allocate (w_q(3*nat,nqmesh) )
  w_q=sign(sqrt(abs(w2_q)),w2_q) * 0.002_dp / Ha_to_eV !meV to Ry


!  Energy range for Eliashberg is read from intw.in
!  with Ry units throughout.
!  Step in energy for integrals and phonon DOS:
  omega_step = (omega_fin-omega_ini)/real(nomega-1,dp)

  write(*,'(A)') '|    .... done reading                              |'
  write(*,'(A)') '|================================================== |'


!================================================================================
! Calculate a2F integral
! TODO duda spin index in a2F??
!================================================================================


  !phonon arrays
  allocate( dyn_qint(3*nat,3*nat), u_qint(3*nat,3*nat), w2_qint(3*nat), w_qint(3*nat) )
  ! ep elements
  allocate (gep_int(3*nat) )
  ! integrals
  allocate (alpha2F(nspin,nspin,3*nat,nomega), lambda(3*nat))

  ! output file for a2F
  unit_a2f=find_free_unit()
  file_a2f=trim(mesh_dir)//trim(prefix)//trim('_a2F_interp.dat')
  open(unit_a2f, file=file_a2f, status='unknown')
  write(unit_a2f,'(A)') '#omega(Ry)  alpha2F(total)    alpha2F(1:nmode)'

  alpha2F = 0.0_dp
  lambda = 0.0_dp

  write(*,'(A)') '!     Calculating a2F...                            !'

! ik, ikp indices implicitly contain the FS sheet index, i.e. the band indices ib, ib'
! to be selected, so instead of iterating over nkpt_tr_tot, I separate over sheets

ik = 0
do ish = 1, nfs_sheets_tot
      ib = nfs_sheet(ish) !band index for k
do iks = 1, nkpt_tr_ibz(ish)
      ik = ik + 1   !k-index over nkpt_tr_tot in the Irreducible BZ

      kpoint = kpts_tr(:,ik)  ! this is cartesians x 2pi/alat. Transform to cryst.
      call cryst_to_cart (1, kpoint, at, -1)

  ikp = 0
  do ishp = 1, nfs_sheets_tot
      ibp = nfs_sheet(ishp) !band index for k'
  do iksp = 1, nkpt_tr(ishp)
      ikp = ikp + 1   !k'-index over nkpt_tr_tot

      kpoint_p = kpts_tr(:,ikp)  ! this is cartesians x 2pi/alat. Transform to cryst.
      call cryst_to_cart (1, kpoint_p, at, -1)

      qpoint = kpoint_p-kpoint

      !interpolate phonon:
      call dyn_interp_1q(qpoint, dyn_qint)
      call dyn_diagonalize_1q(3*nat, dyn_qint, u_qint, w2_qint)
      w_qint=sign(sqrt(abs(w2_qint)),w2_qint)
      !u_qint contains the polarization vector (each column is a mode)
      !phonon frequency in meV units: pass to Ry
      w_qint= w_qint*0.002_dp / Ha_to_eV

      !Generate ep elements with the added mass factor and polarization vector
      !component for each mode and sum over modes.
      !In the eigenvector matrix u_qint, each column is a mode. Rows are individual atom displacements.
      !The matrix elements g are given for atom displacements and I want them as a
      !function of modes to compare with QE a2F.dos* files

      ! band indices to be used are the FS sheet indices given at the beginning of the loop.
      ! I select only those ones for the interpolated elements

      do is=1,nspin
      do js=1,nspin

         gep_int = cmplx_0
         do imode = 1,3*nat
            if (w_qint(imode)>eps_8) then !stable mode frequency
               do iat=1,3*nat
                  k=((iat-1)/3)+1  !atom index
                  rfacq=2.0_dp * w_qint(imode) * amass(ityp(k)) * pmass
                  gep_int(imode) = gep_int(imode)  +  &
                       aep_mat_el(ikp,ik,is,js,iat) * u_qint(iat,imode) / sqrt(rfacq)
               end do
            end if
         end do

        ! Sum over modes and weight with areas and velocities v_k * v_k' for contribution of k',k pair to a2F.
        ! Mind, the areas were in (tpi/alat)**2 units.

        dsk_vk_2 = kpts_tr_area(ik) * kpts_tr_area(ikp) * (tpi/alat)**4 / (vabsk_tr(ik) * vabsk_tr(ikp) )
        dsk_vk_2 = dsk_vk_2 * factor_area_ibz(ish) ! adds weight of this irreducible wedge (IBZ) to the integral
        do iomega = 1, nomega  !frequencies
           omega = omega_ini + omega_step*real(iomega-1,dp)
           do imode = 1,3*nat
              rfacq = smeared_delta(omega-w_qint(imode), osmear_q)  !smear omega(q)
              alpha2F(is,js,imode,iomega) = alpha2F(is,js,imode,iomega) + &
                      (abs(gep_int(imode)))**2 * rfacq * dsk_vk_2
           end do !imode
        end do !iomega

      end do !spin js
      end do !spin is

    end do !k' (kpoint)
    end do !k' (sheet)

end do !k  (kpoint)
end do !k  (sheet)

  ! divide by N(E_F) and 1BZ volume
  ! TODO:
  ! DUDA con la normalizacion... para variar...
  ! DUDA con el spin... 
  alpha2F=alpha2F/dosef/vol1bz

  !  print out by spin blocks
  do is=1,nspin 
  do js=1,nspin 

     write(unit_a2f,*) '#spins = ', is,js

     ! print a2F(omega)
     omega = omega_ini
     do iomega=1, nomega  !frequencies
        write(unit_a2f,'(e16.4,2x,15e16.4)') omega,  sum(alpha2F(is,js,:,iomega)), &
                alpha2F(is,js,:,iomega)
        omega = omega + omega_step
     end do

     ! calculate lambda 
     write(unit_a2f,*) '# lambda = 2 \int dw a2F(w)/w for each mode'
     do imode = 1,3*nat
        lambda(imode) = 0.0_dp
        omega = omega_ini + omega_step
        do iomega=2,nomega !skip omega=0 value
           lambda(imode) = lambda(imode) + alpha2F(is,js,imode,iomega) / omega
           omega = omega + omega_step
        end do
        lambda(imode) = lambda(imode) * 2.0_dp * omega_step
        write(unit_a2f,*) '#imode = ', imode, ', lambda = ', lambda(imode)
     end do
     write(unit_a2f,*) '#total lambda = ', sum(lambda)

     write(unit_a2f,*) ' '
     write(unit_a2f,*) ' '

   end do    !js
   end do    !is

  close(unit_a2f)

!================================================================================
!   Finish
!================================================================================

    deallocate(aep_mat_el)

    stop     
end program a2F_on_trFS
