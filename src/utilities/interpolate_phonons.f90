! MBR June 2024


program interpolatephononstau
        ! Like interpolate_phonons.f90, but using atom-pair-adapter WS vectors

    use kinds, only: dp
    use intw_useful_constants, only: cmplx_0, Ha_to_Ry, Ha_to_eV
    use intw_utility, only: find_free_unit, find_k_1BZ_and_G,  &
            generate_kmesh, cryst_to_cart, smeared_delta, &
            generate_and_allocate_kpath
    use intw_input_parameters, only: read_input, read_cards, &
            nk1, nk2, nk3, mesh_dir, prefix, &
            nq1, nq2, nq3, ph_dir, fc_mat, &
            read_cards, exist_qpath, nqpath, nqspecial, qspecial, exist_kpath, &
            nq1_dosph, nq2_dosph, nq3_dosph, &
            nomega, omega_ini, omega_fin, osmear_q, read_for_dynmat
    use intw_reading, only: read_parameters_data_file_xml, set_num_bands, &
                            nat, nspin, bg, at, alat, tpiba
    use intw_ph, only: nqmesh, qmesh, read_ph_information_xml
    use intw_ph_interpolate


!================================================================================
!       Declare the variables
!================================================================================
    implicit none

    character(256) :: fc_file_name, phband_file_name
    logical        :: read_status
    integer :: iq, nq, ir, iq1, iq2, iq3, iomega, imode
    integer :: ph_unit, dos_unit
    real(dp) :: omega_step, omega, rfacq
    real(dp) :: qpoint(3), qpoint_cart(3), rcart(3)
    real(dp) , allocatable :: qpath(:,:), dqpath(:), dosph(:,:)
    real(dp) , allocatable :: qmesh_int(:,:), w2_qint(:), w_qint(:)
    complex(dp) , allocatable :: dyn_qint(:,:), u_qint(:,:)

!--------------------------------------------------------------------------------
!================================================================================
!       Talk to the user
!================================================================================
write(*,'(A)') '====================================================='
write(*,'(A)') '|         program phonons interpolate               |'
write(*,'(A)') '|        ---------------------------------          |'
write(*,'(A)') '====================================================='
write(*,'(A)') '|    waiting for input file...                      |'

!================================================================================
!       read the input file
!       Read in the necessary information from standard input
!================================================================================
    call read_input(read_status)
 
    if (read_status ) then
       stop
    end if

    ! generate q-list for phonon bands plot with special points in Q_PATH
    call read_cards
    if ( .not. exist_qpath) then
       write(*,*)' Q_PATH not found. Phonon bands/DOS cannot be interpolated. Stopping.'
       stop
    end if        

!================================================================================
!       read the parameters from the SCF QE calculation
!================================================================================
    call read_parameters_data_file_xml()  

!================================================================================
!   Build qpoint path to plot bands.
!   The nqpath number of points from the input might fluctuate.
!================================================================================

    call generate_and_allocate_kpath (at, bg, tpiba, nqpath, nqspecial, qspecial, qpath, dqpath)
    !write(*,*) nqpath
    !do iq=1,nqpath
    !   print *,  qpath(:,iq), dqpath
    !end do

!================================================================================
!   Phonon meshes
!================================================================================

    nqmesh=nq1*nq2*nq3
    allocate(qmesh(3,nqmesh))
    call generate_kmesh(qmesh, nq1,nq2,nq3)

    ! Generate Wigner-Seitz mesh with atom position correction
    call allocate_and_build_ws_irvec_qtau()

    write(*,*) 'Meshes generated'

  !================================================================================
  !       Read all the information about phonons and pseudopotentials
  !================================================================================
  call read_ph_information_xml()

  !================================================================================
  ! Read the force constant matrix from this file.
  ! Calculate dynamical matrices in qmesh and transform to real space
  !================================================================================
  !
  ! two options to get dyn_q:
  !
  print *,' get dyn_q matrices'
  if (read_for_dynmat == 'fc' ) then ! read force constants
          call allocate_and_build_dyn_qmesh(fc_mat)
  else if (read_for_dynmat == 'dynq' ) then ! read dyn files        
          call allocate_and_build_dyn_qmesh2()
  end if        

  ! do iq=1,nqmesh
  !   ! print"(i,3f,x,100f)", iq, qmesh(:,iq), sign(sqrt(abs(w2_q(:,iq))), w2_q(:,iq)) ! a.u.
  !   print"(i,3f,x,100f)", iq, qmesh(:,iq), sign(sqrt(abs(w2_q(:,iq))), w2_q(:,iq))*Ha_to_eV*1000.0_dp ! meV
  !   ! print"(i,3f,x,100f)", iq, qmesh(:,iq), sign(sqrt(abs(w2_q(:,iq))), w2_q(:,iq))*Ha_to_eV*8065.610_dp ! cm-1
  !   ! print"(i,3f,x,100f)", iq, qmesh(:,iq), sign(sqrt(abs(w2_q(:,iq))), w2_q(:,iq))*6579.6839205_dp ! THz
  ! end do

  !
  ! transform to R space with atom position correction
  print *,' transform to R space'
  call dyn_q_to_dyn_rtau()
  
  ! test decay of dyn_r elements with distance
!  do ir=1,nrpts_q
!     rcart = real(irvec_q(:,ir),dp)
!     call cryst_to_cart (1, rcart, at, 1)
!     rcart = rcart * alat  ! bohr units
!     write(520,'(i5,f16.6,8e16.4)') ir,  sqrt ( sum(rcart*rcart) ), &
!             abs(dyn_r(1,1,ir)), abs(dyn_r(1,2,ir)), abs(dyn_r(1,4,ir)), abs(dyn_r(1,5,ir))
!  end do

  !================================================================================
  ! Interpolate on qpath for bands plot
  !================================================================================

  print *,' interpolate phonons'

  allocate( dyn_qint(3*nat,3*nat), u_qint(3*nat,3*nat), w2_qint(3*nat), w_qint(3*nat) )

  phband_file_name = trim(mesh_dir)//trim(prefix)//trim(".qbnd_int")
  ph_unit = find_free_unit()
  open(ph_unit,file=phband_file_name,status='unknown')
  write(ph_unit,'(A)') '# q-point   omega(imode=1)[meV]  omega(2)[meV]   omega(3)[meV] ...'
  do iq=1,nqpath
     qpoint = qpath(:,iq)
     call dyn_interp_1q_tau(qpoint, dyn_qint)
     call dyn_diagonalize_1q(3*nat, dyn_qint, u_qint, w2_qint) ! freqs are given in a.u
     w_qint = sign(sqrt(abs(w2_qint)),w2_qint) * Ha_to_eV*1000.0_dp
     write(ph_unit,'(20e14.6)') dqpath(iq), w_qint
  end do
  close(ph_unit)
  write(*,*)' Phonon bands interpolation finished and written to file', phband_file_name

  !================================================================================
  ! Interpolate on fine q-grid for bands plot
  ! Parameters of DOS plot from namelist /DOS_ph/
  !================================================================================

  omega_step = (omega_fin-omega_ini)/real(nomega-1,dp)
  allocate (dosph(3*nat,nomega))
  dosph = 0.0_dp

  ! fine q-grid
  do iq1 = 1, nq1_dosph
     qpoint(1) = real(iq1-1,dp) / real(nq1_dosph,dp)
  do iq2 = 1, nq2_dosph
     qpoint(2) = real(iq2-1,dp) / real(nq2_dosph,dp)
  do iq3 = 1, nq3_dosph
     qpoint(3) = real(iq3-1,dp) / real(nq3_dosph,dp)

        ! Interpolate frequency in qpoint
        call dyn_interp_1q_tau(qpoint, dyn_qint)
        call dyn_diagonalize_1q(3*nat, dyn_qint, u_qint, w2_qint)
        w_qint=sign(sqrt(abs(w2_qint)),w2_qint) 
        !phonon frequency in a.u.: pass to Ry
        w_qint= w_qint*Ha_to_Ry
        ! Smear omega(q) for DOS (gaussian)
        do imode = 1,3*nat
           do iomega=1, nomega  !frequencies in Ry
              omega = omega_ini + omega_step*real(iomega-1,dp)
              rfacq = smeared_delta(omega-w_qint(imode), osmear_q)  
              dosph(imode,iomega) = dosph(imode,iomega) + rfacq
           end do   
       end do    

   end do
   end do
   end do
   dosph = dosph / real(nq1_dosph*nq2_dosph*nq3_dosph,dp) ! normalize for Nq points

   ! write DOS to file
   phband_file_name = trim(mesh_dir)//trim(prefix)//trim(".qdos_int")
   ph_unit = find_free_unit()
   open(ph_unit,file=phband_file_name,status='unknown')
   write(ph_unit,'(A)') '# omega[Ry]   PDOS(imode=1)  PDOS(imode=2)  PDOS(imode=3) ...'
   do iomega=1,nomega
       omega = omega_ini + omega_step*real(iomega-1,dp)
       write(ph_unit,'(20e14.6)') omega, dosph(:,iomega)
   end do
   close(ph_unit)

   write(*,'(A)') '|  DOS sum test:                                    |' 
   write(*,*)'       DOS integral (trapeze) = ', omega_step*sum(dosph(:,:)) 
   write(*,*)'       Number of modes =', 3*nat
   write(*,*)' Phonon DOS interpolation finished and written to file', phband_file_name


!================================================================================
!       clean up and finish
!================================================================================

    write(*,'(A)') '====================================================='
    write(*,'(A)') '|               end program phonons                 |'
    write(*,'(A)') '|        ---------------------------------          |'
    write(*,'(A)') '====================================================='



end program interpolatephononstau
