! MBR February 2024

! April: TODO implement a kpath / qpath generator

program interpolatephonons
        ! Uses nnkp and ham_r files to write an interpolated band structure.
        ! Additionally a .qpt with the path for band plotting is needed.
        ! In principle, only those 4 files are needed by this utility,
        ! with additional details for DOS given here explicitely.

    use kinds, only: dp
    use intw_useful_constants, only: cmplx_0
    use intw_utility, only: find_free_unit, find_k_1BZ_and_G, switch_indices, &
            generate_kmesh, cryst_to_cart
    use intw_input_parameters, only: nk1, nk2, nk3, mesh_dir, prefix, read_input, &
            nq1, nq2, nq3, ph_dir, fc_mat
    use intw_reading, only: read_parameters_data_file_xml, set_num_bands, &
                            nat, nspin, bg, at, alat
    use intw_ph, only: nqmesh, qmesh, read_ph_information_xml
    use intw_ph_interpolate


!================================================================================
!       Declare the variables
!================================================================================
    implicit none

    character(256) :: fc_file_name, qband_file_name, phband_file_name
    logical        :: read_status
    integer :: iq, nq, ir
    integer :: nqmesh_int, nq1_int, nq2_int, nq3_int
    integer :: ph_unit, qpt_unit
    real(dp) :: qpoint(3), qpoint_cart(3), rcart(3)
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


!================================================================================
!       read the parameters from the SCF QE calculation
!================================================================================
    call read_parameters_data_file_xml()  

!================================================================================
!   Phonon mesh and interpolation meshes
!================================================================================

    nqmesh=nq1*nq2*nq3
    allocate(qmesh(3,nqmesh))
    call generate_kmesh(qmesh, nq1,nq2,nq3)

    ! Wigner-Seitz
    call allocate_and_build_ws_irvec_q()

    ! interpolation grid (alternatively, we use special directions to plot phonon dispersion)
    nq1_int = 10*nq1
    nq2_int = 10*nq2
    nq3_int = 10*nq3
    nqmesh_int = nq1_int*nq2_int*nq3_int
    allocate(qmesh_int(3,nqmesh_int))
    call generate_kmesh(qmesh_int,nq1_int,nq2_int,nq3_int)

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
  ! 1. read force constants
!  fc_file_name = trim(mesh_dir)//trim(ph_dir)//trim(fc_mat)
!  call allocate_and_build_dyn_qmesh(fc_file_name)
  ! 2. read dyn files
  call allocate_and_build_dyn_qmesh2()

  ! transform to R space
  call dyn_q_to_dyn_r()

  do iq=1,nqmesh
     print *, iq, qmesh(:,iq), w2_q(:,iq)
  end do

  ! test decay of dyn_r elements with distance
  do ir=1,nrpts_q
     rcart = real(irvec_q(:,ir),dp)
     call cryst_to_cart (1, rcart, at, 1)
     rcart = rcart * alat  ! bohr units
     write(520,'(i5,f16.6,8e16.4)') ir,  sqrt ( sum(rcart*rcart) ), &
             abs(dyn_r(1,1,ir)), abs(dyn_r(1,2,ir)), abs(dyn_r(1,4,ir)), abs(dyn_r(1,5,ir))
  end do

  !================================================================================
  ! Interpolate on fine grid
  !================================================================================

  allocate( dyn_qint(3*nat,3*nat), u_qint(3*nat,3*nat), w2_qint(3*nat), w_qint(3*nat) )

!  do iq=1,nqmesh !nqmesh_int
!     qpoint = qmesh(:,iq) !qmesh_int(:,iq)
!     qpoint_cart=qpoint
!     call cryst_to_cart (1, qpoint_cart, bg, +1)
!     call dyn_interp_1q(qpoint, dyn_qint)
!     call dyn_diagonalize_1q(3*nat, dyn_qint, u_qint, w2_qint)
!     print *, iq, qpoint, qpoint_cart, w2_qint
!  end do

  ! Alternatively, read qpoints from file cu_band.qpt (same as cu_band.kpt from wannphon test)
  qband_file_name = trim(mesh_dir)//trim(prefix)//trim("_band.qpt")
  phband_file_name = trim(mesh_dir)//trim(prefix)//trim("_qband.dat")
  qpt_unit = find_free_unit()
  ph_unit = qpt_unit + 1
  print *, qpt_unit, ph_unit
  open(qpt_unit,file=qband_file_name,status='old')
  open(ph_unit,file=phband_file_name,status='unknown')
  read(qpt_unit,*) nq
  do iq=1,nq
     read(qpt_unit,*) qpoint
     call dyn_interp_1q(qpoint, dyn_qint)
     call dyn_diagonalize_1q(3*nat, dyn_qint, u_qint, w2_qint)
     w_qint=sign(sqrt(abs(w2_qint)),w2_qint)
     !freqs are given in meV units.
     !pass qpoint to cartesian for output.
     qpoint_cart=qpoint
     call cryst_to_cart (1, qpoint_cart, bg, +1)
     write(ph_unit,'(i3,3f10.4, 15e16.6)') iq, qpoint_cart, w_qint
  end do
  close(qpt_unit)
  close(ph_unit)


!================================================================================
!       clean up and finish
!================================================================================

    write(*,'(A)') '====================================================='
    write(*,'(A)') '|               end program phonons                 |'
    write(*,'(A)') '|        ---------------------------------          |'
    write(*,'(A)') '====================================================='



end program interpolatephonons
