program interpolate
        ! Uses nnkp and ham_r files to write an interpolated band structure.
        ! Additionally a .kpt with the path for band plotting is needed.
        ! In principle, only those 4 files are needed by this utility,
        ! with additional details for DOS given here explicitely.

    use kinds, only: dp
    use intw_utility, only: find_free_unit
    use intw_w90_setup, only: interpolate_1k, interpolated_DOS, &
            allocate_and_read_ham_r, allocate_and_read_u_mesh
    use intw_input_parameters, only: nk1, nk2, nk3, mesh_dir, prefix, read_input
    use intw_reading, only: read_parameters_data_file_xml, set_num_bands, &
                            num_bands_intw, num_wann_intw

!================================================================================
!       Declare the variables
!================================================================================
    implicit none

    character(256) :: nnkp_file, kpath_file, bpath_file, dos_file
    logical        :: read_status
    integer :: io_unit_k, io_unit_b, io_unit_d, nkpath, ne, nki1, nki2, nki3, ik, ie, iw
    real(dp) :: kpoint(3), kpoint_dist
    real(dp) :: eini, efin, esmear, estep
    real(dp) , allocatable :: eig_int(:), DOS(:), PDOS(:,:), weights(:)
    complex(dp) , allocatable :: u_int(:,:)


!--------------------------------------------------------------------------------
!================================================================================
!       Talk to the user
!================================================================================
write(*,'(A)') '====================================================='
write(*,'(A)') '|              program interpolate                  |'
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
!
! From this, I only need seednames, actually, 
! until intw.in is modified to include postprocessing options
!================================================================================
    call read_parameters_data_file_xml()  
   
!================================================================================
! Set the number of bands for the calculation
!================================================================================

    ! (this will read num_wann_intw and num_bands_intw dimensions from .nnkp if available)
    call set_num_bands()

    print *, num_wann_intw

!================================================================================
!       read u_mesh, ham_r file from w902intw
!================================================================================

    call allocate_and_read_ham_r()

!================================================================================
!   Test: interpolate bands over path
!================================================================================

    allocate (eig_int(num_wann_intw))
    allocate (weights(num_wann_intw))
    allocate (u_int(num_wann_intw,num_wann_intw))

    io_unit_k = find_free_unit()
    kpath_file = trim(mesh_dir)//trim(prefix)//trim('_band.kpt')
    open(unit=io_unit_k,file=kpath_file,status='old')

    io_unit_b = find_free_unit()
    bpath_file = trim(mesh_dir)//trim(prefix)//trim('.bnd_int')
    open(unit=io_unit_b,file=bpath_file,status='unknown')

    read(io_unit_k,*) nkpath
    do ik = 1,nkpath
       read(io_unit_k,*) kpoint
       call interpolate_1k (kpoint, eig_int, u_int)
       do iw = 1,num_wann_intw
          weights = (abs(u_int(iw,:)))**2   
          write(io_unit_b,'(2i4,20e14.6)') ik, iw, eig_int(iw),weights
       end do
    end do

    close(io_unit_k)
    close(io_unit_b)

    write(*,*)' Bands finished'

!================================================================================
!   Test: DOS over a fine grid
!================================================================================

    ne = 100
    eini = -10.0_dp
    efin = 10.0_dp
    esmear = 0.05_dp
    estep = (efin-eini)/real(ne-1,dp)
    nki1 = 20
    nki2 = 20
    nki3 = 20
    allocate (DOS(ne), PDOS(ne,num_wann_intw))
    call interpolated_DOS (nki1, nki2, nki3, eini, efin, esmear, ne, DOS, PDOS)
    io_unit_d = find_free_unit()
    dos_file = trim(mesh_dir)//trim(prefix)//trim('.dos_int')
    open(unit=io_unit_d,file=dos_file,status='unknown')
    do ie = 1,ne
       write(io_unit_d,'(20e14.6)') eini+(ie-1)*estep, DOS(ie), PDOS(ie,1:num_wann_intw)
    end do
    close(io_unit_d)

    write(*,*)' DOS finished'

!==========================
! Finally make sure allocate_and_read_u_mesh works (I do not need this for interpolation)
!==========================
!!!     call allocate_and_read_u_mesh()


!================================================================================
!       clean up and finish
!================================================================================

    deallocate (eig_int,DOS)

    write(*,'(A)') '====================================================='
    write(*,'(A)') '|               end program interpolate             |'
    write(*,'(A)') '|        ---------------------------------          |'
    write(*,'(A)') '====================================================='


end program interpolate
