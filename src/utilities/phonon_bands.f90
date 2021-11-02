!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       -----------------------------------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program ep_melements 

!  use intw_pseudo

  use intw_useful_constants
  use intw_intw2wannier
  use intw_tests
  use intw_symmetries
  use intw_fft
  use intw_input_parameters

  use intw_band_crossing

  use intw_ph 
 
!================================================================================
!       Declare the variables 
!================================================================================
implicit none

integer       ::       ikpt
integer       ::       i_folder
integer       ::       i, j, k, ik 
integer       ::       switch 
integer       ::       nk_irr , nkmesh

logical       ::       read_status, have_nnkp

logical       ::       k_points_consistent


character(256)::       nnkp_file, method



real(dp),allocatable:: kpoints_QE(:,:), kmesh(:,:)

real(dp)            :: time1, time2

real(dp),allocatable :: kpoints_irr  (:,:) 
complex(dp), allocatable ::dyn_q(:,:)
real (dp), allocatable :: w2(:)
real(dp) ::q_point(3)
integer :: ii,jj
character(len=256) :: fc_file_name
real(dp) :: at_frc(3,3)

integer :: number_of_special_points, nk_vec_path
real(dp),allocatable     :: list_k_path(:,:)
real(dp),allocatable     :: list_x_path(:)
real(dp),allocatable     :: list_xticks(:)
integer :: io_unit

!--------------------------------------------------------------------------------
call get_timing(time1)
!================================================================================
!       Talk to the user
!================================================================================
write(*,20) '====================================================='
write(*,20) '|                  program me                       |'
write(*,20) '|        ---------------------------------          |'
write(*,20) '====================================================='
write(*,20) '|    waiting for input file...                      |'


!================================================================================
!       read the input file 
!       Read in the necessary information from standard input
!================================================================================
call read_input(read_status)

method = trim(intw2W_method)

if (read_status ) then
       stop 
end if        

!================================================================================
!       read the parameters from the SCF QE calculation, in  
!       particular, read in the symmetry matrices!!! 
!================================================================================

call read_parameters_data_file_xml()

!================================================================================
!       Check that $prefix.nnkp is present
!================================================================================
nnkp_file = trim(mesh_dir)//trim(prefix)//".nnkp"

inquire(file=nnkp_file,exist=have_nnkp)

if(.not. have_nnkp) then
   write(*,20)      '**********************************************************'
   write(*,20)      '* Could not find the file '//trim(nnkp_file)
   write(*,20)      '* Did you run W90 -pp $seed to get the parameter file?   ' 
   write(*,20)      '**********************************************************'
stop
end if
write(*,20) '|       - .nnkp file found                          |'
write(*,20) '|           ---------------------------------       |'

!================================================================================
!       read in the kpoints from the QE folders
!================================================================================
! read the koints information
allocate(kpoints_QE(3,nkpoints_QE))
call read_kpoints_data_file_xml(kpoints_QE) 

!================================================================================
! Build the kmesh corresponding to the parameters in the input file 
!================================================================================
nkmesh = nk1*nk2*nk3
allocate(kmesh(3,nkmesh))
call generate_kmesh(kmesh,nk1,nk2,nk3)

!================================================================================
! check that kmesh and nnkp_kpoints are consistent with one another 
!       This insures that the Wannier data is consistent with intw data.
!================================================================================
call intw2W90_check_mesh(nkmesh,kmesh)

write(*,20) '|       - The mesh in the Wannier90 input.win file  |'
write(*,20) '|         and the intw mesh are equal.              |'
write(*,20) '|           ---------------------------------       |'

if     (nspin==1) then
   write(*,20) '|       - The calculation is paramagnetic nspin=1   |'
   write(*,20) '|                                                   |'
   write(*,20) '|           ---------------------------------       |'
elseif (nspin==2) then
   write(*,20) '|       - Spin calculation nspin = 2                |'
   if (noncolin) then
      write(*,20) '|         Non-collinear Spin calculation            |'
   endif
   write(*,20) '|           ---------------------------------       |'
else
   write(*,20) '*****************************************************'
   write(*,20) '* ERROR: Allowed values for nspin are 1 or 2        *'
   write(*,20) '*            program stops.                         *'
   write(*,20) '*****************************************************'
   stop
endif



!================================================================================
!    Find the size of the irreducible set of k-points (IBZ)
!    and check that the number of kpoints corresponds to either
!    a full mesh or the IBZ.
!================================================================================
call find_size_of_irreducible_k_set(nk1,nk2,nk3,kmesh,nk_irr)

!This is only for testing: The result for nk_irr is different in both.
allocate(kpoints_irr(3,nk1*nk2*nk3))
call find_the_irreducible_k_set(nk1,nk2,nk3,kmesh,kpoints_irr,nk_irr)
deallocate(kpoints_irr )

if (nkpoints_QE /= nkmesh .and. nkpoints_QE /= nk_irr) then
        ! the points in the folder are not consistent with the
        ! input file. Stop the program!
        write(*,20) '*****************************************************'
        write(*,20) '*      The number of kpoints present in the QE      *'
        write(*,20) '*      folders are not consistent with a full       *'
        write(*,20) '*      Brillouin Zone or an irreducible Brillouin   *'
        write(*,20) '*      zone! Review your input...                   *'
        write(*,20) '*                   Program stops.                  *'
        write(*,20) '*****************************************************'
        write(*,20) '* debug information:                                *'
        write(*,*) '*        nkpoints_QE = ',nkpoints_QE 
        write(*,*) '*        nkmesh      = ',nkmesh
        write(*,*) '*        nk_irr      = ',nk_irr, nsym, tr_symmetry

        stop
else if (nkpoints_QE == nk_irr) then
        ! The points in the QE folders *could* be a valid choice for the IBZ;
        ! this must be checked!
        full_mesh = .false.
        IBZ       = .true.
else if (nkpoints_QE == nkmesh ) then
        ! The points in the QE folders *could* be consistent with a full mesh;
        ! this must be checked!
        full_mesh = .true.
        IBZ       = .false.
end if


!================================================================================
!      allocate the symmetry arrays 
!      CAREFUL! the subroutine needs to know the global value of "full_mesh",
!      so it is crucial that this allocation occurs AFTER setting full_mesh
!================================================================================
call allocate_symmetry_related_k(nk1,nk2,nk3,nsym)

!================================================================================
!     Compute the indices of the inverse rotation matrices.
!     This will be useful later in the exectution.
!================================================================================
call find_inverse_symmetry_matrices_indices()

!================================================================================
!     Set up the array spin_symmetry_matrices, which contain
!     the matrices which must be used to rotate spin
!================================================================================
call allocate_and_build_spin_symmetry_matrices(nsym)



!================================================================================
!      Fill the symmetry arrays appropriately
!================================================================================

call set_symmetry_relations(nk1, nk2, nk3, nkpoints_QE, kpoints_QE,kmesh, k_points_consistent, &
             QE_folder_sym, sym_G, symlink )



!================================================================================
!       Tell the user what is in the QE folders
!================================================================================
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
        write(*,20)      '**********************************************************'
        write(*,20)      '* The kpoints present in the QE folders are not consistent' 
        write(*,20)      '* with the parameters of the input file!                 ' 
        write(*,20)      '**********************************************************'
        write(*,20) '* debug information:                                *'
        write(*,*) '*        nkpoints_QE = ',nkpoints_QE 
        write(*,*) '*        nkmesh      = ',nkmesh
        write(*,*) '*        nk_irr      = ',nk_irr, nsym, tr_symmetry
        write(*,*) '*        IBZ         = ',IBZ
        stop
end if


!================================================================================
!       Check that the requested calculation is possible
!================================================================================
if (intw2W_fullzone ) then          
   write(*,20) '|       - intw2W_fullzone = .true.                  |'
   if ( full_mesh ) then
        write(*,20) '|         all k-points are explicitely calculated   |'
        write(*,20) '|         no symmetry is assumed.                   |'
        write(*,20) '|             (This is mostly for testing)          |'
        write(*,20) '|           ---------------------------------       |'
   else 
        write(*,20)      '**********************************************************'
        write(*,20)      '* A full mesh is not present in the QE folders!          '
        write(*,20)      '* The requested calculation is impossible.               ' 
        write(*,20)      '*                   program stops.                       ' 
        write(*,20)      '**********************************************************'
        stop
   end if
else 
        write(*,20) '|       - intw2W_fullzone = .false.                 |'
        write(*,20) '|         Symmetries will be utilized.              |'
        write(*,20) '|           ---------------------------------       |'
end if

        write(*,20) '|       - reading pseudopotentials from UPF files     |'
        write(*,20) '|            (defined in .save data files)          |'
        write(*,20) '|                                                   |'
!call readpp
        write(*,20) '|                       PPs OK                     |'
        write(*,20) '|           ---------------------------------       |'

call read_ph_information_xml()

fc_file_name=trim(trim(mesh_dir)//"/"//trim(ph_dir)//"/"//trim(fc_mat))

call readfc ( fc_file_name &
            , nq1, nq2, nq3, nat, alat, at_frc, ntyp, amass )

!================================================================================
!       Get the special points from the file
!================================================================================
if (bands_points /= 0 ) then
        call get_interpolation_special_points( number_of_special_points)
end if

!================================================================================
!       Define the interpolation path     
!================================================================================

nk_vec_path = bands_points

allocate(list_k_path(3,nk_vec_path))
allocate(list_x_path(nk_vec_path))
allocate(list_xticks(number_of_special_points))

call build_list_kpath(number_of_special_points,nk_vec_path,    &
                      list_k_path,list_x_path, list_xticks)

allocate(dyn_q(3*nat,3*nat), w2(3*nat))

io_unit=find_free_unit()
open(unit=io_unit, file=trim(ph_bands_file),status="unknown")

do i=1,nk_vec_path
 q_point(1:3) =  list_k_path(1:3,i)
 call mat_inv_four_t ( q_point, nq1, nq2, nq3,  3*nat , frc , dyn_q )
 call diagonalize_cmat (3*nat, dyn_q,w2)
 write(io_unit, "(i5,10000f18.12)")  i, sqrt(abs(w2))
enddo

close(unit=io_unit)

deallocate(dyn_q, w2)


!================================================================================
!       Finish 
!================================================================================
call get_timing(time2)
write(*,20) '|                     ALL DONE                      |'
write(*,30) '|     total time: ',time2-time1,' seconds            |'
write(*,20) '====================================================='


20 format(A)
30 format(A,F8.2,6X,A)

contains

subroutine diagonalize_cmat (n,a,w)
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

end program ep_melements

