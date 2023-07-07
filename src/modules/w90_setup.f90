!----------------------------------------------------------------------------!
!       intw project.
!
!----------------------------------------------------------------------------!
!
module intw_w90_setup
  !
  !----------------------------------------------------------------------------!
  !
  ! TODO in new version, what is the module that contains num_bands, num_wann?
  ! I think it should be intw_intw2wannier.
  !
  ! The utility W902INTW should:
  ! 1.- call read_nnkp_file in intw_intw2wannier to obtain nnkp_exclude_bands, nnkp_kpoints, etc (see use line)
  ! 2.- call this read_w90_chk
  ! 3.- call allocate_and_build_u_mesh
  ! 4.- write_formatted_u_mesh
  !
  ! Then, any utility using u_mesh can simply use the information obtained by
  ! read_formatted_u_mesh
  !
  ! TODO think of a wannier type. It may be filled up when doing read_formatted_u_mesh.
  !
  !----------------------------------------------------------------------------!
  !
  use kinds, only: dp
  implicit none
  !
  ! variables
  public :: use_disentanglement
  public :: n_wss, n_ws_search     !!! search space for WS vectors
  public :: irvec, nrpts, ndegen   !!! these will substitute w90_hamiltonian: irvec, nrpts, ndegen
  public :: ham_k, ham_r           !!! ham_k, ham_r
  public :: ndimwin, lwindow   !!! DUDA ... probably we do not need these outside the module 
                                 !!! if the information is contained in u_mesh...
  public :: u_mesh
  !
  ! subroutines
  public :: read_w90_chk, allocate_and_build_u_mesh, write_formatted_u_mesh, allocate_and_build_ws_irvec, &
            allocate_and_build_ham_k, allocate_and_build_ham_r
  !
  private
  !
  save
  !
  logical :: use_disentanglement
  logical, allocatable :: lwindow(:,:)
  integer :: nrpts
  integer, allocatable :: ndegen(:), ndimwin(:)
  integer, parameter :: n_wss=27  !! TODO give somewhere as input
  integer, dimension(3) , parameter :: n_ws_search =(/ 1,1,1 /) !! TODO give somewhere as input
  real(kind=dp), allocatable :: irvec(:,:)
  complex(kind=dp), allocatable :: u_mesh(:,:,:), ham_k(:,:,:), ham_r(:,:,:)
  complex(kind=dp), allocatable  :: u_matrix(:,:,:), u_matrix_opt(:,:,:)
  !
contains
  !
  !----------------------------------------------------------------------------!
  subroutine read_w90_chk()
  !----------------------------------------------------------------------------!
  use intw_reading, only: nbands
  use intw_input_parameters, only: mesh_dir, prefix, nk1, nk2, nk3
  use intw_utility, only: find_free_unit
  use intw_intw2wannier, only: nnkp_exclude_bands, nnkp_excluded_bands, nnkp_num_kpoints, &
           nnkp_kpoints, num_bands_intw, num_wann_intw  
  implicit none

  character(20) :: checkpoint
  character(33) :: header
  character(256) :: filename
  integer :: io_unit_chk, nb, nexc, nkpt, nw, nntot
  integer :: n1, n2, n3, i, j, ik, iw, ib
  integer :: exc_bands(nnkp_exclude_bands)
  real(kind=dp) :: dir_latt(3,3), rec_latt(3,3)
  real(kind=dp), allocatable  :: kpts(:,:)

  io_unit_chk = find_free_unit()
  filename = trim(mesh_dir)//trim(prefix)//trim('.chk')
  open(unit=io_unit_chk,file=filename,status='unknown')

  ! TODO add comments in consistency check with existing data from nnkp
  ! (I think that full-zone k-list should be read from there, since 
  ! in the intw2wannier step it was checked that both lists coincide)

  ! .nnkp file should be read before calling this subroutine

  ! bands
  read (io_unit_chk) header
  read (io_unit_chk) nb  ! no. of bands
  !JLB
  !if ( nb .ne. nbands ) then
  if ( nb .ne. num_bands_intw ) then
   write(*,*) "Number of bands in .chk is not the same as nbands-exclude_bands in .nnkp file"
   write(*,*) "Stopping..."
   stop
  end if
  read (io_unit_chk) nexc  ! no. of excluded bands
  if ( nnkp_exclude_bands .ne. nexc) stop
  read (io_unit_chk) exc_bands(1:nexc)  ! excluded bands
  !! TODO check if ( nnkp_excluded_bands .ne. exc_bands ) stop
  
  ! lattice
  read (io_unit_chk) ((dir_latt(i,j), i=1,3), j=1,3)  ! direct
  read (io_unit_chk) ((rec_latt(i,j), i=1,3), j=1,3)  ! reciprocal

  ! k grid
  read (io_unit_chk) nkpt
  if ( nkpt .ne. nnkp_num_kpoints) stop
  read (io_unit_chk) n1, n2, n3
  if ( (n1 .ne. nk1) .or. (n2 .ne. nk2) .or. (n3 .ne. nk3) ) stop
  allocate (kpts(3,nkpt))
  read (io_unit_chk) ((kpts(i,ik), i=1,3), ik=1,nkpt)
  ! TODO check if ( kpts == nnkp_kpoints )

  read (io_unit_chk) nntot ! no. of neighbours, not used
  read (io_unit_chk) nw
  if ( nw .ne. num_wann_intw) then
   write(*,*) "Number of Wannier functions in .chk is not the same as number of projections in .nnkp file"
   write(*,*) "Stopping..."
   stop
  end if
  read (io_unit_chk) checkpoint ! not used
  read (io_unit_chk) use_disentanglement

  ! u_matrix_opt
  
  allocate (lwindow(num_bands_intw,nnkp_num_kpoints), ndimwin(nnkp_num_kpoints))
  allocate (u_matrix_opt(num_bands_intw,num_wann_intw,nnkp_num_kpoints))

  if (use_disentanglement) then
       read (io_unit_chk) ((lwindow(ib,ik),ib=1,num_bands_intw),ik=1,nnkp_num_kpoints)
       read (io_unit_chk) (ndimwin(ik), ik=1,nnkp_num_kpoints)
       read (io_unit_chk) (((u_matrix_opt(ib,iw,ik),ib=1,num_bands_intw),iw=1,num_wann_intw),ik=1,nnkp_num_kpoints)
  end if     

  ! u_matrix
  allocate (u_matrix(num_wann_intw,num_wann_intw,nnkp_num_kpoints))
  read (io_unit_chk) (((u_matrix(ib,iw,ik),ib=1,num_wann_intw),iw=1,num_wann_intw),ik=1,nnkp_num_kpoints)

  ! These are also contained in chk, but we will not use them in INTW:
  !
  ! m_matrix
  ! centers
  ! spreads

  close(io_unit_chk)
  return
  end subroutine read_w90_chk
  !
  !----------------------------------------------------------------------------!
    subroutine allocate_and_build_u_mesh()
    !----------------------------------------------------------------------------!
    ! This is pretty much produce_u_mesh from W90_to_wannier
    !
    ! This subroutine computes the "rotation" matrices which must be applied
    ! to the matrix elements to take them from band space to Wannier space.
    !
    ! If disentanglement is used, u_mesh    = u_matrix_opt x u_matrix
    ! If no disentanglement is used, u_mesh = u_matrix
    !
    ! note that the number of k-points in the coarse mesh, in Wannier90,
    ! is defined in the variable num_kpts. Consistency demands that this number
    ! be equal to nkmesh, the same quantity in the intw program.
    !
    ! TODO explore the possibility of symmetrization...
    !----------------------------------------------------------------------------!
    use intw_useful_constants, only: cmplx_0
    use intw_intw2wannier, only: nnkp_num_kpoints, num_bands_intw, num_wann_intw, nnkp_excluded_bands  
    implicit none

    integer :: i, ik, n1, n2, nb1

    call read_w90_chk ()       ! allocate and read use_disentanglement, lwindow, u_matrix_opt, u_matrix

    allocate (u_mesh(num_bands_intw,num_wann_intw,nnkp_num_kpoints))  !fullzone k-points 
    u_mesh = cmplx_0

    if (use_disentanglement) then     ! DUDA no seria mejor añadir aqui lwindow??? (reordenar m_opt, se puede hacer aqui o
                                      !  cada vez que se use u_mesh con la lista de eigvals...)
      do ik=1, nnkp_num_kpoints
        n1 = 0
        do nb1=1, num_bands_intw
        !!!! if (excluded_bands(nb1)) cycle !! DUDA excluded_bands or lwindow  ?????
          n1 = n1 + 1
          do n2=1,num_wann_intw
            do i=1,num_wann_intw
              u_mesh(n1,n2,ik) = u_mesh(n1,n2,ik) + u_matrix_opt(n1,i,ik)*u_matrix(i,n2,ik)
            enddo !i
          enddo !n2
        enddo !nb1
      enddo !ik
    else
      u_mesh = u_matrix
    endif
    !
    ! deallocate these arrays from W90, which are no longer useful
    if (use_disentanglement) deallocate(u_matrix_opt)
    deallocate(u_matrix)
    !
    return
  end subroutine allocate_and_build_u_mesh


  !
  !----------------------------------------------------------------------------!
  !subroutine read_formatted_u_mesh ()
  !        TODO first take a decision on how will it be printed.
  !        Create a type wannier and store stuff there upon reading
  !----------------------------------------------------------------------------!


  !----------------------------------------------------------------------------!
  subroutine write_formatted_u_mesh (eigenval)
  !----------------------------------------------------------------------------!
  use intw_reading, only: nbands
  use intw_input_parameters, only: mesh_dir, prefix
  use intw_utility, only: find_free_unit
  use intw_intw2wannier, only: nnkp_exclude_bands, nnkp_excluded_bands, &
          nnkp_num_kpoints, num_bands_intw, num_wann_intw, nnkp_kpoints 
  implicit none

  character(256) :: filename
  integer :: i,ik,ib,iw,io_unit_u
  real(kind=dp), intent(in) :: eigenval(num_bands_intw,nnkp_num_kpoints)

  io_unit_u = find_free_unit()
  filename = trim(mesh_dir)//trim(prefix)//trim('.u_mesh')
  open(unit=io_unit_u,file=filename,status='unknown')

  write(io_unit_u,*)'NBANDS'
  write(io_unit_u,*) nbands
  write(io_unit_u,*)'EXCLUDED_BANDS'
  write(io_unit_u,*) nnkp_exclude_bands
  write(io_unit_u,*) (nnkp_excluded_bands(i),i=1,nbands)
  write(io_unit_u,*)'num_bands  num_wann  num_kpt'
  write(io_unit_u,*) num_bands_intw, num_wann_intw, nnkp_num_kpoints
  write(io_unit_u,*)'USE_DISENTANGLEMENT'
  write(io_unit_u,*) use_disentanglement

  do ik = 1,nnkp_num_kpoints
     write(io_unit_u,*) 'K-POINT', ik
     write(io_unit_u,"(3f16.10)") (nnkp_kpoints(i,ik), i=1,3 )
     write(io_unit_u,*) 'LWINDOW'
     write(io_unit_u,*) (lwindow(ib,ik), ib=1,num_bands_intw)
     write(io_unit_u,*) 'EIGENVALUES'
     write(io_unit_u,"(5es18.8,/)") (eigenval(ib,ik), ib=1,num_bands_intw)
     write(io_unit_u,*) 'u_mesh'
     do ib = 1,num_bands_intw
        write(io_unit_u,"(5es18.8,/)")  (u_mesh(ib,iw,ik), iw=1,num_wann_intw)
     end do
  end do

  close(io_unit_u)
  end subroutine write_formatted_u_mesh
  !
  !----------------------------------------------------------------------------!
  subroutine allocate_and_build_ws_irvec (nk1,nk2,nk3)
  !----------------------------------------------------------------------------!
  !  Calculate real-space Wigner-Seitz lattice vectors
  !----------------------------------------------------------------------------!
  !
  use intw_reading, only: alat, at
  use intw_useful_constants, only: eps_6
  use intw_utility, only: generate_kmesh
  !
  implicit none
  !
  integer, intent(in) :: nk1,nk2,nk3
  !
  integer :: ik, nws, i,j,k,l
  integer :: permu(n_wss), rdeg_ws_max(nk1*nk2*nk3) 
  real(kind=dp) :: kmesh(3,nk1*nk2*nk3)
  real(kind=dp) :: r_cryst(3,n_wss), r_length(n_wss), r_cart(3), r_ws_max(3,n_wss,nk1*nk2*nk3)
  !
  call generate_kmesh (kmesh,nk1,nk2,nk3)
  !
  nws = 0  ! total number of WS vectors
  do ik = 1,nk1*nk2*nk3
      l = 0
      do i = -n_ws_search(1),n_ws_search(1)
      do j = -n_ws_search(2),n_ws_search(2)
      do k = -n_ws_search(3),n_ws_search(3)
         l = l + 1
         ! generate equivalents to point kmesh(:,ik)
         r_cryst(1,l) =  kmesh(1,ik) + real(i,dp)
         r_cryst(2,l) =  kmesh(2,ik) + real(j,dp)
         r_cryst(3,l) =  kmesh(3,ik) + real(k,dp)
         r_cart = r_cryst(:,l)
         !R-vector from crystallographic to cartesian
         call cryst_to_cart (1, r_cryst, at, 1)
         r_cart = r_cart * alat  ! bohr units
         r_length(l) = sqrt ( sum(r_cart*r_cart) )
     end do   
     end do   
     end do   
     ! order by ascending length
     call HPSORT_real(n_wss,r_length,permu)
     ! store first vector (shortest)
     r_ws_max(:,1,ik) = r_cryst(:,permu(1))
     rdeg_ws_max(ik) = 1
     nws = nws + 1
     ! detect degeneracies and store vectors if degenerate
     do l = 2,n_wss
        if ( abs(r_length(l) - r_length(l-1))<eps_6) then
           r_ws_max(:,l,ik) = r_cryst(:,permu(l))
           rdeg_ws_max(ik) = rdeg_ws_max(ik) + 1
           nws = nws + 1
        else
           exit
        end if
     end do  
  end do
  !
  ! Data for Wannier
  !
  ! build WS kpoint list and degeneracies
  nrpts = nws
  allocate ( irvec(3,nrpts) )
  allocate ( ndegen(nrpts) )
  ! allocate and build WS kpoint list and degeneracies (irvec, ndegen)
  i = 0
  do ik = 1,nk1*nk2*nk3
     do l = 1, rdeg_ws_max(ik)
         i = i + 1
         irvec(:,i) = r_ws_max(:,l,ik) 
         ndegen(i) = rdeg_ws_max(ik)
     end do
  end do

  return
  end subroutine allocate_and_build_ws_irvec

  !
  !----------------------------------------------------------------------------!
  subroutine allocate_and_build_ham_k (eigenval)
          ! eigenval has already removed excl_bands
  !----------------------------------------------------------------------------!
  !  Construct ham_k with u_matrices on the fullzone k-mesh
  !----------------------------------------------------------------------------!
  !
  use intw_useful_constants, only: cmplx_0
  use intw_intw2wannier, only: nnkp_num_kpoints, num_bands_intw, num_wann_intw, nnkp_kpoints  
  !
  implicit none
  !
  real(kind=dp) , intent(in) :: eigenval(num_bands_intw,nnkp_num_kpoints)
  !
  integer :: ik, iw, jw, nwin
  complex(kind=dp) :: cterm
  ! 
  allocate (ham_k(num_wann_intw,num_wann_intw,nnkp_num_kpoints))
  ham_k = cmplx_0
  do ik = 1,nnkp_num_kpoints
     nwin = ndimwin(ik)
     do iw = 1,num_wann_intw
       do jw = 1,iw 
           if (use_disentanglement) then
                   !  DUDA ver mi comentario sobre lwindow en allocate_and_build_u_mesh...
                   !  excl_bands should have been already excluded.
              ! Pick up eigenvalues inside the opt window for this k.
              ! They correspond to the first items in u_matrix_opt.
              cterm = sum( conjg(u_mesh(1:nwin,iw,ik)) * pack(eigenval(:,ik),lwindow(:,ik)) * u_mesh(1:nwin,jw,ik) )
           else
              cterm = sum( conjg(u_mesh(:,iw,ik)) * eigenval(:,ik) * u_mesh(:,jw,ik) )
           end if
           ham_k(iw,jw,ik) = ham_k(iw,jw,ik) + cterm
           ! force hermiticity
           ham_k(jw,iw,ik) = ham_k(jw,iw,ik) + cterm
       end do
     end do
  end do
  !
  return
  end subroutine allocate_and_build_ham_k

  !!TODO we may need another routine allocate_and_build_ham_k_for1k_only. 
  ! Will need to find the k in the list first
  
  !
  !----------------------------------------------------------------------------!
  subroutine allocate_and_build_ham_r ()
  !----------------------------------------------------------------------------!
  !  Construct ham_r by Fourier transform of ham_k, k in the fullzone k-mesh
  !----------------------------------------------------------------------------!
  !
  use intw_reading, only: alat, bg
  use intw_useful_constants, only: tpi, cmplx_0, cmplx_i
  use intw_intw2wannier, only: nnkp_num_kpoints, num_bands_intw, num_wann_intw, nnkp_kpoints  
  !
  implicit none
  !
  integer :: ik, i
  complex(kind=dp) :: phasefac
  !
  allocate (ham_r(num_wann_intw, num_wann_intw,nrpts))
  ham_r = cmplx_0
  do i = 1,nrpts
     do ik = 1,nnkp_num_kpoints
          phasefac = exp( -cmplx_i*tpi*sum( nnkp_kpoints(:,ik)*real(irvec(:,i) ) ) )
          ham_r(:,:,i) = ham_r(:,:,i) + phasefac*ham_k(:,:,ik)
     enddo
  enddo
  ham_r = ham_r/real(nnkp_num_kpoints,dp) 

  return
  end subroutine allocate_and_build_ham_r
!----------------------------------------------------------------------------!
end module intw_w90_setup
!----------------------------------------------------------------------------!


!!! Finally, TODO modules for interpolation utilities
