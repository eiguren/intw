!
! Copyright (C) 2024 INTW group
!
! This file is part of INTW.
!
! INTW is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! INTW is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <https://www.gnu.org/licenses/>.
!
module intw_w90_setup
  !
  !----------------------------------------------------------------------------!
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
  !----------------------------------------------------------------------------!
  !
  use kinds, only: dp

  implicit none
  !
  ! variables
  public :: use_disentanglement
  public :: n_wss, n_ws_search     ! search space for WS vectors
  public :: irvec, nrpts, ndegen   ! these will substitute w90_hamiltonian: irvec, nrpts, ndegen
  public :: ham_k, ham_r           ! ham_k, ham_r
  public :: ndimwin, lwindow       ! NOTE probably we do not need these outside the module
  public :: u_mesh                 ! unitary matrices
  public :: eigenval_intw          ! coarse-mesh eigenvalues used in wannier (read from .eig)

  ! subroutines
  public :: read_w90_chk, allocate_and_build_u_mesh, write_formatted_u_mesh, allocate_and_build_ws_irvec, &
            allocate_and_build_ham_k, allocate_and_build_ham_r, read_eig, write_ham_r, &
            allocate_and_read_ham_r, allocate_and_read_u_mesh, &
            interpolate_1k, interpolated_DOS, &
            wann_rotate_matrix, wann_fourier_1index, wann_IFT_1index, wann_FT_1index_1k
  !
  private
  !
  save
  !
  logical :: use_disentanglement
  logical, allocatable :: lwindow(:,:)
  integer :: nrpts
  integer, allocatable :: ndegen(:), ndimwin(:)
  ! NOTE n_wss and n_ws_search  could be given somewhere as input
  integer, parameter :: n_wss=27
  integer, dimension(3) , parameter :: n_ws_search =(/ 1,1,1 /)
  !real(kind=dp), allocatable :: irvec(:,:)
  integer, allocatable :: irvec(:,:)
  real(kind=dp), allocatable :: eigenval_intw(:,:)
  complex(kind=dp), allocatable :: u_mesh(:,:,:), ham_k(:,:,:), ham_r(:,:,:)
  complex(kind=dp), allocatable  :: u_matrix(:,:,:), u_matrix_opt(:,:,:)
  !
contains
  !
  !----------------------------------------------------------------------------!
  subroutine read_w90_chk()
  !----------------------------------------------------------------------------!
  ! This subroutine assumes set_num_bands() and read_nnkp_file()
  ! have been called before
  !----------------------------------------------------------------------------!
  use intw_reading, only: nbands, num_bands_intw, num_wann_intw
  use intw_input_parameters, only: mesh_dir, prefix, nk1, nk2, nk3
  use intw_utility, only: find_free_unit
  use intw_useful_constants, only: eps_8
  use intw_intw2wannier, only: nnkp_exclude_bands, nnkp_num_kpoints, nnkp_kpoints
  implicit none

  character(20) :: checkpoint
  character(33) :: header
  character(256) :: filename
  integer :: io_unit_chk, nb, nexc, nkpt, nw, nntot
  integer :: n1, n2, n3, i, j, ik, iw, ib
  integer :: exc_bands(nnkp_exclude_bands)
  real(kind=dp) :: dir_latt(3,3), rec_latt(3,3), omega_invariant, kp(3)
  real(kind=dp), allocatable  :: kpts(:,:)

  io_unit_chk = find_free_unit()
  filename = trim(mesh_dir)//trim(prefix)//trim('.chk')
  open(unit=io_unit_chk,file=filename,status='old',form='unformatted')

  ! .nnkp file should be read before calling this subroutine

  ! bands
  read (io_unit_chk) header
  read (io_unit_chk) nb  ! no. of bands
  !JLB
  !if ( nb .ne. nbands ) then
  if ( nb .ne. num_bands_intw ) then
     write(*,*) "Number of bands in .chk is not the same as num_bands_intw = nbands-exclude_bands in .nnkp file"
     write(*,*) "Stopping..."
     stop
  end if
  read (io_unit_chk) nexc  ! no. of excluded bands
  if ( nnkp_exclude_bands .ne. nexc) then
     write(*,*) " nexc in .chk is not the same as nnkp_exclude_bands in .nnkp file"
     write(*,*) "Stopping..."
     stop
  end if
  read (io_unit_chk) exc_bands(1:nexc)  ! excluded bands

  ! lattice
  read (io_unit_chk) ((dir_latt(i,j), i=1,3), j=1,3)  ! direct
  read (io_unit_chk) ((rec_latt(i,j), i=1,3), j=1,3)  ! reciprocal

  ! k grid
  read (io_unit_chk) nkpt
  if ( nkpt .ne. nnkp_num_kpoints) then
     write(*,*) "Number of k-points in .chk is not the same as in .nnkp file"
     write(*,*) "Stopping..."
     stop
  end if
  read (io_unit_chk) n1, n2, n3
  if ( (n1 .ne. nk1) .or. (n2 .ne. nk2) .or. (n3 .ne. nk3) ) then
     write(*,*) " nk1, nk2, nk3 in .chk is not the same as in INTW"
     write(*,*) "Stopping..."
     stop
  end if
  allocate (kpts(3,nkpt))
  read (io_unit_chk) ((kpts(i,ik), i=1,3), ik=1,nkpt)
  ! Check if kpts == nnkp_kpoints.
  ! Note that in intw2W90 we already checked that nnkp_kpoints
  ! coincide with the full-zone points provided by generate_kmesh
  do ik = 1,nkpt
     kp=kpts(:,ik)-nnkp_kpoints(:,ik)
     if ( sqrt(sum(kp*kp)) > eps_8 ) then
        write(*,*) "kpts list in .chk different from nnkp_kpoints (and INTW mesh)"
        write(*,*) "Stopping..."
        stop
     end if
  end do

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
       read (io_unit_chk) omega_invariant     ! not used
       read (io_unit_chk) ((lwindow(ib,ik),ib=1,num_bands_intw),ik=1,nnkp_num_kpoints)
       read (io_unit_chk) (ndimwin(ik), ik=1,nnkp_num_kpoints)
       read (io_unit_chk) (((u_matrix_opt(ib,iw,ik),ib=1,num_bands_intw),iw=1,num_wann_intw),ik=1,nnkp_num_kpoints)
  end if

  ! u_matrix
  allocate (u_matrix(num_wann_intw,num_wann_intw,nnkp_num_kpoints))
  read (io_unit_chk) (((u_matrix(ib,iw,ik),ib=1,num_wann_intw),iw=1,num_wann_intw),ik=1,nnkp_num_kpoints)

  ! These are also contained in .chk but, as for now, we do not use them in INTW:
  !
  ! m_matrix
  ! centers
  ! spreads

  close(io_unit_chk)
  return
  end subroutine read_w90_chk
  !
  !----------------------------------------------------------------------------!
  subroutine read_eig(eigenval)
  !----------------------------------------------------------------------------!
  ! Allocates and reads the eigenvalues for the relevant num_bands_intw from the .eig file    !
  !----------------------------------------------------------------------------!
  use intw_input_parameters, only: mesh_dir, prefix
  use intw_utility, only: find_free_unit
  use intw_reading, only: num_bands_intw
  use intw_intw2wannier, only: nnkp_num_kpoints

  ! I/O
  real(kind=dp), optional, intent(out) :: eigenval(num_bands_intw,nnkp_num_kpoints)
  ! Local variables
  character(256) :: filename
  integer :: io_unit, ik, ib, i, j

  allocate (eigenval_intw(num_bands_intw,nnkp_num_kpoints))

  io_unit=find_free_unit()
  filename = trim(mesh_dir)//trim(prefix)//trim('.eig')
  open(unit=io_unit,file=filename,form='formatted',status='old')
  !
  do ik=1,nnkp_num_kpoints
    do ib=1,num_bands_intw
      !
      read(io_unit,*) i, j, eigenval_intw(ib,ik)
      !
    enddo
  enddo

  close(io_unit)

  if (present(eigenval)) eigenval = eigenval_intw

  end subroutine
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
  !----------------------------------------------------------------------------!
  use intw_useful_constants, only: cmplx_0
  use intw_reading, only: num_bands_intw, num_wann_intw
  use intw_intw2wannier, only: nnkp_num_kpoints
  implicit none

  integer :: i, ik, n1, n2, nb1

  call read_w90_chk ()       ! allocate and read use_disentanglement, lwindow, u_matrix_opt, u_matrix

  allocate (u_mesh(num_bands_intw,num_wann_intw,nnkp_num_kpoints))  !fullzone k-points
  u_mesh = cmplx_0

  if (use_disentanglement) then

    do ik=1, nnkp_num_kpoints
       n1 = 0
       do nb1=1, num_bands_intw

          ! MBR-JBL: decided to use lwindow as in w90, instead of exclude_bands as in previus INTW version
          if (.not. lwindow(nb1, ik)) cycle

          ! At each k, bands are reordered so that
          ! u_matrix_opt for bands within the disentanglement window are filled first
          ! and the rest are just padded with zeros
          n1 = n1 + 1
          do n2=1,num_wann_intw
            do i=1,num_wann_intw
               !u_mesh(n1,n2,ik) = u_mesh(n1,n2,ik) + u_matrix_opt(n1,i,ik)*u_matrix(i,n2,ik)
               !MBR 10/05/24 correction. It should be:
               u_mesh(nb1,n2,ik) = u_mesh(nb1,n2,ik) + u_matrix_opt(n1,i,ik)*u_matrix(i,n2,ik)
            enddo !i
          enddo !n2
       enddo !nb1
    enddo !ik

  else
    u_mesh = u_matrix

  endif
  !
  !deallocate arrays which are no longer useful
  !if (use_disentanglement) deallocate(u_matrix_opt)
  !deallocate(u_matrix)
  !
  return
  end subroutine allocate_and_build_u_mesh



  !----------------------------------------------------------------------------!
  subroutine write_formatted_u_mesh()
  !----------------------------------------------------------------------------!
  use intw_reading, only: nbands, num_bands_intw, num_wann_intw
  use intw_input_parameters, only: mesh_dir, prefix
  use intw_utility, only: find_free_unit
  use intw_intw2wannier, only: nnkp_exclude_bands, nnkp_excluded_bands, &
                               nnkp_num_kpoints, nnkp_kpoints
  implicit none

  character(256) :: filename
  integer :: i,ik,ib,iw,io_unit_u

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
     !
     ! MBR 20/12/24 change format: remove blank lines
     write(io_unit_u,*) 'K-POINT', ik
     write(io_unit_u,"(3f16.10)") (nnkp_kpoints(i,ik), i=1,3 )
     write(io_unit_u,*) 'LWINDOW'
     write(io_unit_u,*) (lwindow(ib,ik), ib=1,num_bands_intw)
     write(io_unit_u,*) 'EIGENVALUES'
     write(io_unit_u,"(5es18.8)") (eigenval_intw(ib,ik), ib=1,num_bands_intw)
     write(io_unit_u,*) 'u_mesh'
     do ib = 1,num_bands_intw
        write(io_unit_u,"(5es18.8)")  (u_mesh(ib,iw,ik), iw=1,num_wann_intw)
     end do
     !
     if (use_disentanglement) then
       write(io_unit_u,*) 'u_matrix_opt'
       do ib = 1,num_bands_intw
          write(io_unit_u,"(5es18.8)")  (u_matrix_opt(ib,iw,ik), iw=1,num_wann_intw)
       end do
       write(io_unit_u,*) 'u_matrix'
       do ib = 1,num_wann_intw
          write(io_unit_u,"(5es18.8)")  (u_matrix(ib,iw,ik), iw=1,num_wann_intw)
       end do
     end if
     !
  end do

  close(io_unit_u)
  end subroutine write_formatted_u_mesh
  !
  !
  subroutine allocate_and_build_ws_irvec(nk1,nk2,nk3)
  !----------------------------------------------------------------------------!
  !  Calculate real-space Wigner-Seitz lattice vectors
  !----------------------------------------------------------------------------!
  !
  use intw_reading, only:  at, alat
  use intw_useful_constants, only: eps_8
  use intw_utility, only: generate_kmesh, cryst_to_cart
  !
  implicit none
  !
  integer, intent(in) :: nk1,nk2,nk3
  !
  logical :: in_ws
  integer :: ik, nboundary, i,j,k,l, l0,l1
  integer :: r_cryst_int(3), Rs(3,n_wss), ndegen_ws(nk1*nk2*nk3*n_wss), irvec_ws(3,nk1*nk2*nk3*n_wss)
  real(kind=dp) :: kmesh(3,nk1*nk2*nk3)
  real(kind=dp) :: r_cryst(3), r_length_l, r_length_l1, r_cart(3)
  !
  call generate_kmesh (kmesh,nk1,nk2,nk3)
  !
  ! generate superlattice replica vectors search mesh
  l = 0
  do i = -n_ws_search(1),n_ws_search(1)
  do j = -n_ws_search(2),n_ws_search(2)
  do k = -n_ws_search(3),n_ws_search(3)
         l = l + 1
         Rs(:,l) = (/ i,j,k /)
         if (i == 0 .and. j == 0 .and. k == 0) l0=l ! Origin O
  end do
  end do
  end do
  !
  nrpts = 0  ! total number of WS vectors
  do ik = 1,nk1*nk2*nk3
     !
     do l = 1, n_wss
        ! r-R(l), where for r-supercell-vector I use a conventional cell mesh of size nk1, nk2, nk3
        ! and R(l) runs over replicas
        r_cryst = ( kmesh(:,ik) - real(Rs(:,l),dp) ) * real( (/ nk1, nk2, nk3 /), dp)
        r_cryst_int = nint(r_cryst)
        !R-vector from crystallographic to cartesian
        r_cart = r_cryst
        call cryst_to_cart (1, r_cart, at, 1)
        r_length_l =  alat * sqrt ( sum(r_cart*r_cart) )  ! distance of r-R(l) to O (cartesian, bohr)
        !
        ! r-R(l) is in the WS if its distance to O is shorter than its
        ! distance to any other O' origin.
        ! If it is equidistant, it lies on the boundary and is degenerate.
        in_ws = .true.
        nboundary = 1
        !
        ! Loop over origins O' given by R(l1)
        do l1 = 1, n_wss
          ! r-R(l)-R(l1)
          r_cryst = ( kmesh(:,ik) - real(Rs(:,l)+Rs(:,l1),dp) ) * real( (/ nk1, nk2, nk3 /), dp)
          r_cart = r_cryst
          call cryst_to_cart (1, r_cart, at, 1)
          r_length_l1 = alat * sqrt ( sum(r_cart*r_cart) )  ! distance of r-R(l) to O' (cartesian, bohr)
          ! compare distances leaving a eps_8 gap
          if ( r_length_l > r_length_l1 + eps_8 .and. l1/=l0) then ! not in the WS => remove vector from list
                  in_ws = .false.
                  exit
          else if ( abs(r_length_l-r_length_l1)<=eps_8 .and. l1/=l0) then ! on the boundary => add degeneracy
                  nboundary = nboundary + 1
          end if
        end do
        !
        ! store r-R(l) and its degeneracy if it is inside WS
        if (in_ws) then
             nrpts=nrpts+1
             irvec_ws(:,nrpts) = r_cryst_int
             ndegen_ws(nrpts) = nboundary
        end if
        !
     end do
  end do ! ik
  !
  ! Data for Wannier: WS kpoint list and degeneracies.
  ! Simply dismiss the array at sites >nrpts, which have not been used
  allocate ( irvec(3,nrpts) )
  allocate ( ndegen(nrpts) )
  ndegen = ndegen_ws(1:nrpts)
  do i=1,3
       irvec(i,:) = irvec_ws(i,1:nrpts)
  end do

  !JLB test
  !write (*, '(I5)') nrpts
  !write (*, '(15I5)') (ndegen(i), i=1, nrpts)
  !do i = 1, nrpts
  !      write (*, '(i5,3I6)') i, irvec(:, i)
  !end do

  return
  end subroutine allocate_and_build_ws_irvec

  !
  !----------------------------------------------------------------------------!
  subroutine allocate_and_build_ham_k()
  !----------------------------------------------------------------------------!
  !  Construct ham_k with u_matrices on the fullzone k-mesh
  !----------------------------------------------------------------------------!
  !
  use intw_useful_constants, only: cmplx_0
  use intw_reading, only: num_wann_intw
  use intw_intw2wannier, only: nnkp_num_kpoints
  !
  implicit none
  !
  !
  integer :: ik, iw, jw, nwin
  complex(kind=dp) :: cterm
  !
  allocate (ham_k(num_wann_intw,num_wann_intw,nnkp_num_kpoints))
  ham_k = cmplx_0
  do ik = 1,nnkp_num_kpoints
     nwin = ndimwin(ik)
     do jw = 1,num_wann_intw
       do iw = 1,jw
           ! JLB: Needs to be checked. I think lwindow has been incorporated in u_mesh building so here just multiply.
           if (use_disentanglement) then
                   !  According to comment on lwindow in allocate_and_build_u_mesh,
                   !  excl_bands should have been already excluded.
              ! Pick up eigenvalues inside the opt window for this k.
              ! They correspond to the first items in u_matrix_opt.
              cterm = sum( conjg(u_mesh(1:nwin,iw,ik)) * pack(eigenval_intw(:,ik),lwindow(:,ik)) * u_mesh(1:nwin,jw,ik) )
           else
              cterm = sum( conjg(u_mesh(:,iw,ik)) * eigenval_intw(:,ik) * u_mesh(:,jw,ik) )
           end if
           !ham_k(iw,jw,ik) = ham_k(iw,jw,ik) + cterm
           !! force hermiticity
           !ham_k(jw,iw,ik) = ham_k(jw,iw,ik) + cterm
           !JLB (also changed order of jw, iw loops above):
           ham_k(iw,jw,ik) = cterm
           if(iw .lt. jw) ham_k(jw,iw,ik) = conjg(cterm)
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
  subroutine allocate_and_build_ham_r()
  !----------------------------------------------------------------------------!
  !  Construct ham_r by Fourier transform of ham_k, k in the fullzone k-mesh
  !----------------------------------------------------------------------------!
  !
  use intw_reading, only: num_wann_intw
  use intw_useful_constants, only: tpi, cmplx_0, cmplx_i
  use intw_intw2wannier, only: nnkp_num_kpoints, nnkp_kpoints
  !
  implicit none
  !
  integer :: ik, i, ib, jb
  complex(kind=dp) :: phasefac
  !
  allocate (ham_r(num_wann_intw, num_wann_intw,nrpts))
  ham_r = cmplx_0
  do i = 1,nrpts
     do ik = 1,nnkp_num_kpoints
       do ib = 1, num_wann_intw
         do jb = 1, num_wann_intw
           phasefac = exp( -cmplx_i*tpi*sum( nnkp_kpoints(:,ik)*irvec(:,i) ) )
           ham_r(ib,jb,i) = ham_r(ib,jb,i) + phasefac*ham_k(ib,jb,ik)
         end do
       end do
     enddo
  enddo
  ham_r = ham_r/real(nnkp_num_kpoints,dp)

  return
  end subroutine allocate_and_build_ham_r
  !
  !----------------------------------------------------------------------------!
  subroutine write_ham_r()
  !----------------------------------------------------------------------------!
  !  Write ham_r to file, to check decay (and compare to w90 if needed)
  !----------------------------------------------------------------------------!
  !
  use intw_utility, only: find_free_unit
  use intw_input_parameters, only: mesh_dir, prefix
  use intw_reading, only:  num_wann_intw
  use intw_intw2wannier, only : generate_header
  !
  implicit none
  !
  character(256) :: filename, header
  integer :: io_unit, ir, in, jn
  !
  io_unit = find_free_unit()
  filename = trim(mesh_dir)//trim(prefix)//trim('_hr_intw.dat')
  open (io_unit, file=filename, form='formatted', status='unknown')
  call generate_header(' ',header)
  write (io_unit, *) trim(header)
  write (io_unit, *) num_wann_intw
  write (io_unit, *) nrpts
  write (io_unit, '(15I5)') (ndegen(ir), ir=1, nrpts)
  do ir = 1, nrpts
    do in = 1, num_wann_intw
      do jn = 1, num_wann_intw
        ! JLB: I would increase the writing accuracy, here just same as w90 for comparison
        write (io_unit, '(5I5,2F12.6)') irvec(:, ir), jn, in, ham_r(jn, in, ir)
      end do
    end do
  end do

  close (io_unit)
  return
  end subroutine write_ham_r
  !

  ! MBR new routines 28/08/23:
  ! allocate_and_read_u_mesh and _ham_r: reads previously generated and
  ! printed u_mesh and _ham_r by utility W902intw.

  !
  !----------------------------------------------------------------------------!
  subroutine allocate_and_read_ham_r()
  !
  !----------------------------------------------------------------------------!
  !  MBR:
  !  This reads: nrpts, ndegen, irvec, ham_r
  !  (num_bands_intw, num_wann_intw have been read previously from nnkp using set_numbands)
  !  from formatted file mesh_dir/prefix_hr_intw.dat
  !----------------------------------------------------------------------------!
  !
  use intw_useful_constants, only: cmplx_0
  use intw_utility, only: find_free_unit
  use intw_reading, only: num_wann_intw
  use intw_input_parameters, only: mesh_dir, prefix
  !
  implicit none
  !
  character(256) :: filename, header
  integer :: io_unit, ir, in, jn, i, j
  !
  ! open file and read dimensions
  !
  io_unit = find_free_unit()
  filename = trim(mesh_dir)//trim(prefix)//trim('_hr_intw.dat')
  open (io_unit, file=filename, form='formatted', status='old')
  read (io_unit, *) header
  read (io_unit, *) ir
  if (ir .ne. num_wann_intw) then
      write(*,*) 'num_wann_intw in ', filename,' does not coincide with nnkp value. Stopping.'
      stop
  end if
  read (io_unit, *) nrpts
  !
  ! allocate arrays
  ! careful, in case they had been allocated by previous call to build ham_r
  ! note: if irvec has been calculated before, reading the _hr file will
  ! overwrite it
  !
  if (.not. allocated(ndegen)) &
     allocate (ndegen(nrpts))
  if (.not. allocated(irvec)) &
     allocate (irvec(3,nrpts) )
  if (.not. allocated(ham_r)) &
     allocate (ham_r(num_wann_intw, num_wann_intw,nrpts))
  ndegen = 0
  irvec = 0
  ham_r = cmplx_0
  !
  ! read irvec, ham_r
  !
  read (io_unit, '(15I5)') (ndegen(ir), ir=1, nrpts)
  do ir = 1, nrpts
    do in = 1, num_wann_intw
      do jn = 1, num_wann_intw
        read (io_unit, '(5I5,2F12.6)') irvec(:, ir), j, i, ham_r(jn, in, ir)
      end do
    end do
  end do
  !
  close (io_unit)
  !
  return
  end subroutine allocate_and_read_ham_r


  subroutine allocate_and_read_u_mesh
  !
  !----------------------------------------------------------------------------!
  !  MBR:
  !  This reads all dimensions in the W90 problem
  !  (except for nbands, num_bands_intw, num_wann_intw have been read previously
  !  from nnkp using set_numbands)
  !  which include: nnkp_num_kpoints, nnkp_exclude_bands,
  !  and quantities use_disentanglement, nnkp_excluded_bands, &
  !  nnkp_kpoints, eigenval, lwindow, u_mesh, etc.
  !  from formatted file mesh_dir/prefix_u_mesh.dat
  !
  !----------------------------------------------------------------------------!
  !
  use intw_reading, only: nbands, num_bands_intw, num_wann_intw
  use intw_input_parameters, only: mesh_dir, prefix
  use intw_utility, only: find_free_unit
  use intw_intw2wannier, only: nnkp_exclude_bands, nnkp_excluded_bands, &
                               nnkp_num_kpoints, nnkp_kpoints
  implicit none

  character(256) :: filename, varname
  integer :: i,ik,ib,iw,io_unit_u

  io_unit_u = find_free_unit()
  filename = trim(mesh_dir)//trim(prefix)//trim('.u_mesh')
  open(unit=io_unit_u,file=filename,status='old')

  read(io_unit_u,*) varname !'NBANDS'
  read(io_unit_u,*) ib
  if (ib .ne. nbands) then
     write(*,*) 'nbands in ', filename,' does not coincide with value in QE. Stopping.'
     stop
  end if
  read(io_unit_u,*) varname !'EXCLUDED_BANDS'
  read(io_unit_u,*) nnkp_exclude_bands
  if (.not. allocated(nnkp_excluded_bands)) &
     allocate (nnkp_excluded_bands(nbands))
  read(io_unit_u,*) (nnkp_excluded_bands(i),i=1,nbands)
  read(io_unit_u,*) varname !'num_bands  num_wann  num_kpt'
  read(io_unit_u,*) ib, iw, nnkp_num_kpoints
  if ( ib .ne. num_bands_intw ) then
     write(*,*) 'num_bands_intw in ', filename,' does not coincide with nnkp value. Stopping.'
     stop
  end if
  if ( iw .ne. num_wann_intw ) then
     write(*,*) 'num_wann_intw in ', filename,' does not coincide with nnkp value. Stopping.'
     stop
  end if
  read(io_unit_u,*) varname ! 'USE_DISENTANGLEMENT'
  read(io_unit_u,*) use_disentanglement
  !
  ! allocate quantities
  ! (careful in case allocate_and_build had been called first)
  !
  if (.not. allocated(nnkp_kpoints)) &
     allocate (nnkp_kpoints(3,nnkp_num_kpoints))
  if (.not. allocated(lwindow)) &
     allocate (lwindow(num_bands_intw,nnkp_num_kpoints))
  if (.not. allocated(eigenval_intw)) &
     allocate (eigenval_intw(num_bands_intw,nnkp_num_kpoints))
  if (.not. allocated(u_mesh)) &
     allocate (u_mesh(num_bands_intw,num_wann_intw,nnkp_num_kpoints))
  if (use_disentanglement) then
     if (.not. allocated(u_matrix_opt)) &
        allocate (u_matrix_opt(num_bands_intw,num_wann_intw,nnkp_num_kpoints))
     if (.not. allocated(u_matrix)) &
        allocate (u_matrix(num_wann_intw,num_wann_intw,nnkp_num_kpoints))
  end if
  !
  do ik = 1,nnkp_num_kpoints
     !
     ! MBR 20/12/24 change format: remove blank lines
     read (io_unit_u,*) varname ! 'K-POINT', ik
     read(io_unit_u,"(3f16.10)") (nnkp_kpoints(i,ik), i=1,3 )
     read(io_unit_u,*) varname ! 'LWINDOW'
     read(io_unit_u,*) (lwindow(ib,ik), ib=1,num_bands_intw)
     read(io_unit_u,*) varname ! 'EIGENVALUES'
     read(io_unit_u,"(5es18.8)") (eigenval_intw(ib,ik), ib=1,num_bands_intw)
     read(io_unit_u,*) varname ! 'u_mesh'
     do ib = 1,num_bands_intw
        read(io_unit_u,"(5es18.8)")  (u_mesh(ib,iw,ik), iw=1,num_wann_intw)
     end do
     !
     if (use_disentanglement) then
       read(io_unit_u,*) varname ! 'u_matrix_opt'
       do ib = 1,num_bands_intw
          read(io_unit_u,"(5es18.8)")  (u_matrix_opt(ib,iw,ik), iw=1,num_wann_intw)
       end do
       read(io_unit_u,*) varname ! 'u_matrix'
       do ib = 1,num_wann_intw
          read(io_unit_u,"(5es18.8)")  (u_matrix(ib,iw,ik), iw=1,num_wann_intw)
       end do
     end if
     !
  end do
  !
  close(io_unit_u)
  !
  end subroutine allocate_and_read_u_mesh



!!! Finally, modules for interpolation utilities

  subroutine interpolate_1k (kpoint, eig_int, u_interp)
  !
  !----------------------------------------------------------------------------!
  !  MBR:
  !  Given a kpoint in crystal coordinates, it interpolates the eigenvalues
  !  using already available variables eigenval_intw, ham_r, ir_vec, ndegen, nrpts.
  !
  !  Interpolation matrix, which can be used to obtain weights for fatbands,
  !  is an optional output
  !
  !----------------------------------------------------------------------------!
  !
  use intw_useful_constants, only: cmplx_0, cmplx_i, tpi
  use intw_reading, only: num_wann_intw
  !
  implicit none
  !
  integer :: ir,i,j
  real(kind=dp), intent(in) :: kpoint(3)
  real(kind=dp), intent(out) :: eig_int(num_wann_intw)
  complex(kind=dp) :: cfac
  complex(kind=dp) :: ham_pack_1k((num_wann_intw*(num_wann_intw+1))/2)
  complex(kind=dp) :: u_pass(num_wann_intw,num_wann_intw)
  complex(kind=dp) , optional, intent(out) :: u_interp(num_wann_intw,num_wann_intw)
  ! for ZHPEVX we use the same dimensions as W90 does
  integer :: neig_found, info
  integer :: iwork(5*num_wann_intw), ifail(num_wann_intw)
  real(kind=dp) :: rwork(7*num_wann_intw)
  complex(kind=dp) :: cwork(2*num_wann_intw)
  !
  external :: ZHPEVX
  !
  ! generate ham_k directly packed at this point
  !
  ham_pack_1k = cmplx_0
  do ir = 1,nrpts
     cfac = exp(cmplx_i*tpi*dot_product(kpoint(:),irvec(:,ir)))/real(ndegen(ir),dp)
     do j=1,num_wann_intw
       do i=1,j
          ham_pack_1k(i+((j-1)*j)/2) = ham_pack_1k(i+((j-1)*j)/2) + cfac * ham_r(i,j,ir)
       end do
     end do
  end do
  !
  ! diagonalize
  ! Note: actual Wannier interpolation matrix will be upass^dagger
  !
  eig_int=0.0_dp
  call ZHPEVX('V', 'A', 'U', num_wann_intw, ham_pack_1k, &
              0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
              neig_found, eig_int, u_pass, num_wann_intw, &
              cwork, rwork, iwork, ifail, info)
  !
  if (info < 0) then
     write(*,*) 'Wrong argument in ZHPEVX. Stopping.'
     stop
  else if (info > 0) then
     write(*,*) 'ZHPEVX failed. Stopping.'
     stop
  end if
  !
  ! Interpolation matrix each eigenergy on the WFs
  !
  if (present(u_interp)) then
      u_interp = conjg(transpose(u_pass))
  end if
  !
  end subroutine interpolate_1k



  subroutine interpolated_DOS (nik1, nik2, nik3, eini, efin, esmear, ne, DOS, PDOS)
  !
  !----------------------------------------------------------------------------!
  !  MBR:
  !
  !  Calculate DOS using a fine grid nik1 x nik2 x nik3
  !  and a lorentzian smearing.
  !  Optionally, write PDOS using weights from u_interp
  !
  !----------------------------------------------------------------------------!
  !
  use intw_useful_constants, only: cmplx_0, cmplx_i, tpi
  use intw_reading, only: num_wann_intw
  !
  implicit none
  !
  integer :: ik1, ik2, ik3, ie, iw
  integer, intent(in) :: nik1, nik2, nik3, ne
  real(kind=dp), intent(in) :: eini, efin, esmear
  real(kind=dp), intent(out) :: DOS(ne)
  real(kind=dp), optional, intent(out) :: PDOS(ne,num_wann_intw)
  real(kind=dp) :: ener, estep, lorentz
  real(kind=dp) :: kpoint(3), eig_int(num_wann_intw)
  complex(kind=dp) :: u_interp(num_wann_intw,num_wann_intw)

  estep=(efin-eini)/real(ne-1,dp)
  DOS = 0.0_dp
  if (present(PDOS)) PDOS = 0.0_dp

  ! construct fine grid of kpoints, interpolate 1 by 1 and add
  ! contribution to DOS(e)

  do ik1 = 1,nik1
     kpoint(1) = real(ik1-1,dp) / real(nik1,dp)
  do ik2 = 1,nik2
     kpoint(2) = real(ik2-1,dp) / real(nik2,dp)
  do ik3 = 1,nik3
     kpoint(3) = real(ik3-1,dp) / real(nik3,dp)
     !
     call interpolate_1k (kpoint, eig_int, u_interp)
     do ie = 1,ne
        ener = eini + (ie-1)*estep
        do iw = 1,num_wann_intw
          lorentz = 1.0_dp / ((ener-eig_int(iw))**2+esmear**2)
          lorentz = lorentz * 0.5_dp*esmear/tpi
          DOS(ie) = DOS(ie) + lorentz
          if (present(PDOS)) PDOS(ie,:) = PDOS(ie,:) + lorentz*(abs(u_interp(iw,:)))**2
        end do
     end do
     !
  end do
  end do
  end do
  DOS = DOS / real(nik1 * nik2 * nik3, dp)  ! normalize for Nk points
  PDOS = PDOS / real(nik1 * nik2 * nik3, dp)  ! normalize for Nk points
  !
  end subroutine interpolated_DOS



  subroutine wann_rotate_matrix (ik1, ik2, matin, matout)
  !
  !----------------------------------------------------------------------------!
  !  MBR 9/1/24:
  !  Given a matrix 2x2 of elements on the coarse k-mesh, this routine
  !  rotates it to the Wannier gauge:
  !  Mat_out(ik1,ik2) =  U^dagger(ik1) * Mat_in(ik1,ik2) * U(ik2)
  !  where U is u_mesh, i.e. it already accounts for the entanglement option,
  !  taking the elements from band space to Wannier space.
  !
  !  It is assumed that the u_mesh has been built or read already.
  !----------------------------------------------------------------------------!
  !
  use intw_useful_constants, only: cmplx_0
  use intw_reading, only: num_bands_intw, num_wann_intw
  !
  implicit none
  !
  integer , intent(in) :: ik1,ik2
  complex(kind=dp) , intent(in) :: matin(num_bands_intw,num_bands_intw)
  complex(kind=dp) , intent(out) :: matout(num_wann_intw,num_wann_intw)
  !
  integer :: iw, jw, ib, jb
  !
  matout=cmplx_0
  do iw=1,num_wann_intw
  do jw=1,num_wann_intw
      do ib=1,num_bands_intw
      do jb=1,num_bands_intw
         matout(iw,jw) = matout(iw,jw) + conjg(u_mesh(ib,iw,ik1)) * matin(ib,jb) * u_mesh(jb,jw,ik2)
      end do
      end do
  end do
  end do
  !
  return
end subroutine wann_rotate_matrix


subroutine wann_fourier_1index (matL, matk, switch, sig)
  !
  !----------------------------------------------------------------------------!
  !  MBR 9/1/24:
  !  Given a matrix of elements to be interpolated, this routine acts on one of
  !  the indices, which contains the input dataset (complex).
  !  Depending on the switch, it does the (inverse) Fourier transform between the
  !  coarse k-mesh vectors and the direct lattice vectors L (Wigner-Seitz).
  !  The inverse FT on the matrix elements rotated to the Wannier gauge
  !  is what will be used in a parent routine or utility for interpolation.

  !  switch = 1 does L --> k ( FT)
  !  switch =-1 does k --> L (IFT)

  !  The optional sign "sig" choses the sign to be used in the exponential.
  !  If not specified, the usual (I)FT convention is used  (same as switch).
  !  It is useful when you want to (I)FT matrix elements on the bra/ket sides
  !  with the corresponding sign.
  !----------------------------------------------------------------------------!
  !
  use intw_useful_constants, only: cmplx_0, cmplx_i, tpi
  use intw_reading, only:  num_wann_intw
  use intw_intw2wannier, only: nnkp_num_kpoints, nnkp_kpoints
  !
  implicit none
  !
  integer , intent(in) :: switch
  integer , optional, intent(in) :: sig
  complex(kind=dp) , intent(inout) :: matk (num_wann_intw,num_wann_intw, nnkp_num_kpoints), &
        matL (num_wann_intw,num_wann_intw, nrpts)
  !
  integer :: ir, ik
  real :: signo
  complex(kind=dp) :: fac
  !
  if (switch .eq. -1) then  ! IFT
     !
     signo = -1.0_dp
     if (present(sig)) signo = real(sig,dp)
     !
     matL = cmplx_0
     do ir = 1,nrpts
        do ik = 1, nnkp_num_kpoints
            fac = exp(signo * cmplx_i*tpi*dot_product(nnkp_kpoints(:,ik),irvec(:,ir)))
            matL(:,:,ir) = matL(:,:,ir) + fac*matk(:,:,ik)
        end do
     end do
     matL = matL / real(nnkp_num_kpoints,dp)
     !
  else if (switch .eq. 1) then  ! FT
     !
     signo = 1.0_dp
     if (present(sig)) signo = real(sig,dp)
     !
     matk = cmplx_0
     do ik = 1, nnkp_num_kpoints
        do ir = 1,nrpts
            fac = exp(signo * cmplx_i*tpi*dot_product(nnkp_kpoints(:,ik),irvec(:,ir)))/real(ndegen(ir),dp)
            matk(:,:,ik) = matk(:,:,ik) + fac*matL(:,:,ir)
        end do
     end do
     !
  else
     write(*,*)' Error in switch of wann_fourier_1index. Stopping'
     stop
  end if
  !
  return
  end subroutine wann_fourier_1index

subroutine wann_IFT_1index (nks, kpoints, matk, matL)
  !
  !----------------------------------------------------------------------------!
  !  MBR 9/1/24:
  !  This does the same as wann_fourier_1index with switch -1,
  !  but for a given kpoints list
  !----------------------------------------------------------------------------!
  !
  use intw_useful_constants, only: cmplx_0, cmplx_i, tpi
  use intw_reading, only:  num_wann_intw
  !
  implicit none
  !
  integer , intent(in) :: nks
  real(kind=dp) , intent(in) :: kpoints(3,nks)
  complex(kind=dp) , intent(in) :: matk (num_wann_intw,num_wann_intw, nks)
  complex(kind=dp) , intent(out) :: matL (num_wann_intw,num_wann_intw, nrpts)
  !
  integer :: ir, ik
  complex(kind=dp) :: fac

     matL = cmplx_0
     do ir = 1,nrpts
        do ik = 1, nks
            fac = exp(-cmplx_i*tpi*dot_product(kpoints(:,ik),irvec(:,ir)))
            matL(:,:,ir) = matL(:,:,ir) + fac*matk(:,:,ik)
        end do
     end do
     matL = matL / real(nks,dp)

  return
end subroutine wann_IFT_1index


subroutine wann_FT_1index_1k (kpoint, matL, matk)
  !
  !----------------------------------------------------------------------------!
  !  MBR 9/1/24:
  !  As above, for one index of a matrix of elements in the direct lattice grid, and given any kpoint
  !  (this is the interpolation step) this calculates
  !   sum_L e^ikL mat(L) / degen(L)
  !----------------------------------------------------------------------------!
  !
  !
  use intw_useful_constants, only: cmplx_0, cmplx_i, tpi
  use intw_reading, only:  num_wann_intw
  !
  implicit none
  !
  real(kind=dp) , intent(in) :: kpoint(3)
  complex(kind=dp) , intent(in) :: matL(num_wann_intw,num_wann_intw, nrpts)
  complex(kind=dp) , intent(out) :: matk(num_wann_intw,num_wann_intw)
  !
  integer :: ir
  complex(kind=dp) :: fac
  !
  matk = cmplx_0
  do ir = 1,nrpts
      fac = exp(cmplx_i*tpi*dot_product(kpoint(:),irvec(:,ir)))/real(ndegen(ir),dp)
      matk(:,:) = matk(:,:) + fac*matL(:,:,ir)
  end do
  !
  return
end subroutine wann_FT_1index_1k


!----------------------------------------------------------------------------!
end module intw_w90_setup
!----------------------------------------------------------------------------!
