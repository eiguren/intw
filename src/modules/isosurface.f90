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
module triFS_isosurface

  !------------------------------------------------------------------!
  ! Module that contains all the necessary variables and subroutines !
  ! to read Wannier hamiltonian and create symmetric isosurface.     !
  !------------------------------------------------------------------!

  use kinds, only: dp

  implicit none

  ! Variables to be read from Tetrahedralization
  integer, public :: nnode, ntetra, nibz_nodes, nibz_faces ! Total number of tetrahedral nodes and faces (both small and big)
  real(dp), allocatable, public :: ncoord(:,:) ! Coordinates of small tetrahedra nodes
  integer, allocatable, public :: index_tetra(:,:) ! Indices of small tetrahedra
  real(dp), allocatable, public :: vibz_nodes(:,:) ! Coordinates of big tetrahedra nodes
  integer, allocatable, public :: ibz_faces(:,:) ! face indices of big tetrahedra
  real(dp), allocatable, public :: vec_face(:,:,:) ! edge vectors of big tetrahedra

  ! Variables of isosurface on IBZ
  integer, allocatable, public :: nvert(:), ntri(:) ! Output triangles and vertices of isosurface on IBZ
  real(dp), allocatable, public :: vert_coord(:,:,:) ! Output triangle vertex coordinates
  real(dp), allocatable, public :: vert_veloc(:,:,:) ! Electron velocity at vertex
  integer, allocatable, public :: vert_index(:,:,:) ! Indices of vertices on output triangles

  ! Variables of full rotated isosurface
  integer, allocatable, public :: nvert_rot(:), ntri_rot(:)
  real(dp), allocatable, public :: vert_coord_rot(:,:,:)
  real(dp), allocatable, public :: vert_veloc_rot(:,:,:)
  integer, allocatable, public :: vert_index_rot(:,:,:)


  public :: read_tetrahedra, preserve_orientation, tetra_cut, calculate_energy_sym, calculate_energy, velocity_sym, &
            velocity, calculate_star_TR, tetranodes_by_SplusG, create_isosurface_IBZ, velocity_on_IBZ, &
            rotate_IBZ_mesh, write_full_isosurface, write_IBZ_isosurface, DOS_isosurface, vertices_area
  private

contains

  subroutine read_tetrahedra()

    use intw_utility, only: find_free_unit

    implicit none

    integer :: inode, itetra, dmns, dummy1, dummy2
    integer :: io_unit, i, j


    write(*,'("| - Reading tetrahedra...                           |")')

    io_unit = find_free_unit()
    open(unit=io_unit, action="read", file="Tetrahedralized_IBZ.node", status="unknown")
    read(unit=io_unit, fmt=*) nnode, dmns, dummy1, dummy2
    allocate(ncoord(3,nnode))
    do i = 1, nnode
      read(unit=io_unit, fmt=*) inode, ( ncoord(j,i), j = 1, 3 )
    end do
    close(unit=io_unit)

    io_unit = find_free_unit()
    open(unit=io_unit, action="read", file="Tetrahedralized_IBZ.ele", status="unknown")
    read(unit=io_unit, fmt=*) ntetra, dummy1, dummy2
    allocate(index_tetra(4,ntetra))
    do i = 1, ntetra
      read(unit=io_unit, fmt=*) itetra, ( index_tetra(j,i), j = 1, 4 )
    end do
    close(unit=io_unit)

    !! JL - optional for plotting
    !! Write Tetrahedrized IBZ in off format
    io_unit = find_free_unit()
    open(unit=io_unit, action="write", file="Tetrahedralized_IBZ.off", status="unknown")
    write(unit=io_unit, fmt="(a)") "OFF"
    write(unit=io_unit, fmt="(3I6)") nnode, 4*ntetra, 0
    write(unit=io_unit, fmt=*)
    ! Vertices of tetrahedra
    do i = 1, nnode
      write(unit=io_unit, fmt="(3f12.6)") ( ncoord(j,i), j = 1, 3 )
    end do
    ! Triangular faces of tetrahedra
    do i = 1, ntetra
      write(unit=io_unit, fmt="(4I6)") 3, index_tetra(1,i), index_tetra(3,i), index_tetra(2,i)
      write(unit=io_unit, fmt="(4I6)") 3, index_tetra(1,i), index_tetra(4,i), index_tetra(2,i)
      write(unit=io_unit, fmt="(4I6)") 3, index_tetra(1,i), index_tetra(3,i), index_tetra(4,i)
      write(unit=io_unit, fmt="(4I6)") 3, index_tetra(2,i), index_tetra(3,i), index_tetra(4,i)
    end do
    close(unit=io_unit)
    !! JL - optional for plotting

    io_unit = find_free_unit()
    open(unit=io_unit, file="IBZ.off", action="read", status="unknown")
    ! Read number of irreducible tetrahedra faces
    read(unit=io_unit, fmt=*)
    read(unit=io_unit, fmt='(3I6)') nibz_nodes, nibz_faces, dummy1
    read(unit=io_unit, fmt=*)
    ! Read vectors defining the face and compute normal (the normal will be the third vector of the face)
    allocate(vibz_nodes(3, nibz_nodes), ibz_faces(3, nibz_faces))
    do i = 1, nibz_nodes
      read(unit=io_unit, fmt='(3f12.6)') ( vibz_nodes(j,i), j = 1, 3 )
    end do
    do i = 1, nibz_faces
      read(unit=io_unit, fmt='(4I6)') dummy1, ( ibz_faces(j,i), j = 1, 3 )
    end do
    allocate(vec_face(3, 3, nibz_faces))
    do i = 1, nibz_faces
      vec_face(1:3,1,i) = vibz_nodes(1:3,ibz_faces(2,i)+1) - vibz_nodes(1:3,ibz_faces(1,i)+1)
      vec_face(1:3,2,i) = vibz_nodes(1:3,ibz_faces(3,i)+1) - vibz_nodes(1:3,ibz_faces(1,i)+1)
      vec_face(1:3,3,i) = vibz_nodes(1:3,ibz_faces(1,i)+1)
    end do
    deallocate(vibz_nodes)
    close(unit=io_unit)

  end subroutine read_tetrahedra


  subroutine preserve_orientation(triangle, v)
    !ensures trangle being in the same orientation as v

    use intw_matrix_vector, only: cross

    implicit none

    real(dp), intent(inout) :: triangle(3,3), v(3)

    real(dp) :: tv(3,3)


    tv = triangle

    if ( dot_product( cross(tv(:,3)-tv(:,1), tv(2,:)-tv(:,1)), v) > 0.0_dp ) then
      triangle(1:3,1) = tv(:,1)
      triangle(1:3,2) = tv(:,2)
      triangle(1:3,3) = tv(:,3)
    else
      triangle(1:3,1) = tv(:,1)
      triangle(1:3,2) = tv(:,3)
      triangle(1:3,3) = tv(:,2)
    end if

  end subroutine preserve_orientation


  subroutine tetra_cut(ef, t, e_, nt, vt)

    use intw_utility, only: hpsort_real
    use intw_matrix_vector, only: norma, cross


    implicit none

    real(dp), intent(in)  :: ef, t(3,4), e_(4)
    integer, intent (out) :: nt
    real(dp), intent(out) :: vt(3,3,2)

    !local
    real(dp):: e(4), tv1(3), tv2(3), tv3(3), tv4(3), area1, area2, area3, area4
    integer :: ind(4), i
    integer, parameter :: n = 4
    real(dp) :: vt_(3,3,2), increasing_e(3)


    nt = 1
    vt = -100.0

    e = e_

    call hpsort_real(n, e, ind)

    increasing_e(1:3) = t(1:3,ind(4)) - t(1:3,ind(1))

    if ((ef<e(1)).or.(ef>e(4))) then
      nt = 0
      return
    else if ((ef>e(1)).and.(ef<e(2))) then
      nt = 1
      vt(1:3,1,1) = t(1:3,ind(1)) + (ef-e(1))/(e(2)-e(1)) * (t(1:3,ind(2))- (t(1:3,ind(1))))
      vt(1:3,2,1) = t(1:3,ind(1)) + (ef-e(1))/(e(3)-e(1)) * (t(1:3,ind(3))- (t(1:3,ind(1))))
      vt(1:3,3,1) = t(1:3,ind(1)) + (ef-e(1))/(e(4)-e(1)) * (t(1:3,ind(4))- (t(1:3,ind(1))))

    else if ((ef>e(1)).and.(ef>e(2)).and.(ef<e(3))) then

      nt = 2

      tv1(1:3) = t(1:3,ind(1)) + (ef-e(1))/(e(3)-e(1)) * (t(1:3,ind(3))- (t(1:3,ind(1))))
      tv2(1:3) = t(1:3,ind(1)) + (ef-e(1))/(e(4)-e(1)) * (t(1:3,ind(4))- (t(1:3,ind(1))))
      tv3(1:3) = t(1:3,ind(2)) + (ef-e(2))/(e(4)-e(2)) * (t(1:3,ind(4))- (t(1:3,ind(2))))
      tv4(1:3) = t(1:3,ind(2)) + (ef-e(2))/(e(3)-e(2)) * (t(1:3,ind(3))- (t(1:3,ind(2))))

      area1 = norma(cross(tv2(1:3)-tv1(1:3), tv3(1:3)-tv1(1:3)))
      area2 = norma(cross(tv3(1:3)-tv1(1:3), tv4(1:3)-tv1(1:3)))

      area3 = norma(cross(tv2(1:3)-tv4(1:3), tv1(1:3)-tv4(1:3)))
      area4 = norma(cross(tv2(1:3)-tv4(1:3), tv3(1:3)-tv4(1:3)))

      if ( abs(area1-area2)<= abs(area3-area4)) then

        vt(1:3,1,1) = tv1(1:3)
        vt(1:3,2,1) = tv3(1:3)
        vt(1:3,3,1) = tv2(1:3)

        vt(1:3,1,2) = tv1(1:3)
        vt(1:3,2,2) = tv4(1:3)
        vt(1:3,3,2) = tv3(1:3)

      else

        vt(1:3,1,1) = tv1(1:3)
        vt(1:3,2,1) = tv4(1:3)
        vt(1:3,3,1) = tv2(1:3)

        vt(1:3,1,2) = tv4(1:3)
        vt(1:3,2,2) = tv3(1:3)
        vt(1:3,3,2) = tv2(1:3)

      endif

    else if ((ef>e(1)).and.(ef>e(2)).and.(ef>e(3)).and.(ef<e(4))) then

      nt = 1
      vt(1:3,1,1) = t(1:3,ind(4)) + (ef-e(4))/(e(2)-e(4)) * (t(1:3,ind(2))- (t(1:3,ind(4))))
      vt(1:3,2,1) = t(1:3,ind(4)) + (ef-e(4))/(e(3)-e(4)) * (t(1:3,ind(3))- (t(1:3,ind(4))))
      vt(1:3,3,1) = t(1:3,ind(4)) + (ef-e(4))/(e(1)-e(4)) * (t(1:3,ind(1))- (t(1:3,ind(4))))

    else
      write(*, *) "error in tetra_cut. Case not treated:", ef, e
    end if

    vt_ = vt

    do i = 1, nt ! here we order the vertives accoring to the orientation.
                 ! if nt = 0 nothing is done below
      call preserve_orientation(vt(:,:,i), increasing_e )
    enddo

  end subroutine tetra_cut


  subroutine calculate_energy_sym(nsym, s, TR_sym, n_bnd, nrpts, ndegen, irvec, ham_r, kvec_int, eig_mean)

    implicit none

    integer, intent(in) :: nsym, n_bnd, nrpts, ndegen(nrpts), irvec(3,nrpts)
    integer, intent(in) :: s(3,3,nsym)
    logical, intent(in) :: TR_sym
    complex(dp), intent(in) :: ham_r(n_bnd,n_bnd,nrpts)
    real(dp), intent(in ) :: kvec_int(3)
    real(dp), intent(out) :: eig_mean(n_bnd)

    real(dp) :: eig_int(n_bnd)
    integer :: istar, nstar, symop(96,2)
    real(dp) :: vstar(3,96)


    call calculate_star_TR(nsym, s, TR_sym, mod(kvec_int(:), 1.0_dp), vstar, nstar, symop)

    eig_mean = 0.0_dp
    do istar = 1, nstar
      call calculate_energy(n_bnd, nrpts, ndegen, irvec, ham_r, vstar(1:3,istar), eig_int)
      eig_mean = eig_mean + eig_int
    enddo

    eig_mean = eig_mean/nstar

  end subroutine calculate_energy_sym


  subroutine calculate_energy(n_bnd, nrpts, ndegen, irvec, ham_r, kvec_int, eig_int)

    use intw_useful_constants, only: tpi, cmplx_i, cmplx_0

    implicit none

    external :: ZHPEVX

    integer, intent(in) :: n_bnd, nrpts, ndegen(nrpts), irvec(3,nrpts)
    complex(dp), intent(in) :: ham_r(n_bnd,n_bnd,nrpts)
    real(dp), intent(in) :: kvec_int(3)
    real(dp), intent(out) :: eig_int(n_bnd)

    complex(dp) :: ham_kprm(n_bnd,n_bnd), u_dagger(n_bnd,n_bnd), ham_pack(n_bnd*(n_bnd+1)/2), fac
    real(dp) :: rdotk
    integer :: loop_rpt, nfound
    complex(dp) :: cwork(2*n_bnd)
    real(dp) :: rwork(7*n_bnd)
    integer :: iwork(5*n_bnd), ifail(n_bnd), info
    integer :: i, j


    ham_kprm = cmplx_0
    do loop_rpt = 1, nrpts
      rdotk    = tpi*dot_product(kvec_int, irvec(:,loop_rpt))
      fac      = exp(cmplx_i*rdotk)/real(ndegen(loop_rpt), dp)
      ham_kprm = ham_kprm+fac*ham_r(:,:,loop_rpt)
    end do
    ! Diagonalise H_k (->basis of eigenstates)
    do j = 1, n_bnd
      do i = 1, j
        ham_pack(i+((j-1)*j)/2) = ham_kprm(i,j)
      enddo
    enddo

    ! Diagonalizing routine from lapack.
    call ZHPEVX('V', 'A', 'U', n_bnd, ham_pack, &
                0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
                nfound, eig_int, &
                u_dagger, &
                n_bnd, cwork, rwork, iwork, ifail, info)

  end subroutine calculate_energy


  subroutine velocity_sym(nrpts, irvec, ndegen, alat, ag, bg, nsym, s, TR_sym, n_bnd, ibnd, ham_r, v_mean, kvec_int)

    use intw_matrix_vector, only: ainv

    implicit none

    integer, intent(in) :: nrpts, irvec(3,nrpts), ndegen(nrpts)
    real(dp), intent(in) :: alat, ag(3,3), bg(3,3)
    integer, intent(in) :: nsym, s(3,3,nsym)
    logical, intent(in) :: TR_sym
    integer, intent(in) :: n_bnd, ibnd
    complex(dp), intent(in) :: ham_r(n_bnd,n_bnd,nrpts)
    real(dp), intent(in) :: kvec_int(3)
    real(dp), intent(out) :: v_mean(3)

    real(dp), parameter :: evau = 27.21138602_dp ! evau = 27.2113845
    integer :: istar, nstar, symmop(96,2) ! symmop(48)
    real(dp) :: v(3), so(3,3), kstar(3,96) ! kstar(3,48)


    !call calculate_star_r (mod(kvec_int(:), 1.0_dp), kstar, nstar, symmop)
    call calculate_star_TR(nsym, s, TR_sym, mod(kvec_int(:), 1.0_dp), kstar, nstar, symmop)

    v_mean = 0.d0
    do istar = 1, nstar
      call velocity(nrpts, irvec, ndegen, alat, ag, bg, n_bnd, ibnd, ham_r, v, kstar(:,istar))
      !v = matmul(transpose(at), v)
      !v = matmul(ainv(bg), v)
      !v = matmul(ag, v)
      v = matmul(transpose(ag), v)

      if(symmop(istar, 1).gt.0) then
        so = real(s(:,:,symmop(istar,1)), dp)
        v = matmul(ainv (so), v)

      else
        so = real(s(:,:,symmop(istar,2)), dp)
        v = -matmul(ainv (so), v)
      end if
      v = matmul(bg, v)
      v_mean = v_mean + v

    enddo
    v_mean = v_mean/nstar

  end subroutine velocity_sym


  subroutine velocity(nrpts, irvec, ndegen, alat, ag, bg, n_bnd, ibnd, ham_r, v, k)

    use intw_useful_constants, only: tpi, cmplx_i, cmplx_0
    use intw_matrix_vector, only: ainv

    implicit none

    external :: ZHPEVX

    integer, intent(in) :: n_bnd, ibnd, nrpts, irvec(3,nrpts), ndegen(nrpts)
    real(dp), intent(in) :: alat, ag(3,3), bg(3,3)
    complex(dp), intent(in) :: ham_r(n_bnd,n_bnd,nrpts)
    real(dp), intent(in) :: k(3)
    real(dp), intent(out) :: v(3)

    real(dp) :: eig(n_bnd)
    integer :: i, j, loop_rpt, nfound
    complex(dp) :: ham_kprm_x(n_bnd,n_bnd), ham_kprm_y(n_bnd,n_bnd), &
                   ham_kprm_z(n_bnd,n_bnd), ham_kprm(n_bnd,n_bnd), fac
    complex(dp) :: u_dagger(n_bnd,n_bnd), ham_pack(n_bnd*(n_bnd+1)/2)
    real(dp)    :: rdotk, rcart(3), kcart(3)
    complex(dp) :: cwork(2*n_bnd)
    real(dp) :: rwork(7*n_bnd)
    integer :: iwork(5*n_bnd), ifail(n_bnd), info
    real(dp), parameter :: evau = 27.21138602_dp ! evau = 27.2113845_dp


    ham_kprm_x = cmplx_0
    ham_kprm_y = cmplx_0
    ham_kprm_z = cmplx_0
    ham_kprm = cmplx_0

    do loop_rpt = 1, nrpts

      rdotk = tpi*dot_product(k, irvec(:,loop_rpt))

      !rcart(:) = matmul(transpose(ag), irvec(:,loop_rpt))
      !rcart(:) = matmul(ainv(bg), irvec(:,loop_rpt))
      rcart(:) = matmul(ag, irvec(:,loop_rpt))

      kcart(:) = matmul(bg, k)

      fac = exp(cmplx_i*rdotk)/real(ndegen(loop_rpt), dp)
      ham_kprm = ham_kprm +fac*ham_r(:,:,loop_rpt)

      fac = cmplx_i*rcart(1)*exp(cmplx_i*rdotk)/real(ndegen(loop_rpt), dp)
      ham_kprm_x = ham_kprm_x+fac*ham_r(:,:,loop_rpt)

      fac = cmplx_i*rcart(2)*exp(cmplx_i*rdotk)/real(ndegen(loop_rpt), dp)
      ham_kprm_y = ham_kprm_y+fac*ham_r(:,:,loop_rpt)

      fac = cmplx_i*rcart(3)*exp(cmplx_i*rdotk)/real(ndegen(loop_rpt), dp)
      ham_kprm_z = ham_kprm_z+fac*ham_r(:,:,loop_rpt)

    end do

    do j = 1, n_bnd
      do i = 1, j
        ham_pack(i+((j-1)*j)/2) = ham_kprm(i,j)
      enddo
    enddo

    call ZHPEVX('V', 'A', 'U', n_bnd, ham_pack, &
                0.0_dp, 0.0_dp, 0, 0, -1.0_dp, &
                nfound, eig, &
                u_dagger, &
                n_bnd, cwork, rwork, iwork, ifail, info)

    v = 0.0
    do i = 1, n_bnd
      do j = 1, n_bnd
        v(1) = v(1) + real(conjg( u_dagger(i,ibnd) ) * ham_kprm_x (i,j) * ( u_dagger(j,ibnd)))
        v(2) = v(2) + real(conjg( u_dagger(i,ibnd) ) * ham_kprm_y (i,j) * ( u_dagger(j,ibnd)))
        v(3) = v(3) + real(conjg( u_dagger(i,ibnd) ) * ham_kprm_z (i,j) * ( u_dagger(j,ibnd)))
      enddo
    enddo

    v = v*alat/evau

  end subroutine velocity


  subroutine calculate_star_TR(nsym, s, TR_symmetry, v, vstar, nstar, symop)

    implicit none

    integer, intent(in) :: nsym
    integer, intent(in) :: s(3,3,nsym)
    logical, intent(in) :: TR_symmetry
    real(dp), intent(in) :: v(3)
    real(dp), intent(out) :: vstar(3,96)
    integer, intent(out) :: nstar, symop(96,2)

    integer :: isym, i
    real(dp) :: vrot(3)


    symop = 0
    nstar = 1
    vstar(1:3,nstar) = v(1:3)
    symop(1,1) = 1

    sym_loop: do isym = 1, nsym

      vrot(:) = matmul(dble(s(:,:,isym)), v(:))

      do i = 1, nstar
        if ( sum(abs(vrot(:)- vstar(1:3,i))) < 10E-5 ) then
          cycle sym_loop
        end if
      enddo

      nstar = nstar + 1
      vstar(1:3,nstar) = vrot(1:3)
      symop(nstar,1) = isym

    enddo sym_loop

    if(TR_symmetry) then
      TR_sym_loop: do isym = 1, nsym
        vrot(:) = matmul(dble(s(:,:,isym)), -v(:)) ! -k

        do i = 1, nstar
          if ( sum(abs(vrot(:)- vstar(1:3,i))) < 10E-5 ) then
            cycle TR_sym_loop
          end if
        enddo

        nstar = nstar + 1
        vstar(1:3,nstar) = vrot(1:3)
        symop(nstar,2) = isym

      enddo TR_sym_loop
    end if

  end subroutine calculate_star_TR


  subroutine tetranodes_by_SplusG(icoord, jcoord, bg, nsym, s, TR_sym, epsvert, related)

    use intw_matrix_vector, only: ainv, norma

    implicit none

    ! I/O
    real(dp), intent(in) :: icoord(3), jcoord(3)
    real(dp), intent(in) :: bg(3,3)
    integer, intent(in) :: nsym, s(3,3,nsym)
    logical, intent(in) :: TR_sym
    real(dp), intent(in) :: epsvert
    logical, intent(out) :: related

    ! Local
    integer :: isym, ig, jg, kg
    real(dp) :: rot_inode(3), inode_crys(3)


    ! initalize
    related = .false.

    ! inode coord in crystal coords
    inode_crys = matmul(ainv(bg), icoord(:))
    do isym = 1, nsym
      do ig = -1, 1
        do jg = -1, 1
          do kg = -1, 1
            !
            if(ig.eq.0 .and. jg.eq.0 .and. kg.eq.0) cycle
            !
            rot_inode(1:3) = matmul(bg, &
                                        matmul(dble(s(:,:,isym)), inode_crys) &
                                        + real((/ig, jg, kg/), dp))
            !
            ! check with jnode
            if (norma(rot_inode(:)-jcoord(:)) .lt. epsvert ) then
              print*, ""
              print*, isym, ig, jg, kg, "symmetry"
              related = .true.
              return
            end if
            ! TR symmetry
            if (TR_sym) then
              !
              rot_inode(1:3) = matmul(bg, &
                                        -matmul( dble(s(:,:,isym)), inode_crys ) &
                                        + real((/ig, jg, kg/), dp))
              !
              ! check with jnode
              if (norma(rot_inode(:)-jcoord(:)) .lt. epsvert ) then
                related = .true.
                return
              end if
              !
            end if ! TR_sym
            !
          end do ! kg
        end do ! jg
      end do ! ig
    end do ! isym

  end subroutine tetranodes_by_SplusG


  subroutine create_isosurface_IBZ(eiso, n_bnd, nrpts, ndegen, irvec, ham_r, alat, ag, bg, nsym, s, TR_sym, verbose, epsvert)

    use intw_matrix_vector, only: ainv, norma

    implicit none

    ! I/O
    real(dp), intent(in) :: eiso ! Eenrgy of isosurface
    integer, intent(in) :: n_bnd ! Numer of bands on hr
    integer, intent(in) :: nrpts, ndegen(nrpts) ! Number of R vectors and their degeneracies on hr
    integer, intent(in) :: irvec(3,nrpts) ! R vectors
    complex(dp), intent(in) :: ham_r(n_bnd,n_bnd,nrpts) ! Hamiltonian in Wannier representation
    real(dp), intent(in) :: alat, ag(3,3), bg(3,3) ! Crystal and reciprocal lattice vectors on columns
    integer, intent(in) :: nsym ! Number of symmetry operations
    integer, intent(in) :: s(3,3,nsym) ! Symmetry-operation matrices
    logical, intent(in) :: TR_sym ! Time-reversal symmetry or not
    logical, intent(in) :: verbose ! .true. if info files are to be written
    real(dp), intent(in) :: epsvert ! Parameter to detect duplicated vertices

    ! Local
    integer :: itet, inod, ibnd, itri, j, iv, ivert, rdcd_vert, mid_nvert, tot_nvert, tot_ntri, ntr
    real(dp) :: tetracoord(3,4), vcoord(3,30000,n_bnd), etetra(ntetra,4,n_bnd), eig_int(n_bnd)
    real(dp) :: kvec_int(3), et(4), vtr(3,3,2)
    integer :: vindex(3,10000,n_bnd)
    ! integer :: jtet, jnod
    ! logical :: SplusG_node
    ! logical, allocatable :: node_done(:,:)


    write(*,'("| - Creating isosurface at ",F14.8," eV...     |")') eiso

    ! Calculate energies
    !allocate(node_done(4,ntetra))
    !node_done(:,:) = .false.
    do itet = 1, ntetra
      inod_loop: do inod = 1, 4
        !
        !if(node_done(inod,itet)) cycle
        !
        tetracoord(:,inod) = ncoord(:,index_tetra(inod,itet)+1)
        kvec_int(:) = matmul(ainv(bg), tetracoord(:,inod)) ! Tetra node in crystal coordinates
        !
        call calculate_energy_sym(nsym, s, TR_sym, n_bnd, nrpts, ndegen, irvec, ham_r, kvec_int, eig_int)
        !call calculate_energy(n_bnd, nrpts, ndegen, irvec, ham_r, kvec_int, eig_int)
        etetra(itet,inod,1:n_bnd) = eig_int(1:n_bnd)
        !
        !! Check if this node has S+G related node
        !if (ANY(abs(matmul(ainv(bg), tetracoord(:,inod))-1.0_dp).lt.epsvert)) then !only check vertices at border of BZ
        !  do jtet = 1, ntetra
        !    do jnod = 1, 4
        !      !
        !      if (ALL(abs(matmul(ainv(bg), ncoord(:,index_tetra(jnod, jtet)+1))-1.0_dp).gt.epsvert)) cycle
        !      if(inod.eq.jnod .and. itet.eq.jtet) cycle
        !      if(node_done(jnod, jtet)) cycle
        !      !
        !      call tetranodes_by_SplusG(tetracoord(:,inod), ncoord(:,index_tetra(jnod,jtet)+1), bg, nsym, s, TR_sym, epsvert, SplusG_node)
        !      if(SplusG_node) then ! {jnod,jtet} is pair of {inod,itet}
        !        print*, itet, inod, jtet, jnod, "related!"
        !        etetra(jtet, jnod,1:n_bnd) = etetra(itet,inod,1:n_bnd)
        !        node_done(jnod,jtet) = .true.
        !        !cycle inod_loop
        !      end if
        !    end do
        !  end do
        !  !
        !end if
        !
      end do inod_loop
    end do
    !deallocate(node_done)

    ! Create isosurface
    allocate(nvert(n_bnd), ntri(n_bnd))
    nvert = 0
    ntri = 0
    mid_nvert = 0
    do ibnd = 1, n_bnd
      do itet = 1, ntetra
        et(1:4) = etetra(itet,1:4,ibnd)
        do inod = 1, 4
          tetracoord(:,inod) = ncoord(:,index_tetra(inod,itet)+1)
        end do
        ntr = 0
        vtr = 0
        call tetra_cut(eiso, tetracoord, et, ntr, vtr)
        do itri = 1, ntr
          ntri(ibnd) = ntri(ibnd)+1
          do iv = 1, 3
            nvert(ibnd) = nvert(ibnd)+1
            vcoord(1:3,nvert(ibnd),ibnd) = vtr(1:3,iv,itri)
            vindex(iv,ntri(ibnd),ibnd) = nvert(ibnd)+mid_nvert
          end do
        end do
      end do
      mid_nvert = mid_nvert+nvert(ibnd)
    end do

    if(verbose) then
      write(*,'("|   Number of triangules and vertices:              |")')
      do ibnd = 1, n_bnd
        write(*,'("|     band ",I4,": ",I6,I6,"                       |")') ibnd, ntri(ibnd), nvert(ibnd)
      end do
    endif

    ! Total number of triangles and vertices
    tot_ntri = 0
    tot_nvert = 0
    do ibnd = 1, n_bnd
      tot_ntri = tot_ntri + ntri(ibnd)
      tot_nvert = tot_nvert + nvert(ibnd)
    end do

    if(verbose) write(*,'("|         ---------------------------------         |")')

    ! Clean list
    write(*,'("| - Cleaning list of triangles and vertices...      |")')
    !
    allocate(vert_index(3,tot_ntri,n_bnd))
    allocate(vert_coord(3,tot_nvert,n_bnd))
    vert_coord(:,:,:) = 0.0_dp
    vert_index(:,:,:) = 0
    rdcd_vert = 0
    ivert = 0
    mid_nvert = 0
    do ibnd = 1, n_bnd
      if(ntri(ibnd).eq.0) cycle
      vert_coord(:,1,ibnd) = vcoord(:,1,ibnd)
      ivert = 0
      rdcd_vert = 0
      do itri = 1, ntri(ibnd)
        v_loop: do iv = 1, 3
          ivert = ivert+1
          do j = 1, rdcd_vert
            if(norma(vcoord(:,ivert,ibnd)-vert_coord(:,j,ibnd)).lt.epsvert) then
              vert_index(iv,itri,ibnd) = j
              cycle v_loop
            end if
          end do
          !
          rdcd_vert = rdcd_vert+1
          vert_index(iv,itri,ibnd) = rdcd_vert
          vert_coord(:,rdcd_vert,ibnd) = vcoord(:,ivert,ibnd)
          !
        end do v_loop
      end do ! itri
      nvert(ibnd) = rdcd_vert
    end do ! ibnd

    if(verbose) then
      write(*,'("|   Number of triangles and vertices:               |")')
      do ibnd = 1, n_bnd
        if(ntri(ibnd).eq.0) cycle
        write(*,'("|     band ",I4,": ",I6,I6,"                       |")') ibnd, ntri(ibnd), nvert(ibnd)
      end do
    endif

    ! Compute velocity on IBZ
    allocate(vert_veloc(3, tot_nvert, n_bnd))
    vert_veloc(:,:,:) = 0.0_dp
    call velocity_on_IBZ(n_bnd, nvert, vert_coord, nrpts, irvec, ndegen, alat, ag, bg, nsym, s, TR_sym, ham_r, vert_veloc)

    ! Total number of triangles and vertices
    tot_ntri = 0
    tot_nvert = 0
    do ibnd = 1, n_bnd
      tot_ntri = tot_ntri + ntri(ibnd)
      tot_nvert = tot_nvert + nvert(ibnd)
    end do

  end subroutine create_isosurface_IBZ


  subroutine velocity_on_IBZ(n_bnd, n_vert, v_coord, nrpts, irvec, ndegen, alat, ag, bg, nsym, s, TR_sym, ham_r, v_veloc)

    use intw_matrix_vector, only: ainv

    implicit none

    ! I/O
    integer, intent(in) :: n_bnd, n_vert(n_bnd)
    real(dp), intent(in) :: v_coord(:,:,:)
    integer, intent(in) :: nrpts, irvec(nrpts), ndegen(nrpts)
    real(dp), intent(in) :: alat, ag(3,3), bg(3,3)
    integer, intent(in) :: nsym, s(3,3,nsym)
    logical, intent(in) :: TR_sym
    complex(dp), intent(in) :: ham_r(n_bnd,n_bnd,nrpts)
    real(dp), intent(out) :: v_veloc(:,:,:)

    ! Local
    integer :: ibnd, iv
    real(dp) :: vcoord_crys(3), v_k(3)


    v_veloc = 0.0_dp
    do ibnd = 1, n_bnd
      !
      if(n_vert(ibnd).eq.0) cycle
      !
      do iv = 1, n_vert(ibnd)
        v_k = 0.0_dp
        vcoord_crys(:) = matmul(ainv(bg), v_coord(:,iv,ibnd))
        call velocity_sym(nrpts, irvec, ndegen, alat, ag, bg, nsym, s, TR_sym, n_bnd, ibnd, ham_r, v_k, vcoord_crys)
        v_veloc(:,iv,ibnd) = v_k(:)
      end do ! iv
    end do ! ibnd

  end subroutine velocity_on_IBZ


  subroutine rotate_IBZ_mesh(bg, nsym, s, n_vert, n_tri, v_coord, v_veloc, v_index, nstar, symop, epsvert, n_vert_rot, n_tri_rot, v_coord_rot, v_veloc_rot, v_index_rot)

    use intw_matrix_vector, only: ainv, norma

    implicit none

    ! I/O variables
    real(dp), intent(in) :: bg(3,3)
    integer, intent(in) :: nsym, s(3,3,nsym)
    integer, intent(in) :: n_vert, n_tri
    real(dp), intent(in) :: v_coord(3,n_vert)
    real(dp), intent(in) :: v_veloc(3,n_vert)
    integer, intent(in) :: v_index(3,n_tri)
    integer, intent(in) :: nstar
    integer, intent(in) :: symop(92,2)
    real(dp), intent(in) :: epsvert ! Parameter to detect duplicated vertices
    integer, intent(out) :: n_vert_rot, n_tri_rot
    real(dp), intent(out) :: v_coord_rot(:,:)
    real(dp), intent(out) :: v_veloc_rot(:,:)
    integer, intent(out) :: v_index_rot(:,:)

    ! Local variables
    integer :: iv, ivert, i, j, it, istar, i_sym, nvert_aux
    real(dp) :: v(3), bgi(3,3)
    real(dp) :: vcoord_aux(3,nstar*n_vert), veloc_aux(3,nstar*n_vert)
    integer :: vindex_aux(3,nstar*n_tri)
    ! real(dp) :: s_cart(3,3,96)


    bgi = ainv(bg)

    !! transform symmetry operations to cartesian coords
    !do istar = 1, nstar
    !  if(symop(istar,1).ne.0) then
    !    i_sym = symop(istar,1)
    !    do i = 1, 3
    !      s_cart(:,i,i_sym) = matmul(bg, s(:,i,i_sym))
    !    end do
    !  end if
    !end do

    nvert_aux = 0
    n_tri_rot = 0
    do istar = 1, nstar
      if ( symop(istar,1) /= 0 ) then
        i_sym = symop(istar,1)
        do iv = 1, n_vert
          nvert_aux = nvert_aux + 1
          v(1:3) = matmul(bgi, v_coord(1:3,iv))
          !vcoord_aux(1:3,nvert_aux) = matmul(real(s(:,:,i_sym), dp), v)
          vcoord_aux(1:3,nvert_aux) = matmul(s(:,:,i_sym), v)
          !! Transformation in cartesian coords
          !vcoord_aux(1:3,nvert_aux) = matmul(s_cart(:,:,i_sym), v_coord(1:3,iv))
          !
          v(1:3) = matmul(bgi, v_veloc(1:3,iv))
          !veloc_aux(:,nvert_aux) = matmul(real(s(:,:,i_sym), dp), v)
          veloc_aux(:,nvert_aux) = matmul(s(:,:,i_sym), v)
        end do
        do it = 1, n_tri
          n_tri_rot = n_tri_rot + 1
          vindex_aux(1:3,n_tri_rot) = v_index(1:3,it) + (istar-1)*n_vert
        end do
      else if( symop(istar,2) /= 0 ) then ! TR symmetry
        i_sym = symop(istar,2)
        do iv = 1, n_vert
          nvert_aux = nvert_aux + 1
          v(1:3) = matmul(bgi, v_coord(1:3,iv))
          !vcoord_aux(1:3,nvert_aux) = -matmul(real(s(:,:,i_sym), dp), v)
          vcoord_aux(1:3,nvert_aux) = -matmul(s(:,:,i_sym), v)
          !! Transformation in cartesian coords
          !vcoord_aux(1:3,nvert_aux) = -matmul(s_cart(:,:,i_sym), v_coord(1:3,iv))
          !
          v(1:3) = matmul(bgi, v_veloc(1:3,iv))
          !veloc_aux(:,nvert_aux) = -matmul(real(s(:,:,i_sym), dp), v)
          veloc_aux(:,nvert_aux) = -matmul(s(:,:,i_sym), v)
        end do
        do it = 1, n_tri
          n_tri_rot = n_tri_rot + 1
          vindex_aux(1:3,n_tri_rot) = v_index(1:3,it) + (istar-1)*n_vert
        end do
      end if
    end do

    ! Clean list
    n_vert_rot = 0
    ivert = 0
    v_coord_rot(:,1) = vcoord_aux(:,1)
    do i = 1, n_tri_rot
      v_loop: do iv = 1, 3
        do j = 1, n_vert_rot
          if ( norma(vcoord_aux(:,vindex_aux(iv,i))-v_coord_rot(:,j)) < epsvert ) then
            v_index_rot(iv,i) = j
            cycle v_loop
          end if
        end do
        n_vert_rot = n_vert_rot+1
        v_index_rot(iv,i) = n_vert_rot
        v_coord_rot(:,n_vert_rot) = vcoord_aux(:,vindex_aux(iv,i))
        v_veloc_rot(:,n_vert_rot) = veloc_aux(:,vindex_aux(iv,i))
      end do v_loop
    end do

    ! Transform back to cartesian coordinates
    do j = 1, n_vert_rot
      v_coord_rot(1:3,j) = matmul(bg, v_coord_rot(1:3,j))
      v_veloc_rot(1:3,j) = matmul(bg, v_veloc_rot(1:3,j))
    end do

  end subroutine rotate_IBZ_mesh


  subroutine write_full_isosurface(bg, nsym, s, TR_sym, n_bnd, verbose, epsvert, tag, prefix)

    use intw_utility, only: find_free_unit, int2str
    use intw_matrix_vector, only: norma

    implicit none

    ! I/O
    real(dp), intent(in) :: bg(3,3)
    integer, intent(in) :: nsym
    integer, intent(in) :: s(3,3,nsym)
    logical, intent(in) :: TR_sym
    integer, intent(in) :: n_bnd
    logical, intent(in) :: verbose
    real(dp), intent(in) :: epsvert
    character(len=*), intent(in) :: tag
    character(len=*), optional, intent(in) :: prefix

    ! Local
    character(len=256) :: filename
    integer :: nstar, symop(96,2)
    real(dp) :: vstar(3,96)
    integer :: tot_ntri, mid_nvert, tot_nvert
    integer :: ibnd, iv, i, j, io_unit


    write(*,'("| - Rotating IBZ...                                 |")')

    ! Calculate star
    do iv = 1, nnode
      call calculate_star_TR(nsym, s, TR_sym, ncoord(:,iv), vstar, nstar, symop)
      if(nstar.eq.nsym) exit
    end do

    if(verbose) write(*, '("|   nstar = ",I6,"                   |")') nstar

    ! Rotate mesh
    tot_ntri = sum(ntri(:))
    tot_nvert = sum(nvert(:))
    allocate(vert_index_rot(3,nstar*tot_ntri,n_bnd))
    allocate(vert_coord_rot(3,nstar*tot_nvert,n_bnd))
    allocate(vert_veloc_rot(3,nstar*tot_nvert,n_bnd))
    allocate(ntri_rot(n_bnd), nvert_rot(n_bnd))
    ntri_rot = 0
    nvert_rot = 0
    do ibnd = 1, n_bnd
      if(nvert(ibnd).eq.0) cycle
      call rotate_IBZ_mesh(bg, nsym, s, nvert(ibnd), ntri(ibnd), vert_coord(:,:,ibnd), &
                           vert_veloc(:,:,ibnd), vert_index(:,:,ibnd), nstar, symop, epsvert, &
                           nvert_rot(ibnd), ntri_rot(ibnd), vert_coord_rot(:,:,ibnd), &
                           vert_veloc_rot(:,:,ibnd), vert_index_rot(:,:,ibnd))
    end do

    ! Total number of triangles and vertices
    tot_ntri = 0
    tot_nvert = 0
    do ibnd = 1, n_bnd
      tot_ntri = tot_ntri + ntri_rot(ibnd)
      tot_nvert = tot_nvert + nvert_rot(ibnd)
    end do

    ! Write rotated mesh in OFF format
    write(*,'("| - Writing full FS in OFF format...                |")')

    if (present(prefix)) then
      filename = trim(prefix)//"."//trim(tag)//".off"
    else
      filename = trim(tag)//".off"
    endif
    write(*,'("|   Filename: ",A," |")') filename(1:max(37,len(trim(filename))))

    !
    ! Write full isosurface in OFF format
    io_unit = find_free_unit()
    open(unit=io_unit, action="write", file=filename, status="unknown")
    write(unit=io_unit, fmt="(a)") "OFF"
    write(unit=io_unit, fmt="(3I10)") tot_nvert, tot_ntri, 0
    write(unit=io_unit, fmt=*)
    do ibnd = 1, n_bnd
      do iv = 1, nvert_rot(ibnd)
        write(unit=io_unit, fmt="(4f18.10)") ( vert_coord_rot(j,iv,ibnd), j = 1, 3 )
      end do
    end do
    mid_nvert = 0
    do ibnd = 1, n_bnd
      do i = 1, ntri_rot(ibnd)
        write(unit=io_unit, fmt="(4I6)") 3, ( mid_nvert+vert_index_rot(j,i,ibnd)-1, j = 1, 3 )
      end do
      mid_nvert = mid_nvert + nvert_rot(ibnd)
    end do
    close(unit=io_unit)

    write(*,'("|   Total number of triangles and vertices:         |")')
    write(*,'("|     Triangles: ",I6,"                             |")') tot_ntri
    write(*,'("|     Vertices:  ",I6,"                             |")') tot_nvert

    !
    ! Write band-separated isosurface in OFF format
    do ibnd = 1, n_bnd
      if(nvert(ibnd).eq.0) cycle
      io_unit = find_free_unit()
      if (present(prefix)) then
        filename=trim(prefix)//"."//trim(int2str(ibnd))//"_"//trim(tag)//".off"
      else
        filename=trim(int2str(ibnd))//"_"//trim(tag)//".off"
      endif
      open(unit=io_unit, action="write", file=filename, status="unknown")
      write(unit=io_unit, fmt="(a)") "OFF"
      write(unit=io_unit, fmt="(3I10)") nvert_rot(ibnd), ntri_rot(ibnd), 0
      write(unit=io_unit, fmt=*)
      do iv = 1, nvert_rot(ibnd)
        write(unit=io_unit, fmt="(4f18.10)") ( vert_coord_rot(j,iv,ibnd), j = 1, 3 )
      end do
      do i = 1, ntri_rot(ibnd)
        write(unit=io_unit, fmt="(4I6)") 3, ( vert_index_rot(j,i,ibnd)-1, j = 1, 3)
      end do
      close(unit=io_unit)
    enddo

    !
    ! Write files with Fermi velocity
    write(*, '("| - Writing Fermi velocity files...                 |")')
    do ibnd = 1, n_bnd
      if(nvert(ibnd).eq.0) cycle
      io_unit = find_free_unit()
      if (present(prefix)) then
        filename=trim(prefix)//"."//trim(int2str(ibnd))//"_"//trim(tag)//"_v_k.dat"
      else
        filename=trim(int2str(ibnd))//"_"//trim(tag)//"_v_k.dat"
      endif
      open(unit=io_unit, action="write", file=filename, status="unknown")
      write(unit=io_unit, fmt="(3I10)") nvert_rot(ibnd)
      do iv = 1, nvert_rot(ibnd)
        write(unit=io_unit, fmt="(I6,4F18.10)") iv, ( vert_veloc_rot(j,iv,ibnd), j = 1, 3 ), norma(vert_veloc_rot(:,iv,ibnd))
      end do
      close(unit=io_unit)
    end do

  end subroutine write_full_isosurface


  subroutine write_IBZ_isosurface(tag, n_bnd, velocities, prefix)

    use intw_utility, only: find_free_unit, int2str
    use intw_matrix_vector, only: norma

    implicit none

    ! I/O
    character(len=*), intent(in) :: tag
    integer, intent(in) :: n_bnd
    logical, intent(in) :: velocities
    character(len=*), optional, intent(in) :: prefix

    ! Local
    character(len=256) :: filename
    integer :: ibnd, io_unit, iv, itri, i, j
    integer :: tot_nvert, tot_ntri, mid_nvert


    write(*,'("| - Writing IBZ FS in OFF format...                 |")')

    if (present(prefix)) then
      filename = trim(prefix)//"."//trim(tag)//".off"
    else
      filename = trim(tag)//".off"
    endif

    write(*,'("|   Filename: ",A," |")') filename(1:max(37,len(trim(filename))))

    !
    ! Write IBZ isosurface in OFF format
    io_unit = find_free_unit()
    open(unit=io_unit, action="write", file=filename, status="unknown")
    write(unit=io_unit, fmt="(a)") "OFF"
    write(unit=io_unit, fmt="(3I10)") tot_nvert, tot_ntri, 0
    write(unit=io_unit, fmt=*)
    mid_nvert = 0
    do ibnd = 1, n_bnd
      if(nvert(ibnd).eq.0) cycle
      do iv = 1, nvert(ibnd)
        write(unit=io_unit, fmt="(3F18.10)") ( vert_coord(j,iv,ibnd), j = 1, 3 )
      end do
    end do
    do ibnd = 1, n_bnd
      do itri = 1, ntri(ibnd)
        write(unit=io_unit, fmt="(4I6)") 3, ( mid_nvert+vert_index(j,itri,ibnd)-1, j = 1, 3 )
      end do
      mid_nvert = mid_nvert + nvert(ibnd)
    end do
    close(unit=io_unit)

    !
    ! Write band-separated IBZ isosurface in OFF format
    tot_nvert = 0
    tot_ntri = 0
    do ibnd = 1, n_bnd
      if(nvert(ibnd).eq.0) cycle
      io_unit = find_free_unit()
      if (present(prefix)) then
        filename=trim(prefix)//"."//trim(int2str(ibnd))//"_"//trim(tag)//".off"
      else
        filename=trim(int2str(ibnd))//"_"//trim(tag)//".off"
      endif
      open(unit=io_unit, action="write", file=filename, status="unknown")
      write(unit=io_unit, fmt="(a)") "OFF"
      write(unit=io_unit, fmt="(3I10)") nvert(ibnd), ntri(ibnd), 0
      write(unit=io_unit, fmt=*)
      do iv = 1, nvert(ibnd)
        write(unit=io_unit, fmt="(4f18.10)") ( vert_coord(j,iv,ibnd), j = 1, 3)
        tot_nvert = tot_nvert + 1
      end do
      do i = 1, ntri(ibnd)
        write(unit=io_unit, fmt="(4I6)") 3, ( vert_index(j,i,ibnd)-1, j = 1, 3)
        tot_ntri = tot_ntri + 1
      end do
      close(unit=io_unit)
    end do

    write(*,'("|   Total number of triangles and vertices:         |")')
    write(*,'("|     Triangles: ",I6,"                             |")') tot_ntri
    write(*,'("|     Vertices:  ",I6,"                             |")') tot_nvert

    !
    ! Write file with Fermi velocity
    if (velocities) then
      !
      write(*,'("| - Writing Fermi velocity files...                 |")')
      !
      do ibnd = 1, n_bnd
        if(nvert(ibnd).eq.0) cycle
        io_unit = find_free_unit()
        if (present(prefix)) then
          filename=trim(prefix)//"."//trim(int2str(ibnd))//"_"//trim(tag)//"_v_k.dat"
        else
          filename=trim(int2str(ibnd))//"_"//trim(tag)//"_v_k.dat"
        endif
        open(unit=io_unit, action="write", file=filename, status="unknown")
        write(unit=io_unit, fmt="(3I10)") nvert(ibnd)
        do iv = 1, nvert(ibnd)
          write(unit=io_unit, fmt="(I6,4F18.10)") iv, ( vert_veloc(j,iv,ibnd), j = 1, 3 ), norma(vert_veloc(:,iv,ibnd))
        end do
        close(unit=io_unit)
      end do
      !
    endif

  end subroutine write_IBZ_isosurface


  subroutine DOS_isosurface(alat, vbz, n_bnd, n_vert, n_tri, v_coord, v_index, v_veloc, n_vert_rot, n_tri_rot, v_coord_rot, v_index_rot, v_veloc_rot)

    use intw_utility, only: find_free_unit, int2str
    use intw_useful_constants, only: pi
    use intw_matrix_vector, only: norma

    implicit none

    ! I/O
    real(dp), intent(in) :: alat, vbz
    integer, intent(in) :: n_bnd, n_vert(:), n_tri(:), v_index(:,:,:), n_vert_rot(:), n_tri_rot(:), v_index_rot(:,:,:)
    real(dp), intent(in) :: v_coord(:,:,:), v_veloc(:,:,:), v_coord_rot(:,:,:), v_veloc_rot(:,:,:)

    ! Local
    integer :: ibnd, iv, io_unit
    real(dp) :: ne
    real(dp), allocatable :: vert_area(:,:), ifs_area(:), rot_vert_area(:,:), fs_area(:)


    ! Allocate variables
    allocate(vert_area(maxval(n_vert(:)),n_bnd), ifs_area(n_bnd))
    allocate(rot_vert_area(maxval(n_vert_rot(:)),n_bnd), fs_area(n_bnd))

    io_unit = find_free_unit()
    open(unit=io_unit, action="write", file="DOS.dat", status="unknown")
    ! Loop over all FS sheets
    do ibnd = 1, n_bnd
      !
      if(n_vert(ibnd).eq.0) cycle
      !
      ! Area of each vertex on IFS
      call vertices_area(n_vert(ibnd), n_tri(ibnd), v_coord(:,:,ibnd), v_index(:,:,ibnd), vert_area(:,ibnd))
      ! Area of IFS
      ifs_area = 0.0_dp
      do iv = 1, n_vert(ibnd)
        ifs_area(ibnd) = ifs_area(ibnd) + vert_area(iv,ibnd)
      end do
      write(io_unit, *) "Area of "//trim(int2str(ibnd))//" IFS =", ifs_area(ibnd), "(2*pi/alat)**2"
      write(io_unit, *) "              ", ifs_area(ibnd)*(2.0_dp*pi/alat)**2, "bohr**-2"
      write(io_unit, *)

      ! Area of each vertex on full FS
      call vertices_area(n_vert_rot(ibnd), n_tri_rot(ibnd), v_coord_rot(:,:,ibnd), v_index_rot(:,:,ibnd), rot_vert_area(:,ibnd))
      ! Area of IFS
      fs_area = 0.0_dp
      do iv = 1, n_vert_rot(ibnd)
        fs_area(ibnd) = fs_area(ibnd) + rot_vert_area(iv,ibnd)
      end do
      write(io_unit, *) "Area of full "//trim(int2str(ibnd))//" FS =", fs_area(ibnd), "(2*pi/alat)**2"
      write(io_unit, *) "                  ", fs_area(ibnd)*(2.0_dp*pi/alat)**2, "bohr**-2"
      write(io_unit, *)


      write(io_unit, *) "full FS area / IFS area =", fs_area(ibnd) / ifs_area(ibnd)
      write(io_unit, *)

    end do ! ibnd

    !--------------- Compute Density of States at FS, N(e_F)
    !
    write(io_unit, *)
    write(io_unit, *) "Volume of BZ (bohr**-3) =", (vbz*(2.0_dp*pi/alat)**3)
    write(io_unit, *) ""

    do ibnd = 1, n_bnd
      if(n_vert_rot(ibnd).eq.0) cycle
      ne = 0.0_dp
      do iv = 1, n_vert_rot(ibnd)
        ne = ne + (1.0_dp/norma(v_veloc_rot(1:3,iv,ibnd))) * (rot_vert_area(iv,ibnd)*(2.0_dp*pi/alat)**2)
      end do
      ne = 2.0_dp * ne / (vbz*(2.0_dp*pi/alat)**3)
      write(io_unit, *) "Partial DOS at FS branch"//trim(int2str(ibnd))//" = ", ne, "a.u."
      write(io_unit, *) "                             ", ne/27.21138602_dp, "1/e.V."
      write(io_unit, *) ""
    end do

    ne = 0.0_dp
    do ibnd = 1, n_bnd
      if(n_vert(ibnd).eq.0) cycle
      do iv = 1, n_vert(ibnd)
        ne = ne + (1.0_dp/norma(v_veloc(1:3,iv,ibnd))) * (vert_area(iv,ibnd)*(2.0_dp*pi/alat)**2)
      end do
    end do
    ne = ne * sum(fs_area(:)) / sum(ifs_area(:))
    ne = 2.0_dp * ne / (vbz*(2.0_dp*pi/alat)**3)

    write(io_unit, *) "DOS at FS using irreducible mesh =", ne, "a.u."
    write(io_unit, *) "                                  ", ne/27.21138602_dp, "1/e.V."
    write(io_unit, *)

    ne = 0.0_dp
    do ibnd = 1, n_bnd
      if(n_vert_rot(ibnd).eq.0) cycle
      do iv = 1, n_vert_rot(ibnd)
        ne = ne + (1.0_dp/norma(v_veloc_rot(1:3,iv,ibnd))) * (rot_vert_area(iv,ibnd)*(2.0_dp*pi/alat)**2)
      end do
    end do

    ne = 2.0_dp * ne / (vbz*(2.0_dp*pi/alat)**3)

    write(io_unit, *) "DOS at FS using full rotated mesh =", ne, "a.u."
    write(io_unit, *) "                                   ", ne/27.21138602_dp, "1/e.V."

    close(io_unit)

    write(*,'("|   DOS at FS: ",f16.8, " a.u.                |")') ne
    write(*,'("|              ",f16.8, " 1/eV                |")') ne/27.21138602_dp

    !--------------- Deallocate variables
    !
    deallocate(ifs_area, fs_area, vert_area, rot_vert_area)


  end subroutine DOS_isosurface


  subroutine vertices_area(n_vert, n_tri, v_coord, v_index, v_area)

    use intw_matrix_vector, only: area_vec
    implicit none

    ! I/O variables
    integer, intent(in) :: n_vert, n_tri
    real(dp), intent(in) :: v_coord(3,n_vert)
    integer, intent(in) :: v_index(3,n_tri)
    real(dp), intent(out) :: v_area(n_vert)

    ! Local variables
    integer, parameter :: ntmax = 100
    integer :: iv, jv, it, jt, ivertex
    integer, dimension(n_vert) :: num_of_nghb
    integer, dimension(n_vert,ntmax) :: tri_of_vertex
    real(dp), dimension(3) :: v1, v2, v3
    real(dp) :: tri_area


    ! Find neighbour triangles of each vertex
    num_of_nghb(:) = 0
    do it = 1, n_tri
      do jv = 1, 3
        ivertex = v_index(jv,it)
        num_of_nghb(ivertex) = num_of_nghb(ivertex) + 1
        if(num_of_nghb(ivertex) .le. ntmax) then
          tri_of_vertex(ivertex,num_of_nghb(ivertex)) = it
        else
          write(*, *) "Error", ivertex ; stop
        end if
      end do
    end do

    ! Compute area of each vertex
    v_area = 0.0_dp
    do iv = 1, n_vert
      do jt = 1, num_of_nghb(iv)
        it = tri_of_vertex(iv,jt)
        v1(:) = v_coord(:,v_index(1,it))
        v2(:) = v_coord(:,v_index(2,it))
        v3(:) = v_coord(:,v_index(3,it))
        tri_area = area_vec(v2-v1, v3-v1)

        v_area(iv) = v_area(iv) + tri_area/3.0_dp

      end do
    end do

  end subroutine vertices_area

end module triFS_isosurface
