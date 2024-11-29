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
module triFS_mesh_opt

  ! TODO: Add a description.

  use kinds, only: dp

  ! Variables of isosurface on IBZ
  integer, allocatable, public :: opt_nvert(:), opt_ntri(:) ! Optimized triangles and vertices of isosurface on IBZ
  real(dp), allocatable, public :: opt_vert_coord(:,:,:) ! Optimized triangle vertex coordinates
  integer, allocatable, public :: opt_vert_index(:,:,:) ! Indices of vertices on optimized triangles

  public :: mesh_optimization, newton_rap, check_border, vertices_related_by_SplusG, rotate_to_SplusG, &
            v2edge, edge_collapse, faces_vectors, inner_faces, find_neighbours, find_neighbours_vertex_list, &
            rm_3tri_vrtx, collapse_triangles, collapse_condition, tangential_relaxation, mean_barycenter
  private

contains

  subroutine mesh_optimization(collapse, relax, newton_raphson, collapse_criteria, relax_iter, newton_iter, &
                               relax_vinface, eps_vinface, eps_dupv, verbose, ef, nrpts, irvec, ndegen, ham_r, &
                               alat, ag, bg, nsym, s, TR_sym, num_wann, nfaces_IBZ, faces_IBZ_as_vert, vert_IBZ, ntri, &
                               nvert, v_coord, v_index)

    use triFS_isosurface, only: write_IBZ_isosurface


    implicit none

    ! I/O
    logical, intent(in) :: collapse, relax ! Mesh optimization options
    integer, intent(in) :: newton_raphson
    real(dp), intent(in) :: collapse_criteria
    integer, intent(in) :: relax_iter, newton_iter
    logical, intent(in) :: relax_vinface
    real(dp), intent(in) :: eps_vinface, eps_dupv
    logical, intent(in) :: verbose
    real(dp), intent(in) :: ef ! Energy of isosurface
    integer, intent(in) :: nrpts, irvec(1:3,nrpts), ndegen(nrpts) ! Wannier/Fourier space variables
    complex(dp), intent(in) :: ham_r(num_wann,num_wann,nrpts) ! Hamiltonian in Wannier basis
    real(dp), intent(in) :: alat, ag(1:3,1:3), bg(1:3,1:3) ! real and reciprocal lattice vectors
    integer, intent(in) :: nsym, s(1:3,1:3,nsym) ! symmetry variables
    logical, intent(in) :: TR_sym ! time-reversal symmetry
    integer, intent(in) :: nfaces_IBZ, num_wann
    integer, intent(in) :: faces_IBZ_as_vert(:,:) ! vertex indeces of IBZ big tetra faces
    real(dp), intent(in) :: vert_IBZ(:,:) ! IBZ big tetra vertex coordinates
    integer, intent(inout) :: ntri(:), nvert(:)
    real(dp), dimension(:,:,:), intent(inout) :: v_coord
    integer, dimension(:,:,:), intent(inout) :: v_index

    ! Local
    character(len=25) :: tag_in
    integer :: ibnd


    !!--------------- Newton-Raphson relaxation before improvement (optional)
    !!
    if(newton_raphson.eq.2) then
      call newton_rap(eps_dupv, eps_vinface, ef, nrpts, irvec, ndegen, ham_r, alat, ag, bg, &
                      nsym, s, TR_sym, num_wann, newton_iter, nfaces_IBZ, faces_IBZ_as_vert, &
                      vert_IBZ, ntri, nvert, v_coord)
    end if

    !--------------- Edge-collapse
    !
    if(collapse) then
      allocate(opt_nvert(num_wann),opt_ntri(num_wann))
      allocate(opt_vert_coord(1:3,maxval(nvert(:)),num_wann),opt_vert_index(1:3,maxval(ntri(:)),num_wann))
      opt_nvert = 0
      opt_ntri = 0
      do ibnd = 1, num_wann
        if (nvert(ibnd).eq.0) cycle
        !
        call edge_collapse(collapse_criteria, nvert(ibnd), ntri(ibnd), nfaces_IBZ, &
                           faces_IBZ_as_vert, vert_IBZ, v_coord(:,:,ibnd), v_index(:,:,ibnd), &
                           opt_nvert(ibnd), opt_ntri(ibnd), opt_vert_coord(:,:,ibnd), &
                           opt_vert_index(:,:,ibnd))
        !
      end do
      !
      nvert = opt_nvert
      ntri = opt_ntri
      v_coord = opt_vert_coord
      v_index = opt_vert_index
      !
      deallocate(opt_ntri, opt_nvert, opt_vert_coord, opt_vert_index)
      !
      ! write collapsed isosurface on IBZ
      tag_in = "collapsed_IBZ_FS_tri.off"
      call write_IBZ_isosurface(tag_in, num_wann)
      !
    end if

    !--------------- Tangential relaxation
    !
    if(relax) then
      do ibnd = 1, num_wann
        if(nvert(ibnd).eq.0) cycle
        !
        call tangential_relaxation(relax_vinface, eps_vinface, eps_dupv, relax_iter, verbose, &
                                   nrpts, irvec, ndegen, alat, ag, bg, nsym, s, TR_sym, &
                                   num_wann, ham_r, ibnd, nvert(ibnd), ntri(ibnd), &
                                   nfaces_IBZ, faces_IBZ_as_vert, vert_IBZ, &
                                   v_index(:,:,ibnd), v_coord(:,:,ibnd))
        !
      end do
      ! write relaxed isosurface on IBZ
      tag_in = "relaxed_IBZ_FS_tri.off"
      call write_IBZ_isosurface(tag_in, num_wann)
    end if

    !--------------- Newton-Raphson relaxation after improvement
    !
    if(newton_raphson.ge.1) then
      call newton_rap(eps_dupv, eps_vinface, ef, nrpts, irvec, ndegen, ham_r, alat, ag, bg, &
                      nsym, s, TR_sym, num_wann, newton_iter, nfaces_IBZ, faces_IBZ_as_vert, &
                      vert_IBZ, ntri, nvert, v_coord)
      ! write Newton-relaxed isosurface on IBZ
      tag_in = "newton_IBZ_FS_tri.off"
      call write_IBZ_isosurface(tag_in, num_wann)
    end if


  end subroutine mesh_optimization


  subroutine newton_rap(eps_dupv, eps_vinface, ef, nrpts, irvec, ndegen, ham_r, alat, ag, bg, &
                        nsym, s, TR_sym, num_wann, newton_iter, nibz_faces, faces_ibz, &
                        verts_ibz, in_ntri, in_nvert, v_coord)

    use intw_useful_constants, only: tpi
    use intw_matrix_vector, only: ainv, norma, cross
    use triFS_isosurface, only: nvert, calculate_energy_sym, velocity_sym


    implicit none

    ! I/O
    real(dp), intent(in) :: eps_dupv, eps_vinface
    real(dp), intent(in) :: ef
    integer, intent(in) :: nrpts, irvec(1:3,nrpts), ndegen(nrpts) ! Wannier/Fourier space variables
    complex(dp), intent(in) :: ham_r(num_wann,num_wann,nrpts) ! Hamiltonian in Wannier basis
    real(dp), intent(in) :: alat, ag(1:3,1:3), bg(1:3,1:3) ! real and reciprocal lattice vectors
    integer, intent(in) :: nsym, s(1:3,1:3,nsym) ! symmetry variables
    logical, intent(in) :: TR_sym ! time-reversal symmetry
    integer, intent(in) :: newton_iter, nibz_faces, num_wann
    integer, intent(in) :: faces_ibz(:,:) ! vertex indeces of IBZ big tetra faces
    real(dp), intent(in) :: verts_ibz(:,:) ! IBZ big tetra vertex coordinates
    integer, intent(in) :: in_ntri(:), in_nvert(:)
    real(dp), intent(inout) :: v_coord(:,:,:)

    ! Parameters
    !real(dp), parameter :: eps = 1.0E-6_dp
    real(dp), parameter :: lambda = 0.1_dp
    ! Local
    real(dp), allocatable :: vface(:,:,:)
    real(dp) :: eig_int(num_wann), v_w(3), vcoord_crys(3)
    real(dp) :: unit_u(3), unit_u1(3)
    integer :: i_iter, ibnd, iv, jv, i, j, maxvert
    real(dp) :: norma_vw
    logical, allocatable :: vert_on_border(:,:), vert_on_corner(:)
    logical, allocatable :: has_SplusG(:), SplusG_pair(:,:,:,:,:,:,:)
    logical, allocatable :: vertex_related(:), relaxed_vertex(:)


    write(*,*) "Newton-Raphson relaxation of FS..."

    maxvert = maxval(in_nvert(:))
    allocate(vert_on_border(maxvert,nibz_faces), vert_on_corner(maxvert))
    allocate(has_SplusG(maxvert), SplusG_pair(maxvert,maxvert,nsym,-1:1,-1:1,-1:1,1:2))
    allocate(vertex_related(maxvert), relaxed_vertex(maxvert))

    ! Define vectors of IBZ big tetra faces
    allocate(vface(3,1:3,nibz_faces))
    do i = 1, nibz_faces
      vface(1:3,1,i) = verts_ibz(1:3,faces_ibz(2,i)) - verts_ibz(1:3,faces_ibz(1,i))
      vface(1:3,2,i) = verts_ibz(1:3,faces_ibz(3,i)) - verts_ibz(1:3,faces_ibz(1,i))
      vface(1:3,3,i) = verts_ibz(1:3,faces_ibz(1,i))
    end do

    ! Loop over different FS sheets
    do ibnd = 1, num_wann
      !
      if (nvert(ibnd).eq.0) cycle


      !--------------- Check which vertices are on border of IBZ
      !
      call check_border(in_nvert(ibnd), v_coord(:,:,ibnd), nibz_faces, vface, nsym, s, TR_sym, &
                        eps_vinface, vert_on_border, vert_on_corner)


      !--------------- Detect which vertices have SplusG pairs
      !
      has_SplusG = .false.
      vertex_related = .false.
      SplusG_pair = .false.
      do iv = 1, in_nvert(ibnd)
        !
        if (vertex_related(iv)) cycle ! a pair of this vertex has been already assigned
        !
        do jv = 1, in_nvert(ibnd)
          !
          if(jv.eq.iv) cycle
          !
          call vertices_related_by_SplusG(eps_dupv, v_coord(1:3,iv,ibnd), &
                                          v_coord(1:3,jv,ibnd), bg, nsym, s, TR_sym, &
                                          SplusG_pair(iv,jv,:,:,:,:,:))
          !
          if (ANY(SplusG_pair(iv,jv,:,:,:,:,:))) then
            !
            print*, "Vertices iv and jv related:", iv, jv
            has_SplusG(iv) = .true.
            vertex_related(jv) = .true.
            !
          end if
        end do ! jv
      end do ! iv


      !--------------- NR iterations on vertices
      !
      do i_iter = 1, newton_iter
        !
        print*, "Newton-Raphson relaxing, iter=", i_iter
        relaxed_vertex = .false.
        v_loop :do iv = 1, in_nvert(ibnd)
          !
          if (relaxed_vertex(iv)) cycle
          !
          if(ANY(SplusG_pair(iv,iv,:,:,:,:,:))) cycle ! don't relax vertex which has S+G with itself.
          !
          vcoord_crys(1:3) = matmul(ainv(bg),v_coord(1:3,iv,ibnd)) ! Vertex in crystal coordinates
          !
          ! Compute energy and velocity of vertex
          call calculate_energy_sym (nsym, s, TR_sym, num_wann, nrpts, ndegen, irvec, ham_r, vcoord_crys, eig_int)
          call velocity_sym (nrpts, irvec, ndegen, alat, ag, bg, nsym, s, TR_sym, num_wann, ibnd, ham_r, v_w, vcoord_crys)
          !
          ! When vertex lies on face project velocity to face-plane
          j_loop: do j = 1, nibz_faces
            !
            if (vert_on_border(iv,j) .and. .not.vert_on_corner(iv)) then
              unit_u(:) = cross(vface(:,2,j),vface(:,1,j)) ! normal vector of face
              unit_u(:) = unit_u(:) / norma( unit_u )
              v_w(1:3) = v_w(1:3) - dot_product(v_w,unit_u)*unit_u(1:3)
              exit j_loop
            endif
          end do j_loop
          !
          ! When vert lies on corner between two faces project velocity to intersection vector
          if(vert_on_corner(iv)) then
            !
            !cycle v_loop
            !
            f_loop : do j = 1, nibz_faces
              !
              if(vert_on_border(iv,j)) then
                !
                !unit_u1(:) = cross(vface(:,2,j), vface(:,1,j))
                !unit_u1(:) = unit_u1(:) / norma( unit_u1 )
                !
                do i = 1, nibz_faces
                  !
                  if(i.eq.j) cycle
                  if(vert_on_border(iv,i)) then
                    !
                    ! Project velocity to corner vector, i.e. edge of the intersection between the two faces
                    call v2edge(eps_vinface, iv, nibz_faces, i, j, vface, v_w)
                    !
                    exit f_loop
                    !
                  end if ! vert_on_border2
                end do ! face_loop2
              end if ! vert_on_border1
            end do f_loop
          end if ! vert_on_face
          !
          norma_vw = norma(v_w(1:3)) ! norm of velocity of vertex
          v_w = matmul(ainv(bg), v_w) ! Change velocity to crystal coordinates
          !
          vcoord_crys(1:3) = vcoord_crys(1:3) - lambda * v_w(1:3)/norma_vw**2 * (eig_int (ibnd)-ef)/tpi ! move vertex
          v_coord(1:3,iv,ibnd) = matmul(bg, vcoord_crys(1:3)) ! transform moved vertex to cartesian coordinates
          relaxed_vertex(iv) = .true.
          !
          ! Move any symmetry related vertex accordingly
          if (has_SplusG(iv)) then
            !
            do jv = 1, in_nvert(ibnd)
              if(ANY(SplusG_pair(iv,jv,:,:,:,:,:))) then ! jv is pair of iv
                !
                v_coord(1:3,jv,ibnd) = rotate_to_SplusG(v_coord(1:3,iv,ibnd), bg, nsym, s, SplusG_pair(iv,jv,:,:,:,:,:))
                relaxed_vertex(jv) = .true.
                !
              end if
            end do ! jv
          end if ! has_SplusG(iv)
          !
        enddo v_loop
        !
        ! Set all vertices as not relaxed again for next loop
        relaxed_vertex(:) = .false.
        !
      enddo ! i_iter
      print*, "Newton-Raphson relaxation finished"
      !
    enddo ! ibnd
    !

    ! Deallocate variables
    deallocate(vface)
    deallocate(vert_on_border, vert_on_corner)
    deallocate(has_SplusG, SplusG_pair)
    deallocate(vertex_related, relaxed_vertex)


  end subroutine newton_rap


  subroutine check_border(nvert, vert_coord, nibz_faces, vface, nsym, s, TR_sym, epsvert, vert_on_border, vert_on_corner)

    use intw_matrix_vector, only: ainv, norma, cross

    implicit none

    ! I/O
    integer, intent(in) :: nvert, nibz_faces
    real(dp), intent(in) :: vert_coord(:,:), vface(:,:,:)
    integer, intent(in) :: nsym, s(1:3,1:3,nsym)
    logical, intent(in) :: TR_sym
    real(dp), intent(in) :: epsvert
    logical, intent(inout) :: vert_on_border(:,:), vert_on_corner(:)

    ! Local variables
    logical :: inner_face(nibz_faces)
    integer :: j, iv, k
    real(dp) :: coord(1:3), face(1:3,1:3), coef(1:3)
    !! Parameters
    !real(dp), parameter :: epsvert = 1.0E-5_dp ! Now defined from input


    ! Check which tetra faces are inner face
    call inner_faces(nibz_faces, vface, inner_face)

    vert_on_border = .false.
    vert_on_corner = .false.
    do iv = 1, nvert
      do j = 1, nibz_faces
        ! Do not consider inner faces
        if(inner_face(j)) cycle
        ! Write coordinate vector of the vertex in terms of tetra face vectors
        ! 3rd vector is first node, zero unless the face is BZ border. Later changed to normal vector of face.
        coord(1:3) = vert_coord(1:3,iv)-vface(1:3,3,j)
        face(:,1) = vface(:,1,j)
        face(:,2) = vface(:,2,j)
        face(:,3) = cross(vface(:,2,j), vface(:,1,j))
        face(:,3) = face(:,3) / norma(face(:,3))
        coef(1:3) = matmul(ainv(face(1:3,1:3)), coord(1:3))
        ! Coordinate vector of the vertex is perpendicular to normal vector of face
        if(abs(coef(3)).lt.epsvert) then
          ! Coordinate vector of the vertex lies in the face plane
          if(coef(1).ge.0.0_dp-epsvert .and. coef(1).le.1.0_dp+epsvert &
            .and. coef(2).ge.0.0_dp-epsvert .and. coef(2).le.1.0_dp+epsvert & !) then
            .and. sum(coef(1:2))-1.0_dp.lt.epsvert ) then
              !print*, "Vertex", iv, "is on tetra face", j
              !print*, "vert iv on face j:", iv, j
              vert_on_border(iv,j) = .true.
          end if
        end if
      end do
      k = 0
      do j = 1, nibz_faces
        if(vert_on_border(iv,j)) k = k+1
      end do
      if(k.ge.2) then
        print*, iv, "vert on corner!"
        vert_on_corner(iv) = .true.
      end if
    end do

  end subroutine check_border


  subroutine vertices_related_by_SplusG(epsvert, ivertex, jvertex, bg, nsym, s, TR_sym, relation)

    use intw_matrix_vector, only: ainv, norma

    implicit none

    ! I/O
    real(dp), intent(in) :: epsvert
    real(dp), intent(in) :: ivertex(1:3), jvertex(1:3), bg(1:3,1:3)
    integer, intent(in) :: nsym, s(1:3,1:3,nsym)
    logical, intent(in) :: TR_sym
    logical, intent(out) :: relation(nsym,-1:1,-1:1,-1:1,1:2) ! true if ivertex and jvertex are related

    ! Local
    integer :: isym, ig, jg, kg
    real(dp) :: ivert_crys(1:3), rot_ivertex(1:3)
    !real(dp), parameter :: epsvert = 1.0E-6_dp


    ! Initialize
    relation = .false.

    ! ivertex in crystal coords
    ivert_crys = matmul(ainv(bg), ivertex)

    do isym = 1, nsym
      do ig = -1, 1
        do jg = -1, 1
          do kg = -1, 1
            !
            if(ig.eq.0 .and. jg.eq.0 .and. kg.eq.0) cycle
            !
            rot_ivertex(1:3) = matmul(bg, &
                                        matmul(dble(s(:,:,isym)), ivert_crys ) &
                                        + real((/ig, jg, kg/), dp) )
            !
            ! check with jvertex
            if (norma(rot_ivertex(:)-jvertex(:)) .lt. epsvert ) then
              relation(isym,ig,jg,kg,1) = .true.
              return
            end if
            ! TR symmetry
            if (TR_sym) then
              !
              rot_ivertex(1:3) = matmul(bg, &
                                        -matmul(dble(s(:,:,isym)), ivert_crys ) &
                                        + real((/ig, jg, kg/), dp) )
              !
              ! check with jvertex
              if (norma(rot_ivertex(:)-jvertex(:)) .lt. epsvert ) then
                relation(isym,ig,jg,kg,2) = .true.
                return
              end if
              !
            end if ! TR_sym
            !
          end do ! kg
        end do ! jg
      end do ! ig
    end do ! isym

  end subroutine vertices_related_by_SplusG


  function rotate_to_SplusG(ivertex, bg, nsym, s, which_SplusG)

    use intw_matrix_vector, only: ainv

    implicit none

    ! I/O
    real(dp), intent(in) :: ivertex(1:3), bg(1:3,1:3)
    integer, intent(in) :: nsym, s(1:3,1:3,nsym)
    logical, intent(in) :: which_SplusG(nsym,-1:1,-1:1,-1:1,1:2)

    ! Local
    real(dp), dimension(1:3) :: rotate_to_SplusG
    integer :: isym, ig, jg, kg
    real(dp) :: ivert_crys(1:3)


    ! ivertex in crystal coords
    ivert_crys = matmul(ainv(bg), ivertex)

    do isym = 1, nsym
      do ig = -1, 1
        do jg = -1, 1
          do kg = -1, 1
            !
            if(ig.eq.0 .and. jg.eq.0 .and. kg.eq.0) cycle
            !
            if (which_SplusG(isym,ig,jg,kg,1)) then
              !
              rotate_to_SplusG(1:3) = matmul(bg, &
                                             matmul(dble(s(:,:,isym)), ivert_crys ) &
                                             + real((/ig, jg, kg/), dp) )
              !
            else if (which_SplusG(isym,ig,jg,kg,2)) then
              !
              rotate_to_SplusG(1:3) = matmul(bg, &
                                             -matmul(dble(s(:,:,isym)), ivert_crys ) &
                                             + real((/ig, jg, kg/), dp) )
              !
            end if
            !
          end do ! kg
        end do ! jg
      end do ! ig
    end do ! isym

  end function rotate_to_SplusG


  subroutine v2edge(epsvert, iv, nfaces, iface, jface, vface, v_w)

    use intw_matrix_vector, only: norma, cross

    implicit none

    ! I/O
    real(dp), intent(in) :: epsvert
    integer, intent(in) :: iv, nfaces, iface, jface
    real(dp), dimension(1:3,1:3,nfaces), intent(in) :: vface
    real(dp), dimension(1:3), intent(inout) :: v_w

    !real(dp), parameter :: epsvert = 1.0E-5_dp

    if( sum(abs(cross(vface(:,1,jface), vface(:,1,iface)))) .lt. epsvert ) then
      if(dot_product( v_w(1:3), vface(1:3,1,iface) ) .gt. epsvert ) then
        v_w(1:3) = dot_product( v_w(1:3), vface(1:3,1,iface) ) * vface(1:3,1,iface) / norma(vface(1:3,1,iface))
      else
        v_w(1:3) = dot_product( v_w(1:3), -vface(1:3,1,iface) ) * (-vface(1:3,1,iface)) / norma(-vface(1:3,1,iface))
      end if
    else if( sum(abs(cross(vface(:,1,jface), vface(:,2,iface)))) .lt. epsvert ) then
      if(dot_product( v_w(1:3), vface(1:3,2,iface) ) .gt. 1.0E-6_dp ) then
        v_w(1:3) = dot_product( v_w(1:3), vface(1:3,2,iface) ) * vface(1:3,2,iface) / norma(vface(1:3,2,iface))
      else
        v_w(1:3) = dot_product( v_w(1:3), -vface(1:3,2,iface) ) * (-vface(1:3,2,iface)) / norma(-vface(1:3,2,iface))
      end if
    else if( sum(abs(cross(vface(:,2,jface),vface(:,1,iface)))) .lt. epsvert ) then
      if(dot_product( v_w(1:3), vface(1:3,1,iface) ) .gt. 1.0E-6_dp ) then
        v_w(1:3) = dot_product( v_w(1:3), vface(1:3,1,iface) ) * vface(1:3,1,iface) / norma(vface(1:3,1,iface))
      else
        v_w(1:3) = dot_product( v_w(1:3), -vface(1:3,1,iface) ) * (-vface(1:3,1,iface)) / norma(-vface(1:3,1,iface))
      end if
    else if( sum(abs(cross(vface(:,2,jface), vface(:,2,iface)))) .lt. epsvert ) then
      if(dot_product( v_w(1:3), vface(1:3,2,iface) ) .gt. 1.0E-6_dp ) then
        v_w(1:3) = dot_product( v_w(1:3), vface(1:3,2,iface) ) * vface(1:3,2,iface) / norma(vface(1:3,2,iface))
      else
        v_w(1:3) = dot_product( v_w(1:3), -vface(1:3,2,iface) ) * (-vface(1:3,2,iface)) / norma(-vface(1:3,2,iface))
      end if
    else if( sum(abs(cross(vface(:,2,jface)-vface(:,1,jface), vface(:,1,iface)))) .lt. epsvert ) then
      if(dot_product( v_w(1:3), vface(1:3,1,iface) ) .gt. 1.0E-6_dp ) then
        v_w(1:3) = dot_product( v_w(1:3), vface(1:3,1,iface) ) * vface(1:3,1,iface) / norma(vface(1:3,1,iface))
      else
        v_w(1:3) = dot_product( v_w(1:3), -vface(1:3,1,iface) ) * (-vface(1:3,1,iface)) / norma(-vface(1:3,1,iface))
      end if
    else if( sum(abs(cross(vface(:,2,jface)-vface(:,1,jface), vface(:,2,iface)))) .lt. epsvert ) then
      if(dot_product( v_w(1:3), vface(1:3,2,iface) ) .gt. 1.0E-6_dp ) then
        v_w(1:3) = dot_product( v_w(1:3), vface(1:3,2,iface) ) * vface(1:3,2,iface) / norma(vface(1:3,2,iface))
      else
        v_w(1:3) = dot_product( v_w(1:3), -vface(1:3,2,iface) ) * (-vface(1:3,2,iface)) / norma(-vface(1:3,2,iface))
      end if
    else if( sum(abs(cross(vface(:,1,jface), vface(:,2,iface)-vface(:,1,iface)))) .lt. epsvert ) then
      if(dot_product( v_w(1:3), vface(1:3,1,iface) ) .gt. 1.0E-6_dp ) then
        v_w(1:3) = dot_product( v_w(1:3), vface(1:3,1,jface) ) * vface(1:3,1,jface) / norma(vface(1:3,1,jface))
      else
        v_w(1:3) = dot_product( v_w(1:3), -vface(1:3,1,jface) ) * (-vface(1:3,1,jface)) / norma(-vface(1:3,1,jface))
      end if
    else if( sum(abs(cross(vface(:,2,jface), vface(:,2,iface)-vface(:,1,iface)))) .lt. epsvert ) then
      if(dot_product( v_w(1:3), vface(1:3,1,iface) ) .gt. 1.0E-6_dp ) then
        v_w(1:3) = dot_product( v_w(1:3), vface(1:3,2,jface) ) * vface(1:3,2,jface) / norma(vface(1:3,2,jface))
      else
        v_w(1:3) = dot_product( v_w(1:3), -vface(1:3,2,jface) ) * (-vface(1:3,2,jface)) / norma(-vface(1:3,2,jface))
      end if
    else if( sum(abs(cross(vface(:,2,jface)-vface(:,2,jface), vface(:,2,iface)-vface(:,1,iface)))) .lt. epsvert ) then
      if ( dot_product( v_w(1:3), (vface(:,2,iface)-vface(:,1,iface)) ) .gt. 1.0E-6_dp ) then
        v_w(1:3) = dot_product( v_w(1:3), (vface(:,2,iface)-vface(:,1,iface)) ) * (vface(:,2,iface)-vface(:,1,iface)) / norma((vface(:,2,iface)-vface(:,1,iface)))
      else
        v_w(1:3) = dot_product( v_w(1:3), (-vface(:,2,iface)+vface(:,1,iface)) ) * (-vface(:,2,iface)+vface(:,1,iface)) / norma((-vface(:,2,iface)+vface(:,1,iface)))
      end if
    else
      print*, "Something went wrong!"
      print*, iv, iface, jface
      stop
    end if

  end subroutine v2edge


  subroutine edge_collapse(eps, in_nvert, in_ntri, nibz_faces, faces_ibz, verts_ibz, dif_vert_coord, dif_vert_index, new_nvert, new_ntri, new_dif_vert_coord, new_dif_vert_index)

    use intw_matrix_vector, only: ainv, norma, cross

    implicit none

    ! I/O variables
    real(dp), intent(in) :: eps
    integer, intent(in) :: in_nvert, in_ntri ! input number of triangles, vertices
    integer, intent(in) :: nibz_faces ! number of faces on IBZ
    integer, intent(in) :: faces_ibz(:,:) ! vertex indeces of IBZ big tetra faces
    real(dp), intent(in) :: verts_ibz(:,:) ! IBZ big tetra vertex coordinates
    real(dp), dimension(:,:), intent(inout) :: dif_vert_coord ! input vertex coords which will be changed
    integer, dimension(:,:), intent(inout) :: dif_vert_index ! input triangle indices which will be changed
    integer, intent(out) :: new_nvert, new_ntri ! output number of vertices and triangles after simplification
    real(dp), dimension(:,:), intent(out) :: new_dif_vert_coord ! output vertices coordinates
    integer, dimension(:,:), intent(out) :: new_dif_vert_index ! output triangle vertices

    ! Local variables
    logical, allocatable :: vert_on_border(:,:), vert_on_corner(:), inner_face(:)
    integer, allocatable :: neighbour_triangles(:,:,:)
    real(dp), dimension(3) :: coef,coord
    real(dp), dimension(1:3,1:3) :: face
    real(dp), dimension(:,:,:), allocatable :: vface
    integer :: i, iv, j, ij, it, ntri_aux, nvert_aux, iter, f, k
    integer :: rmvd_vrtx, rmvd_tri1, rmvd_tri2
    ! Parameters
    real(dp), parameter :: epsvert = 1.0E-6_dp
    integer, parameter :: nghmax = 500, tri_elim_iter = 1000, dummy_tri_iter = 1000


    !--------------- Define vectors of IBZ big tetra faces
    !
    allocate(vface(3,3,nibz_faces))
    call faces_vectors(nibz_faces, faces_ibz, verts_ibz, vface)

    !--------------- Check which tetra faces are inner
    !
    allocate(inner_face(nibz_faces))
    call inner_faces(nibz_faces, vface, inner_face)

    !--------------- Collapse edges
    ! I cannot use check_border subr as vertex list is changing constantly...
    !
    allocate(vert_on_border(in_nvert,nibz_faces), vert_on_corner(in_nvert))
    allocate(neighbour_triangles(in_ntri,3,nghmax))
    new_dif_vert_coord = 0.0_dp
    new_dif_vert_index = 0
    new_ntri = in_ntri
    new_nvert = in_nvert
    do iter = 1, tri_elim_iter

      print*, ""
      print*, "Triangle elimination:", iter
      ntri_aux = new_ntri
      nvert_aux = new_nvert

      print*, "Old number of vertices:", nvert_aux
      print*, "Old number of triangles:", ntri_aux


      !--------------- Remove vertices with three neighbours
      !
      ! Check which vertex should be removed
      print*, ""
      do it = 1, dummy_tri_iter
        print*, "Removing dummy triangles iteration", it

        ! Assign neighbours
        neighbour_triangles = 0
        call find_neighbours(ntri_aux, dif_vert_index, neighbour_triangles, nghmax)

        ! Check which vertices are on tetra faces
        vert_on_border = .false.
        vert_on_corner = .false.
        do iv = 1, nvert_aux
          do j = 1, nibz_faces

            !!! JL remove not considering inner faces. I think it makes the difference when S+G symmetries
            ! Do not consider inner faces
            if(inner_face(j)) cycle
            !!! JL

            ! Write coordinate vector of the vertex in terms of tetra face vectors
            ! 3rd vector is first node, zero unless the face is BZ border. Later changed to normal vector of face.
            coord(1:3) = dif_vert_coord(1:3,iv)-vface(1:3,3,j)
            face(:,1) = vface(:,1,j)
            face(:,2) = vface(:,2,j)
            face(:,3) = cross(vface(:,2,j),vface(:,1,j))
            face(:,3) = face(:,3) / norma(face(:,3))
            coef(1:3) = matmul(ainv(face(1:3,1:3)), coord(1:3))
            ! Coordinate vector of the vertex is perpendicular to normal vector of face
            !if(abs(dot_product(face(:,3),coord(:))).lt.epsvert) then
            if(abs(coef(3)).lt.epsvert) then
              ! Coordinate vector of the vertex lies in the face plane
              if(coef(1).ge.0.0_dp-epsvert .and. coef(1).le.1.0_dp+epsvert &
                 .and. coef(2).ge.0.0_dp-epsvert .and. coef(2).le.1.0_dp+epsvert & !) then
                 .and. sum(coef(1:2))-1.0_dp.lt.epsvert ) then
                  !print*, "vert iv on face j:", iv, j
                  vert_on_border(iv,j) = .true.
              end if
            end if
          end do
          k = 0
          do j = 1, nibz_faces
            if(vert_on_border(iv,j)) k = k+1
          end do
          if(k.ge.2) then
            vert_on_corner(iv) = .true.
          end if
        end do

        rmvd_vrtx = 0
        rmvd_tri1 = 0
        rmvd_tri2 = 0
        f = 0
        call rm_3tri_vrtx(ntri_aux, dif_vert_index, neighbour_triangles, nghmax, nibz_faces, &
                          vert_on_border, rmvd_vrtx, rmvd_tri1, rmvd_tri2, f)

        if(f==0) exit
        ! Define new vertex coordinates and lists
        new_nvert = 0
        do iv = 1, nvert_aux
          if(iv==rmvd_vrtx) cycle
          new_nvert = new_nvert+1
          dif_vert_coord(:,new_nvert) = dif_vert_coord(:,iv)
        end do
        new_ntri = 0
        do i = 1, ntri_aux
          if(i==rmvd_tri1) cycle
          if(i==rmvd_tri2) cycle
          new_ntri = new_ntri+1
          dif_vert_index(:,new_ntri) = dif_vert_index(:,i)
          do ij = 1, 3
            if(dif_vert_index(ij,new_ntri).gt.rmvd_vrtx) dif_vert_index(ij,new_ntri) = dif_vert_index(ij,new_ntri)-1
          end do
        end do
        nvert_aux = new_nvert
        ntri_aux = new_ntri
      end do
      print*, "Number of vertices:", nvert_aux
      print*, "Number of triangles:", ntri_aux


      !--------------- Collapse vertices
      !
      ! Assign neighbours
      neighbour_triangles = 0
      call find_neighbours(ntri_aux, dif_vert_index, neighbour_triangles, nghmax)

      ! Collapse triangles and create new mesh
      call collapse_triangles(eps, ntri_aux, nvert_aux, dif_vert_index, dif_vert_coord, &
                              vert_on_border, nibz_faces, nghmax, neighbour_triangles, &
                              new_ntri, new_nvert, new_dif_vert_coord, new_dif_vert_index)


      !--------------- Re-assign vertex and index coordinates
      !
      do iv = 1, new_nvert
        dif_vert_coord(:,iv) = new_dif_vert_coord(:,iv)
      end do
      dif_vert_index(:,:) = new_dif_vert_index(:,:)

      print*, "Old number of triangles", ntri_aux
      print*, "New number of triangles", new_ntri
      print*, "Old number of vertices", nvert_aux
      print*, "New number of vertices", new_nvert

      if(ntri_aux==new_ntri) then
        print*, "No triangle obeys collapse criteria. Mesh simplification is finished."
        exit
      end if

    end do

    !--------------- Deallocate variables
    !
    deallocate(vface, inner_face, neighbour_triangles, vert_on_border, vert_on_corner)


  end subroutine edge_collapse


  subroutine faces_vectors(nibz_faces, faces_ibz, verts_ibz, vface)
    ! Reads four corners of tetra and outputst two vectors of each face,
    ! with the origin at Gamma (always set as the first tetra node by convention)

    implicit none

    ! I/O
    integer, intent(in) :: nibz_faces
    integer, intent(in) :: faces_ibz(:,:) ! vertex indeces of IBZ big tetra faces
    real(dp), intent(in) :: verts_ibz(:,:) ! IBZ big tetra vertex coordinates
    real(dp), intent(out) :: vface(1:3,1:3,nibz_faces)

    ! Local
    integer :: i


    do i = 1, nibz_faces
      vface(1:3,1,i) = verts_ibz(1:3,faces_ibz(2,i)) - verts_ibz(1:3,faces_ibz(1,i))
      vface(1:3,2,i) = verts_ibz(1:3,faces_ibz(3,i)) - verts_ibz(1:3,faces_ibz(1,i))
      vface(1:3,3,i) = verts_ibz(1:3,faces_ibz(1,i))
    end do

  end subroutine faces_vectors


  subroutine inner_faces(nibz_faces, vface, inner_face)
    ! Checks which tetra faces are inner

    use intw_matrix_vector, only: norma

    implicit none

    ! I/O
    integer, intent(in) :: nibz_faces
    real(dp), intent(in) :: vface(:,:,:)
    logical, intent(out) :: inner_face(nibz_faces)

    ! Local
    integer :: i, j
    ! Parameter
    real(dp), parameter :: epsvert = 1.0E-5_dp


    inner_face = .false.
    do i = 1, nibz_faces
      do j = 1, nibz_faces
        if(j.eq.i) cycle
        if((norma(vface(:,1,i)-vface(:,1,j)).lt.epsvert        &
            .and. norma(vface(:,2,i)-vface(:,2,j)).lt.epsvert  &
            .and. norma(vface(:,3,i)).lt.epsvert               &
            .and. norma(vface(:,3,j)).lt.epsvert)              &
        .or. (norma(vface(:,1,i)-vface(:,2,j)).lt.epsvert      &
            .and. norma(vface(:,2,i)-vface(:,1,j)).lt.epsvert  &
            .and. norma(vface(:,3,i)).lt.epsvert               &
            .and. norma(vface(:,3,j)).lt.epsvert)              ) then
          print*, "Face", i, "is inner face"
          inner_face(i) = .true.
          exit
        end if
      end do
    end do

  end subroutine inner_faces


  subroutine find_neighbours(in_ntri, vindex, neighbour_triangles, nghmax)

    implicit none

    ! I/O
    integer, intent(in) :: in_ntri
    integer, dimension(:,:), intent(in) :: vindex
    integer, dimension(:,:,:), intent(inout) :: neighbour_triangles
    integer, intent(in) :: nghmax

    ! Local
    integer :: i, iv, j, jv, ngh


    ! Loop over vertices
    do i = 1, in_ntri
      do iv = 1, 3
        ! Find neighbour triangles of each vertex
        ngh = 0
        do j = 1, in_ntri
          if(j.eq.i) cycle ! Skip the triangle itself
          do jv = 1, 3
            if(ngh.ge.nghmax) then
              print*, "ERROR:Maximum number of neighbours reached"
              stop
            else if(vindex(iv,i).eq.vindex(jv,j)) then
              ngh = ngh+1
              neighbour_triangles(i,iv,ngh) = j ! j is a neighbour triangle of vertex (i,iv)
            end if
          end do
        end do
      end do
    end do

  end subroutine find_neighbours


  subroutine find_neighbours_vertex_list(ntri, nvert, vindex, neighbour_triangles, nghmax)

    implicit none

    integer, intent(in) :: ntri, nvert
    integer, dimension(:,:), intent(in) :: vindex
    integer, dimension(:,:), intent(inout) :: neighbour_triangles
    integer, intent(in) :: nghmax

    !Local variables
    integer :: iv, j, jv, ngh


    !Loop over vertices
    do iv = 1, nvert
      ! Find neighbour triangles of each vertex
      ngh = 0
      do j = 1, ntri
        do jv = 1, 3
          if(ngh.ge.nghmax) then
            print*, "ERROR:Maximum number of neighbours reached"
            stop
          else if(iv.eq.vindex(jv,j)) then
            ngh = ngh+1
            neighbour_triangles(iv,ngh) = j ! j is a neighbour triangle of vertex iv
          end if
        end do
      end do
    end do

  end subroutine find_neighbours_vertex_list


  subroutine rm_3tri_vrtx(ntri, vindex, neighbours, nghmax, faces, border_vert, rmvd_vrtx, rmvd_tri1, rmvd_tri2, n)
    ! Removes problematic vertices which are shared by three triangles

    implicit none

    ! I/O
    integer, intent(in) :: ntri
    integer, dimension(:,:), intent(inout) :: vindex
    integer, dimension(:,:,:), intent(in) :: neighbours
    integer, intent(in) :: nghmax
    integer, intent(in) :: faces
    logical, dimension(:,:), intent(in) :: border_vert
    integer, intent(out) :: rmvd_vrtx, rmvd_tri1, rmvd_tri2, n

    ! Local
    integer :: i, j, iv, jv, k, f


    n = 0
    i_loop: do i = 1, ntri
      iv_loop: do iv = 1, 3
        k = 0
        do j = 1, nghmax
          if(neighbours(i,iv,j)==0) exit
          k = k+1
        end do ! j
        if(k==2) then ! Two neighbours because the triangle itself is skipped in find_neighbours
          do f = 1, faces
            if(border_vert(vindex(iv,i),f)) then
              cycle iv_loop ! If vertex with three neigh is on border ignore
            end if
          end do
          ! Re-assign position of removed vertex. The new vertex position is the vertex of the neighbour triangle which is not shared with the modified triangle
          do jv = 1, 3
            if(vindex(jv,neighbours(i,iv,1)).ne.vindex(1,i) .and. vindex(jv,neighbours(i,iv,1)).ne.vindex(2,i) .and. vindex(jv,neighbours(i,iv,1)).ne.vindex(3,i)) then
              rmvd_vrtx = vindex(iv,i)
              vindex(iv,i) = vindex(jv,neighbours(i,iv,1))
              exit
            end if
          end do ! jv
          rmvd_tri1 = neighbours(i,iv,1)
          rmvd_tri2 = neighbours(i,iv,2)
          n = n+1
          print*, "Vertex", vindex(iv,i), "has three neighbours"
          exit i_loop
        end if
      end do iv_loop
    end do i_loop

    print*, "Number of triangles with three vertex", n

  end subroutine rm_3tri_vrtx


  subroutine collapse_triangles(eps, ntri_old, nvert_old, vindex, vcoord, vert_border, ntetraf, nghmax, neighbour_triangles, &
                                ntri_new, nvert_new, new_vcoord, new_vindex)
    ! Subroutine that performs edge collapse and triangle re-assignation.
    ! Needs clean up...

    use intw_matrix_vector, only: norma, cross

    implicit none

    ! I/O
    real(dp), intent(in) :: eps
    integer, intent(in) :: ntri_old, nvert_old
    integer, dimension(:,:), intent(inout) :: vindex
    real(dp), dimension(:,:), intent(inout) :: vcoord
    logical, dimension(:,:), intent(in) :: vert_border
    integer, intent(in) :: ntetraf
    integer, intent(in) :: nghmax
    integer, dimension(:,:,:), intent(in) :: neighbour_triangles
    integer, intent(out) :: ntri_new, nvert_new
    real(dp), dimension(:,:), intent(inout) :: new_vcoord
    integer, dimension(:,:), intent(inout) :: new_vindex

    !Local variables
    integer :: i, iv, ij, m, mj, j, jv, f, k1, k2, collapsed_edge, eliminated_vertex
    logical, dimension(ntri_old) :: triangle_elimination
    real(dp), dimension(3) :: oldv1, oldv2, oldnormal, newv1, newv2, newnormal
    real(dp), dimension(3,nvert_old) :: old_vcoord, aux_vcoord


    ntri_new = ntri_old
    eliminated_vertex = huge(eliminated_vertex)
    triangle_elimination = .false.
    old_vcoord = vcoord
    aux_vcoord = vcoord
    !
    tri_loop: do i = 1, ntri_old ! Loop over triangles

      ! Check if triangle fulfils collapse condition
      call collapse_condition(vindex(1,i), vindex(2,i), vindex(3,i), vcoord(:,:), eps, collapsed_edge)


      !--------------- Go over all possibilities of collapse...
      !

      !--------------- Edge 1-2 of the triangle is small
      !
      if ( collapsed_edge == 1 ) then

        ! no vertex is on tetra face, move
        if( (ALL(.not.vert_border(vindex(1,i),:)) &
            .and. ALL(.not.vert_border(vindex(2,i),:))) &
            !! or both vertices are on same tetra face
            !.or. ANY(vert_border(vindex(1,i),:) .and. vert_border(vindex(2,i),:)) &
              ) then

          ! loop to check if any vertex is on two tetra-faces, i.e. corner of IBZ
          k1 = 0
          k2 = 0
          do f = 1, ntetraf
            if(vert_border(vindex(1,i),f)) k1 = k1+1
            if(vert_border(vindex(2,i),f)) k2 = k2+1
            if(k1.ge.2 .or. k2.ge.2) cycle tri_loop
          end do

          print*, "Triangle",i,"collapsed from edge",12

          aux_vcoord(:,vindex(2,i)) = ( old_vcoord(:,vindex(1,i)) + old_vcoord(:,vindex(2,i)) ) / 2.d0
          aux_vcoord(:,vindex(1,i)) = ( old_vcoord(:,vindex(1,i)) + old_vcoord(:,vindex(2,i)) ) / 2.d0

          ! Check if collapse would lead to flip of any neighbour triangle of vertex1
          do ij = 1, nghmax

            if(neighbour_triangles(i,1,ij)==0) exit

            oldv1(:) = old_vcoord(:,vindex(2,neighbour_triangles(i,1,ij)))-old_vcoord(:,vindex(1,neighbour_triangles(i,1,ij)))
            oldv2(:) = old_vcoord(:,vindex(3,neighbour_triangles(i,1,ij)))-old_vcoord(:,vindex(1,neighbour_triangles(i,1,ij)))
            oldnormal(:) = cross(oldv2(:),oldv1(:))
            oldnormal(:) = oldnormal(:)/norma(oldnormal)

            newv1(:) = aux_vcoord(:,vindex(2,neighbour_triangles(i,1,ij)))-aux_vcoord(:,vindex(1,neighbour_triangles(i,1,ij)))
            newv2(:) = aux_vcoord(:,vindex(3,neighbour_triangles(i,1,ij)))-aux_vcoord(:,vindex(1,neighbour_triangles(i,1,ij)))
            newnormal(:) = cross(newv2(:),newv1(:))
            newnormal(:) = newnormal(:)/norma(newnormal)

            if(dot_product(oldnormal, newnormal).lt.0.0_dp) then
              print*, "Triangle will not collapse in order to avoid flip of a neighbour triangle"
              cycle tri_loop
            end if

          end do

          ! Check if collapse would lead to flip of any neighbour triangle of vertex2
          do ij = 1, nghmax

            if(neighbour_triangles(i,2,ij)==0) exit

            oldv1(:) = old_vcoord(:,vindex(2,neighbour_triangles(i,2,ij)))-old_vcoord(:,vindex(1,neighbour_triangles(i,2,ij)))
            oldv2(:) = old_vcoord(:,vindex(3,neighbour_triangles(i,2,ij)))-old_vcoord(:,vindex(1,neighbour_triangles(i,2,ij)))
            oldnormal(:) = cross(oldv2(:),oldv1(:))
            oldnormal(:) = oldnormal(:)/norma(oldnormal)

            newv1(:) = aux_vcoord(:,vindex(2,neighbour_triangles(i,2,ij)))-aux_vcoord(:,vindex(1,neighbour_triangles(i,2,ij)))
            newv2(:) = aux_vcoord(:,vindex(3,neighbour_triangles(i,2,ij)))-aux_vcoord(:,vindex(1,neighbour_triangles(i,2,ij)))
            newnormal(:) = cross(newv2(:),newv1(:))
            newnormal(:) = newnormal(:)/norma(newnormal)

            if(dot_product(oldnormal, newnormal).lt.0.0_dp) then
              print*, "Triangle will not collapse in order to avoid flip of a neighbour triangle"
              cycle tri_loop
            end if

          end do

          ! Collapse triangle sharing the removed edge (opposite triangle)
          do ij = 1, nghmax
          j = neighbour_triangles(i,1,ij)
            if(j==0) exit
            do mj = 1, nghmax
              m = neighbour_triangles(i,2,mj)
              if(m==0) exit
              ! Triangle m is neighbour of the two collapsed vertices, i.e. is the opposite triangle
              if(m==j) then
                if(     ANY(vert_border(vindex(1,m),:)) &
                   .or. ANY(vert_border(vindex(2,m),:)) &
                   .or. ANY(vert_border(vindex(3,m),:)) &
                   ) then
                  print*, "Triangle will not be collapsed as its opposite is on a tetra face"
                  cycle tri_loop
                end if
                print*, "Triangle", j, "is opposite triangle of", i, "and will be collapsed"
                triangle_elimination(m) = .true.
                ntri_new = ntri_new - 1
              end if
            end do
          end do
          ! Re-assign coordinate of vertex1
          vcoord(:,vindex(1,i)) = ( old_vcoord(:,vindex(1,i)) + old_vcoord(:,vindex(2,i)) ) / 2.d0
          ! Collapse triangle
          triangle_elimination(i) = .true.
          ntri_new = ntri_new - 1
          ! Remove vertex2, become vertex1
          eliminated_vertex = vindex(2,i)
          ! Re-define index of vertices of neighbour triangles of collapsed vertex2
          do ij = 1, nghmax
            j = neighbour_triangles(i,2,ij)
            if(j==0) exit
            if(.not.triangle_elimination(j)) then
                do jv = 1, 3
                  if( vindex(2,i).eq.vindex(jv,j) ) then
                    print*, "Re-locating vertex", jv, "of triangle", j, "neighbour of", i, 2
                    vindex(jv,j) = vindex(1,i)
                  end if
                end do
            end if
          end do


          ! -- Exit just to remove one triangle for the moment
          exit
          ! --

        else
          cycle
        end if

      !--------------- Edge 1-3 of the triangle is small
      !
      else if ( collapsed_edge == 2 ) then

        ! no vertex is on tetra face, forward
        if( (ALL(.not.vert_border(vindex(1,i),:)) &
            .and. ALL(.not.vert_border(vindex(3,i),:))) &
            !! or both vertices are on same tetra face
            !.or. ANY(vert_border(vindex(1,i),:) .and. vert_border(vindex(3,i),:)) &
              ) then

          ! loop to check if any vertex is on two tetra-faces, i.e. corner of irrBZ
          k1 = 0
          k2 = 0
          do f = 1, ntetraf
            if(vert_border(vindex(1,i),f)) k1 = k1+1
            if(vert_border(vindex(3,i),f)) k2 = k2+1
            if(k1.ge.2 .or. k2.ge.2) cycle tri_loop
          end do


          print*, "Triangle", i, "collapsed from edge", 13

          aux_vcoord(:,vindex(3,i)) = ( old_vcoord(:,vindex(1,i)) + old_vcoord(:,vindex(3,i)) ) / 2.d0
          aux_vcoord(:,vindex(1,i)) = ( old_vcoord(:,vindex(1,i)) + old_vcoord(:,vindex(3,i)) ) / 2.d0

          ! Check if collapse would lead to flip of any neighbour triangle of vertex1
          do ij = 1, nghmax

            if(neighbour_triangles(i,1,ij)==0) exit

            oldv1(:) = old_vcoord(:,vindex(2,neighbour_triangles(i,1,ij)))-old_vcoord(:,vindex(1,neighbour_triangles(i,1,ij)))
            oldv2(:) = old_vcoord(:,vindex(3,neighbour_triangles(i,1,ij)))-old_vcoord(:,vindex(1,neighbour_triangles(i,1,ij)))
            oldnormal(:) = cross(oldv2(:),oldv1(:))
            oldnormal(:) = oldnormal(:)/norma(oldnormal)

            newv1(:) = aux_vcoord(:,vindex(2,neighbour_triangles(i,1,ij)))-aux_vcoord(:,vindex(1,neighbour_triangles(i,1,ij)))
            newv2(:) = aux_vcoord(:,vindex(3,neighbour_triangles(i,1,ij)))-aux_vcoord(:,vindex(1,neighbour_triangles(i,1,ij)))
            newnormal(:) = cross(newv2(:),newv1(:))
            newnormal(:) = newnormal(:)/norma(newnormal)

            if(dot_product(oldnormal, newnormal).lt.0.0_dp) then
              print*, "Triangle will not collapse in order to avoid flip of a neighbour triangle"
              cycle tri_loop
            end if

          end do

          ! Check if collapse would lead to flip of any neighbour triangle of vertex3
          do ij = 1, nghmax

            if(neighbour_triangles(i,3,ij)==0) exit

            oldv1(:) = old_vcoord(:,vindex(2,neighbour_triangles(i,3,ij)))-old_vcoord(:,vindex(1,neighbour_triangles(i,3,ij)))
            oldv2(:) = old_vcoord(:,vindex(3,neighbour_triangles(i,3,ij)))-old_vcoord(:,vindex(1,neighbour_triangles(i,3,ij)))
            oldnormal(:) = cross(oldv2(:),oldv1(:))
            oldnormal(:) = oldnormal(:)/norma(oldnormal)

            newv1(:) = aux_vcoord(:,vindex(2,neighbour_triangles(i,3,ij)))-aux_vcoord(:,vindex(1,neighbour_triangles(i,3,ij)))
            newv2(:) = aux_vcoord(:,vindex(3,neighbour_triangles(i,3,ij)))-aux_vcoord(:,vindex(1,neighbour_triangles(i,3,ij)))
            newnormal(:) = cross(newv2(:),newv1(:))
            newnormal(:) = newnormal(:)/norma(newnormal)

            if(dot_product(oldnormal, newnormal).lt.0.0_dp) then
              print*, "Triangle will not collapse in order to avoid flip of a neighbour triangle"
              print*, dot_product(oldnormal, newnormal)
              cycle tri_loop
            end if

          end do

          ! Collapse triangle sharing the removed edge (opposite triangle)
          do ij = 1, nghmax
            j = neighbour_triangles(i,1,ij)
            if(j==0) exit
            do mj = 1, nghmax
              m = neighbour_triangles(i,3,mj)
              if(m==0) exit
              ! Triangle m is neighbour of the two collapsed vertices, i.e. is the opposite triangle
              if(m==j) then
                if(     ANY(vert_border(vindex(1,m),:)) &
                   .or. ANY(vert_border(vindex(2,m),:)) &
                   .or. ANY(vert_border(vindex(3,m),:)) &
                   ) then
                  print*, "Triangle will not be collapsed as its opposite is on a tetra face"
                  cycle tri_loop
                end if
                print*, "Triangle", j, "is opposite triangle of", i, "and will be collapsed"
                triangle_elimination(m) = .true.
                ntri_new = ntri_new - 1
              end if
            end do
          end do
          ! Re-assign coordinate of vertex1
          vcoord(:,vindex(1,i)) = ( old_vcoord(:,vindex(1,i)) + old_vcoord(:,vindex(3,i)) ) / 2.d0
          ! Collapse triangle
          triangle_elimination(i) = .true.
          ntri_new = ntri_new - 1
          ! Remove vertex3, become vertex1
          eliminated_vertex = vindex(3,i)
          ! Re-define index of vertices of neighbour triangles of collapsed vertex2
          do ij = 1, nghmax
            j = neighbour_triangles(i,3,ij)
            if(j==0) exit
            if(.not.triangle_elimination(j)) then
              do jv = 1, 3
                if( vindex(3,i).eq.vindex(jv,j) ) then
                  print*, "Re-locating vertex", jv, "of triangle", j, "neighbour of", i, 2
                  vindex(jv,j) = vindex(1,i)
                end if
              end do ! jv
            end if
          end do ! ij

          ! -- Exit just to remove one triangle for the moment
          exit
          ! --

        else ! at least one vertex is on tetra face, cycle
          !
          cycle
          !
        end if ! vertex on tetra face or not


      ! Edge 2-3 of the triangle is small
      !
      else if ( collapsed_edge == 3 ) then

        ! no vertex is on tetra face, move
        if( (ALL(.not.vert_border(vindex(2,i),:)) &
            .and. ALL(.not.vert_border(vindex(3,i),:))) &
            !! or both vertices are on same tetra face
            !.or. ANY(vert_border(vindex(3,i),:) .and. vert_border(vindex(2,i),:)) &
              ) then

          ! loop to check if any vertex is on two tetra-faces, i.e. corner of irrBZ
          k1 = 0
          k2 = 0
          do f = 1, ntetraf
            if(vert_border(vindex(2,i),f)) k1 = k1+1
            if(vert_border(vindex(3,i),f)) k2 = k2+1
            if(k1.ge.2 .or. k2.ge.2) cycle tri_loop
          end do

        print*, "Triangle",i,"collapsed from edge",23

        aux_vcoord(:,vindex(3,i)) = ( old_vcoord(:,vindex(2,i)) + old_vcoord(:,vindex(3,i)) ) / 2.d0
        aux_vcoord(:,vindex(2,i)) = ( old_vcoord(:,vindex(2,i)) + old_vcoord(:,vindex(3,i)) ) / 2.d0

          ! Check if collapse would lead to flip of any neighbour triangle of vertex2
          do ij = 1, nghmax

            if(neighbour_triangles(i,2,ij)==0) exit

            oldv1(:) = old_vcoord(:,vindex(2,neighbour_triangles(i,2,ij)))-old_vcoord(:,vindex(1,neighbour_triangles(i,2,ij)))
            oldv2(:) = old_vcoord(:,vindex(3,neighbour_triangles(i,2,ij)))-old_vcoord(:,vindex(1,neighbour_triangles(i,2,ij)))
            oldnormal(:) = cross(oldv2(:),oldv1(:))
            oldnormal(:) = oldnormal(:)/norma(oldnormal)

            newv1(:) = aux_vcoord(:,vindex(2,neighbour_triangles(i,2,ij)))-aux_vcoord(:,vindex(1,neighbour_triangles(i,2,ij)))
            newv2(:) = aux_vcoord(:,vindex(3,neighbour_triangles(i,2,ij)))-aux_vcoord(:,vindex(1,neighbour_triangles(i,2,ij)))
            newnormal(:) = cross(newv2(:),newv1(:))
            newnormal(:) = newnormal(:)/norma(newnormal)

            if(dot_product(oldnormal, newnormal).lt.0.0_dp) then
              print*, "Triangle will not collapse in order to avoid flip of a neighbour triangle"
              cycle tri_loop
            end if

          end do

          ! Check if collapse would lead to flip of any neighbour triangle of vertex3
          do ij = 1, nghmax

            if(neighbour_triangles(i,3,ij)==0) exit

            oldv1(:) = old_vcoord(:,vindex(2,neighbour_triangles(i,3,ij)))-old_vcoord(:,vindex(1,neighbour_triangles(i,3,ij)))
            oldv2(:) = old_vcoord(:,vindex(3,neighbour_triangles(i,3,ij)))-old_vcoord(:,vindex(1,neighbour_triangles(i,3,ij)))
            oldnormal(:) = cross(oldv2(:),oldv1(:))
            oldnormal(:) = oldnormal(:)/norma(oldnormal)

            newv1(:) = aux_vcoord(:,vindex(2,neighbour_triangles(i,3,ij)))-aux_vcoord(:,vindex(1,neighbour_triangles(i,3,ij)))
            newv2(:) = aux_vcoord(:,vindex(3,neighbour_triangles(i,3,ij)))-aux_vcoord(:,vindex(1,neighbour_triangles(i,3,ij)))
            newnormal(:) = cross(newv2(:),newv1(:))
            newnormal(:) = newnormal(:)/norma(newnormal)

            if(dot_product(oldnormal, newnormal).lt.0.0_dp) then
              print*, "Triangle will not collapse in order to avoid flip of a neighbour triangle"
              cycle tri_loop
            end if

          end do

        ! Collapse triangle sharing the removed edge
        do ij = 1, nghmax
          j = neighbour_triangles(i,2,ij)
          if(j==0) exit
            do mj = 1, nghmax
              m = neighbour_triangles(i,3,mj)
              if(m==0) exit
              ! Triangle m is neighbour of the two collapsed vertices, i.e. is the opposite triangle
              if(m==j) then
                if(     ANY(vert_border(vindex(1,m),:)) &
                   .or. ANY(vert_border(vindex(2,m),:)) &
                   .or. ANY(vert_border(vindex(3,m),:)) &
                   ) then
                  print*, "Triangle will not be collapsed as its opposite is on a tetra face"
                  cycle tri_loop
                end if
                print*, "Triangle", j, "is opposite triangle of", i, "and will be collapsed"
                triangle_elimination(m) = .true.
                ntri_new = ntri_new - 1
              end if
            end do
        end do
        ! Re-assign coordinate of vertex2
        vcoord(:,vindex(2,i)) = ( old_vcoord(:,vindex(2,i)) + old_vcoord(:,vindex(3,i)) ) / 2.d0
        ! Collapse triangle
        triangle_elimination(i) = .true.
        ntri_new = ntri_new - 1
        ! Remove vertex3, become vertex2
        eliminated_vertex = vindex(3,i)
        ! Re-define index of vertices of neighbour triangles of collapsed vertex2
        do ij = 1, nghmax
          j = neighbour_triangles(i,3,ij)
          if(j==0) exit
          if(.not.triangle_elimination(j)) then
              do jv = 1, 3
                if( vindex(3,i).eq.vindex(jv,j) ) then
                  print*, "Re-locating vertex", jv, "of triangle", j, "neighbour of", i, 2
                  vindex(jv,j) = vindex(2,i)
                end if
              end do
          end if
        end do

        ! -- Exit just to remove one triangle for the moment
        exit
        ! --

        else ! at least one vertex is on tetra face, cycle
          !
          cycle
          !
        end if ! any vertex at IBZ face or not

      end if ! loop over edge_collapse possibilities

    end do tri_loop


    ! Define new vertex coordinate and lists
    !
    nvert_new = 0
    do iv = 1, nvert_old
      if(iv.eq.eliminated_vertex) cycle
      nvert_new = nvert_new+1
      new_vcoord(:,nvert_new) = vcoord(:,iv)
    end do

    j = 0
    do i = 1, ntri_old
      if(triangle_elimination(i)) cycle
      j = j+1
      new_vindex(:,j) = vindex(:,i)
      do ij = 1, 3
        if(new_vindex(ij,j).gt.eliminated_vertex) new_vindex(ij,j) = new_vindex(ij,j)-1
      end do
    end do

  end subroutine collapse_triangles


  subroutine collapse_condition(vindex1, vindex2, vindex3, vcoord, eps, clpsd_edge)

    use intw_matrix_vector, only: norma

    implicit none

    ! I/O
    integer, intent(in) :: vindex1, vindex2, vindex3
    real(dp), dimension(:,:), intent(in) :: vcoord
    real(dp), intent(in) :: eps
    integer, intent(out) :: clpsd_edge


    clpsd_edge = 0

    if(norma(vcoord(:,vindex1)-vcoord(:,vindex2)) / &
       (norma(vcoord(:,vindex1)-vcoord(:,vindex2))+&
        norma(vcoord(:,vindex1)-vcoord(:,vindex3))+&
        norma(vcoord(:,vindex2)-vcoord(:,vindex3))) .lt. eps) then

      clpsd_edge = 1

    else if(norma(vcoord(:,vindex1)-vcoord(:,vindex3)) / &
       (norma(vcoord(:,vindex1)-vcoord(:,vindex2))+&
        norma(vcoord(:,vindex1)-vcoord(:,vindex3))+&
        norma(vcoord(:,vindex2)-vcoord(:,vindex3))) .lt. eps) then

      clpsd_edge = 2

    else if(norma(vcoord(:,vindex2)-vcoord(:,vindex3)) / &
       (norma(vcoord(:,vindex1)-vcoord(:,vindex2))+&
        norma(vcoord(:,vindex1)-vcoord(:,vindex3))+&
        norma(vcoord(:,vindex2)-vcoord(:,vindex3))) .lt. eps) then

      clpsd_edge = 3

    end if

  end subroutine collapse_condition


  subroutine tangential_relaxation(relax_vinface, eps_vinface, eps_dupv, niter, verbose, nrpts, irvec, ndegen, alat, ag, bg, &
                                   nsym, s, TR_sym, num_wann, ham_r, ibnd, in_nvert, in_ntri, nibz_faces, ibz_faces, verts_ibz, &
                                   v_index, v_coord)

    use intw_matrix_vector, only: ainv, norma, cross
    use triFS_isosurface, only: velocity

    implicit none

    !I/O variables
    logical, intent(in) :: relax_vinface
    real(dp), intent(in) :: eps_vinface, eps_dupv
    integer, intent(in) :: niter
    logical, intent(in) :: verbose !
    integer, intent(in) :: nrpts, irvec(1:3,nrpts), ndegen(nrpts) ! Wannier/Fourier space variables
    real(dp), intent(in) :: alat, ag(1:3,1:3), bg(1:3,1:3) ! real and reciprocal lattice vectors
    integer, intent(in) :: nsym, s(1:3,1:3,nsym) ! symmetry variables
    logical, intent(in) :: TR_sym ! time-reversal symmetry
    integer, intent(in) :: num_wann ! number of Wannier functions on hr
    complex(dp), intent(in) :: ham_r(num_wann,num_wann,nrpts) ! Hamiltonian in Wannier basis
    integer, intent(in) :: ibnd,in_nvert,in_ntri,nibz_faces
    integer, dimension(:,:), intent(in) :: ibz_faces ! IBZ big tetra faces
    real(dp), dimension(:,:), intent(in) :: verts_ibz ! IBZ big tetra vertex coordinates
    integer, dimension(:,:), intent(in) :: v_index
    real(dp), dimension(:,:), intent(inout) :: v_coord

    ! local variables
    integer :: j, iv, jv, it
    integer, dimension(:,:), allocatable :: neighbours
    logical, dimension(:), allocatable :: vert_on_corner
    logical, dimension(:,:), allocatable :: vert_on_border
    real(dp), dimension(:,:), allocatable :: vert_normal
    real(dp), dimension(:,:,:), allocatable :: vface
    real(dp), dimension(1:3) :: meanbary, v_w, vcoord_crys, dvec
    real(dp), dimension(1:3,1:3) :: face
    logical, allocatable :: has_SplusG(:), SplusG_pair(:,:,:,:,:,:,:)
    logical, allocatable :: vertex_related(:), relaxed_vertex(:)
    ! parameters
    integer, parameter :: nghmax = 500
    !real(dp), parameter :: epsvert = 1.0E-6_dp
    real(dp), parameter :: lambda = 0.1_dp ! Parameter controling vertex movement


    write(*,*) "Tangential relaxation... "

    !--------------- Define vectors of IBZ big tetra faces
    !
    allocate(vface(3,3,nibz_faces))
    call faces_vectors(nibz_faces, ibz_faces, verts_ibz, vface)


    !--------------- Find neighbour triangles of each vertex
    !
    allocate(neighbours(in_nvert,nghmax))
    neighbours = 0
    call find_neighbours_vertex_list(in_ntri, in_nvert, v_index, neighbours, nghmax)


    !--------------- Check which vertices are on border of IBZ
    !
    allocate(vert_on_border(in_nvert,nibz_faces), vert_on_corner(in_nvert))
    call check_border(in_nvert, v_coord(:,:), nibz_faces, vface, nsym, s, TR_sym, &
                      eps_vinface, vert_on_border, vert_on_corner)


    !--------------- Detect which vertices have SplusG pairs
    !
    allocate(has_SplusG(in_nvert), SplusG_pair(in_nvert,in_nvert,nsym,-1:1,-1:1,-1:1,1:2))
    allocate(vertex_related(in_nvert), relaxed_vertex(in_nvert))
    has_SplusG = .false.
    vertex_related = .false.
    SplusG_pair = .false.
    do iv = 1, in_nvert
      !
      if (vertex_related(iv)) cycle ! a pair of this vertex has been already assigned
      !
      do jv = 1, in_nvert
        !
        if(jv.eq.iv) cycle
        !
        call vertices_related_by_SplusG(eps_dupv, v_coord(1:3,iv), v_coord(1:3,jv), bg, &
                                        nsym, s, TR_sym, SplusG_pair(iv,jv,:,:,:,:,:))
        !
        if (ANY(SplusG_pair(iv,jv,:,:,:,:,:))) then
          !
          print*, "Vertices iv and jv related:", iv, jv
          has_SplusG(iv) = .true.
          vertex_related(jv) = .true.
          !
        end if
      end do ! jv
    end do ! iv


    !--------------- Vertex relocation
    !
    allocate(vert_normal(3,in_nvert))
    do it = 1, niter
      !
      relaxed_vertex(:) = .false.
      v_loop: do iv = 1, in_nvert
        !
        if (relaxed_vertex(iv)) cycle
        !
        if(ANY(SplusG_pair(iv,iv,:,:,:,:,:))) cycle ! don't relax vertex which has S+G with itself.
        !
        ! Compute normal vector of vertex with the velocity vector
        vcoord_crys(1:3) = matmul(ainv(bg), v_coord(1:3,iv))
        !call velocity_sym(nrpts, irvec, ndegen, ag, bg, nsym, s, TR_sym, num_wann, ibnd, ham_r, v_w, vcoord_crys)
        !! JL - abiadura simetrizatu gabe, bakarrik plano tangentziala lortzeko da.
        v_w = 0.0_dp
        call velocity(nrpts, irvec, ndegen, alat, ag, bg, num_wann, ibnd, ham_r, v_w, vcoord_crys)
        v_w = matmul(transpose(ag), v_w)
        v_w = matmul(bg, v_w)
        !!! JL
        vert_normal(1:3,iv) = v_w(1:3) / norma(v_w(1:3))
        !
        ! Obtain mean value of barycenters of neighbour triangles
        call mean_barycenter(iv, neighbours, nghmax, v_index, v_coord, meanbary)
        !
        ! Define displacement vector
        !
        if (vert_on_corner(iv)) then ! don't move vertices on corners
          dvec(:) = 0.0_dp
        else if (ANY(vert_on_border(iv,:))) then
          if(relax_vinface) then ! restrict movement to the face plane
            do j = 1, nibz_faces
              if (vert_on_border(iv,j)) then
                face(:,1) = vface(:,1,j)
                face(:,2) = vface(:,2,j)
                face(:,3) = cross(vface(:,2,j),vface(:,1,j))
                face(:,3) = face(:,3) / norma(face(:,3))
                !dvec(:) = lambda * ( (meanbary(:)-v_coord(:,iv)) - dot_product(vert_normal(:,iv),(meanbary(:)-v_coord(:,iv)))*vert_normal(:,iv) )
                dvec(:) = ( (meanbary(:)-v_coord(:,iv)) - dot_product(vert_normal(:,iv), (meanbary(:)-v_coord(:,iv)))*vert_normal(:,iv) )
                !write(*,'(6E18.10)'), dvec(:), dot_product(face(:,3),dvec(:))*face(:,3)
                dvec(:) = ( dvec(:) - dot_product(face(:,3), dvec(:))*face(:,3) )
                ! Motelago mugitu aurpegietako bertizeak
                dvec(:) = lambda / 4.0_dp * dvec(:)
              end if ! vert_on_border(iv,j)
            end do !j
          else ! do not relax vertices at faces
            dvec(:) = 0.0_dp
          end if ! relax_vinface
        else ! vertex is not on face, move normally
          dvec(:) = lambda * ( (meanbary(:)-v_coord(:,iv)) - dot_product(vert_normal(:,iv), (meanbary(:)-v_coord(:,iv)))*vert_normal(:,iv) )
        end if
        !
        ! safety condition to avoid errors
        if(sum(abs(dvec(:))).gt.10.0_dp) then
          print*, "Displacement is too big, something went wrong!!"
          dvec = 0.0_dp
        end if
        !
        !! Print displacements if needed
        !if (verbose) then
        !  if(it==1) print*, dvec(:)
        !  if(it==10) print*, dvec(:)
        !  if(it==100) print*, dvec(:)
        !  if(it==200) print*, dvec(:)
        !  if(it==1000) print*, dvec(:)
        !  if(it==10000) print*, dvec(:)
        !  if(it==100000) print*, dvec(:)
        !end if
        !
        ! Move vertex
        v_coord(:,iv) = v_coord(:,iv) + dvec(:)
        relaxed_vertex(iv) = .true.
        !
        ! Move any symmetry related vertex accordingly
        if (has_SplusG(iv)) then
          !
          do jv = 1, in_nvert
            if(ANY(SplusG_pair(iv,jv,:,:,:,:,:))) then ! jv is pair of iv
              !
              v_coord(1:3,jv) = rotate_to_SplusG(v_coord(1:3,iv), bg, nsym, s, SplusG_pair(iv,jv,:,:,:,:,:))
              relaxed_vertex(jv) = .true.
              !
            end if
          end do ! jv
        end if ! has_SplusG(iv)
        !
      end do v_loop

      ! Print iteration if needed
      if(verbose) then
        if(it==1) then
          write(*, *) "  1 iteration completed"
        else if(it==10) then
          write(*, *) "  10 iterations completed"
        else if(it==100) then
          write(*, *) "  100 iterations completed"
        else if(it==200) then
          write(*, *) "  200 iterations completed"
        else if(it==1000) then
          write(*, *) "  1000 iterations completed"
        else if(it==10000) then
          write(*, *) "  10000 iterations completed"
        else if(it==100000) then
          write(*, *) "  100000 iterations completed"
        end if
      end if

    end do ! niter

    write(*, *) "Tangential relaxation finished"

    !--------------- Deallocate variables
    !
    deallocate(vface)
    deallocate(vert_on_border, vert_on_corner)
    deallocate(has_SplusG, SplusG_pair, vertex_related, relaxed_vertex)
    deallocate(vert_normal, neighbours)

  end subroutine tangential_relaxation


  subroutine mean_barycenter(vi, neighbour_triangles, nghmax, vindex, vcoord, mean_bary)

    use intw_matrix_vector, only: cross

    implicit none

    ! I/O
    integer, intent(in) :: vi
    integer, dimension(:,:), intent(in) :: neighbour_triangles
    integer, intent(in) :: nghmax
    integer, dimension(:,:), intent(in) :: vindex
    real(dp), dimension(:,:), intent(in) :: vcoord
    real(dp), dimension(3), intent(out) :: mean_bary

    ! Local
    integer :: ni, ngh
    real(dp) :: tri_area, sum_area
    real(dp), dimension(3) :: bary, v1, v2


    ngh = 0
    mean_bary(:) = 0.0_dp
    sum_area = 0.0_dp
    do ni = 1, nghmax
      if(neighbour_triangles(vi,ni)==0) exit
      ngh = ngh+1
      bary(:) = (vcoord(:,vindex(1,neighbour_triangles(vi,ngh)))   &
                 + vcoord(:,vindex(2,neighbour_triangles(vi,ngh))) &
                 + vcoord(:,vindex(3,neighbour_triangles(vi,ngh)))) / 3.0_dp

      v1(:) = ( vcoord(:,vindex(2,neighbour_triangles(vi,ngh))) - vcoord(:,vindex(1,neighbour_triangles(vi,ngh))) )
      v2(:) = ( vcoord(:,vindex(3,neighbour_triangles(vi,ngh))) - vcoord(:,vindex(1,neighbour_triangles(vi,ngh))) )
      tri_area = sqrt(sum(cross(v1,v2)**2)) / 2.0_dp
      !tri_area = 1.0_dp
      sum_area = sum_area + tri_area

      mean_bary(:) = mean_bary(:) + tri_area*bary(:)
    end do
    mean_bary(:) = mean_bary(:) / sum_area

  end subroutine mean_barycenter

end module triFS_mesh_opt
