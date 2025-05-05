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
module triFS_geometry

  !------------------------------------------------------------------!
  ! Module that contains all the necessary variables and subroutines !
  ! related to triangles and tetrahedra.                             !
  !------------------------------------------------------------------!

  use kinds, only: dp

  public :: wigner_seitz_cell, plot_poly, polyhedra_off, plot_tetra_off, normal_of_tri, &
            tetrasym, rot_tetra, compact_tetra, overlap_tetra, overlap_tri, overlap_edge, overlap_face, &
            overlap_faces_int, equal_tetra, tetraIBZ_2_vert_faces_edges, irr_faces, sort_edges, &
            write_node_list_face, nodes_on_face_and_rotate_to_plane, triangulate_faces, add_nodes_IBZ_volume
  private


contains

  subroutine wigner_seitz_cell(vec, verb)

    use intw_utility, only: find_free_unit
    use intw_matrix_vector, only: ainv, det, norma, cross

    implicit none

    logical, intent(in) :: verb ! write info file or not
    real(dp), dimension(3,3), intent(in) :: vec ! lattice vectors (#vec, comp)

    integer, parameter :: max_vert = 400, &
                          max_faces = 100, &
                          max_ver_inface = 100 !&
    !real(dp), parameter :: eps1 = 1.0d-09
    real(dp), parameter :: eps1 = 1.0E-07_dp
    real(dp), dimension(4,26) :: eqplane !equation of plane: ax+by+cz=d; eqplane = (a, b, c, d)
    real(dp), dimension(3,max_vert) :: vertex !vertices of polihedra (with repetition)
    integer, dimension(3,max_vert) :: vert_index ! indices of the 3 planes producing each vertex
    integer, dimension(max_faces) :: dif_plane_list !list of different planes
    integer, dimension(max_ver_inface,max_faces) :: ind_verts_indifplane ! indexes of vertices in each plane before cleaning up
    integer, dimension(max_ver_inface,max_faces) :: ind_verts_indifplane_final ! indexes of vertices in each plane after cleaning up
    integer, dimension(max_ver_inface) :: num_verts_inplane !number of vertices in a plane before cleaning up
    integer, dimension(max_ver_inface) :: num_verts_inplane_final !number of vertices in a plane after cleaning up
    logical, dimension(max_ver_inface,max_faces) :: active ! related to ind_verts_indifplane
    integer :: num_poi, num_vert, nnplanes, nverplai, nverplai_final, num_plane
    real(dp) :: di, dj, dk, dist1, dist2
    integer :: i, j, k, ii, jj, l, iv, jv, kv, vrtIV, vrtJV, vrtKV
    real(dp) :: mat(3, 3), dd(3), pvert(3), vector1(3), vector2(3), normalvector(3)
    logical :: isnew, isplane, equalIJ, equalIK, equalJK, areequal
    integer :: io_unit1, info_unit


   ! ------------------------------------------------------
   !   vec(i,j) is a 3x3 matrix with the lattice vectors
   !        i = vector number;  b_i,  i = 1,2,3
   !        j = component; j = x,y,z
   ! -------------------------------------------------------

    !-- Check
    !    print*, 'vec'
    !    do i = 1, 3
    !     write(*, '(a,i1,a,3(1x,f14.8))') 'b', i, ':', ( vec(j,i), j = 1, 3 )
    !    end do
    !-- End check

   ! ---------------------------------------------------------------
   !  Definition of the WS cell:
   !  -Find neighbour points from the origin; 26 points in total
   !  -Then construct equations of the planes that go through P_0 and are
   !   perpendicular to OP0, --> OP0·PP0 = 0; (P = x,y,z); P0 = (x0,y0,z0)
   !  -The plane perpendicular to OP0, passing at half distance between
   !   0 and P0 is:
   !      x0.x+y0.y+z0.z = 0.5(x0^2+y0^2+z0^2)
   !       ax + by + cz = d
   !      eqplane = (a,b,c,d)
   ! -----------------------------------------------------------------

    io_unit1 = find_free_unit()
    open(unit=io_unit1, file='polyhedra.dat', status='unknown')
    if (verb) then
      info_unit = find_free_unit()
      open(unit=info_unit, file='info_ws.dat', status='unknown')
    end if

    num_poi = 0
    do i = -1, 1
      do j = -1, 1
        do k = -1, 1
          if ( i /= 0 .or. j /= 0 .or. k /= 0 ) then !disregard the origin
            !di = dble(i)
            !dj = dble(j)
            !dk = dble(k)
            di = i
            dj = j
            dk = k
            num_poi = num_poi+1
            eqplane(1,num_poi) = di*vec(1,1) + dj*vec(2,1) + dk*vec(3,1)
            eqplane(2,num_poi) = di*vec(1,2) + dj*vec(2,2) + dk*vec(3,2)
            eqplane(3,num_poi) = di*vec(1,3) + dj*vec(2,3) + dk*vec(3,3)
            !eqplane(4,num_poi) = 0.5*(norma(eqplane(1:3,num_poi)))**2.0d0
            ! JL
            eqplane(4,num_poi) = 0.5_dp*(norma(eqplane(1:3,num_poi)))**2
          endif
        enddo
      enddo
    enddo

   ! Vertices are the intersections of three planes. Check all possibilities
   ! Solve for X: A X = D
   ! If det(A) = 0, planes are parallel. There is no intersection.
    num_vert = 0
    do i = 1, num_poi-2
      mat(1,1:3) = eqplane(1:3,i)
      dd(1) = eqplane(4,i)

      do j = i+1, num_poi-1
        mat(2,1:3) = eqplane(1:3,j)
        dd(2) = eqplane(4,j)

        do k = j+1, num_poi
          mat(3,1:3) = eqplane(1:3,k)
          dd(3) = eqplane(4,k)

          ! if ( abs(det(mat)) >= 1.d-13 ) then
          ! JL
          if ( abs(det(mat)) >= 1.0E-07_dp ) then
            pvert = matmul(ainv(mat), dd)

            ! Check whether vertex is closer to Gamma than
            ! any other point. If so, it is a new vertex
            isnew = .true.
            l = 1
            do
            if ( .not. ( l<=num_poi .and. isnew ) ) exit
              dist1 = norma(pvert)
              dist2 = norma(pvert-eqplane(1:3,l))
            if ((dist1-dist2) > eps1) isnew = .false.
              l = l+1
            end do

            if (isnew) then
              num_vert = num_vert+1
              vertex(1:3,num_vert) = pvert(1:3)
              vert_index(1,num_vert) = i
              vert_index(2,num_vert) = j
              vert_index(3,num_vert) = k
              ! --- Check: List of vertices closer to Gamma (with repetition)
              if (verb) then
                write(info_unit, "(i4,3f12.6,3i4)") num_vert, vertex(1:3,num_vert), vert_index(1:3,num_vert)
              end if
              ! --- End check
            end if

          end if ! det /=0
        end do ! k
      end do ! j
    end do ! i

   ! The number of different indexes that appear in vert_index(1:3,:)
   ! is the total number of planes a polyhedra is made of --nnplanes--
   ! The 3 indexes of the first vertex are always new and different
    nnplanes = 0
    do i = 1, 3
      nnplanes = nnplanes + 1
      dif_plane_list(nnplanes) = vert_index(i,1)
    end do

    do i = 2, num_vert
      do j = 1, 3
        isnew = .true.
        l = 1
        do
          if ( .not. ( isnew .and. l<=nnplanes ) ) exit
          if (dif_plane_list(l) == vert_index(j,i) ) isnew = .false.
          l = l+1
        end do ! open do

        if (isnew) then
          nnplanes = nnplanes+1
          dif_plane_list(nnplanes) = vert_index(j,i)
        end if

      end do !j
    end do !i
    if (verb) then
      write(*, '(a,i6,a,i3,a)') 'Found', num_vert, ' possible vertices, between ', &
                                nnplanes, ' different planes'
      write(info_unit, *)
      write(info_unit, '(a,i6,a,i3,a)') 'Found', num_vert, ' possible vertices, between ', &
                                        nnplanes, ' different planes'
      write(info_unit, *)

      ! --- Check: List of different planes
      do i = 1, nnplanes
        write(info_unit, *) i, dif_plane_list(i)
      end do
      ! --- End check
    end if

  ! Vertices which have one of the three indexes in common belong to the same plane.
  ! But these vertices might be identical or collinear.
  ! Find first all vertices belonging to each plane. Then remove repetitions and collinear
  ! points.
    active = .false.
    num_plane = 0
    do i = 1, nnplanes
      isplane = .true.
      !Relate each vertex with the list of different planes
      nverplai = 0
      do j = 1, num_vert
        do jj = 1, 3
          if ( vert_index(jj,j) == dif_plane_list(i) ) then
            !vertex j belongs to plane i
            nverplai = nverplai+1
            ind_verts_indifplane(nverplai,i) = j
            active(nverplai,i) = .true.
          end if
        end do !jj
      end do ! j
      num_verts_inplane(i) = nverplai

      ! --- Check: Vertices in each plane
      if (verb) then
        write(info_unit, *)
        write(info_unit, '(a,2i4)') 'Plane #', i, num_verts_inplane(i)
        write(info_unit, '(100i4)') ( ind_verts_indifplane(ii,i), ii = 1, num_verts_inplane(i) )
      end if
      ! --- End check

      ! If total number of vertices in plane < 3, it is not a plane
      if ( num_verts_inplane(i) < 3 ) isplane = .false.
      ! Otherwise remove collinear points and identical points
      if (isplane) then

        do iv = 1, num_verts_inplane(i)-2
          if (active(iv,i)) then
            do jv = iv+1, num_verts_inplane(i)-1
              if (active(jv,i)) then
                do kv = jv+1, num_verts_inplane(i)
                  if (active(kv,i)) then

                    vrtIV = ind_verts_indifplane(iv,i)
                    vrtJV = ind_verts_indifplane(jv,i)
                    vrtKV = ind_verts_indifplane(kv,i)
                    !remove identical points
                    equalIJ = norma( vertex(:,vrtIV)-vertex(:,vrtJV) ) <= eps1
                    equalIK = norma( vertex(:,vrtIV)-vertex(:,vrtKV) ) <= eps1
                    equalJK = norma( vertex(:,vrtJV)-vertex(:,vrtKV) ) <= eps1

                    areequal = .true.
                    if ( equalIJ .and. equalIK ) then
                      ! Three vertices are identical.
                      ! Keep the one with the largest index iv < jv < kv
                      active(iv,i) = .false.
                      active(jv,i) = .false.
                    else if (equalIJ) then
                      active(iv,i) = .false.
                    else if (equalIK) then
                      active(iv,i) = .false.
                    else if (equalJK) then
                      active(jv,i) = .false.
                    else
                      areequal = .false.
                    end if
                    !
                    !
                    if ( .not. areequal) then
                      vector1 = vertex(:,vrtJV)-vertex(:,vrtIV)
                      vector2 = vertex(:,vrtKV)-vertex(:,vrtIV)
                      normalvector = cross(vector1, vector2)

                      if ( norma(normalvector) < eps1 ) then
                        ! vertices are collinear
                        if (verb) write(info_unit, *) ' are collinear'
                        write(*, *) 'Found collinear vertices. Stopping '
                        if (verb) write(info_unit, *) 'Found collinear vertices. Stopping '
                        STOP
                        !Momentuz hau saltatu. Ez dirudi hau gertatu daitekeenik.
                      end if !if (norma(normalvector) < eps1)
                    end if ! not areequal
                    !
                    !

                  end if ! active kv
                end do ! kv
              end if ! active jv
            end do ! jv
          end if ! active iv
        end do ! iv


        !Store all active vertices in plane i (after repetitions removed)
        if (verb) write(info_unit, *)'Vertices without repetition'
        nverplai_final = 0
        do iv = 1, num_verts_inplane(i)
          if (active(iv, i)) then
            nverplai_final = nverplai_final+1
            ind_verts_indifplane_final(nverplai_final,i) = ind_verts_indifplane(iv,i)
            num_verts_inplane_final(i) = nverplai_final
            if (verb) write(info_unit, '(2i4,3f12.8)') nverplai_final, ind_verts_indifplane_final(nverplai_final,i), &
                                                       ( vertex(j,ind_verts_indifplane_final(nverplai_final,i)), j = 1, 3 )
          end if
        end do ! iv

        if (nverplai_final <= 2 .and. verb) write(info_unit, *) '------THIS IS NOT A FACE!'
        if (nverplai_final>2) then
          ! Write raw file (vertices might not be in order)
          if (verb) write(info_unit, *) '......... THIS IS A FACE!'
          num_plane = num_plane +1
          write(io_unit1, *) num_plane, nverplai_final
          do ii = 1, nverplai_final
            write(io_unit1, "(3(x,f18.15))") ( vertex(j,ind_verts_indifplane_final(ii,i)), j = 1, 3 )
          end do
        end if

      end if ! isplane

      ! Write raw file (vertices might not be in order)
      !if ( isplane .and. (nverplai_final>2) ) then
      ! write(info_unit, *) '......... THIS IS A FACE!'
      !  num_plane = num_plane +1
      !  write(io_unit1, *) num_plane, nverplai_final
      !  do ii = 1, nverplai_final
      !   write(io_unit1, "(3(x,f18.15))") ( vertex(j,ind_verts_indifplane_final(ii,i)), j = 1, 3 )
      !  end do
      ! end if

    end do ! i (list of diferent planes)

    if (verb) then
      write(info_unit, *)
      write(info_unit, '(a, i4,a)') 'The BZ has ', num_plane, " outer faces"
      write(*, '(a, i4,a)') 'The BZ has ', num_plane, " outer faces"

      write(info_unit, *) 'WS done'
      close(info_unit)
    end if
    write(*, *) 'WS done'
    write(*, *)
    close(io_unit1)

  end subroutine wigner_seitz_cell


  subroutine plot_poly()

    use intw_utility, only: find_free_unit, hpsort_real
    use intw_matrix_vector, only: norma

    implicit none

    integer, parameter :: nv_max = 200, npl_max = 30
    real(dp) :: vp(1:3,nv_max,npl_max), dist(nv_max)
    integer :: io_unit1, io_unit2, dummy, nvp(npl_max), ierr, ind(nv_max)
    integer :: vp_dif_index(nv_max,npl_max)
    integer :: ndif, npl, nn
    real(dp) :: vp_dif(1:3,nv_max), vnew(1:3)
    integer :: i, j, k
    logical :: done(nv_max,nv_max)


    done = .false.
    io_unit1 = find_free_unit()
    open(io_unit1, file="polyhedra.dat", status="unknown")

    io_unit2 = find_free_unit()
    open(io_unit2, file="polyplot.dat", status="unknown")

    npl = 0
    do
      read(unit=io_unit1, fmt=*, iostat=ierr) dummy, nn
      if (ierr/=0) exit
      npl = npl+1
      nvp(npl) = nn
      do i = 1, nvp(npl) !each plane
        read(io_unit1, *) ( vp(j,i,npl), j = 1, 3 )
      enddo
    enddo

    ndif = 1
    vp_dif(1:3,1) = vp(1:3,1,1)
    do i = 1, npl
      do j = 1, nvp(i)
        !
        vnew(:) = vp(:,j,i)
        do k = 1, ndif
          if ( sum(abs(vnew-vp_dif(1:3,k)))<1.0E-7 ) then
            vp_dif_index(j,i) = k
            cycle
          endif
        enddo
        !
        ndif = ndif+1
        vp_dif(1:3,ndif)  = vnew(1:3)
        vp_dif_index(j,i) = ndif
        !
      enddo
    enddo

    do i = 1, npl
      do j = 1, nvp(i)
        do k = 1, nvp(i)
          dist(k) = norma(vp(:,j,i) - vp(:,k,i))
        enddo
        call hpsort_real(nvp(i), dist(1:nvp(i)), ind)
        do k = 2, 3
          ! if ( .not.done(vp_dif_index(ind(k),i),vp_dif_index(k,i)) ) then
          write(io_unit2, "(100f12.6)") vp(:,j,i), vp(:,ind(k),i) ! matmul(bg, vp(:,j,i)), matmul(bg, vp(:,ind(k),i))
          done(vp_dif_index(ind(k),i),vp_dif_index(k,i)) = .true.
          done(vp_dif_index(k,i),vp_dif_index(ind(k),i)) = .true.
          ! endif
        enddo
      enddo
    enddo
    close(unit=io_unit1)
    close(unit=io_unit2)

  end subroutine plot_poly


  subroutine polyhedra_off()

    use intw_utility, only: find_free_unit
    use intw_matrix_vector, only: norma

    implicit none

    integer  :: nvp(100)
    integer  :: i, j, k, l, iface, ierr
    real(dp) :: vlist(1:3,10000), v(1:3), vaux(1:3), dist2, dist, vvp(1:3,1000)
    real(dp) :: vec1(3), vec2(3), alpha, alpha2
    integer  :: nv
    integer, dimension(100,100)  :: ind
    integer :: io_unit


    io_unit = find_free_unit()
    open(io_unit, file="polyhedra.dat", status="unknown")

    k = 0
    nv = 0
    ind = 0

      do
        read(unit=io_unit, fmt=*, iostat=ierr) iface, nvp(iface)
        if (ierr/=0) exit
        do i = 1, nvp(iface)
          read(io_unit, *) (vvp(j,i), j = 1, 3)
        end do
        l = 1
        do
          if (l==nvp(iface)) exit
          if (l==1) then
            ! Find neighbour vertex to 1 and assign as second vertex
            !dist = huge(dist)
            alpha = 0.0_dp
            vec1(:) = vvp(:,2)-vvp(:,l)
            do i = 3, nvp(iface)
              !dist2 = norma(vvp(:,i)-vvp(:,l))
              !if (dist2<=dist) then
              !  dist = dist2
              !  j = i
              !  vaux(:) = vvp(:,i)
              !end if
              vec2(:) = vvp(:,i)-vvp(:,l)
              alpha2 = acos( dot_product(vec2, vec1) / (norma(vec2)*norma(vec1)+ 1.0E-07_dp)  )
              !print*, alpha, alpha2
              if (alpha2>alpha) then
                alpha = alpha2
                j = i
                vaux(:) = vvp(:,i)
              end if
            end do
            vvp(:,j) = vvp(:,2)
            vvp(:,2) = vaux(:)
            l = l+1
          else
            ! Find next vertex by looking for smallest angle w/r vector joining previous vertices
            alpha = 2.0_dp*acos(-1.0_dp)
            vec1(:) = vvp(:,l)-vvp(:,l-1)
            do i = l+1, nvp(iface)
              vec2(:) = vvp(:,i)-vvp(:,l)
              alpha2 = acos( dot_product(vec2, vec1) / (norma(vec2)*norma(vec1)+ 1.0E-07_dp)  )
              !print*, alpha, alpha2
              if (alpha2<alpha) then
                alpha = alpha2
                j = i
                vaux(:) = vvp(:,i)
              end if
            end do
            vvp(:,j) = vvp(:,l+1)
            vvp(:,l+1) = vaux(:)
            l = l+1
          end if
        end do

        i_loop: do i = 1, nvp(iface) !each plane
          k = k+1
          !read(io_unit, *) ( v(j), j = 1, 3 )
          v(:) = vvp(:, i)
          if (nv==0) then
            nv = nv+1
            vlist(:,nv) = v(:)
            ind(iface,i) = nv
          else
            do j = 1, nv
              if ((sum(abs(v(:)- vlist(:,j))))<1.0E-7_dp) then
                ind(iface,i) = j
                cycle i_loop
              end if
            end do
            nv = nv+1
            vlist(:,nv) = v(:)
            ind(iface,i) = nv
          end if
        enddo i_loop
    enddo

    close(unit=io_unit)

    io_unit = find_free_unit()
    open(io_unit, file="polyhedra.off", status="unknown")

    write(unit=io_unit, fmt="(a)") "OFF"
    write(unit=io_unit, fmt="(3i6)") nv, iface, 0
    write(unit=io_unit, fmt=*)
    do i = 1, nv
      write(unit=io_unit, fmt="(3f18.10)") vlist(:,i)
    end do
    do i = 1, iface
      write(unit=io_unit, fmt="(100I6)") nvp(i), ( ind(i,j)-1, j = 1, nvp(i) )
    end do
    close(io_unit)

  end subroutine polyhedra_off


  subroutine plot_tetra_off(tag_in, num_tetra, vertx_tetra, plot_surface)

    use intw_utility, only: find_free_unit
    use intw_matrix_vector, only: norma

    implicit none

    character(len=25), intent(in)   :: tag_in
    integer, intent(in) :: num_tetra ! Number of tetrahedra to be plotted
    real(dp), dimension(1:3,1:4,num_tetra), intent(in) :: vertx_tetra ! Coordinates of the four vertex of the tetrahedra
    logical, intent(in) :: plot_surface ! .true. if extra file is to be written with only the surface triangles

    ! local variables
    integer :: i, j, k, l, nv, indx(2000), io_unit, io_unit2
    real(dp) :: vlist(1:3,2000)


    l = 0
    nv = 0
    do i = 1, num_tetra
      j_loop: do j = 1, 4
        l = l+1
        !write(unit=io_unit, fmt='(3f12.6)') ( vertx_tetra(k,j,i), k = 1, 3 )
        if (nv==0) then
          nv = nv+1
          vlist(:,nv) = vertx_tetra(:,j,i)
          indx(l) = nv  ! vertex index
          cycle j_loop
        end if
        do k = 1, nv
          if (norma(vertx_tetra(:,j,i)-vlist(:,k))<10E-6) then ! check if vertex already stored
            indx(l) = k
            cycle j_loop
          end if
        end do
        nv = nv+1
        vlist(:,nv) = vertx_tetra(:,j,i)
        indx(l) = nv
      end do j_loop
    end do

    io_unit = find_free_unit()
    open(unit=io_unit, file=TRIM(tag_in), status="unknown")
    write(unit=io_unit, fmt="(a)") "OFF"
    write(unit=io_unit, fmt="(3I6)") nv, 4*num_tetra, 0
    write(unit=io_unit, fmt=*)

    if (plot_surface) then
      io_unit2 = find_free_unit()
      open(unit=io_unit2, file="faces_BZ.off", status="unknown")
      write(unit=io_unit2, fmt="(a)") "OFF"
      write(unit=io_unit2, fmt="(3I6)") nv, num_tetra, 0
      write(unit=io_unit2, fmt=*)
    end if

    ! Vertices of tetrahedra
    do i = 1, nv
      write(unit=io_unit, fmt="(3f12.6)") ( vlist(j,i), j = 1, 3 )
      if (plot_surface) write(unit=io_unit2, fmt="(3f12.6)") ( vlist(j,i), j = 1, 3 )
    end do

    ! Triangular faces of tetrahedra
    l = 1
    do i = 1, num_tetra
      write(unit=io_unit, fmt="(4I6)") 3, indx(l)-1, indx(l+2)-1, indx(l+1)-1
      write(unit=io_unit, fmt="(4I6)") 3, indx(l)-1, indx(l+3)-1, indx(l+1)-1
      write(unit=io_unit, fmt="(4I6)") 3, indx(l)-1, indx(l+2)-1, indx(l+3)-1
      write(unit=io_unit, fmt="(4I6)") 3, indx(l+1)-1, indx(l+2)-1, indx(l+3)-1
      if (plot_surface) write(unit=io_unit2, fmt="(4I6)") 3, indx(l+1)-1, indx(l+2)-1, indx(l+3)-1
      l = l+4
    end do

    close(unit=io_unit)
    if (plot_surface) close(io_unit2)

  end subroutine plot_tetra_off


  function normal_of_tri(vertices)

    use intw_matrix_vector, only: norma, cross

    implicit none

    real(dp), dimension(1:3,1:3), intent(in) :: vertices
    real(dp), dimension(1:3) :: normal_of_tri
    real(dp), dimension(1:3) :: v1, v2


    v1 = vertices(1:3,2) - vertices(1:3,1)
    v2 = vertices(1:3,3) - vertices(1:3,1)

    normal_of_tri = cross(v2, v1) / norma(cross(v2, v1))

  end function normal_of_tri


  subroutine tetrasym(bcell, nsym, s, TR_sym, nt_max, n_b_tetra_all, b_tetra_all, n_b_tetra_irr, b_tetra_irr, tetra_equiv, tetra_symlink)
    ! Subroutine that divides the BZ in tetrahedra, and finds the irreducible ones.

    use intw_utility, only: find_free_unit
    use intw_matrix_vector, only: cross

    implicit none

    real(dp), intent(in) :: bcell(1:3,1:3)
    integer, intent(in) :: nsym
    integer, intent(in) :: s(1:3,1:3,nsym)
    logical, intent(in) :: TR_sym
    integer, intent(in) :: nt_max
    real(dp), intent(out) :: b_tetra_all(1:3,1:4,nt_max )
    integer, intent(out) :: n_b_tetra_all
    real(dp), intent(out) :: b_tetra_irr(1:3,1:4,nt_max)
    integer, intent(out) :: n_b_tetra_irr, tetra_equiv(nt_max), tetra_symlink(nt_max,1:2)

    real(dp) :: b_tetra_aux(1:3,1:4,nt_max )
    integer  ::  n_b_tetra_aux
    real(dp) :: new_t(1:3,1:4), centre(1:3,24)
    integer  :: nvert, nfaces, nedges, iface
    integer  :: i, j, i_sym, it, jt, iv
    real(dp), dimension(:,:), allocatable :: vcoord
    integer, dimension(:,:), allocatable :: vindex
    integer, dimension(:), allocatable :: nvp
    real(dp) :: tirr_vol, t_vol, t1(3), t2(3), t3(3)
    logical :: t_done(nt_max)
    integer :: io_unit


    t_done = .false.

    !---------------------------------------------------------------
    ! Read BZ polyhedra from file
    !---------------------------------------------------------------

    io_unit = find_free_unit()
    !open(io_unit, file="polyhedra.dat", status="unknown")
    open(io_unit, file="polyhedra.off", status="unknown")
    read(unit=io_unit, fmt=*)
    read(unit=io_unit, fmt='(3I6)') nvert, nfaces, nedges
    read(unit=io_unit, fmt=*)

    allocate(vcoord(1:3,1:nvert))
    allocate(vindex(1:24,nfaces))
    allocate(nvp(nfaces))

    n_b_tetra_all = 0
    do iv = 1, nvert
      read(unit=io_unit, fmt="(3f18.10)") vcoord(:,iv)
    end do
    do iface = 1, nfaces
      read(unit=io_unit, fmt="(100I6)") nvp(iface), ( vindex(j,iface), j = 1, nvp(iface) )
    end do
    ! Make index go from 1 to nvert
    vindex = vindex+1

    !---------------------------------------------------------------
    ! Divide the BZ polyhedron in tetrahedra
    !---------------------------------------------------------------

    centre = 0.0_dp
    f_loop: do iface = 1, nfaces !each plane

      do iv = 1, nvp(iface)
        centre(1:3,iface) = centre(1:3,iface) + vcoord(1:3,vindex(iv,iface))
      end do
      centre(1:3,iface) = centre(1:3,iface) / nvp(iface)  ! center of each plane


      do iv = 1, nvp(iface)
        do j = 1, 2
          n_b_tetra_all = n_b_tetra_all + 1

          tetra_equiv  (n_b_tetra_all)   = -1 ! n_b_tetra_all                             ! Initialize equivalence
          tetra_symlink(n_b_tetra_all,1) =  1
          tetra_symlink(n_b_tetra_all,2) =  0

          b_tetra_all(1:3,1,n_b_tetra_all ) = (/0.0d0, 0.0d0, 0.0d0/)
          b_tetra_all(1:3,2,n_b_tetra_all ) = centre(1:3,iface)

          if ( iv==1 .and. j==1 ) then ! special case for first vertex
            b_tetra_all(1:3,3,n_b_tetra_all ) = (vcoord(1:3,vindex(1,iface)) + vcoord(1:3,vindex(nvp(iface),iface)) )/2.0_dp
            b_tetra_all(1:3,4,n_b_tetra_all ) = vcoord(1:3,vindex(1,iface))
          else if ( iv==nvp(iface) .and. j==2 ) then ! special case for last vertex
            b_tetra_all(1:3,3,n_b_tetra_all ) = vcoord(1:3,vindex(nvp(iface),iface))
            b_tetra_all(1:3,4,n_b_tetra_all ) = (vcoord(1:3,vindex(1,iface)) + vcoord(1:3,vindex(nvp(iface),iface)) )/2.0_dp
          else if ( iv/=1 .and. j==1 ) then
            b_tetra_all(1:3,3,n_b_tetra_all ) = (vcoord(1:3,vindex(iv,iface)) + vcoord(1:3,vindex(iv-1,iface)) )/2.0_dp
            b_tetra_all(1:3,4,n_b_tetra_all ) = vcoord(1:3,vindex(iv,iface))
          else if ( iv/=nvp(iface) .and. j==2 ) then
            b_tetra_all(1:3,3,n_b_tetra_all ) = vcoord(1:3,vindex(iv,iface))
            b_tetra_all(1:3,4,n_b_tetra_all ) = (vcoord(1:3,vindex(iv,iface)) + vcoord(1:3,vindex(iv+1,iface)) )/2.0_dp
          end if

        end do ! j
      end do ! iv

    enddo f_loop !planes
    close(io_unit)
    deallocate(vcoord, vindex, nvp)


    ! Filter and eliminate equal t
    n_b_tetra_aux = n_b_tetra_all
    b_tetra_aux   = b_tetra_all

    n_b_tetra_all = 0
    b_tetra_all   = 0

    aux_loop: do it = 1, n_b_tetra_aux

      do i = 1, 4
        new_t(1:3,i) = b_tetra_aux(1:3,i,it )
      end do

      do jt = 1, n_b_tetra_all
        if (equal_tetra(new_t,b_tetra_all(1:3,1:4,jt))) then
          cycle aux_loop
        endif
      end do

      n_b_tetra_all = n_b_tetra_all +1
      do i = 1, 4
        b_tetra_all(1:3,i,n_b_tetra_all ) = new_t(1:3,i)
      end do

    end do aux_loop


    !---------------------------------------------------------------
    ! Find irreducible tetrahedra sampling of the BZ (find the IBZ)
    !---------------------------------------------------------------

    n_b_tetra_irr = 1
    do i = 1, 4
      b_tetra_irr(1:3,i,1) = b_tetra_all(1:3,i,1 )
    enddo
    tetra_equiv(1) = 1
    tetra_symlink(1,1) = 1
    tetra_symlink(1,2) = 0 !

    it_loop: do it = 2, n_b_tetra_all

      do i = 1, 4
        new_t(1:3,i) = b_tetra_all(1:3,i,it )
      enddo

      do i = 1, n_b_tetra_irr
        do i_sym = 1, nsym
          if ( equal_tetra( new_t(1:3,1:4), rot_tetra(bcell,s(:,:,i_sym), b_tetra_irr(1:3,1:4,i)) ) ) then
            tetra_symlink(it,1) = i_sym
            tetra_symlink(it,2) = 0 !
            tetra_equiv(it) = i
            cycle it_loop
          end if
        end do
        !TR_symmetry
        do i_sym = 1, nsym
          if ( TR_sym .and. (equal_tetra ( - new_t(1:3,1:4), rot_tetra(bcell, s(:,:,i_sym), b_tetra_irr(1:3,1:4,i)) )) ) then
            tetra_symlink(it,1) = i_sym
            tetra_symlink(it,2) = 1 !
            tetra_equiv(it) = i
            cycle it_loop
          end if
        end do
      end do

      n_b_tetra_irr = n_b_tetra_irr +1
      tetra_equiv(it) = n_b_tetra_irr
      tetra_symlink(it,1) = 1
      tetra_symlink(it,2) = 0 !

      do i = 1, 4
        b_tetra_irr(1:3,i,n_b_tetra_irr ) = new_t(1:3,i)
      end do
    end do it_loop

    ! Compute volume of IBZ and BZ
    tirr_vol = 0.0_dp
    do it = 1, n_b_tetra_irr
      t1(:) = b_tetra_irr(:,2,it) - b_tetra_irr(:,1,it)
      t2(:) = b_tetra_irr(:,3,it) - b_tetra_irr(:,1,it)
      t3(:) = b_tetra_irr(:,4,it) - b_tetra_irr(:,1,it)
      tirr_vol = tirr_vol + abs(dot_product(t1, cross(t2, t3)))/6.0_dp
    end do
    write(*, *) "Volume of the irreducible BZ (2pi/alat)**3 =", tirr_vol
    ! Compute volume of IBZ
    t_vol = 0.0_dp
    do it = 1, n_b_tetra_all
      t1(:) = b_tetra_all(:,2,it) - b_tetra_all(:,1,it)
      t2(:) = b_tetra_all(:,3,it) - b_tetra_all(:,1,it)
      t3(:) = b_tetra_all(:,4,it) - b_tetra_all(:,1,it)
      t_vol = t_vol + abs(dot_product(t1, cross(t2, t3)))/6.0_dp
    end do
    write(*, *) "Volume of the full BZ (2pi/alat)**3 =", t_vol
    write(*, *) "BZ volume / IBZ volume =", t_vol / tirr_vol
    write(*, *)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Modification to test compact_tetra
    i = 0
    do it = 1, n_b_tetra_all
      if (i==n_b_tetra_irr) exit
      if (tetra_symlink(it,1)/=1) then
        b_tetra_irr(1:3,1:4,tetra_equiv(it)) = b_tetra_all(1:3,1:4,it)
        i = i+1
      end if
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine tetrasym


  function rot_tetra(bcell, symop, t)

    use intw_matrix_vector, only: ainv

    implicit none

    real(dp), intent(in) :: bcell(1:3,1:3)
    integer, intent(in) :: symop(3,3)
    real(dp), intent(in) :: t(1:3,1:4)

    real(dp) :: v(1:3)
    real(dp) :: rot_tetra(1:3,1:4), bcelli(3,3)
    integer :: i


    bcelli = ainv(bcell)

    do i = 1, 4
      v = matmul(bcelli, t(:,i))
      rot_tetra(1:3,i) = matmul(bcell, matmul(dble(symop(:,:)), v))
      ! rot_tetra(1:3,i) = matmul(s_cart(:,:,i_sym), t(:,i))
    end do

  end function rot_tetra


  subroutine compact_tetra(bcell, nsym, s, TR_sym, nt_max, n_b_tetra_irr, b_tetra_irr, n_b_tetra_all, b_tetra_all, tetra_equiv, tetra_symlink)

    implicit none

    real(dp), intent(in)    :: bcell(1:3,1:3)
    integer, intent(in)          :: nsym
    integer, intent(in)          :: s(1:3,1:3,nsym)
    logical, intent(in)          :: TR_sym
    integer, intent(in)          :: nt_max, n_b_tetra_irr, n_b_tetra_all
    real(dp), intent(inout) :: b_tetra_irr(1:3,1:4,nt_max), b_tetra_all(1:3,1:4,nt_max)
    integer, intent(inout)       :: tetra_equiv(nt_max), tetra_symlink(nt_max,1:2)

    !local
    integer       :: i, j, i_sym, nbt
    real(dp) :: b_tetra_irr_l(1:3,1:4,nt_max)
    integer       :: tetra_equiv_l(nt_max), tetra_symlink_l(nt_max,1:2)
    integer       :: overlap(nsym,2) ! 2 Time.R kontatzeko
    real(dp) :: t1(3,4), t2(3,4)
    integer       :: maxindex(2), it


    b_tetra_irr_l = b_tetra_irr
    tetra_equiv_l = tetra_equiv
    tetra_symlink_l = tetra_symlink

    nbt = 1 !lehenegoa ez dugu mugituko eta kontatzen dugu ia eginda dagoela.
    b_tetra_irr(1:3,1:4,1) = b_tetra_irr_l(1:3,1:4,1)

    ! bueltan b_tetra_irr izango do rotatutakoa ...

    do i = 2, n_b_tetra_irr
      overlap(:,1) = 0
      do i_sym = 1, nsym
        t1 = rot_tetra(bcell, s(:,:,i_sym), b_tetra_irr_l(1:3,1:4,i)) ! b_tetra_irr_l rotatu aurretik
        do j = 1, nbt
          t2 = b_tetra_irr(1:3,1:4,j) ! b_tetra_irr rotatuta, berria dagoeneko
          overlap(i_sym,1) = overlap(i_sym,1) + overlap_tetra(t1, t2)
        enddo
      enddo
      if (TR_sym) then
        overlap(:,2) = 0
        do i_sym = 1, nsym
          t1 = rot_tetra(bcell, s(:,:,i_sym), b_tetra_irr(1:3,1:4,i))
          do j = 1, nbt
            t2 = b_tetra_irr_l(1:3,1:4,j)
            overlap(i_sym,2) = overlap(i_sym,2) + overlap_tetra(t1, t2)
          enddo
        enddo
      end if
      maxindex = maxloc(overlap)

      nbt = nbt+1
      if (maxindex(2)==1) then
        i_sym = maxindex(1)
        b_tetra_irr(1:3,1:4,nbt) = rot_tetra(bcell, s(:,:,i_sym), b_tetra_irr_l(1:3,1:4,i))
      else if (maxindex(2)==2) then
        i_sym = maxindex(1)
        b_tetra_irr(1:3,1:4,nbt) = - rot_tetra(bcell, s(:,:,i_sym), b_tetra_irr_l(1:3,1:4,i))
      end if

    enddo !i

      ! Aktualizatu tetra_equiv eta tetra_symlink
      it_loop: do it = 1, n_b_tetra_all

        do i = 1, 4
          t1(1:3,i) = b_tetra_all(1:3,i,it)
        enddo

        do i = 1, n_b_tetra_irr
          do i_sym = 1, nsym
            if ( equal_tetra( t1(1:3,1:4), rot_tetra(bcell, s(:,:,i_sym), b_tetra_irr(1:3,1:4,i)) ) ) then
              tetra_symlink(it,1) = i_sym
              tetra_symlink(it,2) = 0 !
              tetra_equiv(it) = i
              cycle it_loop
            end if

            if ( TR_sym.and.(equal_tetra( - t1(1:3,1:4), rot_tetra(bcell, s(:,:,i_sym), b_tetra_irr(1:3,1:4,i)) )) ) then
              tetra_symlink(it,1) = i_sym
              tetra_symlink(it,2) = 1 !
              tetra_equiv(it) = i
              cycle it_loop
            end if

          end do
        end do
      end do it_loop

  end subroutine compact_tetra


  function overlap_tetra(t1, t2)

    implicit none

    real(dp) :: t1(1:3,1:4)
    real(dp) :: t2(1:3,1:4)

    integer :: overlap_tetra
    integer :: i, j


    overlap_tetra = 0
    do i = 1, 4
      do j = 1, 4
        if (sum(abs(t1(:,i)-t2(:,j)))<1.0E-4_dp) overlap_tetra  = overlap_tetra + 1
      enddo
    enddo

  end function overlap_tetra


  function overlap_tri(t1, t2)

    implicit none

    integer, dimension(1:3), intent(in) :: t1, t2

    integer :: overlap_tri
    integer :: i, j


    if ((size(t1)>3).or.(size(t2)>3)) then
      write(unit=*, fmt=*)"ERROR (overlap_tri). Input are not triangles.", size(t1), size(t2)
      stop
    endif

    overlap_tri = 0
    do i = 1, size(t1)
      do j = 1, size(t2)
        if (t1(i)==t2(j)) then
          overlap_tri = overlap_tri+1
        endif
      enddo
    enddo

  end function overlap_tri


  function overlap_edge(e1, e2)

    implicit none

    integer, intent(in) :: e1(1:2), e2(1:2)

    integer :: overlap_edge
    integer :: i, j


    overlap_edge = 0
    do i = 1, 2
      do j = 1, 2
        if (abs(e1(i)-e2(j))<1.0E-4_dp) overlap_edge = overlap_edge + 1
      enddo
    enddo

  end function overlap_edge


  function overlap_face(f1, f2)

    use intw_matrix_vector, only: norma

    implicit none

    real(dp), intent(in) :: f1(1:3,1:3), f2(1:3,1:3)

    integer :: overlap_face
    integer :: iv, jv


    overlap_face = 0
    do iv = 1, 3
      do jv = 1, 3
        if (norma(f1(:,iv)-f2(:,jv))<1.0E-6_dp) overlap_face = overlap_face + 1
      enddo
    enddo

  end function overlap_face


  function overlap_faces_int(f1, f2)

    implicit none

    integer, intent(in) :: f1(1:3), f2(1:3)

    integer :: overlap_faces_int
    integer :: iv, jv


    overlap_faces_int = 0
    do iv = 1, 3
      do jv = 1, 3
        if (abs(f1(iv)-f2(jv))<1.0E-6_dp) overlap_faces_int  = overlap_faces_int + 1
      enddo
    enddo

  end function overlap_faces_int


  function equal_tetra(t1, t2)

    implicit none

    real(dp) :: t1(1:3,1:4)
    real(dp) :: t2(1:3,1:4)

    logical :: equal_tetra
    integer :: i, j, c


    equal_tetra = .false.

    c = 0
    do i = 1, 4
      do j = 1, 4
        if (sum(abs(t1(:,i)-t2(:,j)))<1.0E-5_dp) c = c + 1
      enddo
    enddo

    if (c==4) then
      equal_tetra = .true.
    endif

  end function equal_tetra


  subroutine tetraIBZ_2_vert_faces_edges(num_tetra, vertx_tetra, verb, nv, nf, ne, vlist, face_as_vert, edge)

    use intw_utility, only: find_free_unit
    use intw_matrix_vector, only: norma

    implicit none

    integer, intent(in) :: num_tetra ! Number of tetrahedra on the IBZ
    real(dp), intent(in) :: vertx_tetra(1:3,1:4,num_tetra) ! Coordinates of the four vertex of each tetrahedra
    logical, intent(in) :: verb ! Output info files or not
    integer, intent(out) :: nv, nf, ne ! Number of vertices, edges and faces
    real(dp), intent(out) :: vlist(1:3,4*num_tetra) ! List of non-repeated vertices
    integer, intent(out) :: face_as_vert(1:3,4*num_tetra) ! List of non-repeated faces (not reduced by S+G symmetry yet)
    integer, intent(out) :: edge(1:2,6*num_tetra) ! List of non-repeated edges

    ! local variables
    integer :: i, j, k, l, e, f, ex, te, tf, fx, fx2
    integer :: io_unit, io_unit2, io_unit3
    integer :: indx(4*num_tetra), edge_aux(1:2,1:6), edge_indx(6*num_tetra), face_aux(1:3,1:4), face_indx(4*num_tetra)


    ! Detect non-repeated vertices
    l = 0
    nv = 0
    do i = 1, num_tetra
      j_loop: do j = 1, 4
        l = l+1 ! counter tetra+vertex
        !write(unit=io_unit, fmt='(3f12.6)') ( vertx_tetra(k,j,i), k = 1, 3 )
        if (nv==0) then
          nv = nv+1
          vlist(:,nv) = vertx_tetra(:,j,i)
          indx(l) = nv ! vertex index
          cycle j_loop
        end if
        do k = 1, nv
          if (norma(vertx_tetra(:,j,i)-vlist(:,k))<10E-6) then ! check if vertex already stored
            indx(l) = k
            cycle j_loop
          end if
        end do
        nv = nv+1
        vlist(:,nv) = vertx_tetra(:,j,i)
        indx(l) = nv
      end do j_loop
    end do

    ! Define edges (non-duplicated) and faces (non-duplicated) of IBZ
    ne = 0
    te = 0 ! counter tetra+edge
    l = 1
    do i = 1, num_tetra
      if (i==1) then ! first tetrahedra, save all the edges.
        edge(1:2,1) = (/ indx(l), indx(l+1) /)
        edge(1:2,2) = (/ indx(l), indx(l+2) /)
        edge(1:2,3) = (/ indx(l), indx(l+3) /)
        edge(1:2,4) = (/ indx(l+1), indx(l+2) /)
        edge(1:2,5) = (/ indx(l+2), indx(l+3) /)
        edge(1:2,6) = (/ indx(l+3), indx(l+1) /)
        ne = ne+6
        !
        do ex = 1, 6
          te = te+1
          edge_indx(te) = te
        end do
      else
        edge_aux(1:2,1) = (/ indx(l), indx(l+1) /)
        edge_aux(1:2,2) = (/ indx(l), indx(l+2) /)
        edge_aux(1:2,3) = (/ indx(l), indx(l+3) /)
        edge_aux(1:2,4) = (/ indx(l+1), indx(l+2) /)
        edge_aux(1:2,5) = (/ indx(l+2), indx(l+3) /)
        edge_aux(1:2,6) = (/ indx(l+3), indx(l+1) /)
        ! Check if any edge has been already stored
        aux_edge_loop: do ex = 1, 6
          te = te+1
          ! Compare edge_aux with all the stored edges
          do e = 1, ne
            if (overlap_edge(edge_aux(:,ex), edge(:,e))==2) then
              edge_indx(te) = e
              cycle aux_edge_loop ! edge_aux already stored, don't save
            end if
          end do
          ne = ne+1
          edge(:,ne) =  edge_aux(:,ex)
          edge_indx(te) = ne ! index of repeated edge list on non-repeated edge list
        end do aux_edge_loop
      end if ! i==1
      l = l+4
    end do ! i

    ! Triangular faces as vertices
    l = 1 ! counter tetra+vert
    tf = 0 ! counter tetra+face
    nf = 1 ! total number of faces
    do i = 1, num_tetra
      if (i==1) then
        face_as_vert(1:3,nf) = (/ indx(l), indx(l+2), indx(l+1)   /)
        face_as_vert(1:3,nf+1) = (/ indx(l), indx(l+3), indx(l+1)   /)
        face_as_vert(1:3,nf+2) = (/ indx(l), indx(l+2), indx(l+3)   /)
        face_as_vert(1:3,nf+3) = (/ indx(l+1), indx(l+2), indx(l+3) /)
        nf = nf+3
        !
        do fx = 1, 4
          tf = tf+1
          face_indx(tf) = tf
        end do
      else
        face_aux(1:3,1) = (/ indx(l), indx(l+2), indx(l+1)   /)
        face_aux(1:3,2) = (/ indx(l), indx(l+3), indx(l+1)   /)
        face_aux(1:3,3) = (/ indx(l), indx(l+2), indx(l+3)   /)
        face_aux(1:3,4) = (/ indx(l+1), indx(l+2), indx(l+3) /)
        ! Check if any edge has been already stored
        aux_face_loop: do fx = 1, 4
          tf = tf+1
          ! Compare edge_aux with all the stored edges
          do f = 1, nf
            if (overlap_faces_int(face_aux(:,fx), face_as_vert(:,f))==3) then
              face_indx(tf) = f
              !cycle aux_face_loop ! face_aux already stored, don't save
              !! Duplicated face! This means face fx and f are inner faces. Remove f.
              !! This will remove inner faces
              do fx2 = f, nf
                face_as_vert(:,fx2) = face_as_vert(:,fx2+1)
              end do
              nf = nf-1
              cycle aux_face_loop
            end if
          end do
          nf = nf+1
          face_as_vert(:,nf) =  face_aux(:,fx)
          face_indx(tf) = nf ! index of repeated face list on non-repeated face list
        end do aux_face_loop
      end if ! i==1
      l = l+4
    end do ! i

    ! Print vertices, edges and faces files if wanted
    if (verb) then
      io_unit = find_free_unit()
      open(unit=io_unit, file='vertices_IBZ.dat', status='unknown')

      write(unit=io_unit, fmt=*) nv
      do i = 1, nv
        write(unit=io_unit, fmt="(3f12.6)") ( vlist(j,i), j = 1, 3 )
      end do

      close(io_unit)

      io_unit2 = find_free_unit()
      open(unit=io_unit2, file='edges_IBZ.dat', status='unknown')

      write(unit=io_unit2, fmt=*) ne
      do e = 1, ne
        write(unit=io_unit2, fmt=*) edge(1:2,e)
      end do

      close(io_unit2)

      io_unit3 = find_free_unit()
      open(unit=io_unit3, file='faces_IBZ.dat', status='unknown')

      write(unit=io_unit3, fmt=*) nf
      do f = 1, nf
        write(unit=io_unit3, fmt=*) face_as_vert(1:3,f)-1
      end do

      close(io_unit3)
    end if

  end subroutine tetraIBZ_2_vert_faces_edges


  subroutine irr_faces(num_tetra, nsym, s, TR_sym, bcell, verb, vlist, nf, face_as_vert, nfaces_irr, face_Gsymlink, face_indx, face_inv_indx)

    use intw_matrix_vector, only: ainv

    implicit none

    ! I/O
    integer, intent(in) :: num_tetra, nsym, nf
    integer, intent(in) :: s(1:3,1:3,nsym)
    logical, intent(in) :: TR_sym
    real(dp), intent(in) :: bcell(1:3,1:3)
    logical, intent(in) :: verb
    real(dp), intent(in) :: vlist(1:3,4*num_tetra) ! List of non-repeated vertices
    integer, intent(in) :: face_as_vert(1:3,4*num_tetra) ! List of non-repeated faces (not reduced by S+G symmetry yet)
    integer, intent(out) :: nfaces_irr ! Number of irreducible faces
    integer, intent(out) :: face_Gsymlink(4*num_tetra,2,-1:1,-1:1,-1:1), face_indx(4*num_tetra), face_inv_indx(4*num_tetra)
    real(dp) :: face_aux(1:3,1:3)
    real(dp) :: bcelli(1:3,1:3), vert_crys(1:3), rot_face(1:3,1:3)

    ! Local
    integer :: iv, fx, f, isym, ig, jg, kg


    ! Detect irreducible faces considering S+G symmetry operations

    bcelli = ainv(bcell)

    face_Gsymlink = 0
    nfaces_irr = 1
    face_indx(1) = 1
    face_inv_indx(1) = 1
    aux_face_loop : do fx = 2, nf
      do iv = 1, 3
        face_aux(1:3,iv) = vlist(:,face_as_vert(iv,fx))
      end do
      ! Compare face with the already stored faces + all the possible S+G rotations
      do f = 1, nfaces_irr
        do isym = 1, nsym
          do ig = -1, 1
            do jg = -1, 1
              do kg = -1, 1
                if ( ig==0 .and. jg==0 .and. kg==0 ) cycle
                do iv = 1, 3
                  vert_crys = matmul(bcelli, vlist(:,face_as_vert(iv,face_inv_indx(f))))
                  rot_face(1:3,iv) = matmul(bcell, matmul(dble(s(:,:,isym)), vert_crys)) & ! rotate with S
                                      + ig*bcell(1:3,1) + jg*bcell(1:3,2) + kg*bcell(1:3,3) ! add G
                end do
                ! check with face_aux
                if ( overlap_face(face_aux, rot_face) == 3 ) then
                  ! Face is related to stored face f by symm.operation isym, no TR, adding (ig,jg,kg) G vector
                  if (verb) write(*, '("Face",I6," is related to face",I6," by G:",3I6)') fx, face_inv_indx(f), ig, jg, kg
                  face_indx(fx) = f
                  face_Gsymlink(fx,1,ig,jg,kg) = isym
                  face_Gsymlink(fx,2,ig,jg,kg) = 0 ! no TR
                  cycle aux_face_loop
                end if
                ! TR symmetry
                if ( TR_sym )  then
                  do iv = 1, 3
                    vert_crys = matmul(bcelli, vlist(:,face_as_vert(iv,face_inv_indx(f))))
                    rot_face(1:3,iv) = - matmul(bcell, matmul(dble(s(:,:,isym)), vert_crys)) & ! rotate with S + TR
                                        + ig*bcell(1:3,1) + jg*bcell(1:3,2) + kg*bcell(1:3,3) ! add G
                  end do
                  if (overlap_face(face_aux, rot_face)==3) then
                    ! Face is related to stored f face by isym symm.operation + TR, plus (ig,jg,kg) G vector
                    if (verb) write(*, '("Face",I6," is related to face",I6," with TR and by G:",3I6)') fx, face_inv_indx(f), ig, jg, kg
                    face_indx(fx) = f
                    face_Gsymlink(fx,1,ig,jg,kg) = isym
                    face_Gsymlink(fx,2,ig,jg,kg) = 1 !
                    cycle aux_face_loop
                  end if
                end if
                !
              end do !kg
            end do !jg
          end do !ig
        end do !isym
      end do !f
      ! Face not related by sym with any stored face
      nfaces_irr = nfaces_irr + 1
      face_indx(fx) = nfaces_irr
      face_inv_indx(nfaces_irr) = fx  ! Relates new stored face with index in old list
    end do aux_face_loop

    if (verb) write(*, *)
    if (verb) write(*, '("Number of faces on IBZ:",I6)') nf
    if (verb) write(*, '("Number of irreducible faces by S+G:",I6)') nfaces_irr
    if (verb) write(*, *)

  end subroutine irr_faces


  subroutine sort_edges(nen, nemx, nein, nemn)

    implicit none

    ! I/O
    integer, dimension(1:3), intent(in) :: nen
    integer, intent(out) :: nemx, nein, nemn

    ! Local
    integer :: i


    ! special case when all edges have same number of nodes
    if ( nen(1)==nen(2) .and. nen(2)==nen(3) ) then
      nemx = 1
      nein = 2
      nemn = 3
      return
    end if

    nemx = maxloc(nen(:), dim=1)
    nemn = minloc(nen(:), dim=1)
    nein = 1
    do i = 1, 3
      if ( i/=nemn .and. i/=nemx ) nein = i
    end do
    if ( nemx==nemn .or. nemx==nein .or. nemn==nein ) then
      write(*, '("ERROR on sorting edges, stopping")')
      write(*, '("ne1, ne2, ne3 :",3I6)') nen(1), nen(2), nen(3)
      write(*, '("nemax, neint, nemin :",3I6)') nemx, nein, nemn
      stop
    end if

  end subroutine sort_edges


  subroutine write_node_list_face(f, tot_nodes, node_coords)

    use intw_utility, only: find_free_unit, int2str

    implicit none

    ! I/O
    integer, intent(in) :: f, tot_nodes
    real(dp), dimension(:,:), intent(in) :: node_coords

    ! Local
    integer :: io_unit, i


    io_unit = find_free_unit()
    open(unit=io_unit, file=trim(int2str(f))//"nodes_on_face.node", status="unknown")
    write(unit=io_unit, fmt='(I12)') tot_nodes
    do i = 1, tot_nodes
      write(io_unit, fmt='(I6, 3E18.10)') i, node_coords(1:3,i)
    end do
    close(unit=io_unit)

  end subroutine write_node_list_face


  subroutine nodes_on_face_and_rotate_to_plane(f, nk1, nk2, nk3, bg, ntetra, face_as_vert, vert_coord, verb, split_edges_nodes, k, theta)

    use intw_utility, only: find_free_unit
    use intw_matrix_vector, only: norma, cross

    implicit none

    ! I/O
    integer, intent(in) :: f ! face number (only needed for writing node file)
    integer, intent(in) :: nk1, nk2, nk3 ! BZ partition from input
    real(dp), intent(in) :: bg(1:3,1:3) ! reciprocal lattice vectors as columns
    integer, intent(in) :: ntetra
    integer, intent(in) :: face_as_vert(1:3) ! list with vertices of face
    real(dp), intent(in) :: vert_coord(1:3,4*ntetra) ! list of vertices coordinates
    logical, intent(in) :: verb ! output edge nodes or not
    real(dp), intent(out) :: split_edges_nodes(1:3,nk1+nk2+nk3,1:3) ! xzy coords. of nodes of split edges
    real(dp), intent(out) :: k(1:3) ! rotation axis from actual 3d to xy plane
    real(dp), intent(out) :: theta ! rotation angle

    ! Local
    integer :: ne(3) ! number of splitted points for each edge
    integer :: tot_nodes_face ! total number of nodes on face
    real(dp) :: node_list(1:3,3*(nk1+nk2+nk3)) ! list with total node coordinates of face
    real(dp) :: nodes_in_plane(1:2,3*(nk1+nk2+nk3)) ! 2d coords of nodes displaced and rotated to xy plane
    real(dp) :: edges_vector(1:3,1:3) ! 3d vectors of edges of face
    real(dp), dimension(1:3,1:3) :: vf ! Coordinates of vertices of face
    real(dp), dimension(1:3) :: normal, xy_plane ! Normal vector of face, normal of xy plane, and rotation axis
    real(dp) :: node_aux(1:3), dist, dist_aux, vec_extra(1:3)
    real(dp), allocatable :: rotated_nodes(:,:) ! Rotated and displaced nodes
    integer :: io_unit, i, j, l, nemax, neint, nemin, ne_extra


    !!------------------------ Split edges in nodes

    edges_vector(1:3,1) = vert_coord(1:3,face_as_vert(2)) - vert_coord(1:3,face_as_vert(1))
    edges_vector(1:3,2) = vert_coord(1:3,face_as_vert(3)) - vert_coord(1:3,face_as_vert(2))
    edges_vector(1:3,3) = vert_coord(1:3,face_as_vert(1)) - vert_coord(1:3,face_as_vert(3))

    ne(1) = NINT( sqrt(nk1 * abs(dot_product(edges_vector(1:3,1), bg(1:3,1))) + &
                       nk2 * abs(dot_product(edges_vector(1:3,1), bg(1:3,2))) + &
                       nk3 * abs(dot_product(edges_vector(1:3,1), bg(1:3,3)))) )
    ne(2) = NINT( sqrt(nk1 * abs(dot_product(edges_vector(1:3,2), bg(1:3,1))) + &
                       nk2 * abs(dot_product(edges_vector(1:3,2), bg(1:3,2))) + &
                       nk3 * abs(dot_product(edges_vector(1:3,2), bg(1:3,3)))) )
    ne(3) = NINT( sqrt(nk1 * abs(dot_product(edges_vector(1:3,3), bg(1:3,1))) + &
                       nk2 * abs(dot_product(edges_vector(1:3,3), bg(1:3,2))) + &
                       nk3 * abs(dot_product(edges_vector(1:3,3), bg(1:3,3)))) )

    ! Convention: last node will be stored as first of the next edge (not to repeat nodes)
    tot_nodes_face = 0
    do i = 1, ne(1)
      split_edges_nodes(1:3,i,1) = vert_coord(1:3,face_as_vert(1)) + edges_vector(1:3,1)*(i-1)/ne(1)
      tot_nodes_face = tot_nodes_face + 1
      node_list(1:3,tot_nodes_face) = split_edges_nodes(1:3,i,1)
    end do
    do i = 1, ne(2)
      split_edges_nodes(1:3,i,2) = vert_coord(1:3,face_as_vert(2)) + edges_vector(1:3,2)*(i-1)/ne(2)
      tot_nodes_face = tot_nodes_face + 1
      node_list(1:3,tot_nodes_face) = split_edges_nodes(1:3,i,2)
    end do
    do i = 1, ne(3)
      split_edges_nodes(1:3,i,3) = vert_coord(1:3,face_as_vert(3)) + edges_vector(1:3,3)*(i-1)/ne(3)
      tot_nodes_face = tot_nodes_face + 1
      node_list(1:3,tot_nodes_face) = split_edges_nodes(1:3,i,3)
    end do


    !!------------------------ Fill face with more nodes

    ! Sort edges with number of nodes
    call sort_edges(ne, nemax, neint, nemin)

    ! Join nodes on most populated edge with nodes of intermediate edge
    if (ne(nemax)>1) then
      do i = 2, ne(nemax) ! first node is always shared by two edges
        !
        !! detect closest node on intermediate edge and define vector in between
        !dist_aux = huge(dist)
        !do j = 2, ne(neint)
        !  dist = norma(split_edges_nodes(1:3,j,neint) - split_edges_nodes(1:3,i,nemax))
        !  if (dist<dist_aux) then
        !    vec_extra(1:3) = split_edges_nodes(1:3,j,neint) - split_edges_nodes(1:3,i,nemax)
        !    dist_aux = dist
        !  end if
        !end do ! j
        ! assign pair-node on intermediate edge "by hand"
        if ( i <= ne(neint) ) then
          vec_extra(1:3) = split_edges_nodes(1:3,ne(neint)-(i-2),neint) - split_edges_nodes(1:3,i,nemax)
        !else if ( i == ne(neint)+1 ) then
        !  vec_extra(1:3) = split_edges_nodes(1:3,1,neint) - split_edges_nodes(1:3,i,nemax)
        else
          cycle
        end if
        ! add extra nodes on line connecting nodes on edges
        ne_extra = NINT( sqrt(nk1 * abs(dot_product(vec_extra(1:3), bg(1:3,1))) + &
                              nk2 * abs(dot_product(vec_extra(1:3), bg(1:3,2))) + &
                              nk3 * abs(dot_product(vec_extra(1:3), bg(1:3,3)))) )
        if (ne_extra>1) then
          do l = 2, ne_extra
            node_list(1:3,tot_nodes_face+l-1) = split_edges_nodes(1:3,i,nemax) + vec_extra(1:3)*(l-1)/ne_extra
          end do
          tot_nodes_face = tot_nodes_face + ne_extra - 1
        end if ! ne_exgtra>1
        !
      end do ! i
    end if

    ! Write nodes on face to file if needed
    if (verb) call write_node_list_face(f, tot_nodes_face, node_list)


    !!------------------------ Rotate nodes to lie at (x,y) plane and displace w/r first node
    !
    xy_plane = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
    ! Vertices of face
    vf(1:3,1) = vert_coord(1:3,face_as_vert(1))
    vf(1:3,2) = vert_coord(1:3,face_as_vert(2))
    vf(1:3,3) = vert_coord(1:3,face_as_vert(3))
    ! Normal vector of face
    normal = normal_of_tri(vf)
    if ( norma(cross(normal, xy_plane)) < 1.0E-6_dp ) then ! plane of face is already xy plane
      write(*, *) "Face is on xy plane"
      write(*, *) "Vectors:", vf(1:3,1)
      write(*, *) "        ", vf(1:3,2)
      write(*, *) "        ", vf(1:3,3)
      write(*, *) "Normal:", normal
      k = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
      theta = 0.0_dp
    else ! define rotation vector and angle
      k = cross(normal, xy_plane) / norma( cross(normal, xy_plane) )
      theta = acos( dot_product(xy_plane, normal) )
    end if
    !
    ! Rodrigues' rotation formula for each node
    allocate(rotated_nodes(1:3,tot_nodes_face))
    do i = 1, tot_nodes_face
      ! displace to origin and rotate to xy plane
      node_aux(1:3) = node_list(1:3,i) - split_edges_nodes(1:3,1,1)
      rotated_nodes(1:3,i) = node_aux(1:3)*cos(theta) + cross(k, node_aux)*sin(theta) + k(1:3)*dot_product(k, node_aux)*(1.d0-cos(theta))
      ! check
      if (abs(rotated_nodes(3,i))>1.0E-10_dp) then
        write(*, *) "Rotation of face failed, something wrong..."
        write(*, *) "z-component of node i:", i
        write(*, *) "Before rotating      :", node_list(1:3,i)
        write(*, *) "After rotating       :", rotated_nodes(1:3,i)
        write(*, *) "Origin node          :", split_edges_nodes(1:3,1,1)
        write(*, *) "EXIT"
        deallocate(rotated_nodes)
        stop
      end if
      nodes_in_plane(1:2,i) = rotated_nodes(1:2,i)
    end do ! i
    deallocate(rotated_nodes)


    !!------------------------ Write node list in poly format to triangulate

    ! Write .poly file to be read by Triangle
    io_unit = find_free_unit()
    open(unit=io_unit, file="face_split_edges.poly", status="unknown")
    write(io_unit, fmt='(4I6)') tot_nodes_face, 2, 0, 0
    do j = 1, tot_nodes_face
      write(io_unit, fmt='(I6, 2F18.10)') j, nodes_in_plane(1:2,j)
    end do
    write(io_unit, fmt='(2I6)') sum(ne(:)), 0
    ! add segments on border of face
    do j = 1, sum(ne(:)) - 1
      write(io_unit, fmt='(3I6)') j, j, j+1
    end do
    write(io_unit, fmt='(3I6)') sum(ne(:)), sum(ne(:)), 1
    write(io_unit, fmt='(I6)') 0
    close(io_unit)

  end subroutine nodes_on_face_and_rotate_to_plane


  subroutine triangulate_faces(ntetra, nfaces, faces_Gsymlink, nsym, s, faces_indx, faces_inv_indx, nk1, nk2, nk3, bcell, &
                               faces_as_vert, vert_coord, verbose)

    use intw_utility, only: find_free_unit, int2str
    use intw_matrix_vector, only: ainv, norma, cross

    implicit none

    ! I/O
    integer, intent(in) :: ntetra, nfaces ! number of tetrahedra and faces on IBZ
    integer, intent(in) :: faces_Gsymlink(4*ntetra,2,-1:1,-1:1,-1:1) ! info about S+G symm of faces
    integer, intent(in) :: nsym ! number of symmetry operations
    integer, intent(in) :: s(1:3,1:3,nsym) ! symmetry operation matrices
    integer, intent(in) :: faces_indx(4*ntetra), faces_inv_indx(4*ntetra) ! faces indexes from full to reduced list
    integer, intent(in) :: nk1, nk2, nk3 ! partition chosen to perform edge-splitting
    real(dp), intent(in) :: bcell(1:3,1:3) ! reciprocal cell vectors on columns
    integer, intent(in) :: faces_as_vert(1:3,4*ntetra) ! vertices indices on faces
    real(dp), intent(in) :: vert_coord(1:3,4*ntetra) ! Coordinates of vertices on IBZ
    logical, intent(in) :: verbose ! logical controlling output verbosity

    ! Local
    integer :: isym, f, iface_irr, iface_full
    integer :: i, j, k, l, ig, jg, kg, iv, nv, it, nt
    integer :: io_unit, io_unit2, dummy_i
    real(dp) :: kvec(1:3), theta, v_aux(1:3)
    logical, allocatable :: triangulated_face(:)
    integer, allocatable :: nvert_face(:), ntri_face(:), triangles_face(:,:,:), indx(:), tri_face_indx(:,:)
    real(dp), allocatable :: split_edges_nodes(:,:,:), nodes_fplane(:,:), vertices_face(:,:,:), vlist(:,:)


    allocate(nvert_face(nfaces), ntri_face(nfaces), triangulated_face(nfaces))
    allocate(vertices_face(1:3,3*(nk1+nk2+nk3),nfaces), triangles_face(1:3,3*(nk1+nk2+nk3),nfaces))
    allocate(split_edges_nodes(1:3,nk1+nk2+nk3,1:3))
    allocate(nodes_fplane(1:2,3*(nk1+nk2+nk3)))
    triangulated_face(:) = .false.
    do f = 1, nfaces
      !
      if (ANY(faces_Gsymlink(f,:,:,:,:)/=0)) then
        !
        if (verbose) write(*, '("Face", I6," will be obtained by symmetry from", I6)') f, faces_indx(f)
        if (verbose) write(*, '("First triangulate ", I6)') faces_indx(f)
        !
        iface_irr = faces_indx(f)
        iface_full = faces_inv_indx(iface_irr)
        call nodes_on_face_and_rotate_to_plane(iface_full, nk1, nk2, nk3, bcell, ntetra, faces_as_vert(1:3,iface_full), &
                                               vert_coord, verbose, split_edges_nodes, kvec, theta)
        !
        CALL EXECUTE_COMMAND_LINE("triangle -g face_split_edges.poly")
        !
        ! Read output face triangulation file, rotate vertices and write to file
        !
        if (verbose) then
          io_unit2 = find_free_unit()
          open(unit=io_unit2, file=trim(int2str(iface_full))//"_triangle_output.off", status='unknown')
        end if
        !
        io_unit = find_free_unit()
        open(unit=io_unit, file='face_split_edges.1.off', status='unknown')
        read(unit=io_unit, fmt=*)
        read(unit=io_unit, fmt=*) nvert_face(iface_full), ntri_face(iface_full)
        if (verbose) write(unit=io_unit2, fmt=*) "OFF"
        if (verbose) write(unit=io_unit2, fmt=*) nvert_face(iface_full), ntri_face(iface_full)
        do iv = 1, nvert_face(iface_full)
          read(unit=io_unit, fmt=*) v_aux(1:3) ! read
          vertices_face(1:3,iv,iface_full) = v_aux(1:3)*cos(-theta) + cross(kvec, v_aux)*sin(-theta) + kvec(1:3)*dot_product(kvec, v_aux)*(1.d0-cos(-theta)) & ! rotate to true 3D plane
                                             + split_edges_nodes(1:3, 1, 1) ! displace to correct origin
          if (verbose) write(unit=io_unit2, fmt='(3F12.6)') vertices_face(1:3,iv,iface_full)
        end do
        do it = 1, ntri_face(iface_full)
          read(unit=io_unit, fmt=*) dummy_i, triangles_face(1:3,it,iface_full)
          if (verbose) write(unit=io_unit2, fmt='(4I6)') dummy_i, triangles_face(1:3,it,iface_full)
        end do
        triangles_face(:,:,iface_full) = triangles_face(:,:,iface_full) + 1 ! start indexes from 1
        close(unit=io_unit)
        if (verbose) close(unit=io_unit2)
        !
        ! Rotate triangulation to faces related by symmetry
        do ig = -1, 1
          do jg = -1, 1
            do kg = -1, 1
              isym = faces_Gsymlink(f,1,ig,jg,kg)
              if (isym==0) cycle
              if (faces_Gsymlink(f,2,ig,jg,kg)==0) then ! no TR
                do iv = 1, nvert_face(iface_full)
                  v_aux = matmul(ainv(bcell), vertices_face(1:3,iv,iface_full))
                  vertices_face(1:3,iv,f) = matmul(bcell, matmul(dble(s(:,:,isym)), v_aux)) & ! rotate with S
                                            + ig*bcell(1:3,1) + jg*bcell(1:3,2) + kg*bcell(1:3,3) ! add G
                end do !iv
              else ! TR
                do iv = 1, nvert_face(iface_full)
                  v_aux = matmul(ainv(bcell), vertices_face(1:3,iv,iface_full))
                  vertices_face(1:3,iv,f) = -matmul(bcell, matmul(dble(s(:,:,isym)), v_aux)) & ! rotate with S
                                            + ig*bcell(1:3,1) + jg*bcell(1:3,2) + kg*bcell(1:3,3) ! add G
                end do !iv
              end if
            end do !kg
          end do !jg
        end do !ig
        ! Equal number and indexes of triangles
        nvert_face(f) = nvert_face(iface_full)
        ntri_face(f) = ntri_face(iface_full)
        triangles_face(:,:,f) = triangles_face(:,:,iface_full)
        ! Set faces as triangulated
        triangulated_face(f) = .true.
        triangulated_face(iface_full) = .true.
        !
      end if ! ANY(faces_Gsymlink)
    end do
    do f = 1, nfaces
      !
      if (triangulated_face(f)) cycle
      !
      if (verbose) write(*, *)
      if (verbose) write(*, '("Face", I6," will be directly triangulated")') f
      !
      call nodes_on_face_and_rotate_to_plane(f, nk1, nk2, nk3, bcell, ntetra, faces_as_vert(1:3,f), vert_coord, verbose, &
                                             split_edges_nodes, kvec, theta)
      !
      CALL EXECUTE_COMMAND_LINE("triangle -pg face_split_edges.poly")
      !
      ! Read output face triangulation file, rotate vertices and write to file
      !
      if (verbose) then
        io_unit2 = find_free_unit()
        open(unit=io_unit2, file=trim(int2str(f))//"_triangle_output.off", status='unknown')
      end if
      !
      io_unit = find_free_unit()
      open(unit=io_unit, file='face_split_edges.1.off', status='unknown')
      read(unit=io_unit, fmt=*)
      read(unit=io_unit, fmt=*) nvert_face(f), ntri_face(f)
      if (verbose) write(unit=io_unit2, fmt=*) "OFF"
      if (verbose) write(unit=io_unit2, fmt=*) nvert_face(f), ntri_face(f)
      do iv = 1, nvert_face(f)
        read(unit=io_unit, fmt=*) v_aux(1:3) ! read
        vertices_face(1:3,iv,f) = v_aux(1:3)*cos(-theta) + cross(kvec, v_aux)*sin(-theta) + kvec(1:3)*dot_product(kvec, v_aux)*(1.d0-cos(-theta)) & ! rotate to true 3D plane
                                  + split_edges_nodes(1:3,1,1) ! displace to correct origin
        if (verbose) write(unit=io_unit2, fmt='(3F12.6)') vertices_face(1:3,iv,f)
      end do
      do it = 1, ntri_face(f)
        read(unit=io_unit, fmt=*) dummy_i, triangles_face(1:3,it,f)
        if (verbose) write(unit=io_unit2, fmt='(4I6)') dummy_i, triangles_face(1:3,it,f)
      end do
      triangles_face(:,:,f) = triangles_face(:,:,f) + 1 ! start indexes from 1
      close(unit=io_unit)
      if (verbose) close(unit=io_unit2)
      !
      triangulated_face(f) = .true.
      !
    end do !f
    deallocate(split_edges_nodes, nodes_fplane)
    !

    ! Create a single OFF file with all the triangulated faces cleaning repeated nodes
    allocate(vlist(1:3,sum(nvert_face(:))), indx(maxval(nvert_face(:))), tri_face_indx(1:3,sum(ntri_face(:))))
    nv = 0
    nt = 0
    do f = 1, nfaces
      l = 0
      !! Do not plot rotated triangulation
      !if (ANY(faces_Gsymlink(f,:,:,:,:)/=0)) cycle
      !! Plot only rotated face
      !if (ALL(faces_Gsymlink(f,:,:,:,:)==0)) cycle
      !
      iv_loop:do iv = 1, nvert_face(f)
        l = l+1
        if (nv==0) then
          nv = nv+1
          vlist(:,nv) = vertices_face(1:3,iv,f)
          indx(l) = nv ! vertex index
          cycle iv_loop
        end if
        do k = 1, nv
          if (norma(vertices_face(1:3,iv,f)-vlist(1:3,k))<1.0E-5_dp) then ! check if vertex already stored
            indx(l) = k
            cycle iv_loop
          end if
        end do
        nv = nv+1
        vlist(:,nv) = vertices_face(:,iv,f)
        indx(l) = nv
      end do iv_loop
      !
      it_loop:do it = 1, ntri_face(f)
        nt = nt+1
        tri_face_indx(1:3,nt) = (/ indx(triangles_face(1,it,f)), indx(triangles_face(2,it,f)), indx(triangles_face(3,it,f)) /)
      end do it_loop
      !
      !!end if
    end do ! nfaces
    ! Write to OFF file
    io_unit = find_free_unit()
    open(unit=io_unit, file="Triangulated_IBZ.off", status="unknown")
    write(unit=io_unit, fmt="(a)") "OFF"
    write(unit=io_unit, fmt="(3I6)") nv, nt, 0
    write(unit=io_unit, fmt=*)
    ! Vertices of triangles
    do i = 1, nv
      write(unit=io_unit, fmt="(3f12.6)") ( vlist(j,i), j = 1, 3 )
    end do
    ! Triangle indexes
    do i = 1, nt
      write(unit=io_unit, fmt="(4I6)") 3, tri_face_indx(1:3,i) - 1
    end do
    close(unit=io_unit)

    ! Write nodes on faces to file if needed
    if (verbose) then
      io_unit = find_free_unit()
      open(unit=io_unit, file="nodes_on_faces.node", status="unknown")
      write(unit=io_unit, fmt="(3I6)") nv
      write(unit=io_unit, fmt=*)
      ! Vertices of triangles
      do i = 1, nv
        write(unit=io_unit, fmt="(I6, 3f12.6)") i, ( vlist(j,i), j = 1, 3 )
      end do
    end if

    close(io_unit)

    deallocate(vlist, indx, tri_face_indx)
    deallocate(nvert_face, ntri_face, triangulated_face, vertices_face, triangles_face)

  end subroutine triangulate_faces


  subroutine add_nodes_IBZ_volume(nk1, nk2, nk3, m, epsface, bg, n_BZ_tetra_irr, BZ_tetra_irr)
    ! Create regular mesh and check which points lie inside IBZ

    use intw_utility, only: find_free_unit
    use intw_matrix_vector, only: ainv, norma, cross

    implicit none

    ! I/O
    integer, intent(in) :: nk1, nk2, nk3, n_BZ_tetra_irr
    real(dp), intent(in) :: m, epsface
    real(dp), intent(in) :: BZ_tetra_irr(1:3,1:4,n_BZ_tetra_irr)
    real(dp), intent(in) :: bg(1:3,1:3)

    ! Local
    integer :: io_unit, i, j, k, l, it, i_min, i_max, j_min, j_max, k_min, k_max, nk_in
    real(dp) :: vcrys(1:3), vcart(1:3), mat(1:3,1:3), coord(1:3), face(1:3,1:3), coef(1:3)
    real(dp), allocatable :: node_coords(:,:)
    ! real(dp), parameter :: epsvert = 1.0E-6_dp
    ! integer, parameter :: m = 1
    ! real(dp), parameter :: m = 1.5_dp


    i_min = -int(m*nk1)
    i_max =  int(m*nk1)
    j_min = -int(m*nk2)
    j_max =  int(m*nk2)
    k_min = -int(m*nk3)
    k_max =  int(m*nk3)

    write(*, "(100i4)") nk1, nk2, nk3
    write(*, "(100i4)") i_min, i_max, j_min, j_max, k_min, k_max

    allocate(node_coords(3,n_BZ_tetra_irr*(i_max-i_min)*(j_max-j_min)*(k_max-k_min)))

    nk_in = 0
    ! do i = i_min, i_max
    !   do j = j_min, j_max
    !     k_loop: do k = k_min, k_max
    do i = 0, 2*int(m*nk1)
      do j = 0, 2*int(m*nk2)
        k_loop: do k = 0, 2*int(m*nk3)
          !
          !vcrys = (/real(i-1, dp)/(2.0_dp*nk1), real(j-1, dp)/(2.0_dp*nk2), real(k-1, dp)/(2.0_dp*nk3)/)
          !vcrys = (/real(i, dp)/(int(m)*nk1), real(j, dp)/(int(m)*nk2), real(k, dp)/(int(m)*nk3)/)
          vcrys = (/ -1.0_dp + real(i, dp)/(int(m*nk1)), -1.0_dp + real(j, dp)/(int(m*nk2)), -1.0_dp + real(k, dp)/(int(m*nk3)) /)
          vcart = matmul(bg, vcrys)
          !
          do it = 1, n_BZ_tetra_irr
            mat = BZ_tetra_irr(1:3,2:4,it)
            coord = matmul(ainv(mat), vcart)
            !
            ! Check if k-point is inside IBZ
            if ( (-epsface<=coord(1)) .and. (-epsface<=coord(2)) .and. (-epsface<=coord(3)) &
                 .and. (coord(1)-1.d0<=epsface) .and. (coord(2)-1.d0<=epsface) &
                 .and. (coord(3)-1.d0<=epsface).and.(sum(coord)-1.d0<=epsface) ) then

                ! Check if k-point is on border of BZ
                if (ANY(coord(:)-epsface<=0.0_dp) .or. abs(sum(coord(:))-1.0_dp)<=epsface) then
                  ! Node lies on BZ border face, discard
                  cycle
                end if

                ! Check if k-point is on any other face of tetrahedra
                do l = 2, 3
                  face(:,1) = BZ_tetra_irr(:,l,it) - BZ_tetra_irr(:,1,it) ! vface(:,1,j)
                  face(:,2) = BZ_tetra_irr(:,l+1,it) - BZ_tetra_irr(:,1,it)
                  face(:,3) = cross(face(:,2), face(:,1))
                  face(:,3) = face(:,3) / norma(face(:,3))
                  coef(1:3) = matmul(ainv(face(1:3,1:3)), coord(1:3))
                  ! Coordinate vector of the vertex is perpendicular to normal vector of face (face and vector coplanar)
                  if ( abs(dot_product(face(:,3), coord(:))) < epsface ) then
                    ! Coordinate vector of the vertex lies within the face
                    if ( coef(1)>=0.0_dp-epsface .and. coef(1)<=1.0_dp+epsface &
                         .and. coef(2)>=0.0_dp-epsface .and. coef(2)<=1.0_dp+epsface & !) then
                         .and. sum(coef(1:2))-1.0_dp<epsface ) then
                      !
                      ! Node lies on IBZ face, discard
                      cycle
                      !
                    end if
                  end if
                end do ! l

                ! Node lies within a IBZ tetra, and not on its faces, so save
                nk_in = nk_in + 1
                node_coords(1:3,nk_in) = vcart(1:3)

                cycle k_loop

            end if ! Check if k-point is inside IBZ
          end do ! it
          !
        end do k_loop ! k
      end do ! j
    end do ! i

    print*, "Number of points on IBZ:", nk_in

    ! Write nodes to .node file
    io_unit = find_free_unit()
    open(unit=io_unit, file="Triangulated_IBZ.a.node", status="unknown")
    write(unit=io_unit, fmt='(I12)') nk_in
    do i = 1, nk_in
      write(io_unit, fmt='(I6, 3E18.10)') i, node_coords(1:3,i)
    end do
    close(unit=io_unit)

    deallocate(node_coords)

  end subroutine add_nodes_IBZ_volume


end module triFS_geometry
