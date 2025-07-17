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
module siesta2intw_fft

  use precision, only: dp

  implicit none

  public :: nG, nGk_max, list_iG, ngk, gvec, gvec_cart, gamma_only, nl, nlm

  public :: compute_fft_info, write_fft_info

  private

  !
  integer :: nG ! number of G vectors inside the rho cutoff (g2cut)
  integer :: nGk_max ! max. number of G vectors inside the wave function cutoff (cutoff)
  integer, allocatable, dimension(:,:) :: list_iG ! list of G vectors for each k point
  integer, allocatable, dimension(:) :: ngk ! number of G vectors inside the wave function cutoff for each k point
  integer, allocatable, dimension(:,:) :: gvec ! G vectors in crytal units
  real(kind=dp), allocatable, dimension(:,:) :: gvec_cart ! G vectors in cartesian units
  logical :: gamma_only ! whether only the Gamma-point is sampled
  integer, allocatable :: nl(:) ! relation between FFT grid and G vectors
  integer, allocatable :: nlm(:) ! relation between FFT grid and -G vectors for gamma_only calculations


contains


  subroutine compute_fft_info(tpiba, bg, cutoff)
    !
    ! Compute all FFT stuff: G vectors, list_iG...
    !

    ! variables
    use units, only : pi
    use writewave, only: nwk, wfk
    use sorting, only: ordix
    use m_ntm, only: ntm
    use siesta_options, only: g2cut
    use siesta2intw_io, only: stdout
    ! functions and subroutines
    use siesta2intw_io, only: find_free_unit
    use siesta2intw_utils, only: ainv

    implicit none

    real(kind=dp), intent(in) :: tpiba, bg(3,3), cutoff

    integer :: ik
    real(kind=dp) :: klist_cryst(3,nwk)
    integer :: n1, n2, n3
    integer :: ig, jg, kg
    real(dp) :: gv(3)
    real(dp), allocatable :: gg_(:), gvec_cart_(:,:)
    integer, allocatable :: gvec_(:,:), perm(:)
    logical :: assigned(ntm(1), ntm(2), ntm(3))


    write(stdout, *) "- Computing FFT info..."

    !
    n1 = (ntm(1)-1)/2
    n2 = (ntm(2)-1)/2
    n3 = (ntm(3)-1)/2

    ! count number of G vectors (nG) inside the rho cutoff
    nG = 0
    iloop1: do ig=-n1,n1
      !
      if ( gamma_only .and. ig < 0) cycle iloop1
      !
      jloop1: do jg=-n2,n2
        !
        if ( gamma_only .and. ig == 0 .and. jg < 0) cycle jloop1
        !
        kloop1: do kg=-n3,n3
          !
          if ( gamma_only .and. ig == 0 .and. jg == 0 .and. kg < 0) cycle kloop1
          !
          gv = tpiba * matmul(bg, (/ig, jg, kg/))
          if (sum(gv**2) <= g2cut) nG = nG + 1
          !
        enddo kloop1 ! kg
        !
      enddo jloop1 ! jg
      !
    enddo iloop1 ! ig

    ! compute G vectors inside rho cutoff
    allocate(gvec_(3,nG), gvec_cart_(3,nG), gg_(nG), perm(nG))
    gvec_ = 0
    gvec_cart_ = 0.0_dp
    gg_ = 0.0_dp
    perm = 0

    nG=0
    iloop2: do ig=-n1,n1
      !
      if ( gamma_only .and. ig < 0) cycle iloop2
      !
      jloop2: do jg=-n2,n2
        !
        if ( gamma_only .and. ig == 0 .and. jg < 0) cycle jloop2
        !
        kloop2: do kg=-n3,n3
          !
          if ( gamma_only .and. ig == 0 .and. jg == 0 .and. kg < 0) cycle kloop2
          !
          gv = tpiba * matmul(bg, (/ig, jg, kg/))
          !
          if (sum(gv**2) <= g2cut) then
            nG = nG + 1
            gg_(nG) = sum(gv**2)
            gvec_cart_(:,nG) = gv
            gvec_(:,nG) = (/ig, jg, kg/)
          end if
          !
        enddo kloop2 ! kg
        !
      enddo jloop2 ! jg
      !
    enddo iloop2 ! ig

    ! sort G vectors by increasing module
    allocate(gvec(3,nG), gvec_cart(3,nG))
    gvec = 0
    gvec_cart = 0.0_dp
    !
    call ordix( gg_, 1, nG, perm)
    !
    do ig=1,nG
      gvec(:,ig) = gvec_(:,perm(ig))
      gvec_cart(:,ig) = gvec_cart_(:,perm(ig))
    enddo
    !
    deallocate(gg_, gvec_cart_, gvec_, perm)

    ! generate nl and nlm
    assigned = .false.
    allocate(nl(nG))
    nl = 0
    if (gamma_only) then
      allocate(nlm(nG))
      nlm = 0
    endif

    do ig = 1, nG

      n1 = modulo( gvec(1, ig), ntm(1) ) + 1
      n2 = modulo( gvec(2, ig), ntm(2) ) + 1
      n3 = modulo( gvec(3, ig), ntm(3) ) + 1

      if ( .not. assigned(n1, n2, n3) ) then

        assigned(n1, n2, n3) = .true.

        ! compute the scalar index corresponding to n1,n2,n3 and
        ! assign it to nl(ig)

        nl(ig) = n1 + (n2-1)*ntm(1) + (n3-1)*ntm(1)*ntm(2)

      else
        write(stdout,*) ' ERROR in generate_nl. FFT mesh too small?            '
        write(stdout,*) '    More than one G-vector in the gvec array are being'
        write(stdout,*) '    assigned to the same FFT triplet (n1,n2,n3);      '
        write(stdout,*) '    this suggests that the FFT mesh (nr1,nr2,nr3) is  '
        write(stdout,*) '    too small.                                        '
        stop
      endif

      if (gamma_only) then
        !
        if ( gvec(1, ig)==-gvec(1, ig) .and. &
             gvec(2, ig)==-gvec(2, ig) .and. &
             gvec(3, ig)==-gvec(3, ig) ) assigned(n1, n2, n3) = .false.
        !
        n1 = modulo( -gvec(1, ig), ntm(1) ) + 1
        n2 = modulo( -gvec(2, ig), ntm(2) ) + 1
        n3 = modulo( -gvec(3, ig), ntm(3) ) + 1
        !
        if ( .not. assigned(n1, n2, n3) ) then
          !
          assigned(n1, n2, n3) = .true.
          !
          ! compute the index corresponding to n1,n2,n3 and
          ! assign it to nlm(ig)
          !
          nlm(ig) = n1 + (n2-1)*ntm(1) + (n3-1)*ntm(1)*ntm(2)
          !
        else
          !
          write(stdout,*) ' ERROR in generate_nl. FFT mesh too small?            '
          write(stdout,*) '    More than one G-vector in the gvec array are being'
          write(stdout,*) '    assigned to the same FFT triplet (n1,n2,n3);      '
          write(stdout,*) '    this suggests that the FFT mesh (nr1,nr2,nr3) is  '
          write(stdout,*) '    too small.                                        '
          stop
          !
        endif
        !
      endif

    end do ! ig

    ! compute list_iG and ngk
    allocate(list_iG(nG,nwk), ngk(nwk))
    list_iG = 0
    ngk = 0
    !
    do ik=1,nwk
       klist_cryst(:,ik) = matmul( ainv(bg), wfk(:,ik)/tpiba)
    enddo
    !
    do ik=1,nwk
      !
      list_iG(:,ik) = -1
      ngk(ik) = 0
      do ig=1,nG
        !
        gv = tpiba * matmul(bg, gvec(:,ig) + klist_cryst(:,ik) )
        if ( sum(gv**2) <= cutoff ) then
          ngk(ik) = ngk(ik) + 1
          list_iG(ngk(ik),ik) = ig
        end if
        !
      end do ! ig
      !
    end do ! ik

    ! max. number of G vectors inside the wave function cutoff
    nGk_max = maxval(ngk)

  end subroutine compute_fft_info


  subroutine write_fft_info()
    !
    ! Write all FFT stuff: G vectors, list_iG...
    !

    ! variables
    use writewave, only: nwk
    use siesta2intw_io, only: stdout, intwdir
    ! functions and subroutines
    use siesta2intw_io, only: find_free_unit

    implicit none

    integer :: io_unit
    CHARACTER(LEN=256) :: datafile
    integer :: ig, ik


    write(stdout, *) "- Writing FFT info..."

    ! write G vectors
    io_unit = find_free_unit()
    datafile = trim(intwdir)//"gvectors.dat"
    open(unit=io_unit,file=datafile, status="unknown", action="write", form="unformatted")
    !
    write(unit=io_unit) nG
    do ig=1,nG
      write(unit=io_unit) gvec(:,ig)
    end do
    close(unit=io_unit)

    ! write list_iG
    io_unit = find_free_unit()
    datafile = trim(intwdir)//"iGlist.dat"
    open(unit=io_unit, file=datafile, status="unknown", action="write", form="unformatted")
    !
    write(unit=io_unit) nGk_max
    do ik=1,nwk
       write(unit=io_unit) ngk(ik)
       write(unit=io_unit) list_iG(1:ngk(ik),ik)
    end do
    close(unit=io_unit)

  end subroutine write_fft_info

end module siesta2intw_fft