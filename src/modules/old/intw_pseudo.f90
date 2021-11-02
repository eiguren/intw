module intw_pseudo
  !
  use iotk_module
  !
  use intw_reading
  !
  use intw_useful_constants
  !
  use intw_input_parameters
  !
  use kinds,     only : dp

  USE intw_pseudo_types
  USE intw_upf_module
  USE intw_radial_grids, ONLY : radial_grid_type, nullify_radial_grid, deallocate_radial_grid

  USE intw_uspp_param,       ONLY : upf

  USE intw_atom,             ONLY :  msh, rgrid


  implicit none

!haritz
!  ! variables
 public :: rcut
!  ! subroutines
 public :: allocate_upfeak, deallocate_upfeak, read_all_pseudo, read_upf_tofile, average_pp
!  !
 private
!haritz

  real(DP), parameter :: rcut = 10.d0
  integer :: nt, ir

  save

contains

  subroutine allocate_upfeak ()

    IF( ALLOCATED( rgrid ) ) THEN
     DO nt = 1, SIZE( rgrid )
        CALL deallocate_radial_grid( rgrid( nt ) )
        CALL nullify_radial_grid( rgrid( nt ) )
     END DO
     DEALLOCATE( rgrid )
     DEALLOCATE( msh )
  END IF

  ALLOCATE( rgrid( ntyp ), msh( ntyp ) )

  DO nt = 1, ntyp
     CALL nullify_radial_grid( rgrid( nt ) )
  END DO

  IF( ALLOCATED( upf ) ) THEN
     DO nt = 1, SIZE( upf )
        CALL deallocate_pseudo_upf( upf( nt ) )
        CALL nullify_pseudo_upf( upf( nt ) )
     END DO
     DEALLOCATE( upf )
  END IF
  !
  ALLOCATE ( upf( ntyp ) )
  !

  end subroutine allocate_upfeak

  subroutine deallocate_upfeak ()
    deallocate (upf)
    deallocate (rgrid)
  end subroutine deallocate_upfeak

  !---------------------------------------------------------------------
  subroutine read_all_pseudo ()
    !---------------------------------------------------------------------

    use intw_reading
    use iotk_module

    !---------------------------------------------------------------------
    !
    !  Read pseudopotentials in the Unified Pseudopotential Format (UPF)
    !
    implicit none
    integer :: is, ios, iunps = 455
    character (len=256) :: file_pseudo
        !
    is = 0

    call allocate_upfeak ()

    do is=1,ntyp
       write(*,"(a,i2,3a)") '|     Reading PP file(',is,')= ', trim(atom_pfile(is)),"   |"
       file_pseudo=trim(trim(mesh_dir)//trim(atom_pfile(is)))

       open(unit=iunps,file=file_pseudo,status='old',form='formatted',iostat=ios)
!haritz
       ! if (ios.ne.0) stop
       if (ios.ne.0) stop "Error reading PP file"
!haritz
       call read_upf_tofile (file_pseudo, upf(is), rgrid(is))

       close (unit=iunps)
    enddo

    do nt=1,ntyp
     do ir = 1, rgrid(nt)%mesh
        if (rgrid(nt)%r(ir) > rcut) then
           msh (nt) = ir
           goto 5
        endif
     enddo
     msh (nt) = rgrid(nt)%mesh
     !
     ! force msh to be odd for simpson integration
     !
5    msh (nt) = 2 * ( (msh (nt) + 1) / 2) - 1

    enddo


20 format(A)

  end  subroutine read_all_pseudo

  !---------------------------------------------------------------------
  subroutine read_upf_tofile (file_pseudo, upf_l, grid_l)
    !---------------------------------------------------------------------
    ! Asier 09/10/2014 Adapted from QE. code.
    !
    !   This small program reads the pseudopotential in the Unified
    !   Pseudopotential Format and writes three files
    !   in a format which can be plotted. The files are:
    !
    ! PWSCF modules
    !
    !
    USE constants, ONLY : fpi
    USE intw_pseudo_types
    USE intw_upf_module
    USE intw_radial_grids, ONLY : radial_grid_type, nullify_radial_grid
    !
    IMPLICIT NONE
    !
    INTEGER :: iunps, ierr
    !
    !haritz
    ! CHARACTER(80),intent(in) :: file_pseudo
    CHARACTER(256),intent(in) :: file_pseudo
    !haritz
    TYPE (pseudo_upf), intent(out) :: upf_l
    TYPE (radial_grid_type), intent(out) :: grid_l
    !
    !     Local variables
    !
    INTEGER :: ios, n, j
    !
    !  nullify objects as soon as they are instantiated

    CALL nullify_pseudo_upf( upf_l )
    CALL nullify_radial_grid( grid_l )

    iunps=20
    OPEN(UNIT=iunps,FILE=file_pseudo,STATUS='old',FORM='formatted', &
         ERR=100, IOSTAT=ios)
100 CALL errore('read_upf_tofile','open error on file '//file_pseudo,ios)

    CALL read_upf(upf_l, grid_l, ierr, unit=iunps)
    !
    IF (ierr /= 0) &
         CALL errore('read_upf_tofile','reading pseudo upf', abs(ierr))
    !
    CLOSE(iunps)
    !
  END subroutine read_upf_tofile


!----------------------------------------------------------------------------
SUBROUTINE average_pp ( ntyp )
  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE intw_atom,             ONLY : rgrid
  USE intw_uspp_param,       ONLY : upf
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ntyp
  !
  INTEGER  :: nt, nb, nbe, ind, ind1, l
  REAL(DP) :: vionl
  !
  !
  DO nt = 1, ntyp
     !
     IF ( upf(nt)%has_so ) THEN
        !
        IF ( upf(nt)%tvanp ) &
             CALL errore( 'setup', 'US j-average not yet implemented', 1 )
        !
        nbe = 0
        !
        DO nb = 1, upf(nt)%nbeta
           !
           nbe = nbe + 1
           !
           IF ( upf(nt)%lll(nb) /= 0 .AND. &
                ABS( upf(nt)%jjj(nb) - upf(nt)%lll(nb) - 0.5D0 ) < 1.D-7 ) &
              nbe = nbe - 1
        END DO
        !
        upf(nt)%nbeta = nbe
        !
        nbe = 0
        !
        DO nb = 1, upf(nt)%nbeta
           !
           nbe = nbe + 1
           !
           l = upf(nt)%lll(nbe)
           !
           IF ( l /= 0 ) THEN
              !
              IF (ABS(upf(nt)%jjj(nbe)-upf(nt)%lll(nbe)+0.5d0) < 1.d-7) THEN
                 IF ( ABS( upf(nt)%jjj(nbe+1)-upf(nt)%lll(nbe+1)-0.5d0 ) &
                      > 1.d-7 ) call errore('setup','wrong beta functions',1)
                 ind=nbe+1
                 ind1=nbe
              ELSE
                 IF (ABS(upf(nt)%jjj(nbe+1)-upf(nt)%lll(nbe+1)+0.5d0) > 1.d-7) &
                      call errore('setup','wrong beta functions',2)
                 ind=nbe
                 ind1=nbe+1
              ENDIF
              !
              vionl = ( ( l + 1.D0 ) * upf(nt)%dion(ind,ind) + &
                   l * upf(nt)%dion(ind1,ind1) ) / ( 2.D0 * l + 1.D0 )
              !
              upf(nt)%beta(1:rgrid(nt)%mesh,nb) = 1.D0 / ( 2.D0 * l + 1.D0 ) * &
                   ( ( l + 1.D0 ) * SQRT( upf(nt)%dion(ind,ind) / vionl ) * &
                   upf(nt)%beta(1:rgrid(nt)%mesh,ind) + &
                   l * SQRT( upf(nt)%dion(ind1,ind1) / vionl ) * &
                   upf(nt)%beta(1:rgrid(nt)%mesh,ind1) )
              !
              upf(nt)%dion(nb,nb) = vionl
              !
              nbe = nbe + 1
                 !
           ELSE
              !
              upf(nt)%beta(1:rgrid(nt)%mesh,nb) = &
                  upf(nt)%beta(1:rgrid(nt)%mesh,nbe)
              !
              upf(nt)%dion(nb,nb) = upf(nt)%dion(nbe,nbe)
              !
           END IF
           !
           upf(nt)%lll(nb)=upf(nt)%lll(nbe)
           !
        END DO
        !
        nbe = 0
        !
        DO nb = 1, upf(nt)%nwfc
           !
           nbe = nbe + 1
           !
           IF ( upf(nt)%lchi(nb) /= 0 .AND. &
                ABS(upf(nt)%jchi(nb)-upf(nt)%lchi(nb)-0.5D0 ) < 1.D-7 ) &
              nbe = nbe - 1
           !
        END DO
        !
        upf(nt)%nwfc = nbe
        !
        nbe = 0
        !
        do nb = 1, upf(nt)%nwfc
           !
           nbe = nbe + 1
           !
           l = upf(nt)%lchi(nbe)
           !
           IF ( l /= 0 ) THEN
              !
              IF (ABS(upf(nt)%jchi(nbe)-upf(nt)%lchi(nbe)+0.5d0) < 1.d-7) THEN
                 IF ( ABS(upf(nt)%jchi(nbe+1)-upf(nt)%lchi(nbe+1)-0.5d0) > &
                      1.d-7) call errore('setup','wrong chi functions',3)
                 ind=nbe+1
                 ind1=nbe
              ELSE
                 IF ( ABS(upf(nt)%jchi(nbe+1)-upf(nt)%lchi(nbe+1)+0.5d0) > &
                      1.d-7) call errore('setup','wrong chi functions',4)
                 ind=nbe
                 ind1=nbe+1
              END IF
              !
              upf(nt)%chi(1:rgrid(nt)%mesh,nb) = &
                 ((l+1.D0) * upf(nt)%chi(1:rgrid(nt)%mesh,ind)+ &
                   l * upf(nt)%chi(1:rgrid(nt)%mesh,ind1)) / ( 2.D0 * l + 1.D0 )
              !
              nbe = nbe + 1
              !
           ELSE
              !
              upf(nt)%chi(1:rgrid(nt)%mesh,nb) = upf(nt)%chi(1:rgrid(nt)%mesh,nbe)
              !
           END IF
           !
           upf(nt)%lchi(nb)= upf(nt)%lchi(nbe)
           !
        END DO
        !
     END IF
     !
     upf(nt)%has_so = .FALSE.
     !
  END DO
  !
END SUBROUTINE average_pp


end module intw_pseudo
