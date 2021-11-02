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
       if (ios.ne.0) stop
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
    CHARACTER(80),intent(in) :: file_pseudo
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



end module intw_pseudo
