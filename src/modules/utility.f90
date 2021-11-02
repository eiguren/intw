!----------------------------------------------------------------------------!
!	intw project.
!
!----------------------------------------------------------------------------!
module intw_utility
!----------------------------------------------------------------------------!
!
!       This module contains useful functions which implement
!       common tasks. It will be VERY useful to call these functions
!       instead of reimplementing them every time they are needed,
!       especially to insure CONSISTENCY.
!
!----------------------------------------------------------------------------!
!haritz
use kinds, only: dp
!haritz
use intw_useful_constants

contains
  subroutine get_timing(time)
  !----------------------------------------------------------------------------!
  !     Timing the code is very important for optimized performance.
  !	However, some experimence indicates that "timing" is not so easy:     
  !	time can overflow the buffer in which it is held, or it can give
  !	crazy numbers for (so far) unknown reasons. Thus, in order to be
  !	able to quickly modify how time is measured, a single subroutine
  !	will be called throuhout the code to measure time. This subroutine
  !	can then easily be modified as the code writer becomes
  !	aware of "better" timing algorithms. 
  !----------------------------------------------------------------------------!
  implicit none

  real(dp),external  :: dsecnd

  real(dp)           :: time

  time = dsecnd()

  end subroutine get_timing
  
 function ainv(a)

    implicit none

    real(dp), dimension(3,3), intent(in)     :: a
    real(dp), dimension(3,3)                 :: ainv

    real(dp), parameter :: eps = 1.0d-13
    real(dp) :: det
    real(dp), dimension(3,3) :: cofactor

    det =   a(1,1)*a(2,2)*a(3,3)  &
         - a(1,1)*a(2,3)*a(3,2)  &
         - a(1,2)*a(2,1)*a(3,3)  &
         + a(1,2)*a(2,3)*a(3,1)  &
         + a(1,3)*a(2,1)*a(3,2)  &
         - a(1,3)*a(2,2)*a(3,1)

    if (abs(det) .le. eps) then
       ainv = 0.0d0
       write(*,*)"ERROR IN CALCULATING MATRIX INVERSE"
       stop
    end if

    cofactor(1,1) = +(a(2,2)*a(3,3)-a(2,3)*a(3,2))
    cofactor(1,2) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
    cofactor(1,3) = +(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    cofactor(2,1) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
    cofactor(2,2) = +(a(1,1)*a(3,3)-a(1,3)*a(3,1))
    cofactor(2,3) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
    cofactor(3,1) = +(a(1,2)*a(2,3)-a(1,3)*a(2,2))
    cofactor(3,2) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
    cofactor(3,3) = +(a(1,1)*a(2,2)-a(1,2)*a(2,1))

    ainv = transpose(cofactor) / det

  end function ainv

 function cmplx_ainv (a)

    implicit none

    complex(dp), dimension(3,3), intent(in)     :: a
    complex(dp), dimension(3,3)                 :: cmplx_ainv

    real(dp), parameter :: eps = 1.0d-13
    complex(dp) :: det
    complex(dp), dimension(3,3) :: cofactor

    det =   a(1,1)*a(2,2)*a(3,3)  &
         - a(1,1)*a(2,3)*a(3,2)  &
         - a(1,2)*a(2,1)*a(3,3)  &
         + a(1,2)*a(2,3)*a(3,1)  &
         + a(1,3)*a(2,1)*a(3,2)  &
         - a(1,3)*a(2,2)*a(3,1)

    if (abs(det) .le. eps) then
       cmplx_ainv = cmplx_0
       write(*,*)"ERROR IN CALCULATING MATRIX INVERSE"
       stop
    end if

    cofactor(1,1) = +(a(2,2)*a(3,3)-a(2,3)*a(3,2))
    cofactor(1,2) = -(a(2,1)*a(3,3)-a(2,3)*a(3,1))
    cofactor(1,3) = +(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    cofactor(2,1) = -(a(1,2)*a(3,3)-a(1,3)*a(3,2))
    cofactor(2,2) = +(a(1,1)*a(3,3)-a(1,3)*a(3,1))
    cofactor(2,3) = -(a(1,1)*a(3,2)-a(1,2)*a(3,1))
    cofactor(3,1) = +(a(1,2)*a(2,3)-a(1,3)*a(2,2))
    cofactor(3,2) = -(a(1,1)*a(2,3)-a(1,3)*a(2,1))
    cofactor(3,3) = +(a(1,1)*a(2,2)-a(1,2)*a(2,1))

    cmplx_ainv = transpose(cofactor) / det

  end function cmplx_ainv

 function cmplx_ainv_2 (a)

    implicit none

    complex(dp), dimension(2,2), intent(in)     :: a
    complex(dp), dimension(2,2)                 :: cmplx_ainv_2

    real(dp), parameter :: eps = 1.0d-13
    complex(dp) :: det
    complex(dp), dimension(2,2) :: cofactor

    det =  a(1,1)*a(2,2)  &
         - a(1,2)*a(2,1)

    if (abs(det) .le. eps) then
       cmplx_ainv_2 = cmplx_0
       write(*,*)"ERROR IN CALCULATING MATRIX INVERSE"
       stop
    end if

    cofactor(1,1) = +a(2,2)
    cofactor(1,2) = -a(2,1)
    cofactor(2,1) = -a(1,2)
    cofactor(2,2) = +a(1,1)

    cmplx_ainv_2 = transpose(cofactor) / det

  end function cmplx_ainv_2

  function conmesurate_and_coarser (nk1, nk2, nk3, nq1, nq2, nq3)
   !Asier&&Idoia 23 06 2014
   logical :: conmesurate_and_coarser
   integer, intent (in) :: nk1, nk2, nk3, nq1, nq2, nq3

   conmesurate_and_coarser=.true.
 
   !check if electron k grid does not contain the phonon q grid.
   if  ((nk1<nq1).or.(nk2<nq2).or.(nk3<nq3))  conmesurate_and_coarser=.false.

   !check if q is contained in k

   if (.not.multiple (nk1, nq1)) conmesurate_and_coarser=.false.
   if (.not.multiple (nk2, nq2)) conmesurate_and_coarser=.false.
   if (.not.multiple (nk3, nq3)) conmesurate_and_coarser=.false.

  end function conmesurate_and_coarser

   function multiple (i1, i2)
   !Asier&&Idoia 23 06 2014
   !Check if two integers are multiple of each other.
   logical :: multiple
   integer, intent (in) :: i1, i2 
 
   multiple=.false.

   if ((i1/i2*i2==i1).or.(i2/i1*i1==i2)) then
      multiple=.true.
   endif 

   end function  multiple 
  subroutine switch_indices(nk_1,nk_2,nk_3,ikpt,i,j,k,switch)
  !----------------------------------------------------------------------------!
  !     This subroutine efficiently computes the scalar kpoint index from
  !     the triplet index, or vice versa. The mesh is described by nk_1,
  !     nk_2, nk_3, which can be the coarse or the smooth mesh.
  !
  !     ikpt is the scalar index, namely a single integer from 1 to nkpoints
  !     which uniquely labels all the kpoints in the full 1BZ.
  !
  !     i,j,k are the triplet index, which labels the kpoints according to
  !     their coordinates in the MP mesh. 
  !
  !     Following Asier's convention, the 3rd triplet index will loop 
  !     fastest. The convention is:
  !
  !             ikpt = (k-1) + (j-1)*nk_3 + (i-1)*nk_2*nk_3 + 1
  !
  !     in general:
  ! 
  !     switch = +1   :  ikpt    =       1+ (k-1) + (j-1)*nk_3 + (i-1)*nk_2*nk_3
  !
  !     switch = -1   :   i      =       1 + (ikpt-1) / (nk_2*nk_3)
  !                       j      =       1 + (ikpt-1) / nk_3-(i-1)*nk_2
  !                       k      =      ikpt -(j-1)*nk_3-(i-1)*nk_2*nk_3
  !----------------------------------------------------------------------------!
  implicit none

  integer        :: i,       j,       k     
  integer        :: i_m_1,   j_m_1,   k_m_1     

  integer        :: nk_1,  nk_2,  nk_3

  integer        :: switch 

  integer        :: ikpt,  ikpt_m_1 
  
  integer        :: n2n3

  if (switch .eq. 1) then

        !       efficient implementation; 3 +/-,  2 */.
        ikpt    =       k + nk_3*(j-1 + (i-1)*nk_2)
  
  else if (switch .eq. -1) then

        ! define useful variables to minimize the number of 
        ! operations. 
        
        !       The scheme below involves 6 +- and 6 */


        n2n3     = nk_2*nk_3
        ikpt_m_1 = ikpt-1


        i_m_1 = ikpt_m_1/n2n3
        i     = i_m_1 + 1

        j_m_1 = ikpt_m_1/nk_3-i_m_1*nk_2
        j     = j_m_1+1

        k     = ikpt - j_m_1*nk_3-i_m_1*n2n3

  end if

  end subroutine switch_indices

  subroutine switch_indices_zyx(nk_1,nk_2,nk_3,ikpt,i,j,k,switch)
  !----------------------------------------------------------------------------!
  !     This subroutine efficiently computes the scalar kpoint index from
  !     the triplet index, or vice versa. The mesh is described by nk_1,
  !     nk_2, nk_3, which can be the coarse or the smooth mesh.
  !
  !     ikpt is the scalar index, namely a single integer from 1 to nkpoints
  !     which uniquely labels all the kpoints in the full 1BZ.
  !
  !     i,j,k are the triplet index, which labels the kpoints according to
  !     their coordinates in the MP mesh. 
  !
  !     Following Asier's convention, the 3rd triplet index will loop 
  !     fastest. The convention is:
  !  
  !     New convention (Asier&Idoia 20 06 2014) for a coherence with our cfftd code:
  !
  !         Old 
  !             ikpt = (k-1) + (j-1)*nk_3 + (i-1)*nk_2*nk_3 + 1
  !         New 20 06 2014 
  !             ikpt = (i-1) + (j-1)*nk_1 + (k-1)*nk_1*nk2 + 1
  !
  !     in general:
  ! 
  !New (20 06 2014):
  !     switch = +1   :  ikpt    =       1+ (i-1) + (j-1)*nk_1 + (k-1)*nk_1*nk_2
  !
  !     switch = -1   :   i      =       1 + (ikpt-1) / (nk_1*nk_2)
  !                       j      =       1 + (ikpt-1) / nk_1-(k-1)*nk_2
  !                       k      =      ikpt -(j-1)*nk_1-(k-1)*nk_1*nk_2

  !Old:
  !     switch = +1   :  ikpt    =       1+ (k-1) + (j-1)*nk_3 + (i-1)*nk_2*nk_3
  !
  !     switch = -1   :   i      =       1 + (ikpt-1) / (nk_2*nk_3)
  !                       j      =       1 + (ikpt-1) / nk_3-(i-1)*nk_2
  !                       k      =      ikpt -(j-1)*nk_3-(i-1)*nk_2*nk_3

  !----------------------------------------------------------------------------!
  implicit none

  integer        :: i,       j,       k     
  integer        :: i_m_1,   j_m_1,   k_m_1     

  integer        :: nk_1,  nk_2,  nk_3

  integer        :: switch 

  integer        :: ikpt,  ikpt_m_1 
  
  integer        :: n1n2

  if (switch .eq. 1) then

         ikpt    =       i + (j-1)*nk_1 + (k-1)*nk_1*nk_2  
  
  else if (switch .eq. -1) then

        ! define useful variables to minimize the number of 
        ! operations. 
        
        !       The scheme below involves 6 +- and 6 */


        n1n2     = nk_1*nk_2
        ikpt_m_1 = ikpt-1


        k_m_1 = ikpt_m_1/n1n2
        k     = k_m_1 + 1

        j_m_1 = ikpt_m_1/nk_1-k_m_1*nk_2
        j     = j_m_1+1

        i     = ikpt - j_m_1*nk_1-k_m_1*n1n2

  end if

  end subroutine switch_indices_zyx

  subroutine generate_kmesh(kmesh,nk_1,nk_2,nk_3)
  !----------------------------------------------------------------------------!
  !     This subroutine builds the array of k vectors corresponding
  !     to MP indices nk_1,nk_2,nk_3, ordered according to their singlet
  !     index.  
  !----------------------------------------------------------------------------!
  implicit none

  real(dp)       :: kmesh(3,nk_1*nk_2*nk_3)
  integer        :: nk_1,  nk_2,  nk_3
  integer        :: i,     j,     k 
  
  integer        :: ikpt , nkmesh
  integer        :: switch 

  ! Build kmesh
  nkmesh = nk_1*nk_2*nk_3

  switch = -1 ! singlet to triplet

  do ikpt = 1, nkmesh
        call switch_indices(nk_1,nk_2,nk_3,ikpt,i,j,k,switch)

        kmesh(1, ikpt) = dble(i-1)/nk_1
        kmesh(2, ikpt) = dble(j-1)/nk_2
        kmesh(3, ikpt) = dble(k-1)/nk_3
  end do


  end subroutine generate_kmesh

  subroutine find_neighbor(kpoint,nk_1,nk_2,nk_3,i_k,j_k,k_k)
  !----------------------------------------------------------------------------!
  !     Given a kpoint = (kx, ky, kz) in crystal coordinates, which lies
  !     inside the 1BZ, namely 0 <= kx,ky,kz < 1, this subroutine returns
  !     the coordinates of the neighbor which lies at the origin of the tricubic
  !     4x4x4 interpolation grid, in the format (i_k,j_k,k_k). 
  !     The crystal coordinates of the neigbhor are given by 
  !                     ((i_k-1)/nk_1, (j_k-1)/nk_2, (k_k-1)/nk_3). 
  !
  !----------------------------------------------------------------------------!
  implicit none

  real(dp)       :: kpoint(3)
  integer        :: nk_1,  nk_2,  nk_3
  
  integer        :: i_k,    j_k,    k_k

  i_k    =       1+nint(kpoint(1)*nk_1)
  j_k    =       1+nint(kpoint(2)*nk_2)
  k_k    =       1+nint(kpoint(3)*nk_3)

  end subroutine find_neighbor

  subroutine find_maximum_index_int(array,sze,i_max)
  !----------------------------------------------------------------------------!
  !     Given an array of integers array(sze), this subroutine returns 
  !     the largest value in the array.                
  !----------------------------------------------------------------------------!
  implicit none

  integer        :: sze,    maximum,   array(sze) 
  
  integer        :: i,  i_max 

  i_max   = 1
  maximum = array(1)

  do i=2,sze
        if ( array(i) > maximum) then
                maximum = array(i)
                i_max   = i 
        end if
  end do

  end subroutine find_maximum_index_int

  subroutine test_qpt_on_fine_mesh(qpt,nk1s,nk2s,nk3s,test_qpt,i_qpt1,i_qpt2,i_qpt3)
  !----------------------------------------------------------------------------!
  !     
  !    This subroutine tests whether the q-point, qpt(3), is of the form 
  !            
  !     qpt(:) = [ i_qpt1-1   i_qpt2-1   i_qpt3-1 ] + (I,J,K)
  !              [ -------- , -------- , -------- ]
  !              [   nk1s       nk2s       nk3s   ]
  !            
  !
  !	where I,J,K are integers, and i_qpt(1,2,3) are integers between 1 and nk(123)s.
  !
  !	It is important to perform this task in a defensive way. 
  !     The fortran internal functions nint, floor, and modulo are DANGEROUS. 
  !     IT IS CRUCIAL TO DO THIS RIGHT ONCE AND FOR ALL. 
  !
  !----------------------------------------------------------------------------!
  implicit none
  
  ! input
  integer        :: nk1s, nk2s, nk3s       ! the fine mesh parameters
  real(dp)       :: qpt(3)                 ! the q-point 
  
  ! output
  logical        :: test_qpt		   ! is the q-point on the fine mesh?
  integer        :: i_qpt1,i_qpt2,i_qpt3    

  real(dp)       :: nqpt_1,nqpt_2,nqpt_3


  ! internal variables
  real(dp)       :: err_1, err_2, err_3   ! the q-point 

  ! assume qpt = (i-1)/n + I + err,    0<= err < 1

  ! find err
  nqpt_1 = dble(nk1s)*qpt(1)
  nqpt_2 = dble(nk2s)*qpt(2)
  nqpt_3 = dble(nk3s)*qpt(3)

  err_1 = (nqpt_1-nint(nqpt_1))/dble(nk1s)
  err_2 = (nqpt_2-nint(nqpt_2))/dble(nk2s)
  err_3 = (nqpt_3-nint(nqpt_3))/dble(nk3s)


  test_qpt =  err_1 <  eps_8 .and.  &
              err_2 <  eps_8 .and.  &
              err_3 <  eps_8

  if (test_qpt) then  
     ! q is on mesh; find its coordinates.

     i_qpt1 = modulo(nint(nqpt_1),nk1s)+1
     i_qpt2 = modulo(nint(nqpt_2),nk2s)+1
     i_qpt3 = modulo(nint(nqpt_3),nk3s)+1

  else
     ! q is not on mesh! throw an error outside the subroutine
     i_qpt1 = 0
     i_qpt2 = 0
     i_qpt3 = 0
  end if

  end subroutine test_qpt_on_fine_mesh

  subroutine find_k_1BZ_and_G(kpoint,nk1,nk2,nk3,i,j,k,kpt_in_1BZ,G)
  !----------------------------------------------------------------------------!
  !     Given a kpoint(3) in crystal coordinates, this subroutine 
  !     generates:
  !            - k_1BZ(3),   with 0 <= k_1BZ(i) < 1 
  !            - G(3),       G(i) an integer
  !            - i,j,k       the triplet coordinates of the point in the mesh.
  !
  !	It it EXTREMELY important to have a subroutine which performs this
  !	task in a defensive way. The fortran internal functions nint, floor, 
  !	and modulo are DANGEROUS. IT IS CRUCIAL TO DO THIS RIGHT ONCE AND FOR
  !     ALL. 
  !
  !	The basic relationship is 
  !
  !		kpoint = k_1BZ + G
  !
  !     Assume 
  !		kpoint(l) =  G(l) + (i_l-1)/nk_l+ epsilon ,   
  !                  with 
  !                            l    = 1, 2, 3
  !                            i_l  = i, j ,k 
  !                            nk_l = nk1, nk2, nk3
  !                            G(l) an integer         
  !                         epsilon numerical noise
  !----------------------------------------------------------------------------!
  implicit none
  
  ! input
  integer        :: nk1, nk2, nk3             ! the mesh parameters
  real(dp)       :: kpoint(3)                 ! the kpoint of interest 
  
  ! output
  integer        :: i, j, k                   ! the triplet coordinates of k_1BZ
  real(dp)       :: kpt_in_1BZ(3)             ! the k point in the 1BZ
  integer        :: G(3)                      ! the translation vector

  ! internal variables
  integer        :: nG_im1, nG_jm1, nG_km1 
  integer        :: im1,    jm1,    km1 


  ! this step kills epsilon 
   nG_im1 = nint(kpoint(1)*dble(nk1))
   nG_jm1 = nint(kpoint(2)*dble(nk2))
   nG_km1 = nint(kpoint(3)*dble(nk3))

  ! this step gets rid of G
  im1     = modulo(nG_im1,nk1)
  jm1     = modulo(nG_jm1,nk2)
  km1     = modulo(nG_km1,nk3)

  ! G can be extracted. This division must be exact.
  G(1)    = (nG_im1 - im1)/nk1
  G(2)    = (nG_jm1 - jm1)/nk2
  G(3)    = (nG_km1 - km1)/nk3

  ! finally we have the triplet coordinates

  i       = im1 + 1
  j       = jm1 + 1
  k       = km1 + 1

  ! compute the k point in the 1BZ
  kpt_in_1BZ(1) = dble(i-1)/dble(nk1)
  kpt_in_1BZ(2) = dble(j-1)/dble(nk2)
  kpt_in_1BZ(3) = dble(k-1)/dble(nk3)

  end subroutine find_k_1BZ_and_G


SUBROUTINE HPSORT(N,RA,P)
!------------------------------------------------------------
! subroutine which performs heap sort on a list of integers
! and also returns an array identifying the permutation
! which sorted the array.
!
! The subroutine was copied from the internet, and slightly
! modified to handle integers and return the permutation
! array. 
!
! Part of the original header follows:
!------------------------------------------------------------

!*****************************************************
!*  Sorts an array RA of length N in ascending order *
!*                by the Heapsort method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table RA                   *
!*          RA	  table to be sorted                 *
!* OUTPUT:                                           *
!*	    RA    table sorted in ascending order    *
!*	    P     table of indices showing transform *
!*                                                   *
!* NOTE: The Heapsort method is a N Log N routine,   *
!*       and can be used for very large arrays.      *
!* ------------------------------------------------- *
!* REFERENCE:                                        *
!*  "NUMERICAL RECIPES by W.H. Press, B.P. Flannery, *
!*   S.A. Teukolsky and W.T. Vetterling, Cambridge   *
!*   University Press, 1986".                        *
!*****************************************************         
  integer N
  integer RA(N)
  integer P(N)

  integer i 
 
  do i=1,N
    P(i) = i
  end do

  L=N/2+1
  IR=N
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
    L=L-1

    RRA=RA(L)
    PP =P(L)

  else
    RRA=RA(IR)
    PP =P (IR)

    RA(IR)=RA(1)
    P(IR) =P(1)

    IR=IR-1
    if(IR.eq.1)then

      RA(1)=RRA
      P (1)=PP 

      return
    end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
  if(J < IR)then
    if(RA(J) < RA(J+1))  J=J+1
  end if
  if(RRA < RA(J))then
    RA(I)=RA(J)
    P (I)=P (J)

    I=J; J=J+J
  else
    J=IR+1
  end if
  goto 20
  end if
  RA(I)=RRA
  P (I)=PP 
  goto 10
END SUBROUTINE HPSORT


    function find_free_unit()
    !--------------------------------------------------------------------------
    ! This function finds a free input/output unit. It is copied from
    ! the QE distribution and slightly modified to suit our needs.
    !--------------------------------------------------------------------------
    !
    implicit none
    !
    integer :: find_free_unit
    integer :: io_unit
    logical :: opnd

    do io_unit = 99, 1, -1

       inquire( unit = io_unit, opened = opnd )
       if ( .not. opnd ) then
          find_free_unit = io_unit
          return
       end if
       !
    end do 
    !
    write(*,*) 'NO free units available!?!' 
    !
    return
    end function find_free_unit


  subroutine create_or_append_to_file(io_unit,filename)
  !-----------------------------------------------------------------------
  !
  !     This routine tests whether the file with name filename
  !     exists, and opens it for appending or creates it correspondingly
  !     on unit io_unit.
  !-----------------------------------------------------------------------

  implicit none
  integer      ::    io_unit
  logical      ::    file_exists
  character(*) ::    filename

  inquire(file=filename, exist=file_exists) 

  if (file_exists) then
     open(unit=io_unit,file=filename,access='append',status='old')
  else
     open(unit=io_unit,file=filename,status='new')
  end if


  end subroutine create_or_append_to_file

  subroutine cryst_to_cart (nvec, vec, trmat, iflag)
  !-----------------------------------------------------------------------
  !
  !     This routine transforms the atomic positions or the k-point
  !     components from crystallographic to cartesian coordinates 
  !     ( iflag=1 ) and viceversa ( iflag=-1 ).
  !     Output cartesian coordinates are stored in the input ('vec') array
  !
  !-----------------------------------------------------------------------
  !
  implicit none
  !
  integer, intent(in) :: nvec, iflag
  ! nvec:  number of vectors (atomic positions or k-points)
  !        to be transformed from crystal to cartesian and vice versa
  ! iflag: gives the direction of the transformation
  real(8), intent(in) :: trmat (3, 3)
  ! trmat: transformation matrix
  ! if iflag=1:
  !    trmat = at ,  basis of the real-space lattice,       for atoms   or
  !          = bg ,  basis of the reciprocal-space lattice, for k-points
  ! if iflag=-1: the opposite
  real(8), intent(inout) :: vec (3, nvec)
  ! coordinates of the vector (atomic positions or k-points) to be
  ! transformed - overwritten on output
  !
  !    local variables
  !
  integer :: nv, kpol
  ! counter on vectors
  ! counter on polarizations
  real(8) :: vau (3)
  ! workspace
  !
  !     Compute the cartesian coordinates of each vectors
  !     (atomic positions or k-points components)
  !
  do nv = 1, nvec
     if (iflag.eq.1) then
        do kpol = 1, 3
           vau (kpol) = trmat (kpol, 1) * vec (1, nv) + trmat (kpol, 2) &
                * vec (2, nv) + trmat (kpol, 3) * vec (3, nv)
        enddo
     else
        do kpol = 1, 3
           vau (kpol) = trmat (1, kpol) * vec (1, nv) + trmat (2, kpol) &
                * vec (2, nv) + trmat (3, kpol) * vec (3, nv)
        enddo
     endif
     do kpol = 1, 3
        vec (kpol, nv) = vau (kpol)
     enddo
  enddo
  !
  return
  end subroutine cryst_to_cart


SUBROUTINE HPSORT_real(N,RA,P)
!------------------------------------------------------------
! subroutine which performs heap sort on a list of real numbers
! and also returns an array identifying the permutation
! which sorted the array.
!
! The subroutine was copied from the internet, and slightly
! modified to handle integers and return the permutation
! array. 
!
! Part of the original header follows:
!------------------------------------------------------------

!*****************************************************
!*  Sorts an array RA of length N in ascending order *
!*                by the Heapsort method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table RA                   *
!*          RA	  table to be sorted                 *
!* OUTPUT:                                           *
!*	    RA    table sorted in ascending order    *
!*	    P     table of indices showing transform *
!*                                                   *
!* NOTE: The Heapsort method is a N Log N routine,   *
!*       and can be used for very large arrays.      *
!* ------------------------------------------------- *
!* REFERENCE:                                        *
!*  "NUMERICAL RECIPES by W.H. Press, B.P. Flannery, *
!*   S.A. Teukolsky and W.T. Vetterling, Cambridge   *
!*   University Press, 1986".                        *
!*****************************************************         
  integer  N
  real(dp) RA(N)
  integer  P(N)

  integer i 
 
  do i=1,N
    P(i) = i
  end do

  L=N/2+1
  IR=N
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
    L=L-1

    RRA=RA(L)
    PP =P(L)

  else
    RRA=RA(IR)
    PP =P (IR)

    RA(IR)=RA(1)
    P(IR) =P(1)

    IR=IR-1
    if(IR.eq.1)then

      RA(1)=RRA
      P (1)=PP 

      return
    end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
  if(J < IR)then
    if(RA(J) < RA(J+1))  J=J+1
  end if
  if(RRA < RA(J))then
    RA(I)=RA(J)
    P (I)=P (J)

    I=J; J=J+J
  else
    J=IR+1
  end if
  goto 20
  end if
  RA(I)=RRA
  P (I)=PP 
  goto 10
END SUBROUTINE HPSORT_real

!haritz
!subroutine check_netcdf(status)
!    use netcdf
!    integer, intent ( in) :: status
!
!    if(status /= nf90_noerr) then
!      print *, trim(nf90_strerror(status))
!      stop "Stopped"
!    end if
!end subroutine check_netcdf
!haritz

  subroutine find_r_in_WS_cell(at,rvec_cryst,nr1,nr2,nr3,rvec_WS_cryst)
  !----------------------------------------------------------------------------!
  !     Given a real space vector in crystal coordinates rvec_cryst(3),
  !	which has coordinates on a nr1 x nr2 x nr3 grid,
  !	this subroutine finds the coordinates of the corresponding vector
  !     in the Wigner Seitz cell, namely the corresponding vector which is
  !	closest to the origin.
  !----------------------------------------------------------------------------!
  implicit none
  
  ! input
  real(dp)       :: at(3,3)                   ! the real space basis vectors
  real(dp)       :: rvec_cryst(3)             ! the real space vector
  integer        :: nr1, nr2, nr3             ! the real mesh parameters
  
  ! output
  real(dp)       :: rvec_WS_cryst(3)   	      ! vector in the WS cell

  ! internal variables
  integer        :: i, j, k                   ! the triplet coordinates of k_1BZ
  integer        :: Rlat(3)                   ! the translation vector


  real(dp)       :: r_UC(3)   	              ! temporary UC vector
  real(dp)       :: r_WS(3)   	              ! temporary WS vector
  real(dp)       :: list_T(3)                 ! translation coefficients 
  real(dp)       :: T1, T2, T3
  integer        :: it1, it2, it3

  real(dp)       :: ai_dot_aj(3,3)            ! inner product of basis vectors
  real(dp)       :: square_norm
  real(dp)       :: minimum_square_norm

  ! First, find the coordinates of the vector in the unit cell coordinates 
  ! Note that the subroutine below was initially designed for reciprocal
  ! space, but the relevant algorithm is the same.
  call find_k_1BZ_and_G(rvec_cryst,nr1,nr2,nr3,i,j,k,r_UC,Rlat)



  ! a_i^alpha = at(i,alpha), i: basis index, alpha: space index
  ai_dot_aj = matmul(at(:,:),transpose(at(:,:)))

  list_T(:) = (/ -one, zero, one /)


  ! initialize
  minimum_square_norm = dot_product(r_UC,matmul(ai_dot_aj(:,:),r_UC))
  rvec_WS_cryst(:)    = r_UC(:)

  do it1 = 1, 3
     T1  = list_T(it1)
     do it2 = 1, 3
        T2  = list_T(it2)
        do it3 = 1, 3
           T3  = list_T(it3)

           ! Define a possible vector in cartesian coordinates as
           ! r_WS = sum_{i=1}^3 rvec_US(i) a_i
           !
           ! its norm^2 is given by r_WS*r_WS

	   r_WS(:)  =  r_UC(:) + (/T1,T2,T3/)

	   square_norm = dot_product(r_WS,matmul(ai_dot_aj(:,:),r_WS))

	   if ( square_norm < minimum_square_norm ) then
		minimum_square_norm  = square_norm 
  		rvec_WS_cryst(:)     = r_WS(:)

           end if


        end do !it3 
     end do !it2    
  end do !it1    

  end subroutine find_r_in_WS_cell

function cmplx_trace (mat)

complex(kind=dp) :: mat(:,:), cmplx_trace
integer :: d1,d2, i

d1 = size (mat(:,1))
d2 = size (mat(1,:))

if (d1.ne.d2) then 
 write(*,*)"Errorea cmplx_trace"
 stop
end if

cmplx_trace = cmplx_0
do i=1,d1
 cmplx_trace = cmplx_trace + mat(i,i)
enddo

end function  cmplx_trace
!---------------------------------------
!***************************************
!---------------------------------------
function weight_ph(x)

  implicit none

  real(dp) :: x,weight_ph
  real(dp),parameter :: x0=0.002d0 !x0=0.0001d0 !x0=0.0004d0 !x0=0.0006d0 !x0=0.003
  real(dp),parameter :: sigma=0.00001d0 !sigma=0.000001d0 !sigma=0.0001d0 !sigma=0.00002d0 !sigma=0.0003

  weight_ph=1/(1+exp((x0-x)/sigma))

end function weight_ph
!---------------------------------------
!***************************************
!---------------------------------------
!
subroutine errore (a,b,i)
character(len=*), intent(in) :: a,b
integer :: i
write(*,*)"ERROR:", a,b,i
stop
end subroutine errore 

!-----------------------------------------------------------------------
subroutine simpson (mesh, func, rab, asum)
  !-----------------------------------------------------------------------
  !
  !     simpson's rule integration. On input:
  !       mesh = mhe number of grid points (should be odd)
  !       func(i)= function to be integrated
  !       rab(i) = r(i) * dr(i)/di * di
  !     For the logarithmic grid not including r=0 :
  !       r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
  !     For the logarithmic grid including r=0 :
  !       r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
  !     Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr 
  !     where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
  !
  use kinds, ONLY: DP
  implicit none
  integer, intent(in) :: mesh
  real(DP), intent(in) :: rab (mesh), func (mesh)
  real(DP), intent(out):: asum
  !
  real(DP) :: f1, f2, f3, r12
  integer :: i
  !
  !     routine assumes that mesh is an odd number so run check
  !     if ( mesh+1 - ( (mesh+1) / 2 ) * 2 .ne. 1 ) then
  !       write(*,*) '***error in subroutine radlg'
  !       write(*,*) 'routine assumes mesh is odd but mesh =',mesh+1
  !       stop
  !     endif
  asum = 0.0d0
  r12 = 1.0d0 / 12.0d0
  f3 = func (1) * rab (1) * r12

  do i = 2, mesh - 1, 2
     f1 = f3
     f2 = func (i) * rab (i) * r12
     f3 = func (i + 1) * rab (i + 1) * r12
     asum = asum + 4.0d0 * f1 + 16.0d0 * f2 + 4.0d0 * f3
  enddo

  return
end subroutine simpson

!---------------------------------------------------------------------
function qe_erf (x)  
  !---------------------------------------------------------------------
  !
  !     Error function - computed from the rational approximations of
  !     W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
  !
  !     for abs(x) le 0.47 erf is calculated directly
  !     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
  !
  use kinds, only : DP
  implicit none  
  real(DP), intent(in) :: x
  real(DP) :: x2, p1 (4), q1 (4)
  !real(DP), external :: qe_erfc  
  real(DP) :: qe_erf
  data p1 / 2.426679552305318E2_DP, 2.197926161829415E1_DP, &
            6.996383488619136_DP,  -3.560984370181538E-2_DP /
  data q1 / 2.150588758698612E2_DP, 9.116490540451490E1_DP, &
            1.508279763040779E1_DP, 1.000000000000000_DP /
  !
  if (abs (x) > 6.0_DP) then  
     !
     !  erf(6)=1-10^(-17) cannot be distinguished from 1
     !
     qe_erf = sign (1.0_DP, x)  
  else  
     if (abs (x)  <= 0.47_DP) then  
        x2 = x**2  
        qe_erf=x *(p1 (1) + x2 * (p1 (2) + x2 * (p1 (3) + x2 * p1 (4) ) ) ) &
                / (q1 (1) + x2 * (q1 (2) + x2 * (q1 (3) + x2 * q1 (4) ) ) )
     else  
        qe_erf = 1.0_DP - qe_erfc (x)  
     endif
  endif
  !
  return  
end function qe_erf

!---------------------------------------------------------------------
function qe_erfc (x)  
  !---------------------------------------------------------------------
  !
  !     erfc(x) = 1-erf(x)  - See comments in erf
  !
  use kinds, only : DP
  implicit none  
  real(DP),intent(in) :: x
  real(DP)            :: qe_erfc
  real(DP) :: ax, x2, xm2, p2 (8), q2 (8), p3 (5), q3 (5), pim1
  !real(DP), external :: qe_erf  
  data p2 / 3.004592610201616E2_DP,  4.519189537118719E2_DP, &
            3.393208167343437E2_DP,  1.529892850469404E2_DP, &
            4.316222722205674E1_DP,  7.211758250883094_DP,   &
            5.641955174789740E-1_DP,-1.368648573827167E-7_DP /
  data q2 / 3.004592609569833E2_DP,  7.909509253278980E2_DP, &
            9.313540948506096E2_DP,  6.389802644656312E2_DP, &
            2.775854447439876E2_DP,  7.700015293522947E1_DP, &
            1.278272731962942E1_DP,  1.000000000000000_DP /
  data p3 /-2.996107077035422E-3_DP,-4.947309106232507E-2_DP, &
           -2.269565935396869E-1_DP,-2.786613086096478E-1_DP, &
           -2.231924597341847E-2_DP /
  data q3 / 1.062092305284679E-2_DP, 1.913089261078298E-1_DP, &
            1.051675107067932_DP,    1.987332018171353_DP,    &
            1.000000000000000_DP /

  data pim1 / 0.56418958354775629_DP /  
  !        ( pim1= sqrt(1/pi) )
  ax = abs (x)  
  if (ax > 26.0_DP) then  
     !
     !  erfc(26.0)=10^(-296); erfc( 9.0)=10^(-37);
     !
     qe_erfc = 0.0_DP  
  elseif (ax > 4.0_DP) then  
     x2 = x**2  
     xm2 = (1.0_DP / ax) **2  
     qe_erfc = (1.0_DP / ax) * exp ( - x2) * (pim1 + xm2 * (p3 (1) &
          + xm2 * (p3 (2) + xm2 * (p3 (3) + xm2 * (p3 (4) + xm2 * p3 (5) &
          ) ) ) ) / (q3 (1) + xm2 * (q3 (2) + xm2 * (q3 (3) + xm2 * &
          (q3 (4) + xm2 * q3 (5) ) ) ) ) )
  elseif (ax > 0.47_DP) then  
     x2 = x**2  
     qe_erfc = exp ( - x2) * (p2 (1) + ax * (p2 (2) + ax * (p2 (3) &
          + ax * (p2 (4) + ax * (p2 (5) + ax * (p2 (6) + ax * (p2 (7) &
          + ax * p2 (8) ) ) ) ) ) ) ) / (q2 (1) + ax * (q2 (2) + ax * &
          (q2 (3) + ax * (q2 (4) + ax * (q2 (5) + ax * (q2 (6) + ax * &
          (q2 (7) + ax * q2 (8) ) ) ) ) ) ) )
  else  
     qe_erfc = 1.0_DP - qe_erf (ax)  
  endif
  !
  ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
  !
  if (x < 0.0_DP) qe_erfc = 2.0_DP - qe_erfc  
  !
  return  
end function qe_erfc

!
end module intw_utility
!
!
!----------------------------------------------------------------------------!

