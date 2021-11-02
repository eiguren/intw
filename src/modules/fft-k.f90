
  subroutine nl_sort(ngm, nr1, nr2, nr3, gvec, nl)

    implicit none
  
    integer, intent(in) :: ngm, nr1, nr2, nr3
    integer, intent(in) :: gvec(1:3, ngm)
    integer, intent(out) :: nl(ngm)

    !local

    integer  :: n1, n2, n3
   
    integer  :: ng

    do ng = 1, ngm

       n1       = modulo(gvec(1,ng), nr1)+1 !Kuboan sartu
       n2       = modulo(gvec(2,ng), nr2)+1
       n3       = modulo(gvec(3,ng), nr3)+1
                
       nl(ng) = n1 + (n2-1)*nr1 + (n3-1)*nr1*nr2 

  end subroutine nl_sortu

  subroutine wfc_from_g_to_r (list_iG,wfc_g, wfc_r)

  implicit none

  complex(dp), intent(in)   :: wfc_g(nG_max)
  complex(dp), intent(out)   :: wfc_r(nr1*nr2*nr3)

  integer       :: i,  iG 
  integer       :: list_iG(nG_max)

  complex(dp), intent(in)   :: wfc_g(nG_max)
  complex(dp)   :: wfc_r(nr1*nr2*nr3)

  ! initialize work array
  wfc_r(:)  =  cmplx_0 

  ! put wfc_g in wfc_r 
  do i=1,nG_max
      ! identify the G vector by its index, as stored in list_iG
      iG = list_iG(i)
      if (iG == 0) exit
      wfc_r(nl(iG)) = wfc_g(i)
  enddo

   call cfftnd(3,(/nr1,nr2,nr3/),1,wfc_r) ! 

  end subroutine wfc_from_g_to_r 

  subroutine wfc_from_r_to_g (list_iG,wfc_r, wfc_g)
  implicit none

  integer       :: i,  iG 
  integer       :: list_iG(nG_max)

  complex(dp),intent(out)   :: wfc_g(nG_max)
  complex(dp),intent(in )   :: wfc_r(nr1*nr2*nr3)

  ! initialize work array
  wfc_g(:)  =  cmplx_0 

  call cfftnd(3,(/nr1,nr2,nr3/),-1,wfc_r) ! 

  do ig=1,ng_max
       wfc_g(ig)=wfc_r(nl(ig))
  enddo 

  end subroutine wfc_from_r_to_g
 
