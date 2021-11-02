program irakurri_geak

integer, parameter :: dp=selected_real_kind(14)
  character(256) :: mesh_dir
  character(256) :: prefix

mesh_dir="./"
prefix="mgb2"

call getngm (mesh_dir, prefix, ngm)


contains

  subroutine getngm (mesh_dir, prefix, ngm)
  use iotk_module

  implicit none
  integer :: io_unit

  character(256) :: mesh_dir
  character(256) :: prefix


  !output: zenbat g maximo
  !local :
  character(256) :: gvec_file
  character(256) :: attr  
  
  integer, intent(out) :: ngm

   gvec_file=trim(trim(mesh_dir)//trim(prefix)//".save/"//"gvectors.dat")

   io_unit  = 31415 
   call iotk_open_read (io_unit,trim(gvec_file))

      call iotk_scan_empty(io_unit,"INFO",attr)
      call iotk_scan_attr(attr,"gvect_number", ngm)
   
   call iotk_close_read (io_unit)


  end subroutine getngm

  !--------------------------------------------------
  !Read global gvec point data size= ngm 
  !---------------------------------------------------
  subroutine getgvec ( mesh_dir, prefix, ngm, gvec )
  use iotk_module


  implicit none
  integer :: io_unit

  character(256) :: mesh_dir
  character(256) :: prefix

  integer, intent(in) :: ngm
  !input 
  ! output 
  !local
  character(256) :: gvec_file
  character(256) :: attr  
 
  !out g bektoreak kristal koordenatuetan
  real(kind=dp), intent(out) :: gvec(1:3, ngm) 

   gvec_file=trim(trim(mesh_dir)//trim(prefix)//".save/"//"gvectors.dat")

   gvec=0
   io_unit  = 3141592 
   call iotk_open_read (io_unit,trim(gvec_file))

      call iotk_scan_empty(io_unit,"INFO",attr)
      CALL iotk_scan_dat  ( io_unit, "g", gvec(1:3,1:ngm))
   
   call iotk_close_read (io_unit)
  end subroutine getgvec

end program irakurri_geak
