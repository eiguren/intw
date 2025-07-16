integer function read_input_test() result(r)

  use kinds, only: dp
  use intw_test_module, only: assert
  use intw_input_parameters, only: read_input, outdir, prefix, nk1, nk2, nk3

  implicit none

  logical :: read_status

  r = 0

  close(unit=5)
  open(unit=5,file='intw.in')

  read_status = .true.
  call read_input(read_status)

  if (read_status) then
    r = 1
    return
  endif

  write(6,*) "outdir:", trim(outdir)
  if (trim(outdir) /= "./qe/") then
    r = 1
    return
  endif

  write(6,*) "prefix:", trim(prefix)
  if (trim(prefix) /= "cu") then
    r = 1
    return
  endif

  write(6,*) "nk1 nk2 nk3:", nk1, nk2, nk3
  if (.not.assert((/nk1, nk2, nk3/),(/8, 8, 8/),0.0_dp)) then
    r = 1
    return
  endif


  return

end function
