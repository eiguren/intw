
# Check Fortran 90 support
if(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
  message(FATAL_ERROR "Fortran compiler does not support Fortran 90")
endif(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)


# Check GNU Fortran compiler version
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 9 )
    message(FATAL_ERROR "GNU Fortran compiler version >= 9 is required!")
  endif()
endif()


# Set Fortran Flags
include(SetFortranFlags)