
# The directory where thw .mod files are stored. By default is empty and the files are stored in the same directory.


set(SOURCES
  all_wfcs.f90
  cfftnd.f90
  eliashberg.f90
  fft.f90
  force_linearity_G.f90
  gk_sort.f90
  init_pp.f90
  kinds.f90
  input_parameters.f90
  intw2wannier.f90
  intw_pseudo.f90
  intw_pseudo_local.f90
  intw_pseudo_non_local.f90
  lambda.f90
  matrix_elements.f90
  ph_module.f90
  reading.f90
  rotaxis.f90
  set_asr.f90
  sph_bes.f90
  symmetries.f90
  useful_constants.f90
  utility.f90
  ylmr2.f90
  fft_interp.f90
  splinelib.f90
  mcf_interpolacion.f90
  w90_setup.f90
  selfen.f90
)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files)

add_library(intw_modules STATIC ${SOURCES})

# With PRIVATE, the include directories are added only for this target.
# Using PUBLIC, the include directories will also be added for the targets that depend on intw_modules.
# As the utilities also need the QE and Wannier *.mod files, I could use PUBLIC, and this will add this include directories also to the utilities, which depend on intw_modules.
# But in this case I prefer to add them explicitly when the utilities targets are declared, that is why I use PRIVATE.



# And with INTERFACE, the include directories are only added to the targets that depend on intw_modules, but will not be added for the target intw_modules itself.
target_include_directories( intw_modules INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/mod_files)


install(TARGETS intw_modules
  ARCHIVE DESTINATION lib
)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
