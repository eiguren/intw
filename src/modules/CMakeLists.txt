
# The directory where thw .mod files are stored. By default is empty and the files are stored in the same directory.


set(SOURCES
  all_wfcs.f90
  cfftnd.f90
  fft.f90
  kinds.f90
  input_parameters.f90
  intw2wannier.f90
  intw_pseudo.f90
  intw_pseudo_local.f90
  intw_pseudo_non_local.f90
  matrix_elements.f90
  matrix_vector.f90
  ph_module.f90
  reading.f90
  symmetries.f90
  useful_constants.f90
  utility.f90
  fft_interp.f90
  w90_setup.f90
  ph_interpolate.f90
  version.f90
)

set_source_files_properties("version.f90" PROPERTIES COMPILE_DEFINITIONS -DVERSION_STRING="${VERSION_STRING}")


set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files)

add_library(intw_modules STATIC ${SOURCES})

# With PRIVATE, the include directories are added only for this target.
# Using PUBLIC, the include directories will also be added for the targets that depend on intw_modules.
# And with INTERFACE, the include directories are only added to the targets that depend on intw_modules, but will not be added for the target intw_modules itself.

target_include_directories(intw_modules INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/mod_files)

if(USE_OPENMP)
  target_link_libraries(intw_modules PUBLIC OpenMP::OpenMP_Fortran)
endif()



install(TARGETS intw_modules
  ARCHIVE DESTINATION lib
)


if(TRIFS)

  set(SOURCES_TRIFS
    geometry.f90
    isosurface.f90
    mesh_opt.f90
  )


  add_library(triFS_modules STATIC ${SOURCES_TRIFS})
  target_link_libraries(triFS_modules PUBLIC intw_modules)
  target_include_directories(triFS_modules INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/mod_files)
  if(USE_OPENMP)
    target_link_libraries(triFS_modules PUBLIC OpenMP::OpenMP_Fortran)
  endif()

  install(TARGETS intw_modules
    ARCHIVE DESTINATION lib
  )

endif()

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
