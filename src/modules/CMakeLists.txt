
# The directory where thw .mod files are stored. By default is empty and the files are stored in the same directory.


set(SOURCES
  multiply_psi_by_dvKB.f90
  multiply_psi_by_vKB.f90
  add_m_nu_lambda_text.f90
  add_q_m_nu.f90
  allocate_nlpot.f90
  allocate_phq.f90
  all_wfcs.f90
  cfftnd.f90
  DM_q2Rp.f90
  DM_Rp2q.f90
  fft.f90
  force_linearity_G.f90
  gaussian.f90
  G_kqw2kqb.f90
  G_kqw2ReRpw.f90
  gk_sort.f90
  G_ReRpw2kqw.f90
  H_kb2Rew.f90
  H_Rew2kb.f90
  imselfen2reselfen_t0.f90
  init_pp.f90
  init_KB_projectors.f90
  dft_fftw.f90
  kinds.f90
  haritz.f90
  input_parameters.f90
  intw2wannier.f90
  intw_pseudo2.f90
  lambda_calculate.f90
  matrix_elements.f90
  ph_module.f90
  phq_init.f90
  reading.f90
  rotaxis.f90
  selfeliash_calculate.f90
  selfen_add_terms.f90
  set_asr.f90
  sph_bes.f90
  setlocq.f90
  sort.f90
  spin_kb2Rw.f90
  spin_orb.f90
  spin_Rw2kb.f90
  symmetries.f90
  useful_constants.f90
  utility.f90
  W90_to_intw.f90
  wigner_seitz.f90
  wigner_seitz_wannier.f90
  ylmr2.f90
  fft_interp.f90
  splinelib.f90
  mcf_interpolacion.f90
  w90_setup.f90
)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod_files)

add_library(intw_modules STATIC ${SOURCES})

# With PRIVATE, the include directories are added only for this target.
# Using PUBLIC, the include directories will also be added for the targets that depend on intw_modules.
# As the utilities also need the QE and Wannier *.mod files, I could use PUBLIC, and this will add this include directories also to the utilities, which depend on intw_modules.
# But in this case I prefer to add them explicitly when the utilities targets are declared, that is why I use PRIVATE.



target_link_libraries( intw_modules PUBLIC W_lib)
# And with INTERFACE, the include directories are only added to the targets that depend on intw_modules, but will not be added for the target intw_modules itself.
target_include_directories( intw_modules INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/mod_files)


install(TARGETS intw_modules
  ARCHIVE DESTINATION lib
)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
