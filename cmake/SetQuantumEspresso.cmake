#############################################################################
#
#############################################################################

# Maybe it is a good idea to ask the user directly for this directory, since
# it can have any name other than build.
set(QE_BUILD_DIR ${QE_HOME}/build)


set(QE_INCLUDE_DIRS
  ${QE_BUILD_DIR}/atomic/mod/qe_atomic
  ${QE_BUILD_DIR}/CPV/mod/qe_cpv
  ${QE_BUILD_DIR}/CPV/mod/qe_cpv_wfdd_exe
  ${QE_BUILD_DIR}/dft-d3/mod/qe_dftd3
  ${QE_BUILD_DIR}/EPW/mod/qe_epw
  ${QE_BUILD_DIR}/external/fox/common/mod/FoX_fsys
  ${QE_BUILD_DIR}/external/fox/dom/mod/FoX_fsys
  ${QE_BUILD_DIR}/external/fox/fsys/mod/FoX_fsys
  ${QE_BUILD_DIR}/external/fox/sax/mod/FoX_fsys
  ${QE_BUILD_DIR}/external/fox/utils/mod/FoX_fsys
  ${QE_BUILD_DIR}/external/fox/wxml/mod/FoX_fsys
  ${QE_BUILD_DIR}/external/mbd/src
  ${QE_BUILD_DIR}/external/mod/qe_devxlib
  ${QE_BUILD_DIR}/external/mod/qe_w90chk2chk_exe
  ${QE_BUILD_DIR}/external/mod/qe_wannier90
  ${QE_BUILD_DIR}/FFTXlib/mod/qe_fftx
  ${QE_BUILD_DIR}/GWW/mod/qe_gww
  ${QE_BUILD_DIR}/GWW/mod/qe_gww_bse
  ${QE_BUILD_DIR}/GWW/mod/qe_gww_pw4gww
  ${QE_BUILD_DIR}/GWW/mod/qe_gww_simple
  ${QE_BUILD_DIR}/GWW/mod/qe_gww_simplebse
  ${QE_BUILD_DIR}/GWW/mod/qe_gww_simpleip
  ${QE_BUILD_DIR}/HP/mod/qe_hp
  ${QE_BUILD_DIR}/KS_Solvers/mod/qe_kssolver_davidson
  ${QE_BUILD_DIR}/KS_Solvers/mod/qe_kssolver_davidsonrci
  ${QE_BUILD_DIR}/LAXlib/mod/qe_lax
  ${QE_BUILD_DIR}/LR_Modules/mod/qe_lr_modules
  ${QE_BUILD_DIR}/Modules/mod/qe_modules
  ${QE_BUILD_DIR}/NEB/mod/qe_neb
  ${QE_BUILD_DIR}/NEB/mod/qe_neb_pathinterpolation_exe
  ${QE_BUILD_DIR}/PHonon/mod/qe_phonon_alpha2f_exe
  ${QE_BUILD_DIR}/PHonon/mod/qe_phonon_gamma
  ${QE_BUILD_DIR}/PHonon/mod/qe_phonon_matdyn_exe
  ${QE_BUILD_DIR}/PHonon/mod/qe_phonon_ph
  ${QE_BUILD_DIR}/PP/mod/qe_pp
  ${QE_BUILD_DIR}/PP/mod/qe_pp_epsilon_exe
  ${QE_BUILD_DIR}/PP/mod/qe_pp_exe
  ${QE_BUILD_DIR}/PP/mod/qe_pp_fermiproj_exe
  ${QE_BUILD_DIR}/PP/mod/qe_pp_fermisurface_exe
  ${QE_BUILD_DIR}/PP/mod/qe_pp_pawplot_exe
  ${QE_BUILD_DIR}/PP/mod/qe_pp_pw2wannier90_exe
  ${QE_BUILD_DIR}/PP/mod/qe_pp_st_dos_exe
  ${QE_BUILD_DIR}/PWCOND/mod/qe_pwcond_exe
  ${QE_BUILD_DIR}/PW/mod/qe_pw
  ${QE_BUILD_DIR}/TDDFPT/mod/qe_tddfpt
  ${QE_BUILD_DIR}/upflib/mod/qe_upflib
  ${QE_BUILD_DIR}/upflib/mod/qe_upflib_upfconv_exe
  ${QE_BUILD_DIR}/UtilXlib/mod/qe_utilx
  ${QE_BUILD_DIR}/XClib/mod/qe_xclib
  ${QE_BUILD_DIR}/XSpectra/mod/qe_xspectra
  ${QE_BUILD_DIR}/XSpectra/mod/qe_xspectra_gipaw
)

set(QE_LIBS
  ${QE_BUILD_DIR}/lib/libqe_pw.a
  ${QE_BUILD_DIR}/lib/libqe_modules.a
  ${QE_BUILD_DIR}/lib/libqe_upflib.a
  ${QE_BUILD_DIR}/lib/libqe_pp.a
  ${QE_BUILD_DIR}/lib/libqe_fftx.a
  ${QE_BUILD_DIR}/lib/libqe_xclib.a
  ${QE_BUILD_DIR}/lib/libqe_kssolver_davidson.a
  ${QE_BUILD_DIR}/lib/libqe_kssolver_cg.a
  ${QE_BUILD_DIR}/lib/libqe_kssolver_ppcg.a
  ${QE_BUILD_DIR}/lib/libqe_kssolver_paro.a
  ${QE_BUILD_DIR}/lib/libqe_kssolver_dense.a
  ${QE_BUILD_DIR}/lib/libqe_dftd3.a
  ${QE_BUILD_DIR}/lib/libmbd.a
  ${QE_BUILD_DIR}/lib/libFoX_dom.a
  ${QE_BUILD_DIR}/lib/libFoX_sax.a
  ${QE_BUILD_DIR}/lib/libFoX_wxml.a
  ${QE_BUILD_DIR}/lib/libFoX_common.a
  ${QE_BUILD_DIR}/lib/libFoX_utils.a
  ${QE_BUILD_DIR}/lib/libFoX_fsys.a
  ${QE_BUILD_DIR}/lib/libqe_libbeef.a
  ${QE_BUILD_DIR}/lib/libqe_lax.a
  ${QE_BUILD_DIR}/lib/libqe_utilx.a
  ${QE_BUILD_DIR}/lib/libqe_clib.a
)

add_library(QE_lib INTERFACE IMPORTED)
target_include_directories( QE_lib INTERFACE ${QE_INCLUDE_DIRS} )
target_link_libraries( QE_lib INTERFACE ${QE_LIBS})
