
include ../make.sys


PROGS = wanndrho_v0.x read_wann.x intw2W90_2.x bands_interpol.x
#ep_melements.x intw2W90.x intw2W90_2.x read_upf_tofile.x generate_irreducible_mesh.x print_wfc_test.x interpolatu_dV_ind_q.x phon_selfen_texture_q.x phon_selfen_texture_no_spin.x phon_selfen_texture_yes_spin.x phon_selfen_texture_q_indtot.x density_calculation.x transitions.x spin_texture.x spin_texture_spinor.x eliashberg_full_no_spin.x eliashberg_full.x lambda_full_nk_no_spin.x  print_bands_ekk_mathematica.x lambda_full_e_no_spin_graph.x lambda_full_graph.x selfen_direct.x selfen_el_full_nk.x eliashberg_full_no_spin_omega_omegap.x selfen_direct_nk.x imselfen2reselfen_eiguren_goiricelaya.x spectral_function_nk.x lambda_full_nk_texture.x ep_melements_K.x lambda_full_nk.x spin_energies.x phon_selfen_texture_q_indtot_AB.x phon_polarizations_Tl.x lambda_full_nk_no_spin_direct.x phonon_linewidths_nuq.x phonon_adiabatic_renormalization_nuq.x phonon_adiabatic_renormalization_dynmat.x lambda_full_ef_no_spin_direct.x lambda_full_irr_ef_no_spin_direct.x lambda_qe_nk_no_spin_direct.x lambda_qe_irr_ef_no_spin_direct.x

##haritz: $(PROGS_HARITZ)

debug: FFLAGS = $(FFLAGS_DEBUG)
debug: clean all

optim: FFLAGS = $(FFLAGS_OPTIM)
optim: clean all

all: $(PROGS_HARITZ) $(PROGS)


clean:
	rm -f *.x *.o *~ *.mod

include make.dep
