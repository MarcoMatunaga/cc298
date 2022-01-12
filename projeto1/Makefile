# Makefile created by mkmf 19.3.0


MK_TEMPLATE = mkmf_template
include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: a.out
mod_derivatives.o: ./mod_derivatives.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./mod_derivatives.f90
mod_diagonalization.o: ./mod_diagonalization.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./mod_diagonalization.f90
mod_flux_vector_splitting.o: ./mod_flux_vector_splitting.f90 mod_vars.o mod_vars_sw.o mod_diagonalization.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./mod_flux_vector_splitting.f90
mod_functions.o: ./mod_functions.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./mod_functions.f90
mod_output_routines.o: ./mod_output_routines.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./mod_output_routines.f90
mod_vars.o: ./mod_vars.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./mod_vars.f90
mod_vars_proj3.o: ./mod_vars_proj3.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./mod_vars_proj3.f90
mod_vars_sw.o: ./mod_vars_sw.f90 mod_vars.o mod_diagonalization.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./mod_vars_sw.f90
src_CFL.o: ./src_CFL.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_CFL.f90
src_artificial_dissipation.o: ./src_artificial_dissipation.f90 mod_vars.o mod_functions.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_artificial_dissipation.f90
src_blktriad.o: ./src_blktriad.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_blktriad.f90
src_blktriad_pc.o: ./src_blktriad_pc.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_blktriad_pc.f90
src_boundary_conditions_curv.o: ./src_boundary_conditions_curv.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_boundary_conditions_curv.f90
src_boundary_proj3.o: ./src_boundary_proj3.f90 mod_vars.o mod_vars_proj3.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_boundary_proj3.f90
src_euler_explicit.o: ./src_euler_explicit.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_euler_explicit.f90
src_fluxes_curvilinear.o: ./src_fluxes_curvilinear.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_fluxes_curvilinear.f90
src_implicit_beam_warming.o: ./src_implicit_beam_warming.f90 mod_vars.o mod_functions.o mod_diagonalization.o mod_derivatives.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_implicit_beam_warming.f90
src_implicit_beam_warming_penta.o: ./src_implicit_beam_warming_penta.f90 mod_vars.o mod_functions.o mod_diagonalization.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_implicit_beam_warming_penta.f90
src_initial_cond_curv.o: ./src_initial_cond_curv.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_initial_cond_curv.f90
src_initial_cond_proj3.o: ./src_initial_cond_proj3.f90 mod_vars.o mod_vars_proj3.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_initial_cond_proj3.f90
src_inv.o: ./src_inv.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_inv.f90
src_jacobian_eta.o: ./src_jacobian_eta.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_jacobian_eta.f90
src_jacobian_ksi.o: ./src_jacobian_ksi.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_jacobian_ksi.f90
src_mesh.o: ./src_mesh.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_mesh.f90
src_mesh_proj3.o: ./src_mesh_proj3.f90 mod_vars.o mod_vars_proj3.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_mesh_proj3.f90
src_metric_terms.o: ./src_metric_terms.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_metric_terms.f90
src_non_linear_dissipation.o: ./src_non_linear_dissipation.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_non_linear_dissipation.f90
src_penta_block.o: ./src_penta_block.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_penta_block.f90
src_projeto1.o: ./src_projeto1.f90 mod_vars.o mod_vars_proj3.o mod_vars_sw.o mod_output_routines.o mod_flux_vector_splitting.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_projeto1.f90
src_pulliam_chausse.o: ./src_pulliam_chausse.f90 mod_vars.o mod_diagonalization.o mod_functions.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_pulliam_chausse.f90
src_pulliam_chausse_block.o: ./src_pulliam_chausse_block.f90 mod_vars.o mod_diagonalization.o mod_functions.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_pulliam_chausse_block.f90
src_pulliam_chausse_block_penta.o: ./src_pulliam_chausse_block_penta.f90 mod_vars.o mod_diagonalization.o mod_functions.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_pulliam_chausse_block_penta.f90
src_residue.o: ./src_residue.f90 mod_vars.o mod_functions.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_residue.f90
src_residue_flux_vector_splitting_1st.o: ./src_residue_flux_vector_splitting_1st.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_residue_flux_vector_splitting_1st.f90
src_residue_flux_vector_splitting_2nd.o: ./src_residue_flux_vector_splitting_2nd.f90 mod_vars.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_residue_flux_vector_splitting_2nd.f90
src_sw_1st.o: ./src_sw_1st.f90 mod_vars.o mod_vars_sw.o mod_diagonalization.o mod_flux_vector_splitting.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_sw_1st.f90
src_sw_1st_a.o: ./src_sw_1st_a.f90 mod_vars.o mod_vars_sw.o mod_diagonalization.o mod_flux_vector_splitting.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_sw_1st_a.f90
src_sw_2nd.o: ./src_sw_2nd.f90 mod_vars.o mod_vars_sw.o mod_diagonalization.o mod_flux_vector_splitting.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_sw_2nd.f90
src_sw_2nd_a.o: ./src_sw_2nd_a.f90 mod_vars.o mod_vars_sw.o mod_diagonalization.o mod_flux_vector_splitting.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_sw_2nd_a.f90
src_thomas_pulliam_chausse.o: ./src_thomas_pulliam_chausse.f90
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./src_thomas_pulliam_chausse.f90
SRC = ./src_sw_2nd.f90 ./src_pulliam_chausse_block_penta.f90 ./src_mesh.f90 ./src_blktriad.f90 ./src_boundary_conditions_curv.f90 ./src_penta_block.f90 ./src_jacobian_ksi.f90 ./src_CFL.f90 ./src_mesh_proj3.f90 ./src_thomas_pulliam_chausse.f90 ./mod_output_routines.f90 ./src_blktriad_pc.f90 ./mod_vars_sw.f90 ./src_initial_cond_proj3.f90 ./src_implicit_beam_warming_penta.f90 ./mod_derivatives.f90 ./src_jacobian_eta.f90 ./src_sw_1st.f90 ./src_euler_explicit.f90 ./mod_diagonalization.f90 ./src_residue_flux_vector_splitting_2nd.f90 ./src_fluxes_curvilinear.f90 ./src_projeto1.f90 ./src_residue.f90 ./src_non_linear_dissipation.f90 ./src_sw_1st_a.f90 ./mod_vars.f90 ./src_implicit_beam_warming.f90 ./src_artificial_dissipation.f90 ./mod_functions.f90 ./mod_flux_vector_splitting.f90 ./mod_vars_proj3.f90 ./src_boundary_proj3.f90 ./src_residue_flux_vector_splitting_1st.f90 ./src_metric_terms.f90 ./src_inv.f90 ./src_sw_2nd_a.f90 ./src_pulliam_chausse.f90 ./src_pulliam_chausse_block.f90 ./src_initial_cond_curv.f90
OBJ = src_sw_2nd.o src_pulliam_chausse_block_penta.o src_mesh.o src_blktriad.o src_boundary_conditions_curv.o src_penta_block.o src_jacobian_ksi.o src_CFL.o src_mesh_proj3.o src_thomas_pulliam_chausse.o mod_output_routines.o src_blktriad_pc.o mod_vars_sw.o src_initial_cond_proj3.o src_implicit_beam_warming_penta.o mod_derivatives.o src_jacobian_eta.o src_sw_1st.o src_euler_explicit.o mod_diagonalization.o src_residue_flux_vector_splitting_2nd.o src_fluxes_curvilinear.o src_projeto1.o src_residue.o src_non_linear_dissipation.o src_sw_1st_a.o mod_vars.o src_implicit_beam_warming.o src_artificial_dissipation.o mod_functions.o mod_flux_vector_splitting.o mod_vars_proj3.o src_boundary_proj3.o src_residue_flux_vector_splitting_1st.o src_metric_terms.o src_inv.o src_sw_2nd_a.o src_pulliam_chausse.o src_pulliam_chausse_block.o src_initial_cond_curv.o
clean: neat
	-rm -f .a.out.cppdefs $(OBJ) a.out

neat:
	-rm -f $(TMPFILES)

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

a.out: $(OBJ)
	$(LD) $(OBJ) -o a.out  $(LDFLAGS)

