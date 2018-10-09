######################
#  Make file
#####################

SHELL = /bin/sh

# GNU
F90 = gfortran

#FFLAGS = -g -fbounds-check 

FFLAGS = -g -fbounds-check -fdefault-real-8
LIBDIR = ../../lapack-3.8.0/ 
LIB = -llapack -lblas

################################
msr1d : common_unit.o error_functions.o \
	vast_kind_param.o y0trns_M.o free_form.o \
	parameters_fe_M.o rdmesh_M.o \
	rdparm_fe_M.o rddelay_M.o\
	assemble_matrix_ss.o assemble_matrix_transient.o \
    cond_corr.o  \
	datainput_fe_M.o \
	density_corr.o   \
	spatial_matrices.o \
	heat_capacity_corr.o driver_solve.o     \
	inter_shape_fcns.o \
	mesh_creation.o steady_state.o \
    solve_soln_steady.o  \
    solve_power_backward_euler.o solve_power_forward_euler.o \
    solve_precursor_forward_euler.o solve_precursor_backward_euler.o \
    transient_forward_euler.o transient_backward_euler.o \
    solve_temperature.o solve_velocity.o \
    reactivity_feedback.o initialize.o \
    write_out_soln.o 
	$(F90) $(FFLAGS) main_driver.f90 *.o -L $(LIBDIR)$(LIB) 

main_driver.o : main_driver.f90
	$(F90) $(FFLAGS)  -c main_driver.f90  

common_unit.o : common_unit.f90
	$(F90) $(FFLAGS) -c common_unit.f90

error_functions.o : error_functions.f90 common_unit.o
	$(F90) $(FFLAGS) -c error_functions.f90

vast_kind_param.o : vast_kind_param.f90
	$(F90) $(FFLAGS) -c vast_kind_param.f90

y0trns_M.o : y0trns_M.f90
	$(F90) $(FFLAGS) -c y0trns_M.f90

free_form.o : free_form_C.f90 common_unit.o error_functions.o vast_kind_param.o y0trns_M.o
	$(F90) $(FFLAGS) -c free_form_C.f90

parameters_fe_M.o : parameters_fe_M.f90
	$(F90) $(FFLAGS) -c parameters_fe_M.f90

rdmesh_M.o : rdmesh_M.f90 free_form.o parameters_fe_M.o 
	$(F90) $(FFLAGS) -c rdmesh_M.f90 

rddelay_M.o : rdmesh_M.f90 free_form.o parameters_fe_M.o 
	$(F90) $(FFLAGS) -c rddelay_M.f90 

rdparm_fe_M.o : rdparm_fe_M.f90 free_form.o parameters_fe_M.o        
	$(F90) $(FFLAGS) -c rdparm_fe_M.f90

assemble_matrix_ss.o: assemble_matrix_ss.f90 parameters_fe_M.o
	$(F90) $(FFLAGS) -c assemble_matrix_ss.f90 

assemble_matrix_transient.o: assemble_matrix_transient.f90 parameters_fe_M.o
	$(F90) $(FFLAGS) -c assemble_matrix_transient.f90 

cond_corr.o: cond_corr.f90       
	$(F90) $(FFLAGS) -c cond_corr.f90

datainput_fe_M.o:  datainput_fe_M.f90 free_form.o rdparm_fe_M.o rddelay_M.o rdmesh_M.o parameters_fe_M.o    
	$(F90) $(FFLAGS) -c datainput_fe_M.f90 

density_corr.o:  density_corr.f90    
	$(F90) $(FFLAGS) -c density_corr.f90

spatial_matrices.o:  spatial_matrices.f90 parameters_fe_M.o datainput_fe_M.o  
	$(F90) $(FFLAGS) -c spatial_matrices.f90

driver_solve.o: driver_solve.f90  parameters_fe_M.o datainput_fe_M.o 
	$(F90) $(FFLAGS) -c driver_solve.f90

inter_shape_fcns.o: inter_shape_fcns.f90 parameters_fe_M.o
	$(F90) $(FFLAGS) -c inter_shape_fcns.f90

initialize.o: initialize.f90 parameters_fe_M.o
	$(F90) $(FFLAGS) -c initialize.f90

heat_capacity_corr.o: heat_capacity_corr.f90 parameters_fe_M.o
	$(F90) $(FFLAGS) -c heat_capacity_corr.f90

mesh_creation.o: mesh_creation.f90 parameters_fe_M.o   
	$(F90) $(FFLAGS) -c mesh_creation.f90

reactivity_feedback.o: reactivity_feedback.f90 parameters_fe_M.o 
	$(F90) $(FFLAGS) -c reactivity_feedback.f90

solve_soln_steady.o:   solve_soln_steady.f90 parameters_fe_M.o 
	$(F90) $(FFLAGS) -c solve_soln_steady.f90

solve_power_forward_euler.o:   solve_power_forward_euler.f90 parameters_fe_M.o 
	$(F90) $(FFLAGS) -c solve_power_forward_euler.f90

solve_power_backward_euler.o:   solve_power_backward_euler.f90 parameters_fe_M.o 
	$(F90) $(FFLAGS) -c solve_power_backward_euler.f90

solve_precursor_backward_euler.o:   solve_precursor_backward_euler.f90 parameters_fe_M.o 
	$(F90) $(FFLAGS) -c solve_precursor_backward_euler.f90

solve_precursor_forward_euler.o:   solve_precursor_forward_euler.f90 parameters_fe_M.o 
	$(F90) $(FFLAGS) -c solve_precursor_forward_euler.f90

steady_state.o : steady_state.f90 parameters_fe_M.o
	$(F90) $(FFLAGS) -c steady_state.f90

solve_temperature.o : solve_temperature.f90 parameters_fe_M.o
	$(F90) $(FFLAGS) -c solve_temperature.f90

solve_velocity.o : solve_velocity.f90 parameters_fe_M.o
	$(F90) $(FFLAGS) -c solve_velocity.f90

transient_forward_euler.o : transient_forward_euler.f90 parameters_fe_M.o
	$(F90) $(FFLAGS) -c transient_forward_euler.f90

transient_backward_euler.o : transient_backward_euler.f90 parameters_fe_M.o
	$(F90) $(FFLAGS) -c transient_backward_euler.f90
 
write_out_soln.o: write_out_soln.f90
	$(F90) $(FFLAGS) -c write_out_soln.f90


clean: 
	rm *.mod *.o
