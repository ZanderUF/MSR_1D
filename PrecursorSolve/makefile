######################
#  Make file
#####################

# GNU

F90 = gfortran

FFLAGS = -g -fbounds-check  -fdefault-real-8
LIBDIR = ../lapack-3.8.0/ 
LIB = -llapack -lblas

#################################
##
#	common_unit.o    \        
#	error_functions.o\ 
#	vast_kind_param.o \ 
#	y0trns_M.o \ 
#	free_form_C.o \  
#	parameters_fe.o  \   
#	rdmesh_M.o   \
#	assemble_matrix.o \ 
#	boundary_cond.o  \   
#	cond_corr.o  \
#	density_corr.o   \
#	element_matrix.o \
#	driver_solve.o     \
#	inter_shape_fcns.o \ 
#	kappa_corr.o     \
#	mesh_creation.o  \ 
#	solve_soln_steady.o     \
#	timestamp_M.o      
################################
main_driver : common_unit.o error_functions.o \
	vast_kind_param.o y0trns_M.o free_form.o \
	parameters_fe.o rdmesh_M.o \
	rdparm_fe_M.o rddelay_M.o\
	assemble_matrix.o assemble_matrix_transient.o \
    boundary_cond.o cond_corr.o  \
	datainput_fe_M.o \
	density_corr.o   \
	element_matrix.o \
	driver_solve.o     \
	inter_shape_fcns.o \
	kappa_corr.o     \
	mesh_creation.o steady_state.o \
    solve_soln_transient.o transient_solve_euler.o \
    solve_power_transient.o \
    solve_soln_steady.o     \
	initialize.o \
    write_out_soln.o \
    timestamp_M.o      
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

parameters_fe.o : parameters_fe.f90
	$(F90) $(FFLAGS) -c parameters_fe.f90

rdmesh_M.o : rdmesh_M.f90 free_form.o parameters_fe.o 
	$(F90) $(FFLAGS) -c rdmesh_M.f90 

rddelay_M.o : rdmesh_M.f90 free_form.o parameters_fe.o 
	$(F90) $(FFLAGS) -c rddelay_M.f90 

rdparm_fe_M.o : rdparm_fe_M.f90 free_form.o parameters_fe.o        
	$(F90) $(FFLAGS) -c rdparm_fe_M.f90

assemble_matrix.o: assemble_matrix.f90 parameters_fe.o
	$(F90) $(FFLAGS) -c assemble_matrix.f90 

assemble_matrix_transient.o: assemble_matrix_transient.f90 parameters_fe.o
	$(F90) $(FFLAGS) -c assemble_matrix_transient.f90 

boundary_cond.o: boundary_cond.f90 parameters_fe.o
	$(F90) $(FFLAGS) -c boundary_cond.f90   

cond_corr.o: cond_corr.f90       
	$(F90) $(FFLAGS) -c cond_corr.f90

datainput_fe_M.o:  datainput_fe_M.f90 free_form.o rdparm_fe_M.o rddelay_M.o rdmesh_M.o parameters_fe.o    
	$(F90) $(FFLAGS) -c datainput_fe_M.f90 

density_corr.o:  density_corr.f90    
	$(F90) $(FFLAGS) -c density_corr.f90

element_matrix.o:  element_matrix.f90 parameters_fe.o datainput_fe_M.o  
	$(F90) $(FFLAGS) -c element_matrix.f90

driver_solve.o: driver_solve.f90  parameters_fe.o datainput_fe_M.o 
	$(F90) $(FFLAGS) -c driver_solve.f90

inter_shape_fcns.o: inter_shape_fcns.f90 parameters_fe.o
	$(F90) $(FFLAGS) -c inter_shape_fcns.f90

initialize.o: initialize.f90 parameters_fe.o
	$(F90) $(FFLAGS) -c initialize.f90

kappa_corr.o: kappa_corr.f90      
	$(F90) $(FFLAGS) -c kappa_corr.f90

mesh_creation.o: mesh_creation.f90 parameters_fe.o   
	$(F90) $(FFLAGS) -c mesh_creation.f90

solve_soln_steady.o:   solve_soln_steady.f90 parameters_fe.o 
	$(F90) $(FFLAGS) -c solve_soln_steady.f90

solve_power_transient.o:   solve_power_transient.f90 parameters_fe.o 
	$(F90) $(FFLAGS) -c solve_power_transient.f90

solve_soln_transient.o:   solve_soln_steady.f90 parameters_fe.o 
	$(F90) $(FFLAGS) -c solve_soln_transient.f90

steady_state.o : steady_state.f90 parameters_fe.o
	$(F90) $(FFLAGS) -c steady_state.f90

transient_solve_euler.o : transient_solve_euler.f90 parameters_fe.o
	$(F90) $(FFLAGS) -c transient_solve_euler.f90
 
write_out_soln.o: write_out_soln.f90
	$(F90) $(FFLAGS) -c write_out_soln.f90

timestamp_M.o: timestamp_M.f90     
	$(F90) $(FFLAGS) -c timestamp_M.f90


clean: 
	rm *.mod *.o
