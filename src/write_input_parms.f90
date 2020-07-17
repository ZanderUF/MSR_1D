! Write out the input parameters and check they are valid 
!
!

subroutine write_out_parms()

     USE global_parameters_M
     USE time_info_M
     USE mesh_info_M
     USE material_info_M
     USE solution_vectors_M
     USE flags_M
     use Mod_SetupOutputFiles
    
    implicit none

!---

!---Start PARM block
    write(outfile_unit, fmt='(a)'), 'Writing PARM block input specifications'
    write(outfile_unit, fmt='(a)'),    '********************************************'
    write(outfile_unit, fmt='(a,L1)'), 'In DEBUG mode? (T/F)                        ',DEBUG
    write(outfile_unit, fmt='(a,L1)'), 'Reading power profile from DIF3D? (T/F)     ',Read_DIF3D 
write(outfile_unit, fmt='(a,es15.3)'), 'Total initial starting power:               ', total_power_initial
    write(outfile_unit, fmt='(a,I6)'), 'Maximum number of nonlinear iterations:     ', max_nl_iter
    write(outfile_unit, fmt='(a,I6)'), 'Number of delayed neutron groups:           ', num_delay_group
    write(outfile_unit, fmt='(a,I6)'), 'Number of fissionable isotopes:             ', num_isotopes
    write(outfile_unit, fmt='(a)'), ' '

!---Start TIME block
    write(outfile_unit, fmt='(a)'),   'Writing TIME block input specifications'
    write(outfile_unit, fmt='(a)'),    '********************************************'
    write(outfile_unit, fmt='(a,L1)'),'Doing time dependent calculation? (T/F)     ',time_solve
    write(outfile_unit, fmt='(a,I2)'),'Time integration method:                    ',&
                                         td_method_type    
    write(outfile_unit, fmt='(a)'),   'Integration method: 0 --> backward Euler ' 
    write(outfile_unit, fmt='(a)'),   'Integration method: 1 --> forward  Euler '
write(outfile_unit, fmt='(a,f10.7)'), 'Taking constant time steps of:              ', delta_t
write(outfile_unit, fmt='(a,f15.3)'), 'Integrating in time up to:                  ', tmax
write(outfile_unit, fmt='(a,f15.3)'), 'Starting calculation at time:               ', t_initial
write(outfile_unit, fmt='(a,f12.4)'), 'Save full solution vectors every:           ',save_time_interval 
    write(outfile_unit, fmt='(a)'), ' '

!---Start PERT block
    write(outfile_unit, fmt='(a)'),  'Writing PERT block input specifications'
    write(outfile_unit, fmt='(a)'),    '********************************************'
    write(outfile_unit, fmt='(a,I2)'),'Feedback method:                            ',feedback_method
    write(outfile_unit, fmt='(a)'),  'Feedback: 0 --> No feedback'
    write(outfile_unit, fmt='(a)'),  'Feedback: 1 -->            '
    write(outfile_unit, fmt='(a)'),  'Feedback: 2 -->            '
    write(outfile_unit, fmt='(a)'),  'Feedback: 3 -->            '
    write(outfile_unit, fmt='(a)'),  'Feedback: 4 -->            '
    write(outfile_unit, fmt='(a)'), ' '

    write(outfile_unit, fmt='(a,L1)'), 'Performing a STEP perturbation? (T/F)          ',step_flag
    write(outfile_unit, fmt='(a,L1)'), 'Performing a RAMP perturbation? (T/F)          ',ramp_flag
    write(outfile_unit, fmt='(a,L1)'), 'Performing ZIG-ZAG perturbation (T/F)          ',zag_flag
    write(outfile_unit, fmt='(a)'), ' '
    write(outfile_unit, fmt='(a,f5.3)'),'Step start time [s]:                           ',step_start_time
    write(outfile_unit, fmt='(a,f5.3)'),'Step end time [s]:                             ',step_end_time
    write(outfile_unit, fmt='(a,f5.3)'),'Ramp start time [s]:                           ',ramp_start_time
    write(outfile_unit, fmt='(a,f5.3)'),'Ramp end time [s]:                             ',ramp_end_time
    write(outfile_unit, fmt='(a)'), ' '
write(outfile_unit, fmt='(a,f10.8)'), 'Reactivity input:                               ', reactivity_input
    write(outfile_unit,fmt='(a,f8.3)'),'Time constant for pump change:                  ',time_constant
write(outfile_unit, fmt='(a,f8.3)'),  'Flow reduction percentage for pump change:      ',flow_reduction_percent
    
write(outfile_unit, fmt='(a,f15.2)'), 'Mass flow rate [g/s]:                           ',mass_flow
    
    write(outfile_unit, fmt='(a,es20.7)'), 'Average nuetron generation time:           ', gen_time
    write(outfile_unit, fmt='(a)'), ' '

!---Start MESH block
    write(outfile_unit, fmt='(a)'),   'Writing MESH block input specifications'
    write(outfile_unit, fmt='(a)'),    '********************************************'
    write(outfile_unit, fmt='(a,I7)'),'Total number of elements in the system:        ', num_elem
    write(outfile_unit, fmt='(a,I6)'),'Number of nodes per element:                   ', nodes_per_elem
write(outfile_unit, fmt='(a,f10.3)'), 'Size of each element [cm]:                     ', elem_size
   write(outfile_unit, fmt='(a,I6)'), 'Inlet plenum starts at element:                ',Fuel_Inlet_Start 
    write(outfile_unit, fmt='(a,I6)'),'Main Core starts at element:                   ',Fuel_Core_Start 
    write(outfile_unit, fmt='(a,I6)'),'Main Core ends at element:                     ', Fuel_Core_End 
    
    write(outfile_unit, fmt='(a,I6)'),'Outlet plenum ends at element:                 ',Fuel_Outlet_End 
    write(outfile_unit, fmt='(a,I6)'),'Heat Exchanger starts at element:              ',Heat_Exchanger_Start 
    write(outfile_unit, fmt='(a,I6)'),'Heat Exchanger ends at element:                ',Heat_Exchanger_End 
 
write(outfile_unit, fmt='(a,f10.3)'),'Cross sectional area of the piping [cm^2]:      ', Area_Pipe
write(outfile_unit, fmt='(a,f10.3)'),'Cross sectional area of the core   [cm^2]:      ', Area_Core
write(outfile_unit, fmt='(a,f10.3)'),'Cross sectional area of heat exchanger [cm^2]:  ', Area_Heat_Exchanger  
      write(outfile_unit, fmt='(a)'), ' '

end subroutine write_out_parms
