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

    implicit none

!---
    if(DEBUG .eqv. .TRUE.) then
        write(outfile_unit, fmt='(a)'), 'In DEBUG mode, going to print out every matrix'
    end if

    write(outfile_unit, fmt='(a)'), ' '
    
    write(outfile_unit, fmt='(a)'), '0 - backward Euler | 1 - forward Euler'
    write(outfile_unit, fmt='(a,I2)'), 'Performing the time integration method: ',&
                                         td_method_type    
    
    if(step_flag .eqv. .TRUE.) then
       write(outfile_unit, fmt='(a)'), 'Performing a STEP perturbation'
    end if
    write(outfile_unit, fmt='(a)'), ' '
    if(ramp_flag .eqv. .TRUE.) then
        write(outfile_unit, fmt='(a)'), 'Performing a RAMP perturbation' 
    end if

    if(zag_flag .eqv. .TRUE.) then
        write(outfile_unit, fmt='(a)'), 'Performing ZIG-ZAG perturbation' 
    end if

    write(outfile_unit, fmt='(a,f10.7)'), 'Taking constant time steps of: ', delta_t

    write(outfile_unit, fmt='(a,f15.3)'), 'Integrating in time up to: ', tmax

    write(outfile_unit, fmt='(a,f15.3)'), 'Starting calculation at time: ', t_initial

    write(outfile_unit, fmt='(a,I7)'), 'Total number of elements in the &
                                             system: ', num_elem

    write(outfile_unit, fmt='(a,I6)'), 'Number of nodes per element: ', nodes_per_elem

    write(outfile_unit, fmt='(a,f10.3)'), 'Size of each element [cm]: ', elem_size

    write(outfile_unit, fmt='(a,I6)'), 'Inlet plenum starts at element: ',Fuel_Inlet_Start 

    write(outfile_unit, fmt='(a,I6)'), 'Main Core starts at element: ',Fuel_Core_Start 

    write(outfile_unit, fmt='(a,I6)'), 'Main Core ends at element: ', Fuel_Core_End 
    
    write(outfile_unit, fmt='(a,I6)'), 'Outlet plenum ends at element: ',Fuel_Outlet_End 
    
    write(outfile_unit, fmt='(a,f10.3)'), 'Cross sectional &
                                        area of the core [cm^2]: ', Area_Core
    
    write(outfile_unit, fmt='(a,f10.3)'), 'Cross sectional &
                                area of the piping [cm^2]: ', Area_Pipe

    write(outfile_unit, fmt='(a,f15.2)'), 'Mass flow rate [g/s]: ',mass_flow

    write(outfile_unit, fmt='(a,es15.3)'), 'Total initial starting power: ', total_power_initial

    write(outfile_unit, fmt='(a,I6)'), 'Maximum number of nonlinear iterations', max_nl_iter


    write(outfile_unit, fmt='(a,I6)'), 'Number of delayed neutron groups: ', num_delay_group

    write(outfile_unit, fmt='(a,I6)'), 'Number of fissionable isotopes: ', num_isotopes

    write(outfile_unit, fmt='(a,es20.7)'), 'Average nuetron generation time: ', gen_time

    write(outfile_unit, fmt='(a,f10.8)'), 'Reactivity input: ', reactivity_input

    write(outfile_unit, fmt='(a,f12.4)'), 'Save the spatial solution vectors at &
                            this interval: ', save_time_interval 
    write(outfile_unit, fmt='(a)'), ' '

end subroutine write_out_parms
