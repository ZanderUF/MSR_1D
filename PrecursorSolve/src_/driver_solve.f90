!*****************************************************************************80
!
!! Main subroutine for the solve 

subroutine driver_solve ( )
  
USE parameters_fe
USE datainput_fe_M
 
implicit none

!---Dummy

!---Local
    integer  ::  status, i,j, n  
    character(len=80) :: current_path
    character(len=20) :: outfile_name
    character(len=20) :: steady_state_soln_file_name
    character(len=20) :: last_time_file_name
    character(len=20) :: power_soln_file_name
    character(len=7) :: input_file

    status = getcwd(current_path)
    outfile_name = 'outfile.txt'
    steady_state_soln_file_name = 'ss_soln.txt'
    last_time_file_name = 'last_t_soln.txt'
    power_soln_file_name = 'power_amp_soln.txt'
    input_file = 'input_t'

!---Open file for writing out debug information
    open (unit=outfile_unit, file=outfile_name,status='unknown',&
    	  form='formatted',position='asis')
    open (unit=soln_outfile_unit, file=steady_state_soln_file_name,&
    	  status='unknown',form='formatted',position='asis')
    open (unit=power_outfile_unit, file=power_soln_file_name,&
    	  status='unknown',form='formatted',position='asis')
    
!---Read in problem parameters here
    call datainput_fe(input_file)

!---Max dimension of the matrices to be computed, including solution vector
    matrix_length = 3*num_elem  

!---Allocate solution vector and global matrices
    allocate(precursor_soln_new(num_isotopes,num_delay_group,num_elem,nodes_per_elem),  &
             power_soln_new(num_elem,nodes_per_elem), &
             temperature_soln_new( num_elem,nodes_per_elem),  &
             density_soln_new( num_elem,nodes_per_elem),  &
             velocity_soln_new( num_elem,nodes_per_elem),  &
             precursor_soln_prev(num_isotopes,num_delay_group,num_elem,nodes_per_elem), &
             power_soln_prev(num_elem,nodes_per_elem), &
             temperature_soln_prev( num_elem,nodes_per_elem), &
             density_soln_prev( num_elem,nodes_per_elem),  &
             velocity_soln_prev( num_elem,nodes_per_elem), &
             spatial_power_fcn( num_elem, nodes_per_elem) )
    allocate(elem_vec_q_final(num_isotopes,num_delay_group,nodes_per_elem)) 
    allocate(elem_vol_int(num_elem,nodes_per_elem))
    allocate(precursor_soln_last_time( num_isotopes,num_delay_group,num_elem,nodes_per_elem),&
             power_soln_last_time(num_elem,nodes_per_elem) )

!---Test integral over volume of element only int_-1^1 f(x)
    elem_vol_int(:,:) = 0

!---Reactor properties
    mass_elem = 100.0/num_elem

!---Starting element for non fuel region, subtract off from end of domain
    non_fuel_start = num_elem - num_elem_external 

!---Initialize matrix entries 
    elem_matrix_A = 0.0

!---Create 1D mesh
    call mesh_creation

!---Steady state solve for temperature 
    call steady_state

    if( td_method_type == 0) then
        write(outfile_unit, fmt=('(a)')) 'Performing forward Euler time integration' 
        print *,'FORWARD'
        !---Transient solve forward Euler method 
        call transient_forward_euler
    end if
    
    if( td_method_type == 1) then
        write(outfile_unit, fmt=('(a)')) 'Performing backward Euler time integration'
        !---Transient solve backward Euler method
        call transient_backward_euler
    end if

end 
