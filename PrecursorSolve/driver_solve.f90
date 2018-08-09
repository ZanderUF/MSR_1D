!*****************************************************************************80
!
!! Main subroutine for the solve 

subroutine driver_solve ( )
  
USE parameters_fe
USE datainput_fe_M
 
implicit none

!---Dummy

!---Local
    integer  ::  i,j, n   ! counter 

!---Open file for writing out debug information
    open (unit=outfile_unit, file="outfile.txt",status='unknown',form='formatted',position='asis')
!---Open file for writing out solution
    open (unit=soln_outfile_unit, file='ss_solution_file.txt',status='unknown',form='formatted',position='asis')
    open (unit=66, file='last_t_solution_file.txt',status='unknown',form='formatted',position='asis')
!---Read in problem parameters here
    call datainput_fe

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
             amplitude_fcn( num_elem, nodes_per_elem) )
    allocate(elem_vec_q_final(num_isotopes,num_delay_group,nodes_per_elem)) 
    allocate(elem_vol_int(num_elem,nodes_per_elem))

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

!---Transient solve Euler method 
!    call transient_solve_euler

end 
